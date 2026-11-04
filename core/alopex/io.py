"""原子文件、基础内容身份及原生进程流；仅清理本次拥有的临时项。"""
from __future__ import annotations

import csv
import datetime as dt
import hashlib
import io
import json
import os
import shutil
import subprocess
import uuid
from collections.abc import Mapping, Sequence
from pathlib import Path


def sha256_bytes(data: bytes) -> str:
    """返回内存字节串的 SHA-256 摘要。"""

    return hashlib.sha256(data).hexdigest()


def sha256_file(path: str | Path, *, chunk_size: int = 4 * 1024 * 1024) -> str:
    """流式计算单个文件的 SHA-256。"""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json_bytes(value: object) -> bytes:
    """用于稳定内容身份的 canonical JSON 编码（排序键、紧凑分隔符）。"""
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def stable_file_identity(
    path: str | Path, *, hash_content: bool = True
) -> dict[str, object]:
    """记录 path/stat 身份；hash_content 时附带内容 SHA-256，并在录制期间侦测文件变化。"""

    source = Path(path)
    before = source.stat()
    identity: dict[str, object] = {
        "path": str(source.resolve()),
        "size_bytes": before.st_size,
        "mtime_ns": before.st_mtime_ns,
    }
    if hash_content:
        identity["sha256"] = sha256_file(source)
        after = source.stat()
        # 等长替换可保留 mtime；inode/ctime 防止把旧哈希绑定到新文件。
        if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (
            after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns,
        ):
            raise RuntimeError(f"file changed while recording identity: {source}")
    return identity


def atomic_write_bytes(path: str | Path, data: bytes, *, mode: int | None = None) -> None:
    """在目标路径旁写本次独占的临时文件，再以一次原子 rename 发布；仅清理本次创建的文件。"""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp.{os.getpid()}.{uuid.uuid4().hex}")
    created = False
    try:
        with open(temporary, "xb", opener=lambda name, flags: os.open(
            name, flags, 0o600 if mode is not None else 0o666,
        )) as handle:
            created = True
            if mode is not None:
                # 在任何内容写入前收紧权限，避免私有清单在临时文件中短暂可读。
                os.fchmod(handle.fileno(), mode)
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    finally:
        if created:
            temporary.unlink(missing_ok=True)


def atomic_write_json(
    path: str | Path,
    payload: Mapping[str, object],
    *,
    mode: int | None = None,
) -> None:
    """以确定性序列化（排序键、ensure_ascii=False）写出 JSON 并原子发布。"""

    data = (
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    atomic_write_bytes(path, data, mode=mode)


def _read_json(path: Path, label: str) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid {label} JSON: {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{label} must contain a JSON object: {path}")
    return value


def _utc_now() -> str:
    return dt.datetime.now(dt.UTC).isoformat()


def _close_native_stream(process: subprocess.Popen) -> None:
    """关闭本次原生 reader；解析失败时限时终止并回收孩子，不覆盖原始解析异常。"""
    process.stdout.close()
    if process.poll() is None:
        process.terminate()
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()


def _safe_relative(text: object, *, label: str) -> Path:
    relative = Path(str(text))
    if (
        not relative.parts
        or relative.is_absolute()
        or ".." in relative.parts
        or relative == Path(".")
    ):
        raise ValueError(f"unsafe {label}: {text!r}")
    return relative


def _remove_path(path: Path) -> None:
    if not path.exists() and not path.is_symlink():
        return
    if path.is_dir() and not path.is_symlink():
        _make_writable(path)
        shutil.rmtree(path)
    else:
        path.unlink(missing_ok=True)


def _fsync_directory(path: Path) -> None:
    try:
        descriptor = os.open(path, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _validate_copy(source: Path, partial: Path, *, is_bam: bool, samtools: str) -> None:
    if source.stat().st_size != partial.stat().st_size:
        raise OSError(f"partial-copy size mismatch: {source} -> {partial}")
    if is_bam:
        subprocess.run(
            [samtools, "quickcheck", "-v", str(partial)],
            check=True,
            stdout=subprocess.DEVNULL,
        )
    if source.suffix == ".json":
        with partial.open("r", encoding="utf-8") as handle:
            json.load(handle)


def publish(
    pairs: list[tuple[Path, Path]],
    *,
    commit_last: Path | None,
    bam_destination: Path | None,
    samtools: str,
) -> None:
    """将文件组写入本次独占 partial 并校验；替换前撤销旧 commit_last，全部完成后最后提交它。"""
    destinations = [destination for _source, destination in pairs]
    if len(destinations) != len(set(destinations)):
        raise ValueError("duplicate destination in atomic publish bundle")
    if commit_last is not None and commit_last not in destinations:
        raise ValueError("--commit-last must match one --file destination")
    if bam_destination is not None and bam_destination not in destinations:
        raise ValueError("--bam-destination must match one --file destination")

    partials: dict[Path, Path] = {}
    try:
        for source, destination in pairs:
            if not source.is_file():
                raise FileNotFoundError(f"publish source is missing: {source}")
            destination.parent.mkdir(parents=True, exist_ok=True)


            partial = destination.with_name(
                f"{destination.name}.partial.{os.getpid()}.{uuid.uuid4().hex}"
            )
            with source.open("rb") as reader, partial.open("xb") as writer:
                partials[destination] = partial
                shutil.copyfileobj(reader, writer, length=8 * 1024 * 1024)
                writer.flush()
                os.fsync(writer.fileno())
            _validate_copy(
                source,
                partial,
                is_bam=bam_destination is not None and destination == bam_destination,
                samtools=samtools,
            )

        ordered = [destination for destination in destinations if destination != commit_last]
        if commit_last is not None:
            if commit_last.exists() or commit_last.is_symlink():
                commit_last.unlink()
            ordered.append(commit_last)

        for destination in ordered:
            os.replace(partials[destination], destination)
            _fsync_directory(destination.parent)
            del partials[destination]
    finally:
        for partial in partials.values():
            try:
                partial.unlink()
            except FileNotFoundError:
                pass


def render_tsv(
    rows: Sequence[Mapping[str, object]], fieldnames: Sequence[str]
) -> str:
    """按当前 schema 序列化 TSV；固定 Unix 换行与确定性输出。"""

    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=list(fieldnames),
        delimiter="\t",
        lineterminator="\n",
        extrasaction="raise",
    )
    writer.writeheader()
    supported = set(fieldnames)
    for row in rows:
        unsupported = sorted(set(row) - supported)
        if unsupported:
            raise ValueError(
                "Manifest row contains unsupported field(s): "
                + ", ".join(str(field) for field in unsupported)
            )
        writer.writerow(
            {
                field: "" if row.get(field) is None else str(row.get(field, ""))
                for field in fieldnames
            }
        )
    return buffer.getvalue()


def write_tsv_atomic(
    path: str | Path,
    rows: Sequence[Mapping[str, object]],
    fieldnames: Sequence[str],
    *,
    mode: int | None = None,
) -> None:
    """用本次独占的临时文件原子发布 manifest；失败保留旧文件并清理自己的临时文件。"""
    atomic_write_bytes(path, render_tsv(rows, fieldnames).encode("utf-8"), mode=mode)


def _make_tree_read_only(root: Path) -> None:


    for path in sorted(root.rglob("*"), key=lambda item: len(item.parts), reverse=True):
        os.chmod(path, 0o755 if path.is_dir() else 0o444)
    os.chmod(root, 0o755)


def _make_writable(path: Path) -> None:
    if path.is_symlink() or not path.is_dir():
        return
    os.chmod(path, 0o755)
    for child in path.rglob("*"):
        if child.is_symlink() or not child.is_dir():
            continue
        try:
            os.chmod(child, 0o755)
        except FileNotFoundError:
            pass


def _remove_stage(stage: Path) -> None:
    if stage.exists():
        _make_writable(stage)
        shutil.rmtree(stage)
