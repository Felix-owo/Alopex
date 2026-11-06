"""环境与reference构建事务、Doctor固定路线和原始输入自检。"""
from __future__ import annotations

import sys

_LOADED_MODULE_CODE = sys._getframe().f_code

import csv
import datetime as dt
import gzip
import io
import json
import os
import re
import shutil
import subprocess
import uuid
from collections.abc import Iterable, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

from .config import (
    CONFIG_TEMPLATE_YAML,
    ENVIRONMENT_IDENTITY_SCHEMA,
    IDENTITY_NAME,
    MANAGED_CONTIGS,
    MANAGED_SPECIES,
    SPIKE_CONTIGS,
    _PRIMARY_CONTIGS,
    _PROTOCOL_FIXED_CONFIG,
    _fai_contigs,
    _read_reference_binding,
    _reference_file_record,
    _reference_file_state,
    apply_config_defaults,
    envs_section,
    envs_section_sha256,
    finalize_project_config,
    load_project_config,
    managed_alignment_contigs,
    primary_analysis_contigs,
    read_environment_identity,
    reference_binding_path,
    reference_component_files,
    reference_sequence_identity,
    validate_downstream_reference,
    validate_reference_component,
)
from .io import (
    atomic_write_bytes,
    atomic_write_json,
    sha256_file,
)
from .paths import PIPELINE_ROOT
from .runtime import effective_cpu_count


def begin_reference_build(reference: Path, output: Path, annotation: Path | None = None) -> None:
    """工具启动前绑定实际FASTA/FAI及可选GTF，供完成边界复核，非额外运行gate。"""
    ref = reference.resolve()
    payload = {"reference": reference_sequence_identity(ref), "fasta_state": list(_reference_file_state(ref))}
    if annotation is not None:
        payload["annotation"] = {"path": str(annotation.resolve()), **_reference_file_record(annotation)}
    atomic_write_json(output, payload)


def seal_reference_component(reference: Path, kind: str, *, build_input: Path | None = None,
                             root: Path | None = None) -> None:
    """在现有构建事务内、完整文件组安装后写最后标记，失败由该事务恢复。"""
    identity = reference_sequence_identity(reference) if kind != "fasta" else None
    if kind != "fasta":
        if build_input is None:
            raise ValueError("reference producer input binding is required")
        before = json.loads(build_input.read_text())
        if before["reference"] != identity or before["fasta_state"] != list(_reference_file_state(reference)):
            raise ValueError("reference changed during component build")
    receipt = reference_binding_path(reference, kind, root=root)
    # fasta重建时旧receipt已由组装事务撤销；只在此边界计算每个新工件的一次完整hash。
    files = {os.path.relpath(p, receipt.parent): _reference_file_record(p)
             for p in reference_component_files(reference, kind, root)}
    if kind.startswith("bismark-"):
        genome = files[os.path.relpath((root or Path(str(reference) + ".bismark")) / "genome.fa", receipt.parent)]
        if genome["sha256"] != identity["fasta"]["sha256"]:
            raise ValueError("Bismark genome.fa differs from source FASTA")
    atomic_write_json(receipt, {"schema_version": 1, "kind": kind, "reference": identity, "files": files})


def seal_downstream_reference_bindings(reference: Path, annotation: Path, tss: Path, cpg: Path,
                                       build_input: Path) -> None:
    before = json.loads(build_input.read_text())
    identity = reference_sequence_identity(reference)
    if before["reference"] != identity or before["fasta_state"] != list(_reference_file_state(reference)):
        raise ValueError("reference changed during downstream reference build")
    annotation_record = _reference_file_record(annotation, before["annotation"])
    for role, path in (("tss", tss), ("single-cpg", cpg)):
        receipt = Path(str(path) + ".reference.json")
        atomic_write_json(receipt, {"schema_version": 1, "kind": role, "reference": identity,
                                   "annotation": annotation_record,
                                   "files": {path.name: _reference_file_record(path)}})


def relocate_downstream_bindings(reference: Path, tss: Path, cpg: Path) -> None:
    """只在持锁Doctor自身rename后刷新ctime，不重复哈希已验证的同inode工件。

    输入/结果在事务期间禁止外部原位改写；该边界与发布硬链接的ctime处理一致，
    不承诺识别任意并发、同inode等长且恢复mtime的外部篡改。
    """
    identity = reference_sequence_identity(reference)
    for path, kind in ((tss, "tss"), (cpg, "single-cpg")):
        receipt = Path(str(path) + ".reference.json")
        payload = _read_reference_binding(receipt, kind)
        if payload["reference"] != identity or set(payload["files"]) != {path.name}:
            raise ValueError("downstream binding changed during publication")
        record = payload["files"][path.name]
        current = list(_reference_file_state(path))
        if current[:4] != record["state"][:4]:
            raise ValueError(f"downstream artifact changed across owned rename: {path}")
        record["state"] = current
        atomic_write_json(receipt, payload)


def import_reference_binding(reference: Path, kind: str, *, root: Path | None = None,
                             artifact: Path | None = None, annotation: Path | None = None) -> None:
    """接收已封存组件：全内容核对成功后仅刷新目标文件系统状态，不重签科学身份。

    必须已有 producer 回执；缺失、错代或任何 digest 不符均拒绝，不能据此批准旧索引。
    跨平台索引格式是否可用另由目标平台的实际小数据比对验收。
    """
    if artifact is None:
        receipt = reference_binding_path(reference, kind, root=root)
        payload = validate_reference_component(reference, kind, root=root)
    else:
        receipt = Path(str(artifact) + ".reference.json")
        payload = validate_downstream_reference(reference, artifact, kind, annotation)
    records = {relative: _reference_file_record(receipt.parent / relative, record)
               for relative, record in payload["files"].items()}
    if artifact is not None and annotation is not None and annotation.is_file():
        payload["annotation"] = _reference_file_record(annotation, payload["annotation"])
    payload["files"] = records
    atomic_write_json(receipt, payload)


def write_envs_section_spec(envs_path: str | Path, role: str, destination: str | Path) -> Path:
    """导出原生安装规格；默认频道由 Doctor 的空 CONDA_CHANNELS 禁用。"""
    from ruamel.yaml import YAML

    section = dict(envs_section(envs_path, role))
    channels = section.get("channels")
    if not isinstance(channels, list) or not channels or any(
        not isinstance(channel, str) or not channel.strip() for channel in channels
    ):
        raise ValueError(f"environment role {role!r} requires explicit channels")
    # nodefaults 是频道策略，不是下载地址；Conda 的新版 YAML 入口会把它当作频道。
    section["channels"] = [channel for channel in channels if channel != "nodefaults"]
    if not section["channels"]:
        raise ValueError(f"environment role {role!r} has no package channel")
    target = Path(destination)
    parser = YAML(typ="safe", pure=True)
    parser.default_flow_style = False
    parser.sort_base_mapping_type_on_output = False
    stream = io.StringIO()
    parser.dump(section, stream)
    atomic_write_bytes(target, stream.getvalue().encode("utf-8"))
    return target


def record_environment_identity(
    *,
    prefix: str | Path,
    role: str,
    spec_name: str,
    spec_sha256: str,
) -> Mapping[str, object]:
    """为一个 Conda prefix 记录最小的受管环境契约（spec 身份为逻辑名+section SHA）。"""

    prefix_path = Path(prefix).expanduser().resolve()
    if not prefix_path.is_dir():
        raise FileNotFoundError(f"environment prefix is missing: {prefix_path}")
    if not str(spec_sha256).strip():
        raise ValueError("spec_sha256 must be a non-empty digest")

    payload: dict[str, object] = {
        "schema_version": ENVIRONMENT_IDENTITY_SCHEMA,
        "complete": True,
        "role": str(role),
        "prefix": str(prefix_path),
        "specification": {
            "path": str(spec_name),
            "sha256": str(spec_sha256),
        },
    }
    atomic_write_json(prefix_path / IDENTITY_NAME, payload, mode=0o444)
    return payload


RELEASE_SCHEMA = 5


RELEASE_NAME = ".dna_pipeline_release.json"


RELEASE_ROLES = ("control", "biscuit", "bismark", "rastair", "notebook")


def finalize_release(release_s: str, release_id: str, expected_subdir: str) -> None:
    """为一个完整构建的 release 写出 0o444 的 release manifest（提交前最后一步）。"""
    release_path = Path(release_s).expanduser().absolute()
    if release_path.name != release_id:
        raise ValueError(f"release id/path mismatch: {release_id} != {release_path.name}")
    release = release_path.resolve(strict=True)
    records: dict[str, dict[str, object]] = {}
    for role in RELEASE_ROLES:
        role_path = release_path / role
        prefix = role_path.resolve(strict=True)
        if release.parent not in prefix.parents:
            raise ValueError(f"release role escaped managed release store: {role} -> {prefix}")
        identity = read_environment_identity(prefix)
        if identity.get("role") != role:
            raise ValueError(f"environment role mismatch: {role}")
        records[role] = {
            "relative_path": role,
            "resolved_prefix": str(prefix),
            "reused": role_path.is_symlink(),
        }
    payload: dict[str, object] = {
        "schema_version": RELEASE_SCHEMA,
        "complete": True,
        "release_id": release_id,
        "created_at_utc": dt.datetime.now(dt.UTC).isoformat(),
        "release_path": str(release),
        "platform": {"conda_subdir": expected_subdir},
        "roles": records,
    }
    atomic_write_json(release / RELEASE_NAME, payload, mode=0o444)


def publish_release(current_s: str, target_s: str) -> None:
    """以临时 symlink + 一次 os.replace 原子切换 conda/current（唯一 COMMIT 点）。"""
    current = Path(current_s).expanduser().absolute()
    target = Path(target_s).expanduser().resolve(strict=True)
    if not target.is_dir():
        raise FileNotFoundError(f"release directory missing: {target}")
    if current.exists() and not current.is_symlink():
        raise ValueError(f"current path is not a symlink: {current}")
    current.parent.mkdir(parents=True, exist_ok=True)
    temporary = current.with_name(f".{current.name}.publish.{os.getpid()}.{uuid.uuid4().hex}")
    created = False
    try:
        os.symlink(str(target), temporary)
        created = True
        os.replace(temporary, current)
    finally:
        if created:
            temporary.unlink(missing_ok=True)


def _check_installed_environment(prefix: Path, role: str, expected_subdir: str, envs_path: str) -> None:
    """只读核对已安装包；不导入目标环境的科学、绘图或 Notebook 模块。"""
    from conda.models.match_spec import MatchSpec
    from conda.models.version import VersionSpec
    from importlib.metadata import Distribution
    from packaging.markers import default_environment
    from packaging.requirements import Requirement

    records: dict[str, dict[str, Any]] = {}
    native = False
    # 包记录在共享文件系统上逐个同步读取会累积 I/O 等待；有界重叠读取，
    # 按原枚举顺序消费全部记录，保留坏 JSON、重复包和未请求包的平台拒绝。
    with ThreadPoolExecutor(max_workers=min(4, effective_cpu_count())) as reader:
        loaded = reader.map(
            lambda path: json.loads(path.read_text(encoding="utf-8")),
            (prefix / "conda-meta").glob("*.json"),
        )
        for record in loaded:
            name = record["name"]
            if name in records:
                raise ValueError(f"{role}: duplicate installed package {name}")
            records[name] = record
            subdir = record.get("subdir")
            if subdir and subdir != "noarch":
                if subdir != expected_subdir:
                    raise ValueError(f"{role}: installed package {name} has platform {subdir}")
                native = True
    if not native:
        raise ValueError(f"{role}: no native package records")
    section = envs_section(envs_path, role)
    python_version = records.get("python", {}).get("version", "")
    python_short = ".".join(python_version.split(".")[:2])
    site = prefix / "lib" / f"python{python_short}" / "site-packages"
    marker_environment = default_environment()
    marker_environment.update(python_version=python_short, python_full_version=python_version)
    normalize = lambda value: re.sub(r"[-_.]+", "-", value).lower()
    # Python 包可能由 pip 安装或覆盖；有效版本来自 dist-info，不由残留 Conda 记录代替。
    python_distributions = {name: name for name in (
        "anndata", "ipykernel", "jupyterlab", "numba", "numpy", "packaging", "pandas", "polars",
        "pyyaml", "scikit-learn", "seaborn", "tqdm", "jsonschema", "snakemake",
        "snakemake-executor-plugin-slurm", "pip", "pysam", "cutadapt", "multiqc", "awscli",
        "leidenalg", "cykhash", "umi_tools",
    )}
    python_distributions.update({"matplotlib-base": "matplotlib", "python-igraph": "igraph"})
    def distribution_for(name):
        candidates = list(Distribution.discover(name=name, path=[str(site)]))
        if len(candidates) != 1:
            raise ValueError(f"{role}: missing or ambiguous Python metadata for {name}")
        distribution = candidates[0]
        if Path(distribution.locate_file("")).resolve() != site.resolve():
            raise ValueError(f"{role}: Python metadata escaped the managed prefix for {name}")
        if normalize(distribution.metadata["Name"] or "") != normalize(name):
            raise ValueError(f"{role}: Python distribution identity mismatch for {name}")
        return distribution
    for dependency in section["dependencies"]:
        if isinstance(dependency, str):
            spec = MatchSpec(dependency)
            record = records.get(spec.name)
            distribution_name = python_distributions.get(spec.name)
            if distribution_name is not None:
                distribution = distribution_for(distribution_name)
                constraint = spec.get_raw_value("version")
                if constraint and not VersionSpec(constraint).match(distribution.version):
                    raise ValueError(f"{role}: effective Python package does not satisfy {dependency}")
                if (spec.get_raw_value("build") is not None or spec.get_raw_value("channel") is not None):
                    if record is None or not spec.match(record):
                        raise ValueError(f"{role}: installed artifact does not satisfy {dependency}")
            elif record is None or not spec.match(record):
                raise ValueError(f"{role}: installed package does not satisfy {dependency}")
        else:
            for item in dependency["pip"]:
                requirement = Requirement(item)
                if requirement.marker and not requirement.marker.evaluate(marker_environment):
                    continue
                distribution = distribution_for(requirement.name)
                if distribution.version not in requirement.specifier:
                    raise ValueError(f"{role}: installed pip package does not satisfy {item}")


def validate_release(current_s: str, expected_subdir: str, envs_path: str, *, installed: bool = False) -> None:
    """只读校验 current 指针、release manifest、role prefix 与 envs section SHA 契约。"""
    current = Path(current_s).expanduser().absolute()
    if not current.is_symlink():
        raise ValueError(f"current is not a symlink: {current}")
    release = current.resolve(strict=True)
    payload = json.loads((release / RELEASE_NAME).read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or payload.get("schema_version") != RELEASE_SCHEMA or payload.get("complete") is not True:
        raise ValueError("release identity is incomplete or unsupported")
    if payload.get("release_id") != release.name:
        raise ValueError("release id/path mismatch")
    if Path(str(payload.get("release_path", ""))).resolve() != release:
        raise ValueError("release path drift")
    platform = payload.get("platform") or {}
    if platform.get("conda_subdir") != expected_subdir:
        raise ValueError(f"release platform mismatch: {platform.get('conda_subdir')} != {expected_subdir}")
    records = payload.get("roles")
    if not isinstance(records, dict) or set(records) != set(RELEASE_ROLES):
        raise ValueError("release roles are invalid")

    for role in RELEASE_ROLES:
        role_path = release / role
        prefix = role_path.resolve(strict=True)
        if release.parent not in prefix.parents:
            raise ValueError(f"{role} prefix escaped managed release store")
        record = records[role]
        if not isinstance(record, dict) or record.get("relative_path") != role:
            raise ValueError(f"{role} release record is invalid")
        recorded_resolved = record.get("resolved_prefix")
        # producer 已保存解析后的绝对路径；再次 resolve 会掩盖原目录被换成链接的漂移。
        if (not isinstance(recorded_resolved, str)
                or not Path(recorded_resolved).is_absolute()
                or Path(recorded_resolved) != prefix):
            raise ValueError(f"{role} resolved prefix drift")
        environment = read_environment_identity(prefix)
        if environment.get("role") != role:
            raise ValueError(f"{role} environment identity role mismatch")
        recorded = (environment.get("specification") or {}).get("sha256")
        if recorded != envs_section_sha256(envs_path, role):
            raise ValueError(f"{role} specification changed; rebuild required")
        if installed:
            _check_installed_environment(prefix, role, expected_subdir, envs_path)


def environment_is_reusable(prefix_s: str, role: str, envs_path: str) -> None:
    """判断既有 prefix 可否复用：身份记录存在、role 匹配且 envs section SHA 未变。"""
    environment = read_environment_identity(prefix_s)
    if environment.get("role") != role:
        raise ValueError(f"{role} environment identity role mismatch")
    if (environment.get("specification") or {}).get("sha256") != envs_section_sha256(
        envs_path, role
    ):
        raise ValueError(f"{role} specification changed")


def _atomic_text_path(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    return path.with_name(f".{path.name}.tmp.{os.getpid()}.{uuid.uuid4().hex}")


def _gtf_attributes(raw: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in raw.split(";"):
        key, separator, value = item.strip().partition(" ")
        if separator:
            attributes[key] = value.strip().strip('"')
    return attributes


def _annotation_release_signature(species: str) -> tuple[str, str]:
    if species == "hg38":
        return "GRCh38", "version 48"
    if species == "mm10":
        return "GRCm38", "version M25"
    raise ValueError(f"unsupported standard species: {species}")


def build_downstream_references(
    *,
    species: str,
    fasta: str | Path,
    annotation_gtf: str | Path,
    tss_bed: str | Path,
    single_cpg_bed: str | Path,
) -> Mapping[str, object]:
    """构建 TSS bins 与 SnapATAC2 兼容的 single-CpG BED。

    CpG BED 为 BED3、0-based half-open，每条 interval 覆盖 CpG 二核苷酸 ``[C, G+1)``，
    因此两条链上报出的 base-resolution 值都会落入同一个 SnapATAC2 peak_file feature。
    """

    reference = Path(fasta).expanduser().resolve()
    fai = Path(f"{reference}.fai")
    annotation = Path(annotation_gtf).expanduser().resolve()
    tss = Path(tss_bed).expanduser().resolve()
    cpg = Path(single_cpg_bed).expanduser().resolve()
    if not reference.is_file() or reference.stat().st_size == 0:
        raise FileNotFoundError(f"reference FASTA is missing or empty: {reference}")
    if not fai.is_file() or fai.stat().st_size == 0:
        raise FileNotFoundError(f"reference FAI is missing or empty: {fai}")
    if not annotation.is_file() or annotation.stat().st_size == 0:
        raise FileNotFoundError(
            f"reference annotation GTF is missing or empty: {annotation}"
        )

    contigs = _fai_contigs(fai)
    host_contigs = [
        (name, length)
        for name, length in contigs
        if name in primary_analysis_contigs(species)
    ]
    contig_lengths = dict(host_contigs)
    expected_build, expected_release = _annotation_release_signature(species)
    header_lines: list[str] = []
    tss_tmp = _atomic_text_path(tss)
    cpg_tmp = _atomic_text_path(cpg)
    binding_input = _atomic_text_path(tss.with_suffix(".build-input.json"))
    begin_reference_build(reference, binding_input, annotation)
    tss_count = 0
    protein_coding_gene_count = 0
    seen_intervals: set[tuple[str, int, int]] = set()

    try:
        with annotation.open("rt", encoding="utf-8") as source, tss_tmp.open(
            "wt", encoding="utf-8", newline=""
        ) as output:
            for line_number, raw in enumerate(source, start=1):
                if raw.startswith("#"):
                    if len(header_lines) < 100:
                        header_lines.append(raw.rstrip("\n"))
                    continue
                fields = raw.rstrip("\n").split("\t")
                if len(fields) != 9:
                    raise ValueError(
                        f"invalid GTF field count on line {line_number}: {annotation}"
                    )
                (
                    chrom,
                    _source,
                    feature,
                    start_s,
                    end_s,
                    _score,
                    strand,
                    _frame,
                    attrs_s,
                ) = fields
                if (
                    feature != "gene"
                    or chrom not in contig_lengths
                    or strand not in {"+", "-"}
                ):
                    continue
                attrs = _gtf_attributes(attrs_s)
                if attrs.get("gene_type") != "protein_coding":
                    continue
                gene_id = attrs.get("gene_id", "").strip()
                if not gene_id:
                    raise ValueError(
                        "protein-coding GTF gene lacks gene_id on line "
                        f"{line_number}: {annotation}"
                    )
                try:
                    gene_start = int(start_s) - 1
                    gene_end = int(end_s)
                except ValueError as exc:
                    raise ValueError(
                        f"invalid GTF coordinate on line {line_number}: {annotation}"
                    ) from exc
                if (
                    gene_start < 0
                    or gene_end <= gene_start
                    or gene_end > contig_lengths[chrom]
                ):
                    raise ValueError(
                        "GTF gene exceeds canonical contig bounds on line "
                        f"{line_number}: {annotation}"
                    )
                protein_coding_gene_count += 1
                transcription_start = gene_start if strand == "+" else gene_end
                for offset in range(-2000, 2000, 20):
                    if strand == "+":
                        bin_start = transcription_start + offset
                        bin_end = bin_start + 20
                    else:
                        bin_end = transcription_start - offset
                        bin_start = bin_end - 20
                    if bin_start < 0 or bin_end > contig_lengths[chrom]:
                        continue
                    interval = (chrom, bin_start, bin_end)
                    if interval in seen_intervals:
                        continue
                    seen_intervals.add(interval)
                    region = f"{chrom}:{bin_start}-{bin_end}"
                    output.write(
                        f"{chrom}\t{bin_start}\t{bin_end}\t{region}\t.\t{strand}\t"
                        f"{gene_id}\t{region}\t{offset}\n"
                    )
                    tss_count += 1

        header = "\n".join(header_lines)
        if expected_build not in header or expected_release not in header:
            raise ValueError(
                f"{species} GTF must be GENCODE {expected_release} on "
                f"{expected_build}: {annotation}"
            )
        if protein_coding_gene_count == 0 or tss_count == 0:
            raise ValueError(f"GTF produced no protein-coding TSS bins: {annotation}")

        cpg_count = 0
        selected = {name for name, _length in host_contigs}
        current = ""
        current_position = 0
        previous_base = ""
        cpg_tmp.parent.mkdir(parents=True, exist_ok=True)
        with reference.open("rt", encoding="utf-8") as source, cpg_tmp.open(
            "wb"
        ) as raw_output:
            with gzip.GzipFile(
                filename="", mode="wb", fileobj=raw_output, mtime=0
            ) as compressed:
                with io.TextIOWrapper(
                    compressed, encoding="utf-8", newline=""
                ) as output:
                    for raw in source:
                        if raw.startswith(">"):
                            current = raw[1:].split(None, 1)[0]
                            current_position = 0
                            previous_base = ""
                            continue
                        sequence = "".join(raw.split()).upper()
                        if not sequence or current not in selected:
                            continue
                        combined = previous_base + sequence
                        combined_start = current_position - (1 if previous_base else 0)
                        match = combined.find("CG")
                        while match >= 0:
                            start = combined_start + match
                            end = start + 2
                            output.write(f"{current}\t{start}\t{end}\n")
                            cpg_count += 1
                            match = combined.find("CG", match + 2)
                        previous_base = sequence[-1]
                        current_position += len(sequence)
        if cpg_count == 0:
            raise ValueError(f"reference produced no CpG sites: {reference}")

        # 换成员前先验证构建输入仍同代；外层Doctor持有可恢复发布事务。
        before = json.loads(binding_input.read_text())
        if before["reference"] != reference_sequence_identity(reference):
            raise ValueError("reference changed during downstream build")
        os.replace(tss_tmp, tss)
        os.replace(cpg_tmp, cpg)
        seal_downstream_reference_bindings(reference, annotation, tss, cpg, binding_input)
    finally:
        for temporary in (tss_tmp, cpg_tmp, binding_input):
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass

    return {
        "species": species,
        "protein_coding_gene_count": protein_coding_gene_count,
        "tss_bin_count": tss_count,
        "single_cpg_count": cpg_count,
        "tss_bed": str(tss),
        "single_cpg_bed": str(cpg),
    }


def append_spikein_contigs(lambda_fasta: Path, puc19_fasta: Path, target_fasta: Path) -> None:
    """把 lambda 与 pUC19 spike-in contig 追加到临时 assembled FASTA 末尾。"""

    def spike_sequence(path: Path) -> str:
        data = "".join(
            line.strip()
            for line in path.read_text().splitlines()
            if line and not line.startswith(">")
        ).upper()
        if not data:
            raise ValueError(f"empty spike FASTA: {path}")
        return data

    with target_fasta.open("a") as handle:
        for name, path in (("lambda", lambda_fasta), ("pUC19", puc19_fasta)):
            sequence = spike_sequence(path)
            handle.write(f">{name}\n")
            for offset in range(0, len(sequence), 60):
                handle.write(sequence[offset : offset + 60] + "\n")


# NCBI GRCh38 analysis-set CEN/PAR/WGS 与 GRC exclusions 的保留序列子集。
# 原始坐标/内容SHA及稳定来源见 ANALYSIS_PARAMETERS；BED0/half-open，长度不变。
_ANALYSIS_MASK_INTERVALS = {
    "hg38": {
        "chr5": (
            (47309184, 49591369),
        ),
        "chr9": (
            (70701801, 70719794),
        ),
        "chr13": (
            (86252979, 86269524),
            (111753855, 111793441),
        ),
        "chr14": (
            (16000000, 16022537),
            (16140627, 16228649),
            (16228749, 16282882),
            (16282982, 16346517),
            (16346617, 16367287),
            (16367387, 16374057),
            (16374157, 16377402),
            (16377502, 16400063),
            (16404448, 17538659),
            (17538759, 17539099),
            (17539199, 17539539),
            (17539639, 17540837),
            (17540937, 18173523),
        ),
        "chr16": (
            (34339329, 34521510),
        ),
        "chr19": (
            (24498980, 24552652),
            (24552752, 24891256),
            (24908689, 27190874),
        ),
        "chr21": (
            (5010000, 5166246),
            (5966593, 6161371),
            (6427258, 6580181),
            (6789085, 6934219),
            (7743700, 7865746),
            (10864560, 10887097),
            (10887197, 10975219),
            (10975319, 11029452),
            (11029552, 11093087),
            (11093187, 11113857),
            (11113957, 11120627),
            (11120727, 11123972),
            (11124072, 11146633),
            (11146733, 12280944),
            (12281044, 12281384),
            (12281484, 12281824),
            (12281924, 12283122),
            (12283222, 12915808),
        ),
        "chr22": (
            (12954788, 12977325),
            (12977425, 12986171),
            (12994027, 13011653),
            (13014130, 13021322),
            (13021422, 13109444),
            (13109544, 13163677),
            (13163777, 13227312),
            (13227412, 13248082),
            (13248182, 13254852),
            (13254952, 13258197),
            (13258297, 13280858),
            (13280958, 13285143),
            (13285243, 14419454),
            (14419554, 14419894),
            (14419994, 14420334),
            (14420434, 14421632),
            (14421732, 15054318),
        ),
        "chrX": (
            (37085008, 37099262),
        ),
        "chrY": (
            (10000, 2781479),
            (56887902, 57217415),
        ),
        "chrUn_KI270752v1": (
            (0, 27745),
        ),
        "chr11_KI270721v1_random": (
            (0, 100316),
        ),
        "chr22_KI270734v1_random": (
            (0, 165050),
        ),
    },
    "mm10": {
        "chrY": (
            (90745844, 91644698),
        ),
    },
}


def _reference_mask_blocks(species: str, fasta: Path, record: tuple[str, int, int, int, int], intervals):
    """按 FAI 定位官方 analysis-set 区间，保留原 FASTA 换行及所有坐标（BED0）。"""
    name, length, offset, line_bases, line_width = record
    if length != dict(managed_alignment_contigs(species)).get(name):
        raise ValueError(f"analysis masking requires standard {species} contig {name}")
    if offset < 0 or line_bases <= 0 or line_width < line_bases:
        raise ValueError("invalid FASTA/FAI geometry for analysis mask")
    with fasta.open("rb") as handle:
        for lo, hi in intervals:
            start = offset + lo // line_bases * line_width + lo % line_bases
            stop = offset + (hi - 1) // line_bases * line_width + (hi - 1) % line_bases + 1
            handle.seek(start)
            data = handle.read(stop - start)
            sequence = data.translate(None, b"\r\n")
            if len(data) != stop - start or len(sequence) != hi - lo or not set(sequence.upper()) <= set(b"ACGTRYSWKMBDHVN"):
                raise ValueError(f"invalid FASTA/FAI at {name}:{lo}-{hi}")
            yield start, data, hi - lo


def mask_analysis_reference(species: str, fasta: Path) -> None:
    """仅修改 Doctor 未提交的候选 FASTA；等长屏蔽官方 analysis-set 区间。

    PAR 的可比对同源信号保留在 X，不能再解释为 X 独有覆盖。mm10 PAR
    含 gap/序列差异，不宣称 X/Y 整段逐碱基相同或全部信号可无损回收。
    https://www.ncbi.nlm.nih.gov/grc/human （GRCh38 PAR1/PAR2）
    https://www.ncbi.nlm.nih.gov/grc/mouse （GRCm38 PAR）
    """
    if species not in _ANALYSIS_MASK_INTERVALS:
        raise ValueError(f"unsupported analysis-set assembly: {species}")
    records = [line.split("\t") for line in Path(str(fasta) + ".fai").read_text().splitlines() if line.strip()]
    by_name = {r[0]: (r[0], *(int(v) for v in r[1:5])) for r in records}
    if len(by_name) != len(records):
        raise ValueError("duplicate contig in candidate FASTA index")
    # 所有目标区间写前验证；只改碱基、不改布局，原 FAI 仍精确有效。
    blocks = []
    for name, intervals in _ANALYSIS_MASK_INTERVALS[species].items():
        if name not in by_name:
            raise ValueError(f"missing analysis-set contig: {name}")
        blocks.extend(_reference_mask_blocks(species, fasta, by_name[name], intervals))
    translation = bytes(c if c in (10, 13) else ord("N") for c in range(256))
    with fasta.open("r+b") as handle:
        for start, data, _length in blocks:
            masked = data.translate(translation)
            if data != masked:
                handle.seek(start)
                handle.write(masked)


def fai_contract_check(mode: str, label: str, fai: Path, fasta: Path) -> bool:
    """校验 source/assembled FASTA 是否满足标准 contig 长度契约。

    只读取很小的 .fai 并通过 FAI 几何 + FASTA 大小的算术检查排除截断文件，
    assembled hg38/mm10 额外只读对应 analysis-set 区间以核对 hard-mask，不扫描全基因组。
    mode=source 允许 UCSC FASTA 中的额外 contig（primary 与
    spike 必须全部存在且长度精确）；mode=assembled 只允许标准 contig 的精确顺序。
    """

    records: list[tuple[str, int, int, int, int]] = []
    seen: set[str] = set()
    for line in fai.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 5 or not fields[0]:
            return False
        name = fields[0]
        if name in seen:
            return False
        seen.add(name)
        try:
            length, offset, line_bases, line_width = (int(v) for v in fields[1:5])
        except ValueError:
            return False
        if length <= 0 or offset < 0 or line_bases <= 0 or line_width < line_bases:
            return False
        records.append((name, length, offset, line_bases, line_width))
    if not records:
        return False


    fasta_size = fasta.stat().st_size
    for _name, length, offset, line_bases, line_width in records:
        full_rows = (length - 1) // line_bases
        final_row_start = offset + full_rows * line_width
        final_row_bases = (length - 1) % line_bases + 1
        if final_row_start + final_row_bases > fasta_size:
            return False

    length_records = [(name, length) for name, length, *_ in records]
    if mode == "source":
        if label in MANAGED_SPECIES:
            observed = dict(length_records)
            return all(
                observed.get(name) == length
                for name, length in _PRIMARY_CONTIGS[label]
            )
        if label in SPIKE_CONTIGS:
            return (
                len(length_records) == 1
                and length_records[0][1] == SPIKE_CONTIGS[label]
            )
        return False
    if mode == "assembled":
        if label not in MANAGED_CONTIGS or length_records != managed_alignment_contigs(label):
            return False
        if label in _ANALYSIS_MASK_INTERVALS:
            by_name = {r[0]: r for r in records}
            try:
                return all(data.translate(None, b"\r\n") == b"N" * length
                           for name, intervals in _ANALYSIS_MASK_INTERVALS[label].items()
                           for _start, data, length in _reference_mask_blocks(label, fasta, by_name[name], intervals))
            except ValueError:
                return False
        return True
    return False


# Doctor 维护：测试路线和 fixture

DEFAULT_ROUTE = "hg38-cabernet-bismark"


TEST_ROUTES: dict[str, dict[str, str]] = {
    "hg38-taps-rastair": {"species": "hg38", "protocol": "taps", "methylation_backend": "rastair", "reference_key": "hg38"},
    "hg38-cabernet-bismark": {
        "species": "hg38",
        "protocol": "cabernet",
        "methylation_backend": "bismark",
        "reference_key": "hg38",
    },
    "mm10-srd-biscuit": {
        "species": "mm10",
        "protocol": "srd",
        "methylation_backend": "biscuit",
        "reference_key": "mm10",
    },
}


def resolve_route(route: str) -> str:
    """返回具体 route id；空值与 auto 应用隐式默认。"""
    name = str(route or "").strip()
    if name == "auto" or not name:
        return DEFAULT_ROUTE
    if name not in TEST_ROUTES:
        raise ValueError(f"unsupported test route: {name}")
    return name


def apply_common_route_config(cfg: dict[str, Any], route_id: str) -> None:
    """把 route 定义的分析子集应用到目标 config（各消费者共享的唯一实现）。

    写出前合并内置默认，保证合成的项目/运行时 config 总是完整可校验。
    """
    meta = TEST_ROUTES[resolve_route(route_id)]
    cfg["species"] = meta["species"]
    references = cfg.setdefault("references", {})
    if not str(references.get(meta["reference_key"], "")).strip():
        raise SystemExit(
            f"ERROR: references.{meta['reference_key']} is required for {meta['protocol']} test mode"
        )
    analysis = cfg.setdefault("analysis", {})
    analysis["protocol"] = meta["protocol"]
    apply_config_defaults(cfg)
    module = cfg[meta["protocol"]]
    if meta["protocol"] not in _PROTOCOL_FIXED_CONFIG:
        module["methylation_backend"] = meta["methylation_backend"]
    if meta["protocol"] in ("cabernet", "srd"):
        module.setdefault("biscuit", {})["generate_snp"] = False
    if meta["methylation_backend"] == "bismark":
        module.setdefault("bismark", {})["library_type"] = "non_directional"


def synthesize_route_config(source: Path, destination: Path, route: str) -> None:
    """把项目 config 与 route 覆盖合成为 test 模式 effective config。"""
    import yaml

    # Doctor 也必须先拒绝重复键/非法源配置，不能先覆盖路线再把错误输入洗成合法配置。
    cfg = load_project_config(source)
    apply_common_route_config(cfg, resolve_route(route))


    runtime = cfg.setdefault("runtime", {})
    runtime["local_executor_cores"] = max(8, int(runtime.get("local_executor_cores", 8)))
    atomic_write_bytes(destination, yaml.safe_dump(
        finalize_project_config(cfg), sort_keys=False, allow_unicode=True,
    ).encode("utf-8"))


BENCHMARK_SCHEMA_VERSION = 2


BENCHMARK_ID = "doctor_hpc_raw_v2"


def _benchmark_fastq_pairs(r1: Path, r2: Path) -> Iterable[tuple[bytes, bytes]]:
    with gzip.open(r1, "rb") as left, gzip.open(r2, "rb") as right:
        while True:
            records = [[handle.readline() for _ in range(4)] for handle in (left, right)]
            if not records[0][0] and not records[1][0]:
                return
            names = []
            for mate, lines in enumerate(records, 1):
                if (not all(lines) or not lines[0].startswith(b"@") or not lines[2].startswith(b"+")
                        or len(lines[1].rstrip()) != len(lines[3].rstrip())):
                    raise ValueError(f"Malformed Doctor FASTQ: {r1 if mate == 1 else r2}")
                name = lines[0].split()[0]
                names.append(re.sub(rb"/[12]$", b"", name))
            if names[0] != names[1]:
                raise ValueError(f"Doctor FASTQ mate names disagree: {names}")
            yield b"".join(records[0]), b"".join(records[1])


def validate_benchmark(root: Path, manifest_path: Path) -> None:
    """核对真实抽样的文件身份、逐记录配对、来源序号与当前 fixture schema；运行时不访问 HPC。"""
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != BENCHMARK_SCHEMA_VERSION or manifest.get("benchmark_id") != BENCHMARK_ID:
        raise ValueError("Unsupported Doctor benchmark schema/id")

    def verify_entry(entry: dict[str, Any]) -> Path:
        name = str(entry["path"])
        path = root / name
        if Path(name).name != name or not path.is_file() or path.is_symlink():
            raise ValueError(f"Invalid Doctor benchmark file: {name}")
        if sha256_file(path) != entry["sha256"]:
            raise ValueError(f"Doctor benchmark SHA256 mismatch: {name}")
        return path

    verify_entry(manifest["barcode_map"])
    total = 0
    for label, dataset in manifest["datasets"].items():
        files = {entry["path"]: verify_entry(entry) for entry in dataset["fastq_files"]}
        r1 = sorted(name for name in files if "_R1_" in name)
        if not r1 or set(files) != set(r1 + [name.replace("_R1_", "_R2_") for name in r1]):
            raise ValueError("Doctor benchmark FASTQ mate inventory is invalid")
        count = 0
        names = set()
        for name in r1:
            for left, _right in _benchmark_fastq_pairs(files[name], files[name.replace("_R1_", "_R2_")]):
                qname = re.sub(rb"/[12]$", b"", left.splitlines()[0].split()[0])
                if qname in names:
                    raise ValueError(f"Repeated source read in Doctor fixture: {label}")
                names.add(qname)
                count += 1
        ordinals = dataset["source_provenance"]["selected_ordinals_1based"]
        if (count != dataset["read_pairs"] or len(ordinals) != count or ordinals != sorted(set(ordinals))
                or not ordinals or ordinals[0] < 1 or ordinals[-1] > dataset["source_provenance"]["scanned_pairs"]):
            raise ValueError(f"Doctor benchmark count/source ordinal mismatch: {label}")
        if dataset["source_provenance"]["raw_records_modified"] is not False:
            raise ValueError("Doctor real FASTQ must preserve the source records")
        total += count
    print(f"Fixed benchmark verified: {BENCHMARK_ID} | real pairs={total} | datasets={len(manifest['datasets'])}")


def _write_benchmark_fastqs(project: Path, stem: str, pairs: Sequence[tuple[str, str, str]]) -> None:
    for mate in (1, 2):
        with (project / "01_raw" / f"{stem}_R{mate}_001.fastq.gz").open("wb") as raw:
            with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as handle:
                for name, left, right in pairs:
                    seq = left if mate == 1 else right
                    handle.write(f"@{name}/{mate}\n{seq}\n+\n{'I' * len(seq)}\n".encode())


def stage_benchmark_supplement(project: Path, reference: Path, label: str, route: str) -> dict[str, object]:
    """为真实抽样补充明确命名的合成对照、重复/CpH 和 SRD 双管 RNA 布局，原始记录不改写。"""
    backend = TEST_ROUTES[route]["methylation_backend"]
    samtools = PIPELINE_ROOT / "conda/current" / backend / "bin/samtools"
    fasta = subprocess.check_output([str(samtools), "faidx", str(reference), "chrM", "lambda", "pUC19"], text=True)
    sequences: dict[str, str] = {}
    for line in fasta.splitlines():
        if line.startswith(">"):
            chrom = line[1:].split()[0]
            sequences[chrom] = ""
        else:
            sequences[chrom] += line.upper()
    with (project / "00_config/Barcode_Map.csv").open() as handle:
        barcodes = list(csv.DictReader(handle))[:2]
    complement = str.maketrans("ACGT", "TGCA")
    pairs = []
    for chrom, seq in sequences.items():
        for index in range(12):
            start = 300 + 53 * index
            fragment = seq[start:start + 180]
            if len(fragment) != 180 or set(fragment) - set("ACGT"):
                raise ValueError(f"Doctor control reference is not callable: {chrom}:{start}")
            molecule = "".join("T" if base == "C" and not (chrom == "pUC19" and seq[start + pos:start + pos + 2] == "CG")
                               else base for pos, base in enumerate(fragment))
            barcode = barcodes[index % 2]["DNA_Barcode"]
            pairs.append((f"synthetic_control_{chrom}_{index}", barcode + "AGATGTGTATAAGAGACAG" + "G" * 9 + molecule[:100],
                          molecule[-100:].translate(complement)[::-1]))
    fragment = sequences["chrM"][2000:2180]
    if len(fragment) != 180 or set(fragment) - set("ACGT"):
        raise ValueError("Doctor unconverted control reference is not callable")
    for barcode in barcodes:
        for duplicate in range(2):
            pairs.append((f"synthetic_high_cph_{barcode['PlateID']}_{duplicate}",
                          barcode["DNA_Barcode"] + "AGATGTGTATAAGAGACAG" + "G" * 9 + fragment[:100],
                          fragment[-100:].translate(complement)[::-1]))
    dna_pairs = len(pairs)
    rna = []
    if TEST_ROUTES[route]["protocol"] == "srd":
        for barcode in barcodes:
            for index in range(12):
                payload = sequences["chrM"][300 + index:400 + index]
                umi = "ACGT" + "ACGT"[index % 4] * 4
                rna.append((f"synthetic_srd_{barcode['PlateID']}_{index}",
                            "GTCTAACGCGTTAC" + barcode["RNA_Barcode"] + umi + "A" * 5 + payload,
                            payload.translate(complement)[::-1]))
        pairs.extend((name + "_dna_tube", left, right) for name, left, right in rna)
        _write_benchmark_fastqs(project, label + "_RNA", rna)
    _write_benchmark_fastqs(project, label + "_L003", pairs)
    return {"kind": "synthetic_supplement", "dna_pairs": dna_pairs, "rna_pairs_per_tube": len(rna),
            "controls": ["chrM", "lambda", "pUC19"], "high_cph_duplicate_pairs": 4,
            "native_srd_library": False, "reference": str(reference)}


def stage_taps_truth(project: Path, reference: Path) -> None:
    """从固定参考窗口生成独立 TAPS 双链/重叠/重复/9 bp gap 真值 FASTQ，不复用转化文库 fixture。"""
    samtools = PIPELINE_ROOT / "conda/current/rastair/bin/samtools"
    text = subprocess.check_output([str(samtools), "faidx", str(reference), "chr1:1000001-1006000", "lambda", "pUC19"], text=True)
    sequences: dict[str, str] = {}
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].split(":")[0]
            sequences[name] = ""
        else:
            sequences[name] += line.upper()
    barcode = "ACGTACGT"
    (project / "00_config/Barcode_Map.csv").write_text("DNA_Barcode,PlateID,Cell_Order\nACGTACGT,TAPS,1\n")
    complement = str.maketrans("ACGTN", "TGCAN")
    pairs = []
    for chrom, seq in sequences.items():
        for index in range(40):
            start = 300 + 17 * index
            fragment = seq[start:start + 180]
            if len(fragment) != 180 or "N" in fragment:
                raise ValueError(f"TAPS truth reference window is not callable: {chrom}:{start}")
            converted = list(fragment)
            modified = chrom == "pUC19" or (chrom == "chr1" and index % 4 < 2)
            ot = index % 2 == 0
            if modified:
                for pos in range(len(fragment)):
                    context = start + pos if ot else start + pos - 1
                    if seq[context:context + 2] == "CG":
                        converted[pos] = "T" if ot else "A"
            molecule = "".join(converted)
            pair = (barcode + "AGATGTGTATAAGAGACAG" + "CCCCCCCCC" + molecule[:100],
                    molecule[-100:].translate(complement)[::-1])
            pairs.append((f"taps_{chrom}_{index}", *pair))
            if index % 10 == 0:
                pairs.append((f"taps_{chrom}_{index}_duplicate", *pair))
    pairs.append(("taps_unmapped", barcode + "AGATGTGTATAAGAGACAG" + "CCCCCCCCC" + "N" * 100, "N" * 100))
    for mate in (1, 2):
        path = project / "01_raw" / f"Taps_R{mate}.fastq.gz"
        with path.open("wb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as gz:
            for name, r1, r2 in pairs:
                seq = r1 if mate == 1 else r2
                gz.write(f"@{name}/{mate}\n{seq}\n+\n{'I' * len(seq)}\n".encode())
    atomic_write_json(project / "04_logs/taps_truth.json", {"fixture": "taps_native_v1", "pairs": len(pairs),
                      "sample": "Taps_TAPS", "lambda": "unmodified", "pUC19": "CpG_modified", "host": "mixed",
                      "fragment_bp": 180, "read_bp": 100, "gap_bp": 9, "duplicate_input_pairs": 12,
                      "reference": str(reference)})


def verify_taps_truth(project: Path) -> None:
    """核对独立 TAPS FASTQ 真值、BAM 保留、对照分母和 sealed 交付的边界。"""
    truth = json.loads((project / "04_logs/taps_truth.json").read_text())
    with (project / "03_results/QC_Results/sample_manifest.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError("TAPS truth must publish exactly one cell")
    row = rows[0]
    if row["protocol"] != "taps" or row["methylation_backend"] != "rastair" or int(row["dna_reads"]) != truth["pairs"]:
        raise ValueError("TAPS truth protocol/demux mismatch")
    if int(row["bam_primary_records"]) != 2 * int(row["trimmed_pairs"]) or min(
        int(row["bam_duplicate_records"]), int(row["duplicate_pairs"])
    ) <= 0:
        raise ValueError("TAPS truth BAM retention/duplicate mismatch")
    if not (float(row["lambda_false_positive_pct"]) == 0.0 and float(row["pUC19_conversion_pct"]) == 100.0):
        raise ValueError("TAPS control truth mismatch")
    if min(int(row["lambda_cpg_observations"]), int(row["pUC19_cpg_observations"])) <= 0:
        raise ValueError("TAPS controls have no observations")
    if row["high_cph_role"] != "not_applicable" or float(row["non_cpg_methylation_pct"]) != 0.0 or int(row["non_cpg_observations"]) <= 0:
        raise ValueError("TAPS truth requires measured zero host CHG/CHH conversion")
    if not row["snp_path"] or not (project / row["snp_path"]).is_file() or not Path(str(project / row["snp_path"]) + ".tbi").is_file():
        raise ValueError("TAPS truth SNP VCF/index is missing")
    for key in ("retained_bam_path", "retained_bam_index"):
        if not row[key] or not (project / row[key]).is_file():
            raise ValueError("TAPS retained BAM/index is missing")
    marker = json.loads((project / "03_results/run_manifest.json").read_text())
    if marker["status"] != "complete" or ".bam" in json.dumps(marker["outputs"]):
        raise ValueError("TAPS sealed inventory must exclude BAM")
    if (project / "02_work/demux").exists():
        raise ValueError("TAPS published demux scratch was not reclaimed")
    print(f"TAPS synthetic truth PASS: {truth['pairs']} pairs; unmethylated lambda / modified pUC19; marked BAM retained")


def stage_route_config(project: Path, reference: Path, route: str) -> None:
    """由 Doctor 在独立项目内统一配置环境已准备好的路线、真实抽样及显式合成补充。"""
    import yaml

    route_id = resolve_route(route)
    project, reference = project.resolve(), reference.resolve()
    for folder in ("00_config", "01_raw", "02_work", "03_results", "04_logs", "05_tmp"):
        (project / folder).mkdir(parents=True, exist_ok=True)
    cfg = yaml.safe_load(CONFIG_TEMPLATE_YAML)
    cfg["references"][TEST_ROUTES[route_id]["reference_key"]] = str(reference)
    apply_common_route_config(cfg, route_id)
    fixture_root = PIPELINE_ROOT / "resources/doctor_benchmark"
    fixture = json.loads((fixture_root / "benchmark.json").read_text())
    cfg["common"]["demux"]["min_matched_read_pairs"] = 10
    manifest = "sample_id\tdna_raw_sample\trna_sample\tspecies\tprotocol\tnotes\n"
    if route_id == "hg38-taps-rastair":
        stage_taps_truth(project, reference)
        manifest += "Taps\tTaps\tTaps_RNA\thg38\ttaps\tindependent_synthetic_truth_with_external_RNA_link\n"
        cfg["retention"]["keep_final_bam"] = True
        provenance = {"kind": "independent_synthetic_chemistry_truth", "real_taps_pairs": 0}
    else:
        label = fixture["routes"][route_id]["dataset"]
        dataset = fixture["datasets"][label]
        for entry in dataset["fastq_files"]:
            shutil.copyfile(fixture_root / entry["path"], project / "01_raw" / entry["path"])
        shutil.copyfile(fixture_root / fixture["barcode_map"]["path"], project / "00_config/Barcode_Map.csv")
        supplement = stage_benchmark_supplement(project, reference, label, route_id)
        rna = label + "_RNA"
        manifest += f"{label}\t{label}\t{rna}\t{cfg['species']}\t{cfg['analysis']['protocol']}\treal_Cabernet_DNA_plus_explicit_synthetic_supplement\n"
        provenance = {"dataset": label, "real_raw_pairs": dataset["read_pairs"],
                      "source_provenance": dataset["source_provenance"], "supplement": supplement}
    cfg["slurm"]["jobs"] = 8
    cfg["runtime"]["local_executor_cores"] = 2
    (project / "00_config/config.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True))
    (project / "00_config/sample_manifest.tsv").write_text(manifest)
    atomic_write_json(project / "04_logs/doctor_fixture.json", {"benchmark_id": BENCHMARK_ID, "route": route_id,
                      "benchmark_sha256": sha256_file(fixture_root / "benchmark.json"), **provenance})


def _benchmark_cpg_metrics(path: Path) -> dict[str, dict[str, int]]:
    from compression import zstd

    metrics: dict[str, dict[str, int]] = {}
    seen = set()
    with zstd.open(path, "rt", encoding="utf-8") as handle:
        if handle.readline() != "chrom\tpos\tmethyl\tunmethyl\n":
            raise ValueError(f"Doctor canonical CpG header mismatch: {path}")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 4:
                raise ValueError(f"Doctor canonical CpG must have four columns: {path}")
            chrom, pos, methyl, unmethyl = fields
            pos, methyl, unmethyl = int(pos), int(methyl), int(unmethyl)
            if min(pos, methyl, unmethyl) < 0 or methyl + unmethyl <= 0 or (chrom, pos) in seen:
                raise ValueError(f"Doctor canonical CpG invalid/duplicate observation: {path}")
            seen.add((chrom, pos))
            counts = metrics.setdefault(chrom, {"rows": 0, "methyl": 0, "unmethyl": 0})
            counts["rows"] += 1
            counts["methyl"] += methyl
            counts["unmethyl"] += unmethyl
    return metrics


def verify_route_outputs(project: Path, benchmark_path: Path, route: str) -> None:
    """核对多 cell 的拆分/过滤守恒、实际宿主和对照 CpG、SRD 保留及密封交付，保存可追溯回归报告。"""
    if route not in TEST_ROUTES:
        raise ValueError(f"unsupported route: {route}")
    benchmark = json.loads(benchmark_path.read_text())
    with (project / "03_results/QC_Results/sample_manifest.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    run_manifest = json.loads((project / "03_results/run_manifest.json").read_text())
    if run_manifest.get("schema_version") != 8 or run_manifest.get("status") != "complete":
        raise ValueError(f"{route}: sealed run manifest is not complete schema 8")
    if TEST_ROUTES[route]["protocol"] in {"cabernet", "taps"}:
        expected_link = "Taps_RNA" if route == "hg38-taps-rastair" else "hg38_cabernet_RNA"
        if any(row.get("rna_sample") != expected_link or row.get("rna_enrichment_route_id") for row in rows):
            raise ValueError(f"{route}: external RNA association was lost or became a local RNA route")
    if not (project / "03_results/QC_Results/multiqc_report.html").is_file():
        raise ValueError("Doctor MultiQC report missing")
    expected = benchmark["routes"][route]
    cpg_metrics = {}
    if route == "hg38-taps-rastair":
        verify_taps_truth(project)
        for row in rows:
            cpg_metrics[row["sample_id"]] = _benchmark_cpg_metrics(project / row["cpg_path"])
    else:
        if len(rows) != len(expected["cells"]) or {row["sample_id"] for row in rows} != set(expected["cells"]):
            raise ValueError(f"{route}: unexpected final manifest cell set")
        active = []
        for row in rows:
            sample = row["sample_id"]
            if (row["protocol"] != TEST_ROUTES[route]["protocol"] or row["species"] != TEST_ROUTES[route]["species"]
                    or row["methylation_backend"] != TEST_ROUTES[route]["methylation_backend"]):
                raise ValueError(f"{sample}: route identity mismatch")
            wanted = {**expected["cells"][sample], **{key: expected[key] for key in ("dna_input_read_pairs", "dna_matched_read_pairs")}}
            for field, value in wanted.items():
                if row.get(field) != str(value):
                    raise ValueError(f"{sample}: demux {field}={row.get(field)!r}, expected {value!r}")
            if row["dna_status"] != "Pass":
                if row["analysis_status"] != "NOT_ANALYZED" or row["cpg_path"]:
                    raise ValueError(f"{sample}: low/zero barcode was analyzed")
                continue
            active.append(row)
            for field, minimum in (("backend_accepted_pairs", expected["minimum_accepted_pairs_per_active_cell"]),
                                   ("final_retained_pairs", expected["minimum_final_pairs_per_active_cell"])):
                if int(row[field]) < minimum:
                    raise ValueError(f"{sample}: {field} below {minimum}")
            trimmed, accepted, rejected, duplicate, postdedup, removed, final = (int(row[field]) for field in
                ("trimmed_pairs", "backend_accepted_pairs", "backend_rejected_pairs", "duplicate_pairs", "postdedup_pairs",
                 "high_cph_removed_pairs", "final_retained_pairs"))
            if trimmed != accepted + rejected or accepted != duplicate + postdedup or postdedup != removed + final:
                raise ValueError(f"{sample}: pair funnel does not close")
            counts = _benchmark_cpg_metrics(project / row["cpg_path"])
            cpg_metrics[sample] = counts
            if sum(value["rows"] for chrom, value in counts.items() if chrom not in {"chrM", "lambda", "pUC19"}) < expected["minimum_host_cpg_rows_per_active_cell"]:
                raise ValueError(f"{sample}: real host genome CpG coverage missing")
            if any(counts.get(chrom, {}).get("rows", 0) == 0 for chrom in ("chrM", "lambda", "pUC19")):
                raise ValueError(f"{sample}: synthetic control CpG coverage missing")
            if route == "mm10-srd-biscuit":
                for field in ("dna_tube_rna_reads", "rna_enrichment_rna_reads"):
                    if int(row[field]) != 12:
                        raise ValueError(f"{sample}: SRD synthetic RNA count mismatch")
                for field in ("dna_tube_rna_fastq_r1", "dna_tube_rna_fastq_r2", "rna_enrichment_rna_fastq_r1", "rna_enrichment_rna_fastq_r2", "rna_bam_path"):
                    if not row[field] or not (project / row[field]).is_file():
                        raise ValueError(f"{sample}: retained SRD {field} missing")
                for prefix in ("dna_tube_rna_fastq", "rna_enrichment_rna_fastq"):
                    if sum(1 for _ in _benchmark_fastq_pairs(project / row[prefix + "_r1"], project / row[prefix + "_r2"])) != 12:
                        raise ValueError(f"{sample}: collected RNA FASTQ pair count mismatch")
                if not (project / (row["rna_bam_path"] + ".bai")).is_file():
                    raise ValueError(f"{sample}: SRD retained BAM index missing")
        for field, minimum in (("duplicate_pairs", expected["minimum_duplicate_pairs_total"]),
                               ("high_cph_removed_pairs", expected["minimum_high_cph_pairs_total"])):
            if sum(int(row[field]) for row in active) < minimum:
                raise ValueError(f"{route}: {field} regression branch was not exercised")
    declared = {row["cpg_path"].removeprefix("03_results/") for row in rows if row["cpg_path"]}
    public = {str(path.relative_to(project / "03_results")) for path in (project / "03_results/CpG").glob("*")}
    if public != declared or ".bam" in json.dumps(run_manifest["outputs"]):
        raise ValueError("Doctor public inventory contains unexpected scientific output/BAM")
    if (project / "02_work/demux").exists():
        raise ValueError("Doctor published demux scratch was not reclaimed")
    report = {"status": "PASS", "route": route, "fixture": json.loads((project / "04_logs/doctor_fixture.json").read_text()),
              "run": run_manifest["run"], "snapshot": json.loads(Path(run_manifest["run"]["snapshot_path"]).read_text()),
              "cells": rows, "cpg_counts": cpg_metrics, "inventory": run_manifest["outputs"]}
    atomic_write_json(project / "04_logs/doctor_verification.json", report)
    print(f"{route}: PASS; {len(declared)} active / {len(rows)} total cells; exact demux, native mapping/filter funnel, CpG and retention verified")
