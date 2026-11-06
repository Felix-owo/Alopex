"""原始输入、manifest、barcode及解复用提交边界。"""
from __future__ import annotations

import sys

_LOADED_MODULE_CODE = sys._getframe().f_code

import stat
import csv
import gzip
import hashlib
import io
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
import time
import uuid
import warnings
from collections.abc import Iterable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from functools import lru_cache
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from types import MappingProxyType
from typing import Any, TextIO, TypedDict, cast

from .config import (
    BarcodeRow,
    CONFIG_TEMPLATE_YAML,
    CellIdentityRow,
    DROPLET_DESIGN_RESOURCE,
    DROPLET_DESIGN_SHA256,
    DemuxIndex,
    DemuxRow,
    DiscoveredSamples,
    MappedDemuxRuns,
    ProjectSampleRow,
    ProtocolPolicy,
    SAMPLE_ID_RE,
    _PROTOCOLS,
    _integer,
    _object,
    _text,
    default_dna_route,
    file_identity,
    get_protocol,
    get_protocol_policy,
    get_raw_route,
    load_project_config,
    logger,
)
from .io import (
    project_temporary_parent,
    atomic_write_bytes,
    _read_json,
    atomic_write_json,
    render_tsv,
    sha256_bytes,
    sha256_file,
    write_tsv_atomic,
)
from .runtime import (
    effective_cpu_count,
)
from .paths import PIPELINE_ROOT


# 输入与解复用：项目和原始输入

FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz")


_FASTQ_PATTERNS = (
    (
        re.compile(r"^(?P<sample>.+)\.R(?P<read>[12])\.raw\.f(?:ast)?q\.gz$", re.IGNORECASE),
        lambda m: "single",
    ),
    (
        re.compile(
            r"^(?P<sample>.+)[_.-]L(?P<lane>[0-9]{1,3})[_.-]R?(?P<read>[12])(?:[_.-](?P<chunk>[0-9]{3}))?\.f(?:ast)?q\.gz$",
            re.IGNORECASE,
        ),
        lambda m: f"L{int(m.group('lane')):03d}_{m.group('chunk') or '001'}",
    ),
    (
        re.compile(r"^(?P<sample>.+)[_.-]R?(?P<read>[12])[_.-](?P<chunk>[0-9]{3})\.f(?:ast)?q\.gz$", re.IGNORECASE),
        lambda m: f"part_{m.group('chunk')}",
    ),
    (
        re.compile(r"^(?P<sample>.+)[_.-]R?(?P<read>[12])\.f(?:ast)?q\.gz$", re.IGNORECASE),
        lambda m: "single",
    ),
)


def _parse_fastq_inventory_name(path: Path) -> tuple[str, str, str] | None:
    for pattern, segment_fn in _FASTQ_PATTERNS:
        match = pattern.match(path.name)
        if match:
            return match.group("sample"), f"R{match.group('read')}", segment_fn(match)
    return None


def _iter_fastq_inventory(
    dir_raw: Path, planned_renames: Mapping[Path, Path] | None = None,
    selected_samples: set[str] | None = None,
) -> Iterable[tuple[Path, str, str, str]]:
    unsupported: list[Path] = []
    raw_names: dict[str, str] = {}
    ignored_count = 0
    ignored_preview: list[str] = []
    for path in sorted(dir_raw.rglob("*")):
        if not path.is_file():
            continue
        if not path.name.lower().endswith(FASTQ_SUFFIXES):
            continue
        if planned_renames:
            path = planned_renames.get(path, path)
        parsed = _parse_fastq_inventory_name(path)
        if selected_samples is not None and (
            parsed is None or parsed[0].removeprefix("样本_") not in selected_samples
        ):
            ignored_count += 1
            if len(ignored_preview) < 5:
                ignored_preview.append(path.name)
            continue
        if parsed is None:
            unsupported.append(path)
            continue
        original_sample, mate, segment = parsed
        if segment.startswith("L000_"):
            raise ValueError(f"FASTQ lane must be within 1-999: {path}")
        sample = original_sample.removeprefix("样本_")
        if SAMPLE_ID_RE.fullmatch(sample) is None:
            raise ValueError(
                "Raw sample names must use letters, digits, '.', '_' and '-' and cannot be '.' or '..': "
                f"{sample!r} from {path} (only a leading '样本_' prefix is removed)"
            )
        previous = raw_names.setdefault(sample, original_sample)
        if previous != original_sample:
            raise ValueError(
                f"FASTQ sample identity collision after removing '样本_': "
                f"{previous!r} and {original_sample!r} both map to {sample!r}; "
                "use distinct sample basenames instead of merging them implicitly"
            )
        yield path, sample, mate, segment

    if ignored_count:
        logger.warning("未登记的 FASTQ 已跳过：%d 个文件；示例：%s", ignored_count, "; ".join(ignored_preview))
    if unsupported:
        preview = "; ".join(str(path) for path in unsupported[:10])
        more = f"; ... (+{len(unsupported) - 10} more)" if len(unsupported) > 10 else ""
        raise ValueError(
            "FASTQ-like file(s) use unsupported names under "
            f"{dir_raw}: {preview}{more}. Supported examples: "
            "sample.R1.raw.fastq.gz, sample_L001_R1_001.fastq.gz, sample_L1_1.fq.gz, "
            "sample_R1_001.fastq.gz, sample_R1.fastq.gz, sample_1.fq.gz, sample.1.fq.gz "
            "(and matching R2/2; .fq.gz is also accepted)."
        )


def discover_samples(
    dir_raw: Path, *, planned_renames: Mapping[Path, Path] | None = None,
    selected_samples: set[str] | None = None,
) -> dict[str, dict[str, list[str]]]:
    """按确定性的 lane/chunk 顺序发现配对 FASTQ。

    支持 R1/R2 或末尾 1/2、1-999 的 lane；仅去掉样本名前缀“样本_”，不改文件名或批次信息。
    lane 按数值配对排序；归一化身份冲突、重复 mate 与缺配对均拒绝；维护预览可校验虚拟改名清单。
    指定样本集合时仅校验该集合的输入，其余 FASTQ 提示后跳过；未指定时完整扫描供 refresh 使用。
    """
    logger.info("Starting sample discovery in: %s", dir_raw)
    if not dir_raw.is_dir():
        raise FileNotFoundError(f"Raw data directory not found: {dir_raw}")
    inventory: dict[str, dict[str, dict[str, Path]]] = {}
    duplicates: list[str] = []
    for path, sample, mate, segment in _iter_fastq_inventory(dir_raw, planned_renames, selected_samples):
        segment_map = inventory.setdefault(sample, {}).setdefault(segment, {})
        if mate in segment_map and segment_map[mate] != path:
            duplicates.append(
                f"duplicate {mate} for raw sample {sample!r}, segment {segment!r}: "
                f"{segment_map[mate]} and {path}"
            )
            continue
        segment_map[mate] = path

    if duplicates:
        raise ValueError("Ambiguous FASTQ inventory: " + "; ".join(duplicates[:10]))

    incomplete: list[str] = []
    samples: dict[str, dict[str, list[str]]] = {}
    for sample, segments in sorted(inventory.items()):
        families = {
            "lane" if key.startswith("L") else "part" if key.startswith("part_") else "single"
            for key in segments
        }
        if len(families) > 1:
            raise ValueError(
                f"Raw sample {sample!r} mixes incompatible FASTQ naming layouts: {sorted(segments)}"
            )
        if "single" in families and len(segments) != 1:
            raise ValueError(f"Raw sample {sample!r} has more than one single-file FASTQ segment")
        r1_files: list[str] = []
        r2_files: list[str] = []
        for segment, mates in sorted(segments.items()):
            missing = [mate for mate in ("R1", "R2") if mate not in mates]
            if missing:
                incomplete.append(
                    f"{sample!r} segment {segment!r} missing {','.join(missing)}; found "
                    + ", ".join(str(mates[m]) for m in sorted(mates))
                )
                continue
            r1_files.append(str(mates["R1"].resolve()))
            r2_files.append(str(mates["R2"].resolve()))
        if r1_files and len(r1_files) == len(r2_files):
            samples[sample] = {"r1": r1_files, "r2": r2_files}

    if incomplete:
        raise FileNotFoundError("Incomplete paired FASTQ input: " + "; ".join(incomplete[:10]))

    preview = ", ".join(
        f"{sample}({len(info['r1'])} pair-file{'s' if len(info['r1']) != 1 else ''})"
        for sample, info in list(sorted(samples.items()))[:5]
    )
    if len(samples) > 5:
        preview += f", ... (+{len(samples) - 5} more)"
    logger.info(
        "Finished sample discovery. Found %d paired raw sample(s)%s",
        len(samples),
        f": {preview}" if preview else "",
    )
    return samples


def _scan_fastq_layout(
    dir_raw: Path,
) -> tuple[dict[str, set[str]], dict[str, dict[Path, dict[str, list[tuple[str, Path]]]]]]:
    families: dict[str, set[str]] = {}
    layouts: dict[str, dict[Path, dict[str, list[tuple[str, Path]]]]] = {}
    for path, sample, mate, segment in _iter_fastq_inventory(dir_raw):
        family = (
            "lane"
            if segment.startswith("L")
            else "part"
            if segment.startswith("part_")
            else "single"
        )
        families.setdefault(sample, set()).add(family)
        layouts.setdefault(sample, {}).setdefault(path.parent, {}).setdefault(
            mate, []
        ).append((family, path))
    return families, layouts


def resolve_single_pair_fastq_collisions(
    dir_raw: Path,
    *,
    apply: bool = True,
) -> list[tuple[Path, Path]]:
    """把散落多个文件夹的同名单式 FASTQ 对改名为 chunk 式 ``{sample}_R{mate}_NNN.fastq.gz``
    （文件夹按字典序），使下一次 ``discover_samples`` 将其合并为同一多 chunk 样本。

    lane/chunk 命名混布、R1/R2 不成对、同 mate 字节级重复副本、目标名已存在均拒绝且不改动文件树；
    ``apply=False`` 只返回计划，正式应用后必须通过全新 ``discover_samples`` 验证，
    失败则按反向顺序回滚全部改名。
    """
    dir_raw = Path(dir_raw)
    sample_families, sample_layouts = _scan_fastq_layout(dir_raw)

    plans: list[tuple[Path, Path]] = []

    def _refuse(reason: str) -> None:
        raise ValueError(
            f"FASTQ collision normalization refused ({reason}); nothing was renamed"
        )

    for sample in sorted(sample_families):
        layout = sample_layouts[sample]

        plain: dict[Path, dict[str, list[Path]]] = {}
        for parent, mates in layout.items():
            for mate, records in mates.items():
                for family, path in records:
                    if family == "single":
                        plain.setdefault(parent, {}).setdefault(mate, []).append(path)
        if not plain:
            continue
        duplicated_mate = any(
            len({parent for parent, mm in plain.items() if mate in mm}) >= 2
            for mate in ("R1", "R2")
        )
        if not duplicated_mate:
            continue

        if sample_families[sample] != {"single"}:
            _refuse(
                f"sample {sample!r} owns files outside the plain single-pair naming "
                "(lane or chunk style); rename the batches manually instead"
            )


        ordered_parents = sorted(plain, key=str)
        for parent in ordered_parents:
            broken = [
                f"{mate}x{len(paths)}"
                for mate, paths in sorted(plain[parent].items())
                if len(paths) != 1
            ] + [
                f"{mate}x0" for mate in ("R1", "R2") if mate not in plain[parent]
            ]
            if broken:
                _refuse(
                    f"folder {parent} holds raw sample {sample!r} with "
                    f"{', '.join(broken)} instead of exactly one R1/R2 pair"
                )

        digests: dict[str, dict[str, str]] = {"R1": {}, "R2": {}}
        for mate in ("R1", "R2"):
            for parent in ordered_parents:
                path = plain[parent][mate][0]
                digests[mate][str(path)] = sha256_file(path, chunk_size=1024 * 1024)
            seen: dict[str, str] = {}
            for owner in sorted(digests[mate], key=str):
                fingerprint = digests[mate][owner]
                if fingerprint in seen:
                    _refuse(
                        f"two {mate} files of sample {sample!r} are byte-identical and look "
                        f"like duplicate copies instead of independent data: {seen[fingerprint]} "
                        f"and {owner}; delete the redundant copy manually"
                    )
                seen[fingerprint] = owner

        for seq, parent in enumerate(ordered_parents, start=1):
            for mate in ("R1", "R2"):
                source = plain[parent][mate][0]
                target = parent / f"{sample}_{mate}_{seq:03d}.fastq.gz"
                if target.exists() or target.is_symlink():
                    _refuse(f"planned target already exists: {target}")
                plans.append((source, target))

    if not plans or not apply:
        return plans

    with _fastq_rename_transaction(dir_raw, plans):
        pass
    return plans


@contextmanager
def _fastq_rename_transaction(
    dir_raw: Path, plans: Sequence[tuple[Path, Path]],
) -> Iterator[None]:
    """改名、重新发现及调用方提交共用回滚边界；保留中断前真实完成的文件操作。"""
    attempted: list[tuple[Path, Path]] = []
    try:
        for source, target in plans:
            if target.exists() or target.is_symlink():
                raise FileExistsError(f"FASTQ rename target already exists: {target}")
            # rename 已成功但尚未返回 Python 时也可能收到 Ctrl-C，须先登记操作。
            attempted.append((source, target))
            source.rename(target)
        if plans:
            try:
                discover_samples(dir_raw)
            except Exception as exc:
                raise RuntimeError(
                    f"Collision renames were rolled back because the merged inventory under "
                    f"{dir_raw} still does not validate: {exc}"
                ) from exc
        yield
    except BaseException:
        for source, target in reversed(attempted):
            if target.exists() or target.is_symlink():
                if source.exists() or source.is_symlink():
                    raise FileExistsError(f"FASTQ rollback source already exists: {source}")
                target.rename(source)
        raise
    if plans:
        logger.info(
            "Normalized %d colliding FASTQ file(s) into chunk-style segments", len(plans)
        )


_BARCODE_RE = re.compile(r"^[ACGT]+$")


_SUPPORTED_BARCODE_LENGTHS = {8, 10}


def _header_map(fieldnames: list[str] | None) -> dict[str, str]:
    headers: dict[str, str] = {}
    for name in fieldnames or []:
        key = (name or "").lstrip("\ufeff").strip().lower()
        if not key or key in headers:
            raise ValueError(f"Barcode_Map.csv has a blank or duplicate column: {name!r}")
        headers[key] = name
    return headers


def _validate_barcode(value: str, *, column: str, row_number: int) -> None:
    if not _BARCODE_RE.fullmatch(value):
        raise ValueError(
            f"Barcode_Map.csv row {row_number} column {column} must contain only "
            f"uppercase protected A/C/G/T bases; found {value!r}"
        )
    if len(value) not in _SUPPORTED_BARCODE_LENGTHS:
        raise ValueError(
            f"Barcode_Map.csv row {row_number} column {column} must be 8 or 10 bp; "
            f"found {len(value)} for {value!r}"
        )


def _audit_one_substitution_neighborhoods(values: list[str], label: str) -> None:
    by_length: dict[int, list[str]] = {}
    for value in values:
        by_length.setdefault(len(value), []).append(value)
    conflicts: list[tuple[str, str, int]] = []
    for barcodes in by_length.values():
        for i, left in enumerate(barcodes):
            for right in barcodes[i + 1 :]:
                distance = sum(a != b for a, b in zip(left, right, strict=True))
                if distance <= 2:
                    conflicts.append((left, right, distance))
                    if len(conflicts) >= 12:
                        break
            if len(conflicts) >= 12:
                break
        if len(conflicts) >= 12:
            break
    if conflicts:
        detail = "; ".join(f"{a}<->{b}:d={d}" for a, b, d in conflicts)
        raise ValueError(
            f"{label} one-substitution correction neighborhoods overlap; "
            f"unique radius-1 correction is not guaranteed: {detail}"
        )


def _barcode_field(
    row: Mapping[str | None, object], headers: Mapping[str, str], key: str
) -> str:
    original = headers.get(key)
    return str(row.get(original, "") or "").strip() if original else ""


def validate_barcode_map_for_protocol(
    path: Path, policy: ProtocolPolicy
) -> list[BarcodeRow]:
    """按当前 protocol 读取并校验带类型的 barcode 行。"""
    path = Path(path)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        headers = _header_map(reader.fieldnames)
        required = ["dna_barcode", "cell_order"] + ([] if policy.name == "droplet" else ["plateid"])
        if policy.requires_rna_barcode:
            required.append("rna_barcode")
        missing = [name for name in required if name not in headers]
        if missing:
            raise ValueError(
                f"Barcode map {path} is incompatible with protocol={policy.name!r}; "
                f"missing column(s): {', '.join(missing)}"
            )

        seen_dna: set[str] = set()
        seen_rna: set[str] = set()
        seen_plate: set[str] = set()
        seen_order: set[str] = set()
        dna_lengths: set[int] = set()
        rna_lengths: set[int] = set()
        usable_rows = 0
        barcode_rows: list[BarcodeRow] = []

        for row_index, raw_row in enumerate(reader, start=2):
            if None in raw_row or any(value is None for value in raw_row.values()):
                raise ValueError(f"Barcode_Map.csv row {row_index} field count differs from its header")
            dna = _barcode_field(raw_row, headers, "dna_barcode")
            rna = _barcode_field(raw_row, headers, "rna_barcode")
            plate = _barcode_field(raw_row, headers, "plateid")
            order = _barcode_field(raw_row, headers, "cell_order")
            if not any((dna, rna, plate, order)):
                continue
            usable_rows += 1
            if not dna or (not plate and policy.name != "droplet") or not order or (policy.requires_rna_barcode and not rna):
                raise ValueError(
                    f"Barcode_Map.csv row {row_index} has empty required field(s) for "
                    f"protocol={policy.name}: DNA_Barcode={dna!r}, RNA_Barcode={rna!r}, "
                    f"PlateID={plate!r}, Cell_Order={order!r}"
                )

            if policy.name == "droplet":
                cell_identity("sample", "droplet", plate, dna)
                if dna not in droplet_design_barcodes():
                    raise ValueError("Droplet called barcode is outside the fixed DD-MET5 design")
                if rna:
                    raise ValueError("Droplet does not support RNA barcodes")
            else:
                _validate_barcode(dna, column="DNA_Barcode", row_number=row_index)
            if plate and SAMPLE_ID_RE.fullmatch(plate) is None:
                raise ValueError(f"Barcode_Map.csv row {row_index} has an unsafe PlateID: {plate!r}")
            dna_lengths.add(len(dna))
            if rna:
                _validate_barcode(rna, column="RNA_Barcode", row_number=row_index)
                rna_lengths.add(len(rna))

            if dna in seen_dna:
                raise ValueError(f"Duplicate DNA_Barcode {dna!r} at row {row_index}")
            if plate and plate in seen_plate:
                raise ValueError(f"Duplicate PlateID {plate!r} at row {row_index}")
            if order in seen_order:
                raise ValueError(f"Duplicate Cell_Order {order!r} at row {row_index}")
            seen_dna.add(dna)
            seen_plate.add(plate)
            seen_order.add(order)

            if policy.requires_rna_barcode:
                if rna in seen_rna:
                    raise ValueError(f"Duplicate RNA_Barcode {rna!r} at row {row_index}")
                seen_rna.add(rna)
            barcode_rows.append(
                {
                    "dna_barcode": dna,
                    "rna_barcode": rna,
                    "plate_id": plate,
                    "cell_order": order,
                }
            )

        if usable_rows == 0:
            raise ValueError(f"Barcode map contains no usable rows: {path}")
        if not dna_lengths.issubset({17} if policy.name == "droplet" else _SUPPORTED_BARCODE_LENGTHS):
            raise ValueError(f"Unsupported DNA barcode lengths: {sorted(dna_lengths)}")
        if policy.requires_rna_barcode and len(rna_lengths) != 1:
            raise ValueError(
                "SRD/dna-rna mode requires one uniform RNA barcode length; "
                f"found {sorted(rna_lengths)}"
            )
        if policy.name != "droplet":
            _audit_one_substitution_neighborhoods(sorted(seen_dna), "DNA_Barcode")
        if policy.requires_rna_barcode:
            _audit_one_substitution_neighborhoods(sorted(seen_rna), "RNA_Barcode")
        return barcode_rows


PROJECT_SAMPLE_MANIFEST_FIELDS = (
    "sample_id",
    "dna_raw_sample",
    "rna_sample",
    "species",
    "protocol",
    "notes",
)


INITIAL_CELL_MANIFEST_FIELDS = (
    "sample_id",
    "project_sample_id",
    "species",
    "protocol",
    "plate_id",
    "cell_order",
    "dna_barcode",
    "rna_barcode",
    "dna_raw_sample",
    "rna_sample",
    "notes",
)


def _require_header(fields: list[str], expected: tuple[str, ...], *, source: Path, label: str) -> None:
    if fields != list(expected):
        raise ValueError(
            f"{label} header/order does not match the current contract: {source}; "
            f"expected={list(expected)!r}, observed={fields!r}"
        )


def read_project_sample_manifest(
    path: Path,
    *,
    species: str,
    protocol: str,
) -> list[ProjectSampleRow]:
    """读取并校验用户维护的项目 sample manifest。"""
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return _parse_project_sample_manifest(handle, path=path, species=species, protocol=protocol)


def _parse_project_sample_manifest(
    handle: TextIO, *, path: Path, species: str, protocol: str
) -> list[ProjectSampleRow]:
    input_columns = get_protocol(protocol).raw_input_columns
    reader = csv.DictReader(handle, delimiter="\t")
    fields = list(reader.fieldnames or [])
    _require_header(fields, PROJECT_SAMPLE_MANIFEST_FIELDS, source=path, label="sample manifest")
    rows: list[ProjectSampleRow] = []
    seen_ids: dict[str, int] = {}
    seen_raw: dict[tuple[str, str], tuple[str, str, int]] = {}
    for raw in reader:
        line_no = reader.line_num
        if None in raw:
            raise ValueError(
                f"sample_manifest.tsv line {line_no} has more fields than its header"
            )
        row = {key: str(raw.get(key, "") or "").strip() for key in PROJECT_SAMPLE_MANIFEST_FIELDS}
        if not any(row.values()):
            continue
        sample_id = row["sample_id"]
        if not sample_id:
            raise ValueError(f"sample_manifest.tsv line {line_no} has empty sample_id")
        if SAMPLE_ID_RE.fullmatch(sample_id) is None:
            raise ValueError(
                f"sample_manifest.tsv line {line_no} sample_id={sample_id!r} contains unsupported characters"
            )
        if not any(row[column] for column in input_columns):
            raise ValueError(
                f"sample_manifest.tsv line {line_no} must define a local FASTQ input: "
                + " or ".join(input_columns)
            )
        if not row["species"]:
            raise ValueError(f"sample_manifest.tsv line {line_no} has empty species")
        if row["species"] != species:
            raise ValueError(
                f"sample_manifest.tsv line {line_no} species={row['species']!r} "
                f"does not match config species={species!r}"
            )
        row_protocol = row["protocol"].lower()
        if not row_protocol:
            raise ValueError(f"sample_manifest.tsv line {line_no} has empty protocol")
        if row_protocol != protocol:
            raise ValueError(
                f"sample_manifest.tsv line {line_no} protocol={row['protocol']!r} does not match {protocol!r}"
            )
        row["protocol"] = row_protocol
        sid = row["sample_id"]
        if sid in seen_ids:
            raise ValueError(
                f"sample_manifest.tsv line {line_no} repeats sample_id={sid!r} from line {seen_ids[sid]}"
            )
        seen_ids[sid] = line_no
        for column in ("dna_raw_sample", "rna_sample"):
            raw_sample = row[column]
            if not raw_sample:
                continue
            if SAMPLE_ID_RE.fullmatch(raw_sample) is None:
                raise ValueError(
                    f"sample_manifest.tsv line {line_no} {column}={raw_sample!r} contains unsupported characters"
                )
            namespace = "Raw sample" if column in input_columns else "External RNA sample"
            key = (namespace, raw_sample)
            if key in seen_raw:
                other_sid, other_col, other_line = seen_raw[key]
                raise ValueError(
                    f"{namespace} {raw_sample!r} is mapped more than once: "
                    f"{other_sid}/{other_col} (line {other_line}) and {sid}/{column} (line {line_no})"
                )
            seen_raw[key] = (sid, column, line_no)
        rows.append(row)
    if not rows:
        raise ValueError(f"sample_manifest.tsv contains no usable rows: {path}")
    return rows


def _read_existing_manifest_loose(path: Path) -> list[ProjectSampleRow]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if not fields:
            return []
        if fields != list(PROJECT_SAMPLE_MANIFEST_FIELDS):
            raise ValueError(
                f"Existing sample manifest {path} header/order does not match the current contract: "
                f"expected={list(PROJECT_SAMPLE_MANIFEST_FIELDS)!r}, observed={fields!r}"
            )
        rows: list[ProjectSampleRow] = []
        for raw in reader:
            if None in raw:
                raise ValueError(f"Existing sample manifest {path} line {reader.line_num} has too many fields")
            row = {field: str(raw.get(field, "") or "").strip() for field in PROJECT_SAMPLE_MANIFEST_FIELDS}
            if any(row.values()):
                rows.append(row)
        return rows


def build_refreshed_rows(
    discovered_samples: Iterable[str],
    existing_rows: Sequence[Mapping[str, str]],
    *,
    species: str,
    protocol: str,
) -> list[dict[str, str]]:
    """由 FASTQ inventory 构建确定性的项目 manifest 行，并提示移除失去本地输入的样本及关联。"""
    protocol = str(protocol).strip().lower()
    if protocol not in _PROTOCOLS:
        raise ValueError(f"Unsupported protocol for manifest refresh: {protocol!r}")
    input_columns = get_protocol(protocol).raw_input_columns
    discovered = sorted({str(sample).strip() for sample in discovered_samples if str(sample).strip()})
    discovered_set = set(discovered)
    rows_by_id: dict[str, dict[str, str]] = {}
    raw_owners: dict[str, tuple[str, str]] = {}
    external_rna: set[str] = set()
    for old in existing_rows:
        current_dna = str(
            old.get("dna_raw_sample", "") or ""
        ).strip()
        current_rna = str(old.get("rna_sample", "") or "").strip()
        if current_rna and "rna_sample" not in input_columns:
            external_rna.add(current_rna)
        kept_dna = current_dna if current_dna in discovered_set else ""
        kept_rna = current_rna if "rna_sample" not in input_columns or current_rna in discovered_set else ""
        if not kept_dna and not ("rna_sample" in input_columns and kept_rna):
            association = (
                f"；外部 RNA 关联 {current_rna!r} 随该行移除"
                if current_rna and "rna_sample" not in input_columns else ""
            )
            warnings.warn(
                f"样本 {old.get('sample_id', '')!r} 已无本地 FASTQ，刷新结果将移除该行{association}。",
                UserWarning,
                stacklevel=2,
            )
            continue
        sample_id = str(old.get("sample_id", "") or "").strip()
        if not sample_id:
            raise ValueError(
                "An existing manifest row that still matches current FASTQs has an empty sample_id; "
                "assign it once before refresh so biological identity is not guessed."
            )
        if sample_id in rows_by_id:
            raise ValueError(f"Existing sample manifest repeats sample_id={sample_id!r}")
        row = {
            "sample_id": sample_id,
            "dna_raw_sample": kept_dna,
            "rna_sample": kept_rna,
            "species": species,
            "protocol": protocol,
            "notes": str(old.get("notes", "") or "").strip(),
        }
        rows_by_id[sample_id] = row
        for column in input_columns:
            raw_sample = row[column]
            if not raw_sample:
                continue
            previous = raw_owners.get(raw_sample)
            if previous is not None:
                raise ValueError(
                    f"Raw sample {raw_sample!r} is mapped more than once: "
                    f"{previous[0]}/{previous[1]} and {sample_id}/{column}"
                )
            raw_owners[raw_sample] = (sample_id, column)

    unmatched = [sample for sample in discovered if sample not in raw_owners]
    if "rna_sample" not in input_columns:
        misplaced = sorted(set(unmatched) & external_rna)
        if misplaced:
            raise ValueError(
                "FASTQ sample(s) declared only as external RNA links are present in DNA 01_raw: "
                + ", ".join(misplaced)
                + ". Keep RNA FASTQs in the RNA project; DNA inputs require an explicit dna_raw_sample mapping."
            )
    if protocol == "srd" and unmatched:
        warnings.warn(
            "New SRD raw sample(s) default to independent DNA routes: "
            + ", ".join(unmatched)
            + ". Move a raw name to rna_sample explicitly if it is an RNA-enrichment tube.",
            UserWarning,
            stacklevel=2,
        )


    for raw_sample in unmatched:
        if raw_sample in rows_by_id:
            raise ValueError(
                f"New raw sample {raw_sample!r} collides with an existing sample_id; "
                "edit sample_manifest.tsv explicitly instead of inferring a pairing"
            )
        rows_by_id[raw_sample] = {
            "sample_id": raw_sample,
            "dna_raw_sample": raw_sample,
            "rna_sample": "",
            "species": species,
            "protocol": protocol,
            "notes": "",
        }
        raw_owners[raw_sample] = (raw_sample, "dna_raw_sample")

    if set(raw_owners) != discovered_set:
        raise AssertionError("Manifest refresh did not assign every raw sample exactly once")
    return sorted(rows_by_id.values(), key=lambda row: row["sample_id"])


def _canonical_project_manifest_rows(
    rows: Sequence[Mapping[str, object]],
) -> list[dict[str, str]]:


    normalized: list[dict[str, str]] = []
    for raw in rows:
        row = {
            field: str(raw.get(field, "") or "").strip()
            for field in PROJECT_SAMPLE_MANIFEST_FIELDS
        }
        row["protocol"] = row["protocol"].lower()
        normalized.append(row)
    return sorted(normalized, key=lambda row: row["sample_id"])


def _resolve_config(project_arg: str | None, config_arg: str | None) -> tuple[Path, Path]:
    cwd = Path.cwd().resolve()
    if config_arg:
        config_path = Path(config_arg).expanduser()
        config_path = (cwd / config_path).resolve() if not config_path.is_absolute() else config_path.resolve()
        if config_path.parent.name != "00_config" or config_path.name != "config.yaml":
            raise ValueError("Alopex only accepts <project>/00_config/config.yaml")
        project = config_path.parent.parent.resolve()
        if project_arg and Path(project_arg).expanduser().resolve() != project:
            raise ValueError("--project and --configfile resolve to different project roots")
    else:
        project = Path(project_arg).expanduser().resolve() if project_arg else cwd
        config_path = project / "00_config" / "config.yaml"

    if not project.is_dir():
        raise FileNotFoundError(f"Project directory does not exist: {project}")
    if not config_path.is_file():
        raise FileNotFoundError(f"Project config not found: {config_path}")
    return project, config_path.resolve()


def refresh_manifest(
    *,
    project_dir: Path,
    config_path: Path,
    dry_run: bool = False,
) -> dict[str, object]:
    """先在内存中校验完整改名与 manifest 计划；预览不写文件，正式写出失败时回滚改名。"""
    config = load_project_config(config_path)

    if not project_dir.is_dir():
        raise FileNotFoundError(f"Project directory does not exist: {project_dir}")

    species = str(config["species"])
    protocol = get_protocol_policy(config).name
    raw_dir = project_dir / "01_raw"
    manifest_path = project_dir / "00_config" / "sample_manifest.tsv"


    fastq_renames = resolve_single_pair_fastq_collisions(raw_dir, apply=False)
    discovered = discover_samples(raw_dir, planned_renames=dict(fastq_renames))
    if not discovered:
        raise ValueError(f"No paired FASTQ samples found under: {raw_dir}")
    existing_rows = _read_existing_manifest_loose(manifest_path)
    rows = build_refreshed_rows(
        discovered.keys(), existing_rows, species=species, protocol=protocol
    )
    if not rows:
        raise ValueError("Manifest refresh produced no rows")


    new_text = render_tsv(rows, PROJECT_SAMPLE_MANIFEST_FIELDS)
    validated = _parse_project_sample_manifest(
        io.StringIO(new_text), path=manifest_path, species=species, protocol=protocol
    )
    map_discovered_samples_to_manifest(discovered, validated, protocol=protocol)
    changed = _canonical_project_manifest_rows(existing_rows) != rows
    backup_path: Path | None = None

    if not dry_run:
        original_manifest = manifest_path.read_bytes() if changed and manifest_path.exists() else None
        original_mode = stat.S_IMODE(manifest_path.stat().st_mode) if original_manifest is not None else None
        with _fastq_rename_transaction(raw_dir, fastq_renames):
            try:
                if changed:
                    if original_manifest is not None and existing_rows:
                        backup_path = manifest_path.with_name(f"{manifest_path.name}.bak")
                        # 不原地截断旧备份，也不沿既有链接写入其他文件。
                        atomic_write_bytes(
                            backup_path, original_manifest,
                            mode=original_mode,
                        )
                    write_tsv_atomic(manifest_path, rows, PROJECT_SAMPLE_MANIFEST_FIELDS, mode=original_mode)
            except BaseException:
                # 替换已经完成后仍可能收到中断；先恢复清单，再由外层逆序恢复文件名。
                if changed:
                    if original_manifest is None:
                        manifest_path.unlink(missing_ok=True)
                    elif not manifest_path.exists() or manifest_path.read_bytes() != original_manifest:
                        atomic_write_bytes(manifest_path, original_manifest, mode=original_mode)
                raise

    return {
        "project_dir": project_dir,
        "config_path": config_path,
        "raw_dir": raw_dir,
        "manifest_path": manifest_path,
        "species": species,
        "protocol": protocol,
        "discovered_count": len(discovered),
        "row_count": len(rows),
        "changed": changed,
        "backup_path": backup_path,
        "fastq_renames": fastq_renames,
        "rows": rows,
        "text": new_text,
    }


def init_project_directory(project_dir: str | Path, pipeline_root: str | Path) -> Path:
    """初始化项目骨架：目录、config 模板、manifest 头、Barcode Map 与 launcher symlink。"""
    if not str(project_dir).strip():
        raise ValueError("empty project path")
    root = Path(pipeline_root).expanduser().resolve()
    project = Path(project_dir).expanduser()
    project.mkdir(parents=True, exist_ok=True)
    project = project.resolve()
    if project == root:
        raise ValueError("the Pipeline installation directory cannot also be a project")
    for name in ("00_config", "01_raw", "02_work", "03_results", "04_logs", "05_tmp"):
        (project / name).mkdir(exist_ok=True)
    config = project / "00_config" / "config.yaml"
    if not config.exists():
        config.write_text(CONFIG_TEMPLATE_YAML, encoding="utf-8")
    manifest = project / "00_config" / "sample_manifest.tsv"
    if not manifest.exists():
        manifest.write_text(
            "sample_id\tdna_raw_sample\trna_sample\tspecies\tprotocol\tnotes\n",
            encoding="utf-8",
        )
    barcode = project / "00_config" / "Barcode_Map.csv"
    if not barcode.exists() and load_project_config(config)["analysis"]["protocol"] != "droplet":
        source = root / "resources" / "Barcode_Map.csv"
        if not source.is_file():
            raise FileNotFoundError(f"Barcode Map template is missing: {source}")
        shutil.copy2(source, barcode)
    launcher = project / "run_pipeline.sh"
    if not launcher.exists() and not launcher.is_symlink():
        launcher.symlink_to(root / "core" / "run_pipeline.sh")
    elif not launcher.is_symlink():
        print(f"WARNING: {launcher} already exists and was not replaced.")
    return project


def print_refresh_result(args, result) -> int:
    """refresh-manifest 的用户可见输出（预览/更新/未变化 + 改名计划）。"""
    action = "Previewed" if args.dry_run else ("Updated" if result["changed"] else "Unchanged")
    print(f"{action} sample manifest")
    print(f"  Project:  {result['project_dir']}")
    print(f"  Config:   {result['config_path']}")
    print(f"  Raw:      {result['raw_dir']}")
    print(f"  Manifest: {result['manifest_path']}")
    print(f"  Protocol: {result['protocol']}")
    print(f"  Species:  {result['species']}")
    print(f"  FASTQ samples: {result['discovered_count']}")
    print(f"  Manifest rows: {result['row_count']}")
    if result["fastq_renames"]:
        state = "Planned renames" if args.dry_run else "Normalized FASTQ pairs"
        print(f"  {state}: {len(result['fastq_renames'])} file(s)")
        for source, target in result["fastq_renames"]:
            rel_source = Path(source).relative_to(result["raw_dir"])
            print(f"    {rel_source} -> {target.name}")
    if result["backup_path"]:
        print(f"  Backup:   {result['backup_path']}")
    if args.dry_run:
        print("\n" + str(result["text"]), end="")
    return 0


# 输入与解复用：细胞识别和提交边界

def demux_source_revision(source_root: str | Path, *, calling_only: bool = False) -> str:
    """覆盖全部Rust源码；calling仅排除独立输出/replay，保留共享spool校验。"""
    root = Path(source_root).expanduser().resolve()
    main_rs = root / "src" / "main.rs"
    paths = [root / "Cargo.toml", root / "Cargo.lock", main_rs,
             root.parent.parent / DROPLET_DESIGN_RESOURCE]
    missing = [path for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "missing Rust source identity file(s): " + ", ".join(map(str, missing))
        )
    rust_sources = sorted((root / "src").rglob("*.rs"))
    # 新增共享源码默认进入两种身份，不以固定文件清单漏掉实际依赖。
    output_only = {"src/output.rs", "src/droplet/replay.rs"}
    paths.remove(main_rs)
    paths.extend(path for path in rust_sources
                 if not calling_only or path.relative_to(root).as_posix() not in output_only)
    digest = hashlib.sha256()
    for path in sorted(paths, key=lambda item: Path(os.path.relpath(item, root)).as_posix()):
        relative = Path(os.path.relpath(path, root)).as_posix().encode()
        data = path.read_bytes()
        digest.update(len(relative).to_bytes(4, "big"))
        digest.update(relative)
        digest.update(len(data).to_bytes(8, "big"))
        digest.update(data)
    return "sha256:" + digest.hexdigest()


def verify_demux_binary_revision(binary_path: str | Path, expected_revision: str) -> None:
    """校验 demux_rs --build-info 内嵌的源码指纹与当前源码完全一致。"""
    binary = Path(binary_path).expanduser().resolve()
    if not binary.is_file():
        raise FileNotFoundError(f"demux binary missing: {binary}")
    completed = subprocess.run(
        [str(binary), "--build-info"],
        check=True,
        capture_output=True,
        text=True,
        timeout=30,
    )
    payload = json.loads(completed.stdout)
    if not isinstance(payload, dict):
        raise ValueError("demux --build-info did not return a JSON object")
    actual = payload.get("source_revision")
    if actual != expected_revision:
        raise ValueError(
            f"demux binary source revision mismatch: {actual!r} != {expected_revision!r}"
        )


DEMUX_REPORT_SCHEMA_VERSION = 3


_REPORT_FIELDS = {
    "schema_version", "build_version", "source_revision", "mode",
    "retention_policy", "retention_threshold", "input_fastq_pairs",
    "input_read_pairs", "usable_read_pairs", "matched_reads", "read_fates", "samples",
}


_FATE_FIELDS = {"dna_assigned", "rna_assigned", "ambiguous", "unmatched", "short"}


_SAMPLE_FIELDS = {
    "sample_name", "cell_sample_id", "plate_id", "cell_order", "dna_barcode",
    "rna_barcode", "dna_read_count", "rna_read_count", "dna_status", "rna_status",
}


_OUTPUT_STATUSES = {"Pass", "Zero_Output", "Low_Reads_Removed"}


class DemuxReadFates(TypedDict):
    dna_assigned: int
    rna_assigned: int
    ambiguous: int
    unmatched: int
    short: int


class DemuxReportSample(TypedDict):
    sample_name: str
    cell_sample_id: str
    plate_id: str
    cell_order: str
    dna_barcode: str
    rna_barcode: str
    dna_read_count: int
    rna_read_count: int
    dna_status: str
    rna_status: str


class DemuxReport(TypedDict):
    schema_version: int
    build_version: str
    source_revision: str
    mode: str
    retention_policy: str
    retention_threshold: int
    input_fastq_pairs: int
    input_read_pairs: int
    usable_read_pairs: int
    matched_reads: int
    read_fates: DemuxReadFates
    samples: list[DemuxReportSample]


def cell_status(reads: int, threshold: int) -> str:
    """demux cell 状态的唯一推导：0 读 Zero_Output、低于阈值 Low_Reads_Removed、其余 Pass。"""
    return "Zero_Output" if reads == 0 else "Low_Reads_Removed" if reads < threshold else "Pass"


def read_demux_report(path: str | Path) -> DemuxReport:
    """完整加载一份 schema-3 report，不为缺失字段做任何推断。"""

    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(f"Demux report is missing: {source}")
    try:
        raw = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"Invalid demux report JSON {source}: {exc}") from exc
    return _parse_demux_report(raw, source)


def _parse_demux_report(raw: object, source: Path) -> DemuxReport:
    report = _object(raw, _REPORT_FIELDS, "Demux report", source)

    schema = _integer(report["schema_version"], "schema_version", source)
    if schema != DEMUX_REPORT_SCHEMA_VERSION:
        raise ValueError(f"Unsupported demux report schema_version={schema}; expected 3: {source}")
    for field in ("build_version", "source_revision", "mode", "retention_policy"):
        _text(report[field], field, source)
    mode = cast(str, report["mode"])
    if mode not in {"dna-only", "dna-only-taps", "dna-rna", "dna-only-droplet"}:
        raise ValueError(f"Demux report mode is unsupported: {mode!r}: {source}")
    if report["retention_policy"] != "matched_read_pairs":
        raise ValueError(f"Demux report retention_policy is unsupported: {source}")

    counts = {
        field: _integer(report[field], field, source)
        for field in (
            "retention_threshold", "input_fastq_pairs", "input_read_pairs",
            "usable_read_pairs", "matched_reads",
        )
    }
    if counts["input_fastq_pairs"] < 1:
        raise ValueError(f"input_fastq_pairs must be at least 1: {source}")
    fate = _object(report["read_fates"], _FATE_FIELDS, "read_fates", source)
    fates = {key: _integer(fate[key], f"read_fates.{key}", source) for key in _FATE_FIELDS}
    if counts["input_read_pairs"] != sum(fates.values()):
        raise ValueError(f"Demux report input/read-fate accounting is inconsistent: {source}")
    if counts["usable_read_pairs"] != counts["input_read_pairs"] - fates["short"]:
        raise ValueError(f"Demux report usable/input/short accounting is inconsistent: {source}")
    if counts["matched_reads"] != fates["dna_assigned"] + fates["rna_assigned"]:
        raise ValueError(f"Demux report matched/read-fate accounting is inconsistent: {source}")
    if mode != "dna-rna" and fates["rna_assigned"]:
        raise ValueError(f"DNA-only demux report has RNA-assigned reads: {source}")

    samples = report["samples"]
    if not isinstance(samples, list) or not samples:
        raise ValueError(f"Demux report samples must be a non-empty array: {source}")
    sample_names: set[str] = set()
    cell_ids: set[str] = set()
    dna_total = rna_total = 0
    threshold = counts["retention_threshold"]
    for index, value in enumerate(samples):
        label = f"samples[{index}]"
        sample = _object(value, _SAMPLE_FIELDS, label, source)
        for field in (
            "sample_name", "cell_sample_id", "plate_id", "cell_order", "dna_barcode",
            "dna_status", "rna_status",
        ):
            _text(sample[field], f"{label}.{field}", source, empty=(field == "plate_id" and mode == "dna-only-droplet"))
        _text(sample["rna_barcode"], f"{label}.rna_barcode", source, empty=True)
        sample_name, cell_id, plate = (
            cast(str, sample[key]) for key in ("sample_name", "cell_sample_id", "plate_id")
        )
        if SAMPLE_ID_RE.fullmatch(sample_name) is None or SAMPLE_ID_RE.fullmatch(cell_id) is None:
            raise ValueError(f"{label} has an invalid sample identity: {source}")
        if cell_id != cell_identity(sample_name, "droplet" if mode == "dna-only-droplet" else "cabernet", plate, sample["dna_barcode"]):
            raise ValueError(f"{label}.cell_sample_id does not match sample_name/plate_id: {source}")
        if cell_id in cell_ids:
            raise ValueError(f"Demux report repeats cell_sample_id={cell_id!r}: {source}")
        cell_ids.add(cell_id)
        sample_names.add(sample_name)

        dna_reads = _integer(sample["dna_read_count"], f"{label}.dna_read_count", source)
        rna_reads = _integer(sample["rna_read_count"], f"{label}.rna_read_count", source)
        dna_status, rna_status = cast(str, sample["dna_status"]), cast(str, sample["rna_status"])
        if dna_status not in _OUTPUT_STATUSES or dna_status != cell_status(dna_reads, threshold):
            raise ValueError(f"{label}.dna_status does not match dna_read_count: {source}")
        if mode != "dna-rna":
            if rna_reads or rna_status != "Not_Used":
                raise ValueError(f"DNA-only {label} has active RNA output fields: {source}")
        elif not sample["rna_barcode"]:
            raise ValueError(f"DNA/RNA {label} has an empty rna_barcode: {source}")
        elif rna_status not in _OUTPUT_STATUSES or rna_status != cell_status(rna_reads, threshold):
            raise ValueError(f"{label}.rna_status does not match rna_read_count: {source}")
        dna_total += dna_reads
        rna_total += rna_reads

    if len(sample_names) != 1:
        raise ValueError(f"Demux report must contain exactly one project sample identity: {source}")
    if dna_total != fates["dna_assigned"] or rna_total != fates["rna_assigned"]:
        raise ValueError(f"Demux report sample/read-fate accounting is inconsistent: {source}")
    return cast(DemuxReport, report)


def _validate_droplet_stage_paths(inputs: Sequence[str | Path], spool: Path,
                                  counts: Path, structure: Path) -> None:
    """最终提交路径也须隔离；原生程序只看得到本次临时输出，无法代替这层核对。"""
    protected: dict[Path, tuple[int, int]] = {}
    def protect_input_entries(path: Path, current: Path) -> Path:
        parts = path.parts
        if path.is_absolute():
            current, parts = Path(path.anchor), parts[1:]
        for index, part in enumerate(parts):
            if part == "..":
                current = current.parent
                continue
            current = current / part
            # 连输入入口的链接/祖先也保护，不能删掉指向外部raw的入口后再尝试读取。
            metadata = current.lstat()
            protected[current] = (metadata.st_dev, metadata.st_ino)
            if stat.S_ISLNK(metadata.st_mode):
                current.resolve(strict=True)  # 保留循环/无效链接拒绝，再逐跳保护中间入口。
                current = protect_input_entries(current.readlink(), current.parent)
                metadata = current.stat()
                protected[current] = (metadata.st_dev, metadata.st_ino)
            if index < len(parts) - 1 and not stat.S_ISDIR(metadata.st_mode):
                raise NotADirectoryError(f"Droplet stage input ancestor must be a directory: {path}")
        return current
    for source in inputs:
        path = Path(source)
        expected = path.resolve(strict=True)
        expected_metadata = expected.stat()
        actual = protect_input_entries(path, Path.cwd())
        actual_metadata = actual.stat()
        if (actual_metadata.st_dev, actual_metadata.st_ino) != (expected_metadata.st_dev, expected_metadata.st_ino):
            raise RuntimeError(f"Droplet stage input changed while resolving: {source}")
        protected[actual] = (actual_metadata.st_dev, actual_metadata.st_ino)

    outputs: list[tuple[Path, bool, tuple[int, int] | None, set[tuple[int, int]]]] = []
    for destination, directory in ((spool, True), (counts, False), (structure, False)):
        absolute = destination.absolute()
        resolved = Path(absolute.anchor).resolve(strict=True)
        parts = absolute.parts[1:]
        for index, part in enumerate(parts):
            if part == "..":
                resolved = resolved.parent
                continue
            resolved = resolved / part
            try:
                metadata = resolved.lstat()
            except FileNotFoundError:
                continue
            last = index == len(parts) - 1
            if stat.S_ISLNK(metadata.st_mode):
                if last:
                    raise ValueError(f"Droplet stage output must not be a symlink: {destination}")
                resolved = resolved.resolve(strict=True)
                metadata = resolved.stat()
            if not last and not stat.S_ISDIR(metadata.st_mode):
                raise ValueError(f"Droplet stage output ancestor must be a directory: {destination}")
        identity = None
        if resolved.exists():
            metadata = resolved.stat()
            identity = (metadata.st_dev, metadata.st_ino)
            valid_type = stat.S_ISDIR(metadata.st_mode) if directory else stat.S_ISREG(metadata.st_mode)
            if not valid_type or (not directory and metadata.st_nlink != 1):
                raise ValueError(f"Droplet stage output must have an exclusive valid type: {destination}")
        if (identity is not None and identity in protected.values()) or any(
                resolved == source or (directory and source.is_relative_to(resolved)) for source in protected):
            raise ValueError(f"Droplet stage output overlaps an input: {destination}")
        parents = set()
        for parent in resolved.parents:
            try:
                metadata = parent.stat()
            except FileNotFoundError:
                continue
            if not stat.S_ISDIR(metadata.st_mode):
                raise ValueError(f"Droplet stage output ancestor must be a directory: {destination}")
            parents.add((metadata.st_dev, metadata.st_ino))
        for other, other_directory, other_identity, other_parents in outputs:
            if (resolved == other or (directory and other.is_relative_to(resolved))
                    or (other_directory and resolved.is_relative_to(other))
                    or (identity is not None and identity == other_identity)
                    or (directory and identity is not None and identity in other_parents)
                    or (other_directory and other_identity is not None and other_identity in parents)):
                raise ValueError(f"Droplet stage outputs overlap: {destination}")
        outputs.append((resolved, directory, identity, parents))

    # 尚不存在的输出没有inode；在同一实际命名空间中试建，才能识别大小写/Unicode等别名。
    groups: dict[tuple[int, int], tuple[Path, list[tuple[Path, bool]]]] = {}
    for resolved, directory, _, _ in outputs:
        parent = resolved.parent
        while not parent.exists():
            parent = parent.parent
        metadata = parent.stat()
        members = groups.setdefault((metadata.st_dev, metadata.st_ino), (parent, []))[1]
        members.append((resolved.relative_to(parent), directory))
    for parent, members in groups.values():
        if len(members) < 2:
            continue
        parent = parent.resolve(strict=True)
        parent_metadata = parent.stat()
        probe = Path(tempfile.mkdtemp(prefix=".alopex_output_names.", dir=parent))
        probe_metadata = probe.stat()
        failure = None
        try:
            planned: list[tuple[Path, bool]] = []
            for suffix, directory in members:
                path = probe / suffix
                path.parent.mkdir(parents=True, exist_ok=True)
                if directory:
                    path.mkdir(exist_ok=True)
                else:
                    with path.open("xb"):
                        pass
                for other, other_directory in planned:
                    if (path.samefile(other)
                            or (directory and any(path.samefile(p) for p in other.parents))
                            or (other_directory and any(other.samefile(p) for p in path.parents))):
                        raise ValueError(f"Droplet stage output names overlap on the target filesystem: {suffix}")
                planned.append((path, directory))
        except BaseException as error:
            failure = error
            raise
        finally:
            try:
                current_parent, current_probe = parent.lstat(), probe.lstat()
                if (stat.S_ISLNK(current_parent.st_mode) or stat.S_ISLNK(current_probe.st_mode)
                        or (current_parent.st_dev, current_parent.st_ino) != (parent_metadata.st_dev, parent_metadata.st_ino)
                        or (current_probe.st_dev, current_probe.st_ino) != (probe_metadata.st_dev, probe_metadata.st_ino)):
                    raise OSError(f"Droplet output namespace probe changed before cleanup: {probe}")
                shutil.rmtree(probe)
            except OSError as cleanup_error:
                if failure is None:
                    raise
                logger.warning("Droplet output namespace probe cleanup failed; preserved %s: %s", probe, cleanup_error)


def stage_droplet_reads(r1: Sequence[str | Path], r2: Sequence[str | Path], spool: str | Path,
                        counts: str | Path, structure: str | Path, *, binary: str | Path,
                        decoder_python: str | Path, threads: int) -> None:
    """配对扫描一次并提交原始条码暂存及完整计数，任一解码或写出失败不提交暂存边界。"""
    if not r1 or len(r1) != len(r2) or threads < 1:
        raise ValueError("Droplet stage requires ordered paired FASTQs and positive threads")
    spool, counts, structure = map(Path, (spool, counts, structure))
    _validate_droplet_stage_paths([*r1, *r2, binary, decoder_python], spool, counts, structure)
    spool.parent.mkdir(parents=True, exist_ok=True)
    counts.parent.mkdir(parents=True, exist_ok=True)
    structure.parent.mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix=".droplet_stage.", dir=project_temporary_parent(spool)))
    owned_root = root.stat()
    preserve_backup = False
    try:
        for destination in (spool, counts, structure):
            if destination.parent.stat().st_dev != root.stat().st_dev:
                raise OSError("Droplet stage outputs and temporary directory must share a filesystem")
        staged = root / "spool"
        command = [str(binary), "droplet-stage", "--spool", str(staged), "--counts-output", str(root / "counts.tsv"),
                   "--count-metrics", str(root / "structure.json"), "--decoder-python", str(decoder_python), "--threads", str(threads)]
        for flag, paths in (("--r1", r1), ("--r2", r2)):
            for path in paths:
                command.extend((flag, str(path)))
        subprocess.run(command, check=True)
        current_root = root.lstat()
        if (current_root.st_dev, current_root.st_ino) != (owned_root.st_dev, owned_root.st_ino):
            raise RuntimeError(f"Droplet stage temporary directory changed: {root}")
        plan = json.loads((staged / "manifest.json").read_text())
        metrics = json.loads((root / "structure.json").read_text())
        if plan["version"] != 1 or plan["input_fastq_pairs"] != len(r1) or plan["pairs"] != metrics["input_reads"]:
            raise ValueError("Droplet staged pair accounting differs from structure counts")
        _validate_droplet_stage_paths([*r1, *r2, binary, decoder_python], spool, counts, structure)
        backup = root / "previous"
        backup.mkdir()
        saved: dict[Path, tuple[Path, tuple[int, int, int, int]]] = {}
        installed: dict[Path, tuple[int, int, int, int]] = {}
        def path_identity(path: Path) -> tuple[int, int, int, int]:
            metadata = path.lstat()
            return metadata.st_dev, metadata.st_ino, metadata.st_size, metadata.st_mtime_ns
        try:
            # 先撤含manifest的spool，再换两个成员；spool最后提交，回滚也最后恢复。
            for name, destination in (("spool", spool), ("counts.tsv", counts), ("structure.json", structure)):
                if destination.exists():
                    previous = backup / name
                    # 先登记再rename，覆盖rename已生效、Python尚未记录时到达的中断。
                    saved[destination] = (previous, path_identity(destination))
                    os.replace(destination, previous)
            for source, destination in ((root / "counts.tsv", counts),
                                        (root / "structure.json", structure), (staged, spool)):
                installed[destination] = path_identity(source)
                os.replace(source, destination)
        except BaseException as error:
            try:
                for destination, identity in reversed(installed.items()):
                    if not destination.exists():
                        continue
                    if path_identity(destination) != identity:
                        raise OSError(f"Droplet stage rollback destination changed: {destination}")
                    if destination.is_dir():
                        shutil.rmtree(destination)
                    else:
                        destination.unlink()
                for destination in (counts, structure, spool):
                    if destination not in saved:
                        continue
                    previous, identity = saved[destination]
                    if not previous.exists() and not previous.is_symlink():
                        if (not destination.exists() or path_identity(destination) != identity):
                            raise OSError(f"Droplet stage backup and original destination changed: {destination}")
                        continue
                    if path_identity(previous) != identity:
                        raise OSError(f"Droplet stage backup changed: {previous}")
                    if destination.exists() or destination.is_symlink():
                        raise FileExistsError(f"Droplet stage rollback destination is occupied: {destination}")
                    os.replace(previous, destination)
            except BaseException as recovery_error:
                preserve_backup = True
                message = f"Droplet stage rollback failed; original backups retained at {root}"
                logger.error(message)
                error.add_note(message)
                raise error from recovery_error
            raise
    finally:
        if not preserve_backup:
            try:
                current_root = root.lstat()
                if (current_root.st_dev, current_root.st_ino) != (owned_root.st_dev, owned_root.st_ino):
                    logger.warning("Refusing to clean replaced Droplet stage directory: %s", root)
                else:
                    shutil.rmtree(root)
            except FileNotFoundError:
                pass
            except OSError as cleanup_error:
                logger.warning("Droplet stage temporary cleanup failed at %s: %s", root, cleanup_error)


def droplet_output_disk_mb(spool: str | Path) -> int:
    """按已提交暂存压缩量预留输出 scratch；给副本和压缩差异留余量，不进入分析配置。"""
    spool = Path(spool)
    plan = json.loads((spool / "manifest.json").read_text())
    total = sum((spool / f"{i:04}.pack").stat().st_size for i in range(plan["shards"]))
    return max(1024, math.ceil(total * 2 / 1024 ** 2) + 512)


def copy_demux_fastq(source: str | Path, destination: str | Path, workers: int = 1) -> dict[str, object]:
    """按完整文件物化 FASTQ；有界复制并尽力释放已消费缓存，失败仍由调用事务回滚。"""
    started = time.monotonic()
    source, destination = Path(source), Path(destination)
    if source.is_symlink() or destination.is_symlink() or source.resolve() == destination.resolve():
        raise ValueError("FASTQ materialization requires distinct real directories")
    destination.mkdir(parents=True, exist_ok=True)
    if any(destination.iterdir()):
        raise ValueError("FASTQ materialization destination must be empty")
    paths = []
    for path in source.rglob("*"):
        info = path.lstat()
        if stat.S_ISDIR(info.st_mode):
            continue
        if not stat.S_ISREG(info.st_mode) or not path.name.endswith(".fastq.gz"):
            raise ValueError(f"Unexpected demux output: {path}")
        paths.append((path, path.relative_to(source), info.st_size))
    paths.sort(key=lambda record: (-record[2], str(record[1])))
    for parent in sorted({relative.parent for _, relative, _ in paths}, key=lambda p: (len(p.parts), str(p))):
        (destination / parent).mkdir(parents=True, exist_ok=True)
    advise = getattr(os, "posix_fadvise", None)
    discard = getattr(os, "POSIX_FADV_DONTNEED", None)
    streamed = callable(advise) and discard is not None
    limit = 4 if streamed else 8
    count = min(limit, effective_cpu_count(workers))
    def release(fd: int, offset: int, length: int) -> None:
        try:
            advise(fd, offset, length, discard)
        except OSError:
            pass
    def one(record: tuple[Path, Path, int]) -> int:
        path, relative, size = record
        target = destination / relative
        if streamed:
            copied = 0
            advised = 0
            buffer = bytearray(1024 * 1024)
            with path.open("rb", buffering=0) as incoming, target.open("xb", buffering=0) as outgoing:
                while n := incoming.readinto(buffer):
                    view = memoryview(buffer)[:n]
                    while view:
                        written = os.write(outgoing.fileno(), view)
                        if written <= 0:
                            raise OSError(f"Incomplete demux FASTQ write: {target}")
                        view = view[written:]
                    copied += n
                    if copied - advised >= 8 * 1024 * 1024:
                        release(incoming.fileno(), advised, copied - advised)
                        release(outgoing.fileno(), advised, copied - advised)
                        advised = copied
                release(incoming.fileno(), advised, 0)
                release(outgoing.fileno(), advised, 0)
            if copied != size:
                raise OSError(f"Incomplete demux FASTQ read: {path}")
        else:
            shutil.copyfile(path, target)
        if target.stat().st_size != size:
            raise OSError(f"Incomplete demux FASTQ copy: {path}")
        return size
    with ThreadPoolExecutor(max_workers=count) as pool:
        copied = sum(pool.map(one, paths))
    return {"files": len(paths), "bytes": copied, "copy_seconds": time.monotonic() - started,
            "copy_workers": count, "copy_strategy": "bounded_stream" if streamed else "native_copy"}


@lru_cache(maxsize=1)
def droplet_design_barcodes() -> frozenset[str]:
    """读取固定内容身份的 DD-MET5 设计条码，供 calling 约束身份，不依赖 RNA 名单。"""
    path = PIPELINE_ROOT / DROPLET_DESIGN_RESOURCE
    data = path.read_bytes()
    if sha256_bytes(data) != DROPLET_DESIGN_SHA256:
        raise ValueError("DD-MET5 design barcode resource SHA256 mismatch")
    return frozenset(gzip.decompress(data).decode("ascii").splitlines())


def call_droplet_cells(
    counts_path: str | Path, output: str | Path, metrics_path: str | Path,
    *, demux_binary: str | Path | None = None, threads: int = 1,
    spool: str | Path | None = None,
) -> None:
    """按固定谷底及质量/共享分子筛查调用，再限定为设计条码并记录被排除的原始身份。"""
    rows = []
    seen = set()
    qualities = {}
    with Path(counts_path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["barcode", "count", "below_q20", "q30"]:
            raise ValueError("Droplet counts require barcode/count/below_q20/q30 columns")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Droplet count row columns disagree with header at {counts_path}:{reader.line_num}")
            barcode, count = row["barcode"], int(row["count"])
            if not re.fullmatch(r"[ACGT]{17}", barcode) or count < 1 or barcode in seen:
                raise ValueError("Invalid or duplicate Droplet barcode count")
            low, high = ([int(n) for n in row[key].split(",")] for key in ("below_q20", "q30"))
            if len(low) != 17 or len(high) != 17 or any(a < 0 or b < 0 or a + b > count for a, b in zip(low, high)):
                raise ValueError("Invalid Droplet per-base quality counts")
            seen.add(barcode)
            if "C" not in barcode and count > 100:
                rows.append((barcode, count))
                qualities[barcode] = (low, high)
    rows.sort(key=lambda row: (-row[1], row[0]))
    metrics = {"method": "log10_count_plus_one_valley_200bins_23smooth_20distance", "eligible_barcodes": len(rows), "called_cells": 0}
    failure = "Droplet calling requires two separated peaks among AGT barcodes with count >100"
    if rows and rows[0][1] != rows[-1][1]:
        values = [math.log10(count + 1) for _, count in rows]
        lower, upper = min(values), max(values)
        width = (upper - lower) / 200
        hist = [0] * 200
        for value in values:
            hist[min(199, int((value - lower) / width))] += 1
        smooth = [sum(hist[max(0, i - 11):min(200, i + 12)]) / 23 for i in range(200)]
        centers = [lower + (i + 0.5) * width for i in range(200)]
        peaks = []
        i = 1
        while i < 199:
            start = i
            while i < 199 and abs(smooth[i + 1] - smooth[start]) <= 1e-12:
                i += 1
            if smooth[start] > smooth[start - 1] + 1e-12 and i < 199 and smooth[i] > smooth[i + 1] + 1e-12:
                peaks.append((start + i) // 2)
            i += 1
        selected = []
        for peak in sorted(peaks, key=lambda p: (-smooth[p], p)):
            if all(abs(peak - other) >= 20 for other in selected):
                selected.append(peak)
        metrics.update(histogram=hist, smoothed=smooth, bin_centers=centers, peaks=selected)
        if len(selected) >= 2:
            left, right = sorted(selected[:2])
            valley = smooth[left + 1:right]
            low = min(valley)
            minima = [i + left + 1 for i, value in enumerate(valley) if abs(value - low) <= 1e-12]
            start = end = int(minima[0])
            for pos in minima[1:]:
                if pos != end + 1:
                    break
                end = int(pos)
            valley_bin = (start + end) // 2
            threshold = float(10 ** centers[valley_bin] - 1)
            called = [(bc, n) for bc, n in rows if n >= threshold]
            metrics.update(valley_bin=valley_bin, threshold=threshold, called_cells=len(called), called_read_pairs=sum(n for _, n in called))
            failure = "" if called and low < min(smooth[left], smooth[right]) - 1e-12 else failure
        else:
            called = []
    else:
        called = []
    metrics["status"] = "failed" if failure else "screening"
    metrics["failure"] = failure
    Path(metrics_path).parent.mkdir(parents=True, exist_ok=True)
    Path(metrics_path).write_text(json.dumps(metrics, indent=2) + "\n")
    if failure:
        Path(output).unlink(missing_ok=True)
        raise ValueError(failure)
    metrics["initial_called_cells"] = len(called)
    metrics["initial_called_read_pairs"] = metrics["called_read_pairs"]
    screen = {"method": "unique_hamming1_quality_complete_parent_scan", "parent_min_ratio": 5,
              "child_min_below_q20_fraction": 0.2, "low_quality_min_enrichment": 5,
              "independent_support": "child Q30 reads at the differing base below original valley threshold",
              "minimum_shared_umis": 3, "minimum_shared_inserts": 3,
              "signature": "exact_umi12_insert32", "child_signature_cap": 4096, "parent_scan": "complete",
              "evaluated": [], "removed_barcodes": 0, "removed_exact_read_pairs": 0}
    metrics["barcode_error_screen"] = screen
    counts = dict(called)
    candidates = []
    for child, count in called:
        neighbors = [(neighbor, pos) for pos, base in enumerate(child) for alternative in "AGT"
                     if alternative != base
                     if counts.get(neighbor := child[:pos] + alternative + child[pos + 1:], 0) > count]
        if not neighbors:
            continue
        decision = {"barcode": child, "count": count, "higher_neighbors": sorted(bc for bc, _ in neighbors),
                    "decision": "retain", "reason": "multiple_higher_neighbors"}
        screen["evaluated"].append(decision)
        if len(neighbors) != 1:
            continue
        parent, pos = neighbors[0]
        low, high = qualities[child]
        parent_low = qualities[parent][0][pos]
        decision.update(parent=parent, position=pos + 1, parent_count=counts[parent],
                        below_q20=low[pos], q30=high[pos], parent_below_q20=parent_low)
        if counts[parent] < 5 * count:
            decision["reason"] = "insufficient_abundance_ratio"
        elif high[pos] >= threshold:
            decision["reason"] = "independent_q30_support"
        elif low[pos] < 0.2 * count or low[pos] * counts[parent] < 5 * parent_low * count:
            decision["reason"] = "insufficient_quality_evidence"
        else:
            decision["reason"] = "pending_molecule_evidence"
            candidates.append(decision)
    try:
        design = droplet_design_barcodes()
        if candidates:
            if spool is None or demux_binary is None or threads < 1:
                raise ValueError("Droplet error screening requires staged reads and demux binary for molecule evidence")
            selected = [{key: r[key] for key in ("barcode", "parent")} for r in candidates]
            with tempfile.TemporaryDirectory(prefix=".droplet_evidence.", dir=project_temporary_parent(metrics_path)) as tmp:
                pair_file, evidence_file = Path(tmp) / "pairs.json", Path(tmp) / "molecules.json"
                pair_file.write_text(json.dumps(selected))
                command = [str(demux_binary), "--mode", "dna-only-droplet", "--threads", str(threads),
                           "--evidence-pairs", str(pair_file), "--evidence-output", str(evidence_file)]
                command.extend(["--droplet-spool", str(spool)])
                subprocess.run(command, check=True)
                evidence = json.loads(evidence_file.read_text())
            if (evidence["signature"] != screen["signature"] or evidence["child_cap"] != screen["child_signature_cap"]
                    or evidence["parent_scan"] != "complete" or set(evidence["pairs"]) != {r["barcode"] for r in selected}):
                raise ValueError("Droplet molecule evidence identity mismatch")
            proposed = set()
            for decision in candidates:
                child, parent = decision["barcode"], decision["parent"]
                values = evidence["pairs"][child]
                shared = values["shared_signatures"]
                sampled, child_reads, parent_reads = (values[key] for key in
                    ("child_sampled_signatures", "child_eligible_reads", "parent_eligible_reads"))
                if (values["parent"] != parent or any(type(n) is not int for n in (sampled, child_reads, parent_reads))
                        or not 0 <= len(shared) <= sampled <= min(evidence["child_cap"], child_reads)
                        or not 0 <= child_reads <= counts[child] or not len(shared) <= parent_reads <= counts[parent]
                        or type(values["child_saturated"]) is not bool
                        or (values["child_saturated"] and (sampled != evidence["child_cap"] or child_reads <= sampled))
                        or any(re.fullmatch(r"[0-9a-f]{22}", key) is None for key in shared) or len(set(shared)) != len(shared)):
                    raise ValueError("Invalid Droplet molecule evidence signatures or pair identity")
                shared = sorted(shared)
                umis = {key[:6] for key in shared}
                inserts = {key[6:] for key in shared}
                decision.update(child_sampled_signatures=sampled, child_eligible_reads=child_reads,
                                child_saturated=values["child_saturated"], parent_eligible_reads=parent_reads,
                                shared_signatures=len(shared), shared_umis=len(umis), shared_inserts=len(inserts),
                                shared_examples=shared[:3])
                decision["reason"] = "insufficient_shared_molecules"
                if len(umis) >= 3 and len(inserts) >= 3:
                    proposed.add(child)
            removed = set()
            for decision in candidates:
                if decision["barcode"] in proposed:
                    if decision["parent"] in proposed:
                        decision["reason"] = "parent_is_error_candidate"
                    else:
                        decision.update(decision="remove", reason="quality_and_shared_molecules")
                        removed.add(decision["barcode"])
            screen.update(removed_barcodes=len(removed), removed_exact_read_pairs=sum(counts[bc] for bc in removed))
            called = [(bc, n) for bc, n in called if bc not in removed]
        screen.update(molecule_checked=len(candidates),
                      insufficient_molecule_evidence=sum(r["reason"] == "insufficient_shared_molecules" for r in candidates),
                      child_signatures_capped=sum(r["child_saturated"] for r in candidates))
        design_called = {bc for bc, _ in called if bc in design}
        excluded = []
        for barcode, count in called:
            if barcode in design:
                continue
            neighbors = sorted(neighbor for pos, base in enumerate(barcode) for alt in "AGT" if alt != base
                               if (neighbor := barcode[:pos] + alt + barcode[pos + 1:]) in design_called)
            excluded.append({"barcode": barcode, "count": count, "called_design_neighbors": neighbors,
                             "read_assignment": "unique_hamming1" if len(neighbors) == 1 else
                                                "ambiguous" if neighbors else "unmatched"})
        metrics["barcode_design_filter"] = {
            "method": "called_intersect_fixed_design_preserve_uncalled_design",
            "resource": DROPLET_DESIGN_RESOURCE, "sha256": DROPLET_DESIGN_SHA256,
            "design_barcodes": len(design), "called_before_design_filter": len(called),
            "excluded_candidates": len(excluded), "excluded_exact_read_pairs": sum(r["count"] for r in excluded),
            "excluded": excluded,
        }
        called = [(bc, n) for bc, n in called if bc in design_called]
        if not called:
            raise ValueError("Droplet calling retained no cells within the fixed DD-MET5 design")
        metrics.update(status="complete", called_cells=len(called), called_read_pairs=sum(n for _, n in called))
    except Exception as exc:
        metrics.update(status="failed", failure=f"Droplet calling failed: {exc}")
        Path(metrics_path).write_text(json.dumps(metrics, indent=2) + "\n")
        Path(output).unlink(missing_ok=True)
        raise
    Path(metrics_path).write_text(json.dumps(metrics, indent=2) + "\n")
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    content = "DNA_Barcode,Cell_Order\n" + "".join(f"{bc},{i}\n" for i, (bc, _) in enumerate(called, 1))
    if not destination.exists() or destination.read_text() != content:
        atomic_write_bytes(destination, content.encode("utf-8"))


def _cell_order_key(value: str) -> tuple[int, int | str]:
    """与 Rust i64 排序一致；非 ASCII 数字或越界值按原字符串排序，不改写身份。"""
    match = re.fullmatch(r"([+-]?)([0-9]+)", value)
    if match:
        # 去前导零只服务数值比较，避免 Python 超长十进制转换限制改变 Rust 的有效整数集合。
        digits = match[2].lstrip("0") or "0"
        if len(digits) <= 19:
            number = int(match[1] + digits)
            if -(1 << 63) <= number < (1 << 63):
                return (0, number)
    return (1, value)


def cell_identity(project: str, protocol: str, plate: str, barcode: str) -> str:
    """按孔板或液滴协议构造唯一 cell ID。"""
    if protocol == "droplet":
        if plate or not re.fullmatch(r"[AGT]{17}", barcode):
            raise ValueError("Droplet requires a 17 bp AGT barcode and empty plate_id")
        return f"{project}_{barcode}"
    if not plate or SAMPLE_ID_RE.fullmatch(plate) is None:
        raise ValueError("Plate protocol requires a valid plate_id")
    return f"{project}_{plate}"


def write_initial_cell_manifest(
    project_manifest: str,
    barcode_map: str,
    out_tsv: str,
    *,
    species: str,
    protocol: str,
    called_cells: Sequence[Sequence[str]] = (),
    project_sample_id: str | None = None,
) -> None:
    """按孔板 map 或 called cells 展开身份，可只写指定的已登记 Droplet 文库。"""
    project_rows = read_project_sample_manifest(
        Path(project_manifest), species=species, protocol=protocol
    )
    if project_sample_id is not None:
        if protocol != "droplet" or project_sample_id not in {r["sample_id"] for r in project_rows}:
            raise ValueError("Per-library cell manifest requires a registered Droplet sample")
        project_rows = [row for row in project_rows if row["sample_id"] == project_sample_id]
    called = dict(called_cells)
    if protocol == "droplet":
        if len(called) != len(called_cells) or set(called) != {r["sample_id"] for r in project_rows}:
            raise ValueError("Called-cell libraries must exactly match project manifest")
    elif called:
        raise ValueError("Called-cell inputs require droplet protocol")
    rows: list[CellIdentityRow] = []
    seen: set[str] = set()
    for project in project_rows:
        project_id = project["sample_id"]
        barcode_rows = validate_barcode_map_for_protocol(Path(called[project_id] if protocol == "droplet" else barcode_map), get_protocol(protocol))
        for barcode in barcode_rows:
            cell_id = cell_identity(project_id, protocol, barcode["plate_id"], barcode["dna_barcode"])
            if SAMPLE_ID_RE.fullmatch(cell_id) is None:
                raise ValueError(
                    f"Canonical cell Sample ID contains unsupported characters: {cell_id!r}. "
                    "Check sample_id and Barcode_Map.csv PlateID."
                )
            if cell_id in seen:
                raise ValueError(f"Canonical cell Sample ID collision: {cell_id!r}")
            seen.add(cell_id)
            rows.append(
                {
                    "sample_id": cell_id,
                    "project_sample_id": project_id,
                    "species": project["species"],
                    "protocol": project["protocol"],
                    "plate_id": barcode["plate_id"],
                    "cell_order": barcode["cell_order"],
                    "dna_barcode": barcode["dna_barcode"],
                    "rna_barcode": barcode["rna_barcode"],
                    "dna_raw_sample": project["dna_raw_sample"],
                    "rna_sample": project["rna_sample"],
                    "notes": project["notes"],
                }
            )

    def sort_key(row: Mapping[str, str]) -> tuple[object, ...]:
        raw_order = str(row.get("cell_order", "") or "").strip()
        return (row.get("project_sample_id", ""), _cell_order_key(raw_order), row.get("plate_id", ""))

    write_tsv_atomic(
        out_tsv,
        sorted(rows, key=sort_key),
        INITIAL_CELL_MANIFEST_FIELDS,
    )


def read_cell_identity_manifest(path: str | Path) -> dict[str, CellIdentityRow]:
    """读取 demux 前的 cell identity 快照并按 cell ID 索引。"""
    source = Path(path)
    with source.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        _require_header(fields, INITIAL_CELL_MANIFEST_FIELDS, source=source, label="initial cell manifest")
        rows: dict[str, CellIdentityRow] = {}
        for line_no, raw in enumerate(reader, start=2):
            if None in raw:
                raise ValueError(
                    f"Initial cell manifest line {line_no} has more fields than its header: {source}"
                )
            row = {key: str(raw.get(key, "") or "").strip() for key in INITIAL_CELL_MANIFEST_FIELDS}
            cell_id = row.get("sample_id", "")
            if not cell_id:
                if any(row.values()):
                    raise ValueError(
                        f"Initial cell manifest line {line_no} has values but no sample_id: {source}"
                    )
                continue
            if SAMPLE_ID_RE.fullmatch(cell_id) is None:
                raise ValueError(
                    f"Initial cell manifest line {line_no} has invalid sample_id={cell_id!r}: {source}"
                )
            if SAMPLE_ID_RE.fullmatch(row.get("project_sample_id", "")) is None:
                raise ValueError(
                    f"Initial cell manifest line {line_no} has invalid project_sample_id="
                    f"{row.get('project_sample_id')!r}: {source}"
                )
            if cell_id != cell_identity(row["project_sample_id"], row["protocol"], row["plate_id"], row["dna_barcode"]):
                raise ValueError(
                    f"Initial cell manifest line {line_no} sample_id does not match "
                    f"project_sample_id/plate_id: {source}"
                )
            if not row["species"] or row["protocol"] not in _PROTOCOLS:
                raise ValueError(
                    f"Initial cell manifest line {line_no} has invalid species/protocol: {source}"
                )
            if not row["cell_order"] or not row["dna_barcode"]:
                raise ValueError(
                    f"Initial cell manifest line {line_no} has incomplete cell identity: {source}"
                )
            if not any(row[column] for column in get_protocol(row["protocol"]).raw_input_columns):
                raise ValueError(
                    f"Initial cell manifest line {line_no} has no raw sample route: {source}"
                )
            if row["protocol"] == "srd" and not row["rna_barcode"]:
                raise ValueError(
                    f"Initial cell manifest line {line_no} has no SRD RNA barcode: {source}"
                )
            if cell_id in rows:
                raise ValueError(f"Initial cell manifest repeats sample_id={cell_id!r} at line {line_no}")
            rows[cell_id] = row
    if not rows:
        raise ValueError(f"Initial cell manifest contains no cells: {source}")
    return rows


def validate_demux_report_against_cell_manifest(
    report: DemuxReport,
    cell_manifest: str,
    project_sample_id: str,
) -> None:
    """要求 Rust demux 输出身份与 initial cell 快照一致。"""
    expected_all = read_cell_identity_manifest(cell_manifest)
    expected = {
        cell_id: row
        for cell_id, row in expected_all.items()
        if row.get("project_sample_id") == project_sample_id
    }
    if not expected:
        raise ValueError(
            f"Initial cell manifest has no rows for project_sample_id={project_sample_id!r}"
        )
    seen: set[str] = set()
    for sample in report["samples"]:
        reported_project = sample["sample_name"]
        plate = sample["plate_id"]
        cell_id = sample["cell_sample_id"]
        if reported_project != project_sample_id:
            raise ValueError(
                f"Demux report project identity mismatch: {reported_project!r} != {project_sample_id!r}"
            )
        row = expected.get(cell_id)
        if row is None:
            raise ValueError(f"Demux report emitted unknown canonical Sample ID: {cell_id!r}")
        checks = {
            "plate_id": plate,
            "cell_order": sample["cell_order"],
            "dna_barcode": sample["dna_barcode"],
            "rna_barcode": sample["rna_barcode"],
        }
        for key, observed in checks.items():
            expected_value = str(row.get(key, "") or "").strip()
            if key == "rna_barcode" and not observed:
                continue
            if observed != expected_value:
                raise ValueError(
                    f"Demux identity mismatch for {cell_id!r}: {key}={observed!r}, "
                    f"initial manifest={expected_value!r}"
                )
        if cell_id in seen:
            raise ValueError(f"Demux report repeats canonical Sample ID: {cell_id!r}")
        seen.add(cell_id)
    missing = sorted(set(expected) - seen)
    extra = sorted(seen - set(expected))
    if missing or extra:
        detail = []
        if missing:
            detail.append("missing=" + ",".join(missing[:20]))
        if extra:
            detail.append("extra=" + ",".join(extra[:20]))
        raise ValueError("Demux report / initial cell manifest mismatch: " + "; ".join(detail))


def manifest_raw_samples(manifest_rows: list[ProjectSampleRow], protocol: str) -> set[str]:
    """提取 manifest 登记的本地物理输入；Cabernet/TAPS 的外部 RNA 关联不参与选样。"""
    return {
        row[column] for row in manifest_rows
        for column in get_protocol(protocol).raw_input_columns if row.get(column)
    }


def map_discovered_samples_to_manifest(
    discovered: DiscoveredSamples,
    manifest_rows: list[ProjectSampleRow],
    *,
    protocol: str,
) -> MappedDemuxRuns:
    """仅将 manifest 登记的 raw FASTQ 映射到 canonical 样本与 route；登记输入缺失时拒绝。"""
    protocol = str(protocol).strip().lower()
    input_columns = get_protocol(protocol).raw_input_columns
    expected_raw: dict[str, dict[str, str]] = {}
    for row in manifest_rows:
        sid = row["sample_id"]
        dna_raw = row.get("dna_raw_sample", "")
        rna = row.get("rna_sample", "")
        if dna_raw:
            route_id = default_dna_route(protocol)
            expected_raw[dna_raw] = {
                "sample_id": sid,
                "route_id": route_id,
                "source_column": "dna_raw_sample",
            }
        if rna and "rna_sample" in input_columns:
            expected_raw[rna] = {
                "sample_id": sid,
                "route_id": "srd_rna_enrichment_tube",
                "source_column": "rna_sample",
            }

    discovered_set = set(discovered)
    expected_set = set(expected_raw)
    missing_fastq = sorted(expected_set - discovered_set)
    if missing_fastq:
        raise ValueError("Raw FASTQ / sample_manifest mismatch: manifest raw sample(s) without paired FASTQ: " + ", ".join(missing_fastq))

    mapped: MappedDemuxRuns = {}
    for raw_sample in sorted(expected_raw):
        identity = expected_raw[raw_sample]
        sid = identity["sample_id"]
        route_id = identity["route_id"]
        run_id = sid if identity["source_column"] == "dna_raw_sample" else f"{sid}__rna"
        if run_id in mapped:
            raise ValueError(f"Canonical demux run ID collision: {run_id!r}")
        mapped[run_id] = {
            "r1": discovered[raw_sample]["r1"],
            "r2": discovered[raw_sample]["r2"],
            "raw_sample": raw_sample,
            "sample_id": sid,
            "route_id": route_id,
        }
    return mapped


DEMUX_MANIFEST_FIELDS = (
    "demux_sample", "project_sample_id", "raw_sample", "pipeline_mode", "route_id", "demux_mode",
    "downstream_dna", "demux_report_schema_version", "demux_build_version", "demux_source_revision",
    "retention_policy", "retention_threshold",
    "rna_output_role", "high_cph_role", "plate_id", "cell_order", "dna_barcode", "rna_barcode",
    "dna_reads", "rna_reads", "input_fastq_pairs", "input_read_pairs", "usable_read_pairs",
    "matched_reads", "unmatched_reads", "ambiguous_reads", "short_reads", "unassigned_rate",
    "dna_r1", "dna_r2", "rna_r1", "rna_r2", "dna_status", "rna_status",
)


DEMUX_CELL_STATUSES = {"Pass", "Zero_Output", "Low_Reads_Removed"}


def _validate_cell_status(
    modality: str,
    reads: int,
    status: str,
    *,
    threshold: int,
    line_no: int,
    path: Path,
) -> None:
    if status not in DEMUX_CELL_STATUSES:
        raise ValueError(
            f"Demux manifest line {line_no} has invalid {modality}_status={status!r}: {path}"
        )
    expected = cell_status(reads, threshold)
    if status != expected:
        raise ValueError(
            f"Demux manifest line {line_no} has {modality}_reads={reads} but "
            f"{modality}_status={status!r}; expected {expected!r}: {path}"
        )


def write_demux_manifest(
    report: DemuxReport,
    out_tsv: str,
    dna_out_dir: str,
    rna_out_dir: str,
    pipeline_mode: str,
    route_id: str,
    downstream_dna: bool,
    raw_sample: str,
    cell_manifest_path: str,
) -> None:
    """把 Rust report 转换为规范的 checkpoint manifest。"""
    rep = report
    report_schema = rep["schema_version"]
    report_project_ids = {sample["sample_name"] for sample in rep["samples"]}
    project_sample_id_for_report = next(iter(report_project_ids))
    validate_demux_report_against_cell_manifest(
        rep, cell_manifest_path, project_sample_id_for_report
    )

    demux_mode = rep["mode"]
    policy = get_protocol(pipeline_mode)
    route = get_raw_route(route_id, policy.name)
    if demux_mode != route.demux_mode:
        raise ValueError(
            f"Demux report mode {demux_mode!r} does not match raw route "
            f"{route.route_id!r} ({route.demux_mode!r})"
        )
    if not isinstance(downstream_dna, bool):
        raise ValueError("downstream_dna must be a boolean")
    if downstream_dna != route.downstream_dna:
        raise ValueError(
            f"Raw route downstream-DNA mismatch for {route.route_id!r}: "
            f"requested={downstream_dna}, expected={route.downstream_dna}"
        )
    if SAMPLE_ID_RE.fullmatch(str(raw_sample)) is None:
        raise ValueError(f"raw_sample has an invalid identity: {raw_sample!r}")
    input_fastq_pairs = rep["input_fastq_pairs"]
    input_pairs = rep["input_read_pairs"]
    usable_pairs = rep["usable_read_pairs"]
    matched_reads = rep["matched_reads"]
    fate = rep["read_fates"]
    unmatched_reads = fate["unmatched"]
    ambiguous_reads = fate["ambiguous"]
    short_reads = fate["short"]
    unassigned_reads = max(input_pairs - matched_reads, 0)
    unmatched_rate = (unassigned_reads / input_pairs) if input_pairs else 0.0
    rna_enabled = demux_mode == "dna-rna"

    rows: list[dict[str, object]] = []
    for sample in rep["samples"]:
        project_sample_id = sample["sample_name"]
        plate_id = sample["plate_id"]
        cell_sample_id = sample["cell_sample_id"]
        dna_bc = sample["dna_barcode"]
        rna_bc = sample["rna_barcode"]
        dna_dir = Path(dna_out_dir) / cell_sample_id
        if rna_enabled and rna_bc:
            rna_dir = Path(rna_out_dir) / cell_sample_id
            rna_r1 = str(rna_dir / f"{cell_sample_id}_R1.fastq.gz")
            rna_r2 = str(rna_dir / f"{cell_sample_id}_R2.fastq.gz")
            rna_status = sample["rna_status"]
        else:
            rna_r1 = rna_r2 = rna_status = ""

        rows.append(
            {
                "demux_sample": cell_sample_id,
                "project_sample_id": project_sample_id,
                "raw_sample": raw_sample,
                "pipeline_mode": policy.name,
                "route_id": route.route_id,
                "demux_mode": demux_mode,
                "downstream_dna": "true" if route.downstream_dna else "false",
                "demux_report_schema_version": report_schema,
                "demux_build_version": rep["build_version"],
                "demux_source_revision": rep["source_revision"],
                "retention_policy": rep["retention_policy"],
                "retention_threshold": rep["retention_threshold"],
                "rna_output_role": route.rna_output_role,
                "high_cph_role": route.high_cph_role,
                "plate_id": plate_id,
                "cell_order": sample["cell_order"],
                "dna_barcode": dna_bc,
                "rna_barcode": rna_bc if rna_enabled else "",
                "dna_reads": sample["dna_read_count"],
                "rna_reads": sample["rna_read_count"] if rna_enabled else 0,
                "input_fastq_pairs": input_fastq_pairs,
                "input_read_pairs": input_pairs,
                "usable_read_pairs": usable_pairs,
                "matched_reads": matched_reads,
                "unmatched_reads": unmatched_reads,
                "ambiguous_reads": ambiguous_reads,
                "short_reads": short_reads,
                "unassigned_rate": f"{unmatched_rate:.6f}",
                "dna_r1": str(dna_dir / f"{cell_sample_id}_R1.fastq.gz"),
                "dna_r2": str(dna_dir / f"{cell_sample_id}_R2.fastq.gz"),
                "rna_r1": rna_r1,
                "rna_r2": rna_r2,
                "dna_status": sample["dna_status"],
                "rna_status": rna_status,
            }
        )

    summary = {field: "" for field in DEMUX_MANIFEST_FIELDS}
    summary.update(
        {
            "demux_sample": "__SUMMARY__",
            "raw_sample": raw_sample,
            "pipeline_mode": policy.name,
            "route_id": route.route_id,
            "demux_mode": demux_mode,
            "downstream_dna": "true" if route.downstream_dna else "false",
            "demux_report_schema_version": report_schema,
            "demux_build_version": rep["build_version"],
            "demux_source_revision": rep["source_revision"],
            "retention_policy": rep["retention_policy"],
            "retention_threshold": rep["retention_threshold"],
            "rna_output_role": route.rna_output_role,
            "high_cph_role": route.high_cph_role,
            "input_fastq_pairs": input_fastq_pairs,
            "input_read_pairs": input_pairs,
            "usable_read_pairs": usable_pairs,
            "matched_reads": matched_reads,
            "unmatched_reads": unmatched_reads,
            "ambiguous_reads": ambiguous_reads,
            "short_reads": short_reads,
            "unassigned_rate": f"{unmatched_rate:.6f}",
        }
    )
    rows.append(summary)
    write_tsv_atomic(out_tsv, rows, DEMUX_MANIFEST_FIELDS)
    read_demux_manifest(out_tsv)


def read_demux_manifest(manifest_path: str) -> list[DemuxRow]:
    """完整读取 demux manifest；任何不一致都以 fail closed 报错收场。"""
    path = Path(manifest_path)
    if not path.is_file():
        raise FileNotFoundError(f"Demux manifest is missing: {path}")
    text = path.read_text(encoding="utf-8")
    if not text:
        raise ValueError(f"Demux manifest is empty: {path}")
    if not text.endswith("\n"):
        raise ValueError(f"Demux manifest has no terminal newline and may be truncated: {path}")
    lines = text.splitlines()
    if not lines:
        raise ValueError(f"Demux manifest is empty: {path}")
    header = lines[0].split("\t")
    if header != list(DEMUX_MANIFEST_FIELDS):
        raise ValueError(
            f"Demux manifest header/order does not match the Alopex contract: {path}; "
            f"expected={list(DEMUX_MANIFEST_FIELDS)!r}, observed={header!r}. "
            "Remove the stale/corrupt manifest and rerun demux."
        )

    rows: list[DemuxRow] = []
    row_lines: list[int] = []
    for line_number, line in enumerate(lines[1:], start=2):
        if not line.strip():
            raise ValueError(f"Demux manifest contains a blank row at line {line_number}: {path}")
        values = line.split("\t")
        if len(values) != len(header):
            raise ValueError(
                f"Demux manifest row {line_number} in {path} has {len(values)} fields; "
                f"expected {len(header)}. Remove the corrupt manifest and rerun demux."
            )
        row = dict(zip(header, values, strict=True))
        rows.append(row)
        row_lines.append(line_number)

    summaries = [idx for idx, row in enumerate(rows) if row["demux_sample"] == "__SUMMARY__"]
    if len(summaries) != 1 or summaries[0] != len(rows) - 1:
        raise ValueError(
            f"Demux manifest must contain exactly one terminal __SUMMARY__ row: {path}"
        )
    cells = rows[:-1]
    if not cells:
        raise ValueError(f"Demux manifest contains no cell rows: {path}")
    summary = rows[-1]

    def nonnegative_int(row: Mapping[str, str], field: str, line_no: int) -> int:
        raw = str(row.get(field, ""))
        if not re.fullmatch(r"[0-9]+", raw):
            raise ValueError(
                f"Demux manifest line {line_no} has invalid non-negative integer "
                f"{field}={raw!r}: {path}"
            )
        return int(raw)

    def unit_float(row: Mapping[str, str], field: str, line_no: int) -> float:
        raw = str(row.get(field, ""))
        try:
            value = float(raw)
        except ValueError as exc:
            raise ValueError(
                f"Demux manifest line {line_no} has invalid numeric {field}={raw!r}: {path}"
            ) from exc
        if not 0.0 <= value <= 1.0:
            raise ValueError(
                f"Demux manifest line {line_no} has out-of-range {field}={raw!r}: {path}"
            )
        return value

    summary_line = row_lines[-1]
    if summary["demux_report_schema_version"] != "3":
        raise ValueError(
            f"Demux manifest {path} uses report schema "
            f"{summary['demux_report_schema_version']!r}; Alopex requires schema 3"
        )
    for field in (
        "project_sample_id", "plate_id", "cell_order", "dna_barcode", "rna_barcode",
        "dna_reads", "rna_reads", "dna_r1", "dna_r2", "rna_r1", "rna_r2",
        "dna_status", "rna_status",
    ):
        if summary[field] != "":
            raise ValueError(
                f"Demux manifest summary field {field} must be empty, found {summary[field]!r}: {path}"
            )

    repeated_fields = (
        "raw_sample", "pipeline_mode", "route_id", "demux_mode", "downstream_dna",
        "demux_report_schema_version", "demux_build_version", "demux_source_revision",
        "retention_policy", "retention_threshold",
        "rna_output_role", "high_cph_role", "input_fastq_pairs", "input_read_pairs",
        "usable_read_pairs", "matched_reads", "unmatched_reads", "ambiguous_reads",
        "short_reads", "unassigned_rate",
    )
    for idx, cell in enumerate(cells):
        for field in repeated_fields:
            if cell[field] != summary[field]:
                raise ValueError(
                    f"Demux manifest line {row_lines[idx]} field {field} differs from summary: {path}"
                )

    if summary["downstream_dna"] not in {"true", "false"}:
        raise ValueError(f"Invalid downstream_dna value in demux manifest: {path}")
    if summary["retention_policy"] != "matched_read_pairs":
        raise ValueError(
            f"Unsupported demux retention_policy={summary['retention_policy']!r}: {path}"
        )

    try:
        route = get_raw_route(summary["route_id"], summary["pipeline_mode"])
    except Exception as exc:
        raise ValueError(f"Invalid protocol/raw route in demux manifest {path}") from exc
    observed_route = (
        summary["demux_mode"], summary["downstream_dna"], summary["rna_output_role"],
        summary["high_cph_role"],
    )
    expected_route = (
        route.demux_mode, "true" if route.downstream_dna else "false",
        route.rna_output_role, route.high_cph_role,
    )
    if observed_route != expected_route:
        raise ValueError(
            f"Demux manifest route metadata is inconsistent: observed={observed_route!r}, "
            f"expected={expected_route!r}: {path}"
        )

    threshold = nonnegative_int(summary, "retention_threshold", summary_line)
    totals = {
        field: nonnegative_int(summary, field, summary_line)
        for field in (
            "input_fastq_pairs", "input_read_pairs", "usable_read_pairs", "matched_reads",
            "unmatched_reads", "ambiguous_reads", "short_reads",
        )
    }
    if totals["input_fastq_pairs"] < 1:
        raise ValueError(f"Demux manifest input_fastq_pairs must be >= 1: {path}")
    if totals["usable_read_pairs"] != totals["input_read_pairs"] - totals["short_reads"]:
        raise ValueError(f"Demux manifest usable/input/short read accounting is inconsistent: {path}")
    if (
        totals["matched_reads"] + totals["unmatched_reads"] + totals["ambiguous_reads"]
        != totals["usable_read_pairs"]
    ):
        raise ValueError(f"Demux manifest terminal read-fate accounting is inconsistent: {path}")
    expected_unassigned = (
        (totals["input_read_pairs"] - totals["matched_reads"]) / totals["input_read_pairs"]
        if totals["input_read_pairs"] else 0.0
    )
    if abs(unit_float(summary, "unassigned_rate", summary_line) - expected_unassigned) > 5.1e-7:
        raise ValueError(f"Demux manifest unassigned_rate is inconsistent with read counts: {path}")

    seen_cells: set[str] = set()
    project_ids: set[str] = set()
    assigned_reads = 0
    for idx, cell in enumerate(cells):
        line_no = row_lines[idx]
        cell_id = cell["demux_sample"]
        if SAMPLE_ID_RE.fullmatch(cell_id) is None or cell_id == "__SUMMARY__":
            raise ValueError(f"Demux manifest line {line_no} has invalid cell ID {cell_id!r}: {path}")
        if cell_id in seen_cells:
            raise ValueError(f"Demux manifest repeats cell ID {cell_id!r}: {path}")
        seen_cells.add(cell_id)
        project_id = cell["project_sample_id"]
        if SAMPLE_ID_RE.fullmatch(project_id) is None:
            raise ValueError(
                f"Demux manifest line {line_no} has invalid project_sample_id={project_id!r}: {path}"
            )
        project_ids.add(project_id)
        if (not cell["plate_id"] and route.protocol != "droplet") or not cell["dna_barcode"] or not cell["dna_r1"] or not cell["dna_r2"]:
            raise ValueError(f"Demux manifest line {line_no} has incomplete cell identity/paths: {path}")
        if not cell["cell_order"].strip():
            raise ValueError(f"Demux manifest line {line_no} has empty cell_order: {path}")
        dna_reads = nonnegative_int(cell, "dna_reads", line_no)
        rna_reads = nonnegative_int(cell, "rna_reads", line_no)
        assigned_reads += dna_reads + rna_reads
        unit_float(cell, "unassigned_rate", line_no)

        _validate_cell_status(
            "dna", dna_reads, cell["dna_status"],
            threshold=threshold, line_no=line_no, path=path,
        )
        if route.demux_mode == "dna-rna":
            if not cell["rna_barcode"] or not cell["rna_r1"] or not cell["rna_r2"]:
                raise ValueError(f"DNA/RNA demux row has incomplete RNA identity/paths at line {line_no}: {path}")
            _validate_cell_status(
                "rna", rna_reads, cell["rna_status"],
                threshold=threshold, line_no=line_no, path=path,
            )
        elif cell["rna_barcode"] or cell["rna_r1"] or cell["rna_r2"] or cell["rna_status"] or rna_reads:
            raise ValueError(f"DNA-only demux row contains RNA output fields at line {line_no}: {path}")

    if len(project_ids) != 1:
        raise ValueError(f"Demux manifest must describe exactly one project sample: {path}")
    if assigned_reads != totals["matched_reads"]:
        raise ValueError(
            f"Demux manifest cell read sum {assigned_reads} != matched_reads "
            f"{totals['matched_reads']}: {path}"
        )
    return rows


def write_run_metadata(
    report_path: str,
    mqc_path: str,
    manifest_path: str,
    dna_dir: str,
    rna_dir: str,
    pipeline_mode: str,
    route_id: str,
    downstream_dna: str,
    raw_sample: str,
    cell_manifest_path: str,
    *, calling_metrics: str = "", structure_metrics: str = "",
) -> None:
    """把 Rust report 一次性转换为 MQC JSON 与 checkpoint manifest。"""
    report = read_demux_report(report_path)
    write_demux_mqc(report, mqc_path, pipeline_mode, route_id, calling_metrics=calling_metrics, structure_metrics=structure_metrics)
    write_demux_manifest(
        report,
        manifest_path,
        dna_dir,
        rna_dir,
        pipeline_mode,
        route_id,
        str(downstream_dna).strip().lower() == "true",
        raw_sample,
        cell_manifest_path,
    )


def collect_rna_fastq(
    manifest_path: str, intermediate_dir: str, marker: str, run_id: str
) -> int:
    """把已提交 demux run 的 RNA FASTQ 硬链接收集到 02_work/RNA/FASTQ 保留区。

    demux scratch（02_work/demux/）在成功发布后会被整树回收；本步骤在 DAG 内
    先行建立同 inode 硬链接，保证 SRD 的 RNA 原始数据不受该清理影响。
    """
    rows = read_demux_manifest(manifest_path)
    collected = 0
    for row in rows:
        if row["demux_sample"] == "__SUMMARY__":
            continue
        if str(row.get("rna_status", "")) != "Pass":
            continue
        route_id = str(row.get("route_id", "")).strip()
        cell_id = str(row.get("demux_sample", "")).strip()
        if not route_id or not cell_id:
            raise ValueError(
                f"Demux manifest row with rna_status=Pass lacks route/sample id: {run_id}"
            )
        for mate in (1, 2):
            source = row.get(f"rna_r{mate}", "")
            if not str(source).strip():
                raise ValueError(
                    f"Demux manifest marks {cell_id!r} RNA as Pass without rna_r{mate}: {run_id}"
                )
            destination = rna_fastq_collect_destination(Path(intermediate_dir), route_id, cell_id, mate)
            destination.parent.mkdir(parents=True, exist_ok=True)
            temporary = project_temporary_parent(destination) / (destination.name + f".tmp.{os.getpid()}.{uuid.uuid4().hex}")
            try:
                os.link(source, temporary)
                os.replace(temporary, destination)
            finally:
                temporary.unlink(missing_ok=True)
            collected += 1

    payload = {"run_id": run_id, "collected_files": collected}
    marker_path = Path(marker)
    atomic_write_json(marker_path, payload)
    print(f"INFO: run={run_id} collected_rna_files={collected}")
    return collected


DEMUX_COMPLETION_SCHEMA = 1


_CONTENT_VERIFIED_KEYS = ("manifest", "initial_cell_manifest")


_EXISTENCE_KEYS = ("report", "multiqc")


def _cell_set_digest(cell_ids: list[str]) -> str:
    payload = "".join(f"{cell_id}\n" for cell_id in sorted(cell_ids)).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def validate_demux_manifest_against_cell_manifest(
    manifest_path: str | Path,
    cell_manifest_path: str | Path,
    project_sample_id: str,
    *,
    rows: list[DemuxRow] | None = None,
) -> None:
    """要求 demux manifest 恰好包含 initial cell manifest 中每个 planned cell 各一次。"""
    expected_all = read_cell_identity_manifest(cell_manifest_path)
    expected = {
        cell_id: row
        for cell_id, row in expected_all.items()
        if row.get("project_sample_id") == project_sample_id
    }
    if not expected:
        raise ValueError(
            f"Initial cell manifest has no cells for project_sample_id={project_sample_id!r}"
        )
    observed_rows = [
        row for row in (rows if rows is not None else read_demux_manifest(str(manifest_path)))
        if row["demux_sample"] != "__SUMMARY__"
    ]
    observed = {row["demux_sample"]: row for row in observed_rows}
    if set(observed) != set(expected):
        missing = sorted(set(expected) - set(observed))
        extra = sorted(set(observed) - set(expected))
        detail = []
        if missing:
            detail.append("missing=" + ",".join(missing[:20]))
        if extra:
            detail.append("extra=" + ",".join(extra[:20]))
        raise ValueError(
            "Demux manifest / initial cell manifest cell-set mismatch: " + "; ".join(detail)
        )
    for cell_id, demux in observed.items():
        planned = expected[cell_id]
        for field in ("project_sample_id", "plate_id", "cell_order", "dna_barcode"):
            if str(demux.get(field, "")) != str(planned.get(field, "")):
                raise ValueError(
                    f"Demux identity mismatch for {cell_id!r}: {field}="
                    f"{demux.get(field)!r}, initial manifest={planned.get(field)!r}"
                )
        if demux.get("rna_barcode") and demux["rna_barcode"] != planned.get("rna_barcode", ""):
            raise ValueError(
                f"Demux identity mismatch for {cell_id!r}: rna_barcode="
                f"{demux['rna_barcode']!r}, initial manifest={planned.get('rna_barcode')!r}"
            )


def _validate_committed_fastqs(
    rows: list[DemuxRow], dna_dir: Path, rna_dir: Path
) -> None:
    dna_root = dna_dir.resolve()
    rna_root = rna_dir.resolve()
    if not dna_root.is_dir() or not rna_root.is_dir():
        raise FileNotFoundError(
            f"Committed demux output directory is missing: DNA={dna_root}, RNA={rna_root}"
        )
    for row in rows:
        if row["demux_sample"] == "__SUMMARY__":
            continue
        for modality, root in (("dna", dna_root), ("rna", rna_root)):
            if row[f"{modality}_status"] != "Pass":
                continue
            for mate in (1, 2):
                fastq = Path(row[f"{modality}_r{mate}"]).resolve()
                try:
                    fastq.relative_to(root)
                except ValueError as exc:
                    raise ValueError(
                        f"Committed {modality.upper()} FASTQ escapes its run directory: {fastq}"
                    ) from exc
                if not fastq.is_file():
                    raise FileNotFoundError(
                        f"Demux manifest marks {row['demux_sample']!r} {modality.upper()} as Pass "
                        f"but FASTQ is missing: {fastq}"
                    )


def _load_committed_completion(
    completion_path: str | Path,
) -> tuple[dict[str, Any], dict[str, Path], str, str]:
    completion_file = Path(completion_path).resolve()
    if not completion_file.is_file():
        raise FileNotFoundError(f"Demux completion manifest is missing: {completion_file}")
    try:
        value = json.loads(completion_file.read_text(encoding="utf-8"))
    except Exception as exc:
        raise ValueError(f"Invalid demux completion manifest {completion_file}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"Demux completion manifest root must be an object: {completion_file}")
    if value.get("schema_version") != DEMUX_COMPLETION_SCHEMA or value.get("status") != "committed":
        raise ValueError(f"Demux generation is not a committed schema-1 generation: {completion_file}")
    generation = str(value.get("generation_id", ""))
    if re.fullmatch(r"[0-9a-f]{32}", generation) is None:
        raise ValueError(f"Demux completion generation_id is invalid: {completion_file}")
    run_id = str(value.get("run_id", ""))
    project_sample_id = str(value.get("project_sample_id", ""))
    if SAMPLE_ID_RE.fullmatch(run_id) is None or SAMPLE_ID_RE.fullmatch(project_sample_id) is None:
        raise ValueError(f"Demux completion has an invalid run/project identity: {completion_file}")

    identities: dict[str, Path] = {}
    for key in (*_CONTENT_VERIFIED_KEYS, *_EXISTENCE_KEYS):
        recorded = value.get(key)
        if not isinstance(recorded, dict):
            raise ValueError(f"Demux completion is missing file identity {key!r}: {completion_file}")
        path = Path(str(recorded.get("path", ""))).resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Committed demux file is missing: {path}")
        identities[key] = path
    for key in ("dna_dir", "rna_dir"):
        directory = Path(str(value.get(key, ""))).resolve()
        if not directory.is_dir():
            raise FileNotFoundError(f"Committed demux output directory is missing: {directory}")
    return value, identities, run_id, project_sample_id


def write_demux_completion(
    out_json: str,
    *,
    run_id: str,
    project_sample_id: str,
    manifest_path: str,
    report_path: str,
    mqc_path: str,
    dna_dir: str,
    rna_dir: str,
    cell_manifest_path: str,
) -> None:
    """为一个完整发布的 demux generation 写出 commit 边界记录。"""
    rows = read_demux_manifest(manifest_path)
    validate_demux_manifest_against_cell_manifest(
        manifest_path, cell_manifest_path, project_sample_id, rows=rows
    )
    _validate_committed_fastqs(rows, Path(dna_dir), Path(rna_dir))
    cell_ids = [row["demux_sample"] for row in rows if row["demux_sample"] != "__SUMMARY__"]
    completion: DemuxCompletion = {
        "schema_version": DEMUX_COMPLETION_SCHEMA,
        "status": "committed",
        "generation_id": uuid.uuid4().hex,
        "run_id": str(run_id),
        "project_sample_id": str(project_sample_id),
        "cell_count": len(cell_ids),
        "cell_ids_sha256": _cell_set_digest(cell_ids),
        "manifest": file_identity(manifest_path),
        "report": file_identity(report_path),
        "multiqc": file_identity(mqc_path),
        "initial_cell_manifest": file_identity(cell_manifest_path),
        "dna_dir": str(Path(dna_dir).resolve()),
        "rna_dir": str(Path(rna_dir).resolve()),
    }
    output = Path(out_json)
    atomic_write_json(output, completion)


def validate_demux_completion(
    completion_path: str | Path,
    *,
    expected_manifest_path: str | Path | None = None,
    expected_cell_manifest_path: str | Path | None = None,
) -> dict[str, Any]:
    """校验一笔 demux generation 的 commit 记录及其全部语义输入。"""
    value, identities, _run_id, project_sample_id = _load_committed_completion(completion_path)
    completion_file = Path(completion_path).resolve()
    for key in _CONTENT_VERIFIED_KEYS:
        observed = file_identity(identities[key])
        if observed != value[key]:
            raise ValueError(
                f"Demux completion file identity changed for {key}: {identities[key]}; "
                "rerun the demux checkpoint"
            )
    if expected_manifest_path is not None and identities["manifest"] != Path(expected_manifest_path).resolve():
        raise ValueError(f"Demux completion points to the wrong manifest: {completion_file}")
    if (
        expected_cell_manifest_path is not None
        and identities["initial_cell_manifest"] != Path(expected_cell_manifest_path).resolve()
    ):
        raise ValueError(f"Demux completion points to the wrong initial cell manifest: {completion_file}")

    rows = read_demux_manifest(str(identities["manifest"]))
    validate_demux_manifest_against_cell_manifest(
        identities["manifest"], identities["initial_cell_manifest"], project_sample_id, rows=rows
    )
    cell_ids = [row["demux_sample"] for row in rows if row["demux_sample"] != "__SUMMARY__"]
    if value.get("cell_count") != len(cell_ids) or value.get("cell_ids_sha256") != _cell_set_digest(cell_ids):
        raise ValueError(f"Demux completion cell set does not match its manifest: {completion_file}")
    _validate_committed_fastqs(
        rows,
        Path(str(value.get("dna_dir", ""))).resolve(),
        Path(str(value.get("rna_dir", ""))).resolve(),
    )
    value["_manifest_path"] = str(identities["manifest"])
    value["_cell_manifest_path"] = str(identities["initial_cell_manifest"])
    value["_rows"] = rows
    return value


def demux_manifest_index(
    state_dir: str | Path,
    run_ids: list[str] | tuple[str, ...],
    *, protocol: str,
) -> DemuxIndex:
    """校验已提交 generation 与当前协议，统一建立 run/cell 索引及通过名单。

    反复展开 DAG 的调用方必须自行负责 memoisation；本模块保持纯校验/索引层，不藏任何状态。
    """
    state_dir = Path(state_dir).resolve()
    protocol = get_protocol(protocol).name
    normalized_runs = tuple(sorted({str(run_id).strip() for run_id in run_ids}))
    if not normalized_runs or any(
        SAMPLE_ID_RE.fullmatch(run_id) is None for run_id in normalized_runs
    ):
        raise ValueError("run_ids must contain valid canonical demux run IDs")
    completions = [
        state_dir / run_id / "completion.json"
        for run_id in normalized_runs
    ]
    missing = [path for path in completions if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Current run is missing demux completion manifest(s): "
            + ", ".join(str(path) for path in missing)
        )

    by_cell: dict[str, tuple[Path, DemuxRow]] = {}
    passing: list[str] = []
    for completion in completions:
        run_id = completion.parent.name
        manifest = completion.parent / "manifest.tsv"
        committed = validate_demux_completion(
            completion, expected_manifest_path=manifest
        )
        if committed["run_id"] != run_id:
            raise ValueError(
                f"Demux completion directory/run_id mismatch: {completion}"
            )
        rows = tuple(row for row in committed["_rows"] if row["demux_sample"] != "__SUMMARY__")
        for row in rows:
            if row["pipeline_mode"] != protocol:
                raise ValueError(f"Demux manifest was produced under a different protocol: {manifest}")
            if row["downstream_dna"] != "true":
                continue
            sample = row["demux_sample"]
            if sample in by_cell:
                raise ValueError(f"Canonical DNA cell {sample!r} appears in more than one committed demux run")
            by_cell[sample] = (manifest, row)
            if row["dna_status"] == "Pass":
                passing.append(sample)
    return DemuxIndex(state_dir, normalized_runs, MappingProxyType(by_cell), tuple(passing))


def read_demux_generation(
    completion_path: str | Path,
    *,
    expected_manifest_path: str | Path,
    expected_cell_manifest_path: str | Path,
) -> dict[str, Any]:
    """按结构读取一份已提交代际，供聚合作业消费，不做二次完整性校验。

    完整性（文件 hash、cell-set digest、FASTQ 存在、manifest↔cell manifest
    语义）由 DAG 展开期的 ``validate_demux_completion`` 单点负责；本函数只保
    证结构可解析、声明的身份键一致、引用文件在读取时仍然存在。聚合作业自身
    消费的行级绑定（cell 归属、route/raw 绑定）仍由调用方执行。
    """
    value, identities, _run_id, _project_sample_id = _load_committed_completion(completion_path)
    completion_file = Path(completion_path).resolve()
    if identities["manifest"] != Path(expected_manifest_path).resolve():
        raise ValueError(f"Demux completion points to the wrong manifest: {completion_file}")
    if identities["initial_cell_manifest"] != Path(expected_cell_manifest_path).resolve():
        raise ValueError(f"Demux completion points to the wrong initial cell manifest: {completion_file}")

    rows = read_demux_manifest(str(identities["manifest"]))
    cell_ids = [row["demux_sample"] for row in rows if row["demux_sample"] != "__SUMMARY__"]
    if value.get("cell_count") != len(cell_ids) or value.get("cell_ids_sha256") != _cell_set_digest(cell_ids):
        raise ValueError(f"Demux completion cell set does not match its manifest: {completion_file}")

    value["_manifest_path"] = str(identities["manifest"])
    value["_cell_manifest_path"] = str(identities["initial_cell_manifest"])
    value["_rows"] = rows
    return value


def resolve_demux_fastq_from_manifest(
    demux_sample: str,
    read_pair: int,
    *,
    manifest_index: DemuxIndex,
) -> Path:
    """只从已提交的 manifest 解析 dna_status=Pass 的 DNA FASTQ。"""
    base = Path(str(demux_sample)).name
    key = f"dna_r{int(read_pair)}"
    match = manifest_index.by_cell.get(base)

    if match is None:
        raise FileNotFoundError(
            f"No current demux manifest row found for DNA sample {base!r}. "
            "The checkpoint manifest is the sole authority; remove stale demux markers "
            "and rerun demultiplexing instead of guessing a path."
        )
    manifest, row = match
    manifest = manifest.resolve()
    try:
        manifest.relative_to(manifest_index.state_dir)
    except ValueError as exc:
        raise ValueError(
            f"Manifest index entry escapes the committed demux state directory: {manifest}"
        ) from exc
    if row.get("dna_status") != "Pass":
        raise ValueError(
            f"Demux sample {base!r} is present in {manifest} but "
            f"dna_status={row.get('dna_status')!r}; it must not enter downstream DNA analysis."
        )
    value = str(row.get(key, "") or "").strip()
    if not value:
        raise ValueError(f"Manifest {manifest} has no {key} path for {base!r}")
    path = Path(value)
    if not path.is_file():
        raise FileNotFoundError(
            f"Manifest {manifest} points to missing FASTQ for {base!r}: {path}. "
            "Remove the stale/corrupt demux outputs and rerun the checkpoint."
        )
    return path


def rna_fastq_collect_destination(
    intermediate_dir: Path, route_id: str, demux_sample: str, mate: int
) -> Path:
    """SRD RNA FASTQ 在 02_work/RNA/FASTQ 保留区的规范落点（收集规则与聚合共用）。"""
    return (
        Path(intermediate_dir)
        / "RNA"
        / "FASTQ"
        / route_id
        / f"{demux_sample}_R{int(mate)}.fastq.gz"
    )


def write_demux_mqc(
    report: DemuxReport, out_json: str, pipeline_mode: str, route_id: str,
    *, calling_metrics: str = "", structure_metrics: str = "",
) -> None:
    """把 demultiplexing 统计转换为 MultiQC custom content。"""
    rep = report
    policy = get_protocol(pipeline_mode)
    route = get_raw_route(route_id, policy.name)
    demux_mode = rep["mode"]
    if demux_mode != route.demux_mode:
        raise ValueError(
            f"Demux report mode {demux_mode!r} does not match raw route "
            f"{route.route_id!r} ({route.demux_mode!r})"
        )
    data = {}
    for sample in rep["samples"]:
        key = sample["cell_sample_id"]
        row = {
            "Protocol": policy.name,
            "Route": route.route_id,
            "PlateID": sample["plate_id"],
            "Cell_Order": sample["cell_order"],
            "DNA_Barcode": sample["dna_barcode"],
            "DNA_Reads": sample["dna_read_count"],
            "DNA_Status": sample["dna_status"],
        }
        if demux_mode == "dna-rna":
            row.update({
                "RNA_Barcode": sample["rna_barcode"],
                "RNA_Reads": sample["rna_read_count"],
                "RNA_Status": sample["rna_status"],
            })
        data[key] = row

    fate = rep["read_fates"]
    summary = {
        "Input_Read_Pairs": rep["input_read_pairs"],
        "Usable_Read_Pairs": rep["usable_read_pairs"],
        "Matched_Read_Pairs": rep["matched_reads"],
        "Short_Read_Pairs": fate["short"],
        "Ambiguous_Read_Pairs": fate["ambiguous"],
        "Unmatched_Read_Pairs": fate["unmatched"],
    }
    if pipeline_mode == "droplet":
        calling = _read_json(Path(calling_metrics), "Droplet calling")
        structure = _read_json(Path(structure_metrics), "Droplet structure")
        if calling.get("status") != "complete" or calling["called_cells"] != len(rep["samples"]):
            raise ValueError("Droplet calling/report cell count mismatch")
        screen = calling["barcode_error_screen"]
        summary.update({"Called_Cells": calling["called_cells"], "Calling_Threshold": calling["threshold"],
                        "Design_Excluded_Candidates": calling["barcode_design_filter"]["excluded_candidates"],
                        "Design_Excluded_Exact_Read_Pairs": calling["barcode_design_filter"]["excluded_exact_read_pairs"],
                        "Barcode_Error_Candidates_Removed": screen["removed_barcodes"],
                        "Barcode_Molecule_Candidates_Checked": screen["molecule_checked"],
                        "Barcode_Insufficient_Molecule_Evidence": screen["insufficient_molecule_evidence"],
                        "Barcode_Child_Signatures_Capped": screen["child_signatures_capped"],
                        "Structured_Reads": structure["structured_reads"], "Structured_UMI_With_N": structure["umi_with_n"]})
    route_key = re.sub(r"[^A-Za-z0-9_.-]+", "_", route.route_id)
    obj = {
        "id": f"demux_rs_{route_key}",
        "section_name": f"Demultiplexing Stats: {route.route_id}",
        "description": (
            f"Protocol={policy.name}; route={route.route_id}; "
            f"read interpreter={demux_mode}; "
            f"retention threshold={rep['retention_threshold']} read pairs. "
            f"Run summary: {summary}"
        ),
        "plot_type": "table",
        "pconfig": {
            "id": f"demux_rs_table_{route_key}",
            "title": f"Demultiplexing Statistics: {route.route_id}",
            "col1_header": "Sample Name",
        },
        "data": data,
    }
    atomic_write_json(out_json, obj)
