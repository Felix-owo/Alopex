"""从 sealed delivery 生成单细胞甲基化 QC、统计缓存和可视化。

Notebook 经 resolve_notebook_paths / run_qc_processor / load_notebook_qc_frame
引导：解析项目与 Pipeline 根、以独立子进程运行 main、加载 QC 表并初始化
绘图主题；计算结果写入项目的 06_downstream/<delivery_id>，输入与公开科学
结果保持只读。
"""
from __future__ import annotations


import ast
import importlib.metadata
import csv
import hashlib
import json
import sys
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

_PIPELINE_ROOT = Path(__file__).resolve().parents[1]
_CORE_ROOT = _PIPELINE_ROOT / "core"
if str(_CORE_ROOT) not in sys.path:
    sys.path.insert(0, str(_CORE_ROOT))

from dna_pipeline import (
    CPG_REPRESENTATION,
    load_run_manifest,
    load_run_snapshot,
    sha256_file,
    snapshot_input_path,
)

FINAL_MANIFEST = Path("03_results/QC_Results/sample_manifest.tsv")

REQUIRED_MANIFEST_FIELDS = {
    "protocol",
    "sample_id",
    "project_sample_id",
    "plate_id",
    "dna_barcode",
    "rna_barcode",
    "dna_raw_sample",
    "dna_reads",
    "dna_status",
    "analysis_status",
    "backend_qc_status",
    "methylation_backend",
    "cutadapt_input_pairs",
    "trimmed_pairs",
    "both_primary_mapped_pairs",
    "one_primary_mapped_pairs",
    "backend_accepted_pairs",
    "backend_rejected_pairs",
    "unmapped_pairs",
    "ambiguous_pairs",
    "no_primary_pairs",
    "duplicate_pairs",
    "postdedup_pairs",
    "high_cph_assessed_pairs",
    "high_cph_flagged_pairs",
    "high_cph_removed_pairs",
    "final_retained_pairs",
    "mapping_policy",
    "native_mapping_policy",
    "native_mapping_unit",
    "native_mapping_pct",
    "dedup_policy",
    "pair_qc_metric_unit",
    "non_cpg_methylation_pct",
    "non_cpg_metric_source",
    "cpg_path",
    "cpg_representation",
}

PAIR_COUNT_FIELDS = (
    "cutadapt_input_pairs",
    "trimmed_pairs",
    "both_primary_mapped_pairs",
    "one_primary_mapped_pairs",
    "backend_accepted_pairs",
    "backend_rejected_pairs",
    "unmapped_pairs",
    "duplicate_pairs",
    "postdedup_pairs",
    "high_cph_assessed_pairs",
    "high_cph_flagged_pairs",
    "high_cph_removed_pairs",
    "final_retained_pairs",
)


@dataclass(frozen=True)
class SealedQCContext:
    """保存同一已发布交付中通过校验的只读输入。"""

    project_dir: Path
    run_manifest_path: Path
    run_manifest: Mapping[str, object]
    run_manifest_sha256: str
    delivery_id: str
    snapshot_path: Path
    snapshot: Mapping[str, object]
    snapshot_id: str
    snapshot_config_path: Path
    snapshot_config_sha256: str
    config: Mapping[str, object]
    pipeline_root: Path
    reference_path: Path
    reference_fai: Path
    cpg_representation: str
    methylation_backend: str
    output_inventory: Mapping[str, Mapping[str, object]]


def _read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [
            {str(key): str(value or "").strip() for key, value in row.items()}
            for row in reader
        ]
        return rows, list(reader.fieldnames or [])


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"Sealed delivery {label} must be a mapping")
    return value


def _safe_inventory(records: object) -> dict[str, Mapping[str, object]]:
    if not isinstance(records, list):
        raise ValueError("Sealed delivery output inventory is malformed")
    result: dict[str, Mapping[str, object]] = {}
    for raw in records:
        record = _mapping(raw, "output record")
        text = str(record.get("path", ""))
        relative = Path(text)
        if (
            not relative.parts
            or relative.is_absolute()
            or ".." in relative.parts
            or relative == Path(".")
        ):
            raise ValueError(f"Sealed delivery contains an unsafe output path: {text!r}")
        key = relative.as_posix()
        if key in result:
            raise ValueError(f"Sealed delivery repeats output path: {key}")
        result[key] = record
    return result


def _snapshot_reference_fai(
    snapshot: Mapping[str, object], reference_path: Path
) -> Path:
    basis = _mapping(snapshot.get("basis"), "snapshot basis")
    reference = _mapping(basis.get("reference"), "snapshot reference")
    recorded_reference = Path(str(reference.get("resolved_path", ""))).expanduser().resolve()
    if recorded_reference != reference_path:
        raise ValueError(
            "run_manifest reference differs from the sealed run snapshot: "
            f"{reference_path} != {recorded_reference}"
        )
    expected = Path(f"{reference_path}.fai").resolve()
    record = _mapping(reference.get("fai"), "snapshot reference FAI")
    if Path(str(record.get("path", ""))).resolve() != expected:
        raise ValueError(f"Sealed run snapshot does not identify exactly one FAI: {expected}")
    if not expected.is_file():
        raise FileNotFoundError(f"Sealed reference FAI is missing: {expected}")
    stat = expected.stat()
    if stat.st_size != record.get("size_bytes") or stat.st_mtime_ns != record.get("mtime_ns"):
        raise ValueError(f"Sealed reference FAI changed after the Pipeline run: {expected}")
    return expected


def _manifest_delivery_id(manifest: Mapping[str, object]) -> str:
    delivery = _mapping(manifest.get("delivery"), "delivery identity")
    delivery_id = str(delivery.get("id", ""))
    if len(delivery_id) != 64 or any(ch not in "0123456789abcdef" for ch in delivery_id):
        raise ValueError("Published Pipeline delivery_id is invalid")
    return delivery_id


def load_sealed_qc_context(
    project_dir: str | Path,
    *,
    expected_pipeline_root: str | Path | None = None,
) -> SealedQCContext:
    """加载完整交付及其不可变的 config/reference 契约。"""
    project = Path(project_dir).expanduser().resolve()
    manifest_path = project / "03_results" / "run_manifest.json"
    manifest = load_run_manifest(
        manifest_path,
        validate_outputs=True,
    )
    if manifest.get("status") != "complete":
        raise ValueError(f"Published Pipeline delivery is not complete: {manifest_path}")
    delivery_id = _manifest_delivery_id(manifest)

    run = _mapping(manifest.get("run"), "run metadata")
    snapshot_id = str(run.get("snapshot_id", ""))
    if snapshot_id != delivery_id:
        raise ValueError("Published delivery_id and run snapshot_id disagree")
    snapshot_path = Path(str(run.get("snapshot_path", ""))).expanduser().resolve()
    snapshot = load_run_snapshot(snapshot_path)
    if snapshot.get("snapshot_id") != snapshot_id:
        raise ValueError("run_manifest and run snapshot content address disagree")
    basis = _mapping(snapshot.get("basis"), "snapshot basis")
    snapshot_project = Path(
        str(_mapping(basis.get("project"), "snapshot project").get("path", ""))
    ).expanduser().resolve()
    if snapshot_project != project:
        raise ValueError(f"Sealed snapshot belongs to another project: {snapshot_project}")

    pipeline = _mapping(manifest.get("pipeline"), "pipeline metadata")
    if pipeline.get("name") != "Alopex":
        raise ValueError("run_manifest does not describe Alopex")
    pipeline_root = Path(str(pipeline.get("source_root", ""))).expanduser().resolve()
    if expected_pipeline_root is not None:
        expected = Path(expected_pipeline_root).expanduser().resolve()
        if pipeline_root != expected:
            raise ValueError(
                "Downstream processor does not come from the delivery-recorded Pipeline: "
                f"{expected} != {pipeline_root}"
            )

    input_copies = _mapping(snapshot.get("input_copies"), "snapshot input copies")
    config_copy = _mapping(input_copies.get("config"), "snapshot config copy")
    snapshot_config_path = snapshot_input_path(snapshot_path, "config")
    config = yaml.safe_load(snapshot_config_path.read_text(encoding="utf-8")) or {}
    if not isinstance(config, Mapping):
        raise ValueError(f"Sealed config root must be a mapping: {snapshot_config_path}")

    analysis = _mapping(config.get("analysis"), "snapshot config analysis")
    run_species = str(run.get("species", "")).strip()
    run_protocol = str(run.get("protocol", "")).strip().lower()
    config_species = str(config.get("species", "")).strip()
    config_protocol = str(analysis.get("protocol", "")).strip().lower()
    config_backend = str(analysis.get("methylation_backend", "")).strip().lower()
    if not run_species or config_species != run_species:
        raise ValueError("Sealed config species disagrees with the run manifest")
    if run_protocol not in {"cabernet", "srd", "taps", "droplet"} or config_protocol != run_protocol:
        raise ValueError("Sealed config protocol disagrees with the run manifest")

    reference_path = Path(str(run.get("reference_path", ""))).expanduser().resolve()
    if not reference_path.is_file():
        raise FileNotFoundError(f"Sealed Pipeline reference is missing: {reference_path}")
    reference_fai = _snapshot_reference_fai(snapshot, reference_path)

    run_config = _mapping(run, "run metadata")
    cpg_representation = str(run_config.get("cpg_representation", ""))
    if cpg_representation != CPG_REPRESENTATION:
        raise ValueError(
            "Unsupported sealed CpG representation: " f"{cpg_representation!r}"
        )
    methylation_backend = str(run_config.get("methylation_backend", "")).strip().lower()
    if methylation_backend not in {"biscuit", "bismark", "rastair"}:
        raise ValueError(f"Unsupported sealed methylation backend: {methylation_backend!r}")
    if (run_protocol == "taps") != (methylation_backend == "rastair"):
        raise ValueError("Sealed TAPS protocol requires Rastair")
    if config_backend != methylation_backend:
        raise ValueError("Sealed config backend disagrees with the run manifest")

    output_inventory = _safe_inventory(manifest.get("outputs"))
    if "QC_Results/sample_manifest.tsv" not in output_inventory:
        raise ValueError("Sealed delivery inventory omits QC_Results/sample_manifest.tsv")
    return SealedQCContext(
        project_dir=project,
        run_manifest_path=manifest_path,
        run_manifest=manifest,
        run_manifest_sha256=sha256_file(manifest_path),
        delivery_id=delivery_id,
        snapshot_path=snapshot_path,
        snapshot=snapshot,
        snapshot_id=snapshot_id,
        snapshot_config_path=snapshot_config_path,
        snapshot_config_sha256=str(config_copy.get("sha256", "")),
        config=config,
        pipeline_root=pipeline_root,
        reference_path=reference_path,
        reference_fai=reference_fai,
        cpg_representation=cpg_representation,
        methylation_backend=methylation_backend,
        output_inventory=output_inventory,
    )


def _inventory_record_for_project_file(
    context: SealedQCContext, value: str | Path, *, label: str
) -> tuple[Path, Mapping[str, object]]:
    relative = Path(str(value))
    if relative.is_absolute() or ".." in relative.parts or relative == Path("."):
        raise ValueError(f"Current Pipeline {label} has an unsafe path: {value!r}")
    path = (context.project_dir / relative).resolve()
    public_root = (context.project_dir / "03_results").resolve()
    try:
        inventory_key = path.relative_to(public_root).as_posix()
    except ValueError as exc:
        raise ValueError(
            f"Current Pipeline {label} escapes the sealed 03_results tree: {path}"
        ) from exc
    record = context.output_inventory.get(inventory_key)
    if record is None:
        raise ValueError(
            f"Current Pipeline {label} is not bound by the delivery inventory: {inventory_key}"
        )
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Current Pipeline {label} is missing or empty: {path}")
    return path, record


def _required_nonnegative_int(row: Mapping[str, str], key: str, sample_id: str) -> int:
    value = str(row.get(key, "") or "").strip()
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ValueError(f"Cell {sample_id!r} has invalid integer {key}={value!r}") from exc
    if parsed < 0:
        raise ValueError(f"Cell {sample_id!r} has negative {key}={parsed}")
    return parsed


def _required_percentage(
    row: Mapping[str, str], key: str, sample_id: str
) -> float:
    value = str(row.get(key, "") or "").strip()
    try:
        parsed = float(value)
    except ValueError as exc:
        raise ValueError(f"Cell {sample_id!r} has invalid percentage {key}={value!r}") from exc
    if parsed != parsed:
        return parsed
    if not 0 <= parsed <= 100:
        raise ValueError(f"Cell {sample_id!r} has out-of-range {key}={parsed}")
    return parsed


def read_active_cells(
    project_dir: str | Path,
    *,
    sealed_context: SealedQCContext | None = None,
) -> list[dict[str, Any]]:
    """按规范 read-pair 契约从任一 backend 返回已完成分析的细胞。"""
    context = sealed_context or load_sealed_qc_context(project_dir)
    project = Path(project_dir).expanduser().resolve()
    if project != context.project_dir:
        raise ValueError("project_dir and sealed_context refer to different projects")
    manifest = project / FINAL_MANIFEST
    _inventory_record_for_project_file(
        context, FINAL_MANIFEST, label="final sample manifest"
    )
    rows, fields = _read_tsv(manifest)
    missing = sorted(REQUIRED_MANIFEST_FIELDS - set(fields))
    if missing:
        raise ValueError(
            "Current sample_manifest.tsv is missing required field(s): "
            + ", ".join(missing)
        )

    active: list[dict[str, Any]] = []
    seen: set[str] = set()
    for row in rows:
        if not (
            row["dna_status"] == "Pass"
            and row["analysis_status"] == "PASS"
            and row["backend_qc_status"] == "PASS"
            and row["cpg_path"]
        ):
            continue
        sample_id = row["sample_id"]
        if row["protocol"] != context.config["analysis"]["protocol"]:
            raise ValueError(f"Cell {sample_id!r} protocol disagrees with sealed config")
        if not sample_id:
            raise ValueError("Active row in sample_manifest.tsv has an empty sample_id")
        if not row["project_sample_id"]:
            raise ValueError(
                f"Active row {sample_id!r} in sample_manifest.tsv has an empty project_sample_id"
            )
        if sample_id in seen:
            raise ValueError(f"Duplicate active sample_id in sample_manifest.tsv: {sample_id}")
        seen.add(sample_id)
        if row["methylation_backend"].strip().lower() != context.methylation_backend:
            raise ValueError(
                f"Cell {sample_id!r} backend disagrees with the sealed run manifest"
            )
        row_representation = str(row.get("cpg_representation", "") or "").strip()
        if row_representation != context.cpg_representation:
            raise ValueError(
                f"Cell {sample_id!r} CpG representation disagrees with the sealed run manifest"
            )

        cell: dict[str, Any] = dict(row)
        cpg_file, cpg_identity = _inventory_record_for_project_file(
            context, row["cpg_path"], label=f"CpG result for {sample_id}"
        )
        cell["cpg_file"] = cpg_file
        cell["cpg_output_identity"] = dict(cpg_identity)
        cell["demux_dna_reads_count"] = _required_nonnegative_int(
            row, "dna_reads", sample_id
        )
        for field in PAIR_COUNT_FIELDS:
            cell[field] = _required_nonnegative_int(row, field, sample_id)
        for field in ("ambiguous_pairs", "no_primary_pairs"):
            value = str(row.get(field, "") or "").strip()
            cell[field] = (
                None if value == "" else _required_nonnegative_int(row, field, sample_id)
            )
        cell["non_cpg_methylation_pct"] = _required_percentage(
            row, "non_cpg_methylation_pct", sample_id
        )
        if cell["pair_qc_metric_unit"] != "read_pairs":
            raise ValueError(f"Cell {sample_id!r} does not use read-pair QC units")
        backend = context.methylation_backend
        expected_policy = {
            "biscuit": "biscuit_any_primary_pair_funnel",
            "bismark": "bismark_unique_concordant",
            "rastair": "bwa_rastair_eligible_pairs",
        }[backend]
        if cell["mapping_policy"] != expected_policy:
            raise ValueError(f"Cell {sample_id!r} has an invalid mapping_policy")
        expected_native_policy = {
            "biscuit": "biscuit_primary_mapq_ge_40",
            "bismark": "bismark_unique_concordant_pairs",
            "rastair": "bwa_proper_primary_any_mate_mapq20",
        }[backend]
        expected_native_unit = {
            "biscuit": "individual_primary_reads_postdedup",
            "bismark": "read_pairs",
            "rastair": "read_pairs",
        }[backend]
        if cell["native_mapping_policy"] != expected_native_policy:
            raise ValueError(f"Cell {sample_id!r} has an invalid native_mapping_policy")
        if cell["native_mapping_unit"] != expected_native_unit:
            raise ValueError(f"Cell {sample_id!r} has an invalid native_mapping_unit")
        native_pct = float(row["native_mapping_pct"])
        if not 0.0 <= native_pct <= 100.0:
            raise ValueError(f"Cell {sample_id!r} has an invalid native_mapping_pct")
        cell["native_mapping_pct"] = native_pct
        expected_dedup_policy = {
            "biscuit": "dupsifter_wgbs_signature_remove_dups",
            "bismark": "umi_tools_directional_physical_r1_ignore_tlen" if row["protocol"] == "droplet" else "bismark_paired_endpoint_orientation",
            "rastair": "samtools_markdup_flag_exclude_in_rastair",
        }[backend]
        if cell["dedup_policy"] != expected_dedup_policy:
            raise ValueError(f"Cell {sample_id!r} has an invalid dedup_policy")
        expected_non_cpg_source = {
            "biscuit": "biscuit_cph_retention_by_read_position",
            "bismark": "bismark_extraction_chg_chh",
            "rastair": "not_measured_taps_cpg_only",
        }[backend]
        if cell["non_cpg_metric_source"] != expected_non_cpg_source:
            raise ValueError(f"Cell {sample_id!r} has an invalid non-CpG metric source")
        if cell["cutadapt_input_pairs"] != cell["demux_dna_reads_count"]:
            raise ValueError(f"Cell {sample_id!r} Cutadapt and demux pair counts disagree")
        if cell["trimmed_pairs"] > cell["cutadapt_input_pairs"]:
            raise ValueError(f"Cell {sample_id!r} trimmed pairs exceed input pairs")
        if (
            cell["both_primary_mapped_pairs"] + cell["one_primary_mapped_pairs"]
            != cell["backend_accepted_pairs"]
        ):
            raise ValueError(f"Cell {sample_id!r} accepted mapping categories disagree")
        if (
            cell["backend_accepted_pairs"] + cell["backend_rejected_pairs"]
            != cell["trimmed_pairs"]
        ):
            raise ValueError(f"Cell {sample_id!r} mapping funnel does not close")
        if backend == "rastair":
            if cell["ambiguous_pairs"] is not None or cell["no_primary_pairs"] is not None:
                raise ValueError(f"Cell {sample_id!r} has invalid TAPS rejection fields")
            rejection_total = cell["unmapped_pairs"] + sum(
                _required_nonnegative_int(row, key, sample_id)
                for key in ("discordant_pairs", "quality_rejected_pairs")
            )
            if row.get("high_cph_role") != "not_applicable" or not np.isnan(cell["non_cpg_methylation_pct"]):
                raise ValueError(f"Cell {sample_id!r} must declare unmeasured TAPS CpH")
            if any(cell[key] for key in ("high_cph_assessed_pairs", "high_cph_flagged_pairs", "high_cph_removed_pairs")):
                raise ValueError(f"Cell {sample_id!r} must not claim TAPS high-CpH filtering")
        elif backend == "bismark":
            if cell["no_primary_pairs"] is not None or cell["ambiguous_pairs"] is None:
                raise ValueError(f"Cell {sample_id!r} has invalid Bismark rejection fields")
            rejection_total = cell["unmapped_pairs"] + cell["ambiguous_pairs"]
        else:
            if cell["ambiguous_pairs"] is not None or cell["no_primary_pairs"] is None:
                raise ValueError(f"Cell {sample_id!r} has invalid BISCUIT rejection fields")
            rejection_total = cell["unmapped_pairs"] + cell["no_primary_pairs"]
        if rejection_total != cell["backend_rejected_pairs"]:
            raise ValueError(f"Cell {sample_id!r} rejection categories do not close")
        if cell["duplicate_pairs"] > cell["backend_accepted_pairs"]:
            raise ValueError(f"Cell {sample_id!r} duplicate pairs exceed accepted pairs")
        if (
            cell["postdedup_pairs"]
            != cell["backend_accepted_pairs"] - cell["duplicate_pairs"]
        ):
            raise ValueError(f"Cell {sample_id!r} post-dedup pair count is inconsistent")
        if not (
            cell["high_cph_removed_pairs"] <= cell["high_cph_flagged_pairs"]
            <= cell["high_cph_assessed_pairs"] <= cell["postdedup_pairs"]
        ):
            raise ValueError(f"Cell {sample_id!r} high-CpH pair counts are inconsistent")
        if (
            cell["final_retained_pairs"]
            != cell["postdedup_pairs"] - cell["high_cph_removed_pairs"]
        ):
            raise ValueError(f"Cell {sample_id!r} final pair count is inconsistent")
        active.append(cell)

    if not active:
        raise ValueError("Final sample manifest contains no completed DNA cells")
    return active


def build_qc_records(
    project_dir: str | Path,
    *,
    sealed_context: SealedQCContext | None = None,
) -> list[dict[str, Any]]:
    """由 sealed 最终 manifest 构建经校验的计算行。"""
    context = sealed_context or load_sealed_qc_context(project_dir)
    records: list[dict[str, Any]] = []
    for cell in read_active_cells(project_dir, sealed_context=context):
        sample_id = str(cell["sample_id"])
        record: dict[str, Any] = {
            "sample_id": sample_id,
            "demux_clone_id": str(cell["project_sample_id"]),
            "demux_plate_id": str(cell["plate_id"]),
            "demux_dna_barcode": str(cell["dna_barcode"]),
            "demux_rna_barcode": str(cell["rna_barcode"]),
            "demux_raw_sample": str(cell["dna_raw_sample"]),
            "demux_dna_reads_count": int(cell["demux_dna_reads_count"]),
            "cpg_file": Path(cell["cpg_file"]),
            "cpg_output_identity": dict(cell["cpg_output_identity"]),
            "methylation_backend": context.methylation_backend,
            "mapping_policy": str(cell["mapping_policy"]),
            "native_mapping_pct": float(cell["native_mapping_pct"]),
            "native_mapping_policy": str(cell["native_mapping_policy"]),
            "native_mapping_unit": str(cell["native_mapping_unit"]),
            "dedup_policy": str(cell["dedup_policy"]),
            "non_cpg_methylation_pct": float(cell["non_cpg_methylation_pct"]),
            "non_cpg_metric_source": str(cell["non_cpg_metric_source"]),
        }
        for field in PAIR_COUNT_FIELDS:
            record[field] = int(cell[field])
        record["ambiguous_pairs"] = cell["ambiguous_pairs"]
        record["no_primary_pairs"] = cell["no_primary_pairs"]
        records.append(record)
    return records


def _output_identity_for_path(
    context: SealedQCContext, path: Path
) -> Mapping[str, object]:
    public = context.project_dir / "03_results"
    try:
        key = path.resolve().relative_to(public.resolve()).as_posix()
    except ValueError as exc:
        raise ValueError(f"Downstream input escapes sealed outputs: {path}") from exc
    record = context.output_inventory.get(key)
    if record is None:
        raise ValueError(f"Downstream input is absent from delivery inventory: {key}")
    return record


def source_identity(
    project_dir: str | Path,
    records: Iterable[Mapping[str, Any]],
    *,
    sealed_context: SealedQCContext | None = None,
    extra_files: Mapping[str, str | Path] | None = None,
    context: Mapping[str, Any] | None = None,
) -> dict[str, object]:
    """返回绑定 delivery、snapshot 与输出哈希的内容身份。"""
    sealed = sealed_context or load_sealed_qc_context(project_dir)
    project = Path(project_dir).expanduser().resolve()
    if project != sealed.project_dir:
        raise ValueError("project_dir and sealed_context refer to different projects")

    output_records: dict[str, Mapping[str, object]] = {
        "QC_Results/sample_manifest.tsv": sealed.output_inventory["QC_Results/sample_manifest.tsv"]
    }
    for record in records:
        path = Path(record["cpg_file"]).resolve()
        identity = _output_identity_for_path(sealed, path)
        output_records[str(identity["path"])] = identity

    extra_entries = []
    for role, raw_path in sorted((extra_files or {}).items()):
        path = Path(raw_path).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Downstream identity input is missing: {path}")
        stat = path.stat()
        extra_entries.append(
            {
                "role": str(role),
                "path": str(path),
                "size_bytes": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
            }
        )

    payload: dict[str, object] = {
        "schema_version": 5,
        "selection": "sealed_delivery_inventory",
        "delivery_id": sealed.delivery_id,
        "snapshot_id": sealed.snapshot_id,
        "run_manifest_sha256": sealed.run_manifest_sha256,
        "snapshot_config_sha256": sealed.snapshot_config_sha256,
        "cpg_representation": sealed.cpg_representation,
        "context": dict(context or {}),
        "outputs": [dict(output_records[key]) for key in sorted(output_records)],
        "extra_files": extra_entries,
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["identity_sha256"] = hashlib.sha256(encoded).hexdigest()
    return payload


import argparse


import multiprocessing
import os
import re
import shutil
import subprocess

import tempfile
import time
import uuid
from concurrent.futures import ProcessPoolExecutor, as_completed

from functools import lru_cache


def _add_cli_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument(
        "--project-dir",
        type=str,
        required=True,
        help="Alopex project root containing 00_config and 03_results.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=0,
        help="Parallel workers; 0 chooses a conservative default.",
    )
    parser.add_argument(
        "--recompute-raw-adata",
        action="store_true",
        help="Rebuild RawAdata.h5ad from sealed CpG inputs.",
    )
    parser.add_argument(
        "--recompute-methyl-stats",
        action="store_true",
        help="Ignore cached methylation statistics.",
    )
    parser.add_argument(
        "--recompute-composition",
        action="store_true",
        help="Ignore cached read-composition statistics.",
    )
    parser.add_argument(
        "--recompute-gini",
        action="store_true",
        help="Ignore cached Gini statistics.",
    )
    parser.add_argument(
        "--recompute-tss",
        action="store_true",
        help="Ignore cached TSS profile statistics.",
    )
    parser.add_argument(
        "--tss-bed",
        type=str,
        default=None,
        help="TSS window BED; defaults to Pipeline resources/<species>_reference/tss.",
    )
    parser.add_argument(
        "--skip-tss",
        action="store_true",
        help="Generate all non-TSS outputs when TSS resources are intentionally unavailable.",
    )
    return parser


def _early_print_help_if_requested() -> None:
    if not any(arg in {"-h", "--help"} for arg in sys.argv[1:]):
        return
    _add_cli_args(argparse.ArgumentParser(description=__doc__)).print_help()
    raise SystemExit(0)


_early_print_help_if_requested()

import numpy as np
import pandas as pd


try:
    import anndata as ad
    import polars as pl
    import snapatac2 as snap
    from numba import njit
    from tqdm import tqdm
except ImportError as _qc_import_error:
    ad = pl = snap = None
    tqdm = None
    _HEAVY_IMPORT_ERROR = _qc_import_error

    def njit(*args, **kwargs):
        def _identity(fn):
            return fn
        if len(args) == 1 and callable(args[0]) and not kwargs:
            return args[0]
        return _identity


@dataclass(frozen=True)
class RunPaths:
    project_dir: Path
    cpg_dir: Path
    qc_results_dir: Path
    cache_dir: Path
    raw_adata: Path
    methyl_stats_cache: Path
    qc_info: Path
    tss_info: Path
    composition_cache: Path
    gini_cache: Path
    input_identity: Path


@dataclass(frozen=True)
class RunConfig:
    genome: str
    ref_fai: Path
    tss_bed: Path | None
    gini_bin_size: int = 500000


class Timer:
    def __init__(self) -> None:
        self._start: dict[str, float] = {}

    def start(self, name: str) -> None:
        self._start[name] = time.perf_counter()

    def stop(self, name: str) -> float | None:
        t0 = self._start.pop(name, None)
        if t0 is None:
            return None
        return time.perf_counter() - t0


def _log(message: str, level: str = "INFO") -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {level:7} {message}", flush=True)


def _pipeline_dir() -> Path:
    return Path(__file__).resolve().parents[1]


def _resolve_project_path(path: Path, project_dir: Path) -> Path:
    if path.is_absolute():
        return path.resolve()
    return (project_dir / path).resolve()


def _project_directory(
    project_dir: Path,
    relative: Path,
    *,
    create: bool,
) -> Path:

    if relative.is_absolute() or ".." in relative.parts or relative == Path("."):
        raise ValueError(f"Unsafe project directory path: {relative}")
    current = project_dir
    for part in relative.parts:
        current = current / part
        if current.is_symlink():
            raise ValueError(f"Project directory path must not contain symlinks: {current}")
        if current.exists() and not current.is_dir():
            raise ValueError(f"Project directory path is not a directory: {current}")
        if create and not current.exists():
            current.mkdir()
    return current


def _resolve_project_dir(path: Path) -> Path:
    path = path.resolve()
    if not (path / "03_results" / "run_manifest.json").is_file():
        raise FileNotFoundError(f"指定目录没有完整的 sealed Alopex delivery: {path}")
    return path


def _resolve_tss_paths(genome: str, *, pipeline_root: Path | None = None) -> Path:
    genome = str(genome).strip()
    base_dir = Path(f"resources/{genome}_reference/tss").expanduser()
    if not base_dir.is_absolute():
        base_dir = (pipeline_root or _pipeline_dir()) / base_dir
    bed = base_dir / f"{genome}_TSS_2000_2000_20.bed"
    return bed


def _has_file(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def build_paths(
    project_dir: Path,
    delivery_id: str,
) -> RunPaths:
    project_dir = project_dir.resolve()
    cpg_dir = _project_directory(
        project_dir, Path("03_results/CpG"), create=False
    )
    qc_results_dir = _project_directory(
        project_dir,
        Path("06_downstream") / delivery_id / "QC_Results",
        create=True,
    )

    cache_dir = _project_directory(
        project_dir, Path("06_downstream") / delivery_id / ".cache", create=True
    )

    return RunPaths(
        project_dir=project_dir,
        cpg_dir=cpg_dir,
        qc_results_dir=qc_results_dir,
        cache_dir=cache_dir,
        raw_adata=qc_results_dir / "RawAdata.h5ad",
        methyl_stats_cache=cache_dir / "Mapping_Stats.csv",
        qc_info=qc_results_dir / "DNAme_QC_Information.csv",
        tss_info=qc_results_dir / "TSS_Profile_Information.csv",
        composition_cache=cache_dir / "CpG_Signal_Composition.csv",
        gini_cache=cache_dir / "Gini_Index.csv",
        input_identity=cache_dir / "DNAme_QC_Input_Identity.json",
    )


def _select_cpg_cache_state(
    cached_mapping: pd.DataFrame | None,
    cached_composition: pd.DataFrame | None,
) -> pd.DataFrame | None:

    if cached_mapping is None or cached_composition is None:
        return None
    return cached_mapping.merge(
        cached_composition.drop(
            columns=["Signal_Composition_Rule", "Sample_CpG_Signals"],
            errors="ignore",
        ),
        on="Sample_ID",
        how="left",
        validate="one_to_one",
    )


def build_config(
    sealed: SealedQCContext,
    *,
    tss_bed_override: Path | None = None,
    skip_tss: bool = False,
) -> RunConfig:
    """解析 reference FAI，并在需要时强制检查共享 TSS 输入。"""
    config = sealed.config
    project_dir = sealed.project_dir
    genome = str(config.get("species") or "").strip()
    if not genome:
        raise ValueError("Alopex config is missing species")
    ref_fai = sealed.reference_fai
    if skip_tss and tss_bed_override is not None:
        raise ValueError("--skip-tss cannot be combined with explicit TSS paths")
    if skip_tss:
        tss_bed = None
    elif tss_bed_override is not None:
        tss_bed = _resolve_project_path(tss_bed_override.expanduser(), project_dir)
    else:
        tss_bed = _resolve_tss_paths(
            genome, pipeline_root=sealed.pipeline_root
        )
    cfg = RunConfig(
        genome=genome,
        ref_fai=ref_fai,
        tss_bed=tss_bed,
    )
    required_paths: list[tuple[str, Path | None]] = [("reference .fai", cfg.ref_fai)]
    if not skip_tss:
        required_paths.append(("TSS BED", cfg.tss_bed))
    missing = [
        f"{label}: {path}"
        for label, path in required_paths
        if path is None or not _has_file(path)
    ]
    if missing:
        raise FileNotFoundError(
            "Alopex downstream Reference/TSS input is incomplete. "
            "Provide --tss-bed or explicitly use --skip-tss:\n"
            + "\n".join(f"  - {value}" for value in missing)
        )
    return cfg


COMPOSITION_COV_COLS = (
    "lambda_cov",
    "puc19_cov",
    "mt_cov",
    "host_cov",
    "total_cov",
    "Sample_CpG_Signals",
)
COMPOSITION_FRACTION_COLS = {
    "lambda_cov": "_f_lambda",
    "puc19_cov": "_f_puc",
    "mt_cov": "_f_mt",
    "host_cov": "_f_host",
}
MAPPING_STAT_COLUMNS = (
    "Sample_ID",
    "Sample_Unique_CpG_Sites",
    "Sample_Meth_CpG_Rate%",
    "Lambda_Unique_CpG_Sites",
    "Lambda_Meth_CpG_Rate%",
    "pUC19_Unique_CpG_Sites",
    "pUC19_Meth_CpG_Rate%",
)
COMPOSITION_CACHE_COLUMNS = (
    "Sample_ID",
    "Signal_Host_Rate%",
    "Signal_Lambda_Rate%",
    "Signal_pUC19_Rate%",
    "Signal_mtDNA_Rate%",
    "Sample_CpG_Signals",
    "Signal_Composition_Rule",
)
TSS_REQUIRED_COLUMNS = (
    "Sample_ID",
    "position_coarse",
    "meth_frac",
    "count",
    "CloneID",
)
GINI_CACHE_COLUMNS = ("Sample_ID", "Gini_Index")
QC_INFO_COLUMNS = (
    "Sample_ID",
    "CloneID",
    "PlateID",
    "Native_Mapping%",
    "Sample_Unique_CpG_Sites",
    "Sample_Meth_CpG_Rate%",
    "Lambda_Unique_CpG_Sites",
    "Lambda_Meth_CpG_Rate%",
    "pUC19_Unique_CpG_Sites",
    "pUC19_Meth_CpG_Rate%",
    "Non_CpG_Methylation%",
    "High_CpH_Flag_Rate%",
    "Signal_Host_Rate%",
    "Signal_Lambda_Rate%",
    "Signal_pUC19_Rate%",
    "Signal_mtDNA_Rate%",
    "Gini_Index",
    "Duplicate_Pair_Rate%",
    "Trim_Retention%",
    "Final_Pair_Yield%",
)
RAW_ADATA_SCHEMA_VERSION = 4
RAW_ADATA_SCHEMA_MODE = "strand_resolved_headered_empty_safe"
RAW_ADATA_ALGORITHM_VERSION = "snapatac2_import_values_0based_v5"
EMPTY_CPG_SENTINEL_CONTIG = "__DNA_PIPELINE_EMPTY__"

MAX_CPG_PROCESS_WORKERS = 32
_THREAD_LIMIT_ENV_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "RAYON_NUM_THREADS",
    "POLARS_MAX_THREADS",
)


@lru_cache(maxsize=4096)
def _classify_seq_name(seq: str) -> tuple[bool, bool, bool, bool]:
    seq = str(seq).strip()
    is_lambda = bool(re.search(r"(?i)lambda", seq))
    is_puc = bool(re.search(r"(?i)puc19", seq))
    is_mt = bool(re.fullmatch(r"(?i)(chrM|MT|chrMT)", seq))
    is_host = seq.startswith("chr") and not (is_lambda or is_puc or is_mt)
    return is_lambda, is_puc, is_mt, is_host


def _initialize_cpg_worker() -> None:

    for variable in _THREAD_LIMIT_ENV_VARS:
        os.environ[variable] = "1"


def _bounded_process_workers(requested: int, task_count: int) -> int:

    if task_count <= 0:
        return 0
    available = max(1, os.cpu_count() or 1)
    wanted = int(requested) if int(requested) > 0 else available
    return max(1, min(wanted, MAX_CPG_PROCESS_WORKERS, task_count))


def _process_map_ordered(function, tasks: list[tuple], worker_count: int) -> list:

    if not tasks:
        return []
    if worker_count <= 0:
        raise ValueError("worker_count must be positive when tasks are present")

    previous_worker_env = {
        variable: os.environ.get(variable)
        for variable in _THREAD_LIMIT_ENV_VARS
    }


    for variable in _THREAD_LIMIT_ENV_VARS:
        os.environ[variable] = "1"

    executor = None
    futures = {}
    try:
        executor = ProcessPoolExecutor(
            max_workers=worker_count,
            mp_context=multiprocessing.get_context("spawn"),
            initializer=_initialize_cpg_worker,
        )
        futures = {
            executor.submit(function, *arguments): index
            for index, arguments in enumerate(tasks)
        }
        rows = [None] * len(tasks)
        for future in as_completed(futures):
            index = futures[future]
            try:
                rows[index] = future.result()
            except BaseException:
                for pending in futures:
                    if pending is not future:
                        pending.cancel()
                raise
        return rows
    except BaseException:
        for pending in futures:
            pending.cancel()
        raise
    finally:


        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
        for variable, value in previous_worker_env.items():
            if value is None:
                os.environ.pop(variable, None)
            else:
                os.environ[variable] = value


def _atomic_write_csv(df: pd.DataFrame, path: Path) -> None:

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            prefix=f".{path.name}.tmp.",
            suffix=".csv",
            dir=str(path.parent),
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            df.to_csv(handle, index=False)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        temporary = None
        _fsync_directory(path.parent)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _fsync_directory(path: Path) -> None:

    directory_fd = os.open(str(Path(path)), os.O_RDONLY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)


def _atomic_write_json(payload: dict, path: Path) -> None:

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            prefix=f".{path.name}.tmp.",
            suffix=".json",
            dir=str(path.parent),
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        temporary = None
        _fsync_directory(path.parent)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _fsync_file(path: Path) -> None:

    file_fd = os.open(str(Path(path)), os.O_RDONLY)
    try:
        os.fsync(file_fd)
    finally:
        os.close(file_fd)


def _raw_adata_stat_binding(path: Path) -> dict:

    path = Path(path)
    stat = path.stat()
    return {
        "name": path.name,
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "device": stat.st_dev,
        "inode": stat.st_ino,
    }


def _content_binding(path: Path) -> dict:
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
    }


def _raw_adata_schema_matches(
    out_file: Path, schema_payload: object, input_signature: dict | None = None
) -> bool:

    if not isinstance(schema_payload, dict) or not Path(out_file).is_file():
        return False
    generation_id = schema_payload.get("generation_id")
    if re.fullmatch(r"[0-9a-f]{32}", str(generation_id or "")) is None:
        return False
    if schema_payload.get("h5ad") != _raw_adata_stat_binding(out_file):
        return False
    if input_signature is not None:
        observed_input = {
            key: value
            for key, value in schema_payload.items()
            if key not in {"generation_id", "h5ad"}
        }
        if observed_input != input_signature:
            return False
    reference_fai = schema_payload.get("reference_fai")
    cpg_files = schema_payload.get("cpg_files")
    if (
        not isinstance(reference_fai, dict)
        or re.fullmatch(r"[0-9a-f]{64}", str(reference_fai.get("sha256", ""))) is None
        or not isinstance(cpg_files, list)
        or any(
            not isinstance(record, dict)
            or re.fullmatch(r"[0-9a-f]{64}", str(record.get("sha256", ""))) is None
            for record in cpg_files
        )
    ):
        return False
    return (
        schema_payload.get("schema_version") == RAW_ADATA_SCHEMA_VERSION
        and schema_payload.get("algorithm_version") == RAW_ADATA_ALGORITHM_VERSION
        and re.fullmatch(r"[0-9a-f]{64}", str(schema_payload.get("delivery_id", "")))
        is not None
        and re.fullmatch(r"[0-9a-f]{64}", str(schema_payload.get("snapshot_id", "")))
        is not None
        and isinstance(schema_payload.get("sample_ids"), list)
        and schema_payload.get("cpg_header") == "snapatac2_import_values_0based"
        and schema_payload.get("empty_cpg_sentinel")
        == f"{EMPTY_CPG_SENTINEL_CONTIG}:0:0:1"
    )


def _publish_raw_adata_generation(
    temporary: Path,
    out_file: Path,
    schema_file: Path,
    input_signature: dict,
    *,
    generation_id: str | None = None,
) -> dict:


    temporary = Path(temporary)
    out_file = Path(out_file)
    schema_file = Path(schema_file)
    if temporary.parent != out_file.parent or schema_file.parent != out_file.parent:
        raise ValueError("RawAdata temporary, output, and schema must share a directory")
    generation_id = generation_id or uuid.uuid4().hex
    if re.fullmatch(r"[0-9a-f]{32}", generation_id) is None:
        raise ValueError("RawAdata generation_id must be 32 lowercase hex characters")

    _fsync_file(temporary)
    schema_file.unlink(missing_ok=True)
    _fsync_directory(out_file.parent)
    os.replace(temporary, out_file)
    _fsync_directory(out_file.parent)

    schema_payload = {
        **input_signature,
        "generation_id": generation_id,
        "h5ad": _raw_adata_stat_binding(out_file),
    }
    _atomic_write_json(schema_payload, schema_file)
    return schema_payload


def _read_chrom_sizes(fai_path: Path) -> dict[str, int]:

    chrom_sizes: dict[str, int] = {}
    with Path(fai_path).open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"{fai_path}:{line_no}: malformed FAI row")
            chrom = fields[0].strip()
            try:
                length = int(fields[1])
            except ValueError as exc:
                raise ValueError(
                    f"{fai_path}:{line_no}: non-integer chromosome length"
                ) from exc
            if not chrom or length <= 0:
                raise ValueError(f"{fai_path}:{line_no}: invalid chromosome entry")
            if chrom == EMPTY_CPG_SENTINEL_CONTIG:
                raise ValueError(
                    f"{fai_path}:{line_no}: reference uses reserved empty-CpG "
                    f"sentinel contig {EMPTY_CPG_SENTINEL_CONTIG}"
                )


            if "." in chrom:
                continue
            if chrom in chrom_sizes:
                raise ValueError(f"{fai_path}:{line_no}: duplicate chromosome {chrom}")
            chrom_sizes[chrom] = length
    if not chrom_sizes:
        raise ValueError(f"Reference FAI contains no usable chromosomes: {fai_path}")
    return chrom_sizes


def _require_exact_sample_ids(
    df: pd.DataFrame, sample_ids: list[str], label: str
) -> pd.DataFrame:
    if "Sample_ID" not in df.columns:
        raise ValueError(f"{label} is missing Sample_ID")
    wanted = [str(value) for value in sample_ids]
    if len(wanted) != len(set(wanted)):
        raise ValueError("Current manifest contains duplicate Sample_ID values")
    out = df.copy()
    out["Sample_ID"] = out["Sample_ID"].astype(str).str.strip()
    observed = out["Sample_ID"].tolist()
    if any(not value for value in observed) or len(observed) != len(set(observed)):
        raise ValueError(f"{label} contains blank or duplicate Sample_ID values")
    if len(observed) != len(wanted) or set(observed) != set(wanted):
        raise ValueError(f"{label} Sample_ID set differs from current manifest")
    return out.set_index("Sample_ID").loc[wanted].reset_index()


def _mapped_numeric(sid: pd.Series, values: pd.Series) -> np.ndarray:
    return (
        pd.to_numeric(sid.map(values), errors="coerce")
        .fillna(0.0)
        .to_numpy(dtype=float)
    )


def prepare_flat_cpg_symlink_dir(
    temporary_root: Path,
    selected_files: list[Path],
) -> tuple[Path, dict[str, str]]:
    temporary_root = Path(temporary_root).resolve()
    flat_dir = temporary_root / "cpg_symlinks"
    flat_dir.mkdir(parents=True, exist_ok=False)

    files = [Path(path).resolve() for path in selected_files]
    for path in files:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Selected CpG result is missing or empty: {path}")

    linkname2paths = {}
    linkname2clean = {}

    for p in files:
        base = p.name
        link_base = base.removesuffix(".cpg.tsv.zst") + ".tsv.zst"
        linkname2paths.setdefault(link_base, []).append(p)

    for link_base, ps in linkname2paths.items():
        if len(ps) != 1:
            raise ValueError(
                f"More than one active CpG file maps to {link_base!r}: "
                + ", ".join(str(path) for path in ps)
            )
        src = ps[0]
        link_path = flat_dir / link_base
        os.symlink(src.resolve(), link_path)

        clean = _clean_cpg_table_name(link_base)

        linkname2clean[link_base] = clean
        if link_base.endswith(".tsv.zst"):
            linkname2clean[link_base[: -len(".zst")]] = clean

    return flat_dir, linkname2clean


def _clean_cpg_table_name(filename: str) -> str:
    for suffix in (".tsv.zst", ".tsv"):
        if filename.endswith(suffix):
            return filename[: -len(suffix)]
    return filename


def load_or_create_raw_adata(
    source_dir: Path,
    out_file: Path,
    fai_path: Path | None,
    linkname2clean: dict,
    *,
    reuse_existing: bool,
    expected_sample_ids: list[str],
    empty_cpg_sample_ids: set[str],
    sealed_delivery_id: str,
    sealed_snapshot_id: str,
    cpg_output_identities: Mapping[str, Mapping[str, object]],
):
    """按 sealed inventory 绑定 CpG 输入并复用或重建 RawAdata，不重复扫描 CpG 计算哈希。"""
    expected = [str(value) for value in expected_sample_ids]
    empty_cpg = {str(value) for value in empty_cpg_sample_ids}
    if len(expected) != len(set(expected)):
        raise ValueError("Current manifest contains duplicate Sample_ID values")
    if empty_cpg - set(expected):
        raise ValueError("Empty-CpG Sample_ID set differs from current manifest")
    if set(cpg_output_identities) != set(expected):
        raise ValueError("Sealed CpG inventory Sample_ID set differs from current manifest")
    for label, value in (
        ("sealed_delivery_id", sealed_delivery_id),
        ("sealed_snapshot_id", sealed_snapshot_id),
    ):
        if re.fullmatch(r"[0-9a-f]{64}", str(value)) is None:
            raise ValueError(f"{label} must be a 64-character content identity")
    schema_file = out_file.with_name(f"{out_file.stem}.input_schema.json")
    cpg_bindings = []
    for path in sorted(source_dir.iterdir()):
        if not path.is_symlink():
            continue
        clean_name = linkname2clean.get(path.name, _clean_cpg_table_name(path.name))
        record = cpg_output_identities[clean_name]
        cpg_bindings.append({
            "sample_id": clean_name,
            "path": str(path.resolve()),
            "size_bytes": record["size_bytes"],
            "sha256": record["sha256"],
        })
    input_signature = {
        "schema_version": RAW_ADATA_SCHEMA_VERSION,
        "algorithm_version": RAW_ADATA_ALGORITHM_VERSION,
        "algorithm_sha256": _algorithm_fingerprints()["raw"],
        "delivery_id": str(sealed_delivery_id),
        "snapshot_id": str(sealed_snapshot_id),
        "cpg_header": "snapatac2_import_values_0based",
        "empty_cpg_sentinel": f"{EMPTY_CPG_SENTINEL_CONTIG}:0:0:1",
        "sample_ids": sorted(expected),
        "empty_cpg_sample_ids": sorted(empty_cpg),
        "reference_fai": _content_binding(fai_path) if fai_path is not None else None,
        "cpg_files": cpg_bindings,
    }
    schema_payload: object = None
    if schema_file.is_file():
        try:
            schema_payload = json.loads(schema_file.read_text(encoding="utf-8"))
        except (OSError, ValueError, TypeError):
            schema_payload = None
    schema_matches = _raw_adata_schema_matches(
        out_file, schema_payload, input_signature
    )

    if reuse_existing and out_file.exists() and out_file.is_file():
        existing = None
        try:
            existing = snap.read(str(out_file), backed="r")
            observed = [str(value) for value in existing.obs_names]
            if len(observed) != len(set(observed)) or set(observed) != set(expected):
                raise ValueError(
                    "RawAdata Sample_ID set differs from current final manifest"
                )
            if schema_matches:
                validated = existing
                existing = None
                return validated
            raise ValueError(
                "RawAdata does not match the current SnapATAC2 0-based four-column input schema"
            )
        except Exception as e:
            _log(f"读取现有 h5ad 失败，尝试重建: {e}", level="WARNING")
        finally:
            if existing is not None:
                existing.close()

    if fai_path is None or not fai_path.exists():
        return None

    chrom_sizes = _read_chrom_sizes(fai_path)

    out_file.parent.mkdir(parents=True, exist_ok=True)
    import_dir = Path(
        tempfile.mkdtemp(prefix=".snapatac-cpg-input.", dir=str(out_file.parent))
    )
    temporary = out_file.with_name(
        f".{out_file.stem}.tmp.{os.getpid()}.{time.time_ns()}{out_file.suffix}"
    )
    snap_data = None
    try:


        for source in sorted(source_dir.iterdir()):
            if not source.is_symlink():
                continue
            sample_id = linkname2clean.get(
                source.name, _clean_cpg_table_name(source.name)
            )
            target = import_dir / source.name
            if sample_id not in empty_cpg:
                os.symlink(source.resolve(), target)
                continue


            target = import_dir / source.name.removesuffix(".zst")
            target.write_text(
                "chrom\tpos\tmethyl\tunmethyl\n"
                f"{EMPTY_CPG_SENTINEL_CONTIG}\t0\t0\t1\n",
                encoding="utf-8",
            )

        snap_data = snap.pp.import_values(
            str(import_dir),
            chrom_sizes=chrom_sizes,
            file=str(temporary),
            backend="hdf5",
            chunk_size=8,
        )
        new_names = [
            linkname2clean.get(
                bn := Path(str(x)).name,
                _clean_cpg_table_name(bn),
            )
            for x in snap_data.obs_names
        ]
        if len(new_names) != len(set(new_names)):
            raise ValueError("Active CpG inputs produced duplicate Sample_ID values")
        if len(new_names) != len(expected) or set(new_names) != set(expected):
            raise ValueError(
                "New RawAdata Sample_ID set differs from current final manifest"
            )
        snap_data.obs_names = new_names
        snap_data.close()
        snap_data = None
        _publish_raw_adata_generation(
            temporary, out_file, schema_file, input_signature
        )
    finally:
        if snap_data is not None:
            try:
                snap_data.close()
            except Exception as exc:
                _log(f"关闭未完成的 h5ad 失败: {exc}", level="WARNING")
        temporary.unlink(missing_ok=True)
        shutil.rmtree(import_dir, ignore_errors=True)

    return snap.read(str(out_file), backed="r")


def _raw_adata_schema_mode(out_file: Path) -> str:

    schema_file = Path(out_file).with_name(f"{Path(out_file).stem}.input_schema.json")
    if not Path(out_file).is_file():
        return "missing"
    if not schema_file.is_file():
        return "invalid_schema"
    try:
        payload = json.loads(schema_file.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return "invalid_schema"
    if _raw_adata_schema_matches(out_file, payload):
        return RAW_ADATA_SCHEMA_MODE
    return "invalid_schema"


def _algorithm_fingerprints() -> dict[str, str]:

    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if node.body and isinstance(node.body[0], ast.Expr) and isinstance(node.body[0].value, ast.Constant) and isinstance(node.body[0].value.value, str):
                node.body.pop(0)
    bindings = {}

    def bind(node, owner=None):
        owner = node if owner is None else owner
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            bindings[node.name] = owner
        elif isinstance(node, (ast.Assign, ast.AnnAssign, ast.Import, ast.ImportFrom)):
            for child in ast.walk(node):
                if isinstance(child, ast.Name) and isinstance(child.ctx, ast.Store):
                    bindings[child.id] = owner
                elif isinstance(child, ast.alias):
                    bindings[child.asname or child.name.split(".")[0]] = owner
        else:
            for child in ast.iter_child_nodes(node):
                bind(child, owner)

    for node in tree.body:
        bind(node)
    versions = {name: importlib.metadata.version(name)
                for name in ("numpy", "pandas", "scipy", "numba", "snapatac2")}
    roots = {
        "input": ("build_qc_records", "source_identity"),
        "raw": ("load_or_create_raw_adata",),
        "cpg": ("compute_methylation_statistics", "compute_cpg_signal_composition"),
        "gini": ("compute_gini_index",),
        "tss": ("compute_tss_profile",),
    }
    fingerprints = {}
    for stage, names in roots.items():
        pending, selected = list(names), {}
        while pending:
            name = pending.pop()
            if name in selected or name not in bindings:
                continue
            node = bindings[name]
            selected[name] = ast.dump(node, include_attributes=False)
            pending.extend(child.id for child in ast.walk(node)
                           if isinstance(child, ast.Name) and isinstance(child.ctx, ast.Load))
        payload = json.dumps({"code": selected, "libraries": versions}, sort_keys=True)
        fingerprints[stage] = hashlib.sha256(payload.encode()).hexdigest()
    return fingerprints


def _cache_identities(identity: Mapping[str, object], cfg, algorithms: Mapping[str, str]) -> dict[str, str]:

    common = {"input": identity["identity_sha256"], "reader": algorithms["input"]}
    payloads = {
        "cpg": {**common, "algorithm": algorithms["cpg"]},
        "gini": {**common, "algorithm": algorithms["gini"], "raw": algorithms["raw"],
                 "bin_size": cfg.gini_bin_size},
        "tss": {**common, "algorithm": algorithms["tss"], "raw": algorithms["raw"],
                "tss_bed": _content_binding(cfg.tss_bed) if cfg.tss_bed is not None else None},
    }
    return {key: hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()
            for key, value in payloads.items()}


def _open_cpg_text(path: Path):

    if path.name.endswith(".zst"):
        executable = shutil.which("zstdcat")
        if executable is None:
            raise FileNotFoundError("zstdcat is required to read CpG .zst files")
        process = subprocess.Popen(
            [executable, str(path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if process.stdout is None:
            process.kill()
            raise RuntimeError(f"Cannot open zstd stream: {path}")
        return process, process.stdout
    raise ValueError(
        f"CpG input must use the current .cpg.tsv.zst contract: {path}"
    )


def _finish_cpg_process(
    process: subprocess.Popen | None, path: Path, stream_error: BaseException | None
) -> None:
    if process is None:
        return
    stderr = ""
    if process.stderr is not None:
        stderr = process.stderr.read()
        process.stderr.close()
    return_code = process.wait()
    if return_code and stream_error is None:
        raise RuntimeError(
            f"zstdcat failed for {path} (exit {return_code}): {stderr.strip()}"
        )


def _cpg_table_is_empty(path: Path) -> bool:

    process, handle = _open_cpg_text(path)
    stream_error: BaseException | None = None
    found_data = False
    try:
        header_seen = False
        for line in handle:
            if not line.strip():
                continue
            if not header_seen:
                if line.rstrip("\n") != "chrom\tpos\tmethyl\tunmethyl":
                    raise ValueError(f"CpG input has an invalid header: {path}")
                header_seen = True
                continue
            found_data = True
            break
        if not header_seen:
            raise ValueError(f"CpG input has no header: {path}")
    except BaseException as exc:
        stream_error = exc
        raise
    finally:
        handle.close()
        if found_data and process is not None:


            if process.poll() is None:
                process.terminate()
            if process.stderr is not None:
                process.stderr.read()
                process.stderr.close()
            process.wait()
        else:
            _finish_cpg_process(process, path, stream_error)
    return not found_data


def _kahan_add(total: float, correction: float, value: float) -> tuple[float, float]:
    adjusted = value - correction
    updated = total + adjusted
    return updated, (updated - total) - adjusted


def _stream_one_cpg_summary(
    path: Path,
    sample_id: str,
    allowed_chroms: frozenset[str],
    allow_empty: bool = False,
) -> dict:

    if not path.is_file():
        raise FileNotFoundError(f"CpG input does not exist: {path}")
    counts = {"host": 0, "lambda": 0, "pUC19": 0}
    fraction_sums = {"host": 0.0, "lambda": 0.0, "pUC19": 0.0}
    fraction_corrections = {"host": 0.0, "lambda": 0.0, "pUC19": 0.0}
    coverage = {
        "lambda_cov": 0,
        "puc19_cov": 0,
        "mt_cov": 0,
        "host_cov": 0,
        "total_cov": 0,
    }
    process, handle = _open_cpg_text(path)
    stream_error: BaseException | None = None
    try:
        for line_no, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if not line.strip():
                continue
            if line_no == 1:
                if fields != ["chrom", "pos", "methyl", "unmethyl"]:
                    raise ValueError(f"{path}:{line_no}: invalid CpG header")
                continue
            if len(fields) != 4:
                raise ValueError(f"{path}:{line_no}: expected 4 tab-separated fields")
            chrom = fields[0].strip()
            try:
                position = int(fields[1])
                methylated = int(fields[2])
                unmethylated = int(fields[3])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_no}: non-integer CpG field") from exc
            total = methylated + unmethylated
            if position < 0 or methylated < 0 or unmethylated < 0 or total <= 0:
                raise ValueError(f"{path}:{line_no}: invalid CpG counts")
            if max(methylated, unmethylated, total) > 65535:
                raise ValueError(
                    f"{path}:{line_no}: CpG counts exceed the SnapATAC2 u16 range"
                )

            if chrom not in allowed_chroms:
                continue
            is_lambda_comp, is_puc_comp, is_mt_comp, is_host_comp = (
                _classify_seq_name(chrom)
            )
            coverage["total_cov"] += total
            if is_lambda_comp:
                coverage["lambda_cov"] += total
            elif is_puc_comp:
                coverage["puc19_cov"] += total
            elif is_mt_comp:
                coverage["mt_cov"] += total
            elif is_host_comp:
                coverage["host_cov"] += total

            is_lambda = chrom == "lambda"
            is_puc = chrom == "pUC19"
            is_host = chrom.startswith("chr") and not re.fullmatch(
                r"(?i)(chrM|MT|chrMT)", chrom
            )
            category = (
                "lambda" if is_lambda else "pUC19" if is_puc else "host" if is_host else None
            )
            if category is None:
                continue
            counts[category] += 1
            fraction = (
                0.0 if methylated == 0 else float(methylated) / float(total)
            )
            fraction_sums[category], fraction_corrections[category] = _kahan_add(
                fraction_sums[category],
                fraction_corrections[category],
                fraction,
            )
    except BaseException as exc:
        stream_error = exc
        raise
    finally:
        handle.close()
        _finish_cpg_process(process, path, stream_error)

    def mean_percent(category: str) -> float:
        count = counts[category]
        return fraction_sums[category] / count * 100.0 if count else float("nan")

    if coverage["total_cov"] <= 0 and not allow_empty:
        raise ValueError(f"CpG input contains no data rows: {path}")

    return {
        "Sample_ID": str(sample_id),
        "Sample_Unique_CpG_Sites": counts["host"],
        "Sample_Meth_CpG_Rate%": mean_percent("host"),
        "Lambda_Unique_CpG_Sites": counts["lambda"],
        "Lambda_Meth_CpG_Rate%": mean_percent("lambda"),
        "pUC19_Unique_CpG_Sites": counts["pUC19"],
        "pUC19_Meth_CpG_Rate%": mean_percent("pUC19"),
        **coverage,
        "Sample_CpG_Signals": coverage["host_cov"],
    }


def compute_methylation_statistics(
    cpg_files: list[Path],
    sample_ids: list[str],
    *,
    allowed_chroms: frozenset[str],
    workers: int,
    return_composition: bool = False,
    empty_cpg_sample_ids: set[str] | None = None,
) -> pd.DataFrame | tuple[pd.DataFrame, pd.DataFrame]:
    if len(cpg_files) != len(sample_ids):
        raise ValueError("CpG file count and Sample_ID count differ")
    if not cpg_files:
        raise ValueError("No CpG files were supplied")
    empty_cpg = {str(value) for value in (empty_cpg_sample_ids or set())}
    if empty_cpg - set(map(str, sample_ids)):
        raise ValueError("Empty-CpG Sample_ID set differs from CpG inputs")
    worker_count = _bounded_process_workers(workers, len(cpg_files))
    _log(
        f"正在从 CpG 表流式提取统计信息 (workers: {worker_count})...",
        level="DEBUG",
    )
    rows = _process_map_ordered(
        _stream_one_cpg_summary,
        [
            (path, sample_id, allowed_chroms, str(sample_id) in empty_cpg)
            for path, sample_id in zip(cpg_files, sample_ids, strict=True)
        ],
        worker_count,
    )
    combined = pd.DataFrame.from_records(rows)
    mapping = combined[list(MAPPING_STAT_COLUMNS)].copy()
    if not return_composition:
        return mapping
    composition = combined[["Sample_ID", *COMPOSITION_COV_COLS]].copy()
    return mapping, composition


def _load_methyl_stats_cache(
    cache_path: Path,
    sample_ids: list[str],
    empty_cpg_sample_ids: set[str],
) -> pd.DataFrame | None:
    if not cache_path.exists():
        return None
    try:
        df = pd.read_csv(cache_path, float_precision="round_trip")
    except Exception:
        return None
    if df.empty or tuple(df.columns) != MAPPING_STAT_COLUMNS:
        return None
    wanted = [str(x) for x in sample_ids]
    try:
        df = _require_exact_sample_ids(df, wanted, "Mapping_Stats cache")
        count_cols = [
            "Sample_Unique_CpG_Sites",
            "Lambda_Unique_CpG_Sites",
            "pUC19_Unique_CpG_Sites",
        ]
        for column in count_cols:
            values = pd.to_numeric(df[column], errors="raise").to_numpy(dtype=float)
            if (
                not np.isfinite(values).all()
                or (values < 0).any()
                or not np.equal(values, np.floor(values)).all()
            ):
                return None
        for prefix in ("Sample", "Lambda", "pUC19"):
            counts = pd.to_numeric(
                df[f"{prefix}_Unique_CpG_Sites"], errors="raise"
            ).to_numpy(dtype=float)
            rates = pd.to_numeric(
                df[f"{prefix}_Meth_CpG_Rate%"], errors="raise"
            ).to_numpy(dtype=float)
            if np.any((counts > 0) & (~np.isfinite(rates))):
                return None
            if np.any((counts == 0) & (~np.isnan(rates))):
                return None
            finite = rates[np.isfinite(rates)]
            if (finite < 0).any() or (finite > 100).any():
                return None
        empty_mask = df["Sample_ID"].astype(str).isin(empty_cpg_sample_ids)
        if not (df.loc[empty_mask, count_cols] == 0).all(axis=None):
            return None
    except Exception:
        return None
    return df


def _validated_composition_cache(
    cache_path: Path,
    sample_ids: list[str],
    empty_cpg_sample_ids: set[str],
) -> pd.DataFrame | None:
    if not cache_path.is_file():
        return None
    try:
        cached = pd.read_csv(cache_path, float_precision="round_trip")
        if tuple(cached.columns) != COMPOSITION_CACHE_COLUMNS:
            return None
        cached = _require_exact_sample_ids(
            cached, sample_ids, "CpG_Signal_Composition cache"
        )
        if not cached["Signal_Composition_Rule"].astype(str).eq(
            "host_chr_prefix"
        ).all():
            return None
        required_finite = [
            "Signal_Host_Rate%",
            "Signal_Lambda_Rate%",
            "Signal_pUC19_Rate%",
            "Signal_mtDNA_Rate%",
            "Sample_CpG_Signals",
        ]
        for column in required_finite:
            values = pd.to_numeric(cached[column], errors="raise").to_numpy(dtype=float)
            if not np.isfinite(values).all() or (values < 0).any():
                return None
            if column.endswith("Rate%") and (values > 100).any():
                return None
            if column == "Sample_CpG_Signals" and not np.equal(
                values, np.floor(values)
            ).all():
                return None
            cached[column] = values
        rate_sum = cached[
            [
                "Signal_Host_Rate%",
                "Signal_Lambda_Rate%",
                "Signal_pUC19_Rate%",
                "Signal_mtDNA_Rate%",
            ]
        ].sum(axis=1)
        if (rate_sum > 100.0 + 1e-9).any():
            return None
        empty_mask = cached["Sample_ID"].astype(str).isin(empty_cpg_sample_ids)
        if not (
            cached.loc[empty_mask, required_finite].to_numpy(dtype=float) == 0
        ).all():
            return None
        return cached
    except Exception:
        return None


def compute_cpg_signal_composition(
    df_qc: pd.DataFrame,
    cache_path: Path,
    composition_stats: pd.DataFrame,
    empty_cpg_sample_ids: set[str],
) -> pd.DataFrame:
    if df_qc.empty:
        return df_qc

    df_qc = _require_exact_sample_ids(
        df_qc, df_qc["Sample_ID"].astype(str).tolist(), "QC table"
    )
    sample_ids = df_qc["Sample_ID"].tolist()
    comp = composition_stats.copy()

    comp = _require_exact_sample_ids(comp, sample_ids, "CpG composition statistics")
    missing_cov = set(COMPOSITION_COV_COLS) - set(comp.columns)
    if missing_cov:
        raise ValueError(
            "CpG composition statistics missing columns: "
            + ", ".join(sorted(missing_cov))
        )
    for column in COMPOSITION_COV_COLS:
        values = pd.to_numeric(comp[column], errors="raise").to_numpy(dtype=float)
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(f"Invalid CpG composition values in {column}")
        comp[column] = values
    empty_mask = comp["Sample_ID"].astype(str).isin(empty_cpg_sample_ids)
    if (comp.loc[~empty_mask, "total_cov"] <= 0).any():
        raise ValueError("Non-empty CpG composition has no coverage")
    if not (comp.loc[empty_mask, COMPOSITION_COV_COLS] == 0).all(axis=None):
        raise ValueError("Empty CpG composition contains coverage")
    classified = comp[["lambda_cov", "puc19_cov", "mt_cov", "host_cov"]].sum(axis=1)
    if (classified > comp["total_cov"]).any():
        raise ValueError("Classified CpG coverage exceeds total coverage")
    if not np.allclose(comp["Sample_CpG_Signals"], comp["host_cov"], rtol=0, atol=0):
        raise ValueError("Sample_CpG_Signals differs from host CpG coverage")

    total_cov = comp["total_cov"].astype(float).replace(0.0, np.nan)
    for source_col, frac_col in COMPOSITION_FRACTION_COLS.items():
        comp[frac_col] = pd.to_numeric(
            comp.get(source_col, 0.0) / total_cov, errors="coerce"
        ).fillna(0.0)
    comp = comp.set_index("Sample_ID", drop=False)

    base = df_qc[["Sample_ID"]].copy()
    base["Sample_ID"] = base["Sample_ID"].astype(str)
    sid = base["Sample_ID"]
    base["Sample_CpG_Signals"] = pd.to_numeric(
        sid.map(comp["Sample_CpG_Signals"]), errors="coerce"
    ).fillna(0.0)
    f_lambda = _mapped_numeric(sid, comp["_f_lambda"])
    f_puc = _mapped_numeric(sid, comp["_f_puc"])
    f_mt = _mapped_numeric(sid, comp["_f_mt"])
    f_host = _mapped_numeric(sid, comp["_f_host"])

    base["Signal_Lambda_Rate%"] = 100.0 * f_lambda
    base["Signal_pUC19_Rate%"] = 100.0 * f_puc
    base["Signal_mtDNA_Rate%"] = 100.0 * f_mt
    base["Signal_Host_Rate%"] = 100.0 * f_host
    base["Signal_Composition_Rule"] = "host_chr_prefix"
    comp_out = base[list(COMPOSITION_CACHE_COLUMNS)].copy()
    _atomic_write_csv(comp_out, cache_path)

    comp_merge = comp_out.drop(
        columns=["Signal_Composition_Rule", "Sample_CpG_Signals"], errors="ignore"
    )
    return df_qc.merge(
        comp_merge, on="Sample_ID", how="left", validate="one_to_one"
    )


@njit
def _gini_from_sorted_positive(x_sorted: np.ndarray) -> float:
    n = x_sorted.size
    if n == 0:
        return float("nan")
    cumx = np.cumsum(x_sorted)
    return float((n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n)


def _aggregate_tss_matrix(
    tss_mat: ad.AnnData,
    pos_map: pd.Series,
    value_name: str,
    agg_expr: str,
    *,
    chunk_size: int = 4,
) -> pl.DataFrame:
    pos_arr = pos_map.reindex(pd.Index(tss_mat.var_names, dtype=str)).to_numpy()
    obs_names = np.asarray(tss_mat.obs_names, dtype=str)
    parts = []
    for matrix, start, end in tss_mat.chunked_X(chunk_size):
        coo = matrix.tocoo(copy=False)
        mask = ~pd.isna(pos_arr[coo.col])
        if not np.any(mask):
            continue
        parts.append(
            pl.DataFrame(
                {
                    "Sample_ID": obs_names[start:end][coo.row[mask]],
                    "position_coarse": pos_arr[coo.col[mask]].astype(
                        np.int32, copy=False
                    ),
                    value_name: coo.data[mask].astype(np.float64, copy=False),
                }
            )
            .group_by(["Sample_ID", "position_coarse"])
            .agg(getattr(pl.col(value_name), agg_expr)())
        )
    if not parts:
        return pl.DataFrame(
            schema={
                "Sample_ID": pl.Utf8,
                "position_coarse": pl.Int32,
                value_name: pl.Float64,
            }
        )
    return (
        pl.concat(parts)
        .group_by(["Sample_ID", "position_coarse"])
        .agg(getattr(pl.col(value_name), agg_expr)())
    )


def _read_tss_offset_map(tss_bed_path: Path) -> pd.Series:
    pos_df = pd.read_csv(tss_bed_path, sep="\t", header=None, comment="#")
    if pos_df.shape[1] < 9:
        raise ValueError(f"TSS BED must have at least 9 columns: {tss_bed_path}")

    positions = pd.to_numeric(pos_df.iloc[:, 8], errors="coerce")
    values = positions.to_numpy(dtype=float)
    if (
        not np.isfinite(values).all()
        or not np.equal(values, np.floor(values)).all()
        or (values < np.iinfo(np.int32).min).any()
        or (values > np.iinfo(np.int32).max).any()
    ):
        raise ValueError(f"TSS BED offset must be a finite 32-bit integer: {tss_bed_path}")

    region_keys = (
        pos_df.iloc[:, 0].astype(str)
        + ":"
        + pos_df.iloc[:, 1].astype(str)
        + "-"
        + pos_df.iloc[:, 2].astype(str)
    )
    if region_keys.duplicated().any():
        raise ValueError(f"TSS BED contains duplicate genomic intervals: {tss_bed_path}")
    return pd.Series(values.astype(np.int32), index=region_keys, name="position_coarse")


def _validated_gini_cache(
    cache_path: Path,
    sample_ids: list[str],
    empty_cpg_sample_ids: set[str],
) -> pd.DataFrame | None:
    if not cache_path.is_file():
        return None
    try:
        cached = pd.read_csv(cache_path, float_precision="round_trip")
        if tuple(cached.columns) != GINI_CACHE_COLUMNS:
            return None
        cached = _require_exact_sample_ids(cached, sample_ids, "Gini cache")
        values = pd.to_numeric(cached["Gini_Index"], errors="raise").to_numpy(
            dtype=float
        )
        if (
            np.isinf(values).any()
            or (values[np.isfinite(values)] < 0).any()
            or (values[np.isfinite(values)] > 1).any()
        ):
            return None
        empty_mask = cached["Sample_ID"].astype(str).isin(empty_cpg_sample_ids)
        if not np.isnan(values[empty_mask]).all() or not np.isfinite(
            values[~empty_mask]
        ).all():
            return None
        cached["Gini_Index"] = values
        return cached
    except Exception:
        return None


def _validated_tss_cache(
    out_csv: Path,
    sample_ids: set[str],
    empty_cpg_sample_ids: set[str],
) -> pd.DataFrame | None:
    if not out_csv.is_file():
        return None
    try:
        prof = pd.read_csv(
            out_csv, converters=dict.fromkeys(("Sample_ID", "CloneID"), str),
            float_precision="round_trip",
        )
        if tuple(prof.columns) != TSS_REQUIRED_COLUMNS:
            return None
        if prof.empty:
            return prof if sample_ids else None
        sample_values = prof["Sample_ID"].astype(str).str.strip()
        observed_samples = set(sample_values)
        if (
            sample_values.eq("").any()
            or observed_samples - sample_ids
            or observed_samples & empty_cpg_sample_ids
        ):
            return None
        clone_values = prof["CloneID"].astype("string").str.strip()
        if clone_values.isna().any() or clone_values.eq("").any():
            return None
        positions = pd.to_numeric(
            prof["position_coarse"], errors="raise"
        ).to_numpy(dtype=float)
        if (
            not np.isfinite(positions).all()
            or not np.equal(positions, np.floor(positions)).all()
        ):
            return None
        key_frame = pd.DataFrame(
            {"Sample_ID": sample_values, "position_coarse": positions}
        )
        if key_frame.duplicated(["Sample_ID", "position_coarse"]).any():
            return None
        fractions = pd.to_numeric(prof["meth_frac"], errors="raise").to_numpy(
            dtype=float
        )
        counts = pd.to_numeric(prof["count"], errors="raise").to_numpy(dtype=float)
        if (
            not np.isfinite(fractions).all()
            or (fractions < 0).any()
            or (fractions > 1).any()
            or not np.isfinite(counts).all()
            or (counts < 0).any()
        ):
            return None
        prof["Sample_ID"] = sample_values
        prof["position_coarse"] = positions.astype(np.int64)
        prof["meth_frac"] = fractions
        prof["count"] = counts
        prof["CloneID"] = clone_values
        return prof
    except Exception:
        return None


def compute_gini_index(
    raw_adata: ad.AnnData,
    bin_size: int,
    cache_path: Path,
    recompute: bool,
    empty_cpg_sample_ids: set[str],
) -> pd.Series:
    expected_sample_ids = [str(value) for value in raw_adata.obs_names]
    if (not recompute) and cache_path.exists():
        cached = _validated_gini_cache(
            cache_path, expected_sample_ids, empty_cpg_sample_ids
        )
        if cached is not None:
            _log("命中 Gini Index 缓存", level="DEBUG")
            return pd.Series(
                cached["Gini_Index"].to_numpy(dtype=float),
                index=cached["Sample_ID"].astype(str),
                name="Gini_Index",
            )
        _log("Gini Index 缓存不符合当前样本/数值契约，重新计算", level="WARNING")

    _log(f"开始计算 Gini Index (Bin Size: {bin_size})", level="INFO")
    temporary = cache_path.with_name(
        f".gini-tiles.{os.getpid()}.{time.time_ns()}.h5ad"
    )
    adata_bin = None
    try:
        adata_bin = snap.pp.add_tile_matrix(
            raw_adata,
            bin_size=bin_size,
            value_type="total",
            summary_type="sum",
            exclude_chroms=[],
            inplace=False,
            file=str(temporary),
            backend="hdf5",
            chunk_size=4,
        )
        out = np.full(adata_bin.n_obs, np.nan, dtype=np.float64)
        for matrix, start, end in tqdm(
            adata_bin.chunked_X(4),
            total=(adata_bin.n_obs + 3) // 4,
            desc="Gini Index",
        ):
            matrix = matrix.tocsr(copy=False)
            for local_index in range(end - start):
                lo = matrix.indptr[local_index]
                hi = matrix.indptr[local_index + 1]
                values = matrix.data[lo:hi]
                if values.size:
                    out[start + local_index] = _gini_from_sorted_positive(
                        np.sort(values.astype(np.float64, copy=False))
                    )
        s = pd.Series(
            out, index=np.asarray(adata_bin.obs_names, dtype=str), name="Gini_Index"
        )
        result = _require_exact_sample_ids(
            pd.DataFrame({"Sample_ID": s.index, "Gini_Index": s.values}),
            expected_sample_ids,
            "Computed Gini statistics",
        )
        values = pd.to_numeric(result["Gini_Index"], errors="raise").to_numpy(
            dtype=float
        )
        if (
            np.isinf(values).any()
            or (values[np.isfinite(values)] < 0).any()
            or (values[np.isfinite(values)] > 1).any()
        ):
            raise ValueError(
                "Computed Gini values must be within [0, 1] when defined"
            )
        empty_mask = result["Sample_ID"].astype(str).isin(empty_cpg_sample_ids)
        if not np.isnan(values[empty_mask]).all() or not np.isfinite(
            values[~empty_mask]
        ).all():
            raise ValueError("Computed Gini missingness differs from empty-CpG cells")
        _atomic_write_csv(result, cache_path)
        return pd.Series(
            values,
            index=result["Sample_ID"].astype(str),
            name="Gini_Index",
        )
    finally:
        if adata_bin is not None:
            adata_bin.close()
        temporary.unlink(missing_ok=True)


def compute_tss_profile(
    raw_adata: ad.AnnData,
    df_qc: pd.DataFrame,
    tss_bed_path: Path,
    out_csv: Path,
    recompute: bool,
    empty_cpg_sample_ids: set[str],
) -> pd.DataFrame:
    expected_sample_ids = set(df_qc["Sample_ID"].astype(str))
    empty_cpg = {str(value) for value in empty_cpg_sample_ids}
    if empty_cpg - expected_sample_ids:
        raise ValueError("Empty-CpG Sample_ID set differs from QC table")
    if (not recompute) and out_csv.exists():
        prof = _validated_tss_cache(out_csv, expected_sample_ids, empty_cpg)
        if prof is not None:
            _log("命中 TSS Profile 缓存", level="DEBUG")
            return prof
        _log("TSS Profile 缓存不符合当前键/数值契约，重新计算", level="WARNING")

    _log("开始计算 TSS Profile", level="INFO")
    if expected_sample_ids == empty_cpg:
        prof = pd.DataFrame(columns=TSS_REQUIRED_COLUMNS)
        _atomic_write_csv(prof, out_csv)
        validated = _validated_tss_cache(out_csv, expected_sample_ids, empty_cpg)
        if validated is None:
            out_csv.unlink(missing_ok=True)
            raise ValueError("Computed empty TSS profile violates the output contract")
        return validated

    pos_series = _read_tss_offset_map(tss_bed_path)

    def make_and_aggregate(
        value_type: str, summary_type: str, value_name: str, aggregate: str
    ) -> pl.DataFrame:
        temporary = out_csv.with_name(
            f".tss-{value_type}.{os.getpid()}.{time.time_ns()}.h5ad"
        )
        matrix = None
        try:
            matrix = snap.pp.make_peak_matrix(
                raw_adata,
                value_type=value_type,
                peak_file=str(tss_bed_path),
                summary_type=summary_type,
                inplace=False,
                file=str(temporary),
                backend="hdf5",
                chunk_size=4,
            )
            return _aggregate_tss_matrix(
                matrix, pos_series, value_name, aggregate, chunk_size=4
            )
        finally:
            if matrix is not None:
                matrix.close()
            temporary.unlink(missing_ok=True)

    frac_pl = make_and_aggregate("fraction", "mean", "meth_frac", "mean")
    count_pl = make_and_aggregate("target", "sum", "count", "sum")

    prof = frac_pl.join(
        count_pl, on=["Sample_ID", "position_coarse"], how="full", coalesce=True
    ).to_pandas()
    prof = prof.merge(
        df_qc[["Sample_ID", "CloneID"]],
        on="Sample_ID",
        how="left",
        validate="many_to_one",
    )
    _atomic_write_csv(prof, out_csv)
    prof = _validated_tss_cache(out_csv, expected_sample_ids, empty_cpg)
    if prof is None:
        out_csv.unlink(missing_ok=True)
        raise ValueError("Computed TSS profile violates the output contract")

    return prof


def _merge_pipeline_qc(
    pipeline_qc: pd.DataFrame,
    df_qc: pd.DataFrame,
) -> pd.DataFrame:
    source_qc = pipeline_qc.copy()
    required_columns = {
        "sample_id",
        "demux_clone_id",
        "methylation_backend",
        "mapping_policy",
        "native_mapping_policy",
        "native_mapping_unit",
        "native_mapping_pct",
        "dedup_policy",
        "cutadapt_input_pairs",
        "trimmed_pairs",
        "both_primary_mapped_pairs",
        "one_primary_mapped_pairs",
        "backend_accepted_pairs",
        "backend_rejected_pairs",
        "duplicate_pairs",
        "postdedup_pairs",
        "high_cph_assessed_pairs",
        "high_cph_flagged_pairs",
        "high_cph_removed_pairs",
        "final_retained_pairs",
        "non_cpg_methylation_pct",
        "non_cpg_metric_source",
    }
    missing = sorted(required_columns - set(source_qc.columns))
    if missing:
        raise ValueError(
            "Pipeline QC records are missing required field(s): "
            + ", ".join(missing)
        )
    source_qc["sample_id"] = source_qc["sample_id"].astype(str).str.strip()
    if source_qc["sample_id"].eq("").any() or source_qc["sample_id"].duplicated().any():
        raise ValueError("Pipeline QC records contain blank or duplicate sample_id")
    clone_ids = source_qc["demux_clone_id"].astype("string").str.strip()
    if clone_ids.isna().any() or clone_ids.eq("").any():
        raise ValueError("Pipeline QC records contain an empty demux_clone_id")

    count_cols = [
        "cutadapt_input_pairs",
        "trimmed_pairs",
        "both_primary_mapped_pairs",
        "one_primary_mapped_pairs",
        "backend_accepted_pairs",
        "backend_rejected_pairs",
        "duplicate_pairs",
        "postdedup_pairs",
        "high_cph_assessed_pairs",
        "high_cph_flagged_pairs",
        "high_cph_removed_pairs",
        "final_retained_pairs",
    ]
    for column in count_cols:
        values = pd.to_numeric(source_qc[column], errors="raise").to_numpy(dtype=float)
        if (
            not np.isfinite(values).all()
            or (values < 0).any()
            or not np.equal(values, np.floor(values)).all()
        ):
            raise ValueError(f"Pipeline QC records contain invalid {column}")
        source_qc[column] = values
    if (source_qc["cutadapt_input_pairs"] <= 0).any():
        raise ValueError("Pipeline QC records contain zero input read pairs")

    cph_pct = pd.to_numeric(
        source_qc["non_cpg_methylation_pct"], errors="raise"
    ).to_numpy(dtype=float)
    if np.isinf(cph_pct).any() or (cph_pct[np.isfinite(cph_pct)] < 0).any() or (
        cph_pct[np.isfinite(cph_pct)] > 100
    ).any():
        raise ValueError("Pipeline QC records contain invalid CpH percentage")
    source_qc["non_cpg_methylation_pct"] = cph_pct

    def percentage(numerator: str, denominator: str) -> pd.Series:
        denom = source_qc[denominator].replace(0, np.nan)
        return (100.0 * source_qc[numerator] / denom).fillna(0.0)

    source_qc["Trim_Retention%"] = percentage(
        "trimmed_pairs", "cutadapt_input_pairs"
    )
    source_qc["Native_Mapping%"] = pd.to_numeric(
        source_qc["native_mapping_pct"], errors="raise"
    )
    source_qc["Duplicate_Pair_Rate%"] = percentage(
        "duplicate_pairs", "backend_accepted_pairs"
    )
    source_qc["Final_Pair_Yield%"] = percentage(
        "final_retained_pairs", "cutadapt_input_pairs"
    )
    source_qc["High_CpH_Flag_Rate%"] = percentage(
        "high_cph_flagged_pairs", "high_cph_assessed_pairs"
    )

    source_qc.loc[source_qc["methylation_backend"] == "rastair", "High_CpH_Flag_Rate%"] = float("nan")

    column_mapping = {
        "demux_clone_id": "_manifest_clone_id",
        "demux_plate_id": "PlateID",
        "non_cpg_methylation_pct": "Non_CpG_Methylation%",
        "Trim_Retention%": "Trim_Retention%",
        "Native_Mapping%": "Native_Mapping%",
        "Duplicate_Pair_Rate%": "Duplicate_Pair_Rate%",
        "Final_Pair_Yield%": "Final_Pair_Yield%",
        "High_CpH_Flag_Rate%": "High_CpH_Flag_Rate%",
    }

    keep_cols = ["sample_id"] + [c for c in column_mapping if c in source_qc.columns]
    qc_subset = source_qc[keep_cols].rename(
        columns={"sample_id": "Sample_ID", **column_mapping}
    )

    if df_qc is None or df_qc.empty:
        raise ValueError("Current CpG mapping statistics are empty")
    expected_ids = qc_subset["Sample_ID"].astype(str).tolist()
    out = _require_exact_sample_ids(
        df_qc, expected_ids, "CpG mapping statistics"
    ).merge(
        qc_subset,
        on="Sample_ID",
        how="left",
        validate="one_to_one",
    )

    out["Sample_ID"] = out["Sample_ID"].astype(str).str.strip()
    out["CloneID"] = out.pop("_manifest_clone_id").astype("string").str.strip()
    return out


def _run_pipeline_impl(
    paths: RunPaths,
    cfg: RunConfig,
    workers: int,
    recompute_raw_adata: bool,
    recompute_methyl_stats: bool,
    recompute_composition: bool,
    recompute_gini: bool,
    recompute_tss: bool,
    *,
    sealed_context: SealedQCContext | None = None,
    _raw_adata_owner: list,
) -> None:

    timer = Timer()
    overall_t0 = time.perf_counter()

    timer.start("current_outputs")
    _log("Step 1/6: validate sealed Pipeline QC inputs")
    sealed = sealed_context or load_sealed_qc_context(
        paths.project_dir, expected_pipeline_root=_pipeline_dir()
    )
    records = build_qc_records(paths.project_dir, sealed_context=sealed)
    extra_identity_files = {"reference_fai": cfg.ref_fai}
    manifest_sample_ids = [str(record["sample_id"]) for record in records]
    cpg_files = [Path(record.pop("cpg_file")) for record in records]
    cpg_output_identities = {
        str(record["sample_id"]): record.pop("cpg_output_identity") for record in records
    }
    pipeline_qc = pd.DataFrame.from_records(records)
    zero_mapped_sample_ids = {
        str(record["sample_id"])
        for record in records
        if int(record["backend_accepted_pairs"]) == 0
    }
    empty_cpg_sample_ids = {
        str(record["sample_id"])
        for record, path in zip(records, cpg_files, strict=True)
        if _cpg_table_is_empty(path)
    }
    if zero_mapped_sample_ids - empty_cpg_sample_ids:
        raise ValueError("Zero-mapped cells must have an empty current CpG result")
    dt = timer.stop("current_outputs")
    if dt is not None:
        _log(f"  done: {dt:.2f}s")

    timer.start("adata")
    _log("Step 2/6: CpG -> AnnData")
    raw_adata = None
    if recompute_raw_adata:
        _log("Forcing RawAdata rebuild (--recompute-raw-adata)")
    if not paths.cpg_dir.is_dir():
        raise FileNotFoundError(f"Sealed CpG directory is missing: {paths.cpg_dir}")
    with tempfile.TemporaryDirectory(prefix="dna-pipeline-cpg-links-") as temporary:
        flat_dir, linkname2clean = prepare_flat_cpg_symlink_dir(
            Path(temporary), cpg_files
        )
        raw_adata = load_or_create_raw_adata(
            flat_dir,
            paths.raw_adata,
            cfg.ref_fai,
            linkname2clean,
            reuse_existing=not recompute_raw_adata,
            expected_sample_ids=manifest_sample_ids,
            empty_cpg_sample_ids=empty_cpg_sample_ids,
            sealed_delivery_id=sealed.delivery_id,
            sealed_snapshot_id=sealed.snapshot_id,
            cpg_output_identities=cpg_output_identities,
        )
    if raw_adata is not None:
        _raw_adata_owner.append(raw_adata)

    raw_adata_mode = _raw_adata_schema_mode(paths.raw_adata)
    identity = source_identity(
        paths.project_dir,
        [
            {**record, "cpg_file": path}
            for record, path in zip(records, cpg_files, strict=True)
        ],
        sealed_context=sealed,
        extra_files=extra_identity_files,
        context={"raw_adata_mode": raw_adata_mode},
    )
    if raw_adata is None or raw_adata_mode != RAW_ADATA_SCHEMA_MODE:
        raise RuntimeError(
            "RawAdata is unavailable or does not use the current SnapATAC2 0-based four-column CpG "
            f"schema (observed: {raw_adata_mode})"
        )
    identity["cache_identities"] = _cache_identities(identity, cfg, _algorithm_fingerprints())
    previous = {}
    if paths.input_identity.is_file():
        try:
            previous = json.loads(paths.input_identity.read_text(encoding="utf-8"))["cache_identities"]
        except (OSError, ValueError, KeyError, TypeError):
            pass
    if not isinstance(previous, dict):
        previous = {}
    changed = {key for key, value in identity["cache_identities"].items() if previous.get(key) != value}
    if changed:
        _log("Downstream cache identities changed: " + ", ".join(sorted(changed)))
    recompute_methyl_stats |= "cpg" in changed
    recompute_composition |= "cpg" in changed
    recompute_gini |= "gini" in changed
    recompute_tss |= "tss" in changed
    dt = timer.stop("adata")
    if dt is not None:
        _log(f"  done: {dt:.2f}s")

    timer.start("methyl_stats")
    _log("Step 3/6: CpG stream -> mapping stats + signal composition")
    df_qc = pd.DataFrame()
    composition_stats: pd.DataFrame | None = None
    sample_ids = manifest_sample_ids
    cached_mapping = (
        None
        if recompute_methyl_stats
        else _load_methyl_stats_cache(
            paths.methyl_stats_cache, sample_ids, empty_cpg_sample_ids
        )
    )
    cached_composition = (
        None
        if recompute_composition
        else _validated_composition_cache(
            paths.composition_cache,
            sample_ids,
            empty_cpg_sample_ids,
        )
    )
    cached_cpg_results = _select_cpg_cache_state(
        cached_mapping, cached_composition
    )
    if cached_cpg_results is not None:
        df_qc = cached_cpg_results
        _log(f"  use caches: {paths.methyl_stats_cache.name}, {paths.composition_cache.name}")
    else:
        if len(records) != len(cpg_files):
            raise ValueError("Active manifest and CpG input count differ")
        cpg_by_sample: dict[str, Path] = {}
        for record, path in zip(records, cpg_files, strict=True):
            sample_id = str(record["sample_id"])
            if sample_id in cpg_by_sample:
                raise ValueError(f"Duplicate active Sample_ID: {sample_id}")
            if not path.is_file():
                raise FileNotFoundError(f"CpG input does not exist: {path}")
            cpg_by_sample[sample_id] = path
        if set(cpg_by_sample) != set(sample_ids):
            raise ValueError("Active manifest Sample_ID set differs from CpG inputs")
        df_qc, composition_stats = compute_methylation_statistics(
            [cpg_by_sample[sid] for sid in sample_ids],
            sample_ids,
            allowed_chroms=frozenset(_read_chrom_sizes(cfg.ref_fai)),
            workers=(workers if workers and workers > 0 else (os.cpu_count() or 1)),
            return_composition=True,
            empty_cpg_sample_ids=empty_cpg_sample_ids,
        )
        _atomic_write_csv(df_qc, paths.methyl_stats_cache)

    dt = timer.stop("methyl_stats")
    if dt is not None:
        _log(f"  done: {dt:.2f}s")

    timer.start("merge")
    _log("Step 4/6: merge pipeline_qc + df_qc")
    if not pipeline_qc.empty:
        df_qc = _merge_pipeline_qc(pipeline_qc, df_qc)

    dt = timer.stop("merge")
    if dt is not None:
        _log(f"  done: {dt:.2f}s")

    timer.start("composition")
    _log("Step 5/6: CpG signal composition (reuse Step 3 stream + cache)")
    if not df_qc.empty and composition_stats is not None:
        df_qc = compute_cpg_signal_composition(
            df_qc=df_qc,
            cache_path=paths.composition_cache,
            composition_stats=composition_stats,
            empty_cpg_sample_ids=empty_cpg_sample_ids,
        )
    dt = timer.stop("composition")
    if dt is not None:
        _log(f"  done: {dt:.2f}s")

    timer.start("gini_tss_density")
    _log("Step 6/6: gini + tss + density")
    tss_profile = None
    if cfg.tss_bed is None:
        if paths.tss_info.exists():
            paths.tss_info.unlink()
            _log(f"Removed stale skipped-TSS output: {paths.tss_info.name}")

    if raw_adata is not None and not df_qc.empty:
        gini = compute_gini_index(
            raw_adata,
            bin_size=cfg.gini_bin_size,
            cache_path=paths.gini_cache,
            recompute=recompute_gini,
            empty_cpg_sample_ids=empty_cpg_sample_ids,
        )
        df_qc["Gini_Index"] = df_qc["Sample_ID"].astype(str).map(gini)
        gini_missing = set(
            df_qc.loc[df_qc["Gini_Index"].isna(), "Sample_ID"].astype(str)
        )
        if gini_missing != empty_cpg_sample_ids:
            raise ValueError("Gini missingness differs from empty-CpG cells")

        if cfg.tss_bed is not None:
            tss_profile = compute_tss_profile(
                raw_adata=raw_adata,
                df_qc=df_qc,
                tss_bed_path=cfg.tss_bed,
                out_csv=paths.tss_info,
                recompute=recompute_tss,
                empty_cpg_sample_ids=empty_cpg_sample_ids,
            )
        else:
            _log("TSS resources intentionally skipped; TSS CSV/PDF will not be generated", "WARNING")

    if cfg.tss_bed is None or (tss_profile is not None and tss_profile.empty):
        plots_dir = _project_directory(
            paths.project_dir,
            paths.qc_results_dir.relative_to(paths.project_dir) / "plots",
            create=False,
        )
        (plots_dir / "TSS_Profile.pdf").unlink(missing_ok=True)

    dt = timer.stop("gini_tss_density")
    if dt is not None:
        _log(f"  done: {dt:.2f}s")

    timer.start("export")
    _log("Export: write CSVs")
    if not df_qc.empty:
        computed_columns = [
            "Sample_ID",
            "CloneID",
            "Native_Mapping%",
            "Duplicate_Pair_Rate%",
            "Final_Pair_Yield%",
            "Gini_Index",
        ]
        missing_computed_columns = [
            column for column in computed_columns if column not in df_qc.columns
        ]
        if missing_computed_columns:
            raise ValueError(
                "Final QC computation is missing required field(s): "
                + ", ".join(missing_computed_columns)
            )
        sample_id_values = df_qc["Sample_ID"].astype(str)
        empty_cpg_mask = sample_id_values.isin(empty_cpg_sample_ids)
        for column in computed_columns[2:]:
            values = pd.to_numeric(df_qc[column], errors="raise").to_numpy(dtype=float)
            nan_mask = (
                empty_cpg_mask if column == "Gini_Index" else None
            )
            if (
                np.isinf(values).any()
                or (values[np.isfinite(values)] < 0).any()
                or (nan_mask is None and not np.isfinite(values).all())
                or (
                    nan_mask is not None
                    and (
                        not np.isnan(values[nan_mask]).all()
                        or not np.isfinite(values[~nan_mask]).all()
                    )
                )
            ):
                raise ValueError(f"Final QC output contains invalid {column}")
        clone_ids = df_qc["CloneID"].astype("string").str.strip()
        if clone_ids.isna().any() or clone_ids.eq("").any():
            raise ValueError("Final QC output contains an empty CloneID")

        missing_qc_info_columns = [
            column for column in QC_INFO_COLUMNS if column not in df_qc.columns
        ]
        if missing_qc_info_columns:
            raise ValueError(
                "DNAme QC information is missing column(s): "
                + ", ".join(missing_qc_info_columns)
            )
        df_qc = df_qc[list(QC_INFO_COLUMNS)]

        df_qc = _require_exact_sample_ids(
            df_qc, manifest_sample_ids, "Final QC output"
        )
        _atomic_write_csv(df_qc, paths.qc_info)
        _atomic_write_json(identity, paths.input_identity)
    dt = timer.stop("export")
    if dt is not None:
        _log(f"  done: {dt:.2f}s")

    _log(f"Total elapsed: {time.perf_counter() - overall_t0:.2f}s")


def run_pipeline(
    paths: RunPaths,
    cfg: RunConfig,
    workers: int,
    recompute_raw_adata: bool,
    recompute_methyl_stats: bool,
    recompute_composition: bool,
    recompute_gini: bool,
    recompute_tss: bool,
    *,
    sealed_context: SealedQCContext | None = None,
) -> None:
    """运行 QC 流程，并保证所有退出路径都关闭自有的 RawAdata 句柄。"""
    raw_adata_owner: list = []
    try:
        _run_pipeline_impl(
            paths,
            cfg,
            workers,
            recompute_raw_adata,
            recompute_methyl_stats,
            recompute_composition,
            recompute_gini,
            recompute_tss,
            sealed_context=sealed_context,
            _raw_adata_owner=raw_adata_owner,
        )
    finally:
        close_error: BaseException | None = None
        for raw_adata in reversed(raw_adata_owner):
            try:
                raw_adata.close()
            except BaseException as exc:
                if close_error is None:
                    close_error = exc
        if close_error is not None and sys.exc_info()[0] is None:
            raise close_error


def parse_args() -> argparse.Namespace:
    return _add_cli_args(argparse.ArgumentParser(description=__doc__)).parse_args()


def _log_config_summary(paths: RunPaths, cfg: RunConfig):
    _log("=" * 60)
    _log("Running configuration")
    _log(f"项目目录 (Project):    {paths.project_dir}")
    _log(f"CpG 输入 (CpG Dir):     {paths.cpg_dir}")
    _log(f"QC 输出 (QC Results):   {paths.qc_results_dir}")
    _log(f"参考基因组 (Genome):   {cfg.genome}")
    _log("-" * 60)
    _log("关键路径检查:")
    _log(
        f"TSS BED 路径: {cfg.tss_bed} "
        f"{'✅' if cfg.tss_bed is not None and cfg.tss_bed.exists() else 'SKIP'}"
    )
    _log(f"Ref FAI 路径: {cfg.ref_fai} {'✅' if cfg.ref_fai.exists() else '❌'}")
    _log("=" * 60)


def main() -> None:
    args = parse_args()
    project_dir = _resolve_project_dir(Path(args.project_dir))
    sealed = load_sealed_qc_context(
        project_dir, expected_pipeline_root=_pipeline_dir()
    )
    paths = build_paths(project_dir, sealed.delivery_id)
    cfg = build_config(
        sealed,
        tss_bed_override=Path(args.tss_bed) if args.tss_bed else None,
        skip_tss=bool(args.skip_tss),
    )

    _log_config_summary(paths, cfg)
    run_pipeline(
        paths=paths,
        cfg=cfg,
        workers=args.workers,
        recompute_raw_adata=bool(args.recompute_raw_adata),
        recompute_methyl_stats=bool(args.recompute_methyl_stats),
        recompute_composition=bool(args.recompute_composition),
        recompute_gini=bool(args.recompute_gini),
        recompute_tss=bool(args.recompute_tss),
        sealed_context=sealed,
    )


if __name__ == "__main__":
    main()


import logging
import math


matplotlib = mcolors = plt = ticker = sns = Line2D = PCA = StandardScaler = None
_PLOT_IMPORT_ERROR = None


def _load_plot_dependencies() -> None:
    global matplotlib, mcolors, plt, ticker, sns, Line2D, PCA, StandardScaler
    global _PLOT_IMPORT_ERROR

    if plt is not None:
        return
    if _PLOT_IMPORT_ERROR is not None:
        raise RuntimeError("Downstream plotting dependencies are unavailable") from _PLOT_IMPORT_ERROR

    logging.getLogger("matplotlib.font_manager").setLevel(logging.WARNING)
    logging.getLogger("fontTools").setLevel(logging.ERROR)
    try:
        import matplotlib as _matplotlib
        import matplotlib.colors as _mcolors
        import matplotlib.pyplot as _plt
        import matplotlib.ticker as _ticker
        import seaborn as _sns
        from matplotlib.lines import Line2D as _Line2D
        from sklearn.decomposition import PCA as _PCA
        from sklearn.preprocessing import StandardScaler as _StandardScaler
    except ImportError as exc:
        _PLOT_IMPORT_ERROR = exc
        raise RuntimeError("Downstream plotting dependencies are unavailable") from exc

    matplotlib = _matplotlib
    mcolors = _mcolors
    plt = _plt
    ticker = _ticker
    sns = _sns
    Line2D = _Line2D
    PCA = _PCA
    StandardScaler = _StandardScaler


class Colors:
    TITLE = "#414956"
    SUBTITLE = "#525B69"
    LABEL = "#646D7A"
    TICK = "#7A8390"
    ANNOT = "#8B93A0"
    SPINE = "#B8C0CB"
    REFLINE = "#D5DBE3"
    MEDIAN = "#CF718C"

NATURE_BASE_PALETTE = [
    "#4FA29C", "#CF718C", "#6F90BC", "#C9A45A",
    "#7EA271", "#9A83BE", "#D08B55", "#8B9097",
    "#498A86", "#C15473", "#567BAD", "#B99244",
    "#6D8E61", "#8368AD", "#C1773D", "#787D84",
    "#74B3AF", "#D69BAC", "#96ABC8", "#CFB787",
    "#9CB494", "#B5A7CC", "#D4A784", "#A7AAAF",
]

COMPOSITION_PALETTE = [
    "#6F8197", "#C49772", "#8B9A84", "#A88998", "#D1A2A6"
]

class FontSizes:
    TICK = 13
    LABEL = 14.5
    TITLE = 15.5
    LEGEND = 11.5
    LEGEND_TITLE = 12.5
    PANEL_LETTER = 18

class PlotParams:
    DPI = 250
    SAVE_DPI = 300
    SPINE_WIDTH = 0.95
    SCATTER_SIZE = 40
    SCATTER_ALPHA = 0.7
    LINE_WIDTH = 2.0

def setup_nature_style() -> None:
    """应用全部 DNAme QC 图共用的绘图主题。"""
    _load_plot_dependencies()
    sns.set_theme(
        style="white",
        context="talk",
        font="sans-serif",
        font_scale=1.14,
        rc={
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans", "sans-serif"],
            "figure.dpi": PlotParams.DPI,
            "savefig.dpi": PlotParams.SAVE_DPI,
            "axes.linewidth": PlotParams.SPINE_WIDTH,
            "xtick.major.width": 0.9,
            "ytick.major.width": 0.9,
            "xtick.major.size": 4.0,
            "ytick.major.size": 4.0,
            "xtick.labelsize": FontSizes.TICK,
            "ytick.labelsize": FontSizes.TICK,
            "axes.labelsize": FontSizes.LABEL,
            "axes.titlesize": FontSizes.TITLE,
            "legend.fontsize": FontSizes.LEGEND,
            "legend.title_fontsize": FontSizes.LEGEND_TITLE,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        },
    )

def format_ax(
    ax: matplotlib.axes.Axes,
    title: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
) -> None:
    if title:
        ax.set_title(title, fontsize=FontSizes.TITLE, pad=15, weight="bold", color=Colors.SUBTITLE)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=FontSizes.LABEL, labelpad=9, color=Colors.LABEL)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=FontSizes.LABEL, labelpad=9, color=Colors.LABEL)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(Colors.SPINE)
    ax.spines["bottom"].set_color(Colors.SPINE)
    ax.spines["left"].set_linewidth(PlotParams.SPINE_WIDTH)
    ax.spines["bottom"].set_linewidth(PlotParams.SPINE_WIDTH)
    ax.grid(False)
    ax.tick_params(axis="both", which="major", labelsize=FontSizes.TICK, colors=Colors.TICK)

def add_panel_letter(ax: matplotlib.axes.Axes, letter: str) -> None:
    ax.text(
        -0.13, 1.065, letter,
        transform=ax.transAxes,
        fontsize=FontSizes.PANEL_LETTER,
        fontweight="bold",
        ha="left",
        va="top",
        color=Colors.SUBTITLE,
        alpha=0.8,
    )


def panel_label(index: int) -> str:
    """生成不限数量的 Excel 式面板标签：A..Z、AA..AZ、BA…"""
    if index < 0:
        raise ValueError("panel index must be non-negative")
    value = index + 1
    chars = []
    while value:
        value, remainder = divmod(value - 1, 26)
        chars.append(chr(ord("A") + remainder))
    return "".join(reversed(chars))

def add_right_legend(
    fig: plt.Figure,
    group_levels: list[str],
    color_map: dict[str, str],
    title: str = "CloneID",
    marker: str = "o",
    linewidth: float = 0.6,
) -> None:
    legend_handles = []
    for grp in group_levels:
        if marker == "o":
            handle = Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                label=str(grp),
                markerfacecolor=color_map.get(grp, "#8B9097"),
                markeredgecolor="white",
                markeredgewidth=linewidth,
                markersize=8.8,
                alpha=0.8,
            )
        else:
            handle = Line2D(
                [0],
                [0],
                color=color_map.get(grp, "#8B9097"),
                linewidth=linewidth,
                alpha=0.8,
                label=str(grp),
            )
        legend_handles.append(handle)

    legend = fig.legend(
        handles=legend_handles,
        labels=[str(g) for g in group_levels],
        title=title,
        loc="center left",
        bbox_to_anchor=(0.925, 0.5),
        frameon=False,
        fontsize=FontSizes.LEGEND,
        title_fontsize=FontSizes.LEGEND_TITLE,
    )
    legend.get_title().set_color(Colors.SUBTITLE)
    for text in legend.get_texts():
        text.set_color(Colors.TICK)

def _normalize_label_value(value: Any) -> str:
    if value is None:
        return ""
    return str(value)


def _common_suffix(labels: list[str]) -> str:
    if not labels:
        return ""
    reversed_labels = [label[::-1] for label in labels]
    reversed_suffix = os.path.commonprefix(reversed_labels)
    return reversed_suffix[::-1]


def get_diff_only_labels(
    labels: list[Any],
    ellipsis: str = "…",
    protect_single_label: bool = True,
) -> tuple[list[str], dict[str, str]]:
    clean_labels = [_normalize_label_value(label) for label in labels]

    if protect_single_label and len(clean_labels) <= 1:
        return clean_labels, {"prefix": "", "suffix": ""}

    common_prefix = os.path.commonprefix(clean_labels)
    common_suffix = _common_suffix(clean_labels)

    min_len = min((len(label) for label in clean_labels), default=0)
    if len(common_prefix) + len(common_suffix) >= min_len:
        common_suffix = ""
        if len(common_prefix) >= min_len:
            common_prefix = ""

    display_labels = []
    used = set()

    for original in clean_labels:
        start = len(common_prefix)
        end = len(original) - len(common_suffix) if common_suffix else len(original)
        diff = original[start:end]
        display = diff if diff else ellipsis


        if display in used:
            display = f"{display} ({original})"

        display_labels.append(display)
        used.add(display)

    return display_labels, {"prefix": common_prefix, "suffix": common_suffix}


def get_diff_only_label_map(
    labels: list[Any],
    ellipsis: str = "…",
) -> tuple[dict[str, str], dict[str, str]]:
    unique_labels = []
    for label in labels:
        clean_label = _normalize_label_value(label)
        if clean_label not in unique_labels:
            unique_labels.append(clean_label)

    display_labels, omitted = get_diff_only_labels(unique_labels, ellipsis=ellipsis)
    return dict(zip(unique_labels, display_labels, strict=True)), omitted


UNKNOWN_COLOR = "#C7CDD6"


def normalize_clone_id_value(value: Any) -> Any:
    if pd.isna(value):
        return np.nan
    clean_value = str(value).strip()
    if clean_value == "":
        return np.nan
    return clean_value


def make_cloneid_diff_labels(
    values: list[Any],
    ellipsis: str = "…",
) -> tuple[dict[str, str], dict[str, str]]:
    clean_values = []
    for value in values:
        clean_value = normalize_clone_id_value(value)
        if isinstance(clean_value, str) and clean_value not in clean_values:
            clean_values.append(clean_value)
    return get_diff_only_label_map(clean_values, ellipsis=ellipsis)


def cloneid_omission_note(omitted: dict[str, str]) -> str:
    parts = []
    if omitted.get("prefix"):
        parts.append(f"prefix '{omitted['prefix']}'")
    if omitted.get("suffix"):
        parts.append(f"suffix '{omitted['suffix']}'")
    return " and ".join(parts)


def add_cloneid_omission_note(
    fig: plt.Figure,
    omitted: dict[str, str],
    *,
    y: float,
) -> None:
    """以统一样式标注 CloneID 显示名中被省略的公共部分。"""

    note = cloneid_omission_note(omitted)
    if note:
        fig.text(
            0.5,
            y,
            f'Note: common CloneID {note} omitted',
            ha='center',
            fontsize=10,
            color=Colors.ANNOT,
            style='italic',
        )


def _panel_omission_note(
    sample_labels_shortened: bool,
    clone_note: str,
    *,
    include_clone_note: bool = True,
) -> str:
    if sample_labels_shortened and clone_note and include_clone_note:
        return (
            f"Note: common CloneID {clone_note} omitted; "
            "Sample_ID labels shortened within each CloneID panel by omitting shared text"
        )
    if sample_labels_shortened:
        return "Note: Sample_ID labels shortened within each CloneID panel by omitting shared text"
    if clone_note and include_clone_note:
        return f"Note: common CloneID {clone_note} omitted"
    return ""


def generate_clone_color_map(unique_clones: list[str]) -> dict[str, str]:
    n_clones = len(unique_clones)
    if n_clones <= 2:
        global_palette = ["#E78AAE", "#63C2C0"]
    else:
        _load_plot_dependencies()
        global_palette = (
            NATURE_BASE_PALETTE[:n_clones]
            if n_clones <= len(NATURE_BASE_PALETTE)
            else sns.color_palette("Set2", n_clones)
        )
    return {unique_clones[i]: global_palette[i] for i in range(n_clones)}


def plot_wgbs_qc_panel(
    df: pd.DataFrame,
    batch_col: str = 'CloneID',
    clone_color_map: dict = None,
    out_path: Path = None,
    figsize: tuple = (16.5, 15.5),
):
    """绘制 WGBS 对照与样本甲基化核心 QC 面板。"""
    if df.empty:
        return
    _load_plot_dependencies()

    plot_df = df.copy()
    if batch_col not in plot_df.columns:
        plot_df[batch_col] = 'All'
    if batch_col == 'CloneID':
        plot_df[batch_col] = plot_df[batch_col].map(normalize_clone_id_value)

    col_map = {
        'Lambda_Meth_CpG_Rate%': 'lambda_rate',
        'pUC19_Meth_CpG_Rate%': 'pUC19_rate',
        'Sample_Meth_CpG_Rate%': 'sample_rate',
        'Lambda_Unique_CpG_Sites': 'lambda_nCG',
        'pUC19_Unique_CpG_Sites': 'pUC19_nCG',
        'Sample_Unique_CpG_Sites': 'sample_nCG',
    }
    rename_dict = {key: value for key, value in col_map.items() if key in plot_df.columns}
    plot_df = plot_df.rename(columns=rename_dict)
    plot_df[batch_col] = plot_df[batch_col].astype('category')
    panel_cfg = [
        ('sample_nCG', 'sample_rate', 'Sample Methylation vs CpG Coverage'),
        ('lambda_rate', 'pUC19_rate', 'Bisulfite Conversion Controls'),
        ('lambda_rate', 'sample_rate', 'Sample Methylation vs Lambda Control'),
        ('pUC19_rate', 'sample_rate', 'Sample Methylation vs pUC19 Control'),
        ('lambda_nCG', 'lambda_rate', 'Lambda Coverage vs Methylation'),
        ('pUC19_nCG', 'pUC19_rate', 'pUC19 Coverage vs Methylation'),
    ]
    available_cols = set(plot_df.columns)
    panel_cfg = [
        (x_col, y_col, title)
        for x_col, y_col, title in panel_cfg
        if x_col in available_cols and y_col in available_cols
    ]
    if not panel_cfg:
        print('Not enough columns to plot WGBS QC panel.')
        return

    n_panels = len(panel_cfg)
    nrows = math.ceil(n_panels / 2)
    group_levels = list(plot_df[batch_col].cat.categories)
    if batch_col == 'CloneID':
        label_map, omitted = make_cloneid_diff_labels(group_levels)
        legend_levels = [label_map.get(g, g) for g in group_levels]
        color_map = (
            {grp: clone_color_map.get(grp, UNKNOWN_COLOR) for grp in group_levels}
            if clone_color_map
            else generate_clone_color_map(group_levels)
        )
        legend_color_map = {label_map.get(g, g): color_map[g] for g in group_levels}
    else:
        label_map, omitted = {}, {'prefix': '', 'suffix': ''}
        legend_levels = group_levels
        color_map = (
            {grp: clone_color_map.get(grp, UNKNOWN_COLOR) for grp in group_levels}
            if clone_color_map
            else generate_clone_color_map(group_levels)
        )
        legend_color_map = color_map

    fig, axes = plt.subplots(nrows, 2, figsize=(figsize[0], 5 * nrows), dpi=PlotParams.DPI)
    if nrows == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    def pretty_label(value):
        labels = {
            'sample_nCG': 'Sample nCpG',
            'sample_rate': 'Sample Rate',
            'lambda_rate': 'Lambda Rate',
            'pUC19_rate': 'pUC19 Rate',
            'lambda_nCG': 'Lambda nCpG',
            'pUC19_nCG': 'pUC19 nCpG',
        }
        return labels.get(value, value.replace('_', ' ').title())

    for i, (x_col, y_col, title) in enumerate(panel_cfg):
        ax = axes[i]
        for group in group_levels:
            subset = plot_df[plot_df[batch_col] == group]
            ax.scatter(
                subset[x_col],
                subset[y_col],
                s=PlotParams.SCATTER_SIZE,
                alpha=PlotParams.SCATTER_ALPHA,
                color=color_map[group],
                edgecolor='white',
                linewidth=0.55,
                zorder=3,
            )
        x_med = plot_df[x_col].median()
        y_med = plot_df[y_col].median()
        ax.axvline(x_med, ls='--', lw=0.95, color=Colors.REFLINE, alpha=0.8, zorder=0)
        ax.axhline(y_med, ls='--', lw=0.95, color=Colors.REFLINE, alpha=0.8, zorder=0)
        ax.text(
            0.02,
            0.975,
            f'{pretty_label(x_col)} Median: {x_med:.3f}\n'
            f'{pretty_label(y_col)} Median: {y_med:.3f}',
            transform=ax.transAxes,
            ha='left',
            va='top',
            fontsize=10.2,
            color=Colors.ANNOT,
        )
        format_ax(
            ax,
            title=title,
            xlabel=pretty_label(x_col),
            ylabel=pretty_label(y_col),
        )
        add_panel_letter(ax, panel_label(i))
        if x_col in ['sample_nCG', 'lambda_nCG'] and (plot_df[x_col] > 0).all():
            if plot_df[x_col].max() / max(plot_df[x_col].min(), 1) > 20:
                ax.set_xscale('log')
        if 'rate' in y_col:
            finite_y = pd.to_numeric(plot_df[y_col], errors='coerce')
            finite_y = finite_y[np.isfinite(finite_y)]
            if not finite_y.empty:
                ymin = finite_y.min()
                ymax = finite_y.max()
                pad = max((ymax - ymin) * 0.1, 0.006)
                ax.set_ylim(ymin - pad, ymax + pad)

    for index in range(len(panel_cfg), len(axes)):
        axes[index].axis('off')
    add_right_legend(fig, legend_levels, legend_color_map, title=batch_col)
    add_cloneid_omission_note(fig, omitted, y=0.015)
    plt.tight_layout(rect=[0, 0.03, 0.885, 0.958], w_pad=2.2, h_pad=2.4)
    if out_path:
        fig.savefig(out_path, bbox_inches='tight')
    plt.show()
    plt.close(fig)


def normalize_clone_id_list(clone_ids: Any) -> list[str]:
    if clone_ids is None:
        return []

    if isinstance(clone_ids, str):
        candidate_values = [clone_ids]
    else:
        try:
            candidate_values = list(clone_ids)
        except TypeError:
            candidate_values = [clone_ids]

    clean_values = []
    for value in candidate_values:
        clean_value = normalize_clone_id_value(value)
        if pd.isna(clean_value):
            continue
        clean_value = str(clean_value)
        if clean_value not in clean_values:
            clean_values.append(clean_value)
    return clean_values


def make_clone_selection_suffix(clone_ids: Any) -> str:
    clone_list = normalize_clone_id_list(clone_ids)
    if not clone_list:
        return 'all'
    if len(clone_list) == 1:
        return clone_list[0].replace('/', '_')
    return f'{len(clone_list)}clones'

def normalize_metric_list(metrics: Any) -> list[str]:
    if metrics is None:
        return []
    if isinstance(metrics, str):
        return [metrics]
    try:
        return [str(metric) for metric in metrics]
    except TypeError:
        return [str(metrics)]


def _format_metric_stat(value) -> str:
    if pd.isna(value):
        return 'NA'
    if value >= 1000000.0:
        return f'{value / 1000000.0:.2f}M'
    if value >= 1000.0:
        return f'{value / 1000.0:.2f}K'
    if value >= 1:
        return f'{value:.2f}'
    return f'{value:.3f}'


def plot_final_qc_dashboard(
    df: pd.DataFrame,
    group_col: str = 'CloneID',
    clone_id: Any = None,
    mode: str = 'distribution',
    metric: str = 'Native_Mapping%',
    clone_color_map: dict = None,
    out_path: Path = None,
    figsize: tuple = (16.5, 12),
    simplify_sample_ids: bool = True,
):
    """按 CloneID 样本或分布模式绘制所选 QC 指标面板。"""
    if df.empty or group_col not in df.columns:
        print(f"Error: column '{group_col}' not found.")
        return
    _load_plot_dependencies()

    plot_df = df.copy()
    if group_col == 'CloneID':
        plot_df[group_col] = plot_df[group_col].map(normalize_clone_id_value)

    plot_df = plot_df.dropna(subset=[group_col]).copy()
    plot_df[group_col] = plot_df[group_col].astype(str)

    all_group_ids = sorted(plot_df[group_col].dropna().astype(str).unique().tolist())

    selected_clone_ids = normalize_clone_id_list(clone_id)

    if selected_clone_ids:
        missing = [cid for cid in selected_clone_ids if cid not in all_group_ids]
        if missing:
            print("Error: some clone_id values were not found:")
            for value in missing:
                print(f"  - {value}")
            print("Available CloneID values:")
            for value in all_group_ids:
                print(f"  - {value}")
            return

        plot_df = plot_df.loc[plot_df[group_col].isin(selected_clone_ids)].copy()
        group_ids = [cid for cid in all_group_ids if cid in selected_clone_ids]
    else:
        group_ids = all_group_ids

    if plot_df.empty:
        print("Error: no data available after CloneID filtering.")
        return

    metric_configs = {
        'Sample_Unique_CpG_Sites': ('Unique CpG Sites', True),
        'Gini_Index': ('Library Uniformity (Gini Index)', False),
        'Non_CpG_Methylation%': ('Non-CpG Methylation (%)', False),
        'Native_Mapping%': ('Native Mapping (%)', False),
        'Final_Pair_Yield%': ('Final Pair Yield (%)', False),
        'Duplicate_Pair_Rate%': ('Duplicate Pair Rate (%)', False),
    }

    if mode == 'distribution':
        label_map, omitted = make_cloneid_diff_labels(group_ids)
        plot_df[f'{group_col}_Display'] = plot_df[group_col].map(label_map).fillna(plot_df[group_col])
        plot_df[f'{group_col}_Display'] = plot_df[f'{group_col}_Display'].astype('category')

        valid_metrics = [m for m in metric_configs if m in plot_df.columns]
        if len(valid_metrics) == 0:
            print('Error: no valid QC metrics found in dataframe.')
            return

        display_levels = [label_map.get(g, g) for g in group_ids]
        if clone_color_map:
            color_map = {
                label_map.get(orig, orig): clone_color_map.get(orig, UNKNOWN_COLOR)
                for orig in group_ids
            }
        else:
            color_map = generate_clone_color_map(display_levels)

        n_cols = 2
        n_rows = math.ceil(len(valid_metrics) / n_cols)
        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=(figsize[0], 5 * n_rows),
            dpi=PlotParams.DPI,
        )
        axes = np.array(axes).reshape(-1)

        for (i, metric_name) in enumerate(valid_metrics):
            ax = axes[i]
            display_name, use_log = metric_configs[metric_name]

            xvals = pd.to_numeric(plot_df[metric_name], errors='coerce')
            panel_df = plot_df.loc[xvals.notna(), [f'{group_col}_Display', metric_name]].copy()
            panel_df[metric_name] = pd.to_numeric(panel_df[metric_name], errors='coerce')
            panel_df = panel_df.dropna(subset=[metric_name])

            if panel_df.empty:
                ax.text(
                    0.5, 0.5, 'No data', transform=ax.transAxes,
                    ha='center', va='center', fontsize=12, color=Colors.ANNOT
                )
                ax.set_axis_off()
                continue

            for display_group in display_levels:
                sub = panel_df[panel_df[f'{group_col}_Display'] == display_group]
                if sub.empty:
                    continue
                vals = sub[metric_name].dropna()
                if len(vals) == 0:
                    continue
                sns.histplot(
                    vals,
                    bins=30,
                    stat='count',
                    element='step',
                    fill=True,
                    alpha=0.18 if not selected_clone_ids or len(selected_clone_ids) > 1 else 0.24,
                    linewidth=1.15 if not selected_clone_ids or len(selected_clone_ids) > 1 else 1.35,
                    color=color_map.get(display_group, UNKNOWN_COLOR),
                    ax=ax,
                    log_scale=(use_log, False),
                )

            median_val = panel_df[metric_name].median()
            ax.axvline(median_val, ls='--', lw=1.0, color=Colors.MEDIAN, alpha=0.9, zorder=5)
            ax.text(
                0.02, 0.975, f'Median: {_format_metric_stat(median_val)}',
                transform=ax.transAxes, ha='left', va='top',
                fontsize=10.2, color=Colors.ANNOT,
            )

            format_ax(
                ax,
                title=display_name,
                xlabel=f"{display_name}{(' (Log Scale)' if use_log else '')}",
                ylabel='Count',
            )
            add_panel_letter(ax, panel_label(i))

            if not use_log:
                ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
                ax.ticklabel_format(style='plain', axis='x')

        for j in range(len(valid_metrics), len(axes)):
            axes[j].axis('off')

        legend_title = group_col if not selected_clone_ids else f'{group_col} selected'
        add_right_legend(fig, display_levels, color_map, title=legend_title, marker='line', linewidth=2.2)

        if selected_clone_ids:
            display_summary = ', '.join(display_levels[:3]) + (' ...' if len(display_levels) > 3 else '')
            fig.suptitle(
                f'Final QC Dashboard - {display_summary}',
                fontsize=FontSizes.TITLE + 3,
                weight='bold',
                color=Colors.TITLE,
                y=0.995,
            )

        if not selected_clone_ids:
            add_cloneid_omission_note(fig, omitted, y=0.015)

        plt.tight_layout(rect=[0, 0.03, 0.885, 0.97], w_pad=2.2, h_pad=2.4)

        if out_path:
            fig.savefig(out_path, bbox_inches='tight')
        plt.show()
        plt.close(fig)
        return

    if mode != 'per_clone_samples':
        print("Error: unsupported mode. Use 'distribution' or 'per_clone_samples'.")
        return

    if 'Sample_ID' not in plot_df.columns:
        print("Error: 'Sample_ID' column is required for mode='per_clone_samples'.")
        return

    if metric not in metric_configs:
        print(f"Error: unsupported metric '{metric}'.")
        print("Supported metrics:")
        for name in metric_configs:
            print(f"  - {name}")
        return

    if metric not in plot_df.columns:
        print(f"Error: metric '{metric}' not found in dataframe.")
        return

    display_name, use_log = metric_configs[metric]
    clones = group_ids
    clone_display_map, omitted = make_cloneid_diff_labels(clones)

    if clone_color_map is None:
        clone_color_map = generate_clone_color_map(clones)

    n_cols = 2
    n_rows = max(1, math.ceil(len(clones) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 4.9 * n_rows), dpi=PlotParams.DPI)
    axes = np.array(axes).reshape(-1)
    sample_label_note_needed = False

    for i, clone in enumerate(clones):
        ax = axes[i]
        sub = plot_df.loc[plot_df[group_col] == clone, ['Sample_ID', metric]].copy()
        sub[metric] = pd.to_numeric(sub[metric], errors='coerce')
        sub = sub.dropna(subset=['Sample_ID', metric])

        if sub.empty:
            ax.axis('off')
            continue

        sub['Sample_ID'] = sub['Sample_ID'].astype(str)

        if simplify_sample_ids:
            sample_label_map, sample_omitted = get_diff_only_label_map(sub['Sample_ID'].tolist())
            sub['Sample_ID_Display'] = sub['Sample_ID'].map(sample_label_map)
            if sample_omitted.get('prefix') or sample_omitted.get('suffix'):
                sample_label_note_needed = True
        else:
            sub['Sample_ID_Display'] = sub['Sample_ID']

        sub = sub.sort_values(metric, ascending=False).reset_index(drop=True)

        bar_color = clone_color_map.get(clone, UNKNOWN_COLOR)
        ax.bar(
            np.arange(len(sub)),
            sub[metric].to_numpy(dtype=float),
            color=bar_color,
            alpha=0.85,
            edgecolor='white',
            linewidth=0.4,
        )

        median_val = sub[metric].median()
        ax.axhline(median_val, ls='--', lw=1.0, color=Colors.MEDIAN, alpha=0.9, zorder=4)
        ax.text(
            0.02, 0.975, f'Median: {_format_metric_stat(median_val)}',
            transform=ax.transAxes, ha='left', va='top',
            fontsize=10.2, color=Colors.ANNOT,
        )

        ax.set_xticks(np.arange(len(sub)))
        n_samples = len(sub)
        if n_samples > 60:
            ax.set_xticklabels([])
        else:
            ax.set_xticklabels(sub['Sample_ID_Display'].tolist(), rotation=90, fontsize=5.8, color=Colors.TICK)

        if use_log and (sub[metric] > 0).all():
            ax.set_yscale('log')
            ylabel = f'{display_name} (Log Scale)'
        else:
            ylabel = display_name

        format_ax(
            ax,
            title=f"{clone_display_map.get(clone, clone)} - {display_name}",
            xlabel='',
            ylabel=ylabel,
        )
        add_panel_letter(ax, panel_label(i))

    for j in range(len(clones), len(axes)):
        axes[j].axis('off')

    note_text = _panel_omission_note(
        sample_label_note_needed,
        cloneid_omission_note(omitted),
        include_clone_note=not selected_clone_ids,
    )

    if note_text:
        fig.text(0.5, 0.01, note_text, ha='center', fontsize=10, color=Colors.ANNOT, style='italic')

    plt.tight_layout(rect=[0, 0.03, 1, 0.98])
    if out_path:
        fig.savefig(out_path, bbox_inches='tight')
    plt.show()
    plt.close(fig)


def plot_plate_metric_map(
    df: pd.DataFrame,
    metric: str = 'Sample_Unique_CpG_Sites',
    clone_color_map: dict = None,
    out_path: Path = None,
):
    """按 CloneID 分面绘制一个 QC 指标的孔板位置热图。"""
    if df.empty or 'PlateID' not in df.columns or 'CloneID' not in df.columns:
        return
    _load_plot_dependencies()

    row_order = list('ABCDEFGH')
    row_map = {letter: index for index, letter in enumerate(row_order)}
    plot_df = df.copy()
    plot_df['CloneID'] = plot_df['CloneID'].map(normalize_clone_id_value)
    unique_for_labels = sorted(plot_df['CloneID'].dropna().unique().tolist())
    clone_display_map, omitted = make_cloneid_diff_labels(unique_for_labels)
    plate_id = plot_df['PlateID'].astype('string').str.strip().str.upper()
    if plate_id.fillna('').eq('').all():
        print('当前细胞无 PlateID，跳过孔板图。')
        return
    valid_plate = plate_id.str.fullmatch(r'[A-H](?:[1-9]|1[0-2])', na=False)
    if (~valid_plate).any():
        invalid = plot_df.loc[~valid_plate]
        examples = []
        for idx, row in invalid.head(10).iterrows():
            sample = str(row.get('Sample_ID', idx))
            value = row.get('PlateID')
            value_label = '<empty>' if pd.isna(value) or not str(value).strip() else str(value)
            examples.append(f'{sample}={value_label}')
        remainder = len(invalid) - len(examples)
        suffix = f'; 另有 {remainder} 行' if remainder > 0 else ''
        print(
            f"⚠️ 忽略 {len(invalid)} 行无效 PlateID（应为 A1-H12）："
            f"{'; '.join(examples)}{suffix}"
        )
    plot_df = plot_df.loc[valid_plate].copy()
    if plot_df.empty:
        print('⚠️ 没有可用于 Plate View 的有效 PlateID')
        return
    plot_df['PlateID'] = plate_id.loc[valid_plate]
    plot_df['row_idx'] = plot_df['PlateID'].str[0].map(row_map)
    plot_df['col_idx'] = plot_df['PlateID'].str[1:].astype(int)
    unique_clones = sorted(plot_df['CloneID'].dropna().unique().tolist())
    if not unique_clones:
        print('⚠️ 没有同时包含有效 PlateID 和 CloneID 的记录；请重新生成 DNAme QC 表')
        return

    num_clones = len(unique_clones)
    n_cols = 2
    n_rows = math.ceil(num_clones / n_cols)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(n_cols * 9, n_rows * 6),
        dpi=PlotParams.DPI,
    )
    if num_clones == 1:
        axes = np.array([axes])
    axes_flat = axes.flatten()

    for index, clone in enumerate(unique_clones):
        ax = axes_flat[index]
        subset = plot_df[plot_df['CloneID'] == clone]
        clone_title = clone_display_map.get(clone, clone)
        ax.set_facecolor('white')
        ax.grid(
            True,
            which='major',
            color='#EAEAEA',
            linestyle='--',
            linewidth=0.8,
            zorder=0,
        )
        if clone_color_map and clone in clone_color_map:
            base_color = clone_color_map[clone]
            cmap = mcolors.LinearSegmentedColormap.from_list(
                f'custom_{clone}', ['#F5F5F5', base_color]
            )
        else:
            cmap = 'YlOrBr'
        ax.scatter(
            subset['col_idx'],
            subset['row_idx'],
            s=850,
            c=subset[metric] if metric in subset.columns else np.zeros(len(subset)),
            cmap=cmap,
            edgecolors='white',
            linewidth=0.8,
            alpha=0.85,
            zorder=2,
        )
        if metric in subset.columns:
            for _, row in subset.iterrows():
                val = row[metric]
                if pd.isna(val):
                    continue
                if val >= 1000000.0:
                    label = f'{int(val / 1000000.0)}M'
                elif val >= 1000.0:
                    label = f'{int(val / 1000.0)}K'
                else:
                    label = f'{int(val)}'
                ax.text(
                    row['col_idx'],
                    row['row_idx'],
                    label,
                    ha='center',
                    va='center',
                    fontsize=9,
                    color=Colors.TITLE,
                    fontweight='normal',
                    zorder=3,
                )
        ax.set_xticks(range(1, 13))
        ax.set_yticks(range(8))
        ax.set_yticklabels(row_order)
        ax.invert_yaxis()
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color(Colors.SPINE)
            spine.set_linewidth(1.5)
        ax.tick_params(
            axis='both',
            which='major',
            color=Colors.TICK,
            labelcolor=Colors.LABEL,
            labelsize=11,
            pad=6,
        )
        metric_title = metric.replace('_', ' ').title()
        ax.set_title(
            f'{clone_title} - {metric_title}',
            fontsize=14,
            pad=15,
            weight='bold',
            color=Colors.SUBTITLE,
        )
    for index in range(num_clones, len(axes_flat)):
        axes_flat[index].axis('off')
    add_cloneid_omission_note(fig, omitted, y=0.01)
    plt.tight_layout(rect=[0, 0.03, 1, 1])
    if out_path:
        fig.savefig(out_path, bbox_inches='tight')
    plt.show()
    plt.close(fig)


def plot_cpg_signal_composition(
    df: pd.DataFrame,
    out_path: Path = None,
    simplify_sample_ids: bool = True,
):
    """绘制 host、lambda、pUC19 与 mtDNA 的 CpG signal 比例。"""
    if df.empty or 'Sample_ID' not in df.columns or 'CloneID' not in df.columns:
        return

    cols = ['Signal_Host_Rate%', 'Signal_Lambda_Rate%', 'Signal_pUC19_Rate%', 'Signal_mtDNA_Rate%']
    valid_cols = [c for c in cols if c in df.columns]
    if not valid_cols:
        return
    _load_plot_dependencies()

    plot_df = df.copy()
    plot_df['CloneID'] = plot_df['CloneID'].map(normalize_clone_id_value)
    plot_df = plot_df.dropna(subset=['CloneID', 'Sample_ID'])

    colors = COMPOSITION_PALETTE[:len(valid_cols)]
    clones = sorted(plot_df['CloneID'].dropna().unique().tolist())
    if len(clones) == 0:
        return

    clone_display_map, omitted = make_cloneid_diff_labels(clones)
    n_cols = 2
    n_rows = math.ceil(len(clones) / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 4.8 * n_rows), dpi=PlotParams.DPI)
    axes = np.array(axes).reshape(-1)

    sample_label_note_needed = False

    for (i, clone) in enumerate(clones):
        ax = axes[i]
        sub = plot_df.loc[plot_df['CloneID'] == clone, ['Sample_ID', *valid_cols]].copy()

        if sub.empty:
            ax.axis('off')
            continue

        sub['Sample_ID'] = sub['Sample_ID'].astype(str)

        if simplify_sample_ids:
            sample_label_map, sample_omitted = get_diff_only_label_map(sub['Sample_ID'].tolist())
            sub['Sample_ID_Display'] = sub['Sample_ID'].map(sample_label_map)
            if sample_omitted.get('prefix') or sample_omitted.get('suffix'):
                sample_label_note_needed = True
        else:
            sub['Sample_ID_Display'] = sub['Sample_ID']

        sub = sub.set_index('Sample_ID_Display')[valid_cols]
        sub.plot(
            kind='bar',
            stacked=True,
            ax=ax,
            width=0.85,
            color=colors,
            alpha=0.85,
            edgecolor='white',
            linewidth=0.4,
            legend=False,
        )

        ax.set_title(
            f'{clone_display_map.get(clone, clone)} - CpG Signal Composition',
            weight='bold',
            color=Colors.SUBTITLE,
            pad=12,
        )
        ax.set_xlabel('')
        ax.set_ylabel('Percentage of CpG Signals (%)', color=Colors.LABEL, labelpad=10)

        n_samples = sub.shape[0]
        if n_samples > 60:
            ax.set_xticklabels([])
        else:
            ax.tick_params(axis='x', labelrotation=90, labelsize=5.8)
            for tick_label in ax.get_xticklabels():
                tick_label.set_color(Colors.TICK)

        ax.tick_params(axis='x', length=0, colors=Colors.TICK, labelcolor=Colors.TICK)
        ax.tick_params(axis='y', colors=Colors.TICK, labelcolor=Colors.TICK)

        sns.despine(ax=ax)
        ax.spines['left'].set_color(Colors.SPINE)
        ax.spines['bottom'].set_color(Colors.SPINE)
        ax.spines['left'].set_linewidth(1.2)
        ax.spines['bottom'].set_linewidth(1.2)

    for j in range(len(clones), len(axes)):
        axes[j].axis('off')

    legend_labels = [c.replace('Reads_', '').replace('_Rate', '').replace('_', ' ') for c in valid_cols]
    handles = [
        plt.Rectangle((0, 0), 1, 1, fc=colors[k], ec='white', lw=0.4, alpha=0.85)
        for k in range(len(valid_cols))
    ]
    legend = fig.legend(
        handles=handles,
        labels=legend_labels,
        loc='center left',
        bbox_to_anchor=(0.92, 0.5),
        frameon=False,
        fontsize=FontSizes.LEGEND,
    )
    for text in legend.get_texts():
        text.set_color(Colors.TICK)

    note_text = _panel_omission_note(
        sample_label_note_needed, cloneid_omission_note(omitted)
    )

    if note_text:
        fig.text(
            0.5,
            0.01,
            note_text,
            ha='center',
            fontsize=10,
            color=Colors.ANNOT,
            style='italic',
        )

    plt.tight_layout(rect=[0, 0.03, 0.9, 0.98])
    if out_path:
        fig.savefig(out_path, bbox_inches='tight')
    plt.show()
    plt.close(fig)


def plot_tss_profile(
    df_tss: pd.DataFrame | str | Path,
    clone_color_map: dict = None,
    out_path: Path = None,
):
    """绘制 TSS 甲基化比例与信号数，保留 CSV 身份和浮点精度；无可绘数据时清理旧图。"""
    if isinstance(df_tss, (str, Path)):
        df_tss = pd.read_csv(
            df_tss, converters=dict.fromkeys(("Sample_ID", "CloneID"), str),
            float_precision="round_trip",
        )
    required = {'position_coarse', 'meth_frac', 'CloneID'}
    if df_tss.empty or not required.issubset(df_tss.columns):
        if out_path is not None:
            Path(out_path).unlink(missing_ok=True)
        return
    _load_plot_dependencies()

    plot_df = df_tss.copy()
    plot_df['CloneID'] = plot_df['CloneID'].map(normalize_clone_id_value)
    group_levels = sorted(plot_df['CloneID'].dropna().astype(str).unique().tolist())
    if not group_levels:
        if out_path is not None:
            Path(out_path).unlink(missing_ok=True)
        print('No valid CloneID values: TSS profile skipped.')
        return
    clone_display_map, omitted = make_cloneid_diff_labels(group_levels)
    if clone_color_map is None:
        clone_color_map = generate_clone_color_map(group_levels)

    num_clones = len(group_levels)
    n_rows = num_clones
    fig, axes = plt.subplots(
        n_rows,
        2,
        figsize=(16.5, max(5.5, n_rows * 5.2)),
        dpi=PlotParams.DPI,
    )
    if num_clones == 1:
        axes = np.array([axes])
    for index, clone in enumerate(group_levels):
        subset = plot_df[plot_df['CloneID'] == clone]
        clone_title = clone_display_map.get(clone, clone)
        ax1, ax2 = axes[index]
        sns.lineplot(
            data=subset,
            x='position_coarse',
            y='meth_frac',
            color=clone_color_map.get(clone, None),
            ax=ax1,
            linewidth=PlotParams.LINE_WIDTH,
            legend=False,
        )
        format_ax(
            ax1,
            title=f'{clone_title} - TSS Methylation Profile',
            xlabel='Distance to TSS (bp)',
            ylabel='mCpG Fraction',
        )
        if 'count' in subset.columns:
            sns.lineplot(
                data=subset,
                x='position_coarse',
                y='count',
                color=clone_color_map.get(clone, None),
                ax=ax2,
                linewidth=PlotParams.LINE_WIDTH,
                legend=False,
            )
            format_ax(
                ax2,
                title=f'{clone_title} - TSS mCpG Signal Profile',
                xlabel='Distance to TSS (bp)',
                ylabel='mCpG Signal Sum',
            )
        else:
            ax2.axis('off')
    add_cloneid_omission_note(fig, omitted, y=0.01)
    plt.tight_layout(rect=[0, 0.03, 1, 1])
    if out_path:
        plt.savefig(out_path, bbox_inches='tight')
    plt.show()
    plt.close(fig)


CPG_SLOPE_RATIO_MIN_CHANGE = 0.05


def _robust_inlier_mask(
    values: np.ndarray,
    z_thresh: float = 3.5,
    min_n: int = 6,
) -> np.ndarray:

    values = np.asarray(values, dtype=float)
    mask = np.isfinite(values)
    if mask.sum() < min_n:
        return mask

    finite_values = values[mask]
    median = float(np.median(finite_values))
    deviations = np.abs(finite_values - median)
    mad = float(np.median(deviations))
    if not np.isfinite(mad):
        return mask

    if mad < 1e-12:


        keep_finite = np.isclose(
            finite_values,
            median,
            rtol=1e-9,
            atol=1e-12,
        )
    else:
        modified_z = 0.6745 * deviations / mad
        keep_finite = modified_z <= z_thresh

    mask[np.flatnonzero(mask)] = keep_finite
    return mask


def _leave_one_out_slope_ratio_inlier_mask(
    x: np.ndarray,
    y_unique: np.ndarray,
    y_signal: np.ndarray,
    z_thresh: float = 3.5,
    min_n: int = 6,
    min_relative_change: float = CPG_SLOPE_RATIO_MIN_CHANGE,
) -> np.ndarray:

    n = len(x)
    if n < min_n:
        return np.ones(n, dtype=bool)
    if min_relative_change < 0:
        raise ValueError('Minimum slope-ratio change must be non-negative')

    xx = x * x
    xy_unique = x * y_unique
    xy_signal = x * y_signal
    total_xx = np.sum(xx)
    total_unique = np.sum(xy_unique)
    total_signal = np.sum(xy_signal)
    remaining_xx = total_xx - xx
    remaining_unique = total_unique - xy_unique
    remaining_signal = total_signal - xy_signal

    with np.errstate(divide='ignore', invalid='ignore'):
        full_slope_unique = total_unique / total_xx
        full_slope_signal = total_signal / total_xx
        full_slope_ratio = full_slope_unique / full_slope_signal
        slope_unique = remaining_unique / remaining_xx
        slope_signal = remaining_signal / remaining_xx
        leave_one_out_ratio = slope_unique / slope_signal
        log_slope_ratio = np.log(leave_one_out_ratio)
        relative_change = np.abs(leave_one_out_ratio / full_slope_ratio - 1.0)


    robust_inlier = _robust_inlier_mask(
        log_slope_ratio,
        z_thresh=z_thresh,
        min_n=min_n,
    )
    material_change = np.isfinite(relative_change) & (
        relative_change >= min_relative_change
    )
    return robust_inlier | ~material_change


def _fit_cpg_slope_pair(
    x: np.ndarray,
    y_unique: np.ndarray,
    y_signal: np.ndarray,
    z_thresh: float = 3.5,
    min_n: int = 6,
    min_ratio_change: float = CPG_SLOPE_RATIO_MIN_CHANGE,
) -> tuple[float, float, np.ndarray]:

    x = np.asarray(x, dtype=float)
    y_unique = np.asarray(y_unique, dtype=float)
    y_signal = np.asarray(y_signal, dtype=float)
    if not (x.shape == y_unique.shape == y_signal.shape):
        raise ValueError('CpG slope arrays must have identical shapes')
    if min_ratio_change < 0:
        raise ValueError('Minimum slope-ratio change must be non-negative')

    valid = (
        np.isfinite(x)
        & np.isfinite(y_unique)
        & np.isfinite(y_signal)
        & (x > 0)
        & (y_unique > 0)
        & (y_signal > 0)
    )
    inlier_mask = np.zeros(len(x), dtype=bool)
    if valid.sum() < 2:
        return np.nan, np.nan, inlier_mask

    valid_positions = np.flatnonzero(valid)
    x_valid = x[valid]
    unique_valid = y_unique[valid]
    signal_valid = y_signal[valid]
    keep = np.ones(len(x_valid), dtype=bool)

    if len(x_valid) >= min_n:
        with np.errstate(divide='ignore', invalid='ignore'):
            log_point_ratio = np.log(unique_valid / signal_valid)

        ratio_outlier = ~_robust_inlier_mask(
            log_point_ratio, z_thresh, min_n
        )
        ratio_influential = ~_leave_one_out_slope_ratio_inlier_mask(
            x_valid,
            unique_valid,
            signal_valid,
            z_thresh=z_thresh,
            min_n=min_n,
            min_relative_change=min_ratio_change,
        )


        confirmed_outlier = ratio_outlier & ratio_influential
        keep &= ~confirmed_outlier


        if keep.sum() < 2:
            keep = np.ones(len(x_valid), dtype=bool)

    inlier_mask[valid_positions] = keep
    x_fit = x[inlier_mask]
    unique_fit = y_unique[inlier_mask]
    signal_fit = y_signal[inlier_mask]
    denominator = float(np.dot(x_fit, x_fit))
    if not np.isfinite(denominator) or denominator <= 0:
        return np.nan, np.nan, inlier_mask

    slope_unique = float(np.dot(x_fit, unique_fit) / denominator)
    slope_signal = float(np.dot(x_fit, signal_fit) / denominator)
    return slope_unique, slope_signal, inlier_mask


def load_cpg_density_frame(
    project_dir: str | Path,
    qc_info_path: Path,
    composition_cache_path: Path,
    *,
    sealed_context: SealedQCContext | None = None,
) -> pd.DataFrame:
    """从原始来源组装 CpG_Density 输入帧，不依赖 20 列 QC 表之外的衍生列。

    pair 计数来自 sealed final manifest（build_qc_records 唯一加载链），
    Sample_CpG_Signals 来自 composition cache，CloneID 与 unique CpG sites
    来自当前 QC 表；三个来源的样本集合必须与 QC 表完全一致。
    """
    qc = pd.read_csv(
        qc_info_path, converters=dict.fromkeys(QC_INFO_COLUMNS[:3], str),
        float_precision="round_trip",
    )
    required = ("Sample_ID", "CloneID", "Sample_Unique_CpG_Sites")
    missing = [column for column in required if column not in qc.columns]
    if missing:
        raise ValueError(
            "CpG density inputs are missing QC column(s): " + ", ".join(missing)
        )
    qc["Sample_ID"] = qc["Sample_ID"].astype(str).str.strip()

    records = build_qc_records(project_dir, sealed_context=sealed_context)
    pairs = pd.DataFrame(
        {
            "Sample_ID": [str(record["sample_id"]) for record in records],
            "Input_Read_Pairs": [
                int(record["cutadapt_input_pairs"]) for record in records
            ],
            "Final_Retained_Pairs": [
                int(record["final_retained_pairs"]) for record in records
            ],
        }
    )

    if not composition_cache_path.is_file():
        raise FileNotFoundError(
            f"CpG composition cache is absent: {composition_cache_path}"
        )
    composition = pd.read_csv(composition_cache_path, float_precision="round_trip")
    for column in ("Sample_ID", "Sample_CpG_Signals"):
        if column not in composition.columns:
            raise ValueError(f"CpG composition cache is missing column: {column}")
    composition["Sample_ID"] = composition["Sample_ID"].astype(str).str.strip()

    qc_ids = qc["Sample_ID"].tolist()
    for source, label in (
        (pairs, "sealed pair counts"),
        (composition, "composition signals"),
    ):
        if set(source["Sample_ID"].astype(str)) != set(qc_ids):
            raise ValueError(f"CpG density {label} do not cover the QC samples")

    frame = (
        qc[list(required)]
        .merge(pairs, on="Sample_ID", how="left", validate="one_to_one")
        .merge(
            composition[["Sample_ID", "Sample_CpG_Signals"]],
            on="Sample_ID",
            how="left",
            validate="one_to_one",
        )
    )
    for column in (
        "Sample_Unique_CpG_Sites",
        "Input_Read_Pairs",
        "Final_Retained_Pairs",
        "Sample_CpG_Signals",
    ):
        frame[column] = pd.to_numeric(frame[column], errors="raise")
    return frame


def plot_cpg_density(
    df: pd.DataFrame,
    clone_color_map: dict = None,
    out_path: Path = None,
    min_slope_ratio_change: float = CPG_SLOPE_RATIO_MIN_CHANGE,
):
    """绘制 CpG 密度，并在异常点实质影响 slope ratio 时共同剔除。"""
    if min_slope_ratio_change < 0:
        raise ValueError('Minimum slope-ratio change must be non-negative')

    metric_a = 'Sample_Unique_CpG_Sites'
    metric_b = 'Sample_CpG_Signals'

    denominator_groups = [
        ('Input_Read_Pairs', 'Input Pair', 'input pair', 'Input-Pair'),
        ('Final_Retained_Pairs', 'Final Retained Pair', 'final retained pair', 'Final-Pair'),
    ]

    if df.empty or 'CloneID' not in df.columns:
        return
    missing_denominators = [
        x_col for x_col, *_ in denominator_groups if x_col not in df.columns
    ]
    if missing_denominators:
        print(f"Missing required denominator for CpG slopes: {', '.join(missing_denominators)}")
        return
    _load_plot_dependencies()

    plot_df = df.copy()
    plot_df['CloneID'] = plot_df['CloneID'].map(normalize_clone_id_value)
    plot_df = plot_df.dropna(subset=['CloneID'])
    clones = sorted(plot_df['CloneID'].dropna().unique().tolist())
    if not clones:
        return
    clone_display_map, omitted = make_cloneid_diff_labels(clones)
    if clone_color_map is None:
        clone_color_map = generate_clone_color_map(clones)
    clone_color_map = {
        normalize_clone_id_value(key): value
        for key, value in clone_color_map.items()
        if isinstance(normalize_clone_id_value(key), str)
    }

    def calc_slopes_by_clone(x_col: str) -> tuple[pd.Series, pd.Series]:
        required = [x_col, metric_a, metric_b]
        if any(column not in plot_df.columns for column in required):
            return pd.Series(dtype=float), pd.Series(dtype=float)

        unique_rows = []
        signal_rows = []
        for clone in clones:
            subset_columns = required + (
                ['Sample_ID'] if 'Sample_ID' in plot_df.columns else []
            )
            subset = plot_df.loc[
                plot_df['CloneID'] == clone,
                subset_columns,
            ].copy()
            x = pd.to_numeric(subset[x_col], errors='coerce').to_numpy(dtype=float)
            y_unique = pd.to_numeric(
                subset[metric_a], errors='coerce'
            ).to_numpy(dtype=float)
            y_signal = pd.to_numeric(
                subset[metric_b], errors='coerce'
            ).to_numpy(dtype=float)
            slope_unique, slope_signal, inlier_mask = _fit_cpg_slope_pair(
                x,
                y_unique,
                y_signal,
                min_ratio_change=min_slope_ratio_change,
            )
            if not (np.isfinite(slope_unique) and np.isfinite(slope_signal)):
                continue

            unique_rows.append((clone, slope_unique))
            signal_rows.append((clone, slope_signal))

            valid = (
                np.isfinite(x)
                & np.isfinite(y_unique)
                & np.isfinite(y_signal)
                & (x > 0)
                & (y_unique > 0)
                & (y_signal > 0)
            )
            excluded = valid & ~inlier_mask
            if excluded.any():
                excluded_rows = subset.iloc[np.flatnonzero(excluded)]
                if 'Sample_ID' in excluded_rows.columns:
                    labels = excluded_rows['Sample_ID'].astype(str).tolist()
                else:
                    labels = [f'row {index}' for index in excluded_rows.index]
                print(
                    'Excluded CpG slope-ratio influential outlier(s) '
                    f'[{x_col}, CloneID={clone}, '
                    f'min_change={min_slope_ratio_change:.1%}]: '
                    f'{", ".join(labels)}'
                )

        return pd.Series(dict(unique_rows)), pd.Series(dict(signal_rows))

    def _format_slope_value(v: float) -> str:
        if not np.isfinite(v):
            return ''
        absolute = abs(v)
        if absolute != 0 and (absolute >= 10000.0 or absolute < 0.01):
            return f'{v:.2e}'
        return f'{v:.3f}'.rstrip('0').rstrip('.')

    def _annotate_bars_h(ax, bars, values: np.ndarray) -> None:
        if len(values) == 0:
            return
        xmax = float(np.nanmax(values))
        if np.isfinite(xmax) and xmax > 0:
            ax.set_xlim(0, xmax * 1.24)
            offset = xmax * 0.02
        else:
            offset = 0.0
        for bar, value in zip(bars, values, strict=True):
            label = _format_slope_value(float(value))
            if not label:
                continue
            x_position = bar.get_width() + offset
            y_position = bar.get_y() + bar.get_height() / 2.0
            ax.text(
                x_position,
                y_position,
                label,
                ha='left',
                va='center',
                fontsize=9,
                clip_on=False,
            )

    def _plot_slope_panel(
        ax,
        slope: pd.Series,
        title: str,
        xlabel: str,
        letter: str,
    ) -> None:
        values = slope.reindex(clones).dropna()
        if values.empty:
            ax.text(0.5, 0.5, 'Data Missing', transform=ax.transAxes, ha='center', va='center')
            ax.axis('off')
            return
        y = np.arange(len(values))
        bars = ax.barh(
            y,
            values.values,
            color=[clone_color_map.get(clone, UNKNOWN_COLOR) for clone in values.index],
            edgecolor='white',
            linewidth=0.6,
            alpha=0.9,
        )
        ax.set_yticks(y)
        ax.set_yticklabels(
            [clone_display_map.get(clone, clone) for clone in values.index],
            fontsize=8,
        )
        format_ax(ax, title=title, xlabel=xlabel, ylabel='CloneID')
        add_panel_letter(ax, letter)
        _annotate_bars_h(ax, bars, values.values)

    def _plot_ratio_panel(
        ax,
        slope_a: pd.Series,
        slope_b: pd.Series,
        title: str,
        letter: str,
    ) -> None:
        comp_df = pd.DataFrame(
            {
                'Unique': slope_a.reindex(clones),
                'Theoretical': slope_b.reindex(clones),
            }
        ).dropna(subset=['Unique', 'Theoretical'], how='any')
        if not comp_df.empty:
            comp_df = comp_df[
                np.isfinite(comp_df['Unique'])
                & np.isfinite(comp_df['Theoretical'])
                & (comp_df['Theoretical'] != 0)
            ]
        if comp_df.empty:
            ax.text(0.5, 0.5, 'Data Missing', transform=ax.transAxes, ha='center', va='center')
            ax.axis('off')
            return
        ratio = (comp_df['Unique'] / comp_df['Theoretical']).astype(float)
        clone_order = ratio.index.tolist()
        y = np.arange(len(ratio))
        bars = ax.barh(
            y,
            ratio.to_numpy(dtype=float),
            color=[clone_color_map.get(clone, UNKNOWN_COLOR) for clone in clone_order],
            edgecolor='white',
            linewidth=0.6,
            alpha=0.9,
        )
        ax.set_yticks(y)
        ax.set_yticklabels(
            [clone_display_map.get(clone, clone) for clone in clone_order],
            fontsize=8,
        )
        ax.axvline(
            1.0,
            color=Colors.REFLINE,
            linestyle='--',
            linewidth=1.2,
            alpha=0.9,
            zorder=0,
        )
        format_ax(ax, title=title, xlabel='Ratio', ylabel='CloneID')
        add_panel_letter(ax, letter)
        _annotate_bars_h(ax, bars, ratio.to_numpy(dtype=float))

    fig_h = max(28.0, 0.56 * len(clones) + 12.0)
    fig = plt.figure(figsize=(22, fig_h), dpi=PlotParams.DPI)
    gs = fig.add_gridspec(4, 2, height_ratios=[1.0, 1.2, 1.0, 1.2])
    panel_letters = iter('ABCDEF')
    for group_idx, (
        x_col,
        title_denominator,
        xlabel_denominator,
        ratio_prefix,
    ) in enumerate(denominator_groups):
        slope_a, slope_b = calc_slopes_by_clone(x_col)
        row = group_idx * 2
        ax_a = fig.add_subplot(gs[row, 0])
        ax_b = fig.add_subplot(gs[row, 1])
        ax_ratio = fig.add_subplot(gs[row + 1, :])
        _plot_slope_panel(
            ax_a,
            slope_a,
            f'Unique CpG Sites per {title_denominator}',
            f'Unique CpG sites / {xlabel_denominator}',
            next(panel_letters),
        )
        _plot_slope_panel(
            ax_b,
            slope_b,
            f'Host CpG Signals per {title_denominator}',
            f'Host CpG coverage signals / {xlabel_denominator}',
            next(panel_letters),
        )
        _plot_ratio_panel(
            ax_ratio,
            slope_a,
            slope_b,
            f'{ratio_prefix} Normalization: Unique / Coverage-Signal Slope Ratio',
            next(panel_letters),
        )

    legend_levels = [clone_display_map.get(clone, clone) for clone in clones]
    legend_color_map = {
        clone_display_map.get(clone, clone): clone_color_map.get(clone, UNKNOWN_COLOR)
        for clone in clones
    }
    add_right_legend(fig, legend_levels, legend_color_map, title='CloneID')
    add_cloneid_omission_note(fig, omitted, y=0.015)
    plt.tight_layout(rect=[0, 0, 0.885, 0.985], w_pad=2.6, h_pad=2.8)
    if out_path:
        plt.savefig(out_path, bbox_inches='tight')
    plt.show()
    plt.close(fig)


def plot_pca_batch_effect(
    df: pd.DataFrame,
    clone_color_map: dict = None,
    out_path: Path = None,
):
    """将 QC 指标投影到 PCA 空间以检查 CloneID 批次效应。"""
    if df.empty:
        return
    features = [
        'Native_Mapping%',
        'Sample_Unique_CpG_Sites',
        'Sample_Meth_CpG_Rate%',
        'Final_Pair_Yield%',
    ]
    valid_f = [f for f in features if f in df.columns]
    if len(valid_f) >= 3 and len(df) >= 2:
        _load_plot_dependencies()
        plot_df = df.copy()
        plot_df['CloneID'] = plot_df['CloneID'].map(normalize_clone_id_value)
        values = StandardScaler().fit_transform(plot_df[valid_f].fillna(0))
        if np.all(np.var(values, axis=0) <= np.finfo(float).eps):
            print('No variation: PCA skipped.')
            return
        pca = PCA(n_components=2)
        coords = pca.fit_transform(values)
        if not np.all(np.isfinite(pca.explained_variance_ratio_)):
            print('No variation: PCA skipped.')
            return
        group_levels = sorted(
            plot_df['CloneID'].dropna().astype(str).unique().tolist()
        )
        label_map, omitted = make_cloneid_diff_labels(group_levels)
        plot_df['CloneID_Display'] = plot_df['CloneID'].map(label_map)
        if clone_color_map is None:
            clone_color_map = generate_clone_color_map(group_levels)
        display_levels = [label_map.get(group, group) for group in group_levels]
        display_color_map = {
            label_map.get(group, group): clone_color_map.get(group, UNKNOWN_COLOR)
            for group in group_levels
        }
        fig, ax = plt.subplots(figsize=(8, 6), dpi=PlotParams.DPI)
        sns.scatterplot(
            x=coords[:, 0],
            y=coords[:, 1],
            hue=plot_df.get('CloneID_Display'),
            palette=display_color_map,
            ax=ax,
            s=PlotParams.SCATTER_SIZE,
            alpha=PlotParams.SCATTER_ALPHA,
            edgecolor='white',
            linewidth=0.55,
            legend=False,
        )
        format_ax(
            ax,
            title='PCA: Batch Effect Check',
            xlabel=f'PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)',
            ylabel=f'PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)',
        )
        add_right_legend(fig, display_levels, display_color_map, title='CloneID')
        add_cloneid_omission_note(fig, omitted, y=0.015)
        plt.tight_layout(rect=[0, 0.03, 0.885, 0.958], w_pad=2.2, h_pad=2.4)
        if out_path:
            plt.savefig(out_path, bbox_inches='tight')
        plt.show()
        plt.close(fig)
    else:
        print('⚠️ 核心特征不足，无法进行 PCA。')


def plot_final_qc_dashboard_batch(
    df: pd.DataFrame,
    metrics: Any,
    group_col: str = 'CloneID',
    clone_id: Any = None,
    clone_color_map: dict = None,
    out_dir: Path = None,
    filename_prefix: str = 'Final_QC_Dashboard',
    mode: str = 'per_clone_samples',
    figsize: tuple = (16.5, 12),
    simplify_sample_ids: bool = True,
):
    """为指标列表生成一个或多个最终 QC 面板。"""
    metric_list = normalize_metric_list(metrics)
    if not metric_list:
        print("Error: no metric specified.")
        return

    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

    clone_suffix = make_clone_selection_suffix(clone_id)

    for metric_name in metric_list:
        safe_metric = (
            str(metric_name)
            .replace('%', 'pct')
            .replace('/', '_')
            .replace('\\', '_')
            .replace(' ', '_')
            .replace(':', '_')
        )
        out_path = None
        if out_dir is not None:
            out_path = out_dir / f'{filename_prefix}_{mode}_{safe_metric}_{clone_suffix}.pdf'

        plot_final_qc_dashboard(
            df=df,
            group_col=group_col,
            clone_id=clone_id,
            mode=mode,
            metric=metric_name,
            clone_color_map=clone_color_map,
            out_path=out_path,
            figsize=figsize,
            simplify_sample_ids=simplify_sample_ids,
        )


@dataclass(frozen=True)
class NotebookPaths:
    """Notebook 引导所需的项目、Pipeline 与 06_downstream 输出路径。"""

    project_dir: Path
    pipeline_root: Path
    delivery_id: str
    qc_results_dir: Path
    plots_dir: Path
    qc_info: Path
    tss_info: Path
    composition_cache: Path


def _notebook_root_option(value: str | Path | None, environment_name: str) -> str:
    value = value or os.environ.get(environment_name)
    return str(value).strip() if value is not None else ""


def _validate_notebook_pipeline_root(candidate: Path) -> Path:
    candidate = candidate.expanduser().resolve()
    processor = candidate / "downstream" / "downstream_qc.py"
    launcher = candidate / "core" / "run_pipeline.sh"
    if not processor.is_file() or not launcher.is_file():
        raise FileNotFoundError(f"Not a usable Alopex: {candidate}")
    return candidate


def _notebook_project_root(dna_project_root: str | Path | None) -> Path:
    option = _notebook_root_option(dna_project_root, "DNA_PROJECT_ROOT")
    if option:
        candidate = Path(option).expanduser().resolve()
        if not (candidate / "03_results" / "run_manifest.json").is_file():
            raise FileNotFoundError(
                f"DNA_PROJECT_ROOT 没有 sealed delivery: {candidate}"
            )
        return candidate
    start = Path.cwd().resolve()
    for base in (start, *start.parents):
        if (base / "03_results" / "run_manifest.json").is_file():
            return base.resolve()
    raise FileNotFoundError(
        "Cannot find a completed Alopex delivery "
        "(03_results/run_manifest.json); set DNA_PROJECT_ROOT."
    )


def _notebook_pipeline_root(
    dna_pipeline_root: str | Path | None, manifest: Mapping[str, object]
) -> Path:
    if manifest.get("status") != "complete":
        raise FileNotFoundError("run_manifest delivery is not complete")
    pipeline = manifest.get("pipeline")
    if not isinstance(pipeline, dict) or pipeline.get("name") != "Alopex":
        raise FileNotFoundError("run_manifest.pipeline is not Alopex metadata")
    option = _notebook_root_option(dna_pipeline_root, "DNA_PIPELINE_ROOT")
    if option:
        return _validate_notebook_pipeline_root(Path(option))
    source_root = pipeline.get("source_root")
    if not isinstance(source_root, str) or not source_root.strip():
        raise FileNotFoundError("run_manifest.pipeline.source_root is missing")
    return _validate_notebook_pipeline_root(Path(source_root))


def resolve_notebook_paths(
    dna_project_root: str | Path | None = None,
    dna_pipeline_root: str | Path | None = None,
) -> NotebookPaths:
    """按 Notebook 契约解析项目与 Pipeline 根并导出 06_downstream 路径。

    显式填写（参数或同名环境变量）不回退；留空时项目按 Jupyter 当前目录
    向上查找 completed delivery，Pipeline 取 sealed run manifest 记录的
    source_root。downstream_qc.py 自身位于 Pipeline 树内，Notebook 第一格
    的定位 bootstrap 之外的完整校验以本函数为唯一实现。
    """
    project_dir = _notebook_project_root(dna_project_root)
    manifest = load_run_manifest(project_dir, validate_outputs=False)
    pipeline_root = _notebook_pipeline_root(dna_pipeline_root, manifest)
    delivery_id = _manifest_delivery_id(manifest)
    qc_results_dir = project_dir / "06_downstream" / delivery_id / "QC_Results"
    return NotebookPaths(
        project_dir=project_dir,
        pipeline_root=pipeline_root,
        delivery_id=delivery_id,
        qc_results_dir=qc_results_dir,
        plots_dir=qc_results_dir / "plots",
        qc_info=qc_results_dir / "DNAme_QC_Information.csv",
        tss_info=qc_results_dir / "TSS_Profile_Information.csv",
        composition_cache=qc_results_dir.parent / ".cache" / "CpG_Signal_Composition.csv",
    )


def run_qc_processor(
    paths: NotebookPaths,
    *,
    tss_bed: str | Path | None = None,
    skip_tss: bool = False,
    recompute_raw_adata: bool = False,
    recompute_methyl_stats: bool = False,
    recompute_composition: bool = False,
    recompute_gini: bool = False,
    recompute_tss: bool = False,
) -> None:
    """以独立子进程运行 downstream_qc 计算 main 并在失败时终止 Notebook。

    线程环境变量对齐当前 Notebook allocation（SLURM_CPUS_ON_NODE 或全部
    CPU），使 Processor 留在 allocation 内的同时其 native 库可用全部核；
    recompute 开关透传 CLI 的 --recompute-* 重算参数，其中 raw_adata 强制
    从 sealed CpG 重建 RawAdata.h5ad，其余默认按缓存身份自动判定复用。
    """
    if skip_tss and tss_bed is not None:
        raise ValueError("skip_tss cannot be combined with an explicit TSS bed")
    processor = paths.pipeline_root / "downstream" / "downstream_qc.py"
    if not processor.is_file():
        raise FileNotFoundError(f"Downstream QC processor is missing: {processor}")
    worker_count = os.environ.get("SLURM_CPUS_ON_NODE", str(os.cpu_count() or 1))
    command = [
        sys.executable,
        str(processor),
        "--project-dir",
        str(paths.project_dir),
        "--workers",
        worker_count,
    ]
    if tss_bed is not None:
        command.extend(["--tss-bed", str(tss_bed)])
    elif skip_tss:
        command.append("--skip-tss")
    for flag, enabled in (
        ("--recompute-raw-adata", recompute_raw_adata),
        ("--recompute-methyl-stats", recompute_methyl_stats),
        ("--recompute-composition", recompute_composition),
        ("--recompute-gini", recompute_gini),
        ("--recompute-tss", recompute_tss),
    ):
        if enabled:
            command.append(flag)
    environment = os.environ.copy()
    for variable in _THREAD_LIMIT_ENV_VARS:
        environment[variable] = worker_count
    _log(f"启动 QC Processor 子进程: {processor}")
    completed = subprocess.run(command, check=False, env=environment)
    if completed.returncode:
        raise RuntimeError(f"QC Processor failed with exit code {completed.returncode}")


def prepare_single_cpg_adata(
    paths: NotebookPaths,
    df_qc: pd.DataFrame,
    *,
    recompute: bool = False,
) -> Path:
    """从当前 RawAdata 生成 single-CpG 矩阵，按细胞身份整合 QC 并关闭全部句柄。

    使用当前物种的 BED3、fraction/mean 和 chunk_size=500；保留全部细胞。
    已验证的矩阵可复用，仅 QC/HQ 改变时只更新 obs。未知或过期文件须显式
    recompute；新矩阵完成并关闭后才替换正式文件，不关闭其他 Notebook 的句柄。
    """
    if snap is None:
        raise ImportError("Single-CpG generation requires the Alopex notebook environment")
    sealed = load_sealed_qc_context(
        paths.project_dir, expected_pipeline_root=paths.pipeline_root
    )
    if sealed.delivery_id != paths.delivery_id:
        raise ValueError("Delivery changed; rerun the Notebook from its first cell")
    run_paths = build_paths(paths.project_dir, sealed.delivery_id)
    if run_paths.qc_results_dir != paths.qc_results_dir:
        raise ValueError("Notebook QC directory differs from the current delivery")
    sample_ids = [cell["sample_id"] for cell in read_active_cells(
        paths.project_dir, sealed_context=sealed
    )]
    if not sample_ids:
        raise ValueError("Single-CpG generation requires at least one cell")
    if not df_qc.columns.is_unique or not set(QC_INFO_COLUMNS).issubset(df_qc.columns):
        raise ValueError("df_qc must contain the complete QC table with unique column names")
    qc = _require_exact_sample_ids(df_qc, sample_ids, "Notebook QC")
    saved_qc = pd.read_csv(
        paths.qc_info, converters=dict.fromkeys(QC_INFO_COLUMNS[:3], str),
        float_precision="round_trip",
    )
    saved_qc = _require_exact_sample_ids(saved_qc, sample_ids, "Saved QC")
    pd.testing.assert_frame_equal(
        qc[list(QC_INFO_COLUMNS)], saved_qc[list(QC_INFO_COLUMNS)],
        check_dtype=False, check_exact=True,
    )
    raw_path = run_paths.raw_adata
    raw_schema = json.loads(raw_path.with_name("RawAdata.input_schema.json").read_text())
    if (not _raw_adata_schema_matches(raw_path, raw_schema)
            or raw_schema["delivery_id"] != sealed.delivery_id
            or set(raw_schema["sample_ids"]) != set(sample_ids)):
        raise ValueError("RawAdata is stale or unverified; rerun the QC Processor")
    genome = str(sealed.config["species"])
    bed = paths.pipeline_root / "resources" / f"{genome}_reference" / "cpg" / f"{genome}.single_cpg.bed.gz"
    signature = {
        "algorithm": "single_cpg_fraction_mean_v1",
        "snapatac2": importlib.metadata.version("snapatac2"),
        "delivery_id": sealed.delivery_id,
        "raw_generation": raw_schema["generation_id"],
        "raw_h5ad": raw_schema["h5ad"],
        "bed": _content_binding(bed),
        "value_type": "fraction", "summary_type": "mean", "chunk_size": 500,
    }
    output = paths.qc_results_dir / "SingleCpG_Adata.h5ad"
    schema_path = output.with_name("SingleCpG_Adata.input_schema.json")
    cached = None
    if schema_path.is_file():
        try:
            cached = json.loads(schema_path.read_text())
        except (OSError, ValueError):
            pass
    reuse = (output.is_file() and isinstance(cached, dict)
             and cached.get("input") == signature
             and cached.get("h5ad") == _raw_adata_stat_binding(output))
    if output.exists() and not reuse and not recompute:
        raise ValueError(
            f"Existing SingleCpG file is unverified or stale: {output}. "
            "Set RECOMPUTE_SINGLECPG=True to rebuild it explicitly."
        )
    temporary = output.with_name(f".{output.stem}.{uuid.uuid4().hex}.tmp.h5ad")
    raw = matrix = None
    try:
        if reuse and not recompute:
            matrix = snap.read(str(output), backed="r")
        else:
            _log("生成 single-CpG 矩阵：fraction/mean，chunk_size=500")
            raw = snap.read(str(raw_path), backed="r")
            _require_exact_sample_ids(qc, list(raw.obs_names), "RawAdata QC")
            matrix = snap.pp.make_peak_matrix(
                raw, value_type="fraction", summary_type="mean",
                peak_file=str(bed), inplace=False, file=str(temporary), chunk_size=500,
            )
        names = list(matrix.obs_names)
        aligned = _require_exact_sample_ids(qc, names, "SingleCpG QC")
        current = matrix.obs[:]
        merged = current.to_pandas()
        for column in aligned.columns:
            merged[column] = aligned[column].array
        updated = pl.from_pandas(merged)
        if not updated.equals(current):
            if reuse and not recompute:
                matrix.close()
                matrix = None
                matrix = snap.read(str(output), backed="r+")
                schema_path.unlink(missing_ok=True)
            matrix.obs = updated
            matrix.obs_names = names
        matrix.close()
        matrix = None
        if not reuse or recompute:
            _fsync_file(temporary)
            schema_path.unlink(missing_ok=True)
            os.replace(temporary, output)
        _atomic_write_json(
            {"input": signature, "h5ad": _raw_adata_stat_binding(output)}, schema_path
        )
        _log(f"SingleCpG 已就绪，{len(names)} 个细胞，文件已关闭: {output}")
        return output
    finally:
        try:
            if matrix is not None:
                matrix.close()
        finally:
            try:
                if raw is not None:
                    raw.close()
            finally:
                temporary.unlink(missing_ok=True)


def load_notebook_qc_frame(
    paths: NotebookPaths,
) -> tuple[pd.DataFrame, dict[str, str] | None]:
    """加载 QC Processor 产出的 DNAme QC 表并初始化绘图主题与 CloneID 配色。"""
    paths.qc_results_dir.mkdir(parents=True, exist_ok=True)
    paths.plots_dir.mkdir(parents=True, exist_ok=True)
    setup_nature_style()
    if paths.qc_info.exists():
        df_qc = pd.read_csv(
            paths.qc_info, converters=dict.fromkeys(QC_INFO_COLUMNS[:3], str),
            float_precision="round_trip",
        )
        print(f"✅ 已加载 QC 数据: {len(df_qc)} 个样本")
    else:
        df_qc = pd.DataFrame()
        print(f"❌ 找不到 QC 数据文件: {paths.qc_info}")
        print("请先成功运行上一格 QC Processor；该步骤需要 snapatac2、anndata、polars、pandas、numpy。")
    clone_color_map: dict[str, str] | None = None
    if not df_qc.empty and "CloneID" in df_qc.columns:
        df_qc["CloneID"] = df_qc["CloneID"].map(normalize_clone_id_value)
        unique_clones = sorted(df_qc["CloneID"].dropna().unique().tolist())
        clone_color_map = generate_clone_color_map(unique_clones)
    return df_qc, clone_color_map
