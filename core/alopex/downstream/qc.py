"""从sealed交付生成甲基化QC、SingleCpG 矩阵和 1 Mb 整数计数与Notebook接入。"""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import importlib.metadata
import json
import os
import re
import stat as _qc_lock_stat
import sys

_LOADED_MODULE_CODE = sys._getframe().f_code
import tempfile
import threading
import time
import uuid
from collections import deque
from collections.abc import Callable, Iterable, Mapping
from concurrent.futures import ThreadPoolExecutor
import contextlib
from contextlib import contextmanager
from dataclasses import dataclass, field as dataclass_field
from itertools import islice
from pathlib import Path
from typing import Any

from ..config import (
    DD_CONSENSUS,
    DD_PROTOCOLS,
    DROPLET_DESIGN_SHA256,
    CPG_REPRESENTATIONS,
    get_backend,
    load_project_config,
    primary_analysis_contigs,
    read_fasta_index,
    reference_sequence_identity,
    validate_downstream_reference,
    validate_tss_reference,
)
from ..delivery import (load_run_manifest, load_run_snapshot, output_file_state,
                        read_sealed_final_manifest, run_manifest_state, verified_published_output_states)
from ..identity import PythonSourceGraph, capture_loaded_source_codes, assert_loaded_sources_current
from ..io import sha256_file, project_temporary_parent, stat_if_present
from ..paths import CORE_ROOT as _CORE_ROOT, PIPELINE_ROOT as _PIPELINE_ROOT
from ..runtime import effective_cpu_count


FINAL_MANIFEST = Path("03_results/QC_Results/sample_manifest.tsv")
MAX_CPG_WORKERS = 8

QC_FILTER_DEFAULTS = {
    "mapping_threshold": None,
    "lambda_max": 10.0,
    "puc19_min": 90.0,
    "cpg_min": 500000.0,
    "gini_max": 0.5,
}

REQUIRED_MANIFEST_FIELDS = {
    "protocol",
    "sample_id",
    "project_sample_id",
    "plate_id",
    "dna_barcode",
    "cell_order",
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
    marker_state: tuple[int, ...]
    output_states: Mapping[str, tuple[int, ...]]

    def assert_current(self) -> None:
        """确认公开完成凭证与最终样本清单仍属加载时的同一次交付。"""
        _check_delivery_file(self.run_manifest_path, self.marker_state, manifest=True)
        self.assert_output_current("QC_Results/sample_manifest.tsv")

    def assert_output_current(self, key: str) -> None:
        _check_delivery_file(self.project_dir / "03_results" / key, self.output_states[key])


def _check_delivery_file(path: Path, expected: tuple[int, ...], *, manifest: bool = False) -> None:
    try:
        current = run_manifest_state(path) if manifest else output_file_state(path.lstat())
    except OSError as exc:
        raise ValueError(f"Sealed delivery changed: {path}; rerun downstream from the current delivery") from exc
    if current != expected:
        raise ValueError(f"Sealed delivery changed: {path}; rerun downstream from the current delivery")


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"Sealed delivery {label} must be a mapping")
    return value


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
    if reference.get("content") != reference_sequence_identity(reference_path):
        raise ValueError("Sealed reference FASTA/FAI content changed after the Pipeline run")
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
    marker_state = run_manifest_state(manifest_path)
    manifest = load_run_manifest(manifest_path, validate_outputs=False)
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
    snapshot_config_path = snapshot_path.parent / "inputs" / str(config_copy["snapshot_path"])
    config = load_project_config(snapshot_config_path)

    analysis = _mapping(config.get("analysis"), "snapshot config analysis")
    run_species = str(run.get("species", "")).strip()
    run_protocol = str(run.get("protocol", "")).strip()
    config_species = str(config.get("species", "")).strip()
    config_protocol = str(analysis.get("protocol", "")).strip()
    config_backend = get_backend(config)
    if not run_species or config_species != run_species:
        raise ValueError("Sealed config species disagrees with the run manifest")
    if run_protocol not in {"cabernet", "taps", *DD_PROTOCOLS} or config_protocol != run_protocol:
        raise ValueError("Sealed config protocol disagrees with the run manifest")

    if run_protocol in DD_PROTOCOLS:
        if run.get("tag_source") != DD_CONSENSUS[run_protocol].source or run.get("barcode_design_sha256") != DROPLET_DESIGN_SHA256:
            raise ValueError("Sealed DNA tag source or barcode design disagrees with its protocol")

    reference_path = Path(str(run.get("reference_path", ""))).expanduser().resolve()
    if not reference_path.is_file():
        raise FileNotFoundError(f"Sealed Pipeline reference is missing: {reference_path}")
    reference_fai = _snapshot_reference_fai(snapshot, reference_path)

    run_config = _mapping(run, "run metadata")
    cpg_representation = str(run_config.get("cpg_representation", ""))
    coordinate_base = analysis["cpg_coordinate_base"]
    expected_representation = CPG_REPRESENTATIONS[coordinate_base]
    if cpg_representation != expected_representation:
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

    output_states = verified_published_output_states(project, manifest)
    output_inventory = {Path(item["path"]).as_posix(): item for item in manifest["outputs"]}
    if "QC_Results/sample_manifest.tsv" not in output_inventory:
        raise ValueError("Sealed delivery inventory omits QC_Results/sample_manifest.tsv")
    context = SealedQCContext(
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
        marker_state=marker_state,
        output_states=output_states,
    )
    context.assert_current()
    return context


def _inventory_record_for_project_file(
    context: SealedQCContext, value: str | Path, *, label: str
) -> tuple[Path, Mapping[str, object]]:
    relative = Path(str(value))
    if relative.is_absolute() or ".." in relative.parts or relative == Path("."):
        raise ValueError(f"Current Pipeline {label} has an unsafe path: {value!r}")
    path = context.project_dir / relative
    public_root = context.project_dir / "03_results"
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
    context.assert_output_current(inventory_key)
    if record["size_bytes"] == 0:
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
    context.assert_current()
    manifest, record = _inventory_record_for_project_file(
        context, FINAL_MANIFEST, label="final sample manifest"
    )
    reader = read_sealed_final_manifest(manifest, record)
    missing = sorted(REQUIRED_MANIFEST_FIELDS - set(reader.fieldnames or ()))
    if missing:
        raise ValueError(
            "Current sample_manifest.tsv is missing required field(s): "
            + ", ".join(missing)
        )

    active: list[dict[str, Any]] = []
    seen: set[str] = set()
    for raw_row in reader:
        row = {str(key): str(value or "").strip() for key, value in raw_row.items()}
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
            "biscuit": "biscuit_nonsecondary_mapq_ge_40",
            "bismark": "bismark_unique_concordant_pairs",
            "rastair": "bwa_proper_primary_any_mate_mapq20",
        }[backend]
        expected_native_unit = {
            "biscuit": "nonsecondary_alignment_records_postdedup",
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
            "bismark": DD_CONSENSUS[row["protocol"]].dedup_policy if row["protocol"] in DD_PROTOCOLS else "bismark_paired_endpoint_orientation",
            "rastair": "samtools_markdup_flag_exclude_in_rastair",
        }[backend]
        if cell["dedup_policy"] != expected_dedup_policy:
            raise ValueError(f"Cell {sample_id!r} has an invalid dedup_policy")
        expected_non_cpg_source = {
            "biscuit": "biscuit_cph_retention_by_read_position",
            "bismark": DD_CONSENSUS[row["protocol"]].cph_source if row["protocol"] in DD_PROTOCOLS else "bismark_extraction_chg_chh",
            "rastair": "taps_oriented_host_chg_chh_conversion",
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
            if row.get("high_cph_role") != "not_applicable":
                raise ValueError(f"Cell {sample_id!r} has an invalid TAPS high-CpH role")
            modified = sum(_required_nonnegative_int(row, f"non_cpg_{ctx}_mod", sample_id) for ctx in ("chg", "chh"))
            total = modified + sum(_required_nonnegative_int(row, f"non_cpg_{ctx}_unmod", sample_id) for ctx in ("chg", "chh"))
            if total != _required_nonnegative_int(row, "non_cpg_observations", sample_id):
                raise ValueError(f"Cell {sample_id!r} has inconsistent TAPS non-CpG counts")
            expected_pct = 100.0 * modified / total if total else float("nan")
            if not (np.isnan(expected_pct) and np.isnan(cell["non_cpg_methylation_pct"])) and not np.isclose(expected_pct, cell["non_cpg_methylation_pct"], rtol=1e-12, atol=1e-12):
                raise ValueError(f"Cell {sample_id!r} has inconsistent TAPS non-CpG percentage")
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
        if row["protocol"] in DD_PROTOCOLS:
            tag = DD_CONSENSUS[row["protocol"]]
            if (cell["duplicate_pairs"] != 0 or row.get("molecule_qc_metric_unit") != tag.unit
                    or row.get("tag_source") != tag.source):
                raise ValueError(f"Cell {sample_id!r} mixes Droplet pair and molecule counts")
            read_count, umi_count, ambiguous_count = (
                _required_nonnegative_int(row, key, sample_id) for key in
                ("cpg_read_observations", tag.count_field, tag.ambiguous_field)
            )
            if umi_count + ambiguous_count > read_count:
                raise ValueError(f"Cell {sample_id!r} has inconsistent Droplet molecule counts")
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
    context.assert_current()
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
            "cell_order": cell["cell_order"],
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
    try:
        key = path.relative_to(context.project_dir / "03_results").as_posix()
    except ValueError as exc:
        raise ValueError(f"Downstream input escapes sealed outputs: {path}") from exc
    record = context.output_inventory.get(key)
    if record is None:
        raise ValueError(f"Downstream input is absent from delivery inventory: {key}")
    context.assert_output_current(key)
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
    sealed.assert_current()

    output_records: dict[str, Mapping[str, object]] = {
        "QC_Results/sample_manifest.tsv": sealed.output_inventory["QC_Results/sample_manifest.tsv"]
    }
    for record in records:
        path = Path(record["cpg_file"])
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
    sealed.assert_current()
    return payload


def _add_cli_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument(
        "--project-dir",
        type=str,
        required=True,
        help="Alopex project root containing 00_config and 03_results.",
    )
    parser.add_argument("--expected-delivery-id", help=argparse.SUPPRESS)
    parser.add_argument(
        "--workers",
        type=int,
        default=0,
        help=f"CPU budget for CpG processing. 0 uses the available CPUs, with at most {MAX_CPG_WORKERS} cell workers.",
    )
    parser.add_argument(
        "--recompute-cpg-stats",
        action="store_true",
        help="Rebuild the single-pass downstream products from sealed CpG inputs.",
    )
    parser.add_argument(
        "--recompute-gini",
        action="store_true",
        help="Recompute the downstream products, including Gini.",
    )
    parser.add_argument(
        "--recompute-tss",
        action="store_true",
        help="Recompute the downstream products, including TSS.",
    )
    parser.add_argument(
        "--tss-bed",
        type=str,
        default=None,
        help="TSS window BED; defaults to Pipeline resources/<species>_reference/tss.",
    )
    parser.add_argument(
        "--tss-sample-fraction",
        type=float,
        default=0.1,
        help="TSS sampling fraction per CloneID; <10 cells use all, otherwise cap at 100. 1 uses all.",
    )
    parser.add_argument(
        "--tss-sample-seed",
        type=int,
        default=0,
        help="Nonnegative seed for reproducible TSS cell sampling; default: 0.",
    )
    parser.add_argument(
        "--skip-tss",
        action="store_true",
        help="Generate all non-TSS outputs when TSS resources are intentionally unavailable.",
    )
    for name, default in QC_FILTER_DEFAULTS.items():
        parser.add_argument("--" + name.replace("_", "-"), type=float, default=default,
                            help=f"SingleCpG HQ filter threshold (default: {default}).")
    return parser


def _early_print_help_if_requested() -> None:
    if not any(arg in {"-h", "--help"} for arg in sys.argv[1:]):
        return
    _add_cli_args(argparse.ArgumentParser(prog="python -m alopex.downstream.qc", description=__doc__)).print_help()
    raise SystemExit(0)


if __name__ == "__main__":
    _early_print_help_if_requested()

# CLI帮助先于科学依赖加载，计算worker不导入绘图库。
import numpy as np
import pandas as pd


try:
    import anndata as ad
    from backports import zstd
    from . import matrix as matrices, diagnostics
    from numba import njit
except ImportError:
    ad = matrices = diagnostics = None

    def njit(*args, **kwargs):
        def _identity(fn):
            return fn
        if len(args) == 1 and callable(args[0]) and not kwargs:
            return args[0]
        return _identity


_LOADED_SOURCE_CODES = capture_loaded_source_codes()


@dataclass(frozen=True)
class RunPaths:
    project_dir: Path
    cpg_dir: Path
    qc_results_dir: Path
    cache_dir: Path
    single_cpg: Path
    single_schema: Path
    doublet_counts: Path
    qc_info: Path
    tss_info: Path
    composition_cache: Path
    input_identity: Path


@dataclass(frozen=True)
class RunConfig:
    genome: str
    ref_fai: Path
    tss_bed: Path | None
    gini_bin_size: int = 500000
    tss_sample_fraction: float = 0.1
    tss_sample_seed: int = 0
    qc_filter_parameters: Mapping[str, float | None] = dataclass_field(default_factory=QC_FILTER_DEFAULTS.copy)


def _log(message: str, level: str = "INFO") -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {level:7} {message}", flush=True)


def _pipeline_dir() -> Path:
    return _PIPELINE_ROOT


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
        try:
            mode = current.lstat().st_mode
        except FileNotFoundError:
            if create:
                current.mkdir(exist_ok=True)
                mode = current.lstat().st_mode
            else:
                continue
        if _qc_lock_stat.S_ISLNK(mode):
            raise ValueError(f"Project directory path must not contain symlinks: {current}")
        if not _qc_lock_stat.S_ISDIR(mode):
            raise ValueError(f"Project directory path is not a directory: {current}")
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


_QC_WRITE_GUARD = threading.Lock()


_QC_WRITE_LOCKS: dict[Path, tuple[int, int]] = {}


def _qc_locks_after_fork() -> None:
    # fork 子进程不能延长父事务的锁寿命，也不能对共享 open description 执行 unlock。
    for _, descriptor in _QC_WRITE_LOCKS.values():
        os.close(descriptor)
    _QC_WRITE_LOCKS.clear()
    _QC_WRITE_GUARD.release()


os.register_at_fork(before=_QC_WRITE_GUARD.acquire,
                    after_in_parent=_QC_WRITE_GUARD.release,
                    after_in_child=_qc_locks_after_fork)


@contextmanager
def _qc_write_lock(directory: Path):
    """保护下游目录事务；跨进程/线程冲突立即拒绝，同线程嵌套复用锁。

    锁文件保留原 inode，绝不 unlink；异常或进程退出由关闭描述符释放锁。
    上游仍只使用 Snakemake 锁，sealed marker 的前后复验保持独立。
    """
    directory = Path(directory)
    if directory.is_symlink():
        raise ValueError(f"QC transaction directory must not be a symlink: {directory}")
    directory = directory.resolve()
    owner = threading.get_ident()
    process_id = os.getpid()
    busy = f"Downstream QC is busy: {directory}; retry after the current writer finishes"
    with _QC_WRITE_GUARD:
        held = _QC_WRITE_LOCKS.get(directory)
        nested = held is not None and held[0] == owner
        if held is not None and not nested:
            raise RuntimeError(busy)
        if not nested:
            descriptor = os.open(directory / ".qc-write.lock", os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW, 0o600)
            try:
                if not _qc_lock_stat.S_ISREG(os.fstat(descriptor).st_mode):
                    raise ValueError(f"QC lock must be a regular file: {directory}")
                try:
                    fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError as exc:
                    raise RuntimeError(busy) from exc
                _QC_WRITE_LOCKS[directory] = (owner, descriptor)
            except BaseException:
                os.close(descriptor)
                raise
    try:
        yield
    finally:
        if not nested and os.getpid() == process_id:
            with _QC_WRITE_GUARD:
                del _QC_WRITE_LOCKS[directory]
                # close 足以释放 flock；避免任何清理路径删除另一个进程使用的锁文件。
                os.close(descriptor)


def resolve_qc_results_dir(project_dir: Path) -> Path:
    """固定输出到06_downstream；是否可复用由完成身份决定，不由目录名决定。"""
    project_dir = project_dir.resolve()
    root = _project_directory(project_dir, Path("06_downstream"), create=False)
    root.mkdir(exist_ok=True)
    return root


def build_paths(
    project_dir: Path,
) -> RunPaths:
    project_dir = project_dir.resolve()
    cpg_dir = _project_directory(
        project_dir, Path("03_results/CpG"), create=False
    )
    qc_results_dir = resolve_qc_results_dir(project_dir)

    cache_dir = _project_directory(
        project_dir, qc_results_dir.relative_to(project_dir) / ".cache", create=False
    )
    cache_dir.mkdir(exist_ok=True)

    return RunPaths(
        project_dir=project_dir,
        cpg_dir=cpg_dir,
        qc_results_dir=qc_results_dir,
        cache_dir=cache_dir,
        single_cpg=qc_results_dir / "SingleCpG_Adata.h5ad",
        single_schema=qc_results_dir / "SingleCpG_Adata.input_schema.json",
        doublet_counts=qc_results_dir / "DNA_Doublet_1Mb_Counts.npz",
        qc_info=qc_results_dir / "DNAme_QC_Information.csv",
        tss_info=qc_results_dir / "TSS_Profile_Information.csv",
        composition_cache=cache_dir / "CpG_Signal_Composition.csv",
        input_identity=cache_dir / "DNAme_QC_Input_Identity.json",
    )


def build_config(
    sealed: SealedQCContext,
    *,
    tss_bed_override: Path | None = None,
    skip_tss: bool = False,
    tss_sample_fraction: float = 0.1,
    tss_sample_seed: int = 0,
    qc_filter_parameters: Mapping[str, float | None] | None = None,
) -> RunConfig:
    """解析 reference FAI，并在需要时强制检查共享 TSS 输入。"""
    config = sealed.config
    _validate_tss_sampling(tss_sample_fraction, tss_sample_seed)
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
        tss_sample_fraction=tss_sample_fraction,
        tss_sample_seed=tss_sample_seed,
        qc_filter_parameters={**QC_FILTER_DEFAULTS, **(qc_filter_parameters or {})},
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
    if cfg.tss_bed is not None:
        validate_tss_reference(sealed.reference_path, cfg.tss_bed)
    return cfg


COMPOSITION_COV_COLS = (
    "lambda_cov",
    "puc19_cov",
    "mt_cov",
    "host_cov",
    "total_cov",
    "Sample_CpG_Signals",
)


COMPOSITION_RATES = {
    "host_cov": "Signal_Host_Rate%",
    "lambda_cov": "Signal_Lambda_Rate%",
    "puc19_cov": "Signal_pUC19_Rate%",
    "mt_cov": "Signal_mtDNA_Rate%",
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

_QC_SITE_COUNTS = (
    "Sample_Unique_CpG_Sites",
    "Lambda_Unique_CpG_Sites",
    "pUC19_Unique_CpG_Sites",
)


def _read_qc_table(path: Path) -> pd.DataFrame:
    """按 CSV 原文读取位点计数与身份，比例保留既有的缺失值和精度。"""
    frame = pd.read_csv(
        path, converters=dict.fromkeys(QC_INFO_COLUMNS[:3], str),
        dtype=dict.fromkeys(_QC_SITE_COUNTS, str), float_precision="round_trip",
    )
    for column in _QC_SITE_COUNTS:
        frame[column] = _nonnegative_counts(frame[column])
    return frame


def _nonnegative_counts(values: pd.Series) -> pd.Series:
    counts = matrices.exact_integers(values, np.int64)
    if (counts < 0).any():
        raise ValueError(f"Negative counts in {values.name}")
    return counts


_THREAD_LIMIT_ENV_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)


def _chromosome_class(name: str, host_contigs: set[str]) -> str:
    if name in host_contigs:
        return "host"
    if name in {"lambda", "pUC19"}:
        return name
    if name.casefold() in {"chrm", "mt", "chrmt"}:
        return "mtDNA"
    return "auxiliary"


def _bounded_cpg_workers(requested: int, task_count: int) -> int:

    if task_count <= 0:
        return 0
    return min(effective_cpu_count(requested), MAX_CPG_WORKERS, task_count)


def _atomic_write_csv(df: pd.DataFrame, path: Path, *, before_commit=None) -> None:

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
            dir=project_temporary_parent(path),
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            df.to_csv(handle, index=False)
            handle.flush()
            os.fsync(handle.fileno())
        if before_commit is not None:
            before_commit()
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


def _atomic_write_json(payload: dict, path: Path, *, before_commit=None) -> None:

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
            dir=project_temporary_parent(path),
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        if before_commit is not None:
            before_commit()
        os.replace(temporary, path)
        temporary = None
        try:
            _fsync_directory(path.parent)
        except OSError:
            # 完成凭证的目录提交失败，不能由下一次读取批准本次结果。
            path.unlink(missing_ok=True)
            raise
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _fsync_file(path: Path) -> None:

    file_fd = os.open(str(Path(path)), os.O_RDONLY)
    try:
        os.fsync(file_fd)
    finally:
        os.close(file_fd)


def _published_file_binding(path: Path) -> dict:

    path = Path(path)
    stat = path.stat()
    return {
        "name": path.name,
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "device": stat.st_dev,
        "inode": stat.st_ino,
    }


def _reference_guard(*paths: Path | None, before_commit=None) -> Callable[[], None]:
    """绑定本次读取的参考文件状态；仅瞬时复验，不把mtime/ctime写入内容身份。"""
    def state(path):
        stat = path.stat()
        return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns

    observed = [(Path(path), state(Path(path))) for path in paths if path is not None]

    def assert_current():
        if before_commit is not None:
            before_commit()
        for path, expected in observed:
            try:
                current = state(path)
            except OSError as exc:
                raise ValueError(f"Reference input changed during downstream computation: {path}") from exc
            if current != expected:
                raise ValueError(f"Reference input changed during downstream computation: {path}")
    return assert_current


def _content_binding(path: Path) -> dict:
    path = Path(path).resolve()
    assert_current = _reference_guard(path)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    assert_current()
    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
    }


def _invalidate_completion(path: Path) -> None:
    """替换任一缓存之前撤销旧完成身份；失败由下一次 Processor 恢复。"""
    try:
        path.unlink()
    except FileNotFoundError:
        return
    _fsync_directory(path.parent)


def _remove_tss_plot(paths: RunPaths) -> None:
    plots = _project_directory(paths.project_dir, Path("06_downstream/plots"), create=False)
    (plots / "TSS_Profile.pdf").unlink(missing_ok=True)


def _remove_skipped_tss(paths: RunPaths) -> None:
    _remove_tss_plot(paths)
    paths.tss_info.unlink(missing_ok=True)


@dataclass(frozen=True)
class CompletedQCContext:
    """当前 QC、SingleCpG 和适用的 1 Mb 计数的共同完成身份。"""

    identity: Mapping[str, Any]
    sample_ids: tuple[str, ...]
    assert_current: Callable[[], None]


def _single_bed(sealed: SealedQCContext) -> Path:
    genome = str(sealed.config["species"])
    return sealed.pipeline_root / "resources" / f"{genome}_reference" / "cpg" / f"{genome}.single_cpg.bed.gz"


def _require_completed_qc(
    paths: RunPaths, sealed: SealedQCContext, *, algorithms=None, verify_diagnostics=(),
) -> CompletedQCContext:
    assert_loaded_sources_current(_LOADED_SOURCE_CODES)
    try:
        qc_guard = _reference_guard(paths.input_identity, paths.qc_info)
        identity = json.loads(paths.input_identity.read_text())
        cells = read_active_cells(paths.project_dir, sealed_context=sealed)
        sample_ids = tuple(cell["sample_id"] for cell in cells)
        cpg_paths = tuple(cell["cpg_output_identity"]["path"] for cell in cells)

        def assert_current():
            assert_loaded_sources_current(_LOADED_SOURCE_CODES)
            sealed.assert_current()
            for key in cpg_paths:
                sealed.assert_output_current(key)
            try:
                qc_guard()
            except ValueError as exc:
                raise ValueError("QC completion changed; rerun the QC Processor") from exc

        algorithms = _algorithm_fingerprints() if algorithms is None else algorithms
        generation = identity["cpg_generation"]
        expected_artifacts = {paths.composition_cache.name, diagnostics.CHROMOSOME_INFO, diagnostics.FILTER_INFO}
        if identity["tss_sampling"] is not None:
            expected_artifacts.add(paths.tss_info.name)
        if set(identity["diagnostic_artifacts"]) != expected_artifacts:
            raise ValueError("incomplete QC artifacts")
        if (not isinstance(generation, str) or re.fullmatch(r"[0-9a-f]{32}", generation) is None
                or any(identity.get(key) != getattr(sealed, key) for key in
                       ("delivery_id", "snapshot_id", "run_manifest_sha256", "snapshot_config_sha256"))
                or identity["algorithms"] != algorithms
                or identity["qc_table"] != _content_binding(paths.qc_info)):
            raise ValueError("stale QC completion identity")
        assert_current = _reference_guard(paths.single_schema, paths.single_cpg, before_commit=assert_current)
        _require_current_single(paths, identity)
        artifacts = [paths.composition_cache if name == paths.composition_cache.name
                     else paths.qc_results_dir / name for name in sorted(expected_artifacts)]
        assert_current = _reference_guard(*artifacts, before_commit=assert_current)
        for artifact in artifacts:
            if verify_diagnostics is None or artifact.name in verify_diagnostics:
                if _content_binding(artifact) != identity["diagnostic_artifacts"][artifact.name]:
                    raise ValueError(f"stale QC artifact: {artifact.name}")
        if sealed.config["analysis"]["protocol"] in DD_PROTOCOLS:
            assert_current = _reference_guard(paths.doublet_counts, before_commit=assert_current)
            if (identity["doublet_counts"] != _published_file_binding(paths.doublet_counts)
                    or identity["doublet_counts_ctime_ns"] != paths.doublet_counts.lstat().st_ctime_ns):
                raise ValueError("stale 1 Mb integer counts")
        assert_current()
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
        raise ValueError("QC is incomplete or stale; rerun the QC Processor") from exc
    return CompletedQCContext(identity, sample_ids, assert_current)


def _require_current_single(paths, identity):
    schema = json.loads(paths.single_schema.read_text())
    source = schema["input"]
    if (source["cpg_generation"] != identity["cpg_generation"]
            or source["delivery_id"] != identity["delivery_id"]
            or source["source_cpg_representation"] != identity["cpg_representation"]
            or source["cpg_coordinate_base"] != _source_cpg_coordinate_base(identity["cpg_representation"])
            or source.get("cpg_interval_convention") != matrices.SINGLE_CPG_INTERVAL_CONVENTION
            or source["bed_sha256"] != identity["context"]["bed"]
            or source["qc_filter_parameters"] != identity["qc_filter_parameters"]
            or schema["h5ad"] != _published_file_binding(paths.single_cpg)
            or identity["single_cpg_ctime_ns"] != paths.single_cpg.lstat().st_ctime_ns):
        raise ValueError("SingleCpG is stale; rerun the QC Processor")
    return schema


def _read_chrom_sizes(fai_path: Path) -> dict[str, int]:
    return {name: length for name, length, *_ in read_fasta_index(fai_path)}


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


def _source_cpg_coordinate_base(representation: str) -> int:
    for base, name in CPG_REPRESENTATIONS.items():
        if representation == name:
            return base
    raise ValueError(f"Unsupported sealed CpG representation: {representation!r}")


@contextmanager
def _managed_cpg_text(path: Path):
    """每个细胞独占一个流式解码器，完整读取所有帧及校验和。"""
    if path.suffix != ".zst":
        raise ValueError(f"CpG input must use the current .cpg.tsv.zst contract: {path}")
    try:
        with zstd.open(path, "rt", encoding="utf-8") as handle:
            yield handle
    except (zstd.ZstdError, EOFError, UnicodeDecodeError) as exc:
        raise ValueError(f"Cannot decode CpG input {path}: {exc}") from exc


def _algorithm_fingerprints() -> dict[str, str]:

    sources = PythonSourceGraph((_CORE_ROOT,), loaded_sources=_LOADED_SOURCE_CODES)
    versions = {name: importlib.metadata.version(name)
                for name in ("numpy", "pandas", "scipy", "numba", "anndata", "h5py")}
    roots = {
        "input": ("build_qc_records", "source_identity"),
        "products": ("_calculate_cpg_products", "_SinglePassCells", "_merge_pipeline_qc", "_single_cpg_obs",
                     "compute_cpg_signal_composition", "_read_chrom_sizes", "_export_qc"),
    }
    fingerprints = {}
    for stage, names in roots.items():
        selected_roots = tuple("alopex.downstream.qc." + name for name in names)
        libraries = versions
        if stage == "products":
            selected_roots += ("alopex.downstream.diagnostics.filter_reasons",)
            libraries = {**versions, "backports.zstd": importlib.metadata.version("backports.zstd"),
                         "libzstd": zstd.zstd_version}
        selected = sources.dependencies(
            selected_roots,
            exclude=("alopex.downstream.qc._algorithm_fingerprints",),
        )
        payload = json.dumps({"code": selected, "libraries": libraries}, sort_keys=True)
        fingerprints[stage] = hashlib.sha256(payload.encode()).hexdigest()
    return fingerprints


def compute_cpg_signal_composition(
    df_qc: pd.DataFrame,
    cache_path: Path,
    composition_stats: pd.DataFrame,
    empty_cpg_sample_ids: set[str],
    *,
    before_commit=None,
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
        comp[column] = matrices.exact_integers(comp[column], np.int64)
        if (comp[column] < 0).any():
            raise ValueError(f"Negative CpG coverage in {column}")
    empty_mask = comp["Sample_ID"].isin(empty_cpg_sample_ids)
    if (comp.loc[~empty_mask, "total_cov"] == 0).any():
        raise ValueError("Non-empty CpG composition has no coverage")
    if (comp.loc[empty_mask, "total_cov"] != 0).any():
        raise ValueError("Empty CpG composition contains coverage")
    if not comp["Sample_CpG_Signals"].equals(comp["host_cov"]):
        raise ValueError("Sample_CpG_Signals differs from host CpG coverage")

    # 逐类扣减，避免有符号整型加总回绕；只在计算比例时使用浮点数。
    remaining = comp["total_cov"]
    denominator = remaining.replace(0, np.nan)
    base = comp[["Sample_ID", "Sample_CpG_Signals"]].copy()
    for column, rate in COMPOSITION_RATES.items():
        if (comp[column] > remaining).any():
            raise ValueError("Classified CpG coverage exceeds total coverage")
        remaining = remaining - comp[column]
        base[rate] = (comp[column] / denominator).fillna(0.0) * 100.0
    base["Signal_Composition_Rule"] = "primary_nuclear_chromosomes"
    comp_out = base[list(COMPOSITION_CACHE_COLUMNS)].copy()
    _atomic_write_csv(comp_out, cache_path, before_commit=before_commit)

    comp_merge = comp_out.drop(
        columns=["Signal_Composition_Rule", "Sample_CpG_Signals"], errors="ignore"
    )
    return df_qc.merge(
        comp_merge, on="Sample_ID", how="left", validate="one_to_one"
    )


def _gini_from_sorted_positive(x_sorted: np.ndarray) -> float:
    n = x_sorted.size
    if n == 0:
        return float("nan")
    # 有序整数窗的成对差之和。单细胞总覆盖已核验，分子与分母仅在最后相除时舍入。
    denominator = n * int(x_sorted.sum())
    weights = np.arange(1 - n, n, 2, dtype=np.int64)
    if denominator <= np.iinfo(np.int64).max:
        numerator = int(weights @ x_sorted)
    else:
        numerator = sum(int(weight) * int(depth) for weight, depth in zip(weights, x_sorted, strict=True))
    return numerator / denominator


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
            dtype={"position_coarse": str, "count": str}, keep_default_na=False, float_precision="round_trip",
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
        positions = matrices.exact_integers(prof["position_coarse"], np.int32)
        key_frame = pd.DataFrame(
            {"Sample_ID": sample_values, "position_coarse": positions}
        )
        if key_frame.duplicated(["Sample_ID", "position_coarse"]).any():
            return None
        fractions = pd.to_numeric(prof["meth_frac"], errors="raise").to_numpy(
            dtype=float
        )
        counts = matrices.exact_integers(prof["count"], np.int64)
        if (
            not np.isfinite(fractions).all()
            or (fractions < 0).any()
            or (fractions > 1).any()
            or (counts < 0).any()
            or np.any((counts == 0) != (fractions == 0))
        ):
            return None
        prof["Sample_ID"] = sample_values
        prof["position_coarse"] = positions
        prof["meth_frac"] = fractions
        prof["count"] = counts
        prof["CloneID"] = clone_values
        return prof
    except Exception:
        return None


@njit(cache=True, nogil=True)
def _coverage_bins(indices, coverage, chrom_offsets, bin_offsets, bin_size, base, primary):
    """按主核染色体分箱，以整数累加有观测位点的覆盖深度。"""
    values = np.empty(min(len(indices), bin_offsets[-1]), dtype=np.int64)
    chrom = 0
    previous = -1
    used = 0
    for k, index in enumerate(indices):
        while chrom + 1 < len(chrom_offsets) and index >= chrom_offsets[chrom + 1]:
            chrom += 1
        if chrom + 1 >= len(chrom_offsets):
            raise ValueError("CpG index is outside the reference")
        pos = index - chrom_offsets[chrom] - base
        if pos < 0:
            raise ValueError("CpG index is outside the configured coordinate base")
        if not primary[chrom]:
            continue
        target = bin_offsets[chrom] + pos // bin_size
        if target != previous:
            if target < previous:
                raise ValueError("CpG bins are not ordered")
            values[used] = coverage[k]
            used += 1
            previous = target
        else:
            values[used - 1] += coverage[k]
    return values[:used]


def _validate_tss_sampling(fraction: float, seed: int) -> None:
    if (isinstance(fraction, (bool, np.bool_))
            or not isinstance(fraction, (int, float, np.integer, np.floating))
            or not np.isfinite(fraction) or not 0 < fraction <= 1):
        raise ValueError("TSS sample fraction must be finite and in (0, 1]")
    if isinstance(seed, (bool, np.bool_)) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("TSS sample seed must be a nonnegative integer")


def _select_tss_cells(df_qc: pd.DataFrame, fraction: float, seed: int) -> tuple[pd.DataFrame, dict]:
    _validate_tss_sampling(fraction, seed)
    sample_ids = df_qc["Sample_ID"].astype(str).tolist()
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("TSS sampling requires unique Sample_ID values")
    clones = df_qc["CloneID"]
    if clones.isna().any() or clones.astype(str).str.strip().eq("").any():
        raise ValueError("TSS sampling requires nonempty CloneID values")
    selected_ids = []
    clone_counts = []
    for clone, group in df_qc.groupby("CloneID", sort=True):
        ids = group["Sample_ID"].astype(str).tolist()
        count = len(ids) if len(ids) < 10 or fraction == 1 else min(100, int(np.ceil(len(ids) * fraction)))
        ranked = sorted(ids, key=lambda sid: (
            hashlib.sha256(json.dumps([int(seed), sid], ensure_ascii=True).encode()).hexdigest(), sid
        ))
        selected_ids.extend(ranked[:count])
        clone_counts.append({"clone_id": str(clone), "total_cells": len(ids), "selected_cells": count})
    selected_ids.sort()
    selected = df_qc.set_index("Sample_ID", drop=False).loc[selected_ids].reset_index(drop=True)
    return selected, {
        "method": "clone_stratified_sha256_small10_cap100_v1",
        "fraction": float(fraction),
        "seed": int(seed),
        "total_cells": len(sample_ids),
        "selected_cells": len(selected_ids),
        "sample_ids": selected_ids,
        "clones": clone_counts,
    }





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
    source_qc.loc[source_qc["dedup_policy"].isin([tag.dedup_policy for tag in DD_CONSENSUS.values()]), "Duplicate_Pair_Rate%"] = float("nan")
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


def _export_qc(
    paths: RunPaths,
    df_qc: pd.DataFrame,
    pipeline_qc: pd.DataFrame,
    manifest_sample_ids: list[str],
    identity: Mapping[str, object],
    sealed: SealedQCContext,
    *, before_commit=None,
) -> None:
    assert_current = before_commit or sealed.assert_current
    with _qc_write_lock(paths.qc_results_dir):
        if not df_qc.empty:
            missing_qc_info_columns = [
                column for column in QC_INFO_COLUMNS if column not in df_qc.columns
            ]
            if missing_qc_info_columns:
                raise ValueError(
                    "DNAme QC information is missing column(s): "
                    + ", ".join(missing_qc_info_columns)
                )
            computed_columns = [
                "Sample_ID",
                "CloneID",
                "Native_Mapping%",
                "Duplicate_Pair_Rate%",
                "Final_Pair_Yield%",
                "Gini_Index",
            ]
            sample_id_values = df_qc["Sample_ID"].astype(str)
            # Gini只定义于宿主主核覆盖；对照/线粒体/辅助序列有信号仍可无此分母。
            primary_empty_mask = df_qc["Sample_Unique_CpG_Sites"].eq(0)
            umi_consensus_sample_ids = set(
                pipeline_qc.loc[
                    pipeline_qc["dedup_policy"].isin([tag.dedup_policy for tag in DD_CONSENSUS.values()]),
                    "sample_id",
                ].astype(str)
            )
            for column in computed_columns[2:]:
                values = pd.to_numeric(df_qc[column], errors="raise").to_numpy(dtype=float)
                nan_mask = None
                if column == "Gini_Index":
                    nan_mask = primary_empty_mask
                elif column == "Duplicate_Pair_Rate%":
                    nan_mask = sample_id_values.isin(umi_consensus_sample_ids)
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

            df_qc = df_qc[list(QC_INFO_COLUMNS)]

            df_qc = _require_exact_sample_ids(
                df_qc, manifest_sample_ids, "Final QC output"
            )
            def before_csv_commit():
                assert_current()
                _invalidate_completion(paths.input_identity)

            _atomic_write_csv(df_qc, paths.qc_info, before_commit=before_csv_commit)
            identity["qc_table"] = _content_binding(paths.qc_info)
            _atomic_write_json(identity, paths.input_identity, before_commit=assert_current)


@njit(cache=True, nogil=True)
def _summarize_cpg_row(indices, methyl_counts, depths, offsets, categories, window_offsets, base):
    observed = np.zeros(5, np.int64)
    coverage = np.zeros(5, np.int64)
    fractions = np.zeros(5, np.float64)
    corrections = np.zeros(5, np.float64)
    mc = np.zeros(window_offsets[-1], np.int64)
    cov = np.zeros_like(mc)
    chrom = 0
    for i in range(len(indices)):
        pos = indices[i]
        while pos >= offsets[chrom+1]:
            chrom += 1
        depth = depths[i]
        methyl = methyl_counts[i]
        group = categories[chrom]
        observed[group] += 1
        coverage[group] += depth
        y = methyl / depth - corrections[group]
        total = fractions[group] + y
        corrections[group] = (total - fractions[group]) - y
        fractions[group] = total
        if window_offsets[chrom+1] != window_offsets[chrom]:
            window = window_offsets[chrom] + (pos - offsets[chrom] - base) // 1_000_000
            mc[window] += methyl
            cov[window] += depth
    return observed, coverage, fractions, mc, cov


class _SinglePassCells:
    """按 sealed 细胞身份有界并行。每个 worker 只解码一次并独立算完该行。"""

    def __init__(self, sources, sizes, base, representation, workers):
        self.obs_names = pd.Index(sorted(sources))
        self.n_obs = len(self.obs_names)
        self.reference = matrices.cpg_reference(sizes, base)
        self.uns = {"cpg_coordinate_base": base, "source_cpg_representation": representation}
        self.sources = sources
        self.workers = _bounded_cpg_workers(workers, self.n_obs)

    def map_cells(self, executor, calculate):
        def process(sid):
            path = self.sources[sid]
            with _managed_cpg_text(path) as handle:
                row = matrices.read_cpg_counts(handle, path, self.reference, self.uns["cpg_coordinate_base"])
            return calculate(sid, row)

        if self.workers == 1:
            yield from map(process, self.obs_names)
            return
        ids = iter(self.obs_names)
        # 较慢的前序细胞等待写入时，worker 可提前处理下一行，只保留有界的待提交结果。
        pending = deque(executor.submit(process, sid) for sid in islice(ids, 2 * self.workers))
        try:
            while pending:
                result = pending.popleft().result()
                sid = next(ids, None)
                if sid is not None:
                    pending.append(executor.submit(process, sid))
                yield result
        finally:
            for future in pending:
                future.cancel()


def _calculate_cpg_products(reader, cfg, selected, zero_mapped, bed, single_tmp, counts_tmp):
    base = reader.uns["cpg_coordinate_base"]
    reference = reader.reference
    names = reference.reference_seq_name.tolist()
    lengths = reference.reference_bp_length.to_numpy(dtype=np.int64)
    offsets = np.r_[0, np.cumsum(lengths+base)]
    host_contigs = set(primary_analysis_contigs(cfg.genome))
    host_lengths = {name: length for name, length in zip(names, lengths, strict=True) if name in host_contigs}
    classes = [_chromosome_class(name, host_contigs) for name in names]
    category_order = ("host", "lambda", "pUC19", "mtDNA", "auxiliary")
    categories = np.asarray([category_order.index(group) for group in classes], dtype=np.int64)
    primary = categories == 0
    autosomes = host_contigs - {"chrX", "chrY"}
    is_auto = [counts_tmp is not None and name in autosomes
               for name in names]
    windows = np.r_[0, np.cumsum(np.where(is_auto, (lengths+999999)//1000000, 0))]
    bins = np.r_[0, np.cumsum((lengths + cfg.gini_bin_size - 1)//cfg.gini_bin_size)]
    tss_clones = dict(selected[["Sample_ID", "CloneID"]].itertuples(index=False, name=None))
    profiles, chromosome, statistics, gini = [], [], [], []
    mc = np.zeros((reader.n_obs, windows[-1]), dtype=np.int64)
    cov = np.zeros_like(mc)
    empty = set()

    def calculate(sid, row):
        is_empty = not len(row.indices)
        if not is_empty and sid in zero_mapped:
            raise ValueError(f"Zero-mapped cells must have an empty current CpG result: {sid}")
        sites, depth, fractions, cell_mc, cell_cov = _summarize_cpg_row(
            row.indices, row.methyl, row.coverage, offsets, categories, windows, base)
        means = [fractions[i]/sites[i]*100 if sites[i] else np.nan for i in range(3)]
        stats = (sid, sites[0], means[0], sites[1], means[1], sites[2], means[2],
                 depth[1], depth[2], depth[3], depth[0], depth.sum(), depth[0])
        counts = _coverage_bins(row.indices, row.coverage, offsets, bins, cfg.gini_bin_size, base, primary)
        counts.sort()
        chrom = list(diagnostics.chromosome_row(reference, offsets, row, sid, classes))
        profile = None
        if tss_index is not None and sid in tss_clones:
            profile = diagnostics.tss_profile_row(reader, row, sid, tss_clones[sid], tss_index)
        single = matrices.single_cpg_row(reader, row, single_index, n_vars)
        return single, stats, _gini_from_sorted_positive(counts), chrom, profile, cell_mc, cell_cov, is_empty

    with ThreadPoolExecutor(max_workers=reader.workers) as executor:
        tss_index = tss_future = None
        if cfg.tss_bed:
            if reader.workers > 1:
                tss_future = executor.submit(diagnostics.tss_intervals, cfg.tss_bed, reference, host_contigs)
            else:
                tss_index = diagnostics.tss_intervals(cfg.tss_bed, reference, host_contigs)
        with matrices.single_cpg_writer(reader, bed, single_tmp, host_lengths) as (writer, single_index):
            if tss_future is not None:
                tss_index = tss_future.result()
            n_vars = writer.n_vars
            with contextlib.closing(reader.map_cells(executor, calculate)) as cells:
                for start, (single, stats, score, chrom, profile, cell_mc, cell_cov, is_empty) in enumerate(cells):
                    writer.append(single)
                    statistics.append(stats)
                    gini.append(score)
                    chromosome.extend(chrom)
                    if profile is not None:
                        profiles.append(profile)
                    mc[start], cov[start] = cell_mc, cell_cov
                    if is_empty:
                        empty.add(reader.obs_names[start])
    if len(statistics) != reader.n_obs:
        raise RuntimeError("Incomplete single-pass CpG products")
    if counts_tmp is not None:
        regions = [f"{chrom}:{start+base}-{min(start+1000000,int(length))+base}"
                   for chrom, length, enabled in zip(names, lengths, is_auto, strict=True) if enabled
                   for start in range(0, length, 1000000)]
        with counts_tmp.open("xb") as handle:
            np.savez(handle, mc=mc, cov=cov, ids=np.asarray(reader.obs_names, dtype=str),
                     regions=np.asarray(regions, dtype=str))
    stats = pd.DataFrame(statistics, columns=[*MAPPING_STAT_COLUMNS,*COMPOSITION_COV_COLS])
    return stats, pd.Series(gini,index=reader.obs_names), pd.DataFrame(chromosome), diagnostics.tss_profile_frame(profiles), empty


def _single_cpg_obs(frame, pipeline_qc, reasons):
    """将同一次 QC 的筛选结果和 sealed 细胞身份保存到 SingleCpG。"""
    obs = frame[list(QC_INFO_COLUMNS)].merge(reasons, on=list(QC_INFO_COLUMNS[:3]),
                                            validate="one_to_one", sort=False)
    obs = obs.set_index("Sample_ID", drop=False)
    cells = pipeline_qc.set_index("sample_id").loc[obs.index]
    obs["DNA_Barcode"] = cells["demux_dna_barcode"].to_numpy()
    obs["Cell_Order"] = cells["cell_order"].to_numpy()
    obs["sample"] = obs["Sample_ID"]
    obs["clone_id"] = obs["CloneID"]
    obs["nCG"] = obs["Sample_Unique_CpG_Sites"]
    obs["met_rate"] = obs["Sample_Meth_CpG_Rate%"] / 100
    obs.index.name = None
    for name in obs.select_dtypes(include=["string"]):
        obs[name] = obs[name].astype(object)
    return obs


def run_pipeline(
    paths: RunPaths, cfg: RunConfig, workers: int, recompute_cpg_stats: bool,
    recompute_gini: bool, recompute_tss: bool, *, sealed_context: SealedQCContext | None = None,
) -> None:
    """由同一次 CpG 解码生成完整下游，完成身份在所有文件提交后写出。"""
    with _qc_write_lock(paths.qc_results_dir):
        _run_pipeline_impl(paths, cfg, workers, recompute_cpg_stats or recompute_gini or recompute_tss,
                           sealed_context=sealed_context)


def _run_pipeline_impl(paths, cfg, workers, recompute, *, sealed_context=None):
    overall_t0 = time.perf_counter()
    stage_start = overall_t0

    def stage(number, title):
        nonlocal stage_start
        if number != 1:
            _log(f"阶段完成；耗时：{time.perf_counter()-stage_start:.2f}s")
        stage_start = time.perf_counter()
        _log(f"阶段 {number}/6：{title}")

    stage(1, "核对当前完整交付和 QC 输入")
    assert_loaded_sources_current(_LOADED_SOURCE_CODES)
    sealed = sealed_context or load_sealed_qc_context(paths.project_dir, expected_pipeline_root=_pipeline_dir())
    if cfg.genome != sealed.config["species"] or cfg.ref_fai != sealed.reference_fai:
        raise ValueError("QC species and reference FAI must match the sealed delivery")
    bed = _single_bed(sealed)
    check_references = _reference_guard(sealed.reference_path, cfg.ref_fai, bed, cfg.tss_bed,
        Path(str(bed)+".reference.json"), Path(str(cfg.tss_bed)+".reference.json") if cfg.tss_bed else None,
        before_commit=sealed.assert_current)
    single_reference = validate_downstream_reference(sealed.reference_path, bed, "single-cpg")
    tss_reference = validate_tss_reference(sealed.reference_path, cfg.tss_bed) if cfg.tss_bed else None
    records = build_qc_records(paths.project_dir, sealed_context=sealed)
    if not records:
        raise ValueError("QC requires at least one cell")
    pipeline_qc = pd.DataFrame.from_records(records)
    _log(f"QC 输入：{len(records)} 个细胞；物种：{cfg.genome}；CPU预算：{workers}")
    sources = {record["sample_id"]: record["cpg_file"] for record in records}
    def assert_current():
        check_references()
        for record in records:
            sealed.assert_output_current(record["cpg_output_identity"]["path"])
    if len(sources)!=len(records):
        raise ValueError("Sealed CpG samples must be unique")
    ids = sorted(sources)
    selected, sampling = _select_tss_cells(
        pipeline_qc[["sample_id","demux_clone_id"]].set_axis(["Sample_ID","CloneID"],axis=1),
        cfg.tss_sample_fraction, cfg.tss_sample_seed)
    if cfg.tss_bed is None:
        selected = selected.iloc[:0]
    filter_parameters = dict(cfg.qc_filter_parameters)
    filter_uns = {"qc_filter_parameters": json.dumps(filter_parameters, sort_keys=True)}
    filter_output = paths.qc_results_dir / diagnostics.FILTER_INFO
    algorithms = _algorithm_fingerprints()
    identity = source_identity(paths.project_dir, records, sealed_context=sealed,
        extra_files={"reference_fai": cfg.ref_fai},
        context={"bed": single_reference["files"][bed.name]["sha256"],
                 "tss": tss_reference["files"][cfg.tss_bed.name]["sha256"] if tss_reference else None,
                 "gini_bin_size": cfg.gini_bin_size, "tss_sampling": sampling if cfg.tss_bed else None})
    identity["algorithms"] = algorithms
    identity["tss_sampling"] = sampling if cfg.tss_bed else None
    identity["qc_filter_parameters"] = filter_parameters
    droplet = sealed.config["analysis"]["protocol"] in DD_PROTOCOLS
    if not recompute:
        try:
            completed = _require_completed_qc(paths, sealed, algorithms=algorithms, verify_diagnostics=None)
        except ValueError:
            pass
        else:
            if completed.identity["identity_sha256"] == identity["identity_sha256"]:
                completed.assert_current()
                assert_current()
                if not droplet:
                    _invalidate_completion(paths.doublet_counts)
                if cfg.tss_bed is None:
                    _remove_skipped_tss(paths)
                if completed.identity["qc_filter_parameters"] != filter_parameters:
                    frame = _read_qc_table(paths.qc_info)
                    frame = _require_exact_sample_ids(frame, ids, "Filter QC")
                    reasons = diagnostics.filter_reasons(frame, **filter_parameters)
                    obs = _single_cpg_obs(frame, pipeline_qc, reasons)
                    schema = _require_current_single(paths, completed.identity)
                    current = dict(completed.identity, qc_filter_parameters=filter_parameters)
                    completed.assert_current()
                    assert_current()
                    _invalidate_completion(paths.input_identity)
                    _invalidate_completion(paths.single_schema)
                    matrices.update_obs(paths.single_cpg, obs, filter_uns)
                    _fsync_file(paths.single_cpg)
                    _atomic_write_csv(reasons, filter_output, before_commit=assert_current)
                    current["diagnostic_artifacts"][filter_output.name] = _content_binding(filter_output)
                    schema["h5ad"] = _published_file_binding(paths.single_cpg)
                    current["single_cpg_ctime_ns"] = paths.single_cpg.lstat().st_ctime_ns
                    schema["input"]["qc_filter_parameters"] = filter_parameters
                    _atomic_write_json(schema, paths.single_schema, before_commit=assert_current)
                    _atomic_write_json(current, paths.input_identity, before_commit=assert_current)
                    _log("复用已核验缓存中的 CpG 结果，已更新 HQ 筛选注释")
                    _log(f"QC 已完成；总耗时：{time.perf_counter()-overall_t0:.2f}s")
                    return
                completed.assert_current()
                assert_current()
                _log("已复用当前完整 QC、SingleCpG 和适用的 1 Mb 计数")
                _log(f"QC 已完成；总耗时：{time.perf_counter()-overall_t0:.2f}s")
                return

    stage(2, "准备同一交付的 CpG 和参考")
    base = _source_cpg_coordinate_base(sealed.cpg_representation)
    sizes = _read_chrom_sizes(cfg.ref_fai)
    reader = _SinglePassCells(sources, sizes, base, sealed.cpg_representation, workers)
    zero_mapped = {r["sample_id"] for r in records if r["backend_accepted_pairs"] == 0}
    temp_dir = project_temporary_parent(paths.single_cpg)
    with tempfile.TemporaryDirectory(prefix="cpg-products.", dir=temp_dir) as folder:
        work = Path(folder)
        single_tmp = work/paths.single_cpg.name
        counts_tmp = work/paths.doublet_counts.name if droplet else None
        stage(3, "单次解码，汇总并写出 SingleCpG")
        stats, gini, chromosome, tss, empty = _calculate_cpg_products(
            reader, cfg, selected, zero_mapped, bed, single_tmp, counts_tmp)
        stage(4, "按细胞身份合并 Pipeline QC")
        frame = _merge_pipeline_qc(pipeline_qc, stats[list(MAPPING_STAT_COLUMNS)])
        composition = stats[["Sample_ID",*COMPOSITION_COV_COLS]]
        stage(5, "整理 QC、信号组成及诊断表")
        frame = compute_cpg_signal_composition(frame, work/paths.composition_cache.name, composition, empty)
        frame["Gini_Index"] = frame.Sample_ID.map(gini)
        frame = _require_exact_sample_ids(frame, ids, "Single-pass QC")
        chromosome["CloneID"] = chromosome.Sample_ID.map(frame.set_index("Sample_ID").CloneID)
        chrom_output = paths.qc_results_dir/diagnostics.CHROMOSOME_INFO
        _atomic_write_csv(chromosome, work/chrom_output.name)
        if cfg.tss_bed is not None:
            _atomic_write_csv(tss, work/paths.tss_info.name)
        reasons = diagnostics.filter_reasons(frame, **filter_parameters)
        _atomic_write_csv(reasons, work / filter_output.name)
        matrices.update_obs(single_tmp, _single_cpg_obs(frame, pipeline_qc, reasons), filter_uns)
        with contextlib.closing(matrices.CSRReader(single_tmp)) as single:
            if list(single.obs_names) != ids:
                raise ValueError("SingleCpG row order differs from the sealed CpG input")
        stage(6, "提交 SingleCpG、整数计数与 QC 完成身份")
        assert_loaded_sources_current(_LOADED_SOURCE_CODES)
        assert_current()
        _invalidate_completion(paths.input_identity)
        _invalidate_completion(paths.single_schema)
        _fsync_file(single_tmp)
        assert_current()
        os.replace(single_tmp, paths.single_cpg)
        _fsync_directory(paths.single_cpg.parent)
        generation = uuid.uuid4().hex
        identity["cpg_generation"] = generation
        identity["single_cpg_ctime_ns"] = paths.single_cpg.lstat().st_ctime_ns
        identity["empty_cpg_sample_ids"] = sorted(empty)
        if counts_tmp is not None:
            _fsync_file(counts_tmp)
            os.replace(counts_tmp, paths.doublet_counts)
            _fsync_directory(paths.doublet_counts.parent)
            identity["doublet_counts"] = _published_file_binding(paths.doublet_counts)
            identity["doublet_counts_ctime_ns"] = paths.doublet_counts.lstat().st_ctime_ns
        else:
            _invalidate_completion(paths.doublet_counts)
        destinations = [paths.composition_cache, chrom_output, filter_output]
        if cfg.tss_bed is not None:
            destinations.append(paths.tss_info)
            if tss.empty:
                _remove_tss_plot(paths)
        else:
            _remove_skipped_tss(paths)
        for destination in destinations:
            os.replace(work/destination.name, destination)
            _fsync_directory(destination.parent)
        identity["diagnostic_artifacts"] = {path.name: _content_binding(path) for path in destinations}
        _atomic_write_json({"input": {"cpg_generation": generation, "delivery_id": sealed.delivery_id,
                                      "cpg_coordinate_base": base,
                                      "cpg_interval_convention": matrices.SINGLE_CPG_INTERVAL_CONVENTION,
                                      "bed_sha256": identity["context"]["bed"],
                                      "qc_filter_parameters": filter_parameters,
                                      "source_cpg_representation": sealed.cpg_representation},
                            "h5ad": _published_file_binding(paths.single_cpg)},
                           paths.single_schema, before_commit=assert_current)
        _export_qc(paths, frame, pipeline_qc, ids, identity, sealed, before_commit=assert_current)
    _log(f"阶段完成；耗时：{time.perf_counter()-stage_start:.2f}s")
    _log(f"QC 已完成；总耗时：{time.perf_counter()-overall_t0:.2f}s")


def parse_args() -> argparse.Namespace:
    return _add_cli_args(argparse.ArgumentParser(prog="python -m alopex.downstream.qc", description=__doc__)).parse_args()


def _log_config_summary(paths: RunPaths, cfg: RunConfig):
    _log(f"DNA QC 开始；物种：{cfg.genome}；项目：{paths.project_dir}")
    _log(f"输入：{paths.cpg_dir}；输出：{paths.qc_results_dir}")
    _log(f"参考索引：{cfg.ref_fai}；TSS：{cfg.tss_bed if cfg.tss_bed is not None else '跳过'}")


def main() -> None:
    args = parse_args()
    project_dir = _resolve_project_dir(Path(args.project_dir))
    sealed = load_sealed_qc_context(
        project_dir, expected_pipeline_root=_pipeline_dir()
    )
    if args.expected_delivery_id is not None and sealed.delivery_id != args.expected_delivery_id:
        raise ValueError("Delivery changed; rerun the Notebook from its first cell")
    paths = build_paths(project_dir)
    cfg = build_config(
        sealed,
        tss_bed_override=Path(args.tss_bed) if args.tss_bed else None,
        skip_tss=bool(args.skip_tss),
        tss_sample_fraction=args.tss_sample_fraction,
        tss_sample_seed=args.tss_sample_seed,
        qc_filter_parameters={name: getattr(args, name) for name in QC_FILTER_DEFAULTS},
    )

    _log_config_summary(paths, cfg)
    run_pipeline(
        paths=paths,
        cfg=cfg,
        workers=args.workers,
        recompute_cpg_stats=bool(args.recompute_cpg_stats),
        recompute_gini=bool(args.recompute_gini),
        recompute_tss=bool(args.recompute_tss),
        sealed_context=sealed,
    )


def load_cpg_density_frame(
    paths: NotebookPaths,
) -> pd.DataFrame:
    """从原始来源组装 CpG_Density 输入帧，不依赖 20 列 QC 表之外的衍生列。

    pair 计数来自 sealed final manifest（build_qc_records 唯一加载链），
    Sample_CpG_Signals 来自 composition cache，CloneID 与 unique CpG sites
    来自当前 QC 表；三个来源须属于 Notebook 捕获的同一次已完成 QC。
    """
    sealed = load_sealed_qc_context(paths.project_dir, expected_pipeline_root=paths.pipeline_root)
    if sealed.delivery_id != paths.delivery_id:
        raise ValueError("Delivery changed; rerun the Notebook from its first cell")
    run_paths = build_paths(paths.project_dir)
    if run_paths.qc_results_dir != paths.qc_results_dir:
        raise ValueError("Notebook QC directory differs from the current delivery")
    assert_qc_current = _require_completed_qc(
        run_paths, sealed, verify_diagnostics=(run_paths.composition_cache.name,),
    ).assert_current
    qc = _read_qc_table(run_paths.qc_info)
    required = ("Sample_ID", "CloneID", "Sample_Unique_CpG_Sites")
    missing = [column for column in required if column not in qc.columns]
    if missing:
        raise ValueError(
            "CpG density inputs are missing QC column(s): " + ", ".join(missing)
        )
    qc["Sample_ID"] = qc["Sample_ID"].astype(str).str.strip()

    records = build_qc_records(paths.project_dir, sealed_context=sealed)
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

    composition_cache_path = run_paths.composition_cache
    composition = pd.read_csv(
        composition_cache_path, converters={"Sample_ID": str}, dtype={"Sample_CpG_Signals": str},
    )
    for column in ("Sample_ID", "Sample_CpG_Signals"):
        if column not in composition.columns:
            raise ValueError(f"CpG composition cache is missing column: {column}")
    composition["Sample_ID"] = composition["Sample_ID"].astype(str).str.strip()
    composition["Sample_CpG_Signals"] = _nonnegative_counts(composition["Sample_CpG_Signals"])

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
    assert_qc_current()
    return frame


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
    return str(value) if value is not None else ""


def _validate_notebook_pipeline_root(candidate: Path) -> Path:
    candidate = candidate.expanduser().resolve()
    processor = candidate / "core" / "alopex" / "downstream" / "qc.py"
    launcher = candidate / "core" / "run_pipeline.sh"
    if not processor.is_file() or not launcher.is_file():
        raise FileNotFoundError(f"Not a usable Alopex: {candidate}")
    return candidate


def _notebook_project_root(dna_project_root: str | Path | None) -> Path:
    option = _notebook_root_option(dna_project_root, "DNA_PROJECT_ROOT")
    start = Path(option).expanduser().resolve() if option else Path.cwd().resolve()
    for base in ([start] if option else [start, *start.parents]):
        candidate = base / "03_results" / "run_manifest.json"
        state = stat_if_present(candidate)
        if state is None:
            if stat_if_present(base / "00_config/config.yaml") is not None:
                raise FileNotFoundError(f"Alopex project has no sealed delivery: {candidate}")
            continue
        if not _qc_lock_stat.S_ISREG(state.st_mode):
            raise ValueError(f"Delivery manifest must be a regular file: {candidate}")
        return base
    if option:
        raise FileNotFoundError(f"DNA_PROJECT_ROOT 没有 sealed delivery: {start}")
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
    source_root。alopex.downstream.qc 位于 Pipeline 树内，Notebook 第一格
    的定位 bootstrap 之外的完整校验以本函数为唯一实现。
    """
    assert_loaded_sources_current(_LOADED_SOURCE_CODES)
    project_dir = _notebook_project_root(dna_project_root)
    manifest = load_run_manifest(project_dir, validate_outputs=False)
    pipeline_root = _notebook_pipeline_root(dna_pipeline_root, manifest)
    delivery_id = _manifest_delivery_id(manifest)
    qc_results_dir = resolve_qc_results_dir(project_dir)
    return NotebookPaths(
        project_dir=project_dir,
        pipeline_root=pipeline_root,
        delivery_id=delivery_id,
        qc_results_dir=qc_results_dir,
        plots_dir=_project_directory(project_dir, Path("06_downstream/plots"), create=False),
        qc_info=qc_results_dir / "DNAme_QC_Information.csv",
        tss_info=qc_results_dir / "TSS_Profile_Information.csv",
        composition_cache=qc_results_dir / ".cache" / "CpG_Signal_Composition.csv",
    )


def run_qc_processor(
    paths: NotebookPaths,
    *,
    tss_bed: str | Path | None = None,
    skip_tss: bool = False,
    tss_sample_fraction: float = 0.1,
    tss_sample_seed: int = 0,
    recompute_cpg_stats: bool = False,
    recompute_gini: bool = False,
    recompute_tss: bool = False,
    qc_filter_parameters: Mapping[str, float | None] | None = None,
) -> None:
    """以独立子进程运行 alopex.downstream.qc 计算 main 并在失败时终止 Notebook。

    worker 与 native 线程共同受进程 affinity 和 Slurm allocation 的 CPU 预算约束；
    重算开关透传 CLI，强制重新解码 sealed CpG；默认复用同代完整结果。
    TSS 默认每个 CloneID 抽检 10%（不足 10 个全检，最多 100 个）。
    比例设为 1.0 时全量检查，修改参数触发整代重算。
    """
    assert_loaded_sources_current(_LOADED_SOURCE_CODES)
    _validate_tss_sampling(tss_sample_fraction, tss_sample_seed)
    if skip_tss and tss_bed is not None:
        raise ValueError("skip_tss cannot be combined with an explicit TSS bed")
    processor = paths.pipeline_root / "core" / "alopex" / "downstream" / "qc.py"
    if not processor.is_file():
        raise FileNotFoundError(f"Downstream QC processor is missing: {processor}")
    worker_count = str(effective_cpu_count())
    command = [
        sys.executable,
        "-m", "alopex.downstream.qc",
        "--project-dir",
        str(paths.project_dir),
        "--expected-delivery-id",
        paths.delivery_id,
        "--workers",
        worker_count,
        "--tss-sample-fraction",
        str(tss_sample_fraction),
        "--tss-sample-seed",
        str(tss_sample_seed),
    ]
    if tss_bed is not None:
        command.extend(["--tss-bed", str(tss_bed)])
    elif skip_tss:
        command.append("--skip-tss")
    for name, value in (qc_filter_parameters or {}).items():
        if name not in QC_FILTER_DEFAULTS:
            raise ValueError(f"Unknown HQ filter parameter: {name}")
        if value is not None or name != "mapping_threshold":
            command.extend(("--" + name.replace("_", "-"), str(value)))
    for flag, enabled in (
        ("--recompute-cpg-stats", recompute_cpg_stats),
        ("--recompute-gini", recompute_gini),
        ("--recompute-tss", recompute_tss),
    ):
        if enabled:
            command.append(flag)
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(paths.pipeline_root / "core")
    # 使用本 Pipeline 已构建的解码工具，避免继承 Jupyter 服务的其他 Conda 环境。
    environment["PATH"] = str(paths.pipeline_root / "conda/current/control/bin") + os.pathsep + environment.get("PATH", "")
    for variable in _THREAD_LIMIT_ENV_VARS:
        environment[variable] = worker_count
    from ..presentation import run_displayed
    log_dir = _project_directory(paths.project_dir, Path("04_logs/downstream"), create=False)
    exit_code = run_displayed(
        command, mode="qc", title=f"DNA 甲基化 QC · {paths.project_dir.name}",
        log_dir=log_dir, project=paths.project_dir, budget=int(worker_count),
        env=environment, notebook=True,
    )
    if exit_code:
        raise RuntimeError(f"QC Processor failed with exit code {exit_code}; see 04_logs/downstream/")


def load_notebook_qc_frame(
    paths: NotebookPaths,
) -> tuple[pd.DataFrame, dict[str, str] | None]:
    """加载 Processor 同代 QC 与筛选结果并初始化绘图主题与 CloneID 配色。"""
    sealed = load_sealed_qc_context(paths.project_dir, expected_pipeline_root=paths.pipeline_root)
    if sealed.delivery_id != paths.delivery_id:
        raise ValueError("Delivery changed; rerun the Notebook from its first cell")
    run_paths = build_paths(paths.project_dir)
    if run_paths.qc_results_dir != paths.qc_results_dir:
        raise ValueError("Notebook QC directory differs from the current delivery")
    assert_qc_current = _require_completed_qc(
        run_paths, sealed, verify_diagnostics=(diagnostics.FILTER_INFO,),
    ).assert_current
    _project_directory(paths.project_dir, Path("06_downstream/plots"), create=True)
    from .plotting import setup_nature_style, normalize_clone_id_value, generate_clone_color_map

    setup_nature_style()
    df_qc = _read_qc_table(run_paths.qc_info)
    reasons = pd.read_csv(run_paths.qc_results_dir / diagnostics.FILTER_INFO,
                          converters=dict.fromkeys(QC_INFO_COLUMNS[:3], str))
    df_qc = df_qc.merge(reasons, on=list(QC_INFO_COLUMNS[:3]), how="left", validate="one_to_one")
    print(f"✅ 已加载 QC 数据: {len(df_qc)} 个样本")
    clone_color_map: dict[str, str] | None = None
    if not df_qc.empty and "CloneID" in df_qc.columns:
        df_qc["CloneID"] = df_qc["CloneID"].map(normalize_clone_id_value)
        unique_clones = sorted(df_qc["CloneID"].dropna().unique().tolist())
        clone_color_map = generate_clone_color_map(unique_clones)
    assert_qc_current()
    return df_qc, clone_color_map


def load_notebook_tss_frame(paths: NotebookPaths) -> pd.DataFrame:
    """读取当前已完成 QC 的 TSS 缓存，独立重跑绘图单元也须核对交付。"""
    sealed = load_sealed_qc_context(paths.project_dir, expected_pipeline_root=paths.pipeline_root)
    if sealed.delivery_id != paths.delivery_id:
        raise ValueError("Delivery changed; rerun the Notebook from its first cell")
    run_paths = build_paths(paths.project_dir)
    if run_paths.qc_results_dir != paths.qc_results_dir:
        raise ValueError("Notebook QC directory differs from the current delivery")
    completed = _require_completed_qc(run_paths, sealed, verify_diagnostics=(run_paths.tss_info.name,))
    if completed.identity["tss_sampling"] is None:
        raise ValueError("TSS was skipped in the completed QC run; set SKIP_TSS=False and rerun the QC Processor")
    frame = _validated_tss_cache(
        run_paths.tss_info, set(completed.sample_ids),
        set(completed.identity["empty_cpg_sample_ids"]),
    )
    if frame is None:
        raise ValueError("TSS cache is incomplete or invalid; rerun the QC Processor")
    completed.assert_current()
    return frame


def _load_qc_diagnostic(paths, name):
    """读取完成凭证绑定的诊断CSV；不接纳同名旧表或读取期间的替换。"""
    sealed = load_sealed_qc_context(paths.project_dir, expected_pipeline_root=paths.pipeline_root)
    if sealed.delivery_id != paths.delivery_id:
        raise ValueError('Delivery changed; rerun the Notebook from its first cell')
    run_paths = build_paths(paths.project_dir)
    if run_paths.qc_results_dir != paths.qc_results_dir:
        raise ValueError("Notebook QC directory differs from the current delivery")
    completed = _require_completed_qc(run_paths, sealed, verify_diagnostics=(name,))
    path = run_paths.qc_results_dir / name
    frame = pd.read_csv(path, converters={'Sample_ID': str, 'CloneID': str, 'chrom': str}, float_precision='round_trip')
    completed.assert_current()
    return frame


def load_chromosome_qc_frame(paths):
    """加载全部细胞的同代逐染色体计数/甲基化率。"""
    return _load_qc_diagnostic(paths, diagnostics.CHROMOSOME_INFO)


if __name__ == "__main__":
    main()
