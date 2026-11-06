"""汇总、MultiQC、snapshot、sealed交付与发布后回收。"""
from __future__ import annotations

import stat
import csv
import datetime as dt
import glob
import hashlib
import io
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import uuid
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from .config import (
    CellOutputPaths,
    DROPLET_DESIGN_RESOURCE,
    MULTIQC_BISCUIT_SOURCE_SUFFIXES,
    PUBLIC_RESULT_DIR_NAMES,
    ProjectLayout,
    ProjectSampleRow,
    SAMPLE_ID_RE,
    _parse_project_config,
    backend_active_index_paths,
    canonical_raw_route_id,
    cell_artifacts,
    cpg_representation,
    default_dna_route,
    get_backend,
    get_protocol,
    load_project_config,
    logger,
    read_environment_identity,
    reference_identity,
    resolve_bismark_alignment,
    resolve_reference,
    validate_reference_component,
)
from .inputs import (
    INITIAL_CELL_MANIFEST_FIELDS,
    _cell_order_key,
    _parse_project_sample_manifest,
    discover_samples,
    manifest_raw_samples,
    map_discovered_samples_to_manifest,
    read_cell_identity_manifest,
    read_demux_generation,
    rna_fastq_collect_destination,
)
from .io import (
    project_temporary_parent,
    _make_tree_read_only,
    _read_json,
    _remove_path,
    _remove_stage,
    _safe_relative,
    _utc_now,
    atomic_write_bytes,
    atomic_write_json,
    canonical_json_bytes,
    sha256_bytes,
    sha256_file,
    stable_file_identity,
    write_tsv_atomic,
)
from .runtime import (
    MULTIQC_MAX_CELLS,
    effective_cpu_count,
)
from .science import (
    DROPLET_DEDUP_POLICY,
    _bismark_version_text,
    _high_cph_multiqc_payload,
    _parse_bismark_non_cpg_methylation_pct,
    _parse_cph_methylation_pct,
    _parse_cutadapt_metrics,
    _read_biscuit_native_mapping,
    _read_dupsifter_metrics,
    _read_high_cph_summary,
    build_metrics,
)


# 汇总与交付：逐细胞结果和报告

def _write_pair_level_multiqc(rows: list[dict[str, object]], out_json: Path) -> None:
    data: dict[str, dict[str, object]] = {}
    for row in rows:
        if row.get("backend_qc_status") != "PASS":
            continue
        sample = str(row["sample_id"])
        input_pairs = int(row["cutadapt_input_pairs"])
        trimmed_pairs = int(row["trimmed_pairs"])
        accepted_pairs = int(row["backend_accepted_pairs"])
        duplicate_pairs = int(row["duplicate_pairs"])
        postdedup_pairs = int(row["postdedup_pairs"])
        final_pairs = int(row["final_retained_pairs"])
        non_cpg = float(row["non_cpg_methylation_pct"])
        data[sample] = {
            "Backend": row["methylation_backend"],
            "Mapping_Policy": row["mapping_policy"],
            "Input_Pairs": input_pairs,
            "Trimmed_Pairs": trimmed_pairs,
            "Accepted_Pairs": accepted_pairs,
            "PostDedup_Pairs": postdedup_pairs,
            "Final_Retained_Pairs": final_pairs,
            "Trim_Retention_pct": (
                100.0 * trimmed_pairs / input_pairs if input_pairs else 0.0
            ),
            "Native_Mapping_pct": float(row["native_mapping_pct"]),
            "Native_Mapping_Unit": row["native_mapping_unit"],
            "Duplicate_Pair_Rate_pct": (None if row.get("protocol") == "droplet" else
                100.0 * duplicate_pairs / accepted_pairs if accepted_pairs else 0.0
            ),
            "Dedup_Policy": row["dedup_policy"],
            **{key: row[key] for key in ("cpg_read_observations", "cpg_umi_observations", "cpg_ambiguous_umi_sites") if key in row},
            "PostDedup_Retention_pct": (
                100.0 * postdedup_pairs / accepted_pairs if accepted_pairs else 0.0
            ),
            "Final_Pair_Yield_pct": (
                100.0 * final_pairs / input_pairs if input_pairs else 0.0
            ),
            "High_CpH_Flag_Rate_pct": (None if row["methylation_backend"] == "rastair" else 100.0 * float(row["high_cph_fraction"])),
            "High_CpH_Status": row["high_cph_role"],
            **{key: row[key] for key in ("lambda_cpg_observations", "lambda_false_positive_pct", "pUC19_cpg_observations", "pUC19_conversion_pct", "bam_primary_records") if key in row},
            "Non_CpG_Methylation_pct": non_cpg if math.isfinite(non_cpg) else None,
            **{key: row[key] for key in ("non_cpg_observations", "non_cpg_candidate_pairs", "non_cpg_oriented_pairs", "non_cpg_ambiguous_pairs", "snp_sites") if key in row},
        }
    payload = {
        "id": "dna_pipeline_pair_level_qc",
        "section_name": "Alopex pair-level QC",
        "description": (
            "Read-pair retention funnel plus backend-native mapping QC. BISCUIT native "
            "mapping is primary MAPQ>=40 individual reads; Bismark native mapping is "
            "unique concordant read pairs. These mapping percentages use declared units "
            "and must not be treated as the same numerator definition. TAPS uses proper primary pairs "
            "with at least one MAPQ>=20 non-QC-failed mate; duplicate flags are excluded only during calling. "
            "TAPS Non-CpG is uncorrected CHG/CHH conversion in host pairs with unambiguous CpG-based orientation; "
            "SNPs/errors may contribute and unorientable pairs are excluded. "
            "TAPS high-CpH filtering is not applicable; cDNA is not excluded. Lambda assumes unmethylated "
            "DNA (CpG false positives); pUC19 assumes CpG-methylated DNA (conversion); missing coverage is not assessed."
        ),
        "plot_type": "table",
        "pconfig": {
            "id": "dna_pipeline_pair_level_qc_table",
            "title": "Canonical pair-level QC",
            "col1_header": "Cell",
        },
        "data": data,
    }
    atomic_write_json(out_json, payload)


def _write_cutadapt_multiqc(rows: list[dict[str, object]], out_json: Path) -> None:
    data: dict[str, dict[str, object]] = {}
    for row in rows:
        if row.get("backend_qc_status") != "PASS":
            continue
        input_pairs = int(row["cutadapt_input_pairs"])
        trimmed_pairs = int(row["trimmed_pairs"])
        data[str(row["sample_id"])] = {
            "Input_Pairs": input_pairs,
            "Trimmed_Pairs": trimmed_pairs,
            "Trim_Retention_pct": (
                100.0 * trimmed_pairs / input_pairs if input_pairs else 0.0
            ),
            "R1_Adapter_pct": float(row["r1_adapter_pct"]),
            "R2_Adapter_pct": float(row["r2_adapter_pct"]),
        }
    payload = {
        "id": "dna_pipeline_cutadapt_qc",
        "section_name": "Alopex Cutadapt R1/R2 QC",
        "description": (
            "Lightweight trimming QC parsed from the production Cutadapt JSON; "
            "no FASTQ is rescanned for presentation."
        ),
        "plot_type": "table",
        "pconfig": {
            "id": "dna_pipeline_cutadapt_qc_table",
            "title": "Cutadapt R1 / R2 QC",
            "col1_header": "Cell",
        },
        "data": data,
    }
    atomic_write_json(out_json, payload)


@dataclass(frozen=True)
class FinalManifestOptions:
    """保存已校验的汇总选项，不再次读取项目的可变配置。"""

    backend: str
    cpg_coordinate_base: int
    generate_snp: bool
    keep_final_bam: bool
    publish_rna_bam: bool


@dataclass(frozen=True)
class CellMetricJob:
    """向汇总 worker 传递当前 cell 的身份、路径和选项；只包含可序列化数据。"""

    cell_id: str
    protocol: str
    dna_reads: str
    paths: CellOutputPaths
    public_paths: CellOutputPaths
    options: FinalManifestOptions
    bismark_version: str


@dataclass(frozen=True)
class NativeCellMetrics:
    """保留原生指标口径，并用独立 read-pair 计数连接公共汇总。"""

    accepted_pairs: int
    rejected_pairs: int
    duplicate_pairs: int
    postdedup_pairs: int
    fields: dict[str, object]

    def manifest_fields(self) -> dict[str, object]:
        """展开 pair 漏斗与原生指标，不把 read 或 UMI 指标换成 pair 口径。"""
        return {**self.fields, "backend_accepted_pairs": self.accepted_pairs,
                "backend_rejected_pairs": self.rejected_pairs,
                "duplicate_pairs": self.duplicate_pairs, "postdedup_pairs": self.postdedup_pairs}


def write_final_sample_manifest(
    initial_cell_manifest: str,
    intermediate_dir: str,
    staged_results_dir: str,
    out_tsv: str,
    *,
    options: Mapping[str, object],
    published_results_dir: str | None = None,
    bismark: str | None = None,
    pair_qc_multiqc_json: str | None = None,
    cutadapt_qc_multiqc_json: str | None = None,
    workers: int = 1,
) -> None:
    """把已提交的 demux 身份与已启用的 per-cell 结果产物连接成最终 manifest。

    编排五个纯阶段：选项解析 → 期望 route 规划 → 已提交代际合并 →
    per-cell 指标核算 → 排序与落盘；阶段之间只通过显式参数传递数据。
    """
    intermediate = Path(intermediate_dir)
    demux_state = intermediate / "demux" / "state"
    results = Path(staged_results_dir)
    published_results = Path(published_results_dir) if published_results_dir else results
    opts = _parse_final_options(options)
    layout = ProjectLayout(intermediate, results, opts.backend)
    published_layout = ProjectLayout(intermediate, published_results, opts.backend)
    output = Path(out_tsv)

    initial_cells = read_cell_identity_manifest(initial_cell_manifest)
    cells: dict[str, dict[str, object]] = {
        cell_id: dict(row) for cell_id, row in initial_cells.items()
    }


    expected_runs = _plan_expected_runs(cells, demux_state)
    _merge_committed_generations(
        demux_state,
        initial_cell_manifest,
        cells,
        expected_runs,
    )
    _require_declared_routes(cells)
    _apply_cell_metrics(
        cells,
        layout=layout,
        published_layout=published_layout,
        options=opts,
        bismark=bismark,
        workers=workers,
    )
    _emit_final_tsv(cells, output, pair_qc_multiqc_json, cutadapt_qc_multiqc_json)


def _parse_final_options(options: Mapping[str, object]) -> FinalManifestOptions:

    backend = str(options.get("methylation_backend", "")).strip().lower()
    if backend not in {"biscuit", "bismark", "rastair"}:
        raise ValueError("methylation_backend must be biscuit, bismark or rastair")
    option_fields = {
        "generate_snp",
        "keep_final_bam",
        "publish_rna_bam",
        "methylation_backend",
        "cpg_coordinate_base",
    }
    if set(options) != option_fields:
        missing = sorted(option_fields - set(options))
        extra = sorted(set(options) - option_fields)
        details = []
        if missing:
            details.append("missing=" + ",".join(missing))
        if extra:
            details.append("unsupported=" + ",".join(extra))
        raise ValueError("Final manifest options do not match the current contract: " + "; ".join(details))
    active_options = dict(options)
    for key, value in active_options.items():
        if key in {"methylation_backend", "cpg_coordinate_base"}:
            continue
        if not isinstance(value, bool):
            raise ValueError(f"Final manifest option {key} must be a boolean")
    cpg_representation(active_options["cpg_coordinate_base"])
    return FinalManifestOptions(
        backend=backend, cpg_coordinate_base=active_options["cpg_coordinate_base"],
        generate_snp=active_options["generate_snp"], keep_final_bam=active_options["keep_final_bam"],
        publish_rna_bam=active_options["publish_rna_bam"],
    )


def _plan_expected_runs(
    cells: dict[str, dict[str, object]],
    demux_state: Path,
) -> dict[str, dict[str, str]]:


    expected_runs: dict[str, dict[str, str]] = {}
    for record in cells.values():
        project_id = str(record.get("project_sample_id", "") or "").strip()
        protocol = str(record.get("protocol", "") or "").strip().lower()
        dna_raw_sample = str(record.get("dna_raw_sample", "") or "").strip()
        rna_sample = str(record.get("rna_sample", "") or "").strip()
        routes = []
        if dna_raw_sample:
            routes.append((project_id, default_dna_route(protocol), dna_raw_sample))
        if rna_sample and "rna_sample" in get_protocol(protocol).raw_input_columns:
            routes.append((f"{project_id}__rna", "srd_rna_enrichment_tube", rna_sample))
        for run_id, route_id, raw_sample in routes:
            expected = {
                "project_sample_id": project_id,
                "route_id": route_id,
                "raw_sample": raw_sample,
                "protocol": protocol,
            }
            if expected_runs.setdefault(run_id, expected) != expected:
                raise ValueError(f"Demux run_id={run_id!r} has conflicting project sample routes")
    missing_run_files = []
    for run_id in sorted(expected_runs):
        for name in ("manifest.tsv", "completion.json"):
            path = demux_state / run_id / name
            if not path.is_file():
                missing_run_files.append(path)
    if missing_run_files:
        raise FileNotFoundError(
            "Current run is missing committed demux route file(s): "
            + ", ".join(str(path) for path in missing_run_files)
        )
    return expected_runs


def _require_declared_routes(cells: dict[str, dict[str, object]]) -> None:
    for cell_id, record in cells.items():
        protocol = str(record.get("protocol", "") or "").strip().lower()
        if str(record.get("dna_raw_sample", "") or "").strip():
            expected_route = default_dna_route(protocol)
            if record.get("dna_route_id") != expected_route:
                raise ValueError(
                    f"Cell {cell_id!r} declares dna_raw_sample but has no matching committed "
                    f"demux route {expected_route!r}"
                )
        if "rna_sample" in get_protocol(protocol).raw_input_columns and str(record.get("rna_sample", "") or "").strip():
            if record.get("rna_enrichment_route_id") != "srd_rna_enrichment_tube":
                raise ValueError(
                    f"Cell {cell_id!r} declares rna_sample but has no matching committed "
                    "srd_rna_enrichment_tube demux route"
                )


def _merge_committed_generations(
    demux_state: Path,
    initial_cell_manifest: str,
    cells: dict[str, dict[str, object]],
    expected_runs: dict[str, dict[str, str]],
) -> None:
    for run_id, expected_run in sorted(expected_runs.items()):
        manifest = demux_state / run_id / "manifest.tsv"
        completion = demux_state / run_id / "completion.json"

        committed = read_demux_generation(
            completion,
            expected_manifest_path=manifest,
            expected_cell_manifest_path=(demux_state.parent / "calling" / run_id / "cell_manifest.tsv")
            if expected_run["protocol"] == "droplet" else initial_cell_manifest,
        )
        if expected_run["protocol"] == "droplet":
            local_cells = read_cell_identity_manifest(committed["_cell_manifest_path"])
            expected_cells = {cell_id: {key: str(row[key]) for key in INITIAL_CELL_MANIFEST_FIELDS}
                              for cell_id, row in cells.items()
                              if row["project_sample_id"] == expected_run["project_sample_id"]}
            if local_cells != expected_cells:
                raise ValueError(f"Droplet library identity differs from project cell manifest: {run_id!r}")
        if committed.get("run_id") != run_id:
            raise ValueError(
                f"Demux completion {completion} declares run_id={committed.get('run_id')!r}; "
                f"expected {run_id!r}"
            )
        for demux in committed["_rows"]:
            cell_id = str(demux.get("demux_sample", "") or "").strip()
            if not cell_id or cell_id == "__SUMMARY__":
                continue
            if cell_id not in cells:
                raise ValueError(
                    f"Demux manifest {manifest} contains cell {cell_id!r} absent from initial cell manifest"
                )
            record = cells[cell_id]
            project_id = str(demux.get("project_sample_id", "") or "").strip()


            route = canonical_raw_route_id(str(demux.get("route_id", "") or ""))
            raw = str(demux.get("raw_sample", "") or "")
            if (
                project_id != expected_run["project_sample_id"]
                or route != expected_run["route_id"]
                or raw != expected_run["raw_sample"]
            ):
                raise ValueError(
                    f"Demux route {run_id!r} does not match the declared sample route: "
                    f"observed=({project_id!r}, {route!r}, {raw!r}), "
                    f"expected=({expected_run['project_sample_id']!r}, "
                    f"{expected_run['route_id']!r}, {expected_run['raw_sample']!r})"
                )
            if (
                str(demux.get("pipeline_mode", "")).strip().lower() == "srd"
                and str(demux.get("rna_status", "")) == "Pass"
            ):
                intermediate = demux_state.parent.parent
                project_root = intermediate.parent
                collected = [
                    rna_fastq_collect_destination(intermediate, route, cell_id, mate)
                    for mate in (1, 2)
                ]
                missing = [path for path in collected if not path.is_file()]
                if missing:
                    raise FileNotFoundError(
                        "Committed SRD RNA FASTQ is missing from the collection area "
                        "(rna_collect_fastq must run first): "
                        + ", ".join(str(path) for path in missing)
                    )
                r1, r2 = (str(path.relative_to(project_root)) for path in collected)
                if route == "srd_rna_enrichment_tube":
                    record["rna_enrichment_rna_fastq_r1"] = r1
                    record["rna_enrichment_rna_fastq_r2"] = r2
                else:
                    record["dna_tube_rna_fastq_r1"] = r1
                    record["dna_tube_rna_fastq_r2"] = r2
            if route == "srd_rna_enrichment_tube":
                if record.get("rna_enrichment_route_id"):
                    raise ValueError(f"Cell {cell_id!r} has more than one SRD RNA-enrichment demux route")
                record.update({
                    "rna_enrichment_route_id": route,
                    "rna_enrichment_raw_sample": raw,
                    "rna_enrichment_dna_reads": demux.get("dna_reads", "0"),
                    "rna_enrichment_dna_status": demux.get("dna_status", ""),
                    "rna_enrichment_rna_reads": demux.get("rna_reads", "0"),
                    "rna_enrichment_rna_status": demux.get("rna_status", ""),
                    "rna_enrichment_input_read_pairs": demux.get("input_read_pairs", "0"),
                    "rna_enrichment_matched_read_pairs": demux.get("matched_reads", "0"),
                    "rna_enrichment_unassigned_rate": demux.get("unassigned_rate", ""),
                })
            elif str(demux.get("downstream_dna", "")).strip().lower() == "true":
                if record.get("dna_route_id"):
                    raise ValueError(f"Cell {cell_id!r} has more than one downstream-DNA demux route")
                record.update({
                    "dna_route_id": route,
                    "dna_raw_sample": raw,
                    "dna_reads": demux.get("dna_reads", "0"),
                    "dna_status": demux.get("dna_status", ""),
                    "dna_tube_rna_reads": demux.get("rna_reads", "0"),
                    "dna_tube_rna_status": demux.get("rna_status", ""),
                    "dna_input_read_pairs": demux.get("input_read_pairs", "0"),
                    "dna_matched_read_pairs": demux.get("matched_reads", "0"),
                    "dna_unassigned_rate": demux.get("unassigned_rate", ""),
                    "high_cph_role": demux.get("high_cph_role", ""),
                })


def _taps_cell_metrics(
    cell_id: str, paths: Mapping[str, Path], cutadapt: Mapping[str, object],
) -> NativeCellMetrics:
    qc = json.loads(paths["rastair_qc"].read_text())
    if qc["sample"] != cell_id or qc["backend"] != "rastair":
        raise ValueError("TAPS QC identity disagrees with the active cell")
    for key in ("trimmed_pairs", "bam_primary_records", "backend_accepted_pairs", "backend_rejected_pairs",
                "duplicate_pairs", "postdedup_pairs", "unmapped_pairs", "discordant_pairs", "quality_rejected_pairs"):
        if type(qc[key]) is not int or qc[key] < 0:
            raise ValueError(f"TAPS QC requires a nonnegative integer: {key}")
    if qc["trimmed_pairs"] != cutadapt["trimmed_pairs"] or qc["bam_primary_records"] != 2 * qc["trimmed_pairs"]:
        raise ValueError("TAPS BAM and Cutadapt pair counts disagree")
    accepted, rejected, duplicates, postdedup = (qc[k] for k in (
        "backend_accepted_pairs", "backend_rejected_pairs", "duplicate_pairs", "postdedup_pairs"))
    if accepted + rejected != qc["trimmed_pairs"] or postdedup + duplicates != accepted:
        raise ValueError("TAPS pair funnel does not close")
    if sum(qc[k] for k in ("unmapped_pairs", "discordant_pairs", "quality_rejected_pairs")) != rejected:
        raise ValueError("TAPS rejection categories do not close")
    non_cpg = {k: v for k, v in qc.items() if k.startswith("non_cpg_")}
    if non_cpg.get("non_cpg_metric_source") != "taps_oriented_host_chg_chh_conversion":
        raise ValueError("TAPS non-CpG QC source is invalid")
    for key in ("non_cpg_chg_mod", "non_cpg_chg_unmod", "non_cpg_chh_mod", "non_cpg_chh_unmod",
                "non_cpg_observations", "non_cpg_candidate_pairs", "non_cpg_oriented_pairs",
                "non_cpg_ambiguous_pairs", "snp_sites"):
        if type(qc[key]) is not int or qc[key] < 0:
            raise ValueError(f"TAPS QC requires a nonnegative integer: {key}")
    modified = sum(qc[f"non_cpg_{context}_mod"] for context in ("chg", "chh"))
    total = modified + sum(qc[f"non_cpg_{context}_unmod"] for context in ("chg", "chh"))
    if (total != qc["non_cpg_observations"]
            or qc["non_cpg_oriented_pairs"] + qc["non_cpg_ambiguous_pairs"] != qc["non_cpg_candidate_pairs"]
            or not 0 <= qc["non_cpg_candidate_pairs"] <= postdedup):
        raise ValueError("TAPS non-CpG accounting does not close")
    expected_pct = 100.0 * modified / total if total else float("nan")
    if not (math.isnan(expected_pct) and math.isnan(qc["non_cpg_methylation_pct"])) and not math.isclose(expected_pct, qc["non_cpg_methylation_pct"]):
        raise ValueError("TAPS non-CpG percentage disagrees with counts")
    return NativeCellMetrics(accepted, rejected, duplicates, postdedup, {
        **{k: v for k, v in qc.items() if k.startswith(("bam_", "lambda_", "pUC19_"))},
        "snp_sites": qc["snp_sites"], "methylation_signal": "5mC+5hmC",
        "mapping_policy": "bwa_rastair_eligible_pairs", "native_mapping_policy": "bwa_proper_primary_any_mate_mapq20",
        "native_mapping_unit": "read_pairs", "native_mapping_pct": 100.0 * accepted / qc["trimmed_pairs"] if qc["trimmed_pairs"] else 0.0,
        "both_primary_mapped_pairs": accepted, "one_primary_mapped_pairs": 0,
        "unmapped_pairs": qc["unmapped_pairs"], "ambiguous_pairs": "", "no_primary_pairs": "",
        "discordant_pairs": qc["discordant_pairs"], "quality_rejected_pairs": qc["quality_rejected_pairs"],
        "final_retained_pairs": postdedup,
        "dedup_policy": "samtools_markdup_flag_exclude_in_rastair", "pair_qc_metric_unit": "read_pairs",
        "high_cph_assessed_pairs": 0, "high_cph_flagged_pairs": 0, "high_cph_removed_pairs": 0,
        "high_cph_fraction": "nan", "high_cph_role": "not_applicable", **non_cpg,
    })


def _biscuit_cell_metrics(
    cell_id: str, paths: Mapping[str, Path], trimmed_pairs: int,
) -> NativeCellMetrics:
    if not paths["dupsifter_multiqc"].is_file():
        raise FileNotFoundError(f"Active cell is missing Dupsifter QC: {paths['dupsifter_multiqc']}")
    _native, pairs = _read_dupsifter_metrics(paths["dupsifter_multiqc"], cell_id, trimmed_pairs)
    counts = [pairs.pop(key) for key in (
        "backend_accepted_pairs", "backend_rejected_pairs", "duplicate_pairs", "postdedup_pairs")]
    return NativeCellMetrics(*counts, {
        **_read_biscuit_native_mapping(paths["mapq_table"]), **pairs, "ambiguous_pairs": "",
        "mapping_policy": "biscuit_any_primary_pair_funnel",
        "dedup_policy": "dupsifter_wgbs_signature_remove_dups",
        "non_cpg_methylation_pct": _parse_cph_methylation_pct(paths["cph_table"]),
        "non_cpg_metric_source": "biscuit_cph_retention_by_read_position",
    })


def _bismark_cell_metrics(
    job: CellMetricJob, paths: Mapping[str, Path], trimmed_pairs: int,
    summary: Mapping[str, object],
) -> NativeCellMetrics:
    metrics = build_metrics(
        sample=job.cell_id, version=job.bismark_version,
        alignment_report=paths["bismark_align_report"], dedup_report=paths["bismark_dedup_report"],
        filter_summary=summary, protocol=job.protocol,
    )
    total_pairs = int(metrics["total_pairs"])
    if total_pairs != trimmed_pairs:
        raise ValueError(
            f"Bismark report disagrees with Cutadapt output for {job.cell_id!r}: "
            f"bismark={total_pairs}, cutadapt={trimmed_pairs}"
        )
    accepted = int(metrics["accepted_pairs"])
    rejected = int(metrics["unmapped_pairs"]) + int(metrics["ambiguous_pairs"])
    if accepted + rejected != trimmed_pairs:
        raise ValueError(f"Bismark pair funnel is inconsistent for {job.cell_id!r}")
    fields = {
        "native_mapping_policy": "bismark_unique_concordant_pairs",
        "native_mapping_unit": "read_pairs",
        "native_mapping_pct": 100.0 * accepted / trimmed_pairs if trimmed_pairs else 0.0,
        "both_primary_mapped_pairs": accepted, "one_primary_mapped_pairs": 0,
        "unmapped_pairs": int(metrics["unmapped_pairs"]), "ambiguous_pairs": int(metrics["ambiguous_pairs"]),
        "no_primary_pairs": "", "mapping_policy": "bismark_unique_concordant",
        "dedup_policy": DROPLET_DEDUP_POLICY if job.protocol == "droplet" else "bismark_paired_endpoint_orientation",
    }
    if job.protocol == "droplet":
        contexts = metrics["umi_qc"]["contexts"]
        cph_m = sum(contexts[k]["methyl"] for k in ("CHG", "CHH"))
        cph_total = cph_m + sum(contexts[k]["unmethyl"] for k in ("CHG", "CHH"))
        fields.update({
            "non_cpg_methylation_pct": 100.0 * cph_m / cph_total if cph_total else float("nan"),
            "non_cpg_metric_source": "droplet_umi_consensus_chg_chh",
            "molecule_qc_metric_unit": "cytosine_umi_consensus",
            "cpg_read_observations": contexts["CpG"]["read_observations"],
            "cpg_umi_observations": contexts["CpG"]["methyl"] + contexts["CpG"]["unmethyl"],
            "cpg_ambiguous_umi_sites": contexts["CpG"]["ambiguous_umis"],
        })
    else:
        fields.update({
            "non_cpg_methylation_pct": _parse_bismark_non_cpg_methylation_pct(paths["bismark_extract_report"]),
            "non_cpg_metric_source": "bismark_extraction_chg_chh",
        })
    return NativeCellMetrics(accepted, rejected, int(metrics["duplicate_pairs"]),
                             int(metrics["postdedup_pairs"]), fields)


def _filtered_cell_metrics(
    job: CellMetricJob, paths: Mapping[str, Path], trimmed_pairs: int,
) -> NativeCellMetrics:
    if not paths["high_cph_summary"].is_file():
        raise FileNotFoundError(f"Active cell is missing high-CpH summary: {paths['high_cph_summary']}")
    summary = _read_high_cph_summary(paths["high_cph_summary"], job.cell_id, job.options.backend)
    if summary["protocol"] != job.protocol:
        raise ValueError(f"High-CpH summary protocol disagrees with the cell identity: {paths['high_cph_summary']}")
    if not summary["filter_applied"]:
        raise ValueError(f"High-CpH summary must record the fixed filter action: {paths['high_cph_summary']}")
    native = (_biscuit_cell_metrics(job.cell_id, paths, trimmed_pairs) if job.options.backend == "biscuit"
              else _bismark_cell_metrics(job, paths, trimmed_pairs, summary))
    assessed = int(summary["candidate_read_pairs"])
    removed = int(summary["flagged_read_pairs"])
    if assessed + int(summary["excluded_read_pairs"]) > native.postdedup_pairs:
        raise ValueError(f"High-CpH accounting exceeds post-dedup pairs for {job.cell_id!r}")
    native.fields.update({
        "high_cph_assessed_pairs": assessed, "high_cph_flagged_pairs": removed,
        "high_cph_removed_pairs": removed, "final_retained_pairs": native.postdedup_pairs - removed,
        "pair_qc_metric_unit": "read_pairs", "high_cph_fraction": float(summary["high_cph_fraction"]),
    })
    return native


def _metrics_for_cell(job: CellMetricJob) -> tuple[str, dict[str, object]]:
    options = job.options
    artifacts = cell_artifacts(job.paths, job.protocol,
        generate_snp=options.generate_snp or options.backend == "rastair", keep_final_bam=options.keep_final_bam)
    cutadapt = _parse_cutadapt_metrics(artifacts.metrics["cutadapt_json"])
    try:
        demux_pairs = int(job.dna_reads.strip())
    except ValueError as exc:
        raise ValueError(f"Active cell {job.cell_id!r} has invalid dna_reads") from exc
    if cutadapt["cutadapt_input_pairs"] != demux_pairs:
        raise ValueError(
            f"Cutadapt input disagrees with demux DNA pairs for {job.cell_id!r}: "
            f"cutadapt={cutadapt['cutadapt_input_pairs']}, demux={demux_pairs}"
        )
    if options.backend == "rastair":
        native = _taps_cell_metrics(job.cell_id, artifacts.metrics, cutadapt)
    else:
        native = _filtered_cell_metrics(job, artifacts.metrics, int(cutadapt["trimmed_pairs"]))
    parsed_inputs = set(artifacts.metrics.values())
    for path in artifacts.final_inputs:
        if path in parsed_inputs:
            continue
        if not path.is_file():
            raise FileNotFoundError(f"Active cell is missing a declared result or QC file: {path}")
    project = job.public_paths.result_root.parent
    updates = {**cutadapt, **native.manifest_fields(),
               "methylation_backend": options.backend, "backend_qc_status": "PASS", "analysis_status": "PASS",
               "cpg_path": str(job.public_paths.cpg.relative_to(project)),
               "cpg_representation": cpg_representation(options.cpg_coordinate_base),
               "snp_path": str(job.public_paths.snp.relative_to(project)) if options.generate_snp or options.backend == "rastair" else "",
               "rna_bam_path": str(job.paths.rna_bam.relative_to(project)) if options.publish_rna_bam else ""}
    if options.backend == "rastair":
        updates.update({
            "retained_bam_path": str(job.paths.marked_bam.relative_to(project)) if options.keep_final_bam else "",
            "retained_bam_index": str(job.paths.marked_bam_index.relative_to(project)) if options.keep_final_bam else "",
            "retained_bam_role": "taps_all_alignments_marked_duplicates_for_variant_reanalysis" if options.keep_final_bam else "",
        })
    return job.cell_id, updates


def _apply_cell_metrics(
    cells: dict[str, dict[str, object]],
    *,
    layout: ProjectLayout,
    published_layout: ProjectLayout,
    options: FinalManifestOptions,
    bismark: str | None,
    workers: int = 1,
) -> None:
    backend = options.backend
    bismark_version = _bismark_version_text(bismark) if backend == "bismark" else ""
    jobs: list[CellMetricJob] = []
    for cell_id, record in cells.items():
        is_active = str(record.get("dna_status", "")) == "Pass"
        if not is_active:
            record.update(
                {
                    "methylation_backend": backend,
                    "backend_qc_status": "NOT_RUN",
                    "cpg_path": "",
                    "cpg_representation": "",
                    "snp_path": "",
                    "rna_bam_path": "",
                    "analysis_status": "NOT_ANALYZED",
                }
            )
            continue
        jobs.append(CellMetricJob(
            cell_id=cell_id, protocol=str(record["protocol"]).strip().lower(),
            dna_reads=str(record.get("dna_reads", "")), paths=layout.cell(cell_id),
            public_paths=published_layout.cell(cell_id), options=options, bismark_version=bismark_version,
        ))
    if not jobs:
        return
    workers = min(effective_cpu_count(workers), 8, len(jobs))
    if workers == 1:
        for job in jobs:
            cell_id, updates = _metrics_for_cell(job)
            cells[cell_id].update(updates)
        return
    try:
        import multiprocessing
        from concurrent.futures import ProcessPoolExecutor
        from concurrent.futures.process import BrokenProcessPool

        pool_context = multiprocessing.get_context("fork")
        pool = ProcessPoolExecutor(max_workers=workers, mp_context=pool_context)
    except (OSError, ValueError) as exc:
        print(f"WARNING: metric worker pool unavailable ({exc}); falling back to sequential",
              file=sys.stderr)
        pending = [_metrics_for_cell(job) for job in jobs]
    else:
        try:
            with pool:
                pending = list(pool.map(_metrics_for_cell, jobs, chunksize=max(1, min(16, len(jobs) // workers))))
        except BrokenProcessPool as exc:
            print(f"WARNING: metric worker pool unavailable ({exc}); falling back to sequential",
                  file=sys.stderr)
            pending = [_metrics_for_cell(job) for job in jobs]
    for cell_id, updates in pending:
        cells[cell_id].update(updates)


def _emit_final_tsv(
    cells: dict[str, dict[str, object]],
    output: Path,
    pair_qc_multiqc_json: str | None,
    cutadapt_qc_multiqc_json: str | None,
) -> None:
    preferred = [
        "sample_id", "project_sample_id", "species", "protocol", "plate_id", "cell_order",
        "dna_barcode", "rna_barcode", "dna_raw_sample", "rna_sample",
        "dna_route_id", "dna_reads", "dna_status",
        "dna_tube_rna_reads", "dna_tube_rna_status", "dna_input_read_pairs",
        "dna_matched_read_pairs", "dna_unassigned_rate",
        "dna_tube_rna_fastq_r1", "dna_tube_rna_fastq_r2",
        "rna_enrichment_route_id", "rna_enrichment_raw_sample",
        "rna_enrichment_dna_reads", "rna_enrichment_dna_status",
        "rna_enrichment_rna_reads", "rna_enrichment_rna_status",
        "rna_enrichment_input_read_pairs", "rna_enrichment_matched_read_pairs",
        "rna_enrichment_unassigned_rate", "rna_enrichment_rna_fastq_r1",
        "rna_enrichment_rna_fastq_r2", "analysis_status",
        "cutadapt_input_pairs", "trimmed_pairs", "both_primary_mapped_pairs",
        "one_primary_mapped_pairs", "backend_accepted_pairs", "backend_rejected_pairs",
        "unmapped_pairs", "ambiguous_pairs", "no_primary_pairs", "duplicate_pairs",
        "postdedup_pairs", "high_cph_assessed_pairs", "high_cph_flagged_pairs",
        "high_cph_removed_pairs", "final_retained_pairs",
        "mapping_policy", "native_mapping_policy", "native_mapping_unit",
        "native_mapping_pct", "dedup_policy", "pair_qc_metric_unit",
        "r1_adapter_pct", "r2_adapter_pct", "non_cpg_methylation_pct",
        "non_cpg_metric_source", "high_cph_fraction",
        "high_cph_role", "methylation_backend", "backend_qc_status",
        "cpg_path", "cpg_representation", "snp_path",
        "rna_bam_path", "notes",
    ]
    rows = list(cells.values())
    extra = sorted({key for row in rows for key in row} - set(preferred))
    fields = preferred + extra

    def cell_sort_key(item: Mapping[str, object]) -> tuple[object, ...]:
        raw_order = str(item.get("cell_order", "") or "").strip()
        return (
            str(item.get("project_sample_id", "")),
            _cell_order_key(raw_order),
            str(item.get("plate_id", "")),
            str(item.get("sample_id", "")),
        )

    sorted_rows = sorted(rows, key=cell_sort_key)
    if pair_qc_multiqc_json is not None:
        _write_pair_level_multiqc(sorted_rows, Path(pair_qc_multiqc_json))
    if cutadapt_qc_multiqc_json is not None:
        _write_cutadapt_multiqc(sorted_rows, Path(cutadapt_qc_multiqc_json))

    write_tsv_atomic(output, sorted_rows, fields)


def should_skip_biscuit_dup_report(report_path: str) -> str | None:
    """返回 BISCUIT duplicate 报告不适合进入 MultiQC 的原因；安全时返回 None。"""
    path = Path(report_path)
    try:
        text = path.read_text(errors="replace")
    except OSError as exc:
        return f"unreadable file: {exc}"

    required_patterns = {
        "total_reads": r"Number of reads:\s+(\d+)",
        "q40_reads": r"Number of q40-reads:\s+(\d+)",
    }
    values = {}
    for key, pattern in required_patterns.items():
        match = re.search(pattern, text, re.MULTILINE)
        if match is None:
            return f"missing required field for MultiQC biscuit parser: {key}"
        values[key] = int(match.group(1))

    reasons = []
    if values["total_reads"] == 0:
        reasons.append("Number of reads == 0")
    if values["q40_reads"] == 0:
        reasons.append("Number of q40-reads == 0")
    if reasons:
        return ", ".join(reasons)
    return None


def _is_multiqc_source_candidate(path: Path) -> bool:
    name = path.name

    if not path.is_file():
        return False


    large_suffixes = (
        ".fastq.gz", ".fq.gz", ".bam", ".bai", ".cram", ".crai",
        ".vcf.gz", ".vcf.gz.tbi", ".bed.gz", ".bed.gz.tbi",
        ".tsv.zst", ".zst",
    )
    if name.endswith(large_suffixes):
        return False


    if name.endswith("_mqc.json"):
        return True
    if name.endswith(".cutadapt.json"):
        return True
    if name.endswith("_dup_report.txt"):
        return True


    if name.endswith(
        (
            "_PE_report.txt",
            ".deduplication_report.txt",
            "_splitting_report.txt",
            "M-bias.txt",
        )
    ):
        return True


    parts = path.parts
    # 细胞 ID 可含点号；不能把 ID 内部的点号当作报告的复合格式后缀而静默丢弃。
    if tuple(parts[-4:-2]) == ("biscuit", "qc") and path.suffix in {
        ".txt", ".json", ".tsv", ".csv"
    }:
        return True

    return False


def _iter_explicit_multiqc_files(inputs: Iterable[str | Path]) -> Iterable[Path]:
    for raw in inputs:
        source = Path(raw).expanduser()
        if not source.exists():
            raise FileNotFoundError(f"Declared MultiQC source is missing: {source}")
        if not source.is_file():
            raise ValueError(f"Declared MultiQC source is not a file: {source}")
        yield source


def prefetch_multiqc_native_libraries(prefix: Path) -> tuple[int, int]:
    """有界预读当前 MultiQC 的解释器与 Polars 映像，缓解 Lustre 冷页映射等待。

    不改环境或库的选择；预读失败仅告警，原生 MultiQC 的退出码仍决定规则成败。
    """
    seen: set[tuple[int, int]] = set()
    files, read_bytes = 0, 0
    for pattern in (
        "lib/python*/site-packages/_polars_runtime_*/*.so",
        "lib/python*/site-packages/polars/*.so",
        "lib/libpython*.so*",
        "lib/libpython*.dylib",
        # 解释器也可能在加载 Polars 前阻塞；最后读取，使其在 exec 前保持热页。
        "bin/python",
    ):
        for library in sorted(prefix.glob(pattern)):
            try:
                identity = library.stat()
                if not stat.S_ISREG(identity.st_mode):
                    continue
                key = (identity.st_dev, identity.st_ino)
                if key in seen:
                    continue
                seen.add(key)
                with library.open("rb") as handle:
                    while block := handle.read(8 * 1024 * 1024):
                        read_bytes += len(block)
                files += 1
            except OSError as exc:
                logger.warning("MultiQC 原生库预读失败，继续由原生工具校验环境：%s: %s", library, exc)
    return files, read_bytes


def run_multiqc_report(file_list: Path, config_file: Path, output_dir: Path, comment: str) -> int:
    """在独立 backend 进程内，仅以原生默认值及 Pipeline 显式配置生成报告。"""
    import multiqc
    from multiqc.core.update_config import ClConfig

    sources = json.loads(file_list.read_text(encoding="utf-8"))
    if (not isinstance(sources, list) or not sources
            or any(not isinstance(path, str) or not path for path in sources)):
        raise ValueError("MultiQC source list must be a nonempty JSON array of paths")

    # 原生 CLI 无禁用隐式配置的选项；其递归合并还会保留外部改名、过滤和
    # custom_data 设置。只隔离加载入口，不改 HOME、包文件或原生解析/退出行为。
    find_user_files = multiqc.config.find_user_files
    multiqc.config.find_user_files = lambda: None
    try:
        result = multiqc.run(
            # 原生逐行清单无法表示含换行的路径；positional API 又会做glob。
            # JSON保留路径字符，逐个转义后交原生搜索，防止误纳入同名诱饵。
            *(glob.escape(path) for path in sources),
            cfg=ClConfig(
                config_files=[config_file.resolve(strict=True)],
                output_dir=output_dir, filename="multiqc_report.html",
                force=True, report_comment=comment,
            ),
            interactive=False,
        )
        return result.sys_exit_code
    finally:
        multiqc.config.find_user_files = find_user_files


def write_multiqc_file_list(
    inputs: Iterable[str | Path],
    out_list: str,
) -> tuple[int, list[str]]:
    """写出报告进程消费的唯一 JSON 源列表，保留路径字符及完整 DAG 依赖。"""
    unique_sources: dict[str, Path] = {}
    skipped: list[str] = []
    for candidate in _iter_explicit_multiqc_files(inputs):
        if not _is_multiqc_source_candidate(candidate):
            continue
        source = candidate.resolve(strict=True)
        if source.name.endswith("_dup_report.txt"):
            reason = should_skip_biscuit_dup_report(str(source))
            if reason is not None:
                skipped.append(f"{source}\t{reason}")
                continue
        unique_sources[str(source)] = source

    if not unique_sources:
        raise ValueError(
            "The explicit MultiQC source set is empty; refusing to publish an empty report"
        )

    file_list = Path(out_list)
    atomic_write_bytes(file_list, canonical_json_bytes(sorted(unique_sources)) + b"\n")

    logger.info(
        "Prepared explicit MultiQC file list with %d sources; skipped %d dup reports",
        len(unique_sources),
        len(skipped),
    )
    return len(unique_sources), skipped


def select_multiqc_samples(samples: Iterable[str]) -> list[str]:
    """按细胞 ID 的 SHA256 固定选取最多 1000 个细胞，不依赖输入顺序或 Python 随机种子。"""
    ranked = sorted(set(samples), key=lambda sample: (hashlib.sha256(sample.encode("utf-8")).digest(), sample))
    return sorted(ranked[:MULTIQC_MAX_CELLS])


def prepare_sampled_multiqc_sources(
    inputs: Iterable[str | Path], selected_samples: Sequence[str], total_cells: int,
    out_dir: Path,
) -> tuple[list[str], str]:
    """统一展示子集及 Bismark 报告后缀；原生源须已按集合选定，不修改任何源文件。"""
    selected = set(selected_samples)
    sampled = total_cells > len(selected)
    out_dir.mkdir(parents=True, exist_ok=True)
    atomic_write_json(out_dir / "multiqc_config.json", {
        "fn_clean_trim": [],
        # Cutadapt 对路径内任意连字符的 stdin 误判会退回共用 scratch basename。
        # 显式使用本 Pipeline 的细胞报告名，且只去掉各模块的末尾工件后缀；
        # MultiQC 默认的中间子串截断也会破坏合法细胞身份和 custom table 的关联。
        # BISCUIT 会重复清名；每次从原报告名开始，避免再次截断 ID 自身的同名结尾。
        "use_filename_as_sample_name": ["cutadapt", "biscuit"],
        "fn_clean_exts": [{
            "type": "regex", "module": "bismark",
            "pattern": r"(?:_PE_report|_splitting_report|\.deduplication_report|\.M-bias)\.txt$",
        }, {
            "type": "regex", "module": "cutadapt", "pattern": r"\.cutadapt\.json$",
        }, {
            "type": "regex", "module": "biscuit",
            "pattern": r"_(?:" + "|".join(map(re.escape, MULTIQC_BISCUIT_SOURCE_SUFFIXES)) + r")$",
        }],
    })
    sources = []
    for raw in inputs:
        path = Path(raw)
        if sampled and path.name.endswith("_mqc.json"):
            payload = json.loads(path.read_text(encoding="utf-8"))
            data = payload["data"]
            if set(data) - selected:
                payload["data"] = {sample: values for sample, values in data.items() if sample in selected}
                payload["description"] = (
                    str(payload.get("description", ""))
                    + " 细胞表仅展示本报告抽样集合；文库级 Run summary 仍为全库计数。"
                )
                projected = out_dir / (hashlib.sha256(str(path).encode("utf-8")).hexdigest() + "_mqc.json")
                atomic_write_json(projected, payload)
                path = projected
        sources.append(str(path))
    scope = (
        f"MultiQC 逐细胞质控范围：{len(selected):,} / {total_cells:,} 个活跃 DNA 细胞。"
        + ("超过 1,000 个细胞时，按细胞 ID 的 SHA256 固定抽样 1,000 个；各逐细胞模块使用同一集合。"
           if sampled else "未抽样，保留全部细胞的质控输入。")
        + "文库级解复用概要使用全库计数。完整科学结果和最终 sample manifest 保留全部细胞。"
    )
    return sources, scope


def prepare_high_cph_multiqc_sources(
    inputs: Iterable[str | Path], out_json: Path,
) -> list[str]:
    """从保留的过滤摘要重建已清理的展示 JSON，避免报告重建回溯科学规则。"""
    sources: list[str] = []
    aggregate: dict[str, object] = {
        "id": "dna_pipeline_high_cph_rebuilt", "section_name": "High-CpH Read Assessment",
        "plot_type": "table", "data": {},
    }
    for raw in inputs:
        path = Path(raw)
        if not path.name.endswith("_filter_summary.json"):
            sources.append(str(path))
            continue
        summary = json.loads(path.read_text(encoding="utf-8"))
        payload = _high_cph_multiqc_payload(
            summary["sample"], summary, summary["protocol"], summary["high_cph_role"]
        )
        if not aggregate["data"]:
            aggregate = payload
        else:
            aggregate["data"].update(payload["data"])
    atomic_write_json(out_json, aggregate)
    sources.append(str(out_json))
    return sources


# 汇总与交付：运行快照和持锁发布

SNAPSHOT_SCHEMA_VERSION = 11


_SOURCE_ROOTS = ("core",)


_SOURCE_SUFFIXES = {".py", ".sh", ".yaml", ".toml", ".rs", ".lock"}


_SOURCE_EXCLUDED_PARTS = {
    ".git", ".ipynb_checkpoints", ".pytest_cache", ".snakemake", "__pycache__", "conda", "target", "vendor"
}


def _pipeline_source_identity(root: Path) -> Mapping[str, object]:
    records: list[dict[str, object]] = []
    for top in _SOURCE_ROOTS:
        base = root / top
        if not base.exists():
            continue
        for path in sorted(base.rglob("*")):
            if not path.is_file() or any(part in _SOURCE_EXCLUDED_PARTS for part in path.relative_to(base).parts):
                continue
            if path.suffix.lower() not in _SOURCE_SUFFIXES and path.name != "Snakefile" and path.parent != root / "core/reference":
                continue
            identity = stable_file_identity(path)
            records.append({
                "path": str(path.relative_to(root)),
                "size_bytes": identity["size_bytes"],
                "sha256": identity["sha256"],
            })
    design = root / DROPLET_DESIGN_RESOURCE
    if design.is_file():
        identity = stable_file_identity(design)
        records.append({"path": DROPLET_DESIGN_RESOURCE, "size_bytes": identity["size_bytes"],
                        "sha256": identity["sha256"]})
    return {
        "method": "production-source-tree-v1",
        "sha256": sha256_bytes(canonical_json_bytes(records)),
        "file_count": len(records),
    }


def _environment_observation(label: str, prefix: Path) -> Mapping[str, object]:
    resolved = prefix.expanduser().resolve()
    warnings: list[str] = []
    try:
        record: Mapping[str, object] | None = read_environment_identity(resolved)
    except (FileNotFoundError, ValueError, OSError) as exc:
        record = None
        warnings.append(str(exc))
    return {
        "label": label,
        "resolved_prefix": str(resolved),
        "record": record,
        "capture_warnings": warnings,
    }


def _fastq_inventory(
    raw_dir: Path, manifest_rows: list[ProjectSampleRow], *, protocol: str,
) -> list[Mapping[str, object]]:
    discovered = discover_samples(raw_dir, selected_samples=manifest_raw_samples(manifest_rows, protocol))
    map_discovered_samples_to_manifest(discovered, manifest_rows, protocol=protocol)
    if not discovered:
        raise ValueError(f"no paired FASTQ samples found under: {raw_dir}")
    result: list[Mapping[str, object]] = []
    for sample, mates in sorted(discovered.items()):
        r1s, r2s = list(mates["r1"]), list(mates["r2"])
        if len(r1s) != len(r2s):
            raise ValueError(f"FASTQ inventory became unpaired for sample {sample!r}")
        for segment, (r1, r2) in enumerate(zip(r1s, r2s, strict=True), start=1):
            result.append({
                "raw_sample": sample,
                "segment": segment,
                "r1": stable_file_identity(Path(r1), hash_content=False),
                "r2": stable_file_identity(Path(r2), hash_content=False),
            })
    return result


def _executable_identity(path: Path) -> Mapping[str, object]:
    if not path.is_file():
        return {"path": str(path), "present": False}
    observation: dict[str, object] = {
        "present": True,
        **stable_file_identity(path, hash_content=True),
    }
    if path.name == "demux_rs":
        try:
            completed = subprocess.run(
                [str(path), "--build-info"], check=False, capture_output=True,
                text=True, timeout=20,
            )
            observation["build_info_exit_code"] = completed.returncode
            if completed.returncode == 0:
                observation["build_info"] = json.loads(completed.stdout)
            else:
                observation["build_info_error"] = completed.stderr.strip()
        except (OSError, ValueError, subprocess.SubprocessError) as exc:
            observation["build_info_error"] = str(exc)
    return observation


def _backend_reference_observation(
    config: Mapping[str, object], reference: Path
) -> Mapping[str, object]:
    backend = get_backend(dict(config))
    kind = ("bismark-local" if resolve_bismark_alignment(config).local_alignment else "bismark-combined") if backend == "bismark" else backend
    validate_reference_component(reference, kind)

    if backend == "bismark":
        public = Path(f"{reference}.bismark")
        if not public.is_dir():
            raise FileNotFoundError(f"Bismark reference directory is missing or invalid: {public}")
        resolved = public.resolve(strict=True)
        genome = resolved / "genome.fa"
        if not genome.is_file() or genome.stat().st_size == 0:
            raise FileNotFoundError(f"Incomplete Bismark reference generation: {resolved}")
        local_alignment = resolve_bismark_alignment(config).local_alignment
        indexes = backend_active_index_paths(
            reference, "bismark", bismark_local_alignment=local_alignment
        )
        missing = [path for path in indexes if not path.is_file() or path.stat().st_size == 0]
        if missing:
            raise FileNotFoundError(
                "Incomplete Bismark index generation: "
                + ", ".join(str(path) for path in missing[:8])
            )
        return {
            "backend": "bismark",
            "alignment_mode": "local-faithful" if local_alignment else "end-to-end-combined",
            "public_path": str(public.absolute()),
            "resolved_path": str(resolved),
            "genome": stable_file_identity(genome, hash_content=False),
            "indexes": [stable_file_identity(path, hash_content=False) for path in indexes],
        }

    indexes = backend_active_index_paths(reference, backend)
    missing = [path for path in indexes if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError(
            f"Incomplete {backend} reference generation: "
            + ", ".join(str(path) for path in missing)
        )
    return {
        "backend": backend,
        "indexes": [stable_file_identity(path, hash_content=False) for path in indexes],
    }


def _observe_snapshot_input(path: Path) -> tuple[bytes, dict[str, object]]:
    """一次读取小型运行输入；科学解析和内容身份必须使用同一份字节。"""
    before = path.stat()
    content = path.read_bytes()
    after = path.stat()
    if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (
        after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns,
    ):
        raise RuntimeError(f"file changed while reading run input: {path}")
    return content, {
        "path": str(path.resolve()), "size_bytes": len(content),
        "mtime_ns": before.st_mtime_ns, "sha256": sha256_bytes(content),
    }


def _snapshot_basis(
    *,
    pipeline_root: Path,
    project_dir: Path,
    config_path: Path,
    control_env: Path,
    biscuit_env: Path,
    bismark_env: Path,
    rastair_env: Path,
    allow_external_config: bool = False,
) -> Mapping[str, object]:
    config_bytes, config_identity = _observe_snapshot_input(config_path)
    config = _parse_project_config(config_bytes.decode("utf-8"), path=config_path)
    expected_project = config_path.parent.parent.resolve()
    if not allow_external_config and project_dir.resolve() != expected_project:
        raise ValueError(f"project/config mismatch: {project_dir} != {expected_project}")

    sample_manifest = project_dir / "00_config" / "sample_manifest.tsv"
    barcode_map = project_dir / "00_config" / "Barcode_Map.csv"
    for path in (config_path, sample_manifest, *([] if config["analysis"]["protocol"] == "droplet" else [barcode_map])):
        if not path.is_file():
            raise FileNotFoundError(f"run-start input is missing: {path}")

    species = str(config.get("species", "")).strip()
    reference = Path(resolve_reference(config, species, project_dir, pipeline_root=pipeline_root))
    reference_record = reference_identity(
        reference,
        expected_species=species or None,
    )
    control_observation = _environment_observation("control", control_env)
    backend = get_backend(config)
    active_backend_observation = _environment_observation(
        backend, {"biscuit": biscuit_env, "bismark": bismark_env,
                       "rastair": rastair_env}[backend]
    )
    protocol = str(config.get("analysis", {}).get("protocol", "")).strip().lower()
    manifest_bytes, manifest_identity = _observe_snapshot_input(sample_manifest)
    manifest_rows = _parse_project_sample_manifest(
        io.StringIO(manifest_bytes.decode("utf-8-sig"), newline=""), path=sample_manifest,
        species=species, protocol=protocol,
    )
    return {
        "pipeline": {
            "root": str(pipeline_root.resolve()),
            "source": _pipeline_source_identity(pipeline_root),
        },
        "analysis": {
            "species": species,
            "protocol": protocol,
            "methylation_backend": backend,
            "cpg_coordinate_base": config["analysis"]["cpg_coordinate_base"],
        },
        "project": {
            "path": str(project_dir.resolve()),
            "config": config_identity,
            "sample_manifest": manifest_identity,
            **({} if protocol == "droplet" else {"barcode_map": stable_file_identity(barcode_map, hash_content=True)}),
        },
        "raw_fastqs": _fastq_inventory(
            project_dir / "01_raw",
            manifest_rows,
            protocol=protocol,
        ),
        "reference": reference_record,
        "backend_reference": _backend_reference_observation(config, reference),
        "environments": {
            "control": control_observation,
            backend: active_backend_observation,
        },
        "executables": {
            "demux_rs": _executable_identity(
                pipeline_root / "conda" / "rust" / "bin" / "demux_rs"
            ),
        },
    }


def _copy_observed_file(
    source: Path, destination: Path, expected: Mapping[str, object], input_store: Path,
    *, link_parent: Path | None = None,
) -> Mapping[str, object]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    input_store.mkdir(parents=True, exist_ok=True)
    shared = input_store / f"{expected['sha256']}.{destination.name}"
    if not shared.exists():
        with source.open("rb") as src, destination.open("wb") as dst:
            shutil.copyfileobj(src, dst, length=1024 * 1024)
            dst.flush()
            os.fsync(dst.fileno())
        copied = stable_file_identity(destination, hash_content=True)
        for key in ("size_bytes", "sha256"):
            if copied.get(key) != expected.get(key):
                raise RuntimeError(f"file changed while creating run snapshot: {source}")
        destination.chmod(0o444)
        try:
            os.link(destination, shared)
        except FileExistsError:
            pass
        destination.unlink()
    if shared.is_symlink():
        raise ValueError(f"snapshot input store may not contain a symlink: {shared}")
    if shared.stat().st_mode & 0o222:
        raise ValueError(f"snapshot input store must be read-only: {shared}")
    copied = stable_file_identity(shared, hash_content=True)
    for key in ("size_bytes", "sha256"):
        if copied.get(key) != expected.get(key):
            raise RuntimeError(f"snapshot input store identity changed: {shared}")
    destination.symlink_to(os.path.relpath(shared, link_parent if link_parent is not None else destination.parent))
    return {
        "source_path": str(source.resolve()),
        "snapshot_path": str(destination.name),
        "size_bytes": copied["size_bytes"],
        "sha256": copied["sha256"],
    }


def create_run_snapshot(
    *,
    pipeline_root: str | Path,
    project_dir: str | Path,
    config_path: str | Path,
    control_env: str | Path,
    biscuit_env: str | Path,
    bismark_env: str | Path,
    rastair_env: str | Path,
    allow_external_config: bool = False,
) -> Path:
    """按本次 run basis 创建或复用唯一不可变的启动 snapshot。"""
    pipeline = Path(pipeline_root).expanduser().resolve()
    project = Path(project_dir).expanduser().resolve()
    config = Path(config_path).expanduser().resolve()
    basis = _snapshot_basis(
        pipeline_root=pipeline,
        project_dir=project,
        config_path=config,
        control_env=Path(control_env).expanduser().resolve(),
        biscuit_env=Path(biscuit_env).expanduser().resolve(),
        bismark_env=Path(bismark_env).expanduser().resolve(),
        rastair_env=Path(rastair_env).expanduser().resolve(),
        allow_external_config=allow_external_config,
    )
    snapshot_id = sha256_bytes(canonical_json_bytes(basis))
    snapshots_root = project / "04_logs" / "provenance" / "snapshots"
    input_store = snapshots_root.parent / "inputs"
    snapshots_root.mkdir(parents=True, exist_ok=True)
    final = snapshots_root / snapshot_id
    manifest = final / "manifest.json"
    if manifest.is_file():
        load_run_snapshot(manifest)
        return manifest

    stage = Path(tempfile.mkdtemp(prefix=f".{snapshot_id}.", dir=project_temporary_parent(final, project_dir=project)))
    try:
        inputs_dir = stage / "inputs"
        project_basis = basis["project"]
        copies = {
            "config": _copy_observed_file(
                Path(str(project_basis["config"]["path"])), inputs_dir / "config.yaml",
                project_basis["config"], input_store, link_parent=final / "inputs",
            ),
            "sample_manifest": _copy_observed_file(
                Path(str(project_basis["sample_manifest"]["path"])),
                inputs_dir / "sample_manifest.tsv", project_basis["sample_manifest"], input_store, link_parent=final / "inputs",
            ),
            **({} if basis["analysis"]["protocol"] == "droplet" else {"barcode_map": _copy_observed_file(
                Path(str(project_basis["barcode_map"]["path"])),
                inputs_dir / "Barcode_Map.csv", project_basis["barcode_map"], input_store, link_parent=final / "inputs",
            )}),
        }
        payload = {
            "schema_version": SNAPSHOT_SCHEMA_VERSION,
            "snapshot_id": snapshot_id,
            "recorded_at_utc": dt.datetime.now(dt.UTC).isoformat(),
            "basis": basis,
            "input_copies": copies,
        }
        (stage / "manifest.json").write_text(
            json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        _make_tree_read_only(stage)
        try:
            os.replace(stage, final)
        except OSError:
            if not manifest.is_file():
                raise
            _remove_stage(stage)
        load_run_snapshot(manifest)
        return manifest
    except Exception:
        _remove_stage(stage)
        raise


def load_run_snapshot_basis(path: str | Path) -> Mapping[str, object]:
    manifest = Path(path).expanduser().resolve()
    try:
        payload = json.loads(manifest.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid run snapshot JSON: {manifest}: {exc}") from exc
    if not isinstance(payload, dict) or payload.get("schema_version") != SNAPSHOT_SCHEMA_VERSION:
        raise ValueError(f"unsupported run snapshot schema: {manifest}")
    basis = payload.get("basis")
    if not isinstance(basis, dict):
        raise ValueError(f"run snapshot basis is missing: {manifest}")
    expected_id = sha256_bytes(canonical_json_bytes(basis))
    if payload.get("snapshot_id") != expected_id or manifest.parent.name != expected_id:
        raise ValueError(f"run snapshot content address is invalid: {manifest}")
    return payload


def _validate_input_copy(snapshot_root: Path, record: object, expected: object, label: str) -> Path:
    if not isinstance(record, Mapping) or not isinstance(expected, Mapping):
        raise ValueError(f"run snapshot {label} copy is malformed")
    relative = str(record.get("snapshot_path", ""))
    path = snapshot_root / "inputs" / relative
    if Path(relative).name != relative or not path.is_file():
        raise FileNotFoundError(f"run snapshot {label} copy is missing: {path}")
    stat = path.stat()
    if stat.st_size != expected.get("size_bytes") or sha256_file(path) != expected.get("sha256"):
        raise ValueError(f"run snapshot {label} copy changed: {path}")
    if record.get("size_bytes") != expected.get("size_bytes") or record.get("sha256") != expected.get("sha256"):
        raise ValueError(f"run snapshot {label} copy identity disagrees with basis")
    return path


def load_run_snapshot(path: str | Path) -> Mapping[str, object]:
    """加载并校验该唯一不可变 snapshot 及其小型输入副本。"""
    manifest = Path(path).expanduser().resolve()
    payload = load_run_snapshot_basis(manifest)
    project = payload["basis"].get("project")
    copies = payload.get("input_copies")
    if not isinstance(project, Mapping) or not isinstance(copies, Mapping):
        raise ValueError(f"run snapshot input copies are malformed: {manifest}")
    for key in ("config", "sample_manifest", *([] if payload["basis"]["analysis"]["protocol"] == "droplet" else ["barcode_map"])):
        _validate_input_copy(manifest.parent, copies.get(key), project.get(key), key)
    return payload


RUN_MANIFEST_SCHEMA_VERSION = 8


DELIVERY_READY_SCHEMA_VERSION = 3


_DELIVERY_DIR_NAMES = frozenset(PUBLIC_RESULT_DIR_NAMES)


def _delivery_id(value: object) -> str:
    text = str(value)
    if len(text) != 64 or any(character not in "0123456789abcdef" for character in text):
        raise ValueError(f"invalid delivery id: {text!r}")
    return text


def _output_relative_is_owned(relative: Path) -> bool:
    return relative.parts[0] in _DELIVERY_DIR_NAMES


def _check_result_directories(root: Path, relative: Path, *, label: str) -> None:
    directory = root
    for part in (*relative.parts[:-1], None):
        if directory.is_symlink() or (directory.exists() and not directory.is_dir()):
            raise ValueError(f"{label} requires real result directories: {directory}")
        if part is not None:
            directory /= part


def _staged_file_state(path: Path, *, include_ctime: bool = False) -> dict[str, int]:
    """ctime 仅用于本次读取/哈希期间的检查；硬链接改变它，不写入恢复身份。"""
    stat = path.stat()
    state = {"device": stat.st_dev, "inode": stat.st_ino,
             "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if include_ctime:
        state["ctime_ns"] = stat.st_ctime_ns
    return state


def _staged_output_inventory(results: Path) -> list[dict[str, object]]:
    manifest = Path("QC_Results/sample_manifest.tsv")
    declared = {manifest, Path("QC_Results/multiqc_report.html")}
    _check_result_directories(results, manifest, label="staged delivery")
    if (results / manifest).is_symlink():
        raise ValueError(f"staged result may not be a symlink: {results / manifest}")
    manifest_state = _staged_file_state(results / manifest, include_ctime=True)
    with (results / manifest).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "cpg_path", "snp_path"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("final sample manifest is missing current delivery path fields")
        for row in reader:
            sample_id = str(row.get("sample_id") or "").strip()
            paths = CellOutputPaths(results.parent, results, sample_id,
                                    str(row.get("methylation_backend") or ""))
            for field, files in paths.public_result_files.items():
                value = str(row.get(field) or "").strip()
                if not value:
                    continue
                if SAMPLE_ID_RE.fullmatch(sample_id) is None:
                    raise ValueError(f"invalid delivery sample_id: {sample_id!r}")
                relative = files[0].relative_to(results)
                if value != f"03_results/{relative.as_posix()}":
                    raise ValueError(f"final sample manifest has an invalid {field}: {value!r}")
                declared.update(path.relative_to(results) for path in files)

    records: list[dict[str, object]] = []
    for relative in sorted(declared):
        _check_result_directories(results, relative, label="staged delivery")
        path = results / relative
        if path.is_symlink():
            raise ValueError(f"staged result may not be a symlink: {path}")
        if not path.is_file():
            raise FileNotFoundError(f"declared staged result is missing: {path}")
        state = _staged_file_state(path, include_ctime=True)
        if relative == manifest and state != manifest_state:
            raise ValueError(f"staged manifest changed during inventory: {path}")
        digest = sha256_file(path)
        if _staged_file_state(path, include_ctime=True) != state:
            raise ValueError(f"staged output changed while hashing: {path}")
        del state["ctime_ns"]
        records.append(
            {
                "path": relative.as_posix(),
                "size_bytes": state["size_bytes"],
                "sha256": digest,
                "stage_identity": state,
            }
        )
    return records


def _validate_inventory(
    root: Path,
    records: object,
    *,
    label: str,
    staged: bool = False,
) -> list[dict[str, object]]:
    if not isinstance(records, list) or not records:
        raise ValueError(f"{label} inventory is empty or malformed")
    validated: list[dict[str, object]] = []
    seen: set[str] = set()
    for raw in records:
        if not isinstance(raw, Mapping):
            raise ValueError(f"{label} inventory contains a malformed record")
        if (
            type(raw.get("size_bytes")) is not int
            or raw["size_bytes"] < 0
            or re.fullmatch(r"[0-9a-f]{64}", str(raw.get("sha256", ""))) is None
        ):
            raise ValueError(f"{label} inventory contains an invalid size/hash record")
        relative = _safe_relative(raw.get("path", ""), label=f"{label} output path")
        if not _output_relative_is_owned(relative):
            raise ValueError(f"{label} inventory contains an unknown output: {relative}")
        key = relative.as_posix()
        if key in seen:
            raise ValueError(f"{label} inventory repeats output: {key}")
        seen.add(key)
        _check_result_directories(root, relative, label=label)
        path = root / relative
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(f"{label} output is missing: {path}")
        stat = path.stat()
        if stat.st_size != raw.get("size_bytes"):
            raise ValueError(f"{label} output size changed: {path}")
        if staged:
            state = raw.get("stage_identity")
            if (not isinstance(state, dict) or any(type(value) is not int for value in state.values())
                    or state != _staged_file_state(path)):
                raise ValueError(f"staged output identity changed: {path}; regenerate delivery-ready")
        validated.append(dict(raw))
    return validated


def write_delivery_ready(
    out_json: str | Path,
    *,
    staged_results_dir: str | Path,
    run_snapshot_path: str | Path,
) -> Path:
    """按最终 sample manifest 声明的 CpG/SNP 与固定 QC 产物生成一次 inventory。

    未声明的 stable stage 文件不进入本次发布；只保存 sealed inventory 与 run snapshot 的绑定，
    结构上固定的路径（稳定 stage、公开目录）不落盘。
    """
    output = Path(out_json).expanduser().resolve()
    stage_input = Path(staged_results_dir).expanduser().absolute()
    _check_result_directories(stage_input.parent.parent, Path("02_work/results_stage/_"), label="staged delivery")
    stage_results = stage_input.resolve()
    project = stage_results.parent.parent
    expected_publish_root = project / "05_tmp" / "finalize" / output.parent.name
    if output.parent != expected_publish_root.resolve(strict=False):
        raise ValueError(f"delivery-ready path is outside its publication transaction: {output}")
    expected_stage = (project / "02_work" / "results_stage").resolve(strict=False)
    if stage_results != expected_stage:
        raise ValueError(f"biological result stage must be stable: {stage_results} != {expected_stage}")

    snapshot_path = Path(run_snapshot_path).expanduser().resolve()
    snapshot = load_run_snapshot(snapshot_path)
    basis = snapshot.get("basis")
    if not isinstance(basis, Mapping):
        raise ValueError("run snapshot basis is missing")
    project_record = basis.get("project")
    if not isinstance(project_record, Mapping):
        raise ValueError("run snapshot delivery identity is incomplete")
    snapshot_project = Path(str(project_record.get("path", ""))).resolve()
    if snapshot_project != project.resolve():
        raise ValueError("run snapshot project and delivery project disagree")
    delivery_id = _delivery_id(output.parent.name)
    if snapshot.get("snapshot_id") != delivery_id:
        raise ValueError("delivery generation must use the run snapshot id")

    payload: dict[str, object] = {
        "schema_version": DELIVERY_READY_SCHEMA_VERSION,
        "delivery_id": delivery_id,
        "recorded_at_utc": _utc_now(),
        "project_dir": str(project.resolve()),
        "run_snapshot_path": str(snapshot_path),
        "outputs": _staged_output_inventory(stage_results),
    }
    atomic_write_json(output, payload)
    return output


def load_delivery_ready(
    path: str | Path, *, validate_staged: bool = True
) -> dict[str, object]:
    manifest = Path(path).expanduser().resolve()
    payload = _read_json(manifest, "delivery-ready")
    if payload.get("schema_version") != DELIVERY_READY_SCHEMA_VERSION:
        raise ValueError(f"unsupported delivery-ready schema: {manifest}")
    project = Path(str(payload.get("project_dir", ""))).expanduser().resolve()


    stage = project / "02_work" / "results_stage"
    _check_result_directories(project, Path("02_work/results_stage/_"), label="staged delivery")
    delivery_id = _delivery_id(payload.get("delivery_id", ""))
    expected_root = project / "05_tmp" / "finalize" / delivery_id
    if manifest.parent != expected_root.resolve(strict=False):
        raise ValueError("delivery-ready location and delivery_id disagree")
    if validate_staged:
        _validate_inventory(stage, payload.get("outputs"), label="staged delivery", staged=True)
    return payload


def _manifest_path(project: Path) -> Path:
    return project / "03_results" / "run_manifest.json"


def load_run_manifest(
    path_or_project: str | Path,
    *,
    validate_outputs: bool = True,
) -> dict[str, object]:
    path = Path(path_or_project).expanduser().resolve()
    manifest = path if path.name == "run_manifest.json" else _manifest_path(path)
    payload = _read_json(manifest, "run manifest")
    if payload.get("schema_version") != RUN_MANIFEST_SCHEMA_VERSION:
        raise ValueError(f"unsupported run manifest schema: {manifest}")
    project = Path(str(payload.get("project_dir", ""))).expanduser().resolve()
    if manifest != _manifest_path(project).resolve(strict=False):
        raise ValueError("run manifest project and location disagree")
    if validate_outputs:
        _validate_inventory(
            project / "03_results",
            payload.get("outputs"),
            label="published delivery",
        )
    return payload


_PUBLIC_DIRNAME = "03_results"


def _publish_hardlink(stage: Path, public: Path, relative: Path) -> None:
    destination = public / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp = project_temporary_parent(destination) / f".{destination.name}.tmp.{os.getpid()}.{uuid.uuid4().hex}"
    try:
        os.link(stage / relative, tmp)
        os.replace(tmp, destination)
    finally:
        tmp.unlink(missing_ok=True)


def _prune_stale_public_files(public: Path, allowed: set[str]) -> None:
    for path in sorted(public.rglob("*"), reverse=True):
        relative = path.relative_to(public)
        if path.is_symlink() or path.is_file():
            if relative.as_posix() not in allowed:
                _remove_path(path)
        elif path.is_dir() and not any(path.iterdir()):
            _remove_path(path)


def _run_manifest_payload(ready: Mapping[str, object], snapshot: Mapping[str, object]) -> dict[str, object]:
    snapshot_path = Path(str(ready["run_snapshot_path"])).resolve()
    basis = snapshot["basis"]
    if not isinstance(basis, Mapping):
        raise ValueError("run snapshot basis is missing")
    pipeline_record = basis.get("pipeline")
    analysis_record = basis.get("analysis")
    reference_record = basis.get("reference")
    if not all(isinstance(item, Mapping) for item in (pipeline_record, analysis_record, reference_record)):
        raise ValueError("run snapshot provenance identity is incomplete")
    project_record = basis.get("project")
    if (snapshot.get("snapshot_id") != ready["delivery_id"]
            or not isinstance(project_record, Mapping)
            or Path(str(project_record.get("path", ""))).resolve() != Path(str(ready["project_dir"])).resolve()):
        raise ValueError("delivery-ready and run snapshot generation/project disagree")
    return {
        "schema_version": RUN_MANIFEST_SCHEMA_VERSION,
        "project_dir": str(ready["project_dir"]),
        "recorded_at_utc": _utc_now(),
        "status": "complete",
        "delivery": {
            "id": ready["delivery_id"],
        },
        "pipeline": {
            "name": "Alopex",
            "source_root": str(pipeline_record.get("root", "")),
        },
        "run": {
            "protocol": str(analysis_record.get("protocol", "")),
            "species": str(analysis_record.get("species", "")),
            "methylation_backend": str(analysis_record.get("methylation_backend", "")),
            "cpg_representation": cpg_representation(analysis_record["cpg_coordinate_base"]),
            "snapshot_id": snapshot["snapshot_id"],
            "snapshot_path": str(snapshot_path),
            "reference_path": str(Path(str(reference_record.get("resolved_path", ""))).resolve()),
        },
        "outputs": [{key: item[key] for key in ("path", "size_bytes", "sha256")}
                    for item in ready["outputs"]],
    }


def publish_delivery(ready_path: str | Path) -> Path:
    """发布一笔 sealed 结果树到真实目录 03_results/。

    事务语义是 invalidate-first / commit-last：任何公开文件发生变化之前先删除
    旧 run_manifest.json，安装并校验全部公开文件之后才写回新 marker——marker 的
    存在即代表一代完整交付，中断只会留下无 marker 的待重发状态，不会留下
    「混合代结果 + complete marker」。交付内容全部来自 stable stage 的硬链接，
    逐文件原子落地（os.replace）；发布中断不损坏 stable stage，重新运行
    launcher 会直接重走本流程。生产调用必须位于 Snakemake 持有项目锁的
    onsuccess 回调内，锁在发布与 BAM/demux 清理完成后释放。
    """
    ready_manifest = Path(ready_path).expanduser().resolve()
    ready = load_delivery_ready(ready_manifest, validate_staged=False)
    project = Path(str(ready["project_dir"])).resolve()
    public = project / _PUBLIC_DIRNAME
    stage = project / "02_work" / "results_stage"

    outputs = _validate_inventory(stage, ready.get("outputs"), label="staged delivery", staged=True)
    snapshot_path = Path(str(ready["run_snapshot_path"])).resolve()
    snapshot = load_run_snapshot(snapshot_path)
    payload = _run_manifest_payload(ready, snapshot)
    config = load_project_config(
        snapshot_path.parent / "inputs" / snapshot["input_copies"]["config"]["snapshot_path"]
    )
    for record in outputs:
        _check_result_directories(public, Path(record["path"]), label="published delivery")

    marker = _manifest_path(project)
    _remove_path(marker)

    public.mkdir(parents=True, exist_ok=True)

    allowed: set[str] = set()
    for record in outputs:
        relative = _safe_relative(record["path"], label="published output path")
        allowed.add(relative.as_posix())
        _publish_hardlink(stage, public, relative)

    _prune_stale_public_files(public, allowed)
    _validate_inventory(stage, outputs, label="staged delivery", staged=True)
    _validate_inventory(public, outputs, label="published delivery")
    for record in outputs:
        relative = Path(record["path"])
        if not os.path.samefile(stage / relative, public / relative):
            raise ValueError(f"published delivery is not the staged hardlink: {relative}")

    atomic_write_json(marker, payload)

    if not config["retention"]["keep_final_bam"]:
        try:
            for path in _retained_bam_paths(project, config):
                path.unlink(missing_ok=True)
        except (OSError, ValueError) as exc:
            logger.warning("Could not fully reclaim unretained BAMs: %s", exc)

    finalize_root = project / "05_tmp" / "finalize"
    if ready_manifest.parent.parent == finalize_root:
        _remove_path(ready_manifest.parent)
    clean_demux_scratch_after_publish(project)
    return marker


def _retained_bam_paths(project: Path, config: Mapping[str, object]) -> tuple[Path, ...]:
    layout = ProjectLayout(project / "02_work", project / "03_results", get_backend(config))
    for directory in (layout.intermediate_dir, layout.intermediate_dir / layout.backend):
        if directory.is_symlink():
            raise ValueError(f"Refusing to access retained BAMs through a symlink: {directory}")
    paths = []
    with (layout.result_root / "QC_Results/sample_manifest.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"sample_id", "analysis_status"} <= set(reader.fieldnames or ()):
            raise ValueError("Final manifest is missing BAM retention identity/status")
        for row in reader:
            if row["analysis_status"] != "PASS":
                continue
            if SAMPLE_ID_RE.fullmatch(row["sample_id"]) is None:
                raise ValueError(f"Invalid retained BAM sample_id: {row['sample_id']!r}")
            artifacts = cell_artifacts(layout.cell(row["sample_id"]), config["analysis"]["protocol"],
                                       generate_snp=False, keep_final_bam=True)
            paths.extend(artifacts.retained_bams)
    return tuple(paths)


def clean_demux_scratch_after_publish(project_dir: str | Path) -> None:
    """在已完成发布后先清除 demux 提交边界；失败时保留 FASTQ，拒绝经目录符号链接清理。"""
    project = Path(project_dir).expanduser().resolve()
    scratch = project / "02_work" / "demux"
    if not scratch.is_dir():
        return
    if (project / "02_work").is_symlink() or scratch.is_symlink():
        logger.warning("Refusing to reclaim demux scratch through a symlink: %s", scratch)
        return
    try:
        _remove_path(scratch / "state")
    except OSError as exc:
        logger.warning("Demux state could not be removed; preserving FASTQs: %s: %s", scratch, exc)
        return
    try:
        _remove_path(scratch)
    except OSError as exc:
        logger.warning("Could not fully reclaim demux scratch: %s: %s", scratch, exc)


def _unsealed_project_state_paths(project: Path) -> tuple[str, ...]:
    state_paths: list[str] = []
    for relative in ("02_work", "03_results"):
        path = project / relative
        if path.is_dir() and any(path.iterdir()):
            state_paths.append(relative)
    runtime_tmp = project / "05_tmp"
    if runtime_tmp.is_dir():
        state_paths.extend(
            f"05_tmp/{child.name}"
            for child in sorted(runtime_tmp.iterdir(), key=lambda item: item.name)
            if child.name != "xdg_cache"
        )
    return tuple(state_paths)


def project_status(
    *,
    project_dir: str | Path,
    current_snapshot: str | Path | None = None,
) -> dict[str, object]:
    """检查 resumable/public 项目状态，不因 provenance identity 强制重算。

    provenance snapshot 只是审计记录，不决定既有 Snakemake 工作可否复用；
    复用权威是 Snakemake 自身的 DAG、input、params、code、软件环境与 incomplete-output 逻辑。
    complete_current 比较 launcher 本次 snapshot 的内容地址并检查声明的保留 BAM/索引，不重算 basis。
    """
    project = Path(project_dir).expanduser().resolve()
    public_manifest = _manifest_path(project)
    if public_manifest.is_file():
        try:
            manifest = load_run_manifest(public_manifest)
        except (FileNotFoundError, ValueError, OSError) as exc:
            return {
                "status": "resumable",
                "reason": f"public manifest will be replaced: {exc}",
            }
        if manifest.get("status") == "complete":
            if current_snapshot is not None:
                try:
                    snapshot = load_run_snapshot(Path(current_snapshot))
                    if manifest["run"]["snapshot_id"] == snapshot["snapshot_id"]:
                        config_path = (
                            Path(current_snapshot).resolve().parent / "inputs"
                            / snapshot["input_copies"]["config"]["snapshot_path"]
                        )
                        config = load_project_config(config_path)
                        keep_bam = config["retention"]["keep_final_bam"]
                        for path in _retained_bam_paths(project, config):
                            if keep_bam and not path.is_file():
                                return {
                                    "status": "resumable",
                                    "reason": f"requested retained BAM/index is missing: {path}",
                                }
                            if not keep_bam and (path.exists() or path.is_symlink()):
                                return {
                                    "status": "resumable",
                                    "reason": f"unretained BAM/index awaits cleanup under the Snakemake lock: {path}",
                                }
                        return {
                            "status": "complete_current",
                            "manifest": str(public_manifest),
                            "delivery_id": manifest.get("delivery", {}).get("id", ""),
                        }
                except (FileNotFoundError, ValueError, RuntimeError, OSError, KeyError, TypeError):
                    pass
            return {
                "status": "resumable",
                "reason": "completed delivery differs from current inputs/code/environment; reusing compatible Snakemake outputs",
            }
        return {
            "status": "resumable",
            "reason": "public manifest status is not complete; normal execution will repair it",
        }


    if _unsealed_project_state_paths(project):
        return {
            "status": "resumable",
            "reason": "existing Snakemake work has no provenance snapshot; preserving and reusing it",
        }
    return {"status": "new"}
