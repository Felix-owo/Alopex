"""Backend科学处理、原生报告与整数计数；不调用Doctor或交付。"""
from __future__ import annotations

import csv
import gzip
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict
from collections.abc import Mapping, Sequence
from contextlib import ExitStack
from pathlib import Path
from typing import Any, TextIO

from .config import (
    FILTER_SUMMARY_FIELDS,
    FILTER_SUMMARY_SCHEMA_VERSION,
    SAMPLE_ID_RE,
    _integer,
)
from .io import (
    _close_native_stream,
    _read_json,
    atomic_write_json,
)


# 科学处理：共用统计

_COUNT = r"([0-9][0-9,]*)"


def _required_count(text: str, label: str, patterns: Sequence[str]) -> int:
    matches: list[int] = []
    for pattern in patterns:
        for match in re.finditer(pattern, text, re.IGNORECASE | re.MULTILINE):
            matches.append(int(match.group(1).replace(",", "")))
    if not matches:
        raise ValueError(f"Bismark report is missing required field: {label}")
    if len(set(matches)) != 1:
        raise ValueError(f"Bismark report has conflicting values for {label}: {matches}")
    return matches[0]


def _optional_count(text: str, label: str, patterns: Sequence[str], *, default: int = 0) -> int:
    matches: list[int] = []
    for pattern in patterns:
        for match in re.finditer(pattern, text, re.IGNORECASE | re.MULTILINE):
            matches.append(int(match.group(1).replace(",", "")))
    if not matches:
        return default
    if len(set(matches)) != 1:
        raise ValueError(f"Bismark report has conflicting values for {label}: {matches}")
    return matches[0]


def _required_float(text: str, label: str, pattern: str) -> float:
    matches = [
        float(match.group(1))
        for match in re.finditer(pattern, text, re.IGNORECASE | re.MULTILINE)
    ]
    if not matches:
        raise ValueError(f"Bismark report is missing required field: {label}")
    if any(not math.isclose(matches[0], value, abs_tol=1e-9) for value in matches[1:]):
        raise ValueError(f"Bismark report has conflicting values for {label}: {matches}")
    return matches[0]


def _count_bam_pairs(bam: Path, samtools: str) -> int:
    result = subprocess.run(
        [samtools, "view", "-c", str(bam)],
        check=True, text=True, stdout=subprocess.PIPE,
    )
    records = int(result.stdout.strip())
    if records < 0 or records % 2:
        raise ValueError("Bismark processed BAM contains an odd or negative alignment count")
    return records // 2


def _parse_cutadapt_metrics(path: Path) -> dict[str, int | float]:
    payload = _read_json_if_present(path)
    counts = payload.get("read_counts")
    if not isinstance(counts, dict):
        raise ValueError(f"Cutadapt JSON lacks read_counts: {path}")
    values: dict[str, int] = {}
    for key in ("input", "output", "read1_with_adapter", "read2_with_adapter"):
        value = counts.get(key)
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f"Cutadapt JSON has invalid read_counts.{key}: {path}")
        values[key] = value
    if values["input"] <= 0:
        raise ValueError(f"Cutadapt JSON has zero input read pairs: {path}")
    if values["read1_with_adapter"] > values["input"] or values["read2_with_adapter"] > values["input"]:
        raise ValueError(f"Cutadapt adapter count exceeds input read pairs: {path}")
    if values["output"] > values["input"]:
        raise ValueError(f"Cutadapt output exceeds input read pairs: {path}")
    return {
        "cutadapt_input_pairs": values["input"],
        "trimmed_pairs": values["output"],
        "r1_adapter_pct": 100.0 * values["read1_with_adapter"] / values["input"],
        "r2_adapter_pct": 100.0 * values["read2_with_adapter"] / values["input"],
    }


def _read_json_if_present(path: Path) -> dict[str, object]:
    if not path.is_file():
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise ValueError(f"Invalid JSON file: {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"JSON root must be an object: {path}")
    return value


def _read_high_cph_summary(path: Path, sample: str, backend: str) -> dict[str, object]:
    if backend == "bismark":
        return _read_filter_summary(path, sample)
    payload = _read_json_if_present(path)
    common = {
        "schema_version",
        "backend",
        "sample",
        "metric_unit",
        "candidate_read_pairs",
        "flagged_read_pairs",
        "high_cph_fraction",
        "excluded_read_pairs",
        "filter_applied",
        "threshold",
        "protocol",
        "high_cph_role",
    }
    backend_fields = {
        "missing_metric_rows", "observed_metric_rows", "flagged_metric_rows",
        "observed_rows", "excluded_rows", "metric", "comparison",
    }
    expected = common | backend_fields
    if set(payload) != expected:
        raise ValueError(
            f"High-CpH summary schema mismatch for {sample!r}: "
            f"expected={sorted(expected)}, observed={sorted(payload)}"
        )
    if (
        payload["schema_version"] != 2
        or payload["backend"] != backend
        or payload["sample"] != sample
        or payload["metric_unit"] != "read_pairs"
        or payload["protocol"] not in {"cabernet", "srd", "droplet"}
        or payload["high_cph_role"]
        not in {"cdna_contamination", "residual_high_cph"}
    ):
        raise ValueError(f"High-CpH summary identity is invalid: {path}")
    expected_role = {
        "cabernet": "cdna_contamination",
        "srd": "residual_high_cph",
    }[str(payload["protocol"])]
    if payload["high_cph_role"] != expected_role:
        raise ValueError(f"High-CpH summary protocol/role contract is invalid: {path}")
    counts: dict[str, int] = {}
    for field in (
        "candidate_read_pairs",
        "flagged_read_pairs",
        "excluded_read_pairs",
    ):
        value = payload[field]
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f"High-CpH summary has invalid {field}: {path}")
        counts[field] = value
    if counts["flagged_read_pairs"] > counts["candidate_read_pairs"]:
        raise ValueError(f"High-CpH flagged pairs exceed candidate pairs: {path}")
    if not isinstance(payload["filter_applied"], bool):
        raise ValueError(f"High-CpH summary filter_applied is invalid: {path}")
    fraction = payload["high_cph_fraction"]
    if isinstance(fraction, bool) or not isinstance(fraction, (int, float)):
        raise ValueError(f"High-CpH summary fraction is invalid: {path}")
    expected_fraction = (
        counts["flagged_read_pairs"] / counts["candidate_read_pairs"]
        if counts["candidate_read_pairs"]
        else 0.0
    )
    if not math.isclose(float(fraction), expected_fraction, abs_tol=1e-12):
        raise ValueError(f"High-CpH summary fraction is inconsistent: {path}")
    threshold = payload["threshold"]
    if (
        isinstance(threshold, bool)
        or not isinstance(threshold, (int, float))
        or not 0 <= float(threshold) <= 1
    ):
        raise ValueError(f"High-CpH summary threshold is invalid: {path}")
    row_counts: dict[str, int] = {}
    for field in (
        "missing_metric_rows",
        "observed_metric_rows",
        "flagged_metric_rows",
        "observed_rows",
        "excluded_rows",
    ):
        value = payload[field]
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f"High-CpH summary has invalid {field}: {path}")
        row_counts[field] = value
    if (
        row_counts["missing_metric_rows"]
        + row_counts["observed_metric_rows"]
        + row_counts["excluded_rows"]
        != row_counts["observed_rows"]
        or row_counts["flagged_metric_rows"]
        > row_counts["observed_metric_rows"]
        or counts["candidate_read_pairs"] > row_counts["observed_metric_rows"]
        or counts["flagged_read_pairs"] > row_counts["flagged_metric_rows"]
        or counts["excluded_read_pairs"] > row_counts["excluded_rows"]
        or payload["metric"] != "CpH_retention"
        or payload["comparison"] != ">"
    ):
        raise ValueError(f"BISCUIT High-CpH row accounting is inconsistent: {path}")
    return payload


# 科学处理：BISCUIT

BSCONV_P_COLUMNS = 9


def read_pileup_variant_count(path: str | Path, sample: str) -> int:
    """读取当前细胞的 pileup 计数；非法计数或身份不得触发空结果分支。"""
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if (not isinstance(value, dict) or type(value.get("schema_version")) is not int
            or value["schema_version"] != 1 or value.get("sample") != sample):
        raise ValueError(f"BISCUIT pileup stats schema/sample mismatch: {path}")
    for field in ("mapped_reads", "variant_records"):
        if type(value.get(field)) is not int or value[field] < 0:
            raise ValueError(f"BISCUIT pileup stats requires a nonnegative integer {field}: {path}")
    if (type(value.get("empty")) is not bool
            or value["empty"] != (value["variant_records"] == 0)
            or (value["mapped_reads"] == 0 and value["variant_records"] != 0)):
        raise ValueError(f"BISCUIT pileup stats empty/count accounting mismatch: {path}")
    return value["variant_records"]


def open_text(path: Path) -> TextIO:
    if str(path) == "-":
        return sys.stdin
    return path.open("r", encoding="utf-8")


def normalize_read_name(value: str) -> str:

    return str(value).strip().split()[0]


def indexed_excluded_contigs(bam_path: Path, samtools: str, excluded: set[str]) -> list[tuple[str, int]]:
    """从原生索引读取精确 contig 名及全长，供 BED 查询使用。"""
    proc = subprocess.run(
        [samtools, "idxstats", str(bam_path)],
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or f"samtools idxstats failed for {bam_path}")
    contigs = []
    for line in proc.stdout.splitlines():
        fields = line.split("\t")
        if fields and fields[0] in excluded:
            if len(fields) != 4 or not re.fullmatch(r"[0-9]+", fields[1]) or int(fields[1]) < 1:
                raise ValueError(f"Invalid samtools idxstats contig length for {fields[0]!r}: {bam_path}")
            contigs.append((fields[0], int(fields[1])))
    return contigs


def excluded_reads_from_bam(bam_path: Path, samtools: str, excluded: set[str]) -> set[str]:
    """经 samtools idxstats/view 流式收集排除 contig 上的 QNAME，不物化完整 SAM。"""
    contigs = indexed_excluded_contigs(bam_path, samtools, excluded)
    if not contigs:
        return set()

    # CLI region 会把 -h 当选项、chr1:10-20 当区间，连花括号也不能转义所有合法名称。
    # BED 的第一列按名称精确关联；通常覆盖完整 contig 并保留索引查询。
    # 合法的 # 前缀在 BED 中却是注释，此时只能流式按 SAM 的精确 RNAME 选择。
    # stdout 流式读取，stderr 落临时文件，避免双 PIPE 互相等待。
    with ExitStack() as temporary, tempfile.TemporaryFile(mode="w+t", errors="replace") as diagnostics:
        command = [samtools, "view"]
        if not any(name.startswith('#') for name, _ in contigs):
            regions = temporary.enter_context(tempfile.NamedTemporaryFile(mode="w+t", encoding="utf-8", suffix=".bed"))
            for name, length in contigs:
                regions.write(f"{name}\t0\t{length}\n")
            regions.flush()
            command.extend(["-M", "-L", regions.name])
        command.extend(["--", str(bam_path)])
        proc = subprocess.Popen(
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=diagnostics,
            bufsize=1,
        )
        assert proc.stdout is not None
        names = set()
        try:
            for line in proc.stdout:
                if not line:
                    continue
                fields = line.rstrip("\r\n").split("\t", 3)
                if len(fields) < 3:
                    raise ValueError(f"Invalid samtools view record for {bam_path}: {line!r}")
                if fields[2] in excluded:
                    names.add(normalize_read_name(fields[0]))
            returncode = proc.wait()
            if returncode != 0:
                diagnostics.seek(0)
                raise RuntimeError(
                    f"samtools view failed for {bam_path} (exit {returncode}): "
                    f"{diagnostics.read().strip()}"
                )
        finally:
            _close_native_stream(proc)
    return names


def _high_cph_multiqc_payload(
    sample_name: str,
    summary: Mapping[str, object],
    protocol: str,
    high_cph_role: str,
) -> dict[str, object]:
    if summary.get("backend") == "bismark":
        return {
            "id": "dna_pipeline_bismark_nonconversion",
            "section_name": "Bismark non-conversion filtering",
            "plot_type": "table", "data": {sample_name: dict(summary)},
        }
    fraction = float(summary["high_cph_fraction"])
    row = {
        "Candidate_Read_Pairs": int(summary["candidate_read_pairs"]),
        "Flagged_Read_Pairs": int(summary["flagged_read_pairs"]),
        "High_CpH_Fraction": fraction,
        "Excluded_Read_Pairs": int(summary["excluded_read_pairs"]),
        "Threshold": float(summary["threshold"]),
    }
    interpretation = (
        "candidate cDNA contamination"
        if high_cph_role == "cdna_contamination"
        else "residual high-CpH/non-converted DNA reads"
    )
    return {
        "id": "bsconv_filter",
        "section_name": "High-CpH Read Assessment",
        "description": (
            f"Protocol={protocol}; high-CpH reads are interpreted as {interpretation}. "
            "They are excluded from DNA methylation pileup."
        ),
        "plot_type": "table",
        "pconfig": {
            "id": "bsconv_filter_table",
            "title": "High-CpH Read Assessment",
            "col1_header": "Sample",
        },
        "data": {sample_name: row},
    }


def write_multiqc(
    path: Path,
    sample_name: str,
    summary: Mapping[str, object],
    protocol: str,
    high_cph_role: str,
) -> None:
    """把 high-CpH 评估汇总写成 MultiQC custom-content table JSON。"""
    obj = _high_cph_multiqc_payload(sample_name, summary, protocol, high_cph_role)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2) + "\n")


def _parse_excluded_contigs_json(value: str) -> set[str]:
    try:
        raw = json.loads(value)
    except json.JSONDecodeError as exc:
        raise ValueError(f"--excluded-contigs-json is invalid JSON: {exc}") from exc
    if not isinstance(raw, list) or any(
        not isinstance(value, str) or re.fullmatch(r"\S+", value) is None or "\x00" in value
        for value in raw
    ):
        raise ValueError("--excluded-contigs-json must be a JSON array of non-empty literal strings without whitespace or NUL")
    if len(set(raw)) != len(raw):
        raise ValueError("--excluded-contigs-json must not repeat a contig")
    return set(raw)


def _sort_unique_count(source: Path, destination: Path) -> int:
    env = dict(os.environ)
    env["LC_ALL"] = "C"
    fd, temporary = tempfile.mkstemp(prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent)
    os.close(fd)
    temporary_path = Path(temporary)
    try:
        proc = subprocess.run(
            ["sort", "-u", str(source), "-o", str(temporary_path)],
            text=True,
            capture_output=True,
            env=env,
            check=False,
        )
        if proc.returncode != 0:
            raise ValueError(proc.stderr.strip() or f"sort -u failed for {source}")
        with temporary_path.open("r", encoding="utf-8") as handle:
            count = sum(1 for line in handle if line.strip())
        os.replace(temporary_path, destination)
        return count
    finally:
        temporary_path.unlink(missing_ok=True)


def split_pairs(
    source: Path, bio: Path, excluded: Path, excluded_contigs: set[str]
) -> tuple[int, int]:
    """按原始顺序拆分 paired BAM，校验相邻 R1/R2 并在写出时累计两路 pair 数。"""
    import pysam

    counts = [0, 0]
    with pysam.AlignmentFile(str(source), "rb") as reader, \
         pysam.AlignmentFile(str(bio), "wb", template=reader, format_options=[b"level=1"]) as bio_out, \
         pysam.AlignmentFile(str(excluded), "wb", template=reader, format_options=[b"level=1"]) as excluded_out:
        excluded_ids = {reader.get_tid(name) for name in excluded_contigs} - {-1}
        records = reader.fetch(until_eof=True)
        for left in records:
            right = next(records, None)
            if right is None:
                raise ValueError("Odd number of alignment records in paired Bismark BAM")
            if left.query_name != right.query_name:
                raise ValueError("paired Bismark BAM is not mate-adjacent")
            if not (
                left.is_paired and right.is_paired
                and left.is_read1 != left.is_read2
                and right.is_read1 != right.is_read2
                and left.is_read1 == right.is_read2
                and left.is_read2 == right.is_read1
            ):
                raise ValueError(f"Bismark records are not an R1/R2 pair: {left.query_name!r}")
            route = int(bool(excluded_ids.intersection((
                left.reference_id, left.next_reference_id,
                right.reference_id, right.next_reference_id,
            ))))
            destination = excluded_out if route else bio_out
            destination.write(left)
            destination.write(right)
            counts[route] += 1
    return counts[0], counts[1]


def run_bsconv_extraction(args) -> None:
    """主流程：流式解析 bsconv TSV、排除指定 contig、按阈值划分高 retention QNAME，
    外部 sort -u 按 read-pair 去重后写出 high-read-names、summary 与 MultiQC JSON。"""
    if not 0 <= args.threshold <= 1:
        raise ValueError("--threshold must be between 0 and 1")

    excluded = _parse_excluded_contigs_json(args.excluded_contigs_json)
    excluded_read_names: set[str] | None = None
    total = excluded_rows = metric_rows = missing_metric = 0
    high_metric_rows = 0
    args.high_read_names.parent.mkdir(parents=True, exist_ok=True)


    raw_high_path: Path | None = None
    raw_candidate_path: Path | None = None
    raw_excluded_path: Path | None = None
    candidate_unique_path: Path | None = None
    excluded_unique_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=args.high_read_names.name + ".",
            suffix=".unsorted",
            dir=args.high_read_names.parent,
            delete=False,
        ) as raw_high_handle, tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=args.high_read_names.name + ".",
            suffix=".candidate.unsorted",
            dir=args.high_read_names.parent,
            delete=False,
        ) as raw_candidate_handle, tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=args.high_read_names.name + ".",
            suffix=".excluded.unsorted",
            dir=args.high_read_names.parent,
            delete=False,
        ) as raw_excluded_handle:
            raw_high_path = Path(raw_high_handle.name)
            raw_candidate_path = Path(raw_candidate_handle.name)
            raw_excluded_path = Path(raw_excluded_handle.name)
            with open_text(args.bsconv_tsv) as handle:
                reader = csv.reader(handle, delimiter="\t")
                for row in reader:
                    if not row or row[0].startswith("#"):
                        continue
                    if len(row) != BSCONV_P_COLUMNS or not all(
                        re.fullmatch(r"[0-9]+", value) for value in row[:8]
                    ) or not row[8].strip():
                        raise ValueError(
                            f"Unexpected biscuit bsconv -p row at {args.bsconv_tsv}:{reader.line_num} "
                            f"(expect 8 nonnegative integer cytosine counts + QNAME): {row!r}"
                        )

                    read_name = normalize_read_name(row[8])
                    if excluded_read_names is None:
                        excluded_read_names = (
                            excluded_reads_from_bam(
                                args.alignment_bam, args.samtools, excluded
                            )
                            if excluded
                            else set()
                        )
                    total += 1
                    if read_name in excluded_read_names:
                        excluded_rows += 1
                        raw_excluded_handle.write(read_name + "\n")
                        continue


                    retained = int(row[0]) + int(row[2]) + int(row[6])
                    converted = int(row[1]) + int(row[3]) + int(row[7])
                    sites = retained + converted
                    if sites == 0:
                        missing_metric += 1
                        continue
                    metric_rows += 1
                    raw_candidate_handle.write(read_name + "\n")
                    if retained / sites > args.threshold:
                        high_metric_rows += 1
                        raw_high_handle.write(read_name + "\n")


        candidate_unique_path = raw_candidate_path.with_suffix(".unique")
        excluded_unique_path = raw_excluded_path.with_suffix(".unique")
        candidate_read_pairs = _sort_unique_count(raw_candidate_path, candidate_unique_path)
        flagged_read_pairs = _sort_unique_count(raw_high_path, args.high_read_names)
        excluded_read_pairs = _sort_unique_count(raw_excluded_path, excluded_unique_path)
        candidate_unique_path.unlink(missing_ok=True)
        excluded_unique_path.unlink(missing_ok=True)
    finally:
        for temporary in (
            raw_high_path,
            raw_candidate_path,
            raw_excluded_path,
            candidate_unique_path,
            excluded_unique_path,
        ):
            if temporary is not None:
                temporary.unlink(missing_ok=True)

    summary = {
        "schema_version": 2,
        "backend": "biscuit",
        "sample": args.sample_name,
        "metric_unit": "read_pairs",
        "candidate_read_pairs": candidate_read_pairs,
        "flagged_read_pairs": flagged_read_pairs,
        "high_cph_fraction": (
            flagged_read_pairs / candidate_read_pairs if candidate_read_pairs else 0.0
        ),
        "excluded_read_pairs": excluded_read_pairs,
        "filter_applied": True,
        "missing_metric_rows": missing_metric,
        "observed_metric_rows": metric_rows,
        "flagged_metric_rows": high_metric_rows,
        "observed_rows": total,
        "excluded_rows": excluded_rows,
        "threshold": args.threshold,
        "metric": "CpH_retention",
        "comparison": ">",
        "protocol": args.protocol,
        "high_cph_role": args.high_cph_role,
    }
    args.summary_json.parent.mkdir(parents=True, exist_ok=True)
    args.summary_json.write_text(json.dumps(summary, indent=2) + "\n")
    if args.multiqc_json:
        write_multiqc(
            args.multiqc_json,
            args.sample_name,
            summary,
            args.protocol,
            args.high_cph_role,
        )


DUPSIFTER_STAT_FIELDS: tuple[tuple[str, str], ...] = (
    ("number of individual reads processed", "Reads_Processed"),
    ("number of reads with both reads mapped", "Both_Mapped"),
    ("number of reads with only one read mapped to the forward strand", "One_Mapped_Forward"),
    ("number of reads with only one read mapped to the reverse strand", "One_Mapped_Reverse"),
    ("number of reads with both reads marked as duplicates", "Dup_Both"),
    ("number of reads on the forward strand marked as duplicates", "Dup_Forward"),
    ("number of reads on the reverse strand marked as duplicates", "Dup_Reverse"),
    ("number of individual primary-alignment reads", "Primary_Reads"),
    ("number of individual secondary- and supplementary-alignment reads", "Secondary_Supplementary"),
    ("number of reads with no reads mapped", "No_Mapped"),
    ("number of reads with no primary reads", "No_Primary"),
)


DUPSIFTER_MQC_FIELDS = tuple(field for _label, field in DUPSIFTER_STAT_FIELDS)


def write_dupsifter_mqc(
    stat_path: str,
    out_json: str,
    demux_sample: str,
    *,
    metadata: Mapping[str, str],
) -> None:
    """从 DAG 提交的 manifest 行写 Dupsifter custom content。"""
    row = {str(key): str(value or "") for key, value in metadata.items()}
    if row.get("demux_sample") != demux_sample:
        raise ValueError(
            f"Dupsifter metadata sample mismatch: {row.get('demux_sample')!r} "
            f"!= {demux_sample!r}"
        )
    if str(row.get("downstream_dna", "")).strip().lower() != "true":
        raise ValueError(
            f"Dupsifter sample {demux_sample!r} belongs to non-DNA route "
            f"{row.get('route_id')!r}; refusing downstream DNA QC."
        )
    if row.get("dna_status") != "Pass":
        raise ValueError(
            f"Dupsifter sample {demux_sample!r} has dna_status={row.get('dna_status')!r}; "
            "refusing downstream DNA QC."
        )

    plate_id = row.get("plate_id", "")
    cell_order = row.get("cell_order", "")
    dna_bc = row.get("dna_barcode", "")
    raw_sample = row.get("raw_sample", "")
    pipeline_mode = row.get("pipeline_mode", "")
    route_id = row.get("route_id", "")

    metrics = {}
    mapping = dict(DUPSIFTER_STAT_FIELDS)
    stat = Path(stat_path)
    if stat.exists():
        for line_no, line in enumerate(stat.read_text().splitlines(), start=1):
            if "]" in line and ":" in line:
                content = line.split("]", 1)[1].strip()
                if ":" in content:
                    key, value = content.split(":", 1)
                    key = key.strip()
                    value = value.strip().replace(",", "")
                    if key in mapping:
                        try:
                            metrics[mapping[key]] = int(value)
                        except ValueError as exc:
                            raise ValueError(
                                f"Invalid Dupsifter count at {stat}:{line_no}: {key}={value!r}"
                            ) from exc
    row_key = demux_sample
    output_row = {
        "Protocol": pipeline_mode,
        "Route": route_id,
        "PlateID": plate_id,
        "Cell_Order": cell_order,
        "DNA_Barcode": dna_bc,
    }
    output_row.update(metrics)
    obj = {
        "id": "dupsifter",
        "section_name": "Dupsifter Stats",
        "description": f"Duplicate marking statistics for {demux_sample}; raw source={raw_sample}",
        "plot_type": "table",
        "pconfig": {"id": "dupsifter_table", "title": "Dupsifter Statistics", "col1_header": "Sample Name"},
        "data": {row_key: output_row},
    }
    atomic_write_json(out_json, obj)


def _parse_cph_methylation_pct(path: Path) -> float:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0].strip() != "BISCUITqc CpH Retention by Read Position Table":
        raise ValueError(f"Malformed BISCUIT CpH retention title: {path}")
    if len(lines) < 2:
        raise ValueError(f"BISCUIT CpH retention table lacks a header: {path}")
    reader = csv.DictReader(lines[1:], delimiter="\t")
    required = {"ReadInPair", "Position", "Conversion/Retention", "Count"}
    if not required.issubset(set(reader.fieldnames or [])):
        raise ValueError(f"Malformed BISCUIT CpH retention header: {path}")
    seen: set[tuple[str, int, str]] = set()
    retained = total = 0
    for row in reader:
        read = str(row["ReadInPair"]).strip()
        state = str(row["Conversion/Retention"]).strip()
        if read not in {"1", "2"} or state not in {"C", "R"}:
            raise ValueError(f"Malformed BISCUIT CpH retention row: {path}: {row}")
        try:
            position = int(str(row["Position"]).strip())
            count = int(str(row["Count"]).strip())
        except ValueError as exc:
            raise ValueError(f"Invalid BISCUIT CpH count/position: {path}: {row}") from exc
        if position < 0 or count < 0:
            raise ValueError(f"Negative BISCUIT CpH count/position: {path}: {row}")
        identity = (read, position, state)
        if identity in seen:
            raise ValueError(f"Duplicate BISCUIT CpH state: {path}: {row}")
        seen.add(identity)
        total += count
        if state == "R":
            retained += count
    return 100.0 * retained / total if total else math.nan


def _read_dupsifter_metrics(
    path: Path, sample: str, trimmed_pairs: int
) -> tuple[dict[str, int], dict[str, int]]:
    payload = _read_json_if_present(path)
    data = payload.get("data")
    if not isinstance(data, dict) or set(data) != {sample}:
        raise ValueError(f"Dupsifter QC must contain exactly sample {sample!r}: {path}")
    row = data[sample]
    if not isinstance(row, dict):
        raise ValueError(f"Dupsifter QC sample row is invalid: {path}")
    native: dict[str, int] = {}
    for field in DUPSIFTER_MQC_FIELDS:
        value = row.get(field)
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f"Dupsifter QC has invalid {field}: {path}")
        native[field.lower()] = value

    both = native["both_mapped"]
    one = native["one_mapped_forward"] + native["one_mapped_reverse"]
    unmapped = native["no_mapped"]
    no_primary = native["no_primary"]


    accepted = both + one
    rejected = unmapped + no_primary
    duplicate = native["dup_both"] + native["dup_forward"] + native["dup_reverse"]
    if accepted + rejected != trimmed_pairs:
        raise ValueError(
            "Dupsifter pair categories disagree with Cutadapt output: "
            f"accepted({accepted}) + rejected({rejected}) != trimmed({trimmed_pairs})"
        )
    if duplicate > accepted:
        raise ValueError(
            f"Dupsifter duplicate pairs exceed accepted mapped pairs: {path}"
        )
    return native, {
        "both_primary_mapped_pairs": both,
        "one_primary_mapped_pairs": one,
        "backend_accepted_pairs": accepted,
        "backend_rejected_pairs": rejected,
        "unmapped_pairs": unmapped,
        "no_primary_pairs": no_primary,
        "duplicate_pairs": duplicate,
        "postdedup_pairs": accepted - duplicate,
    }


def _read_biscuit_native_mapping(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing BISCUIT MAPQ table: {path}")
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 3 or lines[:2] != ["BISCUITqc Mapping Quality Table", "MapQ\tCount"]:
        raise ValueError(f"Malformed BISCUIT MAPQ table header or missing counts: {path}")
    rows = [line for line in lines[2:] if line.strip()]
    if not rows:
        raise ValueError(f"BISCUIT MAPQ table lacks count rows: {path}")
    optimal = suboptimal = unaligned = 0
    seen_mapq: set[str | int] = set()
    for line in rows:
        fields = line.split()
        if len(fields) != 2:
            raise ValueError(f"Malformed BISCUIT MAPQ row in {path}: {line!r}")
        label, raw_count = fields[0], fields[1].replace(",", "")
        try:
            count = int(raw_count)
        except ValueError as exc:
            raise ValueError(f"Invalid BISCUIT MAPQ count in {path}: {line!r}") from exc
        if count < 0:
            raise ValueError(f"Negative BISCUIT MAPQ count in {path}: {line!r}")
        label_key: str | int = label
        if label != "unmapped":
            try:
                label_key = int(label)
            except ValueError as exc:
                raise ValueError(f"Invalid BISCUIT MAPQ label in {path}: {label!r}") from exc
            if not 0 <= label_key <= 255:
                raise ValueError(f"Invalid BISCUIT MAPQ range in {path}: {label!r}")
        if label_key in seen_mapq:
            raise ValueError(f"Duplicate BISCUIT MAPQ category in {path}: {label!r}")
        seen_mapq.add(label_key)
        if label == "unmapped":
            unaligned += count
        else:
            if label_key >= 40:
                optimal += count
            else:
                suboptimal += count
    total = optimal + suboptimal + unaligned
    return {
        "optimally_aligned_reads": optimal,
        "suboptimally_aligned_reads": suboptimal,
        "unaligned_native_reads": unaligned,
        "native_mapping_policy": "biscuit_primary_mapq_ge_40",
        "native_mapping_unit": "individual_primary_reads_postdedup",
        "native_mapping_pct": 100.0 * optimal / total if total else 0.0,
    }


# 科学处理：Bismark 与 Droplet UMI

def parse_alignment_report(path: Path) -> dict[str, int | float]:
    text = path.read_text(encoding="utf-8", errors="strict")
    total = _required_count(
        text,
        "sequence_pairs_total",
        [rf"^Sequence pairs analysed in total:\s*{_COUNT}\s*$"],
    )
    unique = _required_count(
        text,
        "unique_best_pairs",
        [rf"^Number of paired-end alignments with a unique best hit:\s*{_COUNT}\s*$"],
    )
    no_alignment = _required_count(
        text,
        "pairs_without_alignment",
        [rf"^Sequence pairs with no alignments under any condition:\s*{_COUNT}\s*$"],
    )
    non_unique = _required_count(
        text,
        "pairs_not_unique",
        [rf"^Sequence pairs did not map uniquely:\s*{_COUNT}\s*$"],
    )
    discarded = _optional_count(
        text,
        "pairs_discarded_no_genomic_sequence",
        [
            rf"^Sequence pairs which were discarded because genomic sequence could not be extracted:\s*{_COUNT}\s*$"
        ],
    )
    mapping_efficiency = _required_float(
        text, "mapping_efficiency", r"^Mapping efficiency:\s*([0-9]+(?:\.[0-9]+)?)%\s*$"
    )
    if unique + no_alignment + non_unique != total:
        raise ValueError(
            "Bismark alignment counts are inconsistent: "
            f"unique({unique}) + unaligned({no_alignment}) + "
            f"non_unique({non_unique}) != total({total})"
        )
    if discarded > unique:
        raise ValueError(
            "Bismark discarded-pair count exceeds uniquely aligned pairs: "
            f"discarded={discarded}, unique={unique}"
        )
    expected_efficiency = 0.0 if total == 0 else unique * 100.0 / total
    if not math.isclose(mapping_efficiency, expected_efficiency, abs_tol=0.11):
        raise ValueError(
            "Bismark mapping efficiency disagrees with the exact report counts: "
            f"reported={mapping_efficiency}, expected={expected_efficiency:.6f}"
        )
    return {
        "total": total,
        "unique": unique,
        "no_alignment": no_alignment,
        "non_unique": non_unique,
        "discarded": discarded,
        "mapping_efficiency": mapping_efficiency,
    }


def parse_dedup_report(path: Path) -> dict[str, int]:
    text = path.read_text(encoding="utf-8", errors="strict")
    total = _required_count(
        text,
        "dedup_pairs_analysed",
        [rf"^Total number of alignments analysed in .+?:\s*{_COUNT}(?:\s|$)"],
    )
    removed = _required_count(
        text,
        "duplicate_pairs_removed",
        [rf"^Total number duplicated alignments removed:\s*{_COUNT}(?:\s|$)"],
    )
    leftover = _required_count(
        text,
        "deduplicated_pairs_left",
        [rf"^Total count of deduplicated leftover sequences:\s*{_COUNT}(?:\s|$)"],
    )
    if removed > total or leftover != total - removed:
        raise ValueError(
            "Bismark deduplication counts are inconsistent: "
            f"total={total}, removed={removed}, leftover={leftover}"
        )
    return {"total": total, "removed": removed, "leftover": leftover}


def _read_filter_summary(path: Path, sample: str) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("Bismark non-conversion summary must be an object")
    if set(payload) != FILTER_SUMMARY_FIELDS:
        raise ValueError(
            "Bismark non-conversion summary schema mismatch: "
            f"expected={sorted(FILTER_SUMMARY_FIELDS)}, observed={sorted(payload)}"
        )
    if (
        type(payload["schema_version"]) is not int
        or payload["schema_version"] != FILTER_SUMMARY_SCHEMA_VERSION
        or payload["backend"] != "bismark"
        or payload["metric_unit"] != "read_pairs"
    ):
        raise ValueError("Invalid Bismark non-conversion summary identity")
    if payload["sample"] != sample:
        raise ValueError(
            f"Bismark non-conversion summary sample mismatch: {payload['sample']!r}"
        )
    for field in (
        "candidate_read_pairs",
        "flagged_read_pairs",
        "excluded_read_pairs",
        "final_retained_pairs",
        "minimum_count",
    ):
        if isinstance(payload[field], bool) or not isinstance(payload[field], int):
            raise ValueError(f"Bismark non-conversion summary {field} must be an integer")
    for field in ("high_cph_fraction", "threshold"):
        if isinstance(payload[field], bool) or not isinstance(
            payload[field], (int, float)
        ):
            raise ValueError(f"Bismark non-conversion summary {field} must be numeric")
    candidate = payload["candidate_read_pairs"]
    high_cph = payload["flagged_read_pairs"]
    excluded = payload["excluded_read_pairs"]
    retained = payload["final_retained_pairs"]
    fraction = float(payload["high_cph_fraction"])
    threshold = float(payload["threshold"])
    minimum_count = payload["minimum_count"]
    filter_applied = payload["filter_applied"]
    if candidate < 0 or excluded < 0 or not 0 <= high_cph <= candidate:
        raise ValueError("Invalid Bismark non-conversion summary counts")
    if retained != candidate - high_cph + excluded:
        raise ValueError(
            "Bismark non-conversion retained pairs disagree with the filter funnel: "
            f"retained={retained}, expected={candidate - high_cph + excluded} "
            f"(candidate={candidate}, flagged={high_cph}, excluded={excluded})"
        )
    expected_fraction = 0.0 if candidate == 0 else high_cph / candidate
    if not math.isclose(fraction, expected_fraction, abs_tol=1e-12):
        raise ValueError("Bismark non-conversion summary fraction is inconsistent")
    if not 0 <= threshold <= 1:
        raise ValueError("Bismark non-conversion threshold must be within [0, 1]")
    if minimum_count < 1:
        raise ValueError("Bismark non-conversion minimum_count must be >= 1")
    if not isinstance(filter_applied, bool):
        raise ValueError("Bismark non-conversion filter_applied must be boolean")
    if payload["protocol"] not in {"cabernet", "srd", "droplet"}:
        raise ValueError("Bismark non-conversion protocol is invalid")
    if payload["high_cph_role"] not in {
        "cdna_contamination",
        "residual_high_cph",
    }:
        raise ValueError("Bismark non-conversion high_cph_role is invalid")
    expected_role = {"cabernet": "cdna_contamination", "droplet": "cdna_contamination", "srd": "residual_high_cph"}
    if payload["high_cph_role"] != expected_role[payload["protocol"]]:
        raise ValueError("Bismark non-conversion protocol/role contract is invalid")
    if not filter_applied:
        raise ValueError("Bismark non-conversion summary must record the fixed filter action")
    return payload


DROPLET_DEDUP_POLICY = "per_cytosine_umi_consensus_hamming1"


DROPLET_MIN_MAPQ = 10


DROPLET_MIN_BASEQ = 20


def droplet_umi_consensus(groups: Mapping[str, Sequence[int]]) -> tuple[int, int, int, int]:
    """逐位点按原始丰度合并 Hamming-1 UMI；多数碱基优先，平票比较平均质量，仍平票则不调用。"""
    counts = {umi: int(v[0]) + int(v[1]) for umi, v in groups.items()}
    remaining = {umi: list(v) for umi, v in groups.items()}
    merged = ambiguous = methyl = unmethyl = 0
    for umi in sorted(counts, key=lambda key: (-counts[key], key)):
        if umi not in remaining:
            continue
        value = remaining[umi]
        if counts[umi] > 1 and len(remaining) > 1:
            if len(remaining) < 48:
                neighbors = [key for key in remaining if counts[key] < counts[umi]
                             and sum(a != b for a, b in zip(umi, key)) == 1]
                neighbors.sort(key=lambda key: next((i, "ACGTN".index(b)) for i, (a, b) in enumerate(zip(umi, key)) if a != b))
            else:
                neighbors = [umi[:i] + alternate + umi[i + 1:] for i, base in enumerate(umi)
                             for alternate in "ACGTN" if alternate != base]
            for neighbor in neighbors:
                if neighbor in remaining and counts[neighbor] < counts[umi]:
                    other = remaining.pop(neighbor)
                    value = [a + b for a, b in zip(value, other)]
                    merged += 1
        remaining.pop(umi)
        m, u, qm, qu = value
        if m > u or (m == u and qm > qu):
            methyl += 1
        elif u > m or (m == u and qu > qm):
            unmethyl += 1
        else:
            ambiguous += 1
    return methyl, unmethyl, merged, ambiguous


def call_droplet_cpg(bam: str | Path, reference: str | Path, coverage: str | Path,
                     qc: str | Path, mbias: str | Path, *, sample: str,
                     scratch: str | Path, samtools: str, threads: int = 1) -> None:
    """保留过滤 BAM 的全部 pairs，按参考胞嘧啶/UMI 生成分子共识 coverage、QC 及 read-level M-bias。"""
    import pysam

    if SAMPLE_ID_RE.fullmatch(sample) is None:
        raise ValueError("Invalid Droplet sample identity")
    root = Path(tempfile.mkdtemp(prefix="droplet_consensus_", dir=scratch))
    proc = None
    try:
        converted = root / "conversion.bam"
        ordered = root / "coordinate.bam"
        pairs = umi_n = 0
        seen = set()
        checked_contigs = set()
        with pysam.FastaFile(str(reference)) as fasta:
            reference_lengths = dict(zip(fasta.references, fasta.lengths))
        with pysam.AlignmentFile(str(bam), "rb") as source, pysam.AlignmentFile(str(converted), "wb", template=source) as target:
            iterator = iter(source)
            for first in iterator:
                second = next(iterator, None)
                if second is None or first.query_name != second.query_name or first.query_name in seen:
                    raise ValueError("Droplet calling requires unique, adjacent complete read pairs")
                seen.add(first.query_name)
                if (first.flag, second.flag) not in {(99, 147), (163, 83), (147, 99), (83, 163)}:
                    raise ValueError("Unexpected Bismark paired flags")
                if (not first.has_tag("XG") or not second.has_tag("XG")
                        or first.get_tag("XG") != second.get_tag("XG")):
                    raise ValueError("Droplet mates disagree on Bismark genome conversion")
                if (any(not record.has_tag("XR") for record in (first, second))
                        or {first.get_tag("XR"), second.get_tag("XR")} != {"CT", "GA"}):
                    raise ValueError("Droplet mates require complementary Bismark read conversions")
                if (first.reference_id != second.reference_id
                        or any(record.next_reference_id != mate.reference_id
                               or record.next_reference_start != mate.reference_start
                               for record, mate in ((first, second), (second, first)))):
                    raise ValueError("Droplet mate coordinates disagree with the actual pair")
                umi = first.query_name.rsplit(":", 1)[-1]
                if not re.fullmatch(r"[ACGTN]{12}", umi):
                    raise ValueError("Droplet QNAME must end with a 12 bp UMI")
                for record in (first, second):
                    mate = 1 if record.is_read1 else 2
                    contig = record.reference_name
                    if contig not in checked_contigs:
                        if contig not in reference_lengths or source.get_reference_length(contig) != reference_lengths[contig]:
                            raise ValueError(f"Droplet BAM/reference contig mismatch: {contig}")
                        checked_contigs.add(contig)
                    if (record.reference_end is None
                            or not 0 <= record.reference_start < record.reference_end <= reference_lengths[contig]):
                        raise ValueError(f"Droplet alignment exceeds reference bounds: {contig}")
                    if any(not record.has_tag(tag) for tag in ("XM", "XR", "XG")):
                        raise ValueError("Bismark methylation tags missing")
                    if (record.query_sequence is None or record.query_qualities is None
                            or len(record.get_tag("XM")) != record.query_length):
                        raise ValueError("Droplet calling requires sequence, base qualities and matching XM length")
                    if record.has_tag("UR") and record.get_tag("UR") != umi:
                        raise ValueError("Droplet UR tag disagrees with QNAME")
                    conversion = record.get_tag("XG")
                    if conversion not in {"CT", "GA"}:
                        raise ValueError("Invalid Bismark genome conversion")
                    record.set_tag("UR", umi, value_type="Z")
                    record.set_tag("YI", f"{mate}{'R' if record.is_reverse else 'F'}", value_type="Z")
                    record.is_reverse = conversion == "GA"
                    record.mate_is_reverse = conversion == "GA"
                    target.write(record)
                pairs += 1
                umi_n += int("N" in umi)
        del seen
        subprocess.run([samtools, "sort", "-@", str(max(0, min(threads, 4) - 1)),
                        "-m", "512M", "-T", str(root / "sort"), "-o", str(ordered), str(converted)], check=True)
        converted.unlink()
        contexts = {key: {"methyl": 0, "unmethyl": 0, "read_observations": 0,
                          "umi_groups": 0, "merged_umis": 0, "ambiguous_umis": 0,
                          "strand_sites": 0} for key in ("CpG", "CHG", "CHH")}
        cycles: dict[tuple[str, str, int], list[int]] = {}
        command = [samtools, "mpileup", "-Q", str(DROPLET_MIN_BASEQ), "-q", str(DROPLET_MIN_MAPQ),
                   "-B", "-d", "0", "-f", str(reference), "--output-BP", "--output-extra", "RLEN,UR,YI",
                   "--no-output-ins-mods", "--no-output-ins", "--no-output-ins",
                   "--no-output-del", "--no-output-del", "--no-output-ends", str(ordered)]
        with pysam.FastaFile(str(reference)) as fa, gzip.open(root / "coverage.gz", "wt") as handle:
            proc = subprocess.Popen(command, stdout=subprocess.PIPE, text=True, bufsize=1024 * 1024)
            chrom = None
            sequence = ""
            for line in proc.stdout:
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 10:
                    raise ValueError("Unexpected samtools mpileup columns")
                if int(fields[3]) == 0:
                    continue
                base = fields[2].upper()
                if base not in {"C", "G"}:
                    continue
                if chrom != fields[0]:
                    chrom = fields[0]
                    sequence = fa.fetch(chrom).upper()
                pos = int(fields[1]) - 1
                if not 0 <= pos < len(sequence) or sequence[pos] != base:
                    raise ValueError("Droplet pileup coordinate disagrees with reference")
                context_seq = sequence[pos:pos + 3] if base == "C" else sequence[max(0, pos - 2):pos + 1][::-1].translate(str.maketrans("ACGT", "TGCA"))
                if context_seq.startswith("CG"):
                    context = "CpG"
                elif len(context_seq) == 3 and context_seq[1] in "ACT" and context_seq[2] in "ACGT":
                    context = "CHG" if context_seq[2] == "G" else "CHH"
                else:
                    continue
                bases, qualities = fields[4:6]
                positions, lengths, umis, mates = (value.split(",") for value in fields[6:10])
                if len({len(bases), len(qualities), len(positions), len(lengths), len(umis), len(mates), int(fields[3])}) != 1:
                    raise ValueError(f"Unpaired mpileup observations at {chrom}:{pos + 1}")
                methyl_base, unmethyl_base = (".", "T") if base == "C" else (",", "a")
                groups: dict[str, list[int]] = {}
                for observed, quality, position, length, umi, mate in zip(bases, qualities, positions, lengths, umis, mates):
                    if observed not in (methyl_base, unmethyl_base):
                        continue
                    if not re.fullmatch(r"[ACGTN]{12}", umi) or mate not in {"1F", "1R", "2F", "2R"}:
                        raise ValueError("Invalid Droplet pileup molecule identity")
                    state = int(observed == unmethyl_base)
                    values = groups.setdefault(umi, [0, 0, 0, 0])
                    values[state] += 1
                    values[state + 2] += ord(quality) - 33
                    cycle = int(position) if mate[1] == "F" else int(length) - int(position) + 1
                    if not 1 <= cycle <= int(length):
                        raise ValueError("Invalid Droplet pileup read position")
                    cycles.setdefault((context, mate[0], cycle), [0, 0])[state] += 1
                if not groups:
                    continue
                m, u, merged, ambiguous = droplet_umi_consensus(groups)
                stats = contexts[context]
                stats["read_observations"] += sum(v[0] + v[1] for v in groups.values())
                stats["umi_groups"] += len(groups)
                stats["merged_umis"] += merged
                stats["ambiguous_umis"] += ambiguous
                stats["methyl"] += m
                stats["unmethyl"] += u
                stats["strand_sites"] += int(m + u > 0)
                if context == "CpG" and m + u:
                    handle.write(f"{chrom}\t{pos + 1}\t{pos + 1}\t{100 * m / (m + u):.8f}\t{m}\t{u}\n")
            proc.stdout.close()
            if proc.wait() != 0:
                raise RuntimeError("Droplet samtools mpileup failed")
            proc = None
        with (root / "mbias.txt").open("w") as handle:
            handle.write("Alopex Droplet read-level M-bias; MAPQ>=10, baseQ>=20, before UMI consensus\n\n")
            for context in contexts:
                for mate in ("1", "2"):
                    handle.write(f"{context} context (R{mate})\n====================\nposition\tcount_methylated\tcount_unmethylated\t% methylation\tcoverage\n")
                    for ctx, read, cycle in sorted(cycles):
                        if ctx == context and read == mate:
                            m, u = cycles[(ctx, read, cycle)]
                            handle.write(f"{cycle}\t{m}\t{u}\t{100 * m / (m + u):.8f}\t{m + u}\n")
                    handle.write("\n")
        cpg = contexts["CpG"]
        payload = {"schema_version": 1, "sample": sample, "dedup_policy": DROPLET_DEDUP_POLICY,
                   "min_mapq": DROPLET_MIN_MAPQ, "min_baseq": DROPLET_MIN_BASEQ,
                   "metric_unit": "cytosine_umi_consensus", "processed_pairs": pairs,
                   "umi_with_n_pairs": umi_n, "contexts": contexts,
                   "id": "alopex_droplet_umi_consensus", "section_name": "Droplet UMI consensus",
                   "description": "逐胞嘧啶 UMI 共识；read-level 观测及分子计数不等同 read-pair 去重率。M-bias 为质量过滤后、UMI 共识前的 read-level 统计。",
                   "plot_type": "table", "data": {sample: {"CpG_read_observations": cpg["read_observations"],
                       "CpG_UMI_observations": cpg["methyl"] + cpg["unmethyl"],
                       "CpG_ambiguous_UMIs": cpg["ambiguous_umis"], "CpG_strand_sites": cpg["strand_sites"]}}}
        (root / "qc.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
        read_droplet_umi_qc(root / "qc.json", sample)
        for name, target in (("coverage.gz", coverage), ("mbias.txt", mbias), ("qc.json", qc)):
            os.replace(root / name, target)
    finally:
        if proc is not None:
            _close_native_stream(proc)
        shutil.rmtree(root)


def read_droplet_umi_qc(path: str | Path, sample: str) -> dict[str, Any]:
    """校验逐胞嘧啶 UMI 共识计数、单位与实际参数，禁止作为 read-pair 去重率解释。"""
    data = _read_json(Path(path), "Droplet UMI consensus QC")
    for key, value in {"schema_version": 1, "sample": sample, "dedup_policy": DROPLET_DEDUP_POLICY,
                       "min_mapq": DROPLET_MIN_MAPQ, "min_baseq": DROPLET_MIN_BASEQ,
                       "metric_unit": "cytosine_umi_consensus"}.items():
        if type(data.get(key)) is not type(value) or data.get(key) != value:
            raise ValueError(f"Invalid Droplet UMI QC {key}: {path}")
    for key in ("processed_pairs", "umi_with_n_pairs"):
        data[key] = _integer(data.get(key), key, Path(path))
    if not 0 <= data["umi_with_n_pairs"] <= data["processed_pairs"]:
        raise ValueError("Invalid Droplet UMI pair counts")
    for context in ("CpG", "CHG", "CHH"):
        values = data.get("contexts", {}).get(context, {})
        for key in ("methyl", "unmethyl", "read_observations", "umi_groups", "merged_umis", "ambiguous_umis", "strand_sites"):
            values[key] = _integer(values.get(key), key, Path(path))
        called = values["methyl"] + values["unmethyl"]
        if ((data["processed_pairs"] == 0 and values["read_observations"] != 0)
                or bool(called) != bool(values["strand_sites"])):
            raise ValueError("Droplet UMI observations disagree with input pairs or called sites")
        if not (values["strand_sites"] <= called <= values["umi_groups"] <= values["read_observations"]):
            raise ValueError("Droplet UMI observation accounting is inconsistent")
        if values["umi_groups"] != called + values["merged_umis"] + values["ambiguous_umis"]:
            raise ValueError("Droplet UMI consensus accounting does not close")
    return data


def build_metrics(
    *,
    sample: str,
    version: str,
    alignment_report: Path,
    dedup_report: Path,
    filter_summary: Mapping[str, Any],
    protocol: str = "cabernet",
) -> dict[str, Any]:
    """reconcile 各 pair-level 来源：alignment/dedup 报告与恒定执行的 non-conversion
    已统一校验的过滤摘要（含 processed BAM 的 retained 实测计数）必须守恒，矛盾即抛错。"""
    version_lines = [line.strip() for line in version.splitlines() if line.strip()]
    if not sample or not version_lines:
        raise ValueError("sample and Bismark version are required")
    version_label = next(
        (
            line
            for line in version_lines
            if re.search(
                r"\bbismark\b.*?(?<![A-Za-z0-9])v?3(?:\.[0-9]+)+\b",
                line,
                flags=re.IGNORECASE,
            )
        ),
        "",
    )
    if not version_label:
        raise ValueError(f"Expected a Bismark 3.x version string, observed: {version!r}")
    alignment = parse_alignment_report(alignment_report)
    unique_pairs = int(alignment["unique"])
    discarded_pairs = int(alignment["discarded"])
    accepted_pairs = unique_pairs - discarded_pairs

    umi_qc = read_droplet_umi_qc(dedup_report, sample) if protocol == "droplet" else None
    dedup = {"total": accepted_pairs, "removed": 0, "leftover": accepted_pairs} if umi_qc is not None else parse_dedup_report(dedup_report)
    if dedup["total"] != accepted_pairs:
        raise ValueError(
            "Bismark deduplication input disagrees with usable uniquely aligned pairs: "
            f"analysed={dedup['total']}, expected={accepted_pairs} "
            f"(unique={unique_pairs}, discarded={discarded_pairs})"
        )
    duplicate_pairs = dedup["removed"]
    postdedup_pairs = dedup["leftover"]

    filtered = filter_summary
    candidate_pairs = int(filtered["candidate_read_pairs"])
    excluded_pairs = int(filtered["excluded_read_pairs"])
    flagged_pairs = int(filtered["flagged_read_pairs"])
    if candidate_pairs + excluded_pairs != postdedup_pairs:
        raise ValueError(
            "Bismark filtering summary disagrees with post-dedup pairs: "
            f"candidate({candidate_pairs}) + excluded({excluded_pairs}) "
            f"!= postdedup({postdedup_pairs})"
        )
    removed_pairs = flagged_pairs

    final_retained_pairs = int(filtered["final_retained_pairs"])
    if final_retained_pairs != postdedup_pairs - removed_pairs:
        raise ValueError(
            "Bismark processed-BAM retained pairs disagree with the pair funnel: "
            f"retained={final_retained_pairs}, "
            f"expected={postdedup_pairs - removed_pairs}"
        )

    if umi_qc is not None and umi_qc["processed_pairs"] != final_retained_pairs:
        raise ValueError("Droplet caller input disagrees with non-conversion retained pairs")
    return {
        "umi_qc": umi_qc,
        "total_pairs": int(alignment["total"]),
        "accepted_pairs": accepted_pairs,
        "unmapped_pairs": int(alignment["no_alignment"]) + discarded_pairs,
        "ambiguous_pairs": int(alignment["non_unique"]),
        "duplicate_pairs": duplicate_pairs,
        "postdedup_pairs": postdedup_pairs,
        "final_retained_pairs": final_retained_pairs,
    }


def write_filter_summary(
    *,
    sample: str,
    candidate_pairs: int,
    removed_pairs: int,
    excluded_pairs: int,
    retained_pairs: int,
    percentage: float,
    minimum_count: int,
    protocol: str,
    high_cph_role: str,
    summary: Path,
    multiqc: Path,
) -> None:
    """按规则内实测计数写出 schema v3 filter summary 与 MultiQC table JSON。

    retained 必须等于 candidate - removed + excluded（processed BAM 的实测
    pair 数）；任何计数矛盾都在写出前失败。
    """
    if isinstance(minimum_count, bool) or not isinstance(minimum_count, int) or minimum_count < 1:
        raise SystemExit("filter summary minimum_count must be an integer >= 1")
    if (
        isinstance(percentage, bool)
        or not isinstance(percentage, (int, float))
        or not 0 <= percentage <= 100
    ):
        raise SystemExit("filter summary percentage must be within [0, 100]")
    if min(candidate_pairs, removed_pairs, excluded_pairs, retained_pairs) < 0:
        raise SystemExit("filter summary counts must be non-negative integers")
    if removed_pairs > candidate_pairs:
        raise SystemExit("flagged pairs exceed candidate pairs")
    if retained_pairs != candidate_pairs - removed_pairs + excluded_pairs:
        raise SystemExit(
            "retained pairs disagree with the filter funnel: "
            f"retained={retained_pairs}, "
            f"expected={candidate_pairs - removed_pairs + excluded_pairs}"
        )
    payload = {
        "schema_version": FILTER_SUMMARY_SCHEMA_VERSION,
        "backend": "bismark",
        "sample": sample,
        "metric_unit": "read_pairs",
        "candidate_read_pairs": candidate_pairs,
        "flagged_read_pairs": removed_pairs,
        "high_cph_fraction": (
            removed_pairs / candidate_pairs if candidate_pairs else 0.0
        ),
        "excluded_read_pairs": excluded_pairs,
        "final_retained_pairs": retained_pairs,
        "threshold": percentage / 100.0,
        "minimum_count": minimum_count,
        "filter_applied": True,
        "protocol": protocol,
        "high_cph_role": high_cph_role,
    }
    summary.write_text(
        json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    multiqc.write_text(
        json.dumps(
            _high_cph_multiqc_payload(sample, payload, protocol, high_cph_role),
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n",
        encoding="utf-8",
    )


def run_bismark_pair_split(args) -> int:
    """把 BAM 路由到 bio/excluded，标准输出仅返回两路 pair 计数供当前规则使用。"""
    excluded_contigs = _parse_excluded_contigs_json(args.excluded_contigs_json)
    args.bio.parent.mkdir(parents=True, exist_ok=True)
    args.excluded.parent.mkdir(parents=True, exist_ok=True)
    print(*split_pairs(args.input_bam, args.bio, args.excluded, excluded_contigs))
    return 0


def _parse_bismark_non_cpg_methylation_pct(path: Path) -> float:
    text = path.read_text(encoding="utf-8", errors="strict")
    counts: dict[tuple[str, str], int] = {}
    for state in ("Methylated", "Unmethylated"):
        for context in ("CHG", "CHH"):
            if state == "Methylated":
                patterns = (
                    rf"^(?:Total )?methylated C(?:'|’)?s in {context} context:\s*([0-9][0-9,]*)\s*$",
                    rf"^(?:Total )?methylated cytosines in {context} context:\s*([0-9][0-9,]*)\s*$",
                )
            else:
                patterns = (
                    rf"^(?:Total )?unmethylated C(?:'|’)?s in {context} context:\s*([0-9][0-9,]*)\s*$",
                    rf"^(?:Total )?unmethylated cytosines in {context} context:\s*([0-9][0-9,]*)\s*$",
                    rf"^Total C to T conversions in {context} context:\s*([0-9][0-9,]*)\s*$",
                )
            observed = [
                int(match.group(1).replace(",", ""))
                for pattern in patterns
                for match in re.finditer(pattern, text, re.IGNORECASE | re.MULTILINE)
            ]
            if not observed:
                raise ValueError(
                    f"Bismark extraction report lacks {state} {context} count: {path}"
                )
            if len(set(observed)) != 1:
                raise ValueError(
                    f"Bismark extraction report has conflicting {state} {context} counts: {path}"
                )
            counts[(state, context)] = observed[0]
    methylated = sum(counts[("Methylated", context)] for context in ("CHG", "CHH"))
    total = methylated + sum(
        counts[("Unmethylated", context)] for context in ("CHG", "CHH")
    )
    return 100.0 * methylated / total if total else math.nan


def _bismark_version_text(bismark: str | None) -> str:
    if not bismark:
        return ""
    version = subprocess.run(
        [bismark, "--version"], capture_output=True, text=True, check=True
    )
    return (version.stdout + version.stderr).strip()


# 科学处理：TAPS

TAPS_BWA_INSERT_UPPER_FLOOR = 1000


TAPS_BWA_INSERT_PROBE_PAIRS = 10000


TAPS_BWA_BATCH_BASES = 150000000


def taps_bwa_insert_model(r1: str, r2: str, bwa: str, prefix: str, scratch: str, threads: int) -> str:
    """用 BWA 原生 FR 推断均值与标准差，将上限至少放宽至 1000 bp；证据不足时沿用自动推断。"""
    import pysam

    sampled = 0
    root = Path(scratch)
    paths = [root / "insert_probe_R1.fastq", root / "insert_probe_R2.fastq"]
    log = root / "insert_probe.log"
    try:
        with pysam.FastxFile(r1) as left, pysam.FastxFile(r2) as right, paths[0].open("w") as one, paths[1].open("w") as two:
            for a, b in zip(left, right, strict=True):
                if a.name.removesuffix("/1") != b.name.removesuffix("/2"):
                    raise ValueError("TAPS insert probe requires matching paired FASTQ names")
                one.write(str(a) + "\n")
                two.write(str(b) + "\n")
                sampled += 1
                if sampled == TAPS_BWA_INSERT_PROBE_PAIRS:
                    break
        model = ""
        if sampled:
            with log.open("w") as handle:
                subprocess.run([bwa, "mem", "-Y", "-t", str(max(1, min(4, threads))), prefix, *map(str, paths)],
                               stdout=subprocess.DEVNULL, stderr=handle, check=True)
            native = log.read_text()
            block = re.search(r"analyzing insert size distribution for orientation FR\.\.\.(.*?)(?:analyzing insert size distribution for orientation RF|skip orientation RF|skip orientation RR)", native, re.S)
            if block:
                location = re.search(r"\] mean and std.dev: \(([^,]+), ([^)]+)\)", block[1])
                bounds = re.search(r"proper pairs: \((\d+), (\d+)\)", block[1])
                if location and bounds:
                    mean, deviation = map(float, location.groups())
                    lower, upper = map(int, bounds.groups())
                    if math.isfinite(mean) and math.isfinite(deviation) and mean > 0 and deviation > 0 and lower <= upper:
                        model = f"{location[1]},{location[2]},{max(TAPS_BWA_INSERT_UPPER_FLOOR, upper)},{lower}"
        print(f"TAPS insert probe: pairs={sampled} FR model={model or 'native automatic'}", file=sys.stderr)
        return model
    finally:
        if log.is_file():
            print(log.read_text(), file=sys.stderr, end="")
        for path in [*paths, log]:
            path.unlink(missing_ok=True)


TAPS_MIN_MAPQ = 20


TAPS_MIN_BASEQ = 30


TAPS_INCLUDE_FLAGS = 3


TAPS_EXCLUDE_FLAGS = 3852


def write_taps_alignment_qc(bam: str, samtools: str, scratch: str, sample: str, output: str) -> None:
    """按 QNAME 核算完整 BAM 与 Rastair 可用 pair；仅任一 mate 达到调用条件的 pair 进入漏斗。"""
    import itertools
    import pysam

    counts = dict.fromkeys(("bam_records", "bam_primary_records", "bam_duplicate_records",
                           "bam_unmapped_records", "bam_supplementary_records", "bam_secondary_records",
                           "trimmed_pairs", "backend_accepted_pairs", "backend_rejected_pairs",
                           "unmapped_pairs", "discordant_pairs", "quality_rejected_pairs",
                           "duplicate_pairs", "postdedup_pairs"), 0)
    with pysam.AlignmentFile(bam, "rb") as source:
        header = source.header.to_dict()
        if header.get("HD", {}).get("SO") != "coordinate" or not source.has_index():
            raise ValueError("TAPS BAM must be coordinate sorted and indexed")
        groups = header.get("RG", [])
        group_ids = {rg.get("ID") for rg in groups}
        if (not groups or None in group_ids or "" in group_ids or len(group_ids) != len(groups)
                or any(rg.get("SM") != sample for rg in groups)):
            raise ValueError("TAPS BAM requires unique read groups bound to the active sample")
    process = subprocess.Popen([samtools, "collate", "-O", "-u", "-T", str(Path(scratch) / "qc_collate"), bam], stdout=subprocess.PIPE)
    try:
        with pysam.AlignmentFile(process.stdout, "rb") as source:
            for name, group in itertools.groupby(source, key=lambda r: r.query_name):
                records = list(group)
                primary = [r for r in records if not r.is_secondary and not r.is_supplementary]
                if len(primary) != 2 or {r.is_read1 for r in primary} != {True, False} or any(
                    not r.is_paired or r.is_read1 == r.is_read2 for r in primary
                ):
                    raise ValueError(f"TAPS primary pair is incomplete or duplicated: {name}")
                for r in records:
                    if not r.has_tag("RG") or r.get_tag("RG") not in group_ids:
                        raise ValueError(f"TAPS record has an unbound read group: {name}")
                    counts["bam_records"] += 1
                    counts["bam_primary_records"] += int(not r.is_secondary and not r.is_supplementary)
                    for key, flag in (("duplicate", r.is_duplicate), ("unmapped", r.is_unmapped),
                                      ("supplementary", r.is_supplementary), ("secondary", r.is_secondary)):
                        counts[f"bam_{key}_records"] += int(flag)
                counts["trimmed_pairs"] += 1
                accepted = [r for r in primary if r.flag & TAPS_INCLUDE_FLAGS == TAPS_INCLUDE_FLAGS
                            and not r.flag & (TAPS_EXCLUDE_FLAGS & ~1024) and r.mapping_quality >= TAPS_MIN_MAPQ]
                if accepted:
                    counts["backend_accepted_pairs"] += 1
                    if any(not r.is_duplicate for r in accepted):
                        counts["postdedup_pairs"] += 1
                    else:
                        counts["duplicate_pairs"] += 1
                else:
                    counts["backend_rejected_pairs"] += 1
                    if all(r.is_unmapped for r in primary):
                        counts["unmapped_pairs"] += 1
                    elif not all(r.is_proper_pair and not r.is_unmapped and not r.mate_is_unmapped for r in primary):
                        counts["discordant_pairs"] += 1
                    else:
                        counts["quality_rejected_pairs"] += 1
        if process.wait() != 0:
            raise RuntimeError("samtools collate failed during TAPS pair accounting")
    finally:
        _close_native_stream(process)
    atomic_write_json(output, {"sample": sample, "backend": "rastair", "schema_version": 1,
                               "read_groups": header.get("RG", []), "programs": header.get("PG", []),
                               "reference_sequences": header["SQ"], **counts})


def _taps_pair_bases(reads: list, reference) -> tuple[str, dict[int, tuple[str, int]], str, int, str]:
    eligible = [r for r in reads if r.flag & TAPS_INCLUDE_FLAGS == TAPS_INCLUDE_FLAGS and (not r.flag & TAPS_EXCLUDE_FLAGS) and (r.mapping_quality >= TAPS_MIN_MAPQ)]
    if not eligible:
        return ('', {}, '', 0, '')
    chroms = {r.reference_name for r in eligible}
    if len(chroms) != 1:
        raise ValueError('TAPS proper pair spans different reference contigs')
    chrom = eligible[0].reference_name
    start = max(0, min((r.reference_start for r in eligible)) - 2)
    end = min(reference.get_reference_length(chrom), max((r.reference_end for r in eligible)) + 2)
    seq = reference.fetch(chrom, start, end).upper()
    bases = {}
    evidence = 0
    aligned_reads = []
    candidate = False
    for read in eligible:
        observed = read.query_sequence
        qualities = read.query_qualities
        if observed is None or qualities is None:
            raise ValueError('TAPS QC requires BAM sequence and base qualities')
        observed = observed.upper()
        cigar = read.cigartuples or []
        body = [(i, op, n) for i, (op, n) in enumerate(cigar) if op not in (4, 5)]
        linear = len(body) == 1 and body[0][1] in (0, 7, 8) and all((i in (0, len(cigar) - 1) for i, (op, _) in enumerate(cigar) if op in (4, 5)))
        if linear:
            qstart, qend = (read.query_alignment_start, read.query_alignment_end)
            rstart, rend = (read.reference_start, read.reference_end)
            positions = None
            if not candidate:
                candidate = any((q >= TAPS_MIN_BASEQ for q in qualities[qstart:qend]))
        else:
            positions = read.get_reference_positions(full_length=True)
            qstart = qend = rstart = rend = 0
            if not candidate:
                candidate = any((pos is not None and quality >= TAPS_MIN_BASEQ for pos, quality in zip(positions, qualities)))
        aligned_reads.append((observed, qualities, positions, qstart, qend, rstart, rend))
        for motif, flag in (('TG', 1), ('CA', 2)):
            if evidence & flag:
                continue
            qpos = observed.find(motif)
            while qpos >= 0:
                if linear:
                    contiguous = qstart <= qpos and qpos + 1 < qend
                    pos = rstart + qpos - qstart
                else:
                    pos = positions[qpos]
                    contiguous = pos is not None and positions[qpos + 1] == pos + 1
                if contiguous and qualities[qpos] >= TAPS_MIN_BASEQ and (qualities[qpos + 1] >= TAPS_MIN_BASEQ) and (seq[pos - start:pos - start + 2] == 'CG'):
                    evidence |= flag
                    break
                qpos = observed.find(motif, qpos + 1)
    orientation = 'OT' if evidence == 1 else 'OB' if evidence == 2 else ''
    if not candidate:
        return ('', {}, '', 0, '')
    if orientation:
        target = 'C' if orientation == 'OT' else 'G'
        context_positions = set()
        offset = seq.find(target)
        while offset >= 0:
            if orientation == 'OT':
                valid = offset + 2 < len(seq) and seq[offset + 1] in 'ACT' and (seq[offset + 2] in 'ACGT')
            else:
                valid = offset >= 2 and seq[offset - 1] in 'AGT' and (seq[offset - 2] in 'ACGT')
            if valid:
                context_positions.add(start + offset)
            offset = seq.find(target, offset + 1)
        for observed, qualities, positions, qstart, qend, rstart, rend in aligned_reads:
            if positions is None:
                selected = ((qstart + pos - rstart, pos) for pos in context_positions if rstart <= pos < rend)
            else:
                selected = ((qpos, pos) for qpos, pos in enumerate(positions) if pos in context_positions)
            for qpos, pos in selected:
                quality = qualities[qpos]
                if quality < TAPS_MIN_BASEQ:
                    continue
                base = observed[qpos]
                previous = bases.get(pos)
                if previous is None or quality > previous[1]:
                    bases[pos] = (base, quality)
                elif quality == previous[1] and base != previous[0]:
                    bases[pos] = ('N', quality)
    return (chrom, bases, orientation, start, seq)


class _TapsReferenceWindow:
    def __init__(self, reference):
        self.reference = reference
        self.blocks = OrderedDict()

    def get_reference_length(self, chrom):
        return self.reference.get_reference_length(chrom)

    def fetch(self, chrom, start, end):
        block = start // 262144 * 262144
        bound = min(self.reference.get_reference_length(chrom), block + 262144)
        if end > bound:
            return self.reference.fetch(chrom, start, end)
        key = chrom, block
        if key not in self.blocks:
            self.blocks[key] = self.reference.fetch(chrom, block, bound)
            if len(self.blocks) > 16:
                self.blocks.popitem(last=False)
        else:
            self.blocks.move_to_end(key)
        return self.blocks[key][start - block:end - block]


def taps_non_cpg_metrics(bam: str, reference: str, samtools: str, scratch: str, *, expected_pairs: int) -> dict[str, object]:
    """统计可判链分子的宿主 CHG/CHH 转换率；未知或冲突链不随机分配，重叠 mate 每位置只计一次。

    这是未经 SNP/测序错误校正的原始修饰信号 QC；只用高质量 CpG 转换证据判链，
    指标只代表可判链子集，不能直接作为全基因组真实 mCH 或 cDNA 污染率。
    """
    import itertools
    import pysam

    counts = dict.fromkeys(("non_cpg_candidate_pairs", "non_cpg_oriented_pairs", "non_cpg_ambiguous_pairs",
                           "non_cpg_chg_mod", "non_cpg_chg_unmod", "non_cpg_chh_mod", "non_cpg_chh_unmod"), 0)
    if type(expected_pairs) is not int or expected_pairs < 0:
        raise ValueError("TAPS non-CpG QC requires a nonnegative complete pair count")
    primary_pairs = 0
    process = subprocess.Popen([samtools, "collate", "-f", "-O", "-u", "-T", str(Path(scratch) / "non_cpg_collate"), bam], stdout=subprocess.PIPE)
    try:
        with pysam.FastaFile(reference) as fasta, pysam.AlignmentFile(process.stdout, "rb") as source:
            ref = _TapsReferenceWindow(fasta)
            for name, group in itertools.groupby(source, key=lambda r: r.query_name):
                primary = [r for r in group if not r.is_secondary and not r.is_supplementary]
                if len(primary) != 2 or {r.is_read1 for r in primary} != {True, False} or any(
                    not r.is_paired or r.is_read1 == r.is_read2 for r in primary
                ):
                    raise ValueError(f"TAPS non-CpG QC requires a complete pair: {name}")
                primary_pairs += 1
                chrom, bases, orientation, start, seq = _taps_pair_bases(primary, ref)
                if re.fullmatch(r"chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)", chrom) is None:
                    continue
                counts["non_cpg_candidate_pairs"] += 1
                if not orientation:
                    counts["non_cpg_ambiguous_pairs"] += 1
                    continue
                counts["non_cpg_oriented_pairs"] += 1
                chg_mod = chg_unmod = chh_mod = chh_unmod = 0
                if orientation == "OT":
                    limit = len(seq) - 2
                    for pos, (base, _) in bases.items():
                        if base not in "CT":
                            continue
                        offset = pos - start
                        if offset >= limit or seq[offset] != "C" or seq[offset + 1] not in "ACT":
                            continue
                        third = seq[offset + 2]
                        if third == "G":
                            if base == "T":
                                chg_mod += 1
                            else:
                                chg_unmod += 1
                        elif third in "ACT":
                            if base == "T":
                                chh_mod += 1
                            else:
                                chh_unmod += 1
                else:
                    for pos, (base, _) in bases.items():
                        if base not in "GA":
                            continue
                        offset = pos - start
                        if offset < 2 or seq[offset] != "G" or seq[offset - 1] not in "AGT":
                            continue
                        third = seq[offset - 2]
                        if third == "C":
                            if base == "A":
                                chg_mod += 1
                            else:
                                chg_unmod += 1
                        elif third in "AGT":
                            if base == "A":
                                chh_mod += 1
                            else:
                                chh_unmod += 1
                counts["non_cpg_chg_mod"] += chg_mod
                counts["non_cpg_chg_unmod"] += chg_unmod
                counts["non_cpg_chh_mod"] += chh_mod
                counts["non_cpg_chh_unmod"] += chh_unmod
        if process.wait() != 0:
            raise RuntimeError("samtools collate failed during TAPS non-CpG QC")
        if primary_pairs != expected_pairs:
            raise ValueError("TAPS non-CpG primary pair count disagrees with alignment QC")
    finally:
        _close_native_stream(process)
    modified = counts["non_cpg_chg_mod"] + counts["non_cpg_chh_mod"]
    total = modified + counts["non_cpg_chg_unmod"] + counts["non_cpg_chh_unmod"]
    return {**counts, "non_cpg_observations": total,
            "non_cpg_methylation_pct": 100.0 * modified / total if total else float("nan"),
            "non_cpg_metric_source": "taps_oriented_host_chg_chh_conversion"}


def emit_taps_snps(vcf: str, output: str, sample: str, *, reference: str) -> int:
    """校验原生占位样本及候选位点参考，保留 PASS SNP 注释并绑定当前细胞、建立索引。"""
    import pysam
    import pysam.bcftools

    count = 0
    with tempfile.TemporaryDirectory(prefix=".taps_select.", dir=Path(output).parent) as selected:
        selected_vcf = str(Path(selected) / "variants.vcf")
        pysam.bcftools.view("--no-version", "--no-update", "--include", 'FILTER="PASS" && (GT="alt" || GT="mis")',
                            "--output-type", "v", "--output", selected_vcf, vcf, catch_stdout=False)
        with pysam.VariantFile(selected_vcf) as source, pysam.FastaFile(reference) as ref:
            if list(source.header.samples) != ["sample"]:
                raise ValueError("TAPS SNP VCF requires the native Rastair placeholder sample")
            old_sample = next(iter(source.header.samples))
            reference_lengths = dict(zip(ref.references, ref.lengths))
            reference_bases = _TapsReferenceWindow(ref)
            header_text = str(source.header).rsplit("\t", 1)[0] + "\t" + sample + "\n"
            with pysam.BGZFile(output, "w") as target:
                target.write(header_text.encode())
                for row in source:
                    alleles = row.alleles or ()
                    genotype = row.samples[old_sample].get("GT") or ()
                    if (set(row.filter) != {"PASS"} or len(alleles) < 2
                            or any(len(a) != 1 or a not in "ACGT" for a in alleles)
                            or not any(a is not None and a > 0 for a in genotype)):
                        continue
                    if (row.contig not in reference_lengths or row.start < 0
                            or row.stop > reference_lengths[row.contig]):
                        raise ValueError(f"TAPS SNP coordinate is outside reference FASTA: {row.contig}:{row.pos}")
                    if reference_bases.fetch(row.contig, row.start, row.stop).upper() != row.ref:
                        raise ValueError(f"TAPS SNP REF disagrees with reference FASTA: {row.contig}:{row.pos}")
                    target.write(str(row).encode())
                    count += 1
    with tempfile.TemporaryDirectory(prefix=".taps_snp.", dir=Path(output).parent) as scratch:
        corrected = str(Path(scratch) / "snps.vcf.gz")
        pysam.bcftools.reheader("--fai", f"{reference}.fai", "--output", corrected, output, catch_stdout=False)
        os.replace(corrected, output)
    pysam.tabix_index(output, preset="vcf", force=True)
    return count


def emit_rastair_cpg(bed: str, reference: str, alignment_qc: str, output_qc: str, stream: TextIO,
                     *, bam: str, samtools: str, scratch: str, vcf: str, output_snp: str, sample: str) -> None:
    """输出共用 CpG writer 的整数计数输入；CpG 任一端已调用为变异时排除两端，不使用 beta。"""
    import pysam

    qc = json.loads(Path(alignment_qc).read_text())
    if qc["sample"] != sample:
        raise ValueError("TAPS caller sample differs from alignment QC")
    qc.update(taps_non_cpg_metrics(bam, reference, samtools, scratch, expected_pairs=qc["trimmed_pairs"]))
    qc["snp_sites"] = emit_taps_snps(vcf, output_snp, sample, reference=reference)
    counts = {"native_cpg_rows": 0, "canonical_cpg_rows": 0, "excluded_variant_or_denovo_rows": 0,
              "zero_observation_rows": 0}
    controls = {name: {"mod": 0, "unmod": 0, "sites": 0} for name in ("lambda", "pUC19")}
    required = {"#chr", "start", "end", "unmod", "mod", "coverage", "genotype", "cpg"}
    with pysam.FastaFile(reference) as ref, Path(bed).open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        fieldnames = next(reader, [])
        column = {name: i for i, name in enumerate(fieldnames)}
        if not required.issubset(fieldnames) or len(fieldnames) != len(set(fieldnames)):
            raise ValueError("Rastair BED lacks the native 2.2 header")
        previous = None
        contigs = {name: i for i, name in enumerate(ref.references)}
        pending_key = None
        pending = []
        reference_chrom = ""
        reference_start = reference_end = 0
        reference_sequence = ""

        def flush_cpg() -> None:
            if any(gt != base + "/" + base for _, _, _, _, _, gt, base in pending):
                counts["excluded_variant_or_denovo_rows"] += len(pending)
                return
            for chrom, pos, mod, unmod, depth, gt, base in pending:
                if depth == 0 or mod + unmod == 0:
                    counts["zero_observation_rows"] += 1
                    continue
                stream.write(f"{chrom}\t{pos}\t{pos + 1}\t.\t{mod}\t{unmod}\n")
                counts["canonical_cpg_rows"] += 1
                if chrom in controls:
                    controls[chrom]["mod"] += mod
                    controls[chrom]["unmod"] += unmod
                    controls[chrom]["sites"] += 1

        for row in reader:
            counts["native_cpg_rows"] += 1
            label = f"Rastair BED {bed} line {reader.line_num}"
            if len(row) != len(fieldnames) or any(row[column[key]] == "" for key in required):
                raise ValueError(f"{label} has missing or extra fields")
            origin = row[column["cpg"]]
            if origin not in {"REF", "NEW"}:
                raise ValueError(f"{label} has invalid CpG origin: {origin!r} (requires REF or NEW)")
            chrom = row[column["#chr"]]
            if chrom not in contigs:
                raise ValueError(f"{label} contig {chrom!r} is absent from the reference FASTA")
            numeric = [row[column[k]] for k in ("start", "end", "mod", "unmod", "coverage")]
            if any(re.fullmatch(r"[0-9]+", value) is None for value in numeric):
                raise ValueError(f"{label} has non-integer coordinates/counts (requires nonnegative ASCII integers)")
            pos, end, mod, unmod, depth = map(int, numeric)
            if end != pos + 1 or end > ref.get_reference_length(chrom) or mod + unmod > depth:
                raise ValueError(f"{label}: Invalid Rastair coordinate/counts: {chrom}:{pos}")
            key = (contigs[chrom], pos)
            if previous is not None and key <= previous:
                raise ValueError("Rastair BED is not ordered or contains duplicate positions")
            previous = key
            if origin == "NEW":
                counts["excluded_variant_or_denovo_rows"] += 1
                continue
            if chrom != reference_chrom or not reference_start < pos < reference_end - 1:
                reference_chrom = chrom
                reference_start = max(0, pos - 1)
                reference_end = min(ref.get_reference_length(chrom), pos + 1048577)
                reference_sequence = ref.fetch(chrom, reference_start, reference_end).upper()
            offset = pos - reference_start
            base = reference_sequence[offset:offset + 1]
            start = pos if base == "C" else pos - 1
            cpg_offset = start - reference_start
            if (base not in {"C", "G"} or start < 0
                    or reference_sequence[cpg_offset:cpg_offset + 2] != "CG"):
                raise ValueError(f"Rastair reference CpG disagrees with FASTA: {chrom}:{pos}")
            cpg_key = (chrom, start)
            if cpg_key != pending_key:
                flush_cpg()
                pending = []
                pending_key = cpg_key
            pending.append((chrom, pos, mod, unmod, depth, row[column["genotype"]], base))
        flush_cpg()
    qc.update(counts)
    qc.update({"signal": "5mC+5hmC", "high_cph_status": "not_applicable",
               "min_mapq": TAPS_MIN_MAPQ, "min_baseq": TAPS_MIN_BASEQ,
               "include_flags": TAPS_INCLUDE_FLAGS, "exclude_flags": TAPS_EXCLUDE_FLAGS,
               "guess_read_orientation": True, "mate_overlap_counted_once": True})
    for name, value in controls.items():
        denominator = value["mod"] + value["unmod"]
        qc[f"{name}_cpg_mod"] = value["mod"]
        qc[f"{name}_cpg_observations"] = denominator
        qc[f"{name}_cpg_sites"] = value["sites"]
        metric = "false_positive_pct" if name == "lambda" else "conversion_pct"
        qc[f"{name}_{metric}"] = 100.0 * value["mod"] / denominator if denominator else None
    atomic_write_json(output_qc, qc)
