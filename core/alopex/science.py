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

_LOADED_MODULE_CODE = sys._getframe().f_code
import tempfile
from collections import OrderedDict
from collections.abc import Iterator, Mapping, Sequence
from contextlib import ExitStack, contextmanager, nullcontext
from dataclasses import dataclass
from pathlib import Path
from typing import Any, TextIO

from .config import (
    DD_CONSENSUS,
    DD_PROTOCOLS,
    FILTER_SUMMARY_FIELDS,
    FILTER_SUMMARY_SCHEMA_VERSION,
    SAMPLE_ID_RE,
    get_protocol,
    primary_analysis_contigs,
    _integer,
)
from .io import (
    project_temporary_parent,
    publication_locks,
    require_distinct_paths,
    _close_native_stream,
    _read_json,
    atomic_write_json,
)


# 科学处理：共用统计

def _native_count(token: str, *, grouped: bool = False) -> int:
    """Parse an unsigned ASCII count, optionally accepting grouped thousands."""
    if grouped and "," in token and re.fullmatch(r"[0-9]{1,3}(?:,[0-9]{3})+", token):
        token = token.replace(",", "")
    if not token.isascii() or not token.isdecimal():
        raise ValueError(f"Invalid native count: {token!r}")
    return int(token)


def _bismark_field(text: str, field: str, path: Path) -> str:
    values = {v.strip(" \t") for v in re.findall(rf"^{re.escape(field)}:([^\r\n]*)", text, re.MULTILINE)}
    if len(values) != 1 or "" in values:
        raise ValueError(f"Missing or conflicting Bismark {field} in {path}")
    return values.pop()


def _bismark_report_count(text: str, field: str, path: Path) -> int:
    try:
        return _native_count(_bismark_field(text, field, path))
    except ValueError as exc:
        raise ValueError(f"Invalid Bismark count for {field} in {path}") from exc


def _count_bam_pairs(bam: Path, samtools: str) -> int:
    result = subprocess.run(
        [samtools, "view", "-c", str(bam)],
        check=True, text=True, stdout=subprocess.PIPE,
    )
    records = int(result.stdout.strip())
    if records < 0 or records % 2:
        raise ValueError("Bismark processed BAM contains an odd or negative alignment count")
    return records // 2


def _parse_cutadapt_metrics(path: Path, mates: tuple[str, str], payload: Mapping[str, Any]) -> dict[str, int | float]:
    """核对原生报告与本次拆分的双端路径。交付时 FASTQ 可能已回收。"""
    inputs = payload.get("input")
    if (not isinstance(inputs, dict) or inputs.get("paired") is not True
            or (inputs.get("path1"), inputs.get("path2")) != mates):
        raise ValueError(f"Cutadapt report does not match the committed R1/R2 inputs: {path}")
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
    filtered = counts.get("filtered")
    if not isinstance(filtered, dict) or "too_short" not in filtered:
        raise ValueError(f"Cutadapt JSON lacks read_counts.filtered: {path}")
    too_short = filtered["too_short"]
    if (type(too_short) is not int or too_short < 0 or too_short + values["output"] != values["input"]
            or any(value is not None for key, value in filtered.items() if key != "too_short")):
        raise ValueError(f"Cutadapt filtered and retained pairs disagree with input: {path}")
    # 各协议每端只传一个接头，且每条 read 至多修剪一次。
    for mate in (1, 2):
        adapters = payload.get(f"adapters_read{mate}")
        if (not isinstance(adapters, list) or len(adapters) != 1 or not isinstance(adapters[0], dict)
                or type(adapters[0].get("total_matches")) is not int
                or adapters[0]["total_matches"] != values[f"read{mate}_with_adapter"]):
            raise ValueError(f"Cutadapt R{mate} adapter matches disagree with its read count: {path}")
    return {
        "cutadapt_input_pairs": values["input"],
        "trimmed_pairs": values["output"],
        "r1_adapter_pct": 100.0 * values["read1_with_adapter"] / values["input"],
        "r2_adapter_pct": 100.0 * values["read2_with_adapter"] / values["input"],
    }


def _validate_high_cph_route(backend: str, protocol: str, role: str) -> None:
    protocols = {"biscuit": ("cabernet",), "bismark": ("cabernet", *DD_PROTOCOLS)}
    if protocol not in protocols.get(backend, ()) or role != get_protocol(protocol).high_cph_role:
        raise ValueError(f"Invalid high-CpH protocol/role identity: {backend=}, {protocol=}, {role=}")


def _read_high_cph_summary(path: Path, sample: str, backend: str) -> dict[str, object]:
    if backend == "bismark":
        return _read_filter_summary(path, sample)
    payload = _read_json(path, "high-CpH summary")
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
        type(payload["schema_version"]) is not int
        or payload["schema_version"] != 2
        or payload["backend"] != backend
        or payload["sample"] != sample
        or payload["metric_unit"] != "read_pairs"
    ):
        raise ValueError(f"High-CpH summary identity is invalid: {path}")
    _validate_high_cph_route(backend, payload["protocol"], payload["high_cph_role"])
    candidate = _integer(payload["candidate_read_pairs"], "candidate_read_pairs", path)
    flagged = _integer(payload["flagged_read_pairs"], "flagged_read_pairs", path)
    excluded = _integer(payload["excluded_read_pairs"], "excluded_read_pairs", path)
    if flagged > candidate:
        raise ValueError(f"High-CpH flagged pairs exceed candidate pairs: {path}")
    if payload["filter_applied"] is not True:
        raise ValueError(f"High-CpH summary must record the fixed filter action: {path}")
    fraction = payload["high_cph_fraction"]
    if isinstance(fraction, bool) or not isinstance(fraction, (int, float)):
        raise ValueError(f"High-CpH summary fraction is invalid: {path}")
    expected_fraction = flagged / candidate if candidate else 0.0
    # 两个后端均以整数计数直接相除并写入 JSON，读回的浮点值可精确复算。
    if fraction != expected_fraction:
        raise ValueError(f"High-CpH summary fraction is inconsistent: {path}")
    threshold = payload["threshold"]
    if (
        isinstance(threshold, bool)
        or not isinstance(threshold, (int, float))
        or not 0 <= threshold <= 1
    ):
        raise ValueError(f"High-CpH summary threshold is invalid: {path}")
    metric_rows = _integer(payload["observed_metric_rows"], "observed_metric_rows", path)
    high_rows = _integer(payload["flagged_metric_rows"], "flagged_metric_rows", path)
    excluded_rows = _integer(payload["excluded_rows"], "excluded_rows", path)
    missing = _integer(payload["missing_metric_rows"], "missing_metric_rows", path)
    total = _integer(payload["observed_rows"], "observed_rows", path)
    # 每个类别的 QNAME 来自该类原始行去重。未标记的模板至少各有一行未超阈值。
    for pairs, rows in ((candidate, metric_rows), (flagged, high_rows), (excluded, excluded_rows)):
        if pairs > rows or bool(pairs) != bool(rows):
            raise ValueError(f"BISCUIT High-CpH row accounting is inconsistent: {path}")
    if (
        missing + metric_rows + excluded_rows != total
        or candidate - flagged > metric_rows - high_rows
        or (threshold == 1 and flagged)
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


_BISCUIT_CONTEXT_TAGS = re.compile(rb"(?:^|;)(CX|N5)(?:=([^;]*))?(?=;|$)")
_BISCUIT_FIVE_BASES = re.compile(rb"[ACGTRYSWKMBDHVN]{2}C[ACGTRYSWKMBDHVN]{2}")


def _biscuit_is_cpg(ref: bytes, info: bytes, path: Path, line_no: int) -> bool:
    """按原生 fivenuc_context 核对 CX。负链的 N5 已反向互补，边界 N 对应 CN。"""
    tags = _BISCUIT_CONTEXT_TAGS.findall(info)
    if ref not in (b"C", b"G"):
        if tags:
            raise ValueError(f"BISCUIT VCF {path} line {line_no}: context on a non-cytosine reference")
        return False
    context = dict(tags)
    window = context.get(b"N5", b"")
    if len(tags) != 2 or b"CX" not in context or not _BISCUIT_FIVE_BASES.fullmatch(window):
        raise ValueError(f"BISCUIT VCF {path} line {line_no}: expected one CX and one cytosine-centred N5")
    if b"N" in window:
        expected = b"CN"
    elif window[3:4] == b"G":
        expected = b"CG"
    elif window[4:5] == b"G":
        expected = b"CHG"
    else:
        expected = b"CHH"
    if context[b"CX"] != expected:
        raise ValueError(f"BISCUIT VCF {path} line {line_no}: CX={context[b'CX']!r} disagrees with N5={window!r}")
    return expected == b"CG"


def _biscuit_pileup_records(
    vcf: str | Path, stats: str | Path, sample: str,
) -> Iterator[tuple[bytes, list[bytes] | None]]:
    """流式核对样本、化学体系、上下文与记录数，附上已拆分的 CpG 行供整数提取。"""
    expected = read_pileup_variant_count(stats, sample)
    path = Path(vcf)
    opener = gzip.open if path.suffix == ".gz" else open
    # pileup 由输入 BAM 文件名派生样本列，生产 BAM 固定为 <Sample_ID>.dna.bam。
    header = b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + f"{sample}.dna".encode()
    with opener(path, "rb") as handle:
        chemistry_seen = False
        for line_no, line in enumerate(handle, 1):
            if line.startswith(b"##"):
                if re.match(rb"##chemistry(?:[=\s]|$)", line):
                    chemistry = line.rstrip(b"\r\n").removeprefix(b"##chemistry=")
                    if chemistry_seen or chemistry not in (b"bisulfite", b"em-seq"):
                        raise ValueError(f"BISCUIT VCF {path} line {line_no}: duplicate or incompatible chemistry {chemistry!r}")
                    chemistry_seen = True
                yield line, None
                continue
            if line.rstrip(b"\r\n") != header:
                raise ValueError(f"BISCUIT VCF {path} line {line_no}: expected one native sample {sample}.dna")
            yield line, None
            break
        else:
            raise ValueError(f"BISCUIT VCF {path}: missing native sample header")

        count = 0
        for count, line in enumerate(handle, 1):
            fields = line.rstrip(b"\r\n").split(b"\t")
            if line.startswith(b"#") or len(fields) != 10:
                raise ValueError(f"BISCUIT VCF {path} line {line_no + count}: malformed single-sample record")
            if count > expected:
                raise ValueError(f"BISCUIT VCF {path}: more than {expected} records declared for {sample}")
            is_cpg = _biscuit_is_cpg(fields[3], fields[7], path, line_no + count)
            yield line, fields if is_cpg else None
        if count != expected:
            raise ValueError(f"BISCUIT VCF {path}: expected {expected} records for {sample}, observed {count}")


def biscuit_pileup_lines(vcf: str | Path, stats: str | Path, sample: str) -> Iterator[bytes]:
    """将已校验的原始记录逐行交给 SNP 和 QC 的原生消费者。"""
    for line, _ in _biscuit_pileup_records(vcf, stats, sample):
        yield line


def emit_biscuit_cpg(vcf: str | Path, stream: TextIO, *, stats: str | Path, sample: str) -> None:
    """从单样本原生 VCF 的 SP/CV 输出 CpG 原始整数，不从三位小数 BT 反推。

    BISCUIT 同一过滤循环中，参考 C 的 BSW 转换计入 SP.Y，参考 G 的
    BSC 转换计入 SP.R；CV 为可调用链的 retention+conversion。相反链的
    非信息观测也计入 SP.C/G，因此参考碱基支持仅是 retention 的上界。
    保留原生可调用/SNP/质量/重叠过滤、逐链位置及文件顺序，不合并 CpG dyad。
    """
    for line_no, (_, fields) in enumerate(_biscuit_pileup_records(vcf, stats, sample), 1):
        if fields is None:
            continue
        label = f"BISCUIT VCF {vcf} line {line_no}"
        if not fields[1].isdigit() or int(fields[1]) < 1:
            raise ValueError(f"{label}: invalid CpG position")
        names, values = fields[8].split(b":"), fields[9].split(b":")
        if len(names) != len(set(names)) or len(names) != len(values):
            raise ValueError(f"{label}: malformed or duplicate FORMAT fields")
        calls = dict(zip(names, values))
        # 原生 SNP 干扰可使此行没有 CV；CV=0 同样表示本样本不可调用。
        if b"CV" not in calls:
            continue
        if not calls[b"CV"].isdigit():
            raise ValueError(f"{label}: invalid integer CV")
        coverage = int(calls[b"CV"])
        if not coverage:
            continue
        support = calls.get(b"SP", b"")
        tokens = re.findall(rb"([ACGTYR])([0-9]+)", support)
        if (not tokens or b"".join(base + count for base, count in tokens) != support
                or len({base for base, count in tokens}) != len(tokens)):
            raise ValueError(f"{label}: malformed or duplicate integer SP support")
        counts = {base: int(count) for base, count in tokens}
        conversion = counts.get(b"Y" if fields[3] == b"C" else b"R", 0)
        retention = coverage - conversion
        if not 0 <= retention <= counts.get(fields[3], 0):
            raise ValueError(f"{label}: inconsistent SP/CV counts")
        start = int(fields[1]) - 1
        stream.write(f"{fields[0].decode()}\t{start}\t{start + 1}\t.\t{retention}\t{conversion}\n")


_SAM_QNAME = re.compile(r"[!-?A-~]{1,254}")


def _require_sam_qname(name: str | None, source: str | Path, line: int | None = None) -> str:
    """依赖模板身份的步骤须按完整 SAM QNAME 匹配；单独的 * 表示名称缺失。"""
    if not name or name == "*" or _SAM_QNAME.fullmatch(name) is None:
        location = f"{source}:{line}" if line is not None else str(source)
        raise ValueError(f"Invalid or unavailable SAM QNAME in {location}: {name!r}")
    return name


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
    with ExitStack() as temporary, tempfile.TemporaryFile(mode="w+t", errors="replace", dir=project_temporary_parent(bam_path)) as diagnostics:
        command = [samtools, "view"]
        if not any(name.startswith('#') for name, _ in contigs):
            regions = temporary.enter_context(tempfile.NamedTemporaryFile(mode="w+t", encoding="utf-8", suffix=".bed", dir=project_temporary_parent(bam_path)))
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
                    names.add(_require_sam_qname(fields[0], bam_path))
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
    return {
        "id": "bsconv_filter",
        "section_name": "High-CpH Read Assessment",
        "description": (
            f"Protocol={summary['protocol']}; high-CpH reads are interpreted as candidate cDNA contamination. "
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
    env["TMPDIR"] = str(project_temporary_parent(destination))
    fd, temporary = tempfile.mkstemp(prefix=f".{destination.name}.", suffix=".tmp", dir=project_temporary_parent(destination))
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
            _require_sam_qname(left.query_name, source)
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
    _validate_high_cph_route("biscuit", args.protocol, args.high_cph_role)
    if isinstance(args.threshold, bool) or not isinstance(args.threshold, (int, float)) or not 0 <= args.threshold <= 1:
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
            dir=project_temporary_parent(args.high_read_names),
            delete=False,
        ) as raw_high_handle, tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=args.high_read_names.name + ".",
            suffix=".candidate.unsorted",
            dir=project_temporary_parent(args.high_read_names),
            delete=False,
        ) as raw_candidate_handle, tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=args.high_read_names.name + ".",
            suffix=".excluded.unsorted",
            dir=project_temporary_parent(args.high_read_names),
            delete=False,
        ) as raw_excluded_handle:
            raw_high_path = Path(raw_high_handle.name)
            raw_candidate_path = Path(raw_candidate_handle.name)
            raw_excluded_path = Path(raw_excluded_handle.name)
            with (nullcontext(sys.stdin) if str(args.bsconv_tsv) == "-"
                  else args.bsconv_tsv.open(encoding="utf-8")) as handle:
                # 原生 bsconv 不转义 QNAME；合法双引号属于读段身份，不能按 CSV 解码。
                reader = csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
                for row in reader:
                    if len(row) != BSCONV_P_COLUMNS or not all(
                        re.fullmatch(r"[0-9]+", value) for value in row[:8]
                    ):
                        raise ValueError(
                            f"Unexpected biscuit bsconv -p row at {args.bsconv_tsv}:{reader.line_num} "
                            f"(expect 8 nonnegative integer cytosine counts + QNAME): {row!r}"
                        )

                    read_name = _require_sam_qname(row[8], args.bsconv_tsv, reader.line_num)
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
        atomic_write_json(args.multiqc_json, _high_cph_multiqc_payload(args.sample_name, summary))


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
    cutadapt_report: Path,
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
    for line_no, line in enumerate(stat.read_text().splitlines(), start=1):
        prefix, _, content = line.partition("]")
        key, _, value = content.strip().partition(":")
        field = mapping.get(key.strip())
        if prefix != "[dupsifter" or field is None:
            continue
        try:
            count = _native_count(value.strip(), grouped=True)
        except ValueError as exc:
            raise ValueError(f"Invalid Dupsifter count at {stat}:{line_no}: {key}={value!r}") from exc
        if field in metrics and metrics[field] != count:
            raise ValueError(f"Conflicting Dupsifter count at {stat}:{line_no}: {key}={value!r}")
        metrics[field] = count
    cutadapt = _parse_cutadapt_metrics(cutadapt_report, (row["dna_r1"], row["dna_r2"]),
                                      _read_json(cutadapt_report, "Cutadapt report"))
    _dupsifter_pair_counts(metrics, stat, cutadapt["trimmed_pairs"])
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
        "data": {demux_sample: output_row},
    }
    atomic_write_json(out_json, obj)


def _parse_cph_methylation_pct(path: Path) -> float:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0].strip() != "BISCUITqc CpH Retention by Read Position Table":
        raise ValueError(f"Malformed BISCUIT CpH retention title: {path}")
    if len(lines) < 2:
        raise ValueError(f"BISCUIT CpH retention table lacks a header: {path}")
    reader = csv.DictReader(lines[1:], delimiter="\t", quoting=csv.QUOTE_NONE)
    required = {"ReadInPair", "Position", "Conversion/Retention", "Count"}
    fields = reader.fieldnames or []
    if len(fields) != len(set(fields)) or not required.issubset(fields):
        raise ValueError(f"Malformed BISCUIT CpH retention header: {path}")
    seen: set[tuple[str, int, str]] = set()
    retained = total = 0
    for row in reader:
        if None in row or None in row.values():
            raise ValueError(f"Malformed BISCUIT CpH retention row width at {path}:{reader.line_num + 1}")
        read = row["ReadInPair"].strip()
        state = row["Conversion/Retention"].strip()
        if read not in {"1", "2"} or state not in {"C", "R"}:
            raise ValueError(f"Malformed BISCUIT CpH retention row: {path}: {row}")
        try:
            position = _native_count(row["Position"].strip())
            count = _native_count(row["Count"].strip())
        except ValueError as exc:
            raise ValueError(f"Invalid BISCUIT CpH count/position: {path}: {row}") from exc
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
    payload = _read_json(path, "Dupsifter QC")
    data = payload.get("data")
    if not isinstance(data, dict) or set(data) != {sample}:
        raise ValueError(f"Dupsifter QC must contain exactly sample {sample!r}: {path}")
    row = data[sample]
    if not isinstance(row, dict):
        raise ValueError(f"Dupsifter QC sample row is invalid: {path}")
    return row, _dupsifter_pair_counts(row, path, trimmed_pairs)


def _dupsifter_pair_counts(row: Mapping[str, object], path: Path, trimmed_pairs: int) -> dict[str, int]:
    for field in DUPSIFTER_MQC_FIELDS:
        value = row.get(field)
        if type(value) is not int or value < 0:
            raise ValueError(f"Dupsifter QC has invalid {field}: {path}")
    if row["Reads_Processed"] != row["Primary_Reads"] + row["Secondary_Supplementary"]:
        raise ValueError(f"Dupsifter alignment record counts disagree: {path}")
    for dup, mapped in (("Dup_Both", "Both_Mapped"), ("Dup_Forward", "One_Mapped_Forward"),
                        ("Dup_Reverse", "One_Mapped_Reverse")):
        if row[dup] > row[mapped]:
            raise ValueError(f"Dupsifter {dup} exceeds {mapped}: {path}")
    both = row["Both_Mapped"]
    one = row["One_Mapped_Forward"] + row["One_Mapped_Reverse"]
    unmapped = row["No_Mapped"]
    no_primary = row["No_Primary"]
    accepted = both + one
    rejected = unmapped + no_primary
    # 两端均已比对或均未比对的模板各占两条 primary。单端比对可另含一条未比对 mate。
    min_primary = 2 * (both + unmapped) + one
    if not min_primary <= row["Primary_Reads"] <= min_primary + one or no_primary > row["Secondary_Supplementary"]:
        raise ValueError(f"Dupsifter pair categories disagree with alignment records: {path}")
    duplicate = row["Dup_Both"] + row["Dup_Forward"] + row["Dup_Reverse"]
    if accepted + rejected != trimmed_pairs:
        raise ValueError(
            "Dupsifter pair categories disagree with Cutadapt output: "
            f"accepted({accepted}) + rejected({rejected}) != trimmed({trimmed_pairs})"
        )
    return {
        "both_primary_mapped_pairs": both,
        "one_primary_mapped_pairs": one,
        "backend_accepted_pairs": accepted,
        "backend_rejected_pairs": rejected,
        "unmapped_pairs": unmapped,
        "no_primary_pairs": no_primary,
        "duplicate_pairs": duplicate,
        "postdedup_pairs": accepted - duplicate,
    }


def _biscuit_mapq_counts(path: Path, text: str) -> dict[str, int]:
    """受管版 biscuit qc 为 unmapped 和 MAPQ 0–60 各写一行，包括零计数。"""
    lines = text.splitlines()
    if lines[:2] != ["BISCUITqc Mapping Quality Table", "MapQ\tCount"]:
        raise ValueError(f"Malformed BISCUIT MAPQ table header: {path}")
    counts: dict[str, int] = {}
    for line in lines[2:]:
        fields = line.split()
        if len(fields) != 2:
            raise ValueError(f"Malformed BISCUIT MAPQ row in {path}: {line!r}")
        label, raw_count = fields
        try:
            count = _native_count(raw_count, grouped=True)
        except ValueError as exc:
            raise ValueError(f"Invalid BISCUIT MAPQ count in {path}: {line!r}") from exc
        if label != "unmapped":
            try:
                label = str(_native_count(label))
            except ValueError as exc:
                raise ValueError(f"Invalid BISCUIT MAPQ label in {path}: {label!r}") from exc
        if label in counts:
            raise ValueError(f"Duplicate BISCUIT MAPQ category in {path}: {label!r}")
        counts[label] = count
    if counts.keys() != {"unmapped", *map(str, range(61))}:
        raise ValueError(f"BISCUIT MAPQ table requires unmapped and all bins 0–60: {path}")
    return counts


def _read_biscuit_native_mapping(path: Path) -> dict[str, object]:
    counts = _biscuit_mapq_counts(path, path.read_text(encoding="utf-8"))
    total = sum(counts.values())
    unaligned = counts["unmapped"]
    optimal = sum(counts[str(q)] for q in range(40, 61))
    return {
        "optimally_aligned_reads": optimal,
        "suboptimally_aligned_reads": total - unaligned - optimal,
        "unaligned_native_reads": unaligned,
        "native_mapping_policy": "biscuit_nonsecondary_mapq_ge_40",
        "native_mapping_unit": "nonsecondary_alignment_records_postdedup",
        "native_mapping_pct": 100.0 * optimal / total if total else 0.0,
    }


# 科学处理：Bismark 与 Droplet UMI

def _bismark_trimmed_mates(cutadapt: Mapping[str, Any], path: Path) -> tuple[str, str]:
    """Cutadapt records both output paths verbatim even after job scratch is reclaimed."""
    argv = cutadapt.get("command_line_arguments")
    if not isinstance(argv, list) or not all(isinstance(arg, str) for arg in argv):
        raise ValueError(f"Missing Cutadapt command line for Bismark: {path}")
    if argv[-2:] != [cutadapt["input"]["path1"], cutadapt["input"]["path2"]]:
        raise ValueError(f"Cutadapt command line disagrees with its paired inputs: {path}")
    mates = []
    for flag in ("-o", "-p"):
        if argv.count(flag) != 1:
            raise ValueError(f"Missing or repeated Cutadapt output option {flag}: {path}")
        index = argv.index(flag) + 1
        if index >= len(argv) or not Path(argv[index]).is_absolute():
            raise ValueError(f"Invalid Cutadapt output path for {flag}: {path}")
        mates.append(argv[index])
    if mates[0] == mates[1]:
        raise ValueError(f"Cutadapt paired outputs coincide: {path}")
    return mates[0], mates[1]


def _bismark_cytosine_counts(text: str, path: Path, *, unmethylated_field: str,
                            pairs: int) -> dict[str, tuple[int, int]]:
    """核对报告显示值，保留每种 context 的修饰和未修饰观测数。"""
    total_calls = _bismark_report_count(text, "Total number of C's analysed", path)
    lines = set(text.split("\n"))
    counts = {}
    for context in ("CpG", "CHG", "CHH"):
        label = f"{context} context"
        m = _bismark_report_count(text, f"Total methylated C's in {label}", path)
        u = _bismark_report_count(text, f"{unmethylated_field} in {label}", path)
        counts[context] = (m, u)
        n = m + u
        field = f"C methylated in {label}"
        missing = f"Can't determine percentage of methylated Cs in {label} if value was 0" in lines
        if n:
            valid = not missing and _bismark_field(text, field, path) == f"{m * 100.0 / n:.1f}%"
        else:
            valid = missing and not re.search(rf"^{re.escape(field)}:", text, re.MULTILINE)
        if not valid:
            raise ValueError(f"Bismark {context} percentage disagrees with its counts: {path}")
    if sum(m + u for m, u in counts.values()) != total_calls or (pairs == 0 and total_calls):
        raise ValueError(f"Bismark cytosine total disagrees with its counts: {path}")
    return counts


def parse_alignment_report(path: Path, trimmed_mates: tuple[str, str]) -> dict[str, int | float]:
    text = path.read_bytes().decode("utf-8")
    prefix = f"Bismark report for: {trimmed_mates[0]} and {trimmed_mates[1]} (version: "
    if not text.startswith(prefix):
        raise ValueError(f"Bismark alignment report does not match the Cutadapt outputs: {path}")
    # Paths may themselves contain newlines or option-like text. Only parse the report body.
    text = text[len(prefix):].partition("\n")[2]
    if re.search(r"^Bismark report for: ", text, re.MULTILINE):
        raise ValueError(f"Multiple Bismark alignment reports in {path}")
    total = _bismark_report_count(text, "Sequence pairs analysed in total", path)
    unique = _bismark_report_count(text, "Number of paired-end alignments with a unique best hit", path)
    no_alignment = _bismark_report_count(text, "Sequence pairs with no alignments under any condition", path)
    non_unique = _bismark_report_count(text, "Sequence pairs did not map uniquely", path)
    discarded = _bismark_report_count(
        text, "Sequence pairs which were discarded because genomic sequence could not be extracted", path,
    )
    if unique + no_alignment + non_unique != total or discarded > unique:
        raise ValueError(f"Inconsistent Bismark alignment counts in {path}")
    efficiency = f"{unique * 100.0 / total:.1f}" if total else "0"
    if _bismark_field(text, "Mapping efficiency", path) != efficiency + "%":
        raise ValueError(f"Bismark mapping efficiency disagrees with the exact report counts: {path}")
    _bismark_cytosine_counts(text, path, unmethylated_field="Total unmethylated C's", pairs=unique - discarded)
    return {
        "total": total, "unique": unique, "no_alignment": no_alignment, "non_unique": non_unique,
        "discarded": discarded, "mapping_efficiency": float(efficiency),
    }


def parse_dedup_report(path: Path, input_bam: Path) -> dict[str, int]:
    text = path.read_bytes().decode("utf-8")
    header = f"Total number of alignments analysed in {input_bam}:"
    text = text.lstrip("\n")
    if not text.startswith(header):
        raise ValueError(f"Bismark deduplication report names a different input BAM: {path}")
    value, _, text = text[len(header):].partition("\n")
    if re.search(r"^Total number of alignments analysed in ", text, re.MULTILINE):
        raise ValueError(f"Multiple Bismark deduplication reports in {path}")
    try:
        total = _native_count(value.strip(" \t"))
    except ValueError as exc:
        raise ValueError(f"Invalid Bismark deduplication total in {path}") from exc
    counts = {"total": total}
    for key, field, suffix in (
        ("removed", "Total number duplicated alignments removed", "%)"),
        ("leftover", "Total count of deduplicated leftover sequences", "% of total)"),
    ):
        value = _bismark_field(text, field, path)
        token, _, annotation = value.partition(" ")
        try:
            count = _native_count(token)
        except ValueError as exc:
            raise ValueError(f"Invalid Bismark {field} in {path}") from exc
        percent = count / total * 100.0 if total else 0.0
        if annotation != f"({percent:.2f}{suffix}":
            raise ValueError(f"Bismark deduplication percentage disagrees with {field}: {path}")
        counts[key] = count
    if counts["removed"] + counts["leftover"] != total:
        raise ValueError(f"Inconsistent Bismark deduplication counts in {path}")
    return counts


@dataclass(frozen=True)
class BismarkExtractionReport:
    pairs: int
    counts: dict[str, tuple[int, int]]

    @property
    def non_cpg_pct(self) -> float:
        """合并 CHG、CHH 原始观察，以两种 context 的总覆盖为分母。"""
        chg_m, chg_u = self.counts["CHG"]
        chh_m, chh_u = self.counts["CHH"]
        total = chg_m + chh_m + chg_u + chh_u
        return 100.0 * (chg_m + chh_m) / total if total else math.nan


def parse_bismark_extraction_report(path: Path, sample: str) -> BismarkExtractionReport:
    text = path.read_text(encoding="utf-8")
    bam, _, text = text.partition("\n")
    if bam != f"cell_{sample}.processed.bam":
        raise ValueError(f"Bismark extraction report names a different cell BAM: {path}")
    if _bismark_field(text, "Bismark result file", path) != "paired-end (SAM format)":
        raise ValueError(f"Bismark extraction report is not paired-end: {path}")
    processed = re.findall(r"^Processed ([^\r\n]*) lines in total$", text, re.MULTILINE)
    if len(processed) != 1 or not processed[0].isascii() or not processed[0].isdecimal():
        raise ValueError(f"Invalid Bismark processed pair count in {path}")
    pairs = int(processed[0])
    if _bismark_report_count(text, "Total number of methylation call strings processed", path) != 2 * pairs:
        raise ValueError(f"Bismark methylation call strings disagree with the pair count: {path}")
    counts = _bismark_cytosine_counts(text, path, unmethylated_field="Total C to T conversions", pairs=pairs)
    return BismarkExtractionReport(pairs, counts)


def _bismark_pair_reports(alignment_report: Path, dedup_report: Path | None, sample: str,
                          trimmed_mates: tuple[str, str]) -> tuple[dict[str, int | float], dict[str, int]]:
    alignment = parse_alignment_report(alignment_report, trimmed_mates)
    accepted = int(alignment["unique"]) - int(alignment["discarded"])
    if dedup_report is None:
        return alignment, {"total": accepted, "removed": 0, "leftover": accepted}
    input_bam = Path(trimmed_mates[0]).parent / "bismark" / f"{sample}_pe.bam"
    dedup = parse_dedup_report(dedup_report, input_bam)
    if dedup["total"] != accepted:
        raise ValueError(f"Bismark deduplication input disagrees with usable uniquely aligned pairs: {dedup_report}")
    return alignment, dedup


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
    fraction = payload["high_cph_fraction"]
    threshold = payload["threshold"]
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
    if fraction != expected_fraction:
        raise ValueError("Bismark non-conversion summary fraction is inconsistent")
    if not 0 <= threshold <= 1:
        raise ValueError("Bismark non-conversion threshold must be within [0, 1]")
    if minimum_count < 1:
        raise ValueError("Bismark non-conversion minimum_count must be >= 1")
    if not isinstance(filter_applied, bool):
        raise ValueError("Bismark non-conversion filter_applied must be boolean")
    _validate_high_cph_route("bismark", payload["protocol"], payload["high_cph_role"])
    if not filter_applied:
        raise ValueError("Bismark non-conversion summary must record the fixed filter action")
    return payload




DROPLET_MIN_MAPQ = 10


DROPLET_MIN_BASEQ = 20

# Bismark PE 原生写入 FASTQ R1、R2。键为顺序 FLAG，值为共同 XG 和两端 XR。
# CTOT/CTOB 的 0x40/0x80 与物理端相反。
_BISMARK_PE_CONVERSIONS = {
    (99, 147): ("CT", ("CT", "GA")),  # OT
    (83, 163): ("GA", ("CT", "GA")),  # OB
    (147, 99): ("CT", ("GA", "CT")),  # CTOT
    (163, 83): ("GA", ("GA", "CT")),  # CTOB
}


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


def _droplet_consensus_table(sample: str, protocol: str, cpg: Mapping[str, int]) -> dict[str, Any]:
    """按当前协议和同一份 CpG 整数生成 MultiQC 表。"""
    policy = DD_CONSENSUS[protocol]
    return {
        "id": f"alopex_droplet_{policy.tag}_consensus",
        "section_name": f"{protocol} {policy.unit}",
        "description": f"逐胞嘧啶标签共识，标签来源 {policy.source}。M-bias 按共识前的 read-base 统计。",
        "plot_type": "table",
        "data": {sample: {
            "CpG_read_observations": cpg["read_observations"],
            policy.count_field: cpg["methyl"] + cpg["unmethyl"],
            policy.ambiguous_field: cpg[f"ambiguous_{policy.tag}s"],
            "CpG_strand_sites": cpg["strand_sites"],
        }},
    }


def call_droplet_cpg(bam: str | Path, reference: str | Path, coverage: str | Path,
                     qc: str | Path, mbias: str | Path, *, sample: str,
                     scratch: str | Path, samtools: str, threads: int = 1, protocol: str = "DD-MET5") -> None:
    """保留过滤 BAM 的全部 pairs，按参考胞嘧啶/UMI 生成分子共识 coverage、QC 及 read-level M-bias。"""
    import pysam

    if SAMPLE_ID_RE.fullmatch(sample) is None:
        raise ValueError("Invalid Droplet sample identity")
    policy = DD_CONSENSUS[protocol]
    bam_tag = "UR" if protocol == "DD-MET5" else "ZF"
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
                name = _require_sam_qname(first.query_name, bam)
                if second is None or name != second.query_name or name in seen:
                    raise ValueError("Droplet calling requires unique, adjacent complete read pairs")
                seen.add(name)
                conversions = _BISMARK_PE_CONVERSIONS.get((first.flag, second.flag))
                if conversions is None:
                    raise ValueError("Unexpected Bismark paired flags")
                conversion, read_conversions = conversions
                if (first.reference_id != second.reference_id
                        or any(record.next_reference_id != mate.reference_id
                               or record.next_reference_start != mate.reference_start
                               for record, mate in ((first, second), (second, first)))):
                    raise ValueError("Droplet mate coordinates disagree with the actual pair")
                read_id, _, umi = name.rpartition(":")
                if not read_id or not re.fullmatch(r"[ACGTN]{12}", umi):
                    raise ValueError(
                        f"{protocol} QNAME requires a read identifier followed by a colon and a 12 bp {policy.tag}"
                    )
                for mate, (record, xr) in enumerate(zip((first, second), read_conversions), 1):
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
                    if record.get_tag("XR") != xr or record.get_tag("XG") != conversion:
                        raise ValueError(f"Bismark XR/XG disagree with paired FLAG orientation: {name}")
                    if (record.query_sequence is None or record.query_qualities is None
                            or len(record.get_tag("XM")) != record.query_length):
                        raise ValueError("Droplet calling requires sequence, base qualities and matching XM length")
                    if record.has_tag(bam_tag) and record.get_tag(bam_tag) != umi:
                        raise ValueError(f"{bam_tag} disagrees with the {protocol} QNAME tag")
                    record.set_tag(bam_tag, umi, value_type="Z")
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
                          f"{policy.tag}_groups": 0, f"merged_{policy.tag}s": 0, f"ambiguous_{policy.tag}s": 0,
                          "strand_sites": 0} for key in ("CpG", "CHG", "CHH")}
        cycles: dict[tuple[str, str, int], list[int]] = {}
        command = [samtools, "mpileup", "-Q", str(DROPLET_MIN_BASEQ), "-q", str(DROPLET_MIN_MAPQ),
                   "-B", "-d", "0", "-f", str(reference), "--output-BP", "--output-extra", f"RLEN,{bam_tag},YI",
                   "--no-output-ins-mods", "--no-output-ins", "--no-output-ins",
                   "--no-output-del", "--no-output-del", "--no-output-ends", str(ordered)]
        with pysam.FastaFile(str(reference)) as fa, gzip.open(root / "coverage.gz", "wt") as handle:
            proc = subprocess.Popen(command, stdout=subprocess.PIPE, text=True, bufsize=1024 * 1024)
            reference_bases = _TapsReferenceWindow(fa)
            complement = str.maketrans("ACGT", "TGCA")
            for line in proc.stdout:
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 10:
                    raise ValueError("Unexpected samtools mpileup columns")
                if int(fields[3]) == 0:
                    continue
                base = fields[2].upper()
                if base not in {"C", "G"}:
                    continue
                chrom, pos = fields[0], int(fields[1]) - 1
                if not 0 <= pos < reference_lengths.get(chrom, 0):
                    raise ValueError("Droplet pileup coordinate disagrees with reference")
                if base == "C":
                    context_seq = reference_bases.fetch(chrom, pos, min(pos + 3, reference_lengths[chrom])).upper()
                else:
                    triplet = reference_bases.fetch(chrom, max(0, pos - 2), pos + 1)
                    context_seq = triplet.upper()[::-1].translate(complement)
                if not context_seq.startswith("C"):
                    raise ValueError("Droplet pileup coordinate disagrees with reference")
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
                stats[f"{policy.tag}_groups"] += len(groups)
                stats[f"merged_{policy.tag}s"] += merged
                stats[f"ambiguous_{policy.tag}s"] += ambiguous
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
            handle.write(f"Alopex {protocol} read-level M-bias; MAPQ>=10, baseQ>=20, before {policy.unit}\n\n")
            for context in contexts:
                for mate in ("1", "2"):
                    handle.write(f"{context} context (R{mate})\n====================\nposition\tcount_methylated\tcount_unmethylated\t% methylation\tcoverage\n")
                    for ctx, read, cycle in sorted(cycles):
                        if ctx == context and read == mate:
                            m, u = cycles[(ctx, read, cycle)]
                            handle.write(f"{cycle}\t{m}\t{u}\t{100 * m / (m + u):.8f}\t{m + u}\n")
                    handle.write("\n")
        payload = {"schema_version": 1, "sample": sample, "dedup_policy": policy.dedup_policy,
                   "min_mapq": DROPLET_MIN_MAPQ, "min_baseq": DROPLET_MIN_BASEQ,
                   "metric_unit": policy.unit, "protocol": protocol, "tag_source": policy.source, "processed_pairs": pairs,
                   f"{policy.tag}_with_n_pairs": umi_n, "contexts": contexts,
                   **_droplet_consensus_table(sample, protocol, contexts["CpG"])}
        (root / "qc.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
        read_droplet_umi_qc(root / "qc.json", sample, protocol)
        for name, target in (("coverage.gz", coverage), ("mbias.txt", mbias), ("qc.json", qc)):
            os.replace(root / name, target)
    finally:
        if proc is not None:
            _close_native_stream(proc)
        shutil.rmtree(root)


def read_droplet_umi_qc(path: str | Path, sample: str, protocol: str = "DD-MET5") -> dict[str, Any]:
    """校验逐胞嘧啶 UMI 共识计数、单位与实际参数，禁止作为 read-pair 去重率解释。"""
    policy = DD_CONSENSUS[protocol]
    data = _read_json(Path(path), "Droplet tag consensus QC")
    for key, value in {"schema_version": 1, "sample": sample, "dedup_policy": policy.dedup_policy,
                       "min_mapq": DROPLET_MIN_MAPQ, "min_baseq": DROPLET_MIN_BASEQ,
                       "metric_unit": policy.unit, "tag_source": policy.source, "protocol": protocol}.items():
        if type(data.get(key)) is not type(value) or data.get(key) != value:
            raise ValueError(f"Invalid Droplet UMI QC {key}: {path}")
    for key in ("processed_pairs", f"{policy.tag}_with_n_pairs"):
        data[key] = _integer(data.get(key), key, Path(path))
    if not 0 <= data[f"{policy.tag}_with_n_pairs"] <= data["processed_pairs"]:
        raise ValueError("Invalid Droplet UMI pair counts")
    contexts = data.get("contexts")
    if not isinstance(contexts, dict) or contexts.keys() != {"CpG", "CHG", "CHH"}:
        raise ValueError(f"Invalid Droplet consensus contexts: {path}")
    for context, values in contexts.items():
        if not isinstance(values, dict):
            raise ValueError(f"Invalid Droplet {context} counts: {path}")
        for key in ("methyl", "unmethyl", "read_observations", f"{policy.tag}_groups", f"merged_{policy.tag}s", f"ambiguous_{policy.tag}s", "strand_sites"):
            _integer(values.get(key), key, Path(path))
        called = values["methyl"] + values["unmethyl"]
        observations = values["read_observations"]
        groups = values[f"{policy.tag}_groups"]
        merged = values[f"merged_{policy.tag}s"]
        ambiguous = values[f"ambiguous_{policy.tag}s"]
        if ((not data["processed_pairs"] and observations) or bool(called) != bool(values["strand_sites"])
                or bool(observations) != bool(groups) or bool(groups) != bool(called + ambiguous)):
            raise ValueError(f"Droplet {context} observations disagree with input pairs, groups or called sites: {path}")
        # 未决代表至少有两条原始观测。合并也要求代表的原始丰度严格大于邻居。
        if (values["strand_sites"] > called or groups + ambiguous > observations
                or (merged and observations == groups) or groups != called + merged + ambiguous):
            raise ValueError(f"Droplet {context} consensus accounting does not close: {path}")
    table = _droplet_consensus_table(sample, protocol, contexts["CpG"])
    # Python 的 True 和 1.0 均等于 1，但都不是观测计数。
    if (any(data.get(key) != expected for key, expected in table.items())
            or any(type(count) is not int for count in data["data"][sample].values())):
        raise ValueError(f"Droplet MultiQC table disagrees with {sample} {protocol} CpG counts or units: {path}")
    return data


def build_metrics(
    *,
    sample: str,
    version: str,
    alignment_report: Path,
    trimmed_mates: tuple[str, str],
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
    umi_qc = read_droplet_umi_qc(dedup_report, sample, protocol) if protocol in DD_PROTOCOLS else None
    alignment, dedup = _bismark_pair_reports(
        alignment_report, None if umi_qc is not None else dedup_report, sample, trimmed_mates,
    )
    discarded_pairs = int(alignment["discarded"])
    accepted_pairs = dedup["total"]
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
    _validate_high_cph_route("bismark", protocol, high_cph_role)
    pairs = (candidate_pairs, removed_pairs, excluded_pairs, retained_pairs)
    if any(isinstance(n, bool) or not isinstance(n, int) or n < 0 for n in pairs):
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
            _high_cph_multiqc_payload(sample, payload),
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


@contextmanager
def _taps_collate(bam, samtools, scratch, *, primary_only=False):
    """隔离 samtools 的溢写文件；解析失败时先回收子进程，再清理本次目录。"""
    with tempfile.TemporaryDirectory(prefix=".taps_collate.", dir=scratch) as directory:
        command = [samtools, "collate"]
        if primary_only:
            command.append("-f")
        command.extend(["-O", "-u", "-T", str(Path(directory) / "reads"), bam])
        process = subprocess.Popen(command, stdout=subprocess.PIPE)
        try:
            yield process.stdout
            if process.wait() != 0:
                raise RuntimeError(f"samtools collate failed for TAPS BAM: {bam}")
        finally:
            _close_native_stream(process)


def _taps_read_groups(groups, sample):
    by_id = {group.get("ID"): group for group in groups}
    if (not groups or len(by_id) != len(groups) or not all(by_id)
            or any(group.get("SM") != sample for group in groups)):
        raise ValueError(f"TAPS BAM requires unique read groups bound to sample {sample}")
    return by_id


def _taps_primary_pair(records, name, bam, groups):
    _require_sam_qname(name, bam)
    primary = []
    template_group = None
    for read in records:
        try:
            group = read.get_tag("RG")
        except KeyError:
            group = None
        if group not in groups:
            raise ValueError(f"TAPS record has an unbound read group: {name}")
        if template_group is not None and group != template_group:
            raise ValueError(f"TAPS template spans different read groups: {name}")
        template_group = group
        if not read.is_secondary and not read.is_supplementary:
            primary.append(read)
    if (len(primary) != 2 or primary[0].is_read1 == primary[1].is_read1
            or any(not r.is_paired or r.is_read1 == r.is_read2 for r in primary)):
        raise ValueError(f"TAPS primary pair is incomplete or duplicated: {name}")
    return primary


def write_taps_alignment_qc(bam: str, samtools: str, scratch: str, sample: str, output: str) -> None:
    """按 QNAME 核算完整 BAM 与 Rastair 可用 pair；仅任一 mate 达到调用条件的 pair 进入漏斗。"""
    import itertools
    import pysam

    counts = dict.fromkeys(("bam_records", "bam_primary_records", "bam_duplicate_records",
                           "bam_unmapped_records", "bam_supplementary_records", "bam_secondary_records",
                           "trimmed_pairs", "backend_accepted_pairs", "backend_rejected_pairs",
                           "unmapped_pairs", "discordant_pairs", "quality_rejected_pairs",
                           "duplicate_pairs", "postdedup_pairs"), 0)
    require_distinct_paths([Path(bam), Path(f"{bam}.bai"), Path(output)], label="TAPS BAM, index and alignment QC")
    with pysam.AlignmentFile(bam, "rb") as source:
        header = source.header.to_dict()
        if header.get("HD", {}).get("SO") != "coordinate" or not source.has_index():
            raise ValueError("TAPS BAM must be coordinate sorted and indexed")
    with _taps_collate(bam, samtools, scratch) as stream:
        with pysam.AlignmentFile(stream, "rb") as source:
            collated = source.header.to_dict()
            groups = _taps_read_groups(collated.get("RG", []), sample)
            if groups != _taps_read_groups(header.get("RG", []), sample) or collated["SQ"] != header["SQ"]:
                raise ValueError("TAPS BAM header changed during alignment QC")
            for name, group in itertools.groupby(source, key=lambda r: r.query_name):
                records = list(group)
                primary = _taps_primary_pair(records, name, bam, groups)
                for r in records:
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


def taps_non_cpg_metrics(bam: str, reference: str, samtools: str, scratch: str, *, alignment_qc: Mapping, species: str) -> dict[str, object]:
    """统计可判链分子的宿主 CHG/CHH 转换率；未知或冲突链不随机分配，重叠 mate 每位置只计一次。

    这是未经 SNP/测序错误校正的原始修饰信号 QC；只用高质量 CpG 转换证据判链，
    指标只代表可判链子集，不能直接作为全基因组真实 mCH 或 cDNA 污染率。
    """
    import itertools
    import pysam

    counts = dict.fromkeys(("non_cpg_candidate_pairs", "non_cpg_oriented_pairs", "non_cpg_ambiguous_pairs",
                           "non_cpg_chg_mod", "non_cpg_chg_unmod", "non_cpg_chh_mod", "non_cpg_chh_unmod"), 0)
    expected_pairs = alignment_qc["trimmed_pairs"]
    if type(expected_pairs) is not int or expected_pairs < 0:
        raise ValueError("TAPS non-CpG QC requires a nonnegative complete pair count")
    host_contigs = set(primary_analysis_contigs(species))
    primary_pairs = 0
    with _taps_collate(bam, samtools, scratch, primary_only=True) as stream:
        with pysam.FastaFile(reference) as fasta, pysam.AlignmentFile(stream, "rb") as source:
            header = source.header.to_dict()
            sample = alignment_qc["sample"]
            groups = _taps_read_groups(header.get("RG", []), sample)
            if (groups != _taps_read_groups(alignment_qc["read_groups"], sample)
                    or header["SQ"] != alignment_qc["reference_sequences"]):
                raise ValueError(f"TAPS BAM header disagrees with alignment QC for {sample}")
            ref = _TapsReferenceWindow(fasta)
            for name, group in itertools.groupby(source, key=lambda r: r.query_name):
                primary = _taps_primary_pair(group, name, bam, groups)
                primary_pairs += 1
                chrom, bases, orientation, start, seq = _taps_pair_bases(primary, ref)
                if chrom not in host_contigs:
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
    if primary_pairs != expected_pairs:
        raise ValueError("TAPS non-CpG primary pair count disagrees with alignment QC")
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
    with tempfile.TemporaryDirectory(prefix=".taps_select.", dir=project_temporary_parent(output)) as selected:
        selected_vcf = str(Path(selected) / "variants.vcf")
        pysam.bcftools.view("--no-version", "--no-update", "--include", 'FILTER="PASS" && (GT="alt" || GT="mis")',
                            "--output-type", "v", "--output", selected_vcf, vcf, catch_stdout=False)
        with pysam.VariantFile(selected_vcf) as source, pysam.FastaFile(reference) as ref:
            if list(source.header.samples) != ["sample"]:
                raise ValueError("TAPS SNP VCF requires the native Rastair placeholder sample")
            reference_lengths = dict(zip(ref.references, ref.lengths))
            reference_bases = _TapsReferenceWindow(ref)
            header_text = str(source.header).rsplit("\t", 1)[0] + "\t" + sample + "\n"
            raw = str(Path(selected) / "selected.vcf.gz")
            with pysam.BGZFile(raw, "w") as target:
                target.write(header_text.encode())
                for row in source:
                    alleles = row.alleles or ()
                    if (len(alleles) < 2 or any(a not in {"A", "C", "G", "T"} for a in alleles)
                            or "GT" not in row.format):
                        continue
                    # pysam 的 GT 元组会把越界索引和真正缺失都解码成 None；按 VCF 中的数字校验。
                    record = str(row)
                    sample_data = record.rstrip("\n").rsplit("\t", 1)[1]
                    gt = sample_data.split(":")[row.format.keys().index("GT")]
                    called = [int(a) for a in re.split(r"[/|]", gt.lstrip("/|")) if a != "."]
                    if any(a < 0 or a >= len(alleles) for a in called):
                        raise ValueError(f"TAPS SNP GT refers to an undeclared allele: {row.contig}:{row.pos} {gt}")
                    if not any(called):
                        continue
                    if len(set(alleles)) != len(alleles):
                        raise ValueError(f"TAPS SNP has repeated REF/ALT alleles: {row.contig}:{row.pos}")
                    if (row.contig not in reference_lengths or row.start < 0
                            or row.stop > reference_lengths[row.contig]):
                        raise ValueError(f"TAPS SNP coordinate is outside reference FASTA: {row.contig}:{row.pos}")
                    if reference_bases.fetch(row.contig, row.start, row.stop).upper() != row.ref:
                        raise ValueError(f"TAPS SNP REF disagrees with reference FASTA: {row.contig}:{row.pos}")
                    target.write(record.encode())
                    count += 1
        corrected = str(Path(selected) / "snps.vcf.gz")
        pysam.bcftools.reheader("--fai", f"{reference}.fai", "--output", corrected, raw, catch_stdout=False)
        pysam.tabix_index(corrected, preset="vcf")
        index = Path(f"{output}.tbi")
        index.unlink(missing_ok=True)
        os.replace(corrected, output)
        os.replace(f"{corrected}.tbi", index)
    return count


def _emit_rastair_bed(bed: str, reference: str, stream: TextIO) -> dict:
    """写出共用 CpG writer 的整数输入，并统计原生行与两个对照的分母。"""
    import pysam

    counts = {"native_cpg_rows": 0, "canonical_cpg_rows": 0, "excluded_variant_or_denovo_rows": 0,
              "zero_observation_rows": 0}
    controls = {name: {"mod": 0, "unmod": 0, "sites": 0} for name in ("lambda", "pUC19")}
    required = {"#chr", "start", "end", "strand", "unmod", "mod", "no_snp", "snp", "coverage", "genotype", "cpg"}
    with pysam.FastaFile(reference) as ref, Path(bed).open() as handle:
        reader = csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
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
            if any(gt != base + "/" + base for _, _, _, _, gt, base in pending):
                counts["excluded_variant_or_denovo_rows"] += len(pending)
                return
            for chrom, pos, mod, unmod, gt, base in pending:
                if mod + unmod == 0:
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
            strand, genotype = row[column["strand"]], row[column["genotype"]]
            if strand not in {"+", "-", "."} or re.fullmatch(r"[ACGTN]/[ACGTN]", genotype) is None:
                raise ValueError(f"{label} has invalid Rastair strand or genotype: {strand!r}, {genotype!r}")
            chrom = row[column["#chr"]]
            if chrom not in contigs:
                raise ValueError(f"{label} contig {chrom!r} is absent from the reference FASTA")
            numeric = [row[column[k]] for k in ("start", "end", "mod", "unmod", "no_snp", "snp", "coverage")]
            if any(re.fullmatch(r"[0-9]+", value) is None for value in numeric):
                raise ValueError(f"{label} has non-integer coordinates/counts (requires nonnegative ASCII integers)")
            pos, end, mod, unmod, no_snp, snp, depth = map(int, numeric)
            if end != pos + 1 or end > ref.get_reference_length(chrom) or mod + unmod + no_snp + snp > depth:
                raise ValueError(f"{label}: Invalid Rastair coordinate/counts: {chrom}:{pos}")
            key = (contigs[chrom], pos)
            if previous is not None and key <= previous:
                raise ValueError(f"{label} is not ordered or contains duplicate positions")
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
                raise ValueError(f"{label}: Rastair reference CpG disagrees with FASTA: {chrom}:{pos}")
            if strand != ("+" if base == "C" else "-"):
                raise ValueError(f"{label}: Rastair strand disagrees with reference {base} at {chrom}:{pos}")
            cpg_key = (chrom, start)
            if cpg_key != pending_key:
                flush_cpg()
                pending = []
                pending_key = cpg_key
            pending.append((chrom, pos, mod, unmod, genotype, base))
        flush_cpg()
    for name, value in controls.items():
        denominator = value["mod"] + value["unmod"]
        counts[f"{name}_cpg_mod"] = value["mod"]
        counts[f"{name}_cpg_observations"] = denominator
        counts[f"{name}_cpg_sites"] = value["sites"]
        metric = "false_positive_pct" if name == "lambda" else "conversion_pct"
        counts[f"{name}_{metric}"] = 100.0 * value["mod"] / denominator if denominator else None
    return counts


def emit_rastair_cpg(bed: str, reference: str, alignment_qc: str, output_qc: str, stream: TextIO,
                     *, bam: str, samtools: str, scratch: str, vcf: str, output_snp: str, sample: str, species: str) -> None:
    """按当前 BAM 的细胞身份提交候选 SNP、索引和 QC，CpG 整数行交给共用 writer。"""
    snp_path, qc_path = Path(output_snp), Path(output_qc)
    index = Path(f"{output_snp}.tbi")
    locks = {path.parent.resolve() / ".alopex-rastair-call.lock" for path in (snp_path, qc_path)}
    members = [Path(bed), Path(vcf), Path(bam), Path(f"{bam}.bai"), Path(reference), Path(f"{reference}.fai"),
               Path(alignment_qc), snp_path, index, qc_path, *locks]
    path_label = "TAPS caller inputs and outputs"
    require_distinct_paths(members, label=path_label)
    qc = json.loads(Path(alignment_qc).read_text())
    if qc["sample"] != sample:
        raise ValueError("TAPS caller sample differs from alignment QC")
    qc.update(taps_non_cpg_metrics(bam, reference, samtools, scratch, alignment_qc=qc, species=species))
    qc.update({"signal": "5mC+5hmC", "high_cph_status": "not_applicable",
               "min_mapq": TAPS_MIN_MAPQ, "min_baseq": TAPS_MIN_BASEQ,
               "include_flags": TAPS_INCLUDE_FLAGS, "exclude_flags": TAPS_EXCLUDE_FLAGS,
               "guess_read_orientation": True, "mate_overlap_counted_once": True})
    for parent in {snp_path.parent, qc_path.parent}:
        parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".rastair_call.", dir=project_temporary_parent(snp_path)) as directory:
        snp = str(Path(directory) / "snps.vcf.gz")
        qc["snp_sites"] = emit_taps_snps(vcf, snp, sample, reference=reference)
        qc.update(_emit_rastair_bed(bed, reference, stream))
        stream.flush()
        # 原生文件和输出流均通过后才撤旧 QC；规则再将四件产物一起发布。
        with publication_locks(locks):
            require_distinct_paths(members, label=path_label)
            qc_path.unlink(missing_ok=True)
            os.replace(snp, snp_path)
            os.replace(f"{snp}.tbi", index)
            # 在忽略大小写的卷上，新建的成员可能同时占用了刚撤下的 QC 路径。
            if qc_path.exists():
                qc_path.unlink()
                raise ValueError(f"{path_label} must be distinct: {qc_path}")
            atomic_write_json(qc_path, qc)
