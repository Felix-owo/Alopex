"""资源预算、重试及执行器适配；不依赖科学处理、交付或 Doctor。"""
from __future__ import annotations

import sys

_LOADED_MODULE_CODE = sys._getframe().f_code

import math
import os
import re
from collections.abc import Iterator, Mapping
from dataclasses import dataclass, replace
from pathlib import Path
from types import MappingProxyType


def effective_cpu_count(requested: int = 0) -> int:
    """取宿主、进程 affinity、Slurm allocation 与调用方预算的最小正 CPU 数。"""
    limits = [max(1, os.cpu_count() or 1)]
    if hasattr(os, "sched_getaffinity"):
        try:
            limits.append(max(1, len(os.sched_getaffinity(0))))
        except OSError:
            pass
    process_count = os.process_cpu_count()
    if process_count:
        limits.append(process_count)
    for name in ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE"):
        value = os.environ.get(name, "")
        if value.isdecimal() and int(value) > 0:
            limits.append(int(value))
    if requested > 0:
        limits.append(int(requested))
    return min(limits)


def read_controller_progress(path: str | Path, offset: int = 0) -> tuple[int, str]:
    """有界读取新增进度；半行留待下次，文件缩短时从头读取，长日志行不整体驻留。"""
    progress = ""
    with Path(path).open("rb") as handle:
        size = os.fstat(handle.fileno()).st_size
        handle.seek(offset if 0 <= offset <= size else 0)
        offset = handle.tell()
        line_start, pending = True, ""
        while fragment := handle.readline(64 * 1024):
            if line_start:
                match = re.match(rb"^(\d+) of (\d+) steps \((\d+)%\)", fragment)
                pending = " ".join(part.decode("ascii") for part in match.groups()) if match else ""
            line_start = fragment.endswith(b"\n")
            if line_start:
                offset = handle.tell()
                if pending:
                    progress = pending
    return offset, progress


def disable_snakemake_output_mtime_gate() -> None:
    """禁用 Snakemake 跨节点 output/input mtime 硬失败，保留其余输出与完整性检查。"""

    from snakemake.dag import DAG

    current = getattr(DAG, "check_output_mtime", None)
    if current is None:
        raise RuntimeError("Snakemake DAG.check_output_mtime is unavailable")
    if getattr(current, "_dna_pipeline_nonfatal_output_mtime", False):
        return

    async def accept_output_mtime(_dag: object, _job: object, _outputs: object) -> None:
        return None

    setattr(accept_output_mtime, "_dna_pipeline_nonfatal_output_mtime", True)
    DAG.check_output_mtime = accept_output_mtime


def honor_slurm_no_account(workflow: object) -> None:
    """使空账户配置沿用集群默认值，不触发 executor 的账户猜测。"""
    if not getattr(getattr(workflow, "executor_settings", None), "no_account", False):
        return
    from snakemake_executor_plugin_slurm import Executor

    current = Executor.get_account_arg
    if getattr(current, "_alopex_no_account", False):
        return

    def account_arg(executor: object, job: object) -> Iterator[str]:
        # 当前插件虽接收 no_account，却仍进入猜测；其 sacct 管道还会误作普通 argv。
        # 每次查询当前 workflow，避免同进程的其他显式账户运行受影响。
        if getattr(executor.workflow.executor_settings, "no_account", False):
            yield ""
        else:
            yield from current(executor, job)

    setattr(account_arg, "_alopex_no_account", True)
    Executor.get_account_arg = account_arg


def preserve_executor_cli_paths() -> None:
    """需shell引用的绝对路径使用原生Path序列化，保留换行、反斜线与命令字符。"""
    import snakemake_interface_executor_plugins.utils as executor_utils

    current = executor_utils.format_cli_value
    if getattr(current, "_alopex_literal_paths", False):
        return

    def format_value(value, quote=False, base64_encode=False):
        # 路径已经过Alopex定位；只校正传输类型，不resolve或改变科学文件身份。
        if isinstance(value, str) and value.startswith("/") and quote and not base64_encode:
            value = Path(value)
        return current(value, quote=quote, base64_encode=base64_encode)

    setattr(format_value, "_alopex_literal_paths", True)
    executor_utils.format_cli_value = format_value


STATIC_POLICY = "static_reads"


_TIER_ORDER = ("S", "M", "L", "XL")


STATIC_TIER_MAX_READS = (1_000_000, 4_000_000, 6_000_000, None)


STATIC_TIER_MAX_DEMUX_BYTES = (
    5 * 1024 ** 3,
    50 * 1024 ** 3,
    300 * 1024 ** 3,
    None,
)


STATIC_MAX_ATTEMPTS = 3


STATIC_RETRY_MEMORY_RATIOS = ((3, 2), (2, 1))


BISMARK_CELL_RULE_RUNTIME_MIN = 150


_INTEGER_TEXT = re.compile(r"^[0-9]+$")


def _attempt_memory_mb(first_mem_mb: int, attempt_no: int) -> int:
    attempt_no = _attempt_number(attempt_no)
    memory_mb = int(first_mem_mb)
    if attempt_no == 1:
        return memory_mb
    numerator, denominator = STATIC_RETRY_MEMORY_RATIOS[attempt_no - 2]
    return (memory_mb * numerator + denominator - 1) // denominator


@dataclass(frozen=True)
class ResourceRequest:
    """静态策略选出的一笔不可变调度请求。"""

    policy: str
    rule_key: str
    tier: str
    dna_reads: int | None
    threads: int
    mem_mb: int
    runtime_min: int


@dataclass(frozen=True)
class _StaticRuleRequest:
    threads: tuple[int, ...]
    memory_gib: tuple[int, ...]
    runtime_min: tuple[int, ...]

    def __post_init__(self) -> None:
        lengths = {len(self.threads), len(self.memory_gib), len(self.runtime_min)}
        if len(lengths) != 1 or len(self.threads) != 4:
            raise ValueError("static rule requests must define S/M/L/XL resource tiers")


STATIC_RESOURCE_REQUESTS: Mapping[str, _StaticRuleRequest] = MappingProxyType(
    {
        "demux": _StaticRuleRequest((4, 4, 4, 4), (4, 8, 12, 12), (120, 360, 720, 1440)),
        "droplet_calling": _StaticRuleRequest((32, 32, 32, 32), (16, 16, 16, 16), (1440, 1440, 1440, 1440)),
        "droplet_demux": _StaticRuleRequest((32, 32, 32, 32), (16, 16, 16, 16), (1440, 1440, 1440, 1440)),
        "cutadapt": _StaticRuleRequest((4, 4, 4, 4), (2, 2, 2, 4), (30, 45, 60, 90)),
        "align_sort_dedup": _StaticRuleRequest(
            (16, 16, 16, 16), (24, 32, 40, 40), (120, 180, 240, 240)
        ),


        "rastair_align_mark": _StaticRuleRequest(
            (16, 16, 16, 16), (16, 16, 16, 16), (120, 180, 240, 240)
        ),
        "rastair_call": _StaticRuleRequest(
            (16, 16, 16, 16), (4, 8, 12, 16), (60, 120, 180, 240)
        ),
        "bismark_align_dedup": _StaticRuleRequest(
            (8, 12, 12, 16), (16, 16, 16, 16), (120, 180, 270, 900)
        ),
        "bismark_align_dedup_local": _StaticRuleRequest(
            (12, 12, 16, 16), (24, 24, 32, 32), (120, 180, 270, 900)
        ),
        "high_cph_read_names": _StaticRuleRequest(
            (2, 2, 2, 2), (4, 8, 12, 16), (120, 180, 240, 360)
        ),
        "bsconv_dna": _StaticRuleRequest((2, 2, 2, 2), (4, 8, 12, 16), (120, 180, 240, 360)),
        "high_cph_bam": _StaticRuleRequest(
            (2, 2, 2, 2), (4, 8, 12, 16), (120, 180, 240, 360)
        ),

        "pileup": _StaticRuleRequest((1, 1, 1, 1), (4, 8, 12, 16), (60, 120, 180, 240)),
        "generate_cpg_input": _StaticRuleRequest(
            (1, 1, 1, 1), (2, 4, 4, 8), (60, 120, 180, 240)
        ),


        "bismark_extract": _StaticRuleRequest(
            (8, 8, 8, 8), (16, 16, 16, 16), (60, 90, 150, 240)
        ),
        "generate_snp_input": _StaticRuleRequest(
            (1, 1, 1, 1), (2, 4, 4, 8), (60, 120, 180, 240)
        ),
        "biscuit_qc": _StaticRuleRequest((1, 1, 1, 1), (2, 4, 4, 8), (120, 180, 240, 360)),
    }
)


STATIC_SCRATCH_GIB: Mapping[str, tuple[int, int, int, int]] = MappingProxyType(
    {
        "preprocess_alignment": (16, 32, 64, 64),
        "high_cph_filter": (8, 16, 24, 48),
        "bismark_extract": (8, 16, 24, 48),
    }
)


STATIC_FUSED_RULES: Mapping[str, tuple[str, ...]] = MappingProxyType(
    {
        "high_cph_filter": ("high_cph_read_names", "bsconv_dna"),
    }
)


CONTROLLER_RESOURCE_REQUESTS: Mapping[str, tuple[int, int, int]] = MappingProxyType(
    {
        "initial_cell_manifest": (1, 512, 10),
        "final_sample_manifest": (8, 2048, 30),
        "delivery_ready": (1, 1024, 15),
    }
)


def _read_count(value: object) -> int | None:

    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value if value >= 0 else None
    if isinstance(value, str):
        text = value.strip()
        if _INTEGER_TEXT.fullmatch(text):
            return int(text)
    return None


def size_tier(dna_reads: object) -> str:
    """按 DNA read pairs 判定 S/M/L/XL tier；缺失或非法值为 XL。"""

    reads = _read_count(dna_reads)
    if reads is None:
        return "XL"
    if reads <= STATIC_TIER_MAX_READS[0]:
        return "S"
    if reads <= STATIC_TIER_MAX_READS[1]:
        return "M"
    if reads <= STATIC_TIER_MAX_READS[2]:
        return "L"
    return "XL"


def demux_size_tier(total_bytes: object) -> str:
    """按压缩 FASTQ 总字节判定 demux 的 S/M/L/XL tier。

    ``None``（无 run snapshot 或文件总大小不可得）保守回退 XL。
    """

    total = _read_count(total_bytes)
    if total is None:
        return "XL"
    if total <= STATIC_TIER_MAX_DEMUX_BYTES[0]:
        return "S"
    if total <= STATIC_TIER_MAX_DEMUX_BYTES[1]:
        return "M"
    if total <= STATIC_TIER_MAX_DEMUX_BYTES[2]:
        return "L"
    return "XL"


def demux_resource_request(total_bytes: object, attempt: int = 1) -> ResourceRequest:
    """返回孔板 demux 请求；重试仅提高内存，Droplet 各阶段由其独立入口预算。"""

    attempt_no = _attempt_number(attempt)
    key, row = _rule("demux")
    tier = demux_size_tier(total_bytes)
    tier_index = _tier_index(tier)
    mem_mb = _attempt_memory_mb(row.memory_gib[tier_index] * 1024, attempt_no)
    return ResourceRequest(
        policy=STATIC_POLICY,
        rule_key=key,
        tier=tier,
        dna_reads=None,
        threads=row.threads[tier_index],
        mem_mb=mem_mb,
        runtime_min=row.runtime_min[tier_index],
    )


def droplet_resource_request(total_r1_bytes: object, stage: str, attempt: int = 1) -> ResourceRequest:
    """按紧凑计数的大库斜率与配对缓冲请求 Droplet 单库资源，保留两次内存重试。"""
    if stage not in {"calling", "demux"}:
        raise ValueError("Droplet stage must be calling or demux")
    request = resource_request("droplet_" + stage, attempt=1)
    tier = demux_size_tier(total_r1_bytes)
    base_mem = request.mem_mb
    if stage == "calling":
        total = _read_count(total_r1_bytes)
        if tier == "XL" and total is not None:
            # 保留 F-real-2609SG 全库计数实测斜率、余量及紧凑表换算，逐层向上取整。
            count_mem_mb = math.ceil((2 + 0.048 * (total / 1024 ** 3)) * 1.3 * 1024)
        else:
            count_mem_mb = STATIC_RESOURCE_REQUESTS["demux"].memory_gib[_tier_index(tier)] * 1024
        base_mem = max(base_mem, math.ceil(count_mem_mb / 4) + 4 * 1024)
    return replace(request, tier=tier, mem_mb=_attempt_memory_mb(base_mem, _attempt_number(attempt)))


MULTIQC_MAX_CELLS = 1000


def multiqc_resource_request(source_count: object, attempt: int = 1) -> ResourceRequest:
    """按全量声明依赖数为报告及 Snakemake worker 分配资源，重试仅提高内存。"""

    attempt_no = _attempt_number(attempt)
    count = _read_count(source_count) or 0
    base_mem_mb = math.ceil((512 + 0.075 * count) * 1.5)
    mem_mb = _attempt_memory_mb(max(base_mem_mb, 1024), attempt_no)
    return ResourceRequest(
        policy=STATIC_POLICY + ":multiqc_sources",
        rule_key="multiqc",
        tier="NA",
        dna_reads=None,
        threads=1,
        mem_mb=mem_mb,
        runtime_min=max(60, math.ceil(count / 600)),
    )


def controller_resource_request(rule_key: str, attempt: int = 1) -> ResourceRequest:
    """返回控制器类小任务的固定调度请求（键不存在即硬错误）。"""

    attempt_no = _attempt_number(attempt)
    try:
        threads, mem_mb, runtime_min = CONTROLLER_RESOURCE_REQUESTS[rule_key]
    except KeyError as exc:
        supported = ", ".join(CONTROLLER_RESOURCE_REQUESTS)
        raise ValueError(
            f"unsupported controller rule key {rule_key!r}; expected one of: {supported}"
        ) from exc
    mem = _attempt_memory_mb(int(mem_mb), attempt_no)
    return ResourceRequest(
        policy=STATIC_POLICY + ":controller",
        rule_key=str(rule_key),
        tier="NA",
        dna_reads=None,
        threads=int(threads),
        mem_mb=mem,
        runtime_min=int(runtime_min),
    )


def _attempt_number(attempt: object) -> int:
    if (
        isinstance(attempt, bool)
        or not isinstance(attempt, int)
        or not 1 <= attempt <= STATIC_MAX_ATTEMPTS
    ):
        raise ValueError("attempt must be 1, 2 or 3")
    return attempt


def _rule(rule_key: object) -> tuple[str, _StaticRuleRequest]:
    key = str(rule_key).strip()
    try:
        return key, STATIC_RESOURCE_REQUESTS[key]
    except KeyError as exc:
        supported = ", ".join(STATIC_RESOURCE_REQUESTS)
        raise ValueError(
            f"unsupported resource rule key {key!r}; expected one of: {supported}"
        ) from exc


def _tier_index(tier: str) -> int:
    return _TIER_ORDER.index(tier)


def resource_request(
    rule_key: str, dna_reads: object = None, attempt: int = 1
) -> ResourceRequest:
    """返回一条 rule/cell 的确定性调度请求。"""

    attempt_no = _attempt_number(attempt)
    key, row = _rule(rule_key)
    tier = size_tier(dna_reads)
    tier_index = _tier_index(tier)
    mem_mb = _attempt_memory_mb(row.memory_gib[tier_index] * 1024, attempt_no)
    return ResourceRequest(
        policy=STATIC_POLICY,
        rule_key=key,
        tier=tier,
        dna_reads=_read_count(dna_reads),
        threads=row.threads[tier_index],
        mem_mb=mem_mb,
        runtime_min=row.runtime_min[tier_index],
    )


def resource_threads(rule_key: str, dna_reads: object = None) -> int:
    """返回 ``rule_key`` 按 tier 选定的确定性线程数。"""

    _, row = _rule(rule_key)
    tier = size_tier(dna_reads)
    return row.threads[_tier_index(tier)]


def resource_runtime_min(
    rule_key: str, dna_reads: object, attempt: int = 1
) -> int:
    """返回 walltime 分钟数；retry 不改变 walltime。"""

    return resource_request(rule_key, dna_reads, attempt).runtime_min


def combined_runtime_min(
    rule_keys: tuple[str, ...], dna_reads: object, attempt: int = 1
) -> int:
    """合并同一 fused job 内多个静态 stage 的时限。"""

    _attempt_number(attempt)
    if not rule_keys:
        raise ValueError("rule_keys must contain at least one resource rule")
    return sum(resource_runtime_min(key, dna_reads, attempt) for key in rule_keys)


_BACKEND_STAGE_RULES: Mapping[tuple[str, str], str] = MappingProxyType(
    {
        ("biscuit", "align"): "align_sort_dedup",
        ("rastair", "align"): "rastair_align_mark",
        ("rastair", "call"): "rastair_call",
        ("bismark", "align"): "bismark_align_dedup",
        ("bismark", "extract"): "bismark_extract",
    }
)


_LOCAL_BACKEND_THREADS: Mapping[str, tuple[int, int, int, int]] = MappingProxyType(
    {
        "align_sort_dedup": (8, 12, 16, 16),
        "rastair_align_mark": (8, 12, 16, 16),
        "rastair_call": (4, 4, 4, 4),
        "bismark_align_dedup": (8, 12, 16, 16),
        "bismark_align_dedup_local": (12, 12, 16, 16),
        "bismark_extract": (8, 8, 8, 8),
    }
)


def backend_resource_request(
    backend: str,
    stage: str,
    dna_reads: object = None,
    *,
    attempt: int = 1,
    executor: str = "slurm",
    local_cores: object = None,
    bismark_local_alignment: bool = False,
) -> ResourceRequest:
    """返回结合 backend/executor/sample tier 的调度请求。

    memory/runtime 始终绑定 backend 的读数 tier；只有 local CPU 分配会再被配置的工作站预算截断。
    """

    backend_key = str(backend).strip().lower()
    stage_key = str(stage).strip().lower()
    pair = (backend_key, stage_key)
    try:
        rule_key = _BACKEND_STAGE_RULES[pair]
    except KeyError as exc:
        supported = ", ".join(f"{b}/{s}" for b, s in _BACKEND_STAGE_RULES)
        raise ValueError(
            f"unsupported backend resource stage {backend_key!r}/{stage_key!r}; "
            f"expected one of: {supported}"
        ) from exc

    if pair == ("bismark", "align") and bool(bismark_local_alignment):
        rule_key = "bismark_align_dedup_local"

    base = resource_request(rule_key, dna_reads, attempt=attempt)
    mode = str(executor or "slurm").strip().lower()
    if mode == "slurm":
        return base
    if mode != "local":
        raise ValueError("executor must be 'local' or 'slurm'")

    cores = _read_count(local_cores)
    if cores is None or cores < 1:
        return base
    preferred = _LOCAL_BACKEND_THREADS[rule_key][_TIER_ORDER.index(base.tier)]
    local_job_cap = cores if pair == ("rastair", "call") else max(1, min(16, cores // 2 if cores > 1 else 1))
    threads = max(1, min(preferred, local_job_cap, cores))
    return ResourceRequest(
        policy=f"{STATIC_POLICY}:{backend_key}:{mode}",
        rule_key=base.rule_key,
        tier=base.tier,
        dna_reads=base.dna_reads,
        threads=threads,
        mem_mb=base.mem_mb,
        runtime_min=base.runtime_min,
    )


def scratch_mb(scratch_group: str, dna_reads: object) -> int:
    """返回 fused job 固定的项目 scratch 容量（MiB）。"""

    key = str(scratch_group).strip()
    try:
        values = STATIC_SCRATCH_GIB[key]
    except KeyError as exc:
        supported = ", ".join(STATIC_SCRATCH_GIB)
        raise ValueError(
            f"unsupported scratch resource group {key!r}; expected one of: {supported}"
        ) from exc
    return int(values[_TIER_ORDER.index(size_tier(dna_reads))] * 1024)
