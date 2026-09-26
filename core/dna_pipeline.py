"""Alopex 核心 Python 单模块：config 唯一加载链、协议/backend 策略、
manifest 家族、run snapshot 与发布事务、reference/FASTQ 身份、Droplet 全库计数与分块解复用、静态资源与路径契约、
Snakemake bootstrap mtime 策略、环境身份与 Conda release 状态机，以及统一内部 CLI。

Shell / Snakefile 经 ``python -m dna_pipeline <子命令>`` 调用；子命令按
项目/运行状态/阶段写出/Doctor benchmark/环境/reference 分组。
"""
from __future__ import annotations

import argparse
import bisect
import copy
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import logging
import math
import os
import re
import shlex
import shutil
import struct
import subprocess
import sys
import tempfile
import uuid
import warnings
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from types import MappingProxyType
from typing import Any, TextIO, TypedDict, cast

DROPLET_DESIGN_RESOURCE = "resources/droplet_ME5_U3CB_methylation.txt.gz"
DROPLET_DESIGN_SHA256 = "67981bee6acda8db275257a54e727d1e0a0bf4b884c10d8e761d04eaa9b53588"


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


def demux_source_revision(source_root: str | Path, *, count_only: bool = False) -> str:
    """计算 Rust 身份；calling 身份排除不参与计数和分子复核的配对输出模块。"""
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
    if rust_sources != [main_rs]:
        raise ValueError(
            "demux_rs source tree must contain only src/main.rs: "
            + ", ".join(str(path) for path in rust_sources)
        )
    digest = hashlib.sha256()
    for path in sorted(paths, key=lambda item: Path(os.path.relpath(item, root)).as_posix()):
        relative = Path(os.path.relpath(path, root)).as_posix().encode()
        data = path.read_bytes()
        if count_only and path == main_rs:
            before, output_marker, output_and_rest = data.partition(b"\nmod output {\n")
            _, droplet_marker, after = output_and_rest.partition(b"\nmod droplet {\n")
            if not output_marker or not droplet_marker:
                raise ValueError("cannot locate Rust output and droplet module boundaries")
            data = before + droplet_marker + after
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
        if (before.st_size, before.st_mtime_ns) != (
            after.st_size,
            after.st_mtime_ns,
        ):
            raise RuntimeError(f"file changed while recording identity: {source}")
    return identity


def atomic_write_bytes(path: str | Path, data: bytes, *, mode: int | None = None) -> None:
    """在目标路径旁写临时文件，再以一次原子 rename 发布。"""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp.{os.getpid()}")
    try:
        with temporary.open("wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        if mode is not None:
            os.chmod(temporary, mode)
        os.replace(temporary, destination)
    finally:
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


ENVIRONMENT_IDENTITY_SCHEMA = 8
IDENTITY_NAME = ".dna_pipeline_identity.json"
ENVS_SPEC_NAME = "envs.yaml"


def envs_section(envs_path: str | Path, role: str) -> dict[str, object]:
    """从五角色 envs.yaml 读取指定角色的 conda specification section。"""
    import yaml as _yaml

    payload = _yaml.safe_load(Path(envs_path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or payload.get(role) is None:
        raise ValueError(f"envs.yaml has no section for role {role!r}: {envs_path}")
    section = payload[role]
    if not isinstance(section, dict):
        raise ValueError(f"envs.yaml section {role!r} must be a mapping")
    return section


def envs_section_sha256(envs_path: str | Path, role: str) -> str:
    """角色 section 的 canonical SHA256（环境身份与复用判定的唯一口径）。"""
    section = envs_section(envs_path, role)
    return sha256_bytes(canonical_json_bytes(section))


def envs_spec_logical_name(role: str) -> str:
    """返回环境 spec 的逻辑名，与环境安装位置无关。"""
    return f"{ENVS_SPEC_NAME}#{role}"


def write_envs_section_spec(envs_path: str | Path, role: str, destination: str | Path) -> Path:
    """把角色 section 抽取为独立 conda spec 文件（conda env create 的输入）。"""
    import yaml as _yaml

    section = envs_section(envs_path, role)
    target = Path(destination)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        _yaml.safe_dump(section, sort_keys=False, allow_unicode=True), encoding="utf-8"
    )
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


def read_environment_identity(prefix: str | Path) -> Mapping[str, object]:
    """读取并校验一笔已记录的受管环境身份。"""

    prefix_path = Path(prefix).expanduser().resolve()
    identity_path = prefix_path / IDENTITY_NAME
    try:
        payload = json.loads(identity_path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"environment record is missing: {identity_path}") from exc
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid environment record JSON: {identity_path}: {exc}") from exc
    if not isinstance(payload, dict) or payload.get("schema_version") != ENVIRONMENT_IDENTITY_SCHEMA:
        raise ValueError(f"unsupported environment record: {identity_path}")
    if payload.get("complete") is not True:
        raise ValueError(f"environment record is incomplete: {identity_path}")
    if Path(str(payload.get("prefix", ""))).resolve() != prefix_path:
        raise ValueError(f"environment record prefix mismatch: {identity_path}")
    if not str(payload.get("role", "")).strip():
        raise ValueError(f"environment record role is missing: {identity_path}")
    specification = payload.get("specification")
    if not isinstance(specification, Mapping) or not str(specification.get("sha256", "")).strip():
        raise ValueError(f"environment specification identity is missing: {identity_path}")
    return payload


def canonical_json_bytes(value: object) -> bytes:
    """用于稳定内容身份的 canonical JSON 编码（排序键、紧凑分隔符）。"""
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


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
    temporary = current.with_name(f".{current.name}.publish.{os.getpid()}")
    try:
        temporary.unlink(missing_ok=True)
        os.symlink(str(target), temporary)
        os.replace(temporary, current)
    finally:
        temporary.unlink(missing_ok=True)


def validate_release(current_s: str, expected_subdir: str, envs_path: str) -> None:
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
        recorded_resolved = str(record.get("resolved_prefix", "")).strip()
        if recorded_resolved and Path(recorded_resolved).resolve() != prefix:
            raise ValueError(f"{role} resolved prefix drift")
        environment = read_environment_identity(prefix)
        if environment.get("role") != role:
            raise ValueError(f"{role} environment identity role mismatch")
        recorded = (environment.get("specification") or {}).get("sha256")
        if recorded != envs_section_sha256(envs_path, role):
            raise ValueError(f"{role} specification changed; rebuild required")


def environment_is_reusable(prefix_s: str, role: str, envs_path: str) -> None:
    """判断既有 prefix 可否复用：身份记录存在、role 匹配且 envs section SHA 未变。"""
    environment = read_environment_identity(prefix_s)
    if environment.get("role") != role:
        raise ValueError(f"{role} environment identity role mismatch")
    if (environment.get("specification") or {}).get("sha256") != envs_section_sha256(
        envs_path, role
    ):
        raise ValueError(f"{role} specification changed")


CONFIG_SCHEMA_YAML = """$schema: "https://json-schema.org/draft/2020-12/schema"
type: object
additionalProperties: false
required:
  - schema_version
  - species
  - references
  - analysis
  - demux
  - biscuit
  - bismark
  - runtime
  - slurm
  - retention
properties:
  schema_version:
    description: 配置 schema 版本号，当前固定为 2。
    type: integer
    const: 2

  species:
    description: 当前项目的物种标识，决定默认 reference 选择与结果目录命名。
    type: string
    minLength: 1
    pattern: '^[A-Za-z0-9_.-]+$'

  references:
    description: 物种 -> 标准 reference FASTA 路径映射；由 core/doctor.sh reference 构建与校验。
    type: object
    minProperties: 1
    propertyNames:
      pattern: '^[A-Za-z0-9_.-]+$'
    additionalProperties:
      type: string
      minLength: 1

  analysis:
    description: 实验协议与甲基化分析后端的全局选择。
    type: object
    additionalProperties: false
    required: [protocol, methylation_backend]
    properties:
      protocol:
        description: 建库协议（srd / cabernet / taps / droplet）；manifest 的 rna_sample 在 srd 为本地 RNA tube，在 cabernet/taps 为可选外部 RNA raw 关联，在 droplet 必须留空。
        type: string
        enum: [srd, cabernet, taps, droplet]
      methylation_backend:
        description: 甲基化分析后端（biscuit / bismark / rastair），决定 alignment 与 CpG 生产路线。
        type: string
        enum: [biscuit, bismark, rastair]


  demux:
    description: Rust 解复用器（demux_rs）参数；省略时使用内置默认。
    type: object
    additionalProperties: false
    required: [min_matched_read_pairs, dna_w_spacer_len]
    properties:
      min_matched_read_pairs:
        description: 单 cell barcode 匹配 read pairs 的最低门槛，低于该值的 cell 不产出结果。
        type: integer
        minimum: 0
      dna_w_spacer_len:
        description: DNA 读段布局中 barcode 与 ME motif 间 W spacer 碱基数，须与实际建库一致。
        type: integer
        minimum: 0
        maximum: 20



  high_cph:
    description: Cabernet/SRD/Droplet 的 high-CpH / non-conversion 处理参数，固定统计并剔除；TAPS 不适用，无用户开关；省略时使用内置默认。
    type: object
    additionalProperties: false
    required:
      - excluded_contigs
    properties:
      excluded_contigs:
        description: 不参与 conversion 判断的 control/线粒体 contig 名称列表。
        type: array
        minItems: 0
        uniqueItems: true
        items:
          type: string
          pattern: '^[^,\\s]+$'

  biscuit:
    description: BISCUIT 后端专属参数；仅 methylation_backend=biscuit 时生效。
    type: object
    additionalProperties: false
    required:
      - library_type
      - high_cph_retention_threshold
      - generate_snp
    properties:
      library_type:
        description: BISCUIT 建库方向模式，决定允许的 strand model。
        type: string
        enum: [directional, non_directional]
      high_cph_retention_threshold:
        description: read pair CpH retention 严格大于该值时判为 high-CpH；越低越严格。
        type: number
        minimum: 0
        maximum: 1
      generate_snp:
        description: BISCUIT-only SNP 输出开关；非 BISCUIT 后端必须保持 false。
        type: boolean

  bismark:
    description: Bismark 后端专属参数；仅 methylation_backend=bismark 时生效。
    type: object
    additionalProperties: false
    required:
      - library_type
      - local_alignment
      - non_conversion_percentage_cutoff
      - non_conversion_minimum_count
    properties:
      library_type:
        description: Bismark 建库方向模式，决定搜索的 bisulfite strand space。
        type: string
        enum: [directional, non_directional, pbat]
      local_alignment:
        description: false 用 combined-index end-to-end（默认）、true 用 faithful CT/GA + --local；两模式互斥。
        type: boolean
      non_conversion_percentage_cutoff:
        description: 任一 mate 的 non-CG retention 达到该百分比且满足 minimum_count 时移除整对；越低越严格。
        type: integer
        minimum: 0
        maximum: 100
      non_conversion_minimum_count:
        description: 启用百分比判断所需的最低 informative non-CG C 数量。
        type: integer
        minimum: 1

  runtime:
    description: 与生物学结果无关的执行器运行时容错参数。
    type: object
    additionalProperties: false
    required: [local_executor_cores, latency_wait_seconds]
    properties:
      local_executor_cores:
        description: local executor 的总 CPU 数（snakemake --cores），只影响吞吐不改变结果。
        type: integer
        minimum: 1
      latency_wait_seconds:
        description: snakemake --latency-wait 值，等待输出文件在文件系统可见的秒数。
        type: integer
        minimum: 0

  slurm:
    description: Slurm executor 调度参数；仅 Linux/HPC 使用。
    type: object
    additionalProperties: false
    required: [jobs, controller_cores, account, partition, qos, node_tmpdir]
    properties:
      jobs:
        description: 同时占用的最大 SLURM 作业数，含 Droplet 全库 calling、demux 块和后续 cell 作业。
        type: integer
        minimum: 1
      controller_cores:
        description: Snakemake controller 本地进程核数（--local-cores）。
        type: integer
        minimum: 1
      account:
        description: SLURM 记账账户；null 表示不指定账户。
        oneOf:
          - {type: string, minLength: 1}
          - {type: 'null'}
      partition:
        description: 提交作业的默认分区；null 沿用集群默认分区。
        oneOf:
          - {type: string, minLength: 1}
          - {type: 'null'}
      qos:
        description: 提交作业的 QoS；字符串随作业提交，null 沿用集群默认 QoS。
        oneOf:
          - {type: string, minLength: 1}
          - {type: 'null'}
      node_tmpdir:
        description: 计算节点本地 scratch 目录；null 表示优先使用作业自带的 $SLURM_TMPDIR。
        oneOf:
          - type: string
            pattern: '^/.*'
            not: {const: '/'}
            allOf:
              - not: {pattern: '(^|/)\\.\\.?(/|$)'}
          - type: 'null'

  retention:
    description: 恢复边界策略；控制昂贵中间文件的保留面。
    type: object
    additionalProperties: false
    required: [keep_final_bam]
    properties:
      keep_final_bam:
        description: true 时保留 02_work/ 最终 BAM；TAPS 保留全部记录的 marked BAM+BAI 供变异重分析，其余后端保留过滤后 DNA BAM；默认 false。
        type: boolean

allOf:
  - if:
      properties:
        analysis:
          properties:
            methylation_backend: {enum: [bismark, rastair]}
          required: [methylation_backend]
        biscuit:
          properties:
            generate_snp: {const: true}
          required: [generate_snp]
      required: [analysis, biscuit]
    then: false
  - if:
      properties:
        analysis:
          properties:
            protocol: {const: taps}
    then:
      properties:
        analysis:
          properties:
            methylation_backend: {const: rastair}
  - if:
      properties:
        analysis:
          properties:
            methylation_backend: {const: rastair}
    then:
      properties:
        analysis:
          properties:
            protocol: {const: taps}
"""


CONFIG_TEMPLATE_YAML = """schema_version: 2

species: hg38                    # 标准 hg38/mm10 reference 由 Doctor 构建

references:
  hg38: resources/hg38_reference/hg38.primary.lambda.puc19.fa  # primary genome + lambda + pUC19；相对 Pipeline 根目录
  mm10: resources/mm10_reference/mm10.primary.lambda.puc19.fa  # primary genome + lambda + pUC19；也接受绝对路径

analysis:
  protocol: cabernet              # cabernet / srd / taps / droplet；manifest 的 rna_sample 在 srd 为本地输入，cabernet/taps 为外部 RNA raw 关联，droplet 必须留空
  methylation_backend: bismark    # cabernet/srd 选 biscuit/bismark；droplet 仅 bismark；taps 仅 rastair

demux:
  min_matched_read_pairs: 10   # 单 cell barcode 匹配 read pairs 低于该值时不产出结果；0 关闭门槛
  dna_w_spacer_len: 0          # DNA 建库布局中 barcode 与 ME motif 间 W spacer 碱基数，须与实际建库一致

biscuit:
  library_type: non_directional   # directional / non_directional
  high_cph_retention_threshold: 0.7   # read pair CpH retention 严格大于该值判为 high-CpH；越低越严格
  generate_snp: false             # BISCUIT-only SNP 输出开关；true 时额外发布 snps.bed.gz，非 BISCUIT 后端禁止开启

bismark:
  library_type: non_directional   # directional / non_directional / pbat
  local_alignment: false          # false: combined-index end-to-end（默认）；true: faithful CT/GA + --local
  non_conversion_percentage_cutoff: 70    # 任一 mate 的 non-CG retention 达标且满足 minimum_count 时移除整对；越低越严格
  non_conversion_minimum_count: 5         # 至少这么多个 informative non-CG C 才启用百分比判断，防短 CpH 误杀

high_cph:
  excluded_contigs: [pUC19, lambda, chrM]  # 这些对照不参与 high-CpH 筛除；名称须匹配参考 contig，TAPS 不适用

runtime:
  local_executor_cores: 32        # local executor 的总 CPU 数（snakemake --cores）；影响吞吐不改变结果
  latency_wait_seconds: 120       # snakemake --latency-wait：等待输出文件在文件系统可见的秒数

retention:
  keep_final_bam: false          # 保留 02_work/ 最终 BAM；TAPS 为全部记录的 marked BAM+BAI，供后续 SNP/SNV

slurm:
  jobs: 30                        # SLURM 全局并发上限，含 Droplet 全库 calling、demux 块和后续 cell 作业
  controller_cores: 2             # Snakemake controller 本地进程的 CPU 数（--local-cores）
  account: null                  # 仅在集群要求时填写账户
  partition: null                 # 提交作业的默认分区；null 沿用集群默认分区
  qos: huge                       # 提交作业的默认 QoS；设为 null 时沿用集群默认 QoS
  node_tmpdir: /tmp               # 优先 SLURM_TMPDIR；此路径不可写时回退项目 05_tmp/
"""


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
        "demux_index": _StaticRuleRequest((32, 32, 32, 32), (16, 16, 16, 16), (1440, 1440, 1440, 1440)),
        "demux_chunk": _StaticRuleRequest((6, 6, 6, 6), (8, 8, 8, 8), (120, 120, 120, 120)),
        "demux_merge": _StaticRuleRequest((8, 8, 8, 8), (16, 16, 16, 16), (1440, 1440, 1440, 1440)),
        "cutadapt": _StaticRuleRequest((4, 4, 4, 4), (2, 2, 2, 4), (30, 45, 60, 90)),
        "align_sort_dedup": _StaticRuleRequest(
            (16, 16, 16, 16), (24, 32, 40, 40), (120, 180, 240, 240)
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
        "multiqc": _StaticRuleRequest((1, 1, 1, 1), (4, 4, 4, 4), (240, 240, 240, 240)),
    }
)


STATIC_SCRATCH_GIB: Mapping[str, tuple[int, int, int, int]] = MappingProxyType(
    {
        "preprocess_alignment": (16, 32, 64, 64),
        "high_cph_filter": (8, 16, 24, 48),
        "bismark_extract": (8, 16, 24, 48),
        "demux_chunk": (24, 24, 24, 24),
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
        "final_sample_manifest": (1, 2048, 30),
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


def demux_resource_request(total_bytes: object, attempt: int = 1, count_only: bool = False) -> ResourceRequest:
    """返回孔板 demux 或 Droplet calling 的调度请求，两次重试仅把内存提高至首次的 1.5 倍、2 倍。

    ``count_only=True``（droplet 细胞调用的全条码计数）按实测斜率外推：峰值内存 ≈
    2 GiB + 0.048×输入GiB（F-real-2609SG，619 GiB R1 实测 29.5 GiB、12.04B reads），
    请求值乘 1.3 余量；runtime 按 1.1 min/GiB。``count_only=False`` 用表内保守值，
    XL 与 L 内存同档。字节数不可得时沿用表内 XL 保底值；Droplet 配对阶段使用
    demux_index/demux_chunk/demux_merge 三条独立固定资源策略。
    """

    attempt_no = _attempt_number(attempt)
    key, row = _rule("demux")
    tier = demux_size_tier(total_bytes)
    tier_index = _tier_index(tier)
    total = _read_count(total_bytes)
    if count_only and tier_index == 3 and total is not None:
        input_gib = total / 1024 ** 3
        base_mem_mb = math.ceil((2 + 0.048 * input_gib) * 1.3 * 1024)
        runtime_min = max(720, math.ceil(1.1 * input_gib))
    else:
        base_mem_mb = row.memory_gib[tier_index] * 1024
        runtime_min = row.runtime_min[tier_index]
    mem_mb = _attempt_memory_mb(base_mem_mb, attempt_no)
    return ResourceRequest(
        policy=STATIC_POLICY + (":count" if count_only else ""),
        rule_key=key,
        tier=tier,
        dna_reads=None,
        threads=row.threads[tier_index],
        mem_mb=mem_mb,
        runtime_min=runtime_min,
    )


def multiqc_resource_request(source_count: object, attempt: int = 1) -> ResourceRequest:
    """按 MultiQC 源文件数缩放的调度请求，两次重试仅把内存分别提高至首次的 1.5 倍、2 倍。

    实测锚点（F-real-2609SG）：52,816 个源、13,203 样本，峰值 RSS 4.35 GiB、
    file-list 模式 14 min / 目录发现模式 17 min。请求内存 = 1.5×(512 MiB + 0.075 MiB/源)；
    runtime = max(60, ceil(源数/600)) min。
    """

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
        ("rastair", "align"): "align_sort_dedup",
        ("bismark", "align"): "bismark_align_dedup",
        ("bismark", "extract"): "bismark_extract",
    }
)

_LOCAL_BACKEND_THREADS: Mapping[str, tuple[int, int, int, int]] = MappingProxyType(
    {
        "align_sort_dedup": (8, 12, 16, 16),
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
    local_job_cap = max(1, min(16, cores // 2 if cores > 1 else 1))
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
    """返回 fused job 固定的 node scratch 容量（MiB）。"""

    key = str(scratch_group).strip()
    try:
        values = STATIC_SCRATCH_GIB[key]
    except KeyError as exc:
        supported = ", ".join(STATIC_SCRATCH_GIB)
        raise ValueError(
            f"unsupported scratch resource group {key!r}; expected one of: {supported}"
        ) from exc
    return int(values[_TIER_ORDER.index(size_tier(dna_reads))] * 1024)


PUBLIC_RESULT_DIR_NAMES: tuple[str, ...] = ("CpG", "QC_Results", "SNP")


@dataclass(frozen=True)
class CellOutputPaths:
    """解析一个 cell 的中间文件路径与所选结果根目录下的公开输出路径。"""

    intermediate_dir: Path
    result_root: Path
    sample_id: str
    backend: str

    @property
    def backend_dir(self) -> Path:
        return self.intermediate_dir / self.backend

    @property
    def cpg(self) -> Path:
        """SnapATAC2 import_values 消费的公开 CpG 输入（pos 为 0-based）。"""
        return self.result_root / "CpG" / f"{self.sample_id}.cpg.tsv.zst"

    @property
    def snp(self) -> Path:
        return self.result_root / "SNP" / f"{self.sample_id}.snps.bed.gz"

    @property
    def snp_index(self) -> Path:
        return Path(f"{self.snp}.tbi")

    @property
    def marked_bam(self) -> Path:
        """TAPS 保留全部记录的坐标排序 BAM；重复仅标记，不删除。"""
        return self.backend_dir / f"{self.sample_id}.marked.bam"

    @property
    def marked_bam_index(self) -> Path:
        return Path(f"{self.marked_bam}.bai")

    @property
    def rastair_qc(self) -> Path:
        return self.backend_dir / f"{self.sample_id}.rastair_qc.json"

    @property
    def dna_bam(self) -> Path:
        """Biscuit high-CpH 过滤后的 DNA BAM（temp 中间件，pileup 消费）。"""
        return self.backend_dir / f"{self.sample_id}.dna.bam"

    @property
    def dna_bam_index(self) -> Path:
        return Path(f"{self.dna_bam}.bai")

    @property
    def retained_bam_files(self) -> tuple[Path, ...]:
        """返回 keep_final_bam 要求保留且由当前 backend 实际生成的 BAM/索引。"""
        if self.backend == "rastair":
            return self.marked_bam, self.marked_bam_index
        if self.backend == "biscuit":
            return self.dna_bam, self.dna_bam_index
        return (self.processed_bam,)

    @property
    def rna_bam(self) -> Path:
        """SRD 协议恒定保留的被筛除 BAM（RNA 部分原始数据，位于 02_work/RNA）。"""
        return self.intermediate_dir / "RNA" / "BAM" / f"{self.sample_id}.rna.bam"

    @property
    def rna_bam_index(self) -> Path:
        return Path(f"{self.rna_bam}.bai")

    @property
    def post_alignment_bam(self) -> Path:
        return self.backend_dir / f"{self.sample_id}.post_alignment.bam"

    @property
    def post_alignment_bam_index(self) -> Path:
        return Path(f"{self.post_alignment_bam}.bai")


    @property
    def high_cph_summary(self) -> Path:
        suffix = (
            "bsconv_filter_summary.json"
            if self.backend == "biscuit"
            else "nonconversion_filter_summary.json"
        )
        return self.backend_dir / f"{self.sample_id}.{suffix}"

    @property
    def high_cph_multiqc(self) -> Path:
        suffix = (
            "bsconv_filter_mqc.json"
            if self.backend == "biscuit"
            else "nonconversion_filter_mqc.json"
        )
        return self.intermediate_dir / "multiqc" / self.backend / f"{self.sample_id}.{suffix}"

    @property
    def cutadapt_json(self) -> Path:
        return self.intermediate_dir / "cutadapt" / f"{self.sample_id}.cutadapt.json"

    @property
    def dupsifter_multiqc(self) -> Path:
        return self.intermediate_dir / "biscuit" / f"{self.sample_id}.dupsifter_mqc.json"

    @property
    def high_cph_read_names(self) -> Path:
        return self.backend_dir / f"{self.sample_id}.high_cph_reads.txt"

    @property
    def processed_bam(self) -> Path:
        return self.backend_dir / f"{self.sample_id}.processed.bam"

    @property
    def removed_bam(self) -> Path:
        return self.backend_dir / f"{self.sample_id}.nonconversion_removed.bam"


    @property
    def biscuit_qc_dir(self) -> Path:
        return self.intermediate_dir / "biscuit" / "qc" / self.sample_id

    def biscuit_qc_source(self, suffix: str) -> Path:
        return self.biscuit_qc_dir / f"{self.sample_id}_{suffix}"

    def biscuit_qc_sources(self, suffixes: Iterable[str]) -> tuple[Path, ...]:
        return tuple(self.biscuit_qc_source(suffix) for suffix in suffixes)

    @property
    def bismark_dir(self) -> Path:
        return self.intermediate_dir / "bismark" / self.sample_id

    @property
    def bismark_align_report(self) -> Path:
        return self.bismark_dir / f"{self.sample_id}_PE_report.txt"

    @property
    def bismark_dedup_report(self) -> Path:
        return self.bismark_dir / f"{self.sample_id}.deduplication_report.txt"

    @property
    def droplet_umi_qc(self) -> Path:
        return self.bismark_dir / f"{self.sample_id}.umi_dedup.json"

    @property
    def bismark_extract_report(self) -> Path:
        return self.bismark_dir / f"{self.sample_id}_splitting_report.txt"

    @property
    def bismark_mbias(self) -> Path:
        return self.bismark_dir / f"{self.sample_id}.M-bias.txt"

    @property
    def bismark_cov(self) -> Path:
        return self.bismark_dir / f"{self.sample_id}.bismark.cov.gz"


@dataclass(frozen=True)
class ProjectLayout:
    """在不查询实时 workflow 状态的前提下解析输出路径。"""

    intermediate_dir: Path
    result_root: Path
    backend: str

    def __init__(
        self,
        intermediate_dir: str | Path,
        result_root: str | Path,
        backend: str,
    ) -> None:
        if backend not in {"biscuit", "bismark", "rastair"}:
            raise ValueError(f"Unsupported methylation backend: {backend!r}")
        object.__setattr__(self, "intermediate_dir", Path(intermediate_dir))
        object.__setattr__(self, "result_root", Path(result_root))
        object.__setattr__(self, "backend", backend)

    def cell(self, sample_id: str) -> CellOutputPaths:
        return CellOutputPaths(
            self.intermediate_dir,
            self.result_root,
            str(sample_id),
            self.backend,
        )


_COUNT = r"([0-9][0-9,]*)"


FILTER_SUMMARY_SCHEMA_VERSION = 3
FILTER_SUMMARY_FIELDS = frozenset(
    {
        "schema_version",
        "backend",
        "sample",
        "metric_unit",
        "candidate_read_pairs",
        "flagged_read_pairs",
        "high_cph_fraction",
        "excluded_read_pairs",
        "final_retained_pairs",
        "threshold",
        "minimum_count",
        "filter_applied",
        "protocol",
        "high_cph_role",
    }
)


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


DROPLET_DEDUP_POLICY = "umi_tools_directional_physical_r1_ignore_tlen"


def deduplicate_droplet_bam(bam: str | Path, output: str | Path, qc: str | Path, *, sample: str, scratch: str | Path, umi_tools: str, threads: int = 1) -> None:
    """按物理 R1 与 UMI 去重，并恢复 Bismark flags、配对顺序及甲基化标签。"""
    import pysam
    from importlib.metadata import version

    if SAMPLE_ID_RE.fullmatch(sample) is None:
        raise ValueError("Invalid Droplet UMI sample identity")
    if version("umi_tools") != "1.1.6":
        raise ValueError("Droplet dedup requires umi_tools 1.1.6")
    root = Path(tempfile.mkdtemp(prefix="droplet_umi_", dir=scratch))
    total = kept = umi_n = 0
    seen = set()
    try:
        normalized = root / "normalized.bam"
        with pysam.AlignmentFile(str(bam), "rb") as source, pysam.AlignmentFile(str(normalized), "wb", template=source) as target:
            iterator = iter(source)
            for r1 in iterator:
                r2 = next(iterator, None)
                if r2 is None or r1.query_name != r2.query_name or r1.query_name in seen:
                    raise ValueError("Droplet BAM requires unique, adjacent complete physical read pairs")
                seen.add(r1.query_name)
                if (r1.flag, r2.flag) not in {(99, 147), (163, 83), (147, 99), (83, 163)}:
                    raise ValueError("Unexpected Bismark paired flags")
                umi = r1.query_name.rsplit(":", 1)[-1]
                if not re.fullmatch(r"[ACGTN]{12}", umi):
                    raise ValueError(f"Droplet QNAME must end with a 12 bp UMI: {r1.query_name!r}")
                for index, record in enumerate((r1, r2)):
                    if any(not record.has_tag(tag) for tag in ("XM", "XR", "XG")) or record.has_tag("ZF"):
                        raise ValueError("Bismark methylation tags missing or reserved ZF tag present")
                    record.set_tag("ZF", record.flag, value_type="i")
                    record.flag = (record.flag & ~(64 | 128)) | (64 if index == 0 else 128)
                    record.set_tag("UR", umi, value_type="Z")
                    target.write(record)
                total += 1
                umi_n += int("N" in umi)
        del seen
        coordinate = root / "coordinate.bam"
        selected = root / "selected.bam"
        ordered = root / "ordered.bam"
        pysam.sort("-@", str(max(1, threads)), "-o", str(coordinate), str(normalized))
        pysam.index(str(coordinate))
        if total:
            subprocess.run([umi_tools, "dedup", "--stdin", str(coordinate), "--stdout", str(selected),
                            "--paired", "--ignore-tlen", "--extract-umi-method", "tag", "--umi-tag", "UR",
                            "--method", "directional", "--edit-distance-threshold", "1", "--random-seed", "1"], check=True)
        else:
            shutil.copyfile(coordinate, selected)
        pysam.sort("-n", "-@", str(max(1, threads)), "-o", str(ordered), str(selected))
        with pysam.AlignmentFile(str(ordered), "rb") as source, pysam.AlignmentFile(str(output), "wb", template=source) as target:
            iterator = iter(source)
            for r1 in iterator:
                r2 = next(iterator, None)
                if r2 is None or r1.query_name != r2.query_name or not r1.is_read1 or not r2.is_read2:
                    raise ValueError("UMI dedup did not preserve complete physical pairs")
                for record in (r1, r2):
                    record.flag = record.get_tag("ZF")
                    record.set_tag("ZF", None)
                    target.write(record)
                kept += 1
        if kept > total:
            raise ValueError("UMI pair accounting is inconsistent")
        Path(qc).write_text(json.dumps({"schema_version": 1, "sample": sample, "dedup_policy": DROPLET_DEDUP_POLICY,
            "umi_tools_version": "1.1.6", "method": "directional", "edit_distance": 1, "random_seed": 1,
            "input_pairs": total, "retained_pairs": kept, "removed_pairs": total - kept,
            "umi_with_n_pairs": umi_n}, indent=2) + "\n")
    finally:
        shutil.rmtree(root)


def read_droplet_umi_qc(path: str | Path, sample: str) -> dict[str, int]:
    """校验 UMI 去重参数与 pair 守恒，返回统一的去重计数。"""
    data = _read_json(Path(path), "Droplet UMI QC")
    for key, value in {"schema_version": 1, "sample": sample, "dedup_policy": DROPLET_DEDUP_POLICY, "umi_tools_version": "1.1.6", "method": "directional", "edit_distance": 1, "random_seed": 1}.items():
        if data.get(key) != value:
            raise ValueError(f"Invalid Droplet UMI QC {key}: {path}")
    numbers = {key: _integer(data.get(key), key, Path(path)) for key in ("input_pairs", "retained_pairs", "removed_pairs", "umi_with_n_pairs")}
    if numbers["retained_pairs"] + numbers["removed_pairs"] != numbers["input_pairs"] or numbers["umi_with_n_pairs"] > numbers["input_pairs"]:
        raise ValueError(f"Droplet UMI QC pair accounting mismatch: {path}")
    return {"total": numbers["input_pairs"], "removed": numbers["removed_pairs"], "leftover": numbers["retained_pairs"]}


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

    dedup = read_droplet_umi_qc(dedup_report, sample) if protocol == "droplet" else parse_dedup_report(dedup_report)
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

    return {
        "total_pairs": int(alignment["total"]),
        "accepted_pairs": accepted_pairs,
        "unmapped_pairs": int(alignment["no_alignment"]) + discarded_pairs,
        "ambiguous_pairs": int(alignment["non_unique"]),
        "duplicate_pairs": duplicate_pairs,
        "postdedup_pairs": postdedup_pairs,
        "final_retained_pairs": final_retained_pairs,
    }


def _count_bam_pairs(bam: Path, samtools: str) -> int:
    result = subprocess.run(
        [samtools, "view", "-c", str(bam)],
        check=True, text=True, stdout=subprocess.PIPE,
    )
    records = int(result.stdout.strip())
    if records < 0 or records % 2:
        raise ValueError("Bismark processed BAM contains an odd or negative alignment count")
    return records // 2


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


logger = logging.getLogger("dna_pipeline")
SAMPLE_ID_RE = re.compile(r"^(?!\.{1,2}$)[A-Za-z0-9_.-]+$")
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


def resolve_fai_contigs(
    requested: Sequence[object],
    fai_path: str | Path,
    *,
    field: str = "high_cph.excluded_contigs",
) -> tuple[str, ...]:
    """把配置的 contig 解析为 FAI 中的权威精确名称。

    精确拼写优先；唯一的大小写不敏感匹配可被项目接受，
    歧义或缺失的值在作业启动前即失败。
    """
    path = Path(fai_path)
    if not path.is_file():
        raise FileNotFoundError(f"Reference FASTA index is missing: {path}")
    contigs: list[str] = []
    for line_no, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        name = line.split("\t", 1)[0].strip()
        if not name:
            raise ValueError(f"Reference FAI line {line_no} has an empty contig name: {path}")
        contigs.append(name)
    if not contigs:
        raise ValueError(f"Reference FAI contains no contigs: {path}")
    if len(set(contigs)) != len(contigs):
        raise ValueError(f"Reference FAI repeats a contig name: {path}")

    exact = set(contigs)
    folded: dict[str, list[str]] = {}
    for contig in contigs:
        folded.setdefault(contig.casefold(), []).append(contig)
    resolved: list[str] = []
    for raw in requested:
        value = str(raw or "").strip()
        if not value:
            raise ValueError(f"{field} must not contain an empty contig name")
        if value in exact:
            canonical = value
        else:
            matches = folded.get(value.casefold(), [])
            if not matches:
                raise ValueError(
                    f"{field} contains {value!r}, which is absent from reference FAI {path}"
                )
            if len(matches) != 1:
                raise ValueError(
                    f"{field} value {value!r} is ambiguous in reference FAI {path}: {matches}"
                )
            canonical = matches[0]
        if canonical in resolved:
            raise ValueError(
                f"{field} resolves more than once to reference contig {canonical!r}"
            )
        resolved.append(canonical)
    return tuple(resolved)


def resolve_reference(
    config: Mapping[str, object],
    species: str,
    project_dir: Path,
    *,
    pipeline_root: Path,
) -> str:
    """解析某物种配置的 reference FASTA 路径。"""
    refs = config.get("references", {})
    if not isinstance(refs, Mapping) or species not in refs:
        raise ValueError(
            f"Species {species!r} is not configured under references. "
            f"Add references.{species}: /path/to/reference.fa to 00_config/config.yaml."
        )

    ref_path = str(refs[species] or "").strip()
    if not ref_path:
        raise ValueError(f"Reference path for species {species!r} is empty")

    configured = Path(ref_path).expanduser()
    pipeline_root = Path(pipeline_root).expanduser().resolve()
    candidates = (
        [configured.resolve()]
        if configured.is_absolute()
        else [(pipeline_root / configured).resolve(), (project_dir / configured).resolve()]
    )
    seen = set()
    unique = []
    for candidate in candidates:
        key = str(candidate)
        if key not in seen:
            unique.append(candidate)
            seen.add(key)
    for candidate in unique:
        if candidate.is_file():
            return str(candidate)

    checked = "\n".join(f"  - {candidate}" for candidate in unique)
    setup_script = pipeline_root / "core" / "doctor.sh"
    if species in {"hg38", "mm10"}:
        hint = (
            "For standard hg38/mm10 references, place hg38.fa, mm10.fa, "
            "lambda.fa and pUC19.fa under resources/reference_source/, then run:\n"
            f"  bash \"{setup_script}\" reference"
        )
    else:
        hint = (
            "Custom species are user-managed. Assemble the final FASTA and build "
            "samtools FAI and the selected backend indexes manually, "
            "then point references.<species> to that FASTA."
        )
    raise FileNotFoundError(
        f"Reference genome for species {species!r} was not found.\n"
        f"Configured path: {ref_path}\nChecked:\n{checked}\n\n"
        f"{hint}"
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

    def _sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()

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
                digests[mate][str(path)] = _sha256(path)
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

    _apply_fastq_renames(dir_raw, plans)
    return plans


def _apply_fastq_renames(dir_raw: Path, plans: Sequence[tuple[Path, Path]]) -> None:
    applied: list[tuple[Path, Path]] = []

    def _rollback() -> None:
        for source, target in reversed(applied):
            target.rename(source)

    try:
        for source, target in plans:
            if target.exists() or target.is_symlink():
                raise FileExistsError(f"FASTQ rename target already exists: {target}")
            source.rename(target)
            applied.append((source, target))
    except OSError:
        _rollback()
        raise
    try:
        discover_samples(dir_raw)
    except Exception as exc:
        _rollback()
        raise RuntimeError(
            f"Collision renames were rolled back because the merged inventory under "
            f"{dir_raw} still does not validate: {exc}"
        ) from exc
    logger.info(
        "Normalized %d colliding FASTQ file(s) into chunk-style segments", len(applied)
    )


MANAGED_SPECIES = frozenset({"hg38", "mm10"})


_PRIMARY_CONTIGS: Mapping[str, tuple[tuple[str, int], ...]] = {
    "hg38": (
        ("chr1", 248956422), ("chr2", 242193529), ("chr3", 198295559),
        ("chr4", 190214555), ("chr5", 181538259), ("chr6", 170805979),
        ("chr7", 159345973), ("chr8", 145138636), ("chr9", 138394717),
        ("chr10", 133797422), ("chr11", 135086622), ("chr12", 133275309),
        ("chr13", 114364328), ("chr14", 107043718), ("chr15", 101991189),
        ("chr16", 90338345), ("chr17", 83257441), ("chr18", 80373285),
        ("chr19", 58617616), ("chr20", 64444167), ("chr21", 46709983),
        ("chr22", 50818468), ("chrX", 156040895), ("chrY", 57227415),
        ("chrM", 16569),
    ),
    "mm10": (
        ("chr1", 195471971), ("chr2", 182113224), ("chr3", 160039680),
        ("chr4", 156508116), ("chr5", 151834684), ("chr6", 149736546),
        ("chr7", 145441459), ("chr8", 129401213), ("chr9", 124595110),
        ("chr10", 130694993), ("chr11", 122082543), ("chr12", 120129022),
        ("chr13", 120421639), ("chr14", 124902244), ("chr15", 104043685),
        ("chr16", 98207768), ("chr17", 94987271), ("chr18", 90702639),
        ("chr19", 61431566), ("chrX", 171031299), ("chrY", 91744698),
        ("chrM", 16299),
    ),
}
SPIKE_CONTIGS: Mapping[str, int] = {"lambda": 48502, "pUC19": 2686}
MANAGED_CONTIGS = {
    species: [*(primaries), *SPIKE_CONTIGS.items()]
    for species, primaries in _PRIMARY_CONTIGS.items()
}


def _fai_contigs(path: Path) -> list[tuple[str, int]]:

    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except UnicodeError as exc:
        raise ValueError(f"reference FAI is not valid text: {path}") from exc

    records: list[tuple[str, int]] = []
    seen: set[str] = set()
    for line_number, line in enumerate(lines, start=1):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 5 or not fields[0]:
            raise ValueError(f"invalid reference FAI line {line_number}: {path}")
        name = fields[0]
        if name in seen:
            raise ValueError(f"duplicate reference FAI contig {name!r}: {path}")
        try:
            length, offset, line_bases, line_width = map(int, fields[1:5])
        except ValueError as exc:
            raise ValueError(
                f"invalid numeric reference FAI line {line_number}: {path}"
            ) from exc
        if length <= 0 or offset < 0 or line_bases <= 0 or line_width < line_bases:
            raise ValueError(f"invalid reference FAI geometry on line {line_number}: {path}")
        seen.add(name)
        records.append((name, length))
    if not records:
        raise ValueError(f"reference FAI contains no contigs: {path}")
    return records


def _atomic_text_path(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    return path.with_name(f".{path.name}.tmp.{os.getpid()}")


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
        if name not in {"lambda", "pUC19"}
    ]
    contig_lengths = dict(host_contigs)
    expected_build, expected_release = _annotation_release_signature(species)
    header_lines: list[str] = []
    tss_tmp = _atomic_text_path(tss)
    cpg_tmp = _atomic_text_path(cpg)
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

        os.replace(tss_tmp, tss)
        os.replace(cpg_tmp, cpg)
    finally:
        for temporary in (tss_tmp, cpg_tmp):
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


def reference_identity(
    reference: str | Path,
    *,
    expected_species: str | None = None,
) -> Mapping[str, object]:
    """返回所有 backend 共用的 FASTA + FAI 身份。

    backend 专属 index 单独另行校验记录；未选中的 backend 绝不会成为重算或 run-start gate。
    """
    public_ref = Path(reference).expanduser().absolute()
    if not public_ref.is_file() or public_ref.stat().st_size == 0:
        raise FileNotFoundError(f"reference FASTA is missing or empty: {public_ref}")
    resolved_ref = public_ref.resolve()
    fai = Path(f"{resolved_ref}.fai")
    if not fai.is_file() or fai.stat().st_size == 0:
        raise FileNotFoundError(f"reference FAI is missing or empty: {fai}")
    observed = _fai_contigs(fai)
    if expected_species in MANAGED_SPECIES and observed != MANAGED_CONTIGS[expected_species]:
        raise ValueError(
            f"{expected_species} reference contig/order/length contract mismatch: {resolved_ref}"
        )
    ref_stat = resolved_ref.stat()
    return {
        "public_path": str(public_ref),
        "resolved_path": str(resolved_ref),
        "size_bytes": ref_stat.st_size,
        "mtime_ns": ref_stat.st_mtime_ns,
        "fai": stable_file_identity(fai, hash_content=False),
    }


BWA_INDEX_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")

BISCUIT_INDEX_SUFFIXES = (
    ".bis.pac", ".bis.amb", ".bis.ann", ".par.bwt", ".par.sa", ".dau.bwt", ".dau.sa",
)

BISCUIT_TRANSIENT_SUFFIXES = (".par.pac", ".dau.pac")


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


def fai_contract_check(mode: str, label: str, fai: Path, fasta: Path) -> bool:
    """校验 source/assembled FASTA 是否满足标准 contig 长度契约。

    只读取很小的 .fai 并通过 FAI 几何 + FASTA 大小的算术检查排除截断文件，
    不扫描序列本身。mode=source 允许 UCSC FASTA 中的额外 contig（primary 与
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
        return label in MANAGED_CONTIGS and length_records == MANAGED_CONTIGS[label]
    return False


def backend_active_index_paths(
    reference: str | Path,
    backend: str,
    *,
    bismark_local_alignment: bool = False,
    bismark_root: str | Path | None = None,
) -> list[Path]:
    """枚举当前 backend/模式必需的完整 index 家族；构建期可指定候选 Bismark 根目录。

    Bowtie2 以 .1.bt2 优先、其次 .1.bt2l 选择六文件家族；缺失成员也必须进入清单，不能被扫描漏掉。
    """
    reference = Path(reference)
    if backend == "rastair":
        return [Path(f"{reference}.bwa/genome{suffix}") for suffix in BWA_INDEX_SUFFIXES]
    if backend == "biscuit":
        return [Path(f"{reference}{suffix}") for suffix in BISCUIT_INDEX_SUFFIXES]
    if backend == "bismark":
        index_dir = Path(bismark_root or f"{reference}.bismark") / "Bisulfite_Genome"
        prefixes = (
            [index_dir / "CT_conversion" / "BS_CT", index_dir / "GA_conversion" / "BS_GA"]
            if bismark_local_alignment
            else [index_dir / "Combined" / "BS_combined"]
        )
        paths = []
        for prefix in prefixes:
            extension = "bt2l" if not Path(f"{prefix}.1.bt2").is_file() and Path(f"{prefix}.1.bt2l").is_file() else "bt2"
            paths.extend(Path(f"{prefix}.{part}.{extension}") for part in ("1", "2", "3", "4", "rev.1", "rev.2"))
        return paths
    raise ValueError(f"Unsupported methylation backend: {backend!r}")


SUPPORTED_BISMARK_LIBRARY_TYPES = ("directional", "non_directional", "pbat")
SUPPORTED_BISCUIT_LIBRARY_TYPES = ("directional", "non_directional")


CONFIG_SECTION_DEFAULTS: dict[str, dict[str, object]] = {
    "demux": {
        "min_matched_read_pairs": 10,
        "dna_w_spacer_len": 0,
    },
    "high_cph": {
        "excluded_contigs": ["pUC19", "lambda", "chrM"],
    },
    "biscuit": {
        "library_type": "non_directional",
        "high_cph_retention_threshold": 0.7,
        "generate_snp": False,
    },
    "bismark": {
        "library_type": "non_directional",
        "local_alignment": False,
        "non_conversion_percentage_cutoff": 70,
        "non_conversion_minimum_count": 5,
    },
    "runtime": {
        "local_executor_cores": 32,
        "latency_wait_seconds": 120,
    },
    "slurm": {
        "jobs": 30,
        "controller_cores": 2,
        "account": None,
        "partition": None,
        "qos": "huge",
        "node_tmpdir": "/tmp",
    },
    "retention": {
        "keep_final_bam": False,
    },
}


def apply_config_defaults(config: dict[str, object]) -> dict[str, object]:
    """把内置默认合并进项目 config；缺 section 补 section，缺键补键，不改显式值。"""
    for section, defaults in CONFIG_SECTION_DEFAULTS.items():
        if section not in config:
            config[section] = copy.deepcopy(defaults)
            continue
        current = config[section]
        if not isinstance(current, dict):
            raise ValueError(f"config section {section!r} must be a mapping")
        for key, value in defaults.items():
            current.setdefault(key, copy.deepcopy(value))
    return config


def finalize_project_config(config: Mapping[str, object]) -> dict[str, object]:
    """合并内置默认并完成 schema 与跨字段语义校验（项目 config 的唯一加载链）。"""
    import yaml
    from jsonschema import validate as validate_jsonschema

    merged = apply_config_defaults(copy.deepcopy(dict(config)))
    schema = yaml.safe_load(CONFIG_SCHEMA_YAML)
    validate_jsonschema(instance=merged, schema=schema)
    validate_config_semantics(merged)
    return merged


def load_project_config(path: Path | str) -> dict[str, object]:
    """读取项目 config YAML 并走唯一 finalize 加载链。"""
    import yaml

    config = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    if not isinstance(config, dict):
        raise ValueError(f"project config root must be a mapping: {path}")
    return finalize_project_config(config)


@dataclass(frozen=True)
class BackendPolicy:
    """单一 methylation backend 的静态能力契约。"""

    name: str
    supports_snp: bool
    uses_dupsifter: bool


_BACKEND_POLICIES = {
    "biscuit": BackendPolicy("biscuit", supports_snp=True, uses_dupsifter=True),
    "bismark": BackendPolicy("bismark", supports_snp=False, uses_dupsifter=False),
    "rastair": BackendPolicy("rastair", supports_snp=False, uses_dupsifter=False),
}


def get_backend(config: dict[str, Any]) -> BackendPolicy:
    """返回选定的 backend，并拒绝该 backend 不支持的功能组合。"""
    analysis = config.get("analysis")
    if not isinstance(analysis, dict):
        raise ValueError("analysis must be a mapping")
    name = str(analysis.get("methylation_backend", "")).strip().lower()
    if name not in _BACKEND_POLICIES:
        raise ValueError(
            f"Unsupported methylation backend {name!r}; choose biscuit, bismark or rastair"
        )
    if (analysis.get("protocol") == "taps") != (name == "rastair"):
        raise ValueError("protocol taps requires rastair; Cabernet/SRD require biscuit or bismark")
    policy = _BACKEND_POLICIES[name]
    biscuit = config.get("biscuit") or {}
    if bool(biscuit.get("generate_snp")) and not policy.supports_snp:
        raise ValueError(
            "biscuit.generate_snp=true is currently supported only by the Biscuit backend"
        )
    return policy


def bismark_library_flag(config: dict[str, Any]) -> str:
    """把配置的 library type 翻译为 Bismark CLI flag。"""
    bismark = config.get("bismark")
    if not isinstance(bismark, dict):
        raise ValueError("bismark must be a mapping")
    value = str(bismark.get("library_type", "")).strip().lower()
    if value not in SUPPORTED_BISMARK_LIBRARY_TYPES:
        raise ValueError(
            "bismark.library_type must be directional, non_directional, or pbat"
        )
    return {
        "directional": "",
        "non_directional": "--non_directional",
        "pbat": "--pbat",
    }[value]


def bismark_alignment_instance_count(config: Mapping[str, object]) -> int:
    """返回当前 Bismark 模式使用的 Bowtie2 instance 数，用于划分 rule 总 CPU。"""
    bismark = config.get("bismark")
    if not isinstance(bismark, Mapping):
        raise ValueError("bismark must be a mapping")
    value = str(bismark.get("library_type", "")).strip().lower()
    if value not in SUPPORTED_BISMARK_LIBRARY_TYPES:
        raise ValueError(
            "bismark.library_type must be directional, non_directional, or pbat"
        )
    if not bool(bismark.get("local_alignment", False)):
        return 1
    return 4 if value == "non_directional" else 2


def biscuit_library_mode(config: Mapping[str, object]) -> int:
    """把命名的 paired-end BISCUIT library type 翻译为 ``align -b`` 的取值。"""
    biscuit = config.get("biscuit")
    if not isinstance(biscuit, Mapping):
        raise ValueError("biscuit must be a mapping")
    value = str(biscuit.get("library_type", "")).strip().lower()
    if value not in SUPPORTED_BISCUIT_LIBRARY_TYPES:
        raise ValueError("biscuit.library_type must be directional or non_directional")
    return {"non_directional": 0, "directional": 1}[value]


def validate_config_semantics(config: Mapping[str, object]) -> None:
    """校验 JSON Schema 无法表达的跨字段约束。"""
    high_cph = config.get("high_cph")
    if not isinstance(high_cph, Mapping):
        raise ValueError("high_cph must be a mapping")

    biscuit = config.get("biscuit")
    runtime = config.get("runtime")
    slurm = config.get("slurm")
    if not isinstance(biscuit, Mapping):
        raise ValueError("biscuit must be a mapping")
    if not isinstance(runtime, Mapping):
        raise ValueError("runtime must be a mapping")
    if not isinstance(slurm, Mapping):
        raise ValueError("slurm must be a mapping")
    node_tmpdir = slurm.get("node_tmpdir")
    if node_tmpdir is not None:
        if not isinstance(node_tmpdir, str) or not str(node_tmpdir).strip():
            raise ValueError("slurm.node_tmpdir must be a non-empty path or null")
        if Path(str(node_tmpdir)).expanduser().resolve(strict=False) == Path(
            Path(str(node_tmpdir)).expanduser().resolve(strict=False).anchor
        ):
            raise ValueError("slurm.node_tmpdir must not resolve to the filesystem root")
    for section_name, section, field, minimum, maximum in (
        ("biscuit", biscuit, "high_cph_retention_threshold", 0.0, 1.0),
    ):
        value = section.get(field)
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(float(value))
            or float(value) < minimum
            or (maximum is not None and float(value) > maximum)
        ):
            bounds = f"[{minimum:g}, {maximum:g}]" if maximum is not None else f">= {minimum:g}"
            raise ValueError(f"{section_name}.{field} must be a finite number within {bounds}")
    get_backend(dict(config))
    biscuit_library_mode(config)
    bismark_library_flag(dict(config))
    if config["analysis"]["protocol"] == "droplet":
        if config["analysis"]["methylation_backend"] != "bismark" or config["bismark"]["library_type"] != "non_directional":
            raise ValueError("Droplet requires bismark/non_directional")
        if config["demux"]["dna_w_spacer_len"] != 0:
            raise ValueError("Droplet has a fixed structure; dna_w_spacer_len must be 0")


@dataclass(frozen=True)
class ProtocolPolicy:
    name: str
    demux_mode: str
    requires_rna_barcode: bool
    rna_output_role: str
    high_cph_role: str

    publishes_rna_bam: bool

    @property
    def raw_input_columns(self) -> tuple[str, ...]:
        """返回本项目需要 FASTQ 的 manifest 列；外部 RNA 关联不参与拆分。"""
        return ("dna_raw_sample", "rna_sample") if self.name == "srd" else ("dna_raw_sample",)


@dataclass(frozen=True)
class RawRoutePolicy:
    route_id: str
    protocol: str
    demux_mode: str
    downstream_dna: bool
    rna_output_role: str
    high_cph_role: str


_PROTOCOLS = {
    "droplet": ProtocolPolicy(
        name="droplet", demux_mode="dna-only-droplet", requires_rna_barcode=False,
        rna_output_role="not_used", high_cph_role="cdna_contamination", publishes_rna_bam=False,
    ),
    "taps": ProtocolPolicy(
        name="taps", demux_mode="dna-only-taps", requires_rna_barcode=False,
        rna_output_role="not_used", high_cph_role="not_applicable",
        publishes_rna_bam=False,
    ),
    "cabernet": ProtocolPolicy(
        name="cabernet", demux_mode="dna-only", requires_rna_barcode=False,
        rna_output_role="not_used", high_cph_role="cdna_contamination",
        publishes_rna_bam=False,
    ),
    "srd": ProtocolPolicy(
        name="srd", demux_mode="dna-rna", requires_rna_barcode=True,
        rna_output_role="audit_only", high_cph_role="residual_high_cph",
        publishes_rna_bam=True,
    ),
}

def _raw_route(
    route_id: str,
    protocol: str,
    downstream_dna: bool,
    *,
    rna_output_role: str | None = None,
    high_cph_role: str | None = None,
) -> RawRoutePolicy:
    policy = _PROTOCOLS[protocol]
    return RawRoutePolicy(
        route_id=route_id,
        protocol=protocol,
        downstream_dna=downstream_dna,
        demux_mode=policy.demux_mode,
        rna_output_role=rna_output_role or policy.rna_output_role,
        high_cph_role=high_cph_role or policy.high_cph_role,
    )


_RAW_ROUTES = {
    "droplet_dna": _raw_route("droplet_dna", "droplet", True),
    "taps_dna_nucleus": _raw_route("taps_dna_nucleus", "taps", True),
    "cabernet_dna_nucleus": _raw_route("cabernet_dna_nucleus", "cabernet", True),
    "srd_dna_tube": _raw_route("srd_dna_tube", "srd", True),
    "srd_rna_enrichment_tube": _raw_route(
        "srd_rna_enrichment_tube",
        "srd",
        False,
        rna_output_role="audit_only",
        high_cph_role="not_applicable",
    ),
}

_BARCODE_RE = re.compile(r"^[ACGT]+$")
_SUPPORTED_BARCODE_LENGTHS = {8, 10}


def _canonical_protocol_name(name: str) -> str:
    protocol = str(name or "").strip().lower()
    if protocol not in _PROTOCOLS:
        raise ValueError(
            f"Unsupported pipeline protocol={protocol!r}; expected one of: "
            + ", ".join(sorted(_PROTOCOLS))
        )
    return protocol


def get_protocol_policy(config: Mapping[str, object]) -> ProtocolPolicy:
    protocol = _canonical_protocol_name(str(config["analysis"]["protocol"]))
    return _PROTOCOLS[protocol]


def get_protocol(name: str) -> ProtocolPolicy:
    return _PROTOCOLS[_canonical_protocol_name(name)]


def default_dna_route(protocol: str) -> str:
    """返回协议唯一的下游 DNA 原始样本路线。"""
    return next(r.route_id for r in _RAW_ROUTES.values() if r.protocol == _canonical_protocol_name(protocol) and r.downstream_dna)


def canonical_raw_route_id(route_name: str) -> str:
    return str(route_name or "").strip().lower()


def get_raw_route(route_name: str, protocol_name: str) -> RawRoutePolicy:
    protocol = _canonical_protocol_name(protocol_name)
    route_id = canonical_raw_route_id(route_name)
    route = _RAW_ROUTES.get(route_id)
    if route is None or route.protocol != protocol:
        valid = sorted(r.route_id for r in _RAW_ROUTES.values() if r.protocol == protocol)
        raise ValueError(
            f"Unsupported raw-sample route={route_name!r} for protocol={protocol!r}; "
            f"expected one of: {', '.join(valid)}"
        )
    return route

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


_EXPORTED_SCALARS = {
    "CFG_LOCAL_CORES": "runtime.local_executor_cores",
    "CFG_LATENCY_WAIT": "runtime.latency_wait_seconds",
    "CFG_PROTOCOL": "analysis.protocol",
    "CFG_SPECIES": "species",
    "CFG_METHYLATION_BACKEND": "analysis.methylation_backend",
    "CFG_SLURM_JOBS": "slurm.jobs",
    "CFG_SLURM_CONTROLLER_CORES": "slurm.controller_cores",
    "CFG_SLURM_ACCOUNT": "slurm.account",
    "CFG_SLURM_PARTITION": "slurm.partition",
    "CFG_SLURM_QOS": "slurm.qos",
    "CFG_SLURM_NODE_TMPDIR": "slurm.node_tmpdir",
}


def export_config_scalars(config: Mapping[str, object]) -> dict[str, str]:
    """提取平台路由所需的 config 叶子标量（None 渲染为空串，bool 渲染为 true/false）。"""
    rendered: dict[str, str] = {}
    for var, dotted in _EXPORTED_SCALARS.items():
        value: object = config
        for key in dotted.split("."):
            if not isinstance(value, Mapping) or key not in value:
                raise ValueError(f"required config field is missing: {dotted}")
            value = value[key]
        if value is None:
            rendered[var] = ""
        elif isinstance(value, bool):
            rendered[var] = "true" if value else "false"
        elif isinstance(value, (str, int, float)):
            rendered[var] = str(value)
        else:
            raise ValueError(f"config field must be a scalar: {dotted}")
    return rendered


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
    analysis["methylation_backend"] = meta["methylation_backend"]
    analysis["protocol"] = meta["protocol"]
    cfg.setdefault("biscuit", {})["generate_snp"] = False
    if meta["methylation_backend"] == "bismark":
        cfg.setdefault("bismark", {})["library_type"] = "non_directional"
    apply_config_defaults(cfg)

def synthesize_route_config(source: Path, destination: Path, route: str) -> None:
    """把项目 config 与 route 覆盖合成为 test 模式 effective config。"""
    import yaml

    cfg = yaml.safe_load(source.read_text(encoding="utf-8")) or {}
    apply_common_route_config(cfg, resolve_route(route))


    runtime = cfg.setdefault("runtime", {})
    runtime["local_executor_cores"] = max(8, int(runtime.get("local_executor_cores", 8)))
    destination.write_text(
        yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True), encoding="utf-8"
    )


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
    samtools = Path(__file__).resolve().parent.parent / "conda/current" / backend / "bin/samtools"
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
    samtools = Path(__file__).resolve().parent.parent / "conda/current/rastair/bin/samtools"
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
    if row["high_cph_role"] != "not_applicable" or not math.isnan(float(row["non_cpg_methylation_pct"])):
        raise ValueError("TAPS CpH must be unmeasured")
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
    fixture_root = Path(__file__).resolve().parent.parent / "resources/doctor_benchmark"
    fixture = json.loads((fixture_root / "benchmark.json").read_text())
    cfg["demux"]["min_matched_read_pairs"] = 10
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


SNAPSHOT_SCHEMA_VERSION = 10

_SOURCE_ROOTS = ("core",)
_SOURCE_SUFFIXES = {".py", ".sh", ".yaml", ".toml", ".rs", ".lock"}
_SOURCE_EXCLUDED_PARTS = {
    ".git", ".pytest_cache", ".snakemake", "__pycache__", "conda", "target", "vendor"
}


def _pipeline_source_identity(root: Path) -> Mapping[str, object]:
    records: list[dict[str, object]] = []
    for top in _SOURCE_ROOTS:
        base = root / top
        if not base.exists():
            continue
        for path in sorted(base.rglob("*")):
            if not path.is_file() or any(part in _SOURCE_EXCLUDED_PARTS for part in path.parts):
                continue
            if path.suffix.lower() not in _SOURCE_SUFFIXES and path.name != "Snakefile":
                continue
            records.append({
                "path": str(path.relative_to(root)),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            })
    design = root / DROPLET_DESIGN_RESOURCE
    if design.is_file():
        records.append({"path": DROPLET_DESIGN_RESOURCE, "size_bytes": design.stat().st_size,
                        "sha256": sha256_file(design)})
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
    policy = get_backend(dict(config))
    if policy.name == "bismark":
        public = Path(f"{reference}.bismark")
        if not public.is_dir():
            raise FileNotFoundError(f"Bismark reference directory is missing or invalid: {public}")
        resolved = public.resolve(strict=True)
        genome = resolved / "genome.fa"
        if not genome.is_file() or genome.stat().st_size == 0:
            raise FileNotFoundError(f"Incomplete Bismark reference generation: {resolved}")
        local_alignment = bool(config["bismark"]["local_alignment"])
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

    indexes = backend_active_index_paths(reference, policy.name)
    missing = [path for path in indexes if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError(
            f"Incomplete {policy.name} reference generation: "
            + ", ".join(str(path) for path in missing)
        )
    return {
        "backend": policy.name,
        "indexes": [stable_file_identity(path, hash_content=False) for path in indexes],
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
    config = load_project_config(config_path)
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
        backend.name, {"biscuit": biscuit_env, "bismark": bismark_env,
                       "rastair": rastair_env}[backend.name]
    )
    protocol = str(config.get("analysis", {}).get("protocol", "")).strip().lower()
    return {
        "pipeline": {
            "root": str(pipeline_root.resolve()),
            "source": _pipeline_source_identity(pipeline_root),
        },
        "analysis": {
            "species": species,
            "protocol": protocol,
            "methylation_backend": backend.name,
        },
        "project": {
            "path": str(project_dir.resolve()),
            "config": stable_file_identity(config_path, hash_content=True),
            "sample_manifest": stable_file_identity(sample_manifest, hash_content=True),
            **({} if protocol == "droplet" else {"barcode_map": stable_file_identity(barcode_map, hash_content=True)}),
        },
        "raw_fastqs": _fastq_inventory(
            project_dir / "01_raw",
            read_project_sample_manifest(sample_manifest, species=species, protocol=protocol),
            protocol=protocol,
        ),
        "reference": reference_record,
        "backend_reference": _backend_reference_observation(config, reference),
        "environments": {
            "control": control_observation,
            backend.name: active_backend_observation,
        },
        "executables": {
            "demux_rs": _executable_identity(
                pipeline_root / "conda" / "rust" / "bin" / "demux_rs"
            ),
        },
    }


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


def _copy_observed_file(
    source: Path, destination: Path, expected: Mapping[str, object], input_store: Path
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
    copied = stable_file_identity(shared, hash_content=True)
    for key in ("size_bytes", "sha256"):
        if copied.get(key) != expected.get(key):
            raise RuntimeError(f"snapshot input store identity changed: {shared}")
    destination.symlink_to(os.path.relpath(shared, destination.parent))
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

    stage = Path(tempfile.mkdtemp(prefix=f".{snapshot_id}.", dir=snapshots_root))
    try:
        inputs_dir = stage / "inputs"
        project_basis = basis["project"]
        copies = {
            "config": _copy_observed_file(
                Path(str(project_basis["config"]["path"])), inputs_dir / "config.yaml",
                project_basis["config"], input_store
            ),
            "sample_manifest": _copy_observed_file(
                Path(str(project_basis["sample_manifest"]["path"])),
                inputs_dir / "sample_manifest.tsv", project_basis["sample_manifest"], input_store,
            ),
            **({} if basis["analysis"]["protocol"] == "droplet" else {"barcode_map": _copy_observed_file(
                Path(str(project_basis["barcode_map"]["path"])),
                inputs_dir / "Barcode_Map.csv", project_basis["barcode_map"], input_store,
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


def snapshot_input_path(snapshot: str | Path, key: str) -> Path:
    """从 snapshot 返回 config/sample_manifest/barcode_map 的不可变副本路径。"""
    if key not in {"config", "sample_manifest", "barcode_map"}:
        raise ValueError(f"unsupported snapshot input key: {key}")
    manifest = Path(snapshot).expanduser().resolve()
    payload = load_run_snapshot(manifest)
    root = manifest.parent
    record = payload.get("input_copies", {}).get(key)
    if not isinstance(record, Mapping):
        raise ValueError(f"run snapshot {key} copy is missing")
    return root / "inputs" / str(record["snapshot_path"])


RUN_MANIFEST_SCHEMA_VERSION = 8
DELIVERY_READY_SCHEMA_VERSION = 2
CPG_REPRESENTATION = "snapatac2_import_values_0based_v1"

_DELIVERY_DIR_NAMES = frozenset(PUBLIC_RESULT_DIR_NAMES)


def _utc_now() -> str:
    return dt.datetime.now(dt.UTC).isoformat()


def _read_json(path: Path, label: str) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid {label} JSON: {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{label} must contain a JSON object: {path}")
    return value


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


def _delivery_id(value: object) -> str:
    text = str(value)
    if len(text) != 64 or any(character not in "0123456789abcdef" for character in text):
        raise ValueError(f"invalid delivery id: {text!r}")
    return text


def _remove_path(path: Path) -> None:
    if not path.exists() and not path.is_symlink():
        return
    if path.is_dir() and not path.is_symlink():
        _make_writable(path)
        shutil.rmtree(path)
    else:
        path.unlink(missing_ok=True)


def _output_relative_is_owned(relative: Path) -> bool:
    return relative.parts[0] in _DELIVERY_DIR_NAMES


def _check_result_directories(root: Path, relative: Path, *, label: str) -> None:
    directory = root
    for part in (*relative.parts[:-1], None):
        if directory.is_symlink() or (directory.exists() and not directory.is_dir()):
            raise ValueError(f"{label} requires real result directories: {directory}")
        if part is not None:
            directory /= part


def _staged_output_inventory(results: Path) -> list[dict[str, object]]:
    manifest = Path("QC_Results/sample_manifest.tsv")
    declared = {manifest, Path("QC_Results/multiqc_report.html")}
    _check_result_directories(results, manifest, label="staged delivery")
    if (results / manifest).is_symlink():
        raise ValueError(f"staged result may not be a symlink: {results / manifest}")
    with (results / manifest).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "cpg_path", "snp_path"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("final sample manifest is missing current delivery path fields")
        for row in reader:
            for field, directory, suffix in (
                ("cpg_path", "CpG", ".cpg.tsv.zst"),
                ("snp_path", "SNP", ".snps.bed.gz"),
            ):
                value = str(row.get(field) or "").strip()
                if not value:
                    continue
                sample_id = str(row.get("sample_id") or "").strip()
                if SAMPLE_ID_RE.fullmatch(sample_id) is None:
                    raise ValueError(f"invalid delivery sample_id: {sample_id!r}")
                relative = Path(directory) / f"{sample_id}{suffix}"
                if value != f"03_results/{relative.as_posix()}":
                    raise ValueError(f"final sample manifest has an invalid {field}: {value!r}")
                declared.add(relative)
                if field == "snp_path":
                    declared.add(relative.with_name(relative.name + ".tbi"))

    records: list[dict[str, object]] = []
    for relative in sorted(declared):
        _check_result_directories(results, relative, label="staged delivery")
        path = results / relative
        if path.is_symlink():
            raise ValueError(f"staged result may not be a symlink: {path}")
        if not path.is_file():
            raise FileNotFoundError(f"declared staged result is missing: {path}")
        stat = path.stat()
        records.append(
            {
                "path": relative.as_posix(),
                "size_bytes": stat.st_size,
                "sha256": sha256_file(path),
            }
        )
    return records


def _validate_inventory(
    root: Path,
    records: object,
    *,
    label: str,
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
    stage_results = Path(staged_results_dir).expanduser().resolve()
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


    stage = (project / "02_work" / "results_stage").resolve(strict=False)
    delivery_id = _delivery_id(payload.get("delivery_id", ""))
    expected_root = project / "05_tmp" / "finalize" / delivery_id
    if manifest.parent != expected_root.resolve(strict=False):
        raise ValueError("delivery-ready location and delivery_id disagree")
    if validate_staged:
        _validate_inventory(stage, payload.get("outputs"), label="staged delivery")
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
    tmp = destination.with_name(
        f".{destination.name}.tmp.{os.getpid()}.{uuid.uuid4().hex}"
    )
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


def _run_manifest_payload(ready: Mapping[str, object]) -> dict[str, object]:
    snapshot_path = Path(str(ready["run_snapshot_path"])).resolve()
    snapshot = load_run_snapshot(snapshot_path)
    basis = snapshot["basis"]
    if not isinstance(basis, Mapping):
        raise ValueError("run snapshot basis is missing")
    pipeline_record = basis.get("pipeline")
    analysis_record = basis.get("analysis")
    reference_record = basis.get("reference")
    if not all(isinstance(item, Mapping) for item in (pipeline_record, analysis_record, reference_record)):
        raise ValueError("run snapshot provenance identity is incomplete")
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
            "cpg_representation": CPG_REPRESENTATION,
            "snapshot_id": snapshot["snapshot_id"],
            "snapshot_path": str(snapshot_path),
            "reference_path": str(Path(str(reference_record.get("resolved_path", ""))).resolve()),
        },
        "outputs": [dict(item) for item in ready["outputs"]],
    }


def publish_delivery(ready_path: str | Path) -> Path:
    """发布一笔 sealed 结果树到真实目录 03_results/。

    事务语义是 invalidate-first / commit-last：任何公开文件发生变化之前先删除
    旧 run_manifest.json，安装并校验全部公开文件之后才写回新 marker——marker 的
    存在即代表一代完整交付，中断只会留下无 marker 的待重发状态，不会留下
    「混合代结果 + complete marker」。交付内容全部来自 stable stage 的硬链接，
    逐文件原子落地（os.replace）；发布中断不损坏 stable stage，重新运行
    launcher 会直接重走本流程。生产调用必须位于 Snakemake 持有项目锁的
    onsuccess 回调内，锁在发布与 demux 清理完成后释放。
    """
    ready_manifest = Path(ready_path).expanduser().resolve()
    ready = load_delivery_ready(ready_manifest, validate_staged=False)
    project = Path(str(ready["project_dir"])).resolve()
    public = project / _PUBLIC_DIRNAME
    stage = (project / "02_work" / "results_stage").resolve(strict=False)

    outputs = _validate_inventory(stage, ready.get("outputs"), label="staged delivery")
    payload = _run_manifest_payload(ready)
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
    _validate_inventory(public, outputs, label="published delivery")

    atomic_write_json(marker, payload)

    finalize_root = project / "05_tmp" / "finalize"
    if ready_manifest.parent.parent == finalize_root:
        _remove_path(ready_manifest.parent)
    clean_demux_scratch_after_publish(project)
    return marker


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
                        if config["retention"]["keep_final_bam"]:
                            layout = ProjectLayout(
                                project / "02_work", public_manifest.parent,
                                get_backend(config).name,
                            )
                            with (public_manifest.parent / "QC_Results/sample_manifest.tsv").open() as handle:
                                for row in csv.DictReader(handle, delimiter="\t"):
                                    if row["analysis_status"] != "PASS":
                                        continue
                                    for path in layout.cell(row["sample_id"]).retained_bam_files:
                                        if not path.is_file():
                                            return {
                                                "status": "resumable",
                                                "reason": f"requested retained BAM/index is missing: {path}",
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
    """将文件组 staged 为 .partial 并校验；替换前撤销旧 commit_last，全部完成后最后提交它。"""
    destinations = [destination for _source, destination in pairs]
    if len(destinations) != len(set(destinations)):
        raise ValueError("duplicate destination in atomic publish bundle")
    if commit_last is not None and commit_last not in destinations:
        raise ValueError("--commit-last must match one --file destination")
    if bam_destination is not None and bam_destination not in destinations:
        raise ValueError("--bam-destination must match one --file destination")

    job_token = os.environ.get("SLURM_JOB_ID", "local")
    partials: dict[Path, Path] = {}
    try:
        for source, destination in pairs:
            if not source.is_file():
                raise FileNotFoundError(f"publish source is missing: {source}")
            destination.parent.mkdir(parents=True, exist_ok=True)


            for orphan in destination.parent.glob(f"{destination.name}.partial.*"):
                if orphan.is_file() or orphan.is_symlink():
                    orphan.unlink()

            partial = destination.with_name(
                f"{destination.name}.partial.{job_token}.{os.getpid()}"
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


class DemuxRow(TypedDict):
    demux_sample: str
    project_sample_id: str
    raw_sample: str
    pipeline_mode: str
    route_id: str
    demux_mode: str
    downstream_dna: str
    demux_report_schema_version: str
    demux_build_version: str
    demux_source_revision: str
    retention_policy: str
    retention_threshold: str
    rna_output_role: str
    high_cph_role: str
    plate_id: str
    cell_order: str
    dna_barcode: str
    rna_barcode: str
    dna_reads: str
    rna_reads: str
    input_fastq_pairs: str
    input_read_pairs: str
    usable_read_pairs: str
    matched_reads: str
    unmatched_reads: str
    ambiguous_reads: str
    short_reads: str
    unassigned_rate: str
    dna_r1: str
    dna_r2: str
    rna_r1: str
    rna_r2: str
    dna_status: str
    rna_status: str


class ProjectSampleRow(TypedDict):
    sample_id: str
    dna_raw_sample: str
    rna_sample: str
    species: str
    protocol: str
    notes: str


class CellIdentityRow(ProjectSampleRow):
    project_sample_id: str
    plate_id: str
    cell_order: str
    dna_barcode: str
    rna_barcode: str


class BarcodeRow(TypedDict):
    dna_barcode: str
    rna_barcode: str
    plate_id: str
    cell_order: str


class FileIdentity(TypedDict):
    path: str
    size_bytes: int
    sha256: str


class DemuxCompletion(TypedDict, total=False):
    schema_version: int
    status: str
    generation_id: str
    run_id: str
    project_sample_id: str
    cell_count: int
    cell_ids_sha256: str
    manifest: FileIdentity
    report: FileIdentity
    multiqc: FileIdentity
    initial_cell_manifest: FileIdentity
    dna_dir: str
    rna_dir: str
    _manifest_path: str
    _cell_manifest_path: str
    _rows: list[DemuxRow]


ManifestIndex = dict[str, list[tuple[Path, DemuxRow]]]


class FastqPairFiles(TypedDict):
    r1: list[str]
    r2: list[str]


class MappedDemuxRun(FastqPairFiles):
    raw_sample: str
    sample_id: str
    route_id: str
    source_column: str


DiscoveredSamples = dict[str, FastqPairFiles]
MappedDemuxRuns = dict[str, MappedDemuxRun]


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
) -> None:
    """先在最终路径旁写临时文件，再以一次原子 rename 发布小体积 manifest。"""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    temporary.write_text(render_tsv(rows, fieldnames), encoding="utf-8")
    temporary.replace(output)


def file_identity(path: str | Path) -> FileIdentity:
    source = Path(path).resolve()
    if not source.is_file():
        raise FileNotFoundError(f"Committed demux file is missing: {source}")
    return {
        "path": str(source),
        "size_bytes": source.stat().st_size,
        "sha256": sha256_file(source),
    }


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


def _object(
    value: object, fields: set[str], label: str, source: Path
) -> Mapping[str, object]:
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a JSON object: {source}")
    missing, extra = sorted(fields - value.keys()), sorted(value.keys() - fields)
    if missing or extra:
        detail = "; ".join(
            part for part in (
                "missing=" + ",".join(missing) if missing else "",
                "unsupported=" + ",".join(extra) if extra else "",
            ) if part
        )
        raise ValueError(f"{label} fields do not match schema 3 ({detail}): {source}")
    return value


def _integer(value: object, label: str, source: Path) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{label} must be a non-negative integer: {source}")
    return value


def _text(value: object, label: str, source: Path, *, empty: bool = False) -> str:
    if not isinstance(value, str) or (not empty and not value):
        raise ValueError(f"{label} must be a{' possibly empty' if empty else ' non-empty'} string: {source}")
    return value


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


DEMUX_CHUNK_READ_PAIRS = 50_000_000
DEMUX_CHUNK_RAW_BYTES = 16 * 1024 ** 3


def _rapidgzip_command(*arguments: str) -> list[str]:
    return [sys.executable, "-c", "import rapidgzip; raise SystemExit(rapidgzip.cli())", *map(str, arguments)]


def _read_gzip_line_index(path: Path) -> tuple[list[int], list[int], list[int], int, int]:
    lines, offsets, positions = [], [], []
    with path.open("rb", buffering=8 * 1024**2) as handle:
        if handle.read(20) != b"\0" * 8 + b"gzipindX" + b"\0" * 4:
            raise ValueError(f"invalid gzip line index: {path}")
        count, complete = struct.unpack(">QQ", handle.read(16))
        if count != complete or count > path.stat().st_size // 32:
            raise ValueError(f"incomplete gzip line index: {path}")
        previous_compressed = -1
        for _ in range(count):
            positions.append(handle.tell())
            offset, compressed, bits, window = struct.unpack(">QQII", handle.read(24))
            if bits > 7 or compressed * 8 - bits <= previous_compressed:
                raise ValueError(f"invalid gzip seek point: {path}")
            previous_compressed = compressed * 8 - bits
            handle.seek(window, os.SEEK_CUR)
            line = struct.unpack(">Q", handle.read(8))[0] - 1
            if line < 0 or (lines and (line < lines[-1] or offset < offsets[-1])):
                raise ValueError(f"unordered gzip line index: {path}")
            lines.append(line)
            offsets.append(offset)
        positions.append(handle.tell())
        raw_bytes, newlines = struct.unpack(">QQ", handle.read(16))
        if handle.read(1) or (lines and (newlines < lines[-1] or raw_bytes < offsets[-1])):
            raise ValueError(f"invalid gzip line index footer: {path}")
    if not lines or lines[0] != 0 or offsets[0] != 0:
        raise ValueError(f"gzip line index lacks the first record: {path}")
    lines.append(newlines)
    offsets.append(raw_bytes)
    return lines, offsets, positions, raw_bytes, newlines


def _write_gzip_range_index(source: Path, destination: Path, lines: Sequence[int], positions: Sequence[int],
                            start_pair: int, read_pairs: int) -> None:
    begin = max(0, bisect.bisect_left(lines, start_pair * 4) - 2)
    end = min(len(positions) - 1, bisect.bisect_left(lines, (start_pair + read_pairs) * 4) + 3)
    with source.open("rb", buffering=64 * 1024) as origin, destination.open("xb") as output:
        output.write(origin.read(20))
        count = end - begin + int(begin > 0)
        output.write(struct.pack(">QQ", count, count))
        ranges = ([(positions[0], positions[1])] if begin else []) + [(positions[begin], positions[end])]
        for start, stop in ranges:
            origin.seek(start)
            remaining = stop - start
            while remaining:
                block = origin.read(min(remaining, 1024 ** 2))
                if not block:
                    raise ValueError(f"truncated gzip line index: {source}")
                output.write(block)
                remaining -= len(block)
        origin.seek(positions[-1])
        footer = origin.read(16)
        if len(footer) != 16:
            raise ValueError(f"truncated gzip line index footer: {source}")
        output.write(footer)


def _plan_demux_pair_chunks(data: Sequence[tuple[dict[str, Any], list[int], list[int], list[int]]],
                            destination: Path, source_number: int, chunk_offset: int) -> list[dict[str, Any]]:
    chunks = []
    pairs = data[0][0]["read_pairs"]
    start = 0
    while start < pairs:
        end = min(pairs, start + DEMUX_CHUNK_READ_PAIRS)
        lower = []
        for metadata, lines, offsets, positions in data:
            before = max(0, bisect.bisect_left(lines, start * 4) - 1)
            lower.append(offsets[before])
            limit = bisect.bisect_right(offsets, offsets[before] + DEMUX_CHUNK_RAW_BYTES // 2) - 1
            end = min(end, lines[limit] // 4)
        if end <= start:
            raise ValueError("FASTQ records or gzip index spacing exceed the chunk scratch bound")
        upper = sum(offsets[bisect.bisect_left(lines, end * 4)] - lo
                    for (_, lines, offsets, _), lo in zip(data, lower, strict=True))
        entry = dict(id=f"{chunk_offset + len(chunks):06d}", source=source_number, start_pair=start,
                     read_pairs=end-start, raw_bytes=upper)
        for mate, (metadata, lines, offsets, positions) in enumerate(data, 1):
            name = f"{entry['id']}.R{mate}.gzi"
            _write_gzip_range_index(destination / metadata["index"], destination / name,
                                    lines, positions, start, end - start)
            entry[f"r{mate}_index"] = name
        chunks.append(entry)
        start = end
    for metadata, _, _, _ in data:
        (destination / metadata.pop("index")).unlink()
    return chunks


def index_demux_fastqs(r1: Sequence[str], r2: Sequence[str], output: str | Path, *, threads: int) -> None:
    """顺序读取压缩流并行校验 gzip，裁剪区段索引；不复制或重压缩 FASTQ。"""
    if not r1 or len(r1) != len(r2) or threads < 1:
        raise ValueError("demux index requires equally sized, non-empty R1/R2 lists and positive threads")
    for value in (*r1, *r2):
        source = Path(value)
        if not source.name.lower().endswith((".fastq.gz", ".fq.gz")):
            raise ValueError(f"demux index requires a gzip FASTQ suffix: {source}")
        with source.open("rb") as handle:
            if handle.read(2) != b"\x1f\x8b":
                raise ValueError(f"demux index input lacks gzip magic: {source}")
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.mkdir()
    committed = False
    inputs, chunks = [], []
    decoders = max(2, (int(threads) - 2) // 2)

    def build_index(item: tuple[int, int, str]) -> tuple[dict[str, Any], list[int], list[int], list[int]]:
        number, mate, value = item
        path = Path(value).resolve()
        before = path.stat()
        name = f"source{number:06d}.R{mate}.gzi"
        index = destination / name
        subprocess.run(_rapidgzip_command("-P", str(decoders), "--io-read-method", "sequential", "--verify", "--export-index", str(index),
                       "--index-format", "gztool-with-lines", str(path)), check=True)
        lines, offsets, positions, raw_bytes, newlines = _read_gzip_line_index(index)
        if raw_bytes:
            tail = subprocess.check_output(_rapidgzip_command("-P", "2", "--import-index", str(index),
                "-d", "-c", "--ranges", f"1@{raw_bytes - 1}", str(path)))
            if len(tail) != 1:
                raise ValueError(f"cannot verify final FASTQ line: {path}")
            newlines += int(tail != b"\n")
            lines[-1] = newlines
        if newlines % 4:
            raise ValueError(f"FASTQ does not contain complete four-line records: {path}")
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise ValueError(f"raw FASTQ changed while indexing: {path}")
        return dict(path=str(path), bytes=before.st_size, mtime_ns=before.st_mtime_ns,
                    index=name, raw_bytes=raw_bytes, read_pairs=newlines // 4), lines, offsets, positions

    try:
        with ThreadPoolExecutor(max_workers=2) as pool:
            for number, pair in enumerate(zip(r1, r2, strict=True)):
                data = list(pool.map(build_index, [(number, mate, value) for mate, value in enumerate(pair, 1)]))
                if data[0][0]["read_pairs"] != data[1][0]["read_pairs"]:
                    raise ValueError("paired FASTQ inputs are not in sync")
                inputs.append({"r1": data[0][0], "r2": data[1][0]})
                chunks.extend(_plan_demux_pair_chunks(data, destination, number, len(chunks)))
        if not chunks:
            raise ValueError("demux index found no complete FASTQ records")
        manifest = dict(schema_version=1, input_fastq_pairs=len(inputs), inputs=inputs, chunks=chunks,
                        input_read_pairs=sum(row["read_pairs"] for row in chunks))
        atomic_write_json(destination / "manifest.json", manifest)
        read_demux_chunks(destination)
        committed = True
        print(f"Demux index: {manifest['input_read_pairs']} pairs, {len(chunks)} independent chunks", flush=True)
    finally:
        if not committed:
            shutil.rmtree(destination)


def read_demux_chunks(directory: str | Path) -> dict[str, Any]:
    """加载索引分块清单，核对来源身份、连续区段、scratch 上限和全库配对守恒。"""
    source = Path(directory) / "manifest.json"
    data = json.loads(source.read_text(encoding="utf-8"))
    if not isinstance(data, dict) or data.get("schema_version") != 1:
        raise ValueError(f"invalid demux index manifest: {source}")
    inputs, chunks = data.get("inputs"), data.get("chunks")
    if not isinstance(inputs, list) or not inputs or data.get("input_fastq_pairs") != len(inputs):
        raise ValueError(f"demux index lacks physical FASTQ pairs: {source}")
    for number, pair in enumerate(inputs):
        for mate in (1, 2):
            row = pair.get(f"r{mate}", {})
            if not isinstance(row.get("path"), str) or not Path(row["path"]).is_absolute():
                raise ValueError(f"invalid demux index source: {source}")
            for field in ("bytes", "mtime_ns", "raw_bytes", "read_pairs"):
                _integer(row.get(field), f"input.{field}", source)
        if pair["r1"]["read_pairs"] != pair["r2"]["read_pairs"]:
            raise ValueError(f"unpaired demux index source: {source}")
    if not isinstance(chunks, list) or not chunks:
        raise ValueError(f"demux index has no chunks: {source}")
    covered = [0] * len(inputs)
    previous = 0
    for index, row in enumerate(chunks):
        if not isinstance(row, dict) or row.get("id") != f"{index:06d}":
            raise ValueError(f"demux chunk order/identity is invalid: {source}")
        if any(row.get(f"r{mate}_index") != f"{index:06d}.R{mate}.gzi" for mate in (1, 2)):
            raise ValueError(f"invalid demux chunk index: {source}")
        for field in ("read_pairs", "raw_bytes", "source", "start_pair"):
            _integer(row.get(field), f"chunk.{field}", source)
        number = row["source"]
        if number < previous or number >= len(inputs) or row["start_pair"] != covered[number]:
            raise ValueError(f"demux index has overlapping or missing ranges: {source}")
        if row["raw_bytes"] > DEMUX_CHUNK_RAW_BYTES or row["read_pairs"] > DEMUX_CHUNK_READ_PAIRS:
            raise ValueError(f"demux chunk exceeds its scratch bound: {source}")
        if not row["read_pairs"] and (len(chunks) != 1 or row["raw_bytes"]):
            raise ValueError(f"unexpected empty demux chunk: {source}")
        covered[number] += row["read_pairs"]
        previous = number
    if covered != [pair["r1"]["read_pairs"] for pair in inputs] or _integer(data.get("input_read_pairs"), "input_read_pairs", source) != sum(covered):
        raise ValueError(f"demux index pair accounting mismatch: {source}")
    return data


def _run_indexed_demux(index_dir: str | Path, chunk: str, command: list[str]) -> dict[str, Any]:
    plan = read_demux_chunks(index_dir)
    if not re.fullmatch(r"[0-9]{6}", chunk) or int(chunk) >= len(plan["chunks"]):
        raise ValueError(f"unknown demux chunk: {chunk!r}")
    entry = plan["chunks"][int(chunk)]
    sources = plan["inputs"][entry["source"]]
    readers, read_fds = [], []
    consumer = None
    try:
        command = [*command, "--fastq-stream"]
        for mate in (1, 2):
            row = sources[f"r{mate}"]
            path = Path(row["path"])
            state = path.stat()
            if (state.st_size, state.st_mtime_ns) != (row["bytes"], row["mtime_ns"]):
                raise ValueError(f"raw FASTQ changed since indexing: {path}")
            read_fd, write_fd = os.pipe()
            read_fds.append(read_fd)
            try:
                reader = subprocess.Popen(_rapidgzip_command("-P", "1", "--import-index", str(Path(index_dir) / entry[f"r{mate}_index"]),
                    "-d", "-c", "--ranges", f"{entry['read_pairs']*4}L@{entry['start_pair']*4}L", str(path)), stdout=write_fd)
                readers.append(reader)
            finally:
                os.close(write_fd)
            command.extend((f"--r{mate}", f"/dev/fd/{read_fd}"))
        consumer = subprocess.Popen(command, pass_fds=tuple(read_fds))
        for fd in read_fds:
            os.close(fd)
        read_fds.clear()
        if consumer.wait() != 0:
            raise subprocess.CalledProcessError(consumer.returncode, command)
        for reader in readers:
            if reader.wait() != 0:
                raise subprocess.CalledProcessError(reader.returncode, reader.args)
    finally:
        for fd in read_fds:
            os.close(fd)
        for process in [consumer, *readers]:
            if process is not None and process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
    return entry


def run_demux_chunk(index_dir: str | Path, chunk: str, *, binary: str, barcode_map: str,
                    sample: str, scratch: str, pack: str, metadata: str, threads: int) -> None:
    """在 job scratch 解复用一个配对块，保留非零细胞并打包 gzip；metadata 最后提交。"""
    work = Path(scratch)
    work.mkdir(parents=True, exist_ok=True)
    dna = work / "DNA"
    report_file = work / "report.json"
    if dna.exists() or report_file.exists():
        raise FileExistsError("demux chunk requires an unused job scratch directory")
    command = [str(binary), "--barcode-map", str(barcode_map), "--sample-name=" + sample,
               "--dna-out-dir", str(dna), "--rna-out-dir", str(work / "RNA"),
               "--json-report", str(report_file), "--mode", "dna-only-droplet",
               "--min-matched-read-pairs", "1", "--threads", str(max(1, int(threads) - 4))]
    entry = _run_indexed_demux(index_dir, chunk, command)
    report = read_demux_report(report_file)
    if report["input_read_pairs"] != entry["read_pairs"]:
        raise ValueError(f"demux chunk {chunk} lost or duplicated input pairs")
    destination, commit = Path(pack), Path(metadata)
    destination.parent.mkdir(parents=True, exist_ok=True)
    commit.parent.mkdir(parents=True, exist_ok=True)
    token = ".tmp." + uuid.uuid4().hex
    pack_tmp, metadata_tmp = Path(str(destination) + token), Path(str(commit) + token)
    members = []
    try:
        with pack_tmp.open("xb", buffering=1024**2) as writer:
            for row in report["samples"]:
                member = []
                for mate in (1, 2):
                    offset = writer.tell()
                    if row["dna_read_count"]:
                        path = dna / row["cell_sample_id"] / f"{row['cell_sample_id']}_R{mate}.fastq.gz"
                        with path.open("rb") as reader:
                            shutil.copyfileobj(reader, writer, length=1024**2)
                        if writer.tell() - offset != path.stat().st_size or writer.tell() == offset:
                            raise ValueError(f"incomplete chunk output: {path}")
                    member.extend((offset, writer.tell() - offset))
                members.append(member)
            size = writer.tell()
            writer.flush()
            os.fsync(writer.fileno())
        metadata_tmp.write_text(json.dumps({"schema_version": 1, "chunk": chunk, "pack_bytes": size,
                                           "report": report, "members": members}, separators=(",", ":")) + "\n", encoding="utf-8")
        commit.unlink(missing_ok=True)
        os.replace(pack_tmp, destination)
        os.replace(metadata_tmp, commit)
    finally:
        pack_tmp.unlink(missing_ok=True)
        metadata_tmp.unlink(missing_ok=True)


def merge_demux_chunks(index_dir: str | Path, packs_dir: str | Path, dna_dir: str | Path,
                       report_path: str | Path, *, threshold: int, threads: int) -> None:
    """按细胞组连续读取 pack 并保持原始块序；受句柄预算约束，关闭全部输出后提交报告。"""
    import resource

    soft_limit = resource.getrlimit(resource.RLIMIT_NOFILE)[0]
    fd_budget = max(1, (1024 if soft_limit == resource.RLIM_INFINITY else soft_limit) - 64)
    pack_batch_size = min(64, max(1, fd_budget // 2))
    workers = min(8, int(threads), max(1, (fd_budget - pack_batch_size) // 2))
    cells_per_group = min(32, max(1, (fd_budget - pack_batch_size) // (2 * max(1, workers))))
    plan = read_demux_chunks(index_dir)
    packs, destination = Path(packs_dir), Path(dna_dir)
    if threshold < 0 or threads < 1:
        raise ValueError("demux merge requires non-negative threshold and positive threads")
    destination.mkdir(parents=True, exist_ok=True)
    if any(destination.iterdir()):
        raise FileExistsError(f"demux merge destination is not empty: {destination}")
    aggregate = None
    identity = None
    for entry in plan["chunks"]:
        source = packs / (entry["id"] + ".json")
        data = json.loads(source.read_text(encoding="utf-8"))
        if data.get("schema_version") != 1 or data.get("chunk") != entry["id"]:
            raise ValueError(f"demux chunk metadata identity mismatch: {source}")
        report = _parse_demux_report(data.get("report"), source)
        if report["mode"] != "dna-only-droplet" or report["retention_threshold"] != 1 or report["input_fastq_pairs"] != 1:
            raise ValueError(f"demux chunk used a wrong protocol or per-chunk filter: {source}")
        if report["input_read_pairs"] != entry["read_pairs"]:
            raise ValueError(f"demux chunk read-pair accounting mismatch: {source}")
        fixed = {key: report[key] for key in ("schema_version", "build_version", "source_revision", "mode", "retention_policy")}
        fixed["samples"] = [{key: value for key, value in row.items()
                              if key not in {"dna_read_count", "rna_read_count", "dna_status", "rna_status"}}
                             for row in report["samples"]]
        if identity is not None and identity != fixed:
            raise ValueError(f"demux chunks disagree on implementation or cell identities: {source}")
        identity = fixed
        members = data.get("members")
        if not isinstance(members, list) or len(members) != len(report["samples"]):
            raise ValueError(f"demux pack member count mismatch: {source}")
        offset = 0
        for row, member in zip(report["samples"], members, strict=True):
            if not isinstance(member, list) or len(member) != 4:
                raise ValueError(f"invalid demux pack member: {source}")
            for position, length in (member[:2], member[2:]):
                _integer(position, "member.offset", source)
                _integer(length, "member.length", source)
                if position != offset or bool(length) != bool(row["dna_read_count"]):
                    raise ValueError(f"demux pack offset/count mismatch: {source}")
                offset += length
        if offset != _integer(data.get("pack_bytes"), "pack_bytes", source) or offset != (packs / (entry["id"] + ".pack")).stat().st_size:
            raise ValueError(f"demux pack is truncated or has trailing bytes: {source}")
        if aggregate is None:
            aggregate = copy.deepcopy(report)
        else:
            for field in ("input_read_pairs", "usable_read_pairs", "matched_reads"):
                aggregate[field] += report[field]
            for field in _FATE_FIELDS:
                aggregate["read_fates"][field] += report["read_fates"][field]
            for total, row in zip(aggregate["samples"], report["samples"], strict=True):
                total["dna_read_count"] += row["dna_read_count"]
    assert aggregate is not None
    aggregate["retention_threshold"] = threshold
    aggregate["input_fastq_pairs"] = plan["input_fastq_pairs"]
    for row in aggregate["samples"]:
        row["dna_status"] = cell_status(row["dna_read_count"], threshold)
    _parse_demux_report(aggregate, Path(report_path))
    if aggregate["input_read_pairs"] != plan["input_read_pairs"]:
        raise ValueError("merged demux input pairs do not match the complete index manifest")
    passing = [(index, row["cell_sample_id"]) for index, row in enumerate(aggregate["samples"]) if row["dna_status"] == "Pass"]
    for _, cell in passing:
        (destination / cell).mkdir()
    for batch_start in range(0, len(plan["chunks"]), pack_batch_size):
        batch = plan["chunks"][batch_start:batch_start + pack_batch_size]
        with ExitStack() as opened:
            inputs = []
            for entry in batch:
                metadata = json.loads((packs / (entry["id"] + ".json")).read_text(encoding="utf-8"))
                reader = opened.enter_context((packs / (entry["id"] + ".pack")).open("rb", buffering=0))
                inputs.append((reader.fileno(), metadata["members"]))

            def merge_cells(items: list[tuple[int, str]]) -> None:
                with ExitStack() as outputs:
                    writers = []
                    for index, cell in items:
                        pair = []
                        for mate in (1, 2):
                            target = destination / cell / f"{cell}_R{mate}.fastq.gz"
                            pair.append(outputs.enter_context(target.open(
                                "wb" if batch_start == 0 else "ab", buffering=1024**2)))
                        writers.append((index, cell, pair))
                    for fd, members in inputs:
                        for index, cell, pair in writers:
                            for mate, writer in enumerate(pair):
                                offset, remaining = members[index][mate*2:(mate+1)*2]
                                while remaining:
                                    data = os.pread(fd, min(1024**2, remaining), offset)
                                    if not data:
                                        raise ValueError(f"demux pack became truncated while merging {cell}")
                                    writer.write(data)
                                    offset += len(data)
                                    remaining -= len(data)
                    for _, _, pair in writers:
                        for writer in pair:
                            writer.flush()

            groups = [passing[i:i+cells_per_group] for i in range(0, len(passing), cells_per_group)]
            with ThreadPoolExecutor(max_workers=workers) as pool:
                for _ in pool.map(merge_cells, groups):
                    pass
        print(f"Demux merge: completed blocks {batch_start + 1}-{batch_start + len(batch)} / {len(plan['chunks'])}", flush=True)
    result = Path(report_path)
    result.parent.mkdir(parents=True, exist_ok=True)
    result.write_text(json.dumps(aggregate, indent=2) + "\n", encoding="utf-8")


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
        if protocol == "droplet" and row["rna_sample"]:
            raise ValueError("Droplet does not support RNA association")
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


@lru_cache(maxsize=1)
def droplet_design_barcodes() -> frozenset[str]:
    """读取固定内容身份的 DD-MET5 设计条码，供 calling 约束身份，不依赖 RNA 名单。"""
    path = Path(__file__).resolve().parent.parent / DROPLET_DESIGN_RESOURCE
    data = path.read_bytes()
    if sha256_bytes(data) != DROPLET_DESIGN_SHA256:
        raise ValueError("DD-MET5 design barcode resource SHA256 mismatch")
    return frozenset(gzip.decompress(data).decode("ascii").splitlines())


def call_droplet_cells(
    counts_path: str | Path, output: str | Path, metrics_path: str | Path,
    *, r1: Sequence[str | Path] = (), demux_binary: str | Path | None = None, threads: int = 1,
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
            if not r1 or demux_binary is None or threads < 1:
                raise ValueError("Droplet error screening requires raw R1 and demux binary for molecule evidence")
            selected = [{key: r[key] for key in ("barcode", "parent")} for r in candidates]
            with tempfile.TemporaryDirectory(prefix=".droplet_evidence.", dir=Path(metrics_path).parent) as tmp:
                pair_file, evidence_file = Path(tmp) / "pairs.json", Path(tmp) / "molecules.json"
                pair_file.write_text(json.dumps(selected))
                command = [str(demux_binary), "--mode", "dna-only-droplet", "--threads", str(threads),
                           "--evidence-pairs", str(pair_file), "--evidence-output", str(evidence_file)]
                for path in r1:
                    command.extend(["--r1", str(path)])
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
        temporary = destination.with_name(destination.name + "." + uuid.uuid4().hex + ".tmp")
        temporary.write_text(content)
        os.replace(temporary, destination)


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
) -> None:
    """按孔板 map 或逐文库 called cells 展开规范 cell 身份，写出初始 manifest。"""
    project_rows = read_project_sample_manifest(
        Path(project_manifest), species=species, protocol=protocol
    )
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
        try:
            order: tuple[int, object] = (0, int(raw_order))
        except ValueError:
            order = (1, raw_order)
        return (row.get("project_sample_id", ""), order, row.get("plate_id", ""))

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
    used_run_ids: set[str] = set()
    for raw_sample in sorted(expected_raw):
        identity = expected_raw[raw_sample]
        sid = identity["sample_id"]
        route_id = identity["route_id"]
        run_id = sid if identity["source_column"] == "dna_raw_sample" else f"{sid}__rna"
        if run_id in used_run_ids:
            raise ValueError(f"Canonical demux run ID collision: {run_id!r}")
        used_run_ids.add(run_id)
        mapped[run_id] = {
            "r1": discovered[raw_sample]["r1"],
            "r2": discovered[raw_sample]["r2"],
            "raw_sample": raw_sample,
            "sample_id": sid,
            "route_id": route_id,
            "source_column": identity["source_column"],
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
    *, calling_metrics: str = "", structure_metrics: str = "",
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
            temporary = destination.with_name(destination.name + f".tmp.{os.getpid()}")
            temporary.unlink(missing_ok=True)
            try:
                os.link(source, temporary)
                os.replace(temporary, destination)
            finally:
                temporary.unlink(missing_ok=True)
            collected += 1

    payload = {"run_id": run_id, "collected_files": collected}
    marker_path = Path(marker)
    marker_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = marker_path.with_name(marker_path.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.replace(temporary, marker_path)
    finally:
        temporary.unlink(missing_ok=True)
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
) -> ManifestIndex:
    """逐个校验已提交 generation 并按 cell ID 建立行索引。

    反复展开 DAG 的调用方必须自行负责 memoisation；本模块保持纯校验/索引层，不藏任何状态。
    """
    state_dir = Path(state_dir).resolve()
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

    index: ManifestIndex = {}
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
        for row in committed["_rows"]:
            sample = str(row.get("demux_sample", "") or "").strip()
            if not sample or sample == "__SUMMARY__":
                continue
            index.setdefault(sample, []).append((manifest, row))
    return index


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
    state_dir: str,
    demux_sample: str,
    read_pair: int,
    *,
    manifest_index: ManifestIndex,
) -> Path:
    """只从已提交的 manifest 解析 dna_status=Pass 的 DNA FASTQ。"""
    state_dir = Path(state_dir).resolve()
    base = Path(str(demux_sample)).name
    key = f"dna_r{int(read_pair)}"
    matches = [
        (manifest, row)
        for manifest, row in manifest_index.get(base, [])
        if str(row.get("downstream_dna", "")).strip().lower() == "true"
    ]

    if not matches:
        raise FileNotFoundError(
            f"No current demux manifest row found for DNA sample {base!r}. "
            "The checkpoint manifest is the sole authority; remove stale demux markers "
            "and rerun demultiplexing instead of guessing a path."
        )
    if len(matches) > 1:
        raise ValueError(
            f"DNA sample {base!r} appears in multiple demux manifests: "
            + ", ".join(str(item[0]) for item in matches)
        )

    manifest, row = matches[0]
    manifest = manifest.resolve()
    try:
        manifest.relative_to(state_dir)
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
        if label == "unmapped":
            unaligned += count
        else:
            try:
                mapq = int(label)
            except ValueError as exc:
                raise ValueError(f"Invalid BISCUIT MAPQ label in {path}: {label!r}") from exc
            if mapq >= 40:
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
            "Duplicate_Pair_Rate_pct": (
                100.0 * duplicate_pairs / accepted_pairs if accepted_pairs else 0.0
            ),
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
    layout = ProjectLayout(intermediate, results, opts["backend"])
    published_layout = ProjectLayout(intermediate, published_results, opts["backend"])
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
        backend=opts["backend"],
        publish_rna_bam=opts["publish_rna_bam"],
        generate_snp=opts["generate_snp"],
        keep_final_bam=opts["keep_final_bam"],
        bismark=bismark,
    )
    _emit_final_tsv(cells, output, pair_qc_multiqc_json, cutadapt_qc_multiqc_json)


def _parse_final_options(options: Mapping[str, object]) -> dict[str, object]:

    backend = str(options.get("methylation_backend", "")).strip().lower()
    if backend not in {"biscuit", "bismark", "rastair"}:
        raise ValueError("methylation_backend must be biscuit, bismark or rastair")
    option_fields = {
        "generate_snp",
        "keep_final_bam",
        "publish_rna_bam",
        "methylation_backend",
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
        if key == "methylation_backend":
            continue
        if not isinstance(value, bool):
            raise ValueError(f"Final manifest option {key} must be a boolean")
    return {
        "backend": backend,
        "generate_snp": active_options["generate_snp"],
        "keep_final_bam": active_options["keep_final_bam"],
        "publish_rna_bam": active_options["publish_rna_bam"],
    }


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
            expected_cell_manifest_path=initial_cell_manifest,
        )
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


def _bismark_version_text(bismark: str | None) -> str:
    if not bismark:
        return ""
    version = subprocess.run(
        [bismark, "--version"], capture_output=True, text=True, check=True
    )
    return (version.stdout + version.stderr).strip()


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
        if not any(rg.get("SM") == sample for rg in header.get("RG", [])):
            raise ValueError("TAPS BAM lacks the sample read group")
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
        if process.poll() is None:
            process.terminate()
        process.wait()
    atomic_write_json(output, {"sample": sample, "backend": "rastair", "schema_version": 1,
                               "read_groups": header.get("RG", []), "programs": header.get("PG", []),
                               "reference_sequences": header["SQ"], **counts})


def emit_rastair_cpg(bed: str, reference: str, alignment_qc: str, output_qc: str, stream: TextIO) -> None:
    """输出共用 CpG writer 的整数计数输入；CpG 任一端已调用为变异时排除两端，不使用 beta。"""
    import pysam

    qc = json.loads(Path(alignment_qc).read_text())
    counts = {"native_cpg_rows": 0, "canonical_cpg_rows": 0, "excluded_variant_or_denovo_rows": 0,
              "zero_observation_rows": 0}
    controls = {name: {"mod": 0, "unmod": 0, "sites": 0} for name in ("lambda", "pUC19")}
    required = {"#chr", "start", "end", "unmod", "mod", "coverage", "genotype", "cpg"}
    with pysam.FastaFile(reference) as ref, Path(bed).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not required.issubset(reader.fieldnames or ()) or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("Rastair BED lacks the native 2.2 header")
        previous = None
        contigs = {name: i for i, name in enumerate(ref.references)}
        pending_key = None
        pending = []

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
            if None in row or any(row.get(key) in (None, "") for key in required):
                raise ValueError(f"{label} has missing or extra fields")
            chrom = row["#chr"]
            if chrom not in contigs:
                raise ValueError(f"{label} contig {chrom!r} is absent from the reference FASTA")
            try:
                pos, end, mod, unmod, depth = (int(row[k]) for k in ("start", "end", "mod", "unmod", "coverage"))
            except ValueError as exc:
                raise ValueError(f"{label} has non-integer coordinates/counts") from exc
            if min(pos, mod, unmod, depth) < 0 or end != pos + 1 or mod + unmod > depth:
                raise ValueError(f"Invalid Rastair coordinate/counts: {chrom}:{pos}")
            key = (contigs[chrom], pos)
            if previous is not None and key <= previous:
                raise ValueError("Rastair BED is not ordered or contains duplicate positions")
            previous = key
            if row["cpg"] != "REF":
                counts["excluded_variant_or_denovo_rows"] += 1
                continue
            base = ref.fetch(chrom, pos, pos + 1).upper()
            start = pos if base == "C" else pos - 1
            if base not in {"C", "G"} or start < 0 or ref.fetch(chrom, start, start + 2).upper() != "CG":
                raise ValueError(f"Rastair reference CpG disagrees with FASTA: {chrom}:{pos}")
            cpg_key = (chrom, start)
            if cpg_key != pending_key:
                flush_cpg()
                pending = []
                pending_key = cpg_key
            pending.append((chrom, pos, mod, unmod, depth, row["genotype"], base))
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


def _taps_cell_metrics(job: Mapping[str, object], cutadapt: Mapping[str, object]) -> dict[str, object]:
    paths = job["paths"]
    qc = json.loads(Path(paths["rastair_qc"]).read_text())
    if qc["sample"] != job["cell_id"] or qc["backend"] != "rastair":
        raise ValueError("TAPS QC identity disagrees with the active cell")
    if qc["trimmed_pairs"] != cutadapt["trimmed_pairs"] or qc["bam_primary_records"] != 2 * qc["trimmed_pairs"]:
        raise ValueError("TAPS BAM and Cutadapt pair counts disagree")
    accepted, rejected, duplicates, postdedup = (int(qc[k]) for k in (
        "backend_accepted_pairs", "backend_rejected_pairs", "duplicate_pairs", "postdedup_pairs"))
    if accepted + rejected != qc["trimmed_pairs"] or postdedup + duplicates != accepted:
        raise ValueError("TAPS pair funnel does not close")
    if sum(qc[k] for k in ("unmapped_pairs", "discordant_pairs", "quality_rejected_pairs")) != rejected:
        raise ValueError("TAPS rejection categories do not close")
    if not Path(paths["cpg"]).is_file():
        raise FileNotFoundError(paths["cpg"])
    keep_bam = bool(job["keep_final_bam"])
    if keep_bam:
        for key in ("marked_bam", "marked_bam_index"):
            if not Path(paths[key]).is_file():
                raise FileNotFoundError(paths[key])
    return {**cutadapt, **{k: v for k, v in qc.items() if k.startswith(("bam_", "lambda_", "pUC19_"))},
            "methylation_backend": "rastair", "backend_qc_status": "PASS", "analysis_status": "PASS",
            "cpg_path": job["public"]["cpg"], "cpg_representation": CPG_REPRESENTATION,
            "snp_path": "", "rna_bam_path": "", "methylation_signal": "5mC+5hmC",
            "retained_bam_path": job["public"]["marked_bam"] if keep_bam else "",
            "retained_bam_index": job["public"]["marked_bam_index"] if keep_bam else "",
            "retained_bam_role": "taps_all_alignments_marked_duplicates_for_variant_reanalysis" if keep_bam else "",
            "mapping_policy": "bwa_rastair_eligible_pairs", "native_mapping_policy": "bwa_proper_primary_any_mate_mapq20",
            "native_mapping_unit": "read_pairs", "native_mapping_pct": 100.0 * accepted / qc["trimmed_pairs"] if qc["trimmed_pairs"] else 0.0,
            "both_primary_mapped_pairs": accepted, "one_primary_mapped_pairs": 0,
            "backend_accepted_pairs": accepted, "backend_rejected_pairs": rejected,
            "unmapped_pairs": qc["unmapped_pairs"], "ambiguous_pairs": "", "no_primary_pairs": "",
            "discordant_pairs": qc["discordant_pairs"], "quality_rejected_pairs": qc["quality_rejected_pairs"],
            "duplicate_pairs": duplicates, "postdedup_pairs": postdedup, "final_retained_pairs": postdedup,
            "dedup_policy": "samtools_markdup_flag_exclude_in_rastair", "pair_qc_metric_unit": "read_pairs",
            "high_cph_assessed_pairs": 0, "high_cph_flagged_pairs": 0, "high_cph_removed_pairs": 0,
            "high_cph_fraction": "nan", "high_cph_role": "not_applicable",
            "non_cpg_methylation_pct": "nan", "non_cpg_metric_source": "not_measured_taps_cpg_only"}


def _metrics_for_cell(job: Mapping[str, object]) -> tuple[str, dict[str, object]]:
    cell_id = str(job["cell_id"])
    record = dict(job["record"])
    backend = str(job["backend"])
    publish_rna_bam = bool(job["publish_rna_bam"])
    generate_snp = bool(job["generate_snp"])
    paths = {key: Path(str(value)) for key, value in job["paths"].items()}
    public = job["public"]
    protocol = str(record.get("protocol", "") or "").strip().lower()

    updates: dict[str, object] = {}
    canonical: dict[str, object] = {}
    cutadapt = _parse_cutadapt_metrics(paths["cutadapt_json"])
    try:
        demux_pairs = int(str(record.get("dna_reads", "")).strip())
    except ValueError as exc:
        raise ValueError(f"Active cell {cell_id!r} has invalid dna_reads") from exc
    if cutadapt["cutadapt_input_pairs"] != demux_pairs:
        raise ValueError(
            f"Cutadapt input disagrees with demux DNA pairs for {cell_id!r}: "
            f"cutadapt={cutadapt['cutadapt_input_pairs']}, demux={demux_pairs}"
        )
    trimmed_pairs = int(cutadapt["trimmed_pairs"])

    if backend == "rastair":
        return cell_id, _taps_cell_metrics(job, cutadapt)

    if not paths["high_cph_summary"].is_file():
        raise FileNotFoundError(f"Active cell is missing high-CpH summary: {paths['high_cph_summary']}")
    bsconv = _read_high_cph_summary(paths["high_cph_summary"], cell_id, backend)
    if bsconv["protocol"] != protocol:
        raise ValueError(
            f"High-CpH summary protocol disagrees with the cell identity: {paths['high_cph_summary']}"
        )
    if not bsconv["filter_applied"]:
        raise ValueError(
            f"High-CpH summary must record the fixed filter action: {paths['high_cph_summary']}"
        )

    if backend == "biscuit":
        if not paths["dupsifter_multiqc"].is_file():
            raise FileNotFoundError(
                f"Active cell is missing Dupsifter QC: {paths['dupsifter_multiqc']}"
            )
        _native, pair_counts = _read_dupsifter_metrics(
            paths["dupsifter_multiqc"], cell_id, trimmed_pairs
        )
        alignment_native = _read_biscuit_native_mapping(paths["mapq_table"])
        updates.update(alignment_native)
        canonical.update(pair_counts)
        canonical.update({
            "native_mapping_policy": alignment_native["native_mapping_policy"],
            "native_mapping_unit": alignment_native["native_mapping_unit"],
            "native_mapping_pct": alignment_native["native_mapping_pct"],
        })
        canonical["ambiguous_pairs"] = ""
    else:
        metrics = build_metrics(
            sample=cell_id,
            version=str(job["bismark_version"]),
            alignment_report=paths["bismark_align_report"],
            dedup_report=paths["bismark_dedup_report"],
            filter_summary=bsconv,
            protocol=record["protocol"],
        )
        total_pairs = int(metrics["total_pairs"])
        if total_pairs != trimmed_pairs:
            raise ValueError(
                f"Bismark report disagrees with Cutadapt output for {cell_id!r}: "
                f"bismark={total_pairs}, cutadapt={trimmed_pairs}"
            )
        accepted = int(metrics["accepted_pairs"])
        rejected = int(metrics["unmapped_pairs"]) + int(metrics["ambiguous_pairs"])
        if accepted + rejected != trimmed_pairs:
            raise ValueError(f"Bismark pair funnel is inconsistent for {cell_id!r}")
        canonical.update(
            {
                "native_mapping_policy": "bismark_unique_concordant_pairs",
                "native_mapping_unit": "read_pairs",
                "native_mapping_pct": (100.0 * accepted / trimmed_pairs if trimmed_pairs else 0.0),
                "both_primary_mapped_pairs": accepted,
                "one_primary_mapped_pairs": 0,
                "backend_accepted_pairs": accepted,
                "backend_rejected_pairs": rejected,
                "unmapped_pairs": int(metrics["unmapped_pairs"]),
                "ambiguous_pairs": int(metrics["ambiguous_pairs"]),
                "no_primary_pairs": "",
                "duplicate_pairs": int(metrics["duplicate_pairs"]),
                "postdedup_pairs": int(metrics["postdedup_pairs"]),
            }
        )

    assessed_pairs = int(bsconv["candidate_read_pairs"])
    flagged_pairs = int(bsconv["flagged_read_pairs"])
    excluded_pairs = int(bsconv["excluded_read_pairs"])
    removed_pairs = flagged_pairs
    high_cph_fraction = float(bsconv["high_cph_fraction"])
    postdedup_pairs = int(canonical["postdedup_pairs"])
    if assessed_pairs + excluded_pairs > postdedup_pairs:
        raise ValueError(
            f"High-CpH accounting exceeds post-dedup pairs for {cell_id!r}"
        )
    final_retained_pairs = postdedup_pairs - removed_pairs

    canonical.update(
        {
            "cutadapt_input_pairs": int(cutadapt["cutadapt_input_pairs"]),
            "trimmed_pairs": int(cutadapt["trimmed_pairs"]),
            "high_cph_assessed_pairs": assessed_pairs,
            "high_cph_flagged_pairs": flagged_pairs,
            "high_cph_removed_pairs": removed_pairs,
            "final_retained_pairs": final_retained_pairs,
            "mapping_policy": (
                "biscuit_any_primary_pair_funnel"
                if backend == "biscuit"
                else "bismark_unique_concordant"
            ),
            "dedup_policy": (
                "dupsifter_wgbs_signature_remove_dups"
                if backend == "biscuit"
                else DROPLET_DEDUP_POLICY if record["protocol"] == "droplet" else "bismark_paired_endpoint_orientation"
            ),
            "pair_qc_metric_unit": "read_pairs",
            "r1_adapter_pct": float(cutadapt["r1_adapter_pct"]),
            "r2_adapter_pct": float(cutadapt["r2_adapter_pct"]),
            "high_cph_fraction": high_cph_fraction,
        }
    )
    if backend == "biscuit":
        canonical["non_cpg_methylation_pct"] = _parse_cph_methylation_pct(
            paths["cph_table"]
        )
        canonical["non_cpg_metric_source"] = "biscuit_cph_retention_by_read_position"
    else:
        canonical["non_cpg_methylation_pct"] = (
            _parse_bismark_non_cpg_methylation_pct(paths["bismark_extract_report"])
        )
        canonical["non_cpg_metric_source"] = "bismark_extraction_chg_chh"
    updates.update(canonical)

    required_outputs = [
        (paths["cpg"], "CpG result"),
    ]
    if backend == "biscuit":
        required_outputs.append((paths["biscuit_qc_dir"], "BISCUIT QC directory"))
    else:
        required_outputs.extend(
            [
                (paths["bismark_align_report"], "Bismark alignment report"),
                (paths["bismark_extract_report"], "Bismark extraction report"),
                (paths["bismark_mbias"], "Bismark M-bias report"),
            ]
        )
    for required_path, label in required_outputs:
        if not required_path.exists():
            raise FileNotFoundError(f"Active cell is missing {label}: {required_path}")
    if generate_snp and not paths["snp"].is_file():
        raise FileNotFoundError(f"Active cell is missing enabled SNP result: {paths['snp']}")
    if publish_rna_bam:
        for rna_path, label in (
            (paths["rna_bam"], "SRD RNA BAM"),
            (paths["rna_bam_index"], "SRD RNA BAM index"),
        ):
            if not rna_path.is_file():
                raise FileNotFoundError(f"Active cell is missing {label}: {rna_path}")

    updates["methylation_backend"] = backend
    updates["backend_qc_status"] = "PASS"
    updates["cpg_path"] = str(public["cpg"])
    updates["cpg_representation"] = CPG_REPRESENTATION
    updates["snp_path"] = str(public["snp"]) if generate_snp else ""
    updates["rna_bam_path"] = str(public["rna_bam"]) if publish_rna_bam else ""
    updates["analysis_status"] = "PASS"
    return cell_id, updates


def _apply_cell_metrics(
    cells: dict[str, dict[str, object]],
    *,
    layout: ProjectLayout,
    published_layout: ProjectLayout,
    backend: str,
    publish_rna_bam: bool,
    generate_snp: bool,
    bismark: str | None,
    keep_final_bam: bool,
) -> None:
    published_results = published_layout.result_root
    project_root = published_results.parent
    bismark_version = _bismark_version_text(bismark) if backend == "bismark" else ""
    jobs: list[dict[str, object]] = []
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
        paths = layout.cell(cell_id)
        public_paths = published_layout.cell(cell_id)
        jobs.append(
            {
                "cell_id": cell_id,
                "record": dict(record),
                "backend": backend,
                "publish_rna_bam": publish_rna_bam,
                "generate_snp": generate_snp,
                "keep_final_bam": keep_final_bam,
                "paths": {
                    "rastair_qc": paths.rastair_qc,
                    "marked_bam": paths.marked_bam,
                    "marked_bam_index": paths.marked_bam_index,
                    "cutadapt_json": paths.cutadapt_json,
                    "dupsifter_multiqc": paths.dupsifter_multiqc,
                    "mapq_table": paths.biscuit_qc_source("mapq_table.txt"),
                    "cph_table": paths.biscuit_qc_source("CpHRetentionByReadPos.txt"),
                    "biscuit_qc_dir": paths.biscuit_qc_dir,
                    "high_cph_summary": paths.high_cph_summary,
                    "bismark_align_report": paths.bismark_align_report,
                    "bismark_dedup_report": paths.droplet_umi_qc if record["protocol"] == "droplet" else paths.bismark_dedup_report,
                    "bismark_extract_report": paths.bismark_extract_report,
                    "bismark_mbias": paths.bismark_mbias,
                    "cpg": paths.cpg,
                    "snp": paths.snp,
                    "rna_bam": paths.rna_bam,
                    "rna_bam_index": paths.rna_bam_index,
                },
                "public": {
                    "marked_bam": str(paths.marked_bam.relative_to(project_root)),
                    "marked_bam_index": str(paths.marked_bam_index.relative_to(project_root)),
                    "cpg": str(public_paths.cpg.relative_to(project_root)),
                    "snp": str(public_paths.snp.relative_to(project_root)),
                    "rna_bam": str(paths.rna_bam.relative_to(project_root)),
                },
                "bismark_version": bismark_version,
            }
        )
    if not jobs:
        return
    workers = int(os.environ.get("DNA_PIPELINE_FINAL_WORKERS", "0") or 0)
    if workers <= 0:
        workers = min(8, os.cpu_count() or 1)
    workers = max(1, min(workers, len(jobs)))
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
                pending = list(pool.map(_metrics_for_cell, jobs))
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
        try:
            order_key: tuple[int, object] = (0, int(raw_order))
        except ValueError:
            order_key = (1, raw_order)
        return (
            str(item.get("project_sample_id", "")),
            order_key,
            str(item.get("plate_id", "")),
            str(item.get("sample_id", "")),
        )

    sorted_rows = sorted(rows, key=cell_sort_key)
    if pair_qc_multiqc_json is not None:
        _write_pair_level_multiqc(sorted_rows, Path(pair_qc_multiqc_json))
    if cutadapt_qc_multiqc_json is not None:
        _write_cutadapt_multiqc(sorted_rows, Path(cutadapt_qc_multiqc_json))

    write_tsv_atomic(output, sorted_rows, fields)


MULTIQC_BISCUIT_SOURCE_SUFFIXES = (
    "strand_table.txt",
    "totalReadConversionRate.txt",
    "CpGRetentionByReadPos.txt",
    "mapq_table.txt",
    "totalBaseConversionRate.txt",
    "CpHRetentionByReadPos.txt",
    "isize_table.txt",
    "dup_report.txt",
)


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
    suffixes = "".join(path.suffixes)

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
    if tuple(parts[-4:-2]) == ("biscuit", "qc") and suffixes in {
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


def write_multiqc_file_list(
    inputs: Iterable[str | Path],
    out_list: str,
) -> tuple[int, list[str]]:
    """写出 MultiQC 消费的唯一显式源列表，源依赖由 Snakemake 追踪。"""
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
    file_list.parent.mkdir(parents=True, exist_ok=True)
    file_list_tmp = file_list.with_name(file_list.name + ".tmp")
    file_list_tmp.write_text(
        "".join(f"{path}\n" for path in sorted(unique_sources)),
        encoding="utf-8",
    )
    file_list_tmp.replace(file_list)

    logger.info(
        "Prepared explicit MultiQC file list with %d sources; skipped %d dup reports",
        len(unique_sources),
        len(skipped),
    )
    return len(unique_sources), skipped


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


MANIFEST_FIELDS = PROJECT_SAMPLE_MANIFEST_FIELDS


def _read_existing_manifest_loose(path: Path) -> list[ProjectSampleRow]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if not fields:
            return []
        if fields != list(MANIFEST_FIELDS):
            raise ValueError(
                f"Existing sample manifest {path} header/order does not match the current contract: "
                f"expected={list(MANIFEST_FIELDS)!r}, observed={fields!r}"
            )
        rows: list[ProjectSampleRow] = []
        for raw in reader:
            if None in raw:
                raise ValueError(f"Existing sample manifest {path} line {reader.line_num} has too many fields")
            row = {field: str(raw.get(field, "") or "").strip() for field in MANIFEST_FIELDS}
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
    for old in existing_rows:
        current_dna = str(
            old.get("dna_raw_sample", "") or ""
        ).strip()
        current_rna = str(old.get("rna_sample", "") or "").strip()
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
        external_rna = {row["rna_sample"] for row in rows_by_id.values() if row["rna_sample"]}
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
            for field in MANIFEST_FIELDS
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


    new_text = render_tsv(rows, MANIFEST_FIELDS)
    validated = _parse_project_sample_manifest(
        io.StringIO(new_text), path=manifest_path, species=species, protocol=protocol
    )
    map_discovered_samples_to_manifest(discovered, validated, protocol=protocol)
    changed = _canonical_project_manifest_rows(existing_rows) != rows
    backup_path: Path | None = None

    if not dry_run:
        if fastq_renames:
            _apply_fastq_renames(raw_dir, fastq_renames)
        try:
            if changed:
                if manifest_path.is_file() and existing_rows:
                    backup_path = manifest_path.with_name(f"{manifest_path.name}.bak")
                    shutil.copy2(manifest_path, backup_path)
                write_tsv_atomic(manifest_path, rows, MANIFEST_FIELDS)
        except Exception:
            for source, target in reversed(fastq_renames):
                target.rename(source)
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


BSCONV_P_COLUMNS = 9


def open_text(path: Path) -> TextIO:
    if str(path) == "-":
        return sys.stdin
    return path.open("r", encoding="utf-8")


def normalize_read_name(value: str) -> str:

    return str(value).strip().split()[0]


def indexed_excluded_contigs(bam_path: Path, samtools: str, excluded: set[str]) -> list[str]:
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
            contigs.append(fields[0])
    return contigs


def excluded_reads_from_bam(bam_path: Path, samtools: str, excluded: set[str]) -> set[str]:
    """经 samtools idxstats/view 流式收集排除 contig 上的 QNAME，不物化完整 SAM。"""
    contigs = indexed_excluded_contigs(bam_path, samtools, excluded)
    if not contigs:
        return set()

    proc = subprocess.Popen(
        [samtools, "view", str(bam_path), *contigs],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=1,
    )
    assert proc.stdout is not None
    names = set()
    for line in proc.stdout:
        if not line:
            continue
        names.add(normalize_read_name(line.split("\t", 1)[0]))
    proc.stdout.close()
    stderr = proc.stderr.read() if proc.stderr is not None else ""
    returncode = proc.wait()
    if returncode != 0:
        raise RuntimeError(stderr.strip() or f"samtools view failed for {bam_path}")
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
        not isinstance(value, str) or not value.strip() for value in raw
    ):
        raise ValueError("--excluded-contigs-json must be a JSON array of non-empty strings")
    if len(set(raw)) != len(raw):
        raise ValueError("--excluded-contigs-json must not repeat a contig")
    return {str(value).strip() for value in raw}


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


def run_bismark_pair_split(args) -> int:
    """把 BAM 路由到 bio/excluded，标准输出仅返回两路 pair 计数供当前规则使用。"""
    excluded_contigs = _parse_excluded_contigs_json(args.excluded_contigs_json)
    args.bio.parent.mkdir(parents=True, exist_ok=True)
    args.excluded.parent.mkdir(parents=True, exist_ok=True)
    print(*split_pairs(args.input_bam, args.bio, args.excluded, excluded_contigs))
    return 0


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


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Alopex 核心模块内部 CLI。")
    commands = parser.add_subparsers(dest="command", required=True)


    init = commands.add_parser("init-project", help="初始化项目骨架（目录/模板/launcher symlink）。")
    init.add_argument("path", nargs="?", default=".", help="项目目录；默认当前目录。")
    init.add_argument("--pipeline-root", required=True, help="Pipeline 安装根，用于模板与 symlink。")
    clean_demux = commands.add_parser("clean-demux-scratch", help="完成发布后按提交边界顺序回收 demux scratch。")
    clean_demux.add_argument("--project", required=True)
    scalars = commands.add_parser("config-scalars", help="把平台路由所需的 config 标量输出为可 eval 的 shell 赋值。")
    scalars.add_argument("--config", required=True)
    resolve_tmp = commands.add_parser("resolve-tmpdir", help="把 slurm.node_tmpdir 展开为绝对路径（根目录由语义校验拒绝）。")
    resolve_tmp.add_argument("--path", required=True)
    synthesize = commands.add_parser("synthesize-config", help="合成 test 模式 effective config。")
    synthesize.add_argument("--source", required=True)
    synthesize.add_argument("--destination", required=True)
    synthesize.add_argument("--route", required=True)
    refresh = commands.add_parser("refresh-manifest", help="根据 01_raw/ 刷新 sample manifest。")
    refresh.add_argument("--project", help="项目根目录；默认使用当前目录。")
    refresh.add_argument("--configfile", help="显式指定项目 config；仅接受 <project>/00_config/config.yaml。")
    refresh.add_argument("--dry-run", action="store_true", help="只预览候选 manifest，不写入文件。")


    create = commands.add_parser("snapshot-create", help="创建/复用内容寻址 run snapshot。")
    create.add_argument("--pipeline-root", required=True)
    create.add_argument("--project", required=True)
    create.add_argument("--config", required=True)
    create.add_argument("--control-env", required=True)
    create.add_argument("--biscuit-env", required=True)
    create.add_argument("--bismark-env", required=True)
    create.add_argument("--rastair-env", required=True)
    create.add_argument("--allow-external-config", action="store_true")
    status = commands.add_parser("project-status", help="输出项目交付状态 JSON。")
    status.add_argument("--project", required=True)
    status.add_argument("--current-snapshot")
    ready = commands.add_parser("delivery-ready", help="写出 delivery_ready 发布事务标记。")
    ready.add_argument("--ready", required=True)
    ready.add_argument("--staged-results", required=True)
    ready.add_argument("--run-snapshot", required=True)
    atomic = commands.add_parser("atomic-publish", help="scratch 文件组校验后原子发布。")
    atomic.add_argument("--file", action="append", nargs=2, metavar=("SOURCE", "DESTINATION"),
                        required=True, dest="files")
    atomic.add_argument("--commit-last", default="")
    atomic.add_argument("--bam-destination", default="")
    atomic.add_argument("--samtools", default="samtools")


    taps_qc = commands.add_parser("taps-alignment-qc", help="核算 TAPS BAM 保留与调用集合。")
    for field in ("bam", "samtools", "scratch", "sample", "output"):
        taps_qc.add_argument("--" + field, required=True)
    taps_cpg = commands.add_parser("rastair-cpg", help="转换原生 Rastair 计数并记录对照 QC。")
    for field in ("bed", "reference", "alignment-qc", "output-qc"):
        taps_cpg.add_argument("--" + field, required=True)

    write_initial = commands.add_parser("write-cell-manifest", help="写出初始 cell manifest。")
    write_initial.add_argument("project_manifest")
    write_initial.add_argument("barcode_map")
    write_initial.add_argument("out_tsv")
    write_initial.add_argument("--species", required=True)
    write_initial.add_argument("--protocol", required=True)
    write_initial.add_argument("--called-cells", action="append", type=lambda value: value.split(":", 1), default=[])
    umi = commands.add_parser("dedup-droplet-bam", help="按物理 R1 与 UMI 去重并恢复 Bismark 标志。")
    for field in ("bam", "output", "qc", "scratch", "umi-tools", "sample"):
        umi.add_argument("--" + field, required=True)
    umi.add_argument("--threads", type=int, default=1)
    calling = commands.add_parser("call-droplet-cells", help="按双峰谷底调用 Droplet 细胞。")
    for field in ("counts", "output", "metrics"):
        calling.add_argument("--" + field, required=True)
    calling.add_argument("--r1", action="append", required=True)
    calling.add_argument("--demux-binary", required=True)
    calling.add_argument("--threads", type=int, default=1)
    index_demux = commands.add_parser("index-demux-fastq", help="为原始 Droplet gzip 建立可直接并行读取的区段索引。")
    index_demux.add_argument("--r1", action="append", required=True)
    index_demux.add_argument("--r2", action="append", required=True)
    index_demux.add_argument("--output", required=True)
    index_demux.add_argument("--threads", type=int, required=True)
    chunk_demux = commands.add_parser("demux-chunk", help="在 job scratch 解复用一个块并发布打包结果。")
    for field in ("index-dir", "chunk", "binary", "barcode-map", "sample", "scratch", "pack", "metadata"):
        chunk_demux.add_argument("--" + field, required=True)
    chunk_demux.add_argument("--threads", type=int, required=True)
    merge_demux = commands.add_parser("merge-demux-chunks", help="按块序合并结果并执行全库细胞保留阈值。")
    for field in ("index-dir", "packs-dir", "dna-dir", "report"):
        merge_demux.add_argument("--" + field, required=True)
    merge_demux.add_argument("--threshold", type=int, required=True)
    merge_demux.add_argument("--threads", type=int, required=True)
    write_metadata = commands.add_parser("write-demux-metadata", help="demux 事务内写出 MQC 与 manifest。")
    for arg in ("report_path", "mqc_path", "manifest_path", "dna_dir", "rna_dir",
                "pipeline_mode", "route_id", "downstream_dna", "raw_sample", "cell_manifest_path"):
        write_metadata.add_argument(arg)
    write_metadata.add_argument("--calling-metrics", default="")
    write_metadata.add_argument("--structure-metrics", default="")
    write_completion = commands.add_parser("write-demux-completion", help="demux 事务内写出 completion 提交记录。")
    write_completion.add_argument("out_json")
    write_completion.add_argument("--run-id", required=True)
    write_completion.add_argument("--project-sample-id", required=True)
    write_completion.add_argument("--manifest", required=True)
    write_completion.add_argument("--report", required=True)
    write_completion.add_argument("--mqc", required=True)
    write_completion.add_argument("--dna-dir", required=True)
    write_completion.add_argument("--rna-dir", required=True)
    write_completion.add_argument("--cell-manifest", required=True)
    collect = commands.add_parser("collect-rna-fastq", help="收集 RNA FASTQ 到保留区。")
    collect.add_argument("--manifest", required=True)
    collect.add_argument("--intermediate-dir", required=True)
    collect.add_argument("--marker", required=True)
    collect.add_argument("--run-id", required=True)
    write_final = commands.add_parser("write-sample-manifest", help="写出最终 sample manifest。")
    for arg in ("initial_manifest", "intermediate_dir", "staged_results_dir",
                "published_results_dir", "out_tsv", "pair_qc_json", "cutadapt_qc_json"):
        write_final.add_argument(arg)
    write_final.add_argument("--bismark", required=True)
    write_final.add_argument("--options-json", required=True)
    write_dupsifter = commands.add_parser("write-dupsifter-qc", help="写出 Dupsifter MultiQC 内容。")
    write_dupsifter.add_argument("stat_path")
    write_dupsifter.add_argument("out_json")
    write_dupsifter.add_argument("--sample", required=True)
    write_dupsifter.add_argument("--metadata", required=True, help="Dupsifter 作业 metadata JSON。")
    extract = commands.add_parser("extract-bsconv", help="bsconv TSV 高 CpH 分类提取。")
    extract.add_argument("--bsconv-tsv", required=True, type=Path)
    extract.add_argument("--threshold", required=True, type=float)
    extract.add_argument("--high-read-names", required=True, type=Path)
    extract.add_argument("--summary-json", required=True, type=Path)
    extract.add_argument("--multiqc-json", type=Path)
    extract.add_argument("--sample-name", required=True)
    extract.add_argument("--excluded-contigs-json", required=True,
                         help="JSON array of canonical reference contig names excluded from assessment.")
    extract.add_argument("--alignment-bam", required=True, type=Path)
    extract.add_argument("--samtools", default="samtools")
    extract.add_argument("--protocol", choices=("cabernet", "srd", "droplet"), required=True)
    extract.add_argument("--high-cph-role", choices=("cdna_contamination", "residual_high_cph"),
                         required=True)
    split = commands.add_parser("split-bismark-pairs", help="Bismark BAM 配对路由与计数。")
    split.add_argument("--excluded-contigs-json", required=True)
    split.add_argument("--input-bam", type=Path, required=True)
    split.add_argument("--bio", type=Path, required=True)
    split.add_argument("--excluded", type=Path, required=True)
    demux_identity = commands.add_parser("demux-identity", help="计算源码指纹或校验 demux_rs 编译期身份。")
    demux_identity.add_argument("operation", choices=("revision", "verify-binary"))
    demux_identity.add_argument("path")
    demux_identity.add_argument("expected_revision", nargs="?")
    summary = commands.add_parser("write-filter-summary", help="写出 Bismark filter summary（schema v3）。")
    summary.add_argument("--sample", required=True)
    summary.add_argument("--candidate-pairs", type=int, required=True)
    summary.add_argument("--removed-pairs", type=int, required=True)
    summary.add_argument("--excluded-pairs", type=int, required=True)
    summary.add_argument("--processed-bam", type=Path, required=True)
    summary.add_argument("--samtools", required=True)
    summary.add_argument("--percentage", type=float, required=True)
    summary.add_argument("--minimum-count", type=int, required=True)
    summary.add_argument("--protocol", required=True)
    summary.add_argument("--high-cph-role", required=True)
    summary.add_argument("--summary-path", type=Path, required=True)
    summary.add_argument("--multiqc-path", type=Path, required=True)


    bench_validate = commands.add_parser("benchmark-validate", help="校验 Doctor fixture SHA256 与读数守恒。")
    bench_validate.add_argument("--root", required=True)
    bench_validate.add_argument("--manifest", required=True)
    bench_stage = commands.add_parser("benchmark-stage", help="合成沙箱路线 config。")
    bench_stage.add_argument("--project", required=True)
    bench_stage.add_argument("--reference", required=True)
    bench_stage.add_argument("--route", required=True)
    bench_verify = commands.add_parser("benchmark-verify", help="校验 sealed 输出守恒。")
    bench_verify.add_argument("--project", required=True)
    bench_verify.add_argument("--manifest", required=True)
    bench_verify.add_argument("--route", required=True)
    env_record = commands.add_parser("env-record", help="为一个新 prefix 记录环境身份。")
    env_record.add_argument("prefix")
    env_record.add_argument("--envs", required=True, help="五角色 envs.yaml 路径。")
    env_record.add_argument("--role", required=True)
    env_spec = commands.add_parser("env-spec", help="抽取角色 section 为独立 conda spec 文件。")
    env_spec.add_argument("--envs", required=True)
    env_spec.add_argument("--role", required=True)
    env_spec.add_argument("--destination", required=True)
    env_finalize = commands.add_parser("env-finalize", help="写出 release manifest。")
    env_finalize.add_argument("release")
    env_finalize.add_argument("release_id")
    env_finalize.add_argument("expected_subdir")
    env_publish = commands.add_parser("env-publish", help="原子切换 conda/current。")
    env_publish.add_argument("current")
    env_publish.add_argument("target")
    env_validate = commands.add_parser("env-validate", help="只读校验 current release。")
    env_validate.add_argument("current")
    env_validate.add_argument("expected_subdir")
    env_validate.add_argument("--envs", required=True)
    env_reusable = commands.add_parser("env-reusable", help="判断既有 prefix 可否复用。")
    env_reusable.add_argument("prefix")
    env_reusable.add_argument("--envs", required=True)
    env_reusable.add_argument("--role", required=True)


    fai_check = commands.add_parser("fai-check", help="FAI 契约校验。")
    fai_check.add_argument("--mode", choices=("source", "assembled"), required=True)
    fai_check.add_argument("--label", required=True)
    fai_check.add_argument("--fai", type=Path, required=True)
    fai_check.add_argument("--fasta", type=Path, required=True)
    suffixes = commands.add_parser("index-suffixes", help="列出 BISCUIT index 后缀清单。")
    suffixes.add_argument("--transient", action="store_true")
    indexes = commands.add_parser("index-paths", help="列出选中 backend/模式的完整 index 文件家族。")
    indexes.add_argument("--reference", required=True)
    indexes.add_argument("--backend", choices=("biscuit", "bismark", "rastair"), required=True)
    indexes.add_argument("--bismark-local-alignment", action="store_true")
    indexes.add_argument("--bismark-root")
    build_downstream = commands.add_parser("build-downstream", help="构建下游 TSS/single-CpG reference。")
    build_downstream.add_argument("--species", required=True, choices=sorted(MANAGED_SPECIES))
    build_downstream.add_argument("--fasta", required=True)
    build_downstream.add_argument("--annotation-gtf", required=True)
    build_downstream.add_argument("--tss-bed", required=True)
    build_downstream.add_argument("--single-cpg-bed", required=True)
    spikein = commands.add_parser("append-spikein", help="把 lambda/pUC19 追加到 assembled FASTA。")
    spikein.add_argument("--lambda-fasta", type=Path, required=True)
    spikein.add_argument("--puc19-fasta", type=Path, required=True)
    spikein.add_argument("--target-fasta", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    command = args.command
    try:

        if command == "init-project":
            project = init_project_directory(args.path, args.pipeline_root)
            print(f"Initialized Alopex project: {project}")
            print("Next:")
            print("  1. Put paired FASTQ files under 01_raw/")
            print("  2. ./run_pipeline.sh --refresh-manifest")
            print("  3. ./run_pipeline.sh --dry-run")
            print("  4. ./run_pipeline.sh")
        elif command == "config-scalars":
            config = load_project_config(args.config)
            for var, value in export_config_scalars(config).items():
                print(f"{var}={shlex.quote(value)}")
        elif command == "clean-demux-scratch":
            clean_demux_scratch_after_publish(args.project)
        elif command == "resolve-tmpdir":
            print(Path(args.path).expanduser().resolve(strict=False))
        elif command == "synthesize-config":
            synthesize_route_config(Path(args.source), Path(args.destination), args.route)
        elif command == "refresh-manifest":
            try:
                project, config_path = _resolve_config(args.project, args.configfile)
                result = refresh_manifest(
                    project_dir=project,
                    config_path=config_path,
                    dry_run=bool(args.dry_run),
                )
            except Exception as exc:
                print(f"ERROR: {exc}", file=sys.stderr)
                return 2
            return print_refresh_result(args, result)

        elif command == "snapshot-create":
            print(create_run_snapshot(
                pipeline_root=args.pipeline_root,
                project_dir=args.project,
                config_path=args.config,
                control_env=args.control_env,
                biscuit_env=args.biscuit_env,
                bismark_env=args.bismark_env,
                rastair_env=args.rastair_env,
                allow_external_config=args.allow_external_config,
            ))
        elif command == "project-status":
            result = project_status(
                project_dir=args.project,
                current_snapshot=args.current_snapshot,
            )
            print(json.dumps(result, ensure_ascii=False, sort_keys=True))
        elif command == "delivery-ready":
            run_snapshot = str(args.run_snapshot).strip()
            if not run_snapshot:
                raise ValueError(
                    "publishing requires a launcher-managed run snapshot. "
                    "Direct Snakemake runs stop at 02_work/results_stage; rerun via "
                    "core/run_pipeline.sh to publish."
                )
            write_delivery_ready(
                args.ready,
                staged_results_dir=args.staged_results,
                run_snapshot_path=run_snapshot,
            )
        elif command == "atomic-publish":
            pairs = [(Path(source), Path(destination)) for source, destination in args.files]
            publish(
                pairs,
                commit_last=Path(args.commit_last) if args.commit_last else None,
                bam_destination=Path(args.bam_destination) if args.bam_destination else None,
                samtools=args.samtools,
            )

        elif command == "dedup-droplet-bam":
            deduplicate_droplet_bam(args.bam, args.output, args.qc, sample=args.sample, scratch=args.scratch, umi_tools=args.umi_tools, threads=args.threads)
        elif command == "call-droplet-cells":
            call_droplet_cells(args.counts, args.output, args.metrics, r1=args.r1, demux_binary=args.demux_binary, threads=args.threads)
        elif command == "index-demux-fastq":
            index_demux_fastqs(args.r1, args.r2, args.output, threads=args.threads)
        elif command == "demux-chunk":
            run_demux_chunk(args.index_dir, args.chunk, binary=args.binary, barcode_map=args.barcode_map,
                            sample=args.sample, scratch=args.scratch, pack=args.pack, metadata=args.metadata, threads=args.threads)
        elif command == "merge-demux-chunks":
            merge_demux_chunks(args.index_dir, args.packs_dir, args.dna_dir, args.report,
                               threshold=args.threshold, threads=args.threads)
        elif command == "write-cell-manifest":
            write_initial_cell_manifest(
                args.project_manifest, args.barcode_map, args.out_tsv,
                species=args.species, protocol=args.protocol, called_cells=args.called_cells,
            )
        elif command == "write-demux-metadata":
            write_run_metadata(
                args.report_path, args.mqc_path, args.manifest_path, args.dna_dir,
                args.rna_dir, args.pipeline_mode, args.route_id, args.downstream_dna,
                args.raw_sample, args.cell_manifest_path,
                calling_metrics=args.calling_metrics, structure_metrics=args.structure_metrics,
            )
        elif command == "write-demux-completion":
            write_demux_completion(
                args.out_json,
                run_id=args.run_id,
                project_sample_id=args.project_sample_id,
                manifest_path=args.manifest,
                report_path=args.report,
                mqc_path=args.mqc,
                dna_dir=args.dna_dir,
                rna_dir=args.rna_dir,
                cell_manifest_path=args.cell_manifest,
            )
        elif command == "collect-rna-fastq":
            collect_rna_fastq(args.manifest, args.intermediate_dir, args.marker, args.run_id)
        elif command == "write-sample-manifest":
            write_final_sample_manifest(
                args.initial_manifest, args.intermediate_dir, args.staged_results_dir,
                args.out_tsv,
                options=json.loads(args.options_json),
                published_results_dir=args.published_results_dir,
                bismark=args.bismark,
                pair_qc_multiqc_json=args.pair_qc_json,
                cutadapt_qc_multiqc_json=args.cutadapt_qc_json,
            )
        elif command == "taps-alignment-qc":
            write_taps_alignment_qc(args.bam, args.samtools, args.scratch, args.sample, args.output)
        elif command == "rastair-cpg":
            emit_rastair_cpg(args.bed, args.reference, args.alignment_qc, args.output_qc, sys.stdout)
        elif command == "write-dupsifter-qc":
            write_dupsifter_mqc(
                args.stat_path, args.out_json, args.sample,
                metadata=json.loads(args.metadata),
            )
        elif command == "extract-bsconv":
            run_bsconv_extraction(args)
        elif command == "split-bismark-pairs":
            return run_bismark_pair_split(args)
        elif command == "demux-identity":
            if args.operation == "revision":
                if args.expected_revision is not None:
                    raise ValueError("demux-identity revision accepts only the source path")
                print(demux_source_revision(args.path))
            else:
                if args.expected_revision is None:
                    raise ValueError("demux-identity verify-binary requires the expected revision")
                verify_demux_binary_revision(args.path, args.expected_revision)
        elif command == "write-filter-summary":
            write_filter_summary(
                sample=args.sample,
                candidate_pairs=args.candidate_pairs,
                removed_pairs=args.removed_pairs,
                excluded_pairs=args.excluded_pairs,
                retained_pairs=_count_bam_pairs(args.processed_bam, args.samtools),
                percentage=args.percentage,
                minimum_count=args.minimum_count,
                protocol=args.protocol,
                high_cph_role=args.high_cph_role,
                summary=args.summary_path,
                multiqc=args.multiqc_path,
            )

        elif command == "benchmark-validate":
            validate_benchmark(Path(args.root), Path(args.manifest))
        elif command == "benchmark-stage":
            stage_route_config(Path(args.project), Path(args.reference), args.route)
        elif command == "benchmark-verify":
            verify_route_outputs(Path(args.project), Path(args.manifest), args.route)
        elif command == "env-record":
            record_environment_identity(
                prefix=args.prefix,
                role=args.role,
                spec_name=envs_spec_logical_name(args.role),
                spec_sha256=envs_section_sha256(args.envs, args.role),
            )
        elif command == "env-spec":
            write_envs_section_spec(args.envs, args.role, args.destination)
        elif command == "env-finalize":
            finalize_release(args.release, args.release_id, args.expected_subdir)
        elif command == "env-publish":
            publish_release(args.current, args.target)
        elif command == "env-validate":
            validate_release(args.current, args.expected_subdir, args.envs)
        elif command == "env-reusable":
            environment_is_reusable(args.prefix, args.role, args.envs)

        elif command == "fai-check":
            return 0 if fai_contract_check(args.mode, args.label, args.fai, args.fasta) else 1
        elif command == "index-suffixes":
            suffixes_list = BISCUIT_TRANSIENT_SUFFIXES if args.transient else BISCUIT_INDEX_SUFFIXES
            print("\n".join(suffixes_list))
        elif command == "index-paths":
            print("\n".join(map(str, backend_active_index_paths(
                args.reference, args.backend,
                bismark_local_alignment=args.bismark_local_alignment,
                bismark_root=args.bismark_root,
            ))))
        elif command == "append-spikein":
            append_spikein_contigs(args.lambda_fasta, args.puc19_fasta, args.target_fasta)
        elif command == "build-downstream":
            payload = build_downstream_references(
                species=args.species,
                fasta=args.fasta,
                annotation_gtf=args.annotation_gtf,
                tss_bed=args.tss_bed,
                single_cpg_bed=args.single_cpg_bed,
            )
            print(json.dumps(payload, ensure_ascii=False, sort_keys=True))
        else:
            raise ValueError(f"unsupported dna_pipeline command: {command}")
    except (
        FileNotFoundError,
        ValueError,
        RuntimeError,
        OSError,
        subprocess.SubprocessError,
        json.JSONDecodeError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
