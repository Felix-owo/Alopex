"""配置、协议、数据类型及只读参考契约；不调用交付或构建。"""
from __future__ import annotations

import copy
import json
import logging
import math
import os
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from types import MappingProxyType
from typing import Any, TypedDict

from .io import (
    canonical_json_bytes,
    sha256_bytes,
    sha256_file,
    stable_file_identity,
)
from .paths import CORE_ROOT


logger = logging.getLogger("dna_pipeline")


# 配置与契约

DROPLET_DESIGN_RESOURCE = "resources/droplet_ME5_U3CB_methylation.txt.gz"


DROPLET_DESIGN_SHA256 = "67981bee6acda8db275257a54e727d1e0a0bf4b884c10d8e761d04eaa9b53588"


SAMPLE_ID_RE = re.compile(r"^(?!\.{1,2}$)[A-Za-z0-9_.-]+$")


ENVIRONMENT_IDENTITY_SCHEMA = 8


IDENTITY_NAME = ".dna_pipeline_identity.json"


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


CONFIG_SCHEMA_YAML = r"""$schema: https://json-schema.org/draft/2020-12/schema
type: object
additionalProperties: false
required:
- schema_version
- species
- references
- analysis
- runtime
- slurm
- retention
- common
- cabernet
- srd
- taps
- droplet
properties:
  schema_version:
    description: 配置 schema 版本号，当前固定为 3。
    type: integer
    const: 3
  species:
    description: 当前项目的物种标识，决定默认 reference 选择与结果目录命名。
    type: string
    minLength: 1
    pattern: ^[A-Za-z0-9_.-]+$
  references:
    description: 物种 -> 标准 reference FASTA 路径映射；由 core/doctor.sh reference 构建与校验。
    type: object
    minProperties: 1
    propertyNames:
      pattern: ^[A-Za-z0-9_.-]+$
    additionalProperties:
      type: string
      minLength: 1
  analysis:
    description: 公共协议选择与 CpG 输出坐标。
    type: object
    additionalProperties: false
    required:
    - protocol
    - cpg_coordinate_base
    properties:
      protocol:
        description: 建库协议（srd / cabernet / taps / droplet）；manifest 的 rna_sample 在 srd 为本地 RNA tube，在 cabernet/taps/droplet
          为可选外部 RNA raw 关联。
        type: string
        enum:
        - srd
        - cabernet
        - taps
        - droplet
      cpg_coordinate_base:
        description: 四列 CpG、RawAdata 和下游 SingleCpG 位点的坐标起点；默认 1，设为 0 输出 0-based；不改变 BED、VCF 或 BAM 的原生坐标。
        type: integer
        enum:
        - 0
        - 1
  runtime:
    description: 与生物学结果无关的执行器运行时容错参数。
    type: object
    additionalProperties: false
    required:
    - local_executor_cores
    - latency_wait_seconds
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
    required:
    - jobs
    - controller_cores
    - account
    - partition
    - qos
    - node_tmpdir
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
        - type: string
          minLength: 1
        - type: 'null'
      partition:
        description: 提交作业的默认分区；null 沿用集群默认分区。
        oneOf:
        - type: string
          minLength: 1
        - type: 'null'
      qos:
        description: 提交作业的 QoS；字符串随作业提交，null 沿用集群默认 QoS。
        oneOf:
        - type: string
          minLength: 1
        - type: 'null'
      node_tmpdir:
        description: scratch 候选根；优先容量足够的本地磁盘，Slurm 中自动探测 /tmp 和 /var/tmp；跳过内存盘，无本地盘时警告并回退共享/未知盘。
        oneOf:
        - type: string
          pattern: ^/.*
          not:
            const: /
          allOf:
          - not:
              pattern: (^|/)\.\.?(/|$)
        - type: 'null'
  retention:
    description: 恢复边界策略；控制昂贵中间文件的保留面。
    type: object
    additionalProperties: false
    required:
    - keep_final_bam
    properties:
      keep_final_bam:
        description: true 保留 02_work/ 最终 BAM；false 在持锁成功发布后回收活跃 cell 的最终 BAM/索引，SRD RNA BAM 除外；默认 false。
        type: boolean
  cabernet:
    $ref: '#/$defs/plate_protocol'
  srd:
    $ref: '#/$defs/plate_protocol'
  taps:
    type: object
    additionalProperties: false
    properties:
      demux:
        $ref: '#/$defs/demux'
  droplet:
    type: object
    additionalProperties: false
    properties:
      demux:
        type: object
        additionalProperties: false
        properties:
          min_matched_read_pairs:
            $ref: '#/$defs/demux/properties/min_matched_read_pairs'
      high_cph:
        $ref: '#/$defs/high_cph'
      bismark:
        type: object
        additionalProperties: false
        properties:
          local_alignment:
            $ref: '#/$defs/bismark/properties/local_alignment'
          non_conversion_percentage_cutoff:
            $ref: '#/$defs/bismark/properties/non_conversion_percentage_cutoff'
          non_conversion_minimum_count:
            $ref: '#/$defs/bismark/properties/non_conversion_minimum_count'
  common:
    type: object
    additionalProperties: false
    required:
    - demux
    - high_cph
    - bismark
    - biscuit
    properties:
      demux:
        $ref: '#/$defs/demux'
        required:
        - min_matched_read_pairs
        - dna_w_spacer_len
      high_cph:
        $ref: '#/$defs/high_cph'
        required:
        - excluded_contigs
      bismark:
        $ref: '#/$defs/bismark'
        required:
        - library_type
        - local_alignment
        - non_conversion_percentage_cutoff
        - non_conversion_minimum_count
      biscuit:
        $ref: '#/$defs/biscuit'
        required:
        - library_type
        - high_cph_retention_threshold
        - generate_snp
$defs:
  plate_protocol:
    type: object
    additionalProperties: false
    required:
    - methylation_backend
    properties:
      methylation_backend:
        description: 本协议使用的甲基化后端。
        type: string
        enum:
        - biscuit
        - bismark
      demux:
        $ref: '#/$defs/demux'
      high_cph:
        $ref: '#/$defs/high_cph'
      bismark:
        $ref: '#/$defs/bismark'
      biscuit:
        $ref: '#/$defs/biscuit'
  demux:
    description: Rust 解复用器（demux_rs）参数；省略时使用内置默认。
    type: object
    additionalProperties: false
    properties:
      min_matched_read_pairs:
        description: 单 cell barcode 匹配 read pairs 的最低门槛，低于该值的 cell 不产出结果。
        type: integer
        minimum: 0
        maximum: 18446744073709551615
      dna_w_spacer_len:
        description: DNA 读段布局中 barcode 与 ME motif 间 W spacer 碱基数，须与实际建库一致。
        type: integer
        minimum: 0
        maximum: 20
  high_cph:
    description: Cabernet/SRD/Droplet 的 high-CpH / non-conversion 处理参数，固定统计并剔除；TAPS 不适用，无用户开关；省略时使用内置默认。
    type: object
    additionalProperties: false
    properties:
      excluded_contigs:
        description: 不参与 conversion 判断的 control/线粒体 contig 名称列表。
        type: array
        minItems: 0
        uniqueItems: true
        items:
          type: string
          pattern: ^[^,\s]+$
  biscuit:
    description: BISCUIT 后端专属参数；仅 methylation_backend=biscuit 时生效。
    type: object
    additionalProperties: false
    properties:
      library_type:
        description: BISCUIT 建库方向模式，决定允许的 strand model。
        type: string
        enum:
        - directional
        - non_directional
      high_cph_retention_threshold:
        description: read pair CpH retention 严格大于该值时判为 high-CpH；越低越严格。
        type: number
        minimum: 0
        maximum: 1
      generate_snp:
        description: 选定 BISCUIT 后端时发布 SNP；同一协议选择 Bismark 时必须为 false。
        type: boolean
  bismark:
    description: Bismark 后端专属参数；仅 methylation_backend=bismark 时生效。
    type: object
    additionalProperties: false
    properties:
      library_type:
        description: Bismark 建库方向模式，决定搜索的 bisulfite strand space。
        type: string
        enum:
        - directional
        - non_directional
        - pbat
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
        maximum: 18446744073709551615
"""


CONFIG_TEMPLATE_YAML = """schema_version: 3

species: hg38                    # 标准 hg38/mm10 reference 由 Doctor 构建

references:
  hg38: resources/hg38_reference/hg38.primary.lambda.puc19.fa  # primary genome + lambda + pUC19；相对 Pipeline 根目录
  mm10: resources/mm10_reference/mm10.primary.lambda.puc19.fa  # primary genome + lambda + pUC19；也接受绝对路径

analysis:
  protocol: cabernet              # cabernet / srd / taps / droplet；manifest 的 rna_sample 在 srd 为本地输入，cabernet/taps/droplet 为可选外部 RNA raw 关联
  cpg_coordinate_base: 1         # CpG pos、RawAdata 与 SingleCpG 位点默认 1-based；设为 0 输出 0-based，原生 BED/VCF/BAM 不变

runtime:
  local_executor_cores: 32        # local executor 的总 CPU 数（snakemake --cores）；影响吞吐不改变结果
  latency_wait_seconds: 120       # snakemake --latency-wait：等待输出文件在文件系统可见的秒数

retention:
  keep_final_bam: false          # true 保留最终 BAM/索引；false 在持锁成功发布后回收，SRD RNA BAM 除外

slurm:
  jobs: 30                        # SLURM 全局并发上限，含 Droplet 全库 calling、demux 块和后续 cell 作业
  controller_cores: 2             # Snakemake controller 本地进程的 CPU 数（--local-cores）
  account: null                  # 仅在集群要求时填写账户
  partition: null                 # 提交作业的默认分区；null 沿用集群默认分区
  qos: huge                       # 提交作业的默认 QoS；设为 null 时沿用集群默认 QoS
  node_tmpdir: /tmp               # 本地磁盘优先；Slurm 自动探测 /tmp、/var/tmp；无本地盘时警告回退，跳过内存盘

common:
  demux:
    min_matched_read_pairs: 10   # 单 cell barcode 匹配 read pairs 低于该值时不产出结果；0 关闭门槛
    dna_w_spacer_len: 0          # DNA 建库布局中 barcode 与 ME motif 间 W spacer 碱基数，须与实际建库一致
  high_cph:
    excluded_contigs: [pUC19, lambda, chrM]  # 这些对照不参与 high-CpH 筛除；名称须匹配参考 contig，TAPS 不适用
  bismark:
    library_type: non_directional   # directional / non_directional / pbat
    local_alignment: false          # false: combined-index end-to-end（默认）；true: faithful CT/GA + --local
    non_conversion_percentage_cutoff: 70    # 任一 mate 的 non-CG retention 达标且满足 minimum_count 时移除整对；越低越严格
    non_conversion_minimum_count: 5         # 至少这么多个 informative non-CG C 才启用百分比判断，防短 CpH 误杀
  biscuit:
    library_type: non_directional   # directional / non_directional
    high_cph_retention_threshold: 0.7   # read pair CpH retention 严格大于该值判为 high-CpH；越低越严格
    generate_snp: false             # BISCUIT-only SNP 输出开关；true 时发布 snps.bed.gz；TAPS 默认发布 VCF，不使用此开关

cabernet:
  methylation_backend: bismark  # biscuit / bismark；需要差异参数时在本模块覆盖 common 的同名叶子

srd:
  methylation_backend: bismark  # biscuit / bismark；需要差异参数时在本模块覆盖 common 的同名叶子

taps: {}

droplet: {}
"""


SUPPORTED_BISMARK_LIBRARY_TYPES = ("directional", "non_directional", "pbat")


SUPPORTED_BISCUIT_LIBRARY_TYPES = ("directional", "non_directional")


CONFIG_SECTION_DEFAULTS: dict[str, dict[str, object]] = {
    "analysis": {
        "cpg_coordinate_base": 1,
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
    "common": {
        "demux": {
            "min_matched_read_pairs": 10,
            "dna_w_spacer_len": 0,
        },
        "high_cph": {
            "excluded_contigs": ["pUC19", "lambda", "chrM"],
        },
        "bismark": {
            "library_type": "non_directional",
            "local_alignment": False,
            "non_conversion_percentage_cutoff": 70,
            "non_conversion_minimum_count": 5,
        },
        "biscuit": {
            "library_type": "non_directional",
            "high_cph_retention_threshold": 0.7,
            "generate_snp": False,
        },
    },
    "cabernet": {
        "methylation_backend": "bismark",
    },
    "srd": {
        "methylation_backend": "bismark",
    },
    "taps": {},
    "droplet": {},
}


_PROTOCOL_FIXED_CONFIG = {
    "taps": {"methylation_backend": "rastair"},
    "droplet": {
        "methylation_backend": "bismark",
        "demux": {"dna_w_spacer_len": 0},
        "bismark": {"library_type": "non_directional"},
    },
}


def apply_config_defaults(config: dict[str, object]) -> dict[str, object]:
    """递归补齐公共与协议模块的省略项；保留显式值并隔离各模块的可变默认。"""
    def merge(current, defaults, path):
        if not isinstance(current, dict):
            raise ValueError(f"config section {path!r} must be a mapping")
        for key, value in defaults.items():
            if key not in current:
                current[key] = copy.deepcopy(value)
            elif isinstance(value, dict):
                # YAML anchor 可跨 section 共享映射；补默认值只属于当前路径，
                # 不能向另一协议注入未显式填写（甚至禁止填写）的叶子。
                current[key] = copy.deepcopy(current[key])
                merge(current[key], value, f"{path}.{key}" if path else key)
    merge(config, CONFIG_SECTION_DEFAULTS, "")
    for protocol in ("cabernet", "srd", "taps", "droplet"):
        for key, value in config[protocol].items():
            if key in ("demux", "high_cph", "bismark", "biscuit") and not isinstance(value, dict):
                raise ValueError(f"config section {protocol}.{key!s} must be a mapping")
    return config


def get_protocol_config(config: Mapping[str, Any]) -> dict[str, Any]:
    """合并 common 与协议差异并注入固定协议值；结果与配置独立深拷贝。"""
    analysis = config.get("analysis")
    if not isinstance(analysis, Mapping):
        raise ValueError("analysis must be a mapping")
    protocol = analysis.get("protocol")
    if protocol not in ("cabernet", "srd", "taps", "droplet"):
        raise ValueError(f"Unsupported protocol: {protocol!r}")
    module = config.get(protocol)
    if not isinstance(module, dict) or not isinstance(config.get("common"), dict):
        raise ValueError("common and selected protocol must be mappings")
    groups = ("demux",) if protocol == "taps" else ("demux", "high_cph", "bismark")
    if protocol in ("cabernet", "srd"):
        groups += ("biscuit",)
    effective = {key: copy.deepcopy(value) for key, value in config["common"].items() if key in groups}
    for values in (module, _PROTOCOL_FIXED_CONFIG.get(protocol, {})):
        for key, value in values.items():
            if isinstance(value, dict):
                effective.setdefault(key, {}).update(copy.deepcopy(value))
            else:
                effective[key] = copy.deepcopy(value)
    return effective


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
    return _parse_project_config(Path(path).read_text(encoding="utf-8"), path=path)


def _parse_project_config(text: str, *, path: Path | str) -> dict[str, object]:
    """从一次已观测的配置内容解析，供普通加载与 snapshot 共用。"""
    import yaml

    class UniqueKeyLoader(yaml.SafeLoader):
        def __init__(self, stream):
            super().__init__(stream)
            self.checked_mappings = set()

        def flatten_mapping(self, node):
            # merge 来源可能只被 flatten 而不经过 construct_mapping，必须在展开前验证。
            # 锚点会复用节点；展开后出现的继承/显式覆盖不属于原映射的重复键。
            if node in self.checked_mappings:
                return super().flatten_mapping(node)
            self.checked_mappings.add(node)
            seen = set()
            for key_node, _ in node.value:
                if key_node.tag == "tag:yaml.org,2002:merge":
                    key = "<<"
                else:
                    key = self.construct_object(key_node)
                try:
                    duplicate = key in seen
                    seen.add(key)
                except TypeError as exc:
                    raise ValueError(f"Invalid config mapping key: {path}:{key_node.start_mark.line + 1}") from exc
                if duplicate:
                    raise ValueError(f"Duplicate config key {key!r}: {path}:{key_node.start_mark.line + 1}")
            return super().flatten_mapping(node)

    config = yaml.load(text, Loader=UniqueKeyLoader) or {}
    if not isinstance(config, dict):
        raise ValueError(f"project config root must be a mapping: {path}")
    return finalize_project_config(config)


def get_backend(config: dict[str, Any]) -> str:
    """返回选定的 backend，并拒绝该 backend 不支持的功能组合。"""
    analysis = config.get("analysis")
    if not isinstance(analysis, dict):
        raise ValueError("analysis must be a mapping")
    module = get_protocol_config(config)
    name = str(module.get("methylation_backend", "")).strip().lower()
    if name not in ("biscuit", "bismark", "rastair"):
        raise ValueError(
            f"Unsupported methylation backend {name!r}; choose biscuit, bismark or rastair"
        )
    if (analysis.get("protocol") == "taps") != (name == "rastair"):
        raise ValueError("protocol taps requires rastair; Cabernet/SRD require biscuit or bismark")
    biscuit = module.get("biscuit") or {}
    if bool(biscuit.get("generate_snp")) and name != "biscuit":
        raise ValueError(
            "biscuit.generate_snp=true is currently supported only by the Biscuit backend"
        )
    return name


@dataclass(frozen=True)
class BismarkAlignment:
    """描述所选比对模式的命令参数与并发实例；规则只投影实际消费的叶子值。"""

    local_alignment: bool
    library_flag: str
    instances: int

    @property
    def local_flag(self) -> str:
        return "--local" if self.local_alignment else ""

    @property
    def index_flag(self) -> str:
        return "" if self.local_alignment else "--combined_index"

    @property
    def score_min(self) -> str:
        return "G,20,8" if self.local_alignment else "L,0,-0.2"


def resolve_bismark_alignment(config: Mapping[str, object]) -> BismarkAlignment:
    """统一推导 library flag、模式与 Bowtie2 实例数，不查询 workflow 或文件系统。"""
    bismark = get_protocol_config(config).get("bismark")
    if not isinstance(bismark, Mapping):
        raise ValueError("bismark must be a mapping")
    value = str(bismark.get("library_type", "")).strip().lower()
    if value not in SUPPORTED_BISMARK_LIBRARY_TYPES:
        raise ValueError(
            "bismark.library_type must be directional, non_directional, or pbat"
        )
    local = bool(bismark.get("local_alignment", False))
    return BismarkAlignment(
        local_alignment=local,
        library_flag={"directional": "", "non_directional": "--non_directional", "pbat": "--pbat"}[value],
        instances=(4 if value == "non_directional" else 2) if local else 1,
    )


def cutadapt_adapter_specs(protocol: str) -> tuple[str, str]:
    """按协议返回含最短 overlap 的双端接头，不查询运行状态。"""
    if get_protocol(protocol).name == "droplet":
        return "N{9}CTGTCTCTTATACACATCT;min_overlap=12", "N{9}CTRTCTCTTATACACATCT;min_overlap=12"
    return "CTGTCTCTTATACACATCT;min_overlap=3", "N{9}CTGTCTCTTATACACATCT;min_overlap=12"


def biscuit_library_mode(config: Mapping[str, object]) -> int:
    """把命名的 paired-end BISCUIT library type 翻译为 ``align -b`` 的取值。"""
    biscuit = get_protocol_config(config).get("biscuit")
    if not isinstance(biscuit, Mapping):
        raise ValueError("biscuit must be a mapping")
    value = str(biscuit.get("library_type", "")).strip().lower()
    if value not in SUPPORTED_BISCUIT_LIBRARY_TYPES:
        raise ValueError("biscuit.library_type must be directional or non_directional")
    return {"non_directional": 0, "directional": 1}[value]


def validate_config_semantics(config: Mapping[str, object]) -> None:
    """校验 JSON Schema 无法表达的跨字段约束。"""
    cpg_representation(config["analysis"]["cpg_coordinate_base"])
    # JSON Schema 的 $ 可在末尾换行前匹配；科学身份不能随后被其他 consumer strip 成另一值。
    identifiers = [("species", config["species"])]
    identifiers.extend(("references key", name) for name in config["references"])
    for field, name in identifiers:
        if re.fullmatch(r"[A-Za-z0-9_.-]+", name) is None:
            raise ValueError(f"{field} must contain only ASCII letters, digits, dot, underscore or hyphen: {name!r}")
    for species, path in config["references"].items():
        if "\x00" in path:
            raise ValueError(f"references.{species} must not contain NUL")
    for section in ("common", "cabernet", "srd", "droplet"):
        for name in config[section].get("high_cph", {}).get("excluded_contigs", ()):
            if re.fullmatch(r"[^,\s]+", name) is None or "\x00" in name:
                raise ValueError(f"{section}.high_cph.excluded_contigs must contain literal names without whitespace, comma or NUL: {name!r}")
    backend = get_backend(dict(config))
    slurm = config["slurm"]
    # Bash 变量不能保存 NUL；阻止导出时被截去字符而提交到不同 account/QoS。
    for field in ("account", "partition", "qos"):
        value = slurm[field]
        if value is not None and (not value.strip() or any(ord(char) < 32 or ord(char) == 127 for char in value)):
            raise ValueError(f"slurm.{field} must be a non-blank string without control characters or null")

    node_tmpdir = slurm.get("node_tmpdir")
    if node_tmpdir is not None:
        if not isinstance(node_tmpdir, str) or not str(node_tmpdir).strip():
            raise ValueError("slurm.node_tmpdir must be a non-empty path or null")
        if Path(str(node_tmpdir)).expanduser().resolve(strict=False) == Path(
            Path(str(node_tmpdir)).expanduser().resolve(strict=False).anchor
        ):
            raise ValueError("slurm.node_tmpdir must not resolve to the filesystem root")
    for section in ("common", "cabernet", "srd"):
        biscuit = config[section].get("biscuit", {})
        if "high_cph_retention_threshold" not in biscuit:
            continue
        value = biscuit["high_cph_retention_threshold"]
        if (isinstance(value, bool) or not isinstance(value, (int, float))
                or not math.isfinite(float(value)) or not 0 <= value <= 1):
            raise ValueError(f"{section}.biscuit.high_cph_retention_threshold must be a finite number within [0, 1]")
    if backend == "biscuit":
        biscuit_library_mode(config)
    elif backend == "bismark":
        resolve_bismark_alignment(config)


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


_EXPORTED_SCALARS = {
    "CFG_LOCAL_CORES": "runtime.local_executor_cores",
    "CFG_LATENCY_WAIT": "runtime.latency_wait_seconds",
    "CFG_PROTOCOL": "analysis.protocol",
    "CFG_SPECIES": "species",
    "CFG_SLURM_JOBS": "slurm.jobs",
    "CFG_SLURM_CONTROLLER_CORES": "slurm.controller_cores",
    "CFG_SLURM_ACCOUNT": "slurm.account",
    "CFG_SLURM_PARTITION": "slurm.partition",
    "CFG_SLURM_QOS": "slurm.qos",
    "CFG_SLURM_NODE_TMPDIR": "slurm.node_tmpdir",
}


def export_config_scalars(config: Mapping[str, object]) -> dict[str, str]:
    """提取平台路由所需的 config 叶子标量（None 渲染为空串，bool 渲染为 true/false）。"""
    rendered: dict[str, str] = {"CFG_METHYLATION_BACKEND": get_backend(dict(config))}
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
        elif isinstance(value, float):
            # 当前导出的数字叶子均为整数；保持 schema 对 2.0 等整数值的语义，
            # 只规范 CLI 字面量，不把 2.5 截断成2，也不改配置或 sealed 输入副本。
            if not math.isfinite(value) or not value.is_integer():
                raise ValueError(f"config field must be an integer: {dotted}")
            rendered[var] = str(int(value))
        elif isinstance(value, (str, int)):
            rendered[var] = str(value)
        else:
            raise ValueError(f"config field must be a scalar: {dotted}")
    # default-resources 会求值表达式并展开环境变量；站点名称必须保持字面字符串。
    for field in ("ACCOUNT", "PARTITION"):
        value = rendered[f"CFG_SLURM_{field}"]
        rendered[f"CFG_SLURM_{field}_LITERAL"] = json.dumps(value, ensure_ascii=True).replace("$", r"\u0024")
    return rendered


CPG_REPRESENTATIONS = {
    0: "snapatac2_import_values_0based_v1",
    1: "cpg_four_column_1based_v1",
}


def cpg_representation(coordinate_base: int) -> str:
    """返回四列 CpG 的坐标身份；只接受整数 0 或 1。"""
    if type(coordinate_base) is not int or coordinate_base not in CPG_REPRESENTATIONS:
        raise ValueError("cpg_coordinate_base must be integer 0 or 1")
    return CPG_REPRESENTATIONS[coordinate_base]


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
        """公开四列 CpG 输入；pos 按 analysis.cpg_coordinate_base 写出。"""
        return self.result_root / "CpG" / f"{self.sample_id}.cpg.tsv.zst"

    @property
    def snp(self) -> Path:
        suffix = "vcf.gz" if self.backend == "rastair" else "bed.gz"
        return self.result_root / "SNP" / f"{self.sample_id}.snps.{suffix}"

    @property
    def snp_index(self) -> Path:
        return Path(f"{self.snp}.tbi")

    @property
    def public_result_files(self) -> dict[str, tuple[Path, ...]]:
        """按 final manifest 路径字段列出主文件与必需伴随文件；只定义路径，不判断是否启用。"""
        return {"cpg_path": (self.cpg,), "snp_path": (self.snp, self.snp_index)}

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
        return self.bismark_dir / f"{self.sample_id}.umi_consensus_mqc.json"

    @property
    def droplet_umi_cov(self) -> Path:
        return self.bismark_dir / f"{self.sample_id}.umi.cov.gz"

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


@dataclass(frozen=True)
class CellArtifacts:
    """按消费用途列出所选路线的产物；路径来自 CellOutputPaths，展示文件允许从持久摘要重建。"""

    metrics: Mapping[str, Path]
    deliverables: tuple[Path, ...]
    final_qc: tuple[Path, ...]
    multiqc: tuple[Path, ...]
    retained_bams: tuple[Path, ...]
    rna_bams: tuple[Path, ...]
    high_cph_display: Path | None

    @property
    def final_inputs(self) -> tuple[Path, ...]:
        """列出汇总实际消费的科学结果、持久 QC 和保留 BAM，不包含临时展示文件。"""
        return tuple(dict.fromkeys((*self.deliverables, *self.final_qc, *self.rna_bams)))

    def multiqc_sources(self) -> tuple[Path, ...]:
        """展示 JSON 已回收时只依赖保留的 summary，不回溯 BAM 规则。"""
        display = self.high_cph_display
        if display is not None and not display.is_file():
            return tuple(self.metrics["high_cph_summary"] if source == display else source
                         for source in self.multiqc)
        return self.multiqc


def cell_artifacts(
    paths: CellOutputPaths, protocol: str, *, generate_snp: bool, keep_final_bam: bool,
) -> CellArtifacts:
    """统一提供 DAG、final manifest 与 MultiQC 的逐细胞文件清单，不改变阶段完成边界。"""
    policy = get_protocol(protocol)
    metrics = {"cutadapt_json": paths.cutadapt_json}
    public_files = paths.public_result_files
    deliverables = list(public_files["cpg_path"])
    final_qc = [paths.cutadapt_json]
    multiqc = [paths.cutadapt_json]
    retained = paths.retained_bam_files if keep_final_bam else ()
    final_qc.extend(retained)
    rna = (paths.rna_bam, paths.rna_bam_index) if policy.publishes_rna_bam else ()
    if generate_snp:
        deliverables.extend(public_files["snp_path"])
    if paths.backend == "rastair":
        metrics["rastair_qc"] = paths.rastair_qc
        final_qc.append(paths.rastair_qc)
    else:
        metrics["high_cph_summary"] = paths.high_cph_summary
        final_qc.append(paths.high_cph_summary)
        high_cph_source = paths.high_cph_multiqc
        if paths.backend == "biscuit":
            metrics.update(dupsifter_multiqc=paths.dupsifter_multiqc,
                           mapq_table=paths.biscuit_qc_source("mapq_table.txt"),
                           cph_table=paths.biscuit_qc_source("CpHRetentionByReadPos.txt"))
            final_qc.extend((metrics["dupsifter_multiqc"], metrics["cph_table"]))
            multiqc.extend((metrics["dupsifter_multiqc"], high_cph_source,
                            *paths.biscuit_qc_sources(MULTIQC_BISCUIT_SOURCE_SUFFIXES)))
        else:
            metrics.update(bismark_align_report=paths.bismark_align_report,
                           bismark_extract_report=paths.droplet_umi_qc if policy.name == "droplet" else paths.bismark_extract_report,
                           bismark_dedup_report=paths.droplet_umi_qc if policy.name == "droplet" else paths.bismark_dedup_report)
            reports = [metrics["bismark_align_report"], metrics["bismark_extract_report"], paths.bismark_mbias]
            final_qc.extend(reports)
            multiqc.extend((*reports, high_cph_source))
            if policy.name != "droplet":
                final_qc.append(metrics["bismark_dedup_report"])
                multiqc.append(metrics["bismark_dedup_report"])
    return CellArtifacts(MappingProxyType(metrics), tuple(deliverables), tuple(final_qc),
                         tuple(multiqc), retained, rna,
                         paths.high_cph_multiqc if paths.backend != "rastair" else None)


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


@dataclass(frozen=True)
class DemuxIndex:
    """同一次 completion 校验产生的 run/cell 索引与有序通过名单；缓存由 checkpoint 调用方管理。"""

    state_dir: Path
    run_ids: tuple[str, ...]
    by_cell: Mapping[str, tuple[Path, DemuxRow]]
    passing_cells: tuple[str, ...]


class FastqPairFiles(TypedDict):
    r1: list[str]
    r2: list[str]


class MappedDemuxRun(FastqPairFiles):
    raw_sample: str
    sample_id: str
    route_id: str


DiscoveredSamples = dict[str, FastqPairFiles]


MappedDemuxRuns = dict[str, MappedDemuxRun]


def file_identity(path: str | Path) -> FileIdentity:
    source = Path(path).resolve()
    if not source.is_file():
        raise FileNotFoundError(f"Committed demux file is missing: {source}")
    return {
        "path": str(source),
        "size_bytes": source.stat().st_size,
        "sha256": sha256_file(source),
    }


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


# 配置与契约：参考输入

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


def managed_alignment_contigs(species: str) -> list[tuple[str, int]]:
    """固定assembly的主染色体＋真实辅助/decoy清单；不影响下游主染色体定义。"""
    path = CORE_ROOT / "reference" / f"{species}.alignment_contigs.json"
    records = [(name, int(length)) for name, length in json.loads(path.read_text())]
    return [*records, *SPIKE_CONTIGS.items()]


def primary_analysis_contigs(species: str, *, include_mitochondrial: bool = False) -> tuple[str, ...]:
    return tuple(name for name, _length in _PRIMARY_CONTIGS[species]
                 if include_mitochondrial or name != "chrM")


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

    # YAML 引号内的空格属于路径；裁剪可能静默选中另一份参考基因组。
    ref_path = str(refs[species] or "")
    if not ref_path.strip():
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
            "For standard hg38/mm10 references, place hg38.analysis.fa, mm10.fa, "
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


def _reference_file_state(path: Path) -> tuple[int, ...]:
    st = path.stat()
    return st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns


@lru_cache(maxsize=512)
def _reference_content_hash(path: str, state: tuple[int, ...]) -> str:
    """每个实际文件状态只哈希一次；读取期间变化不得生成有效绑定。"""
    p = Path(path)
    if _reference_file_state(p) != state:
        raise RuntimeError(f"reference file changed before hashing: {p}")
    value = sha256_file(p)
    if _reference_file_state(p) != state:
        raise RuntimeError(f"reference file changed while hashing: {p}")
    return value


def _reference_file_record(path: Path, recorded: Mapping | None = None) -> dict:
    state = _reference_file_state(path)
    if not path.is_file() or state[2] <= 0:
        raise ValueError(f"missing/empty reference artifact: {path}")
    value = (recorded["sha256"] if recorded and list(state) == recorded.get("state")
             else _reference_content_hash(str(path.resolve()), state))
    if recorded and (value != recorded.get("sha256") or state[2] != recorded.get("size_bytes")):
        raise ValueError(f"reference generation content mismatch: {path}")
    if re.fullmatch(r"[0-9a-f]{64}", str(value)) is None:
        raise ValueError(f"invalid reference digest: {path}")
    return {"sha256": value, "size_bytes": state[2], "state": list(state)}


def reference_binding_path(reference: Path, kind: str, *, root: Path | None = None) -> Path:
    """每个独立提交组件各有一个commit-last绑定；不共写全局可变清单。"""
    if kind == "fasta":
        return Path(str(reference) + ".reference.json")
    if kind == "rastair":
        return (root or Path(str(reference) + ".bwa")) / "reference.json"
    if kind == "biscuit":
        return Path(str(reference) + ".biscuit.reference.json")
    if kind in {"bismark-local", "bismark-combined"}:
        return (root or Path(str(reference) + ".bismark")) / (kind + ".reference.json")
    raise ValueError(f"unsupported reference component: {kind}")


def _read_reference_binding(path: Path, kind: str) -> dict:
    payload = json.loads(path.read_text())
    if payload.get("schema_version") != 1 or payload.get("kind") != kind:
        raise ValueError(f"invalid reference component binding: {path}")
    return payload


def _check_bound_artifacts(receipt: Path, payload: Mapping) -> None:
    files = payload.get("files")
    if not isinstance(files, dict) or not files:
        raise ValueError(f"empty reference artifact binding: {receipt}")
    for relative, record in files.items():
        _reference_file_record(receipt.parent / relative, record)


def reference_sequence_identity(reference: str | Path) -> dict:
    """内容身份与文件布局分离：同内容touch不强制科学重算，内容变化拒用旧派生。"""
    ref = Path(reference).resolve()
    receipt = reference_binding_path(ref, "fasta")
    recorded = {}
    if receipt.is_file():
        payload = _read_reference_binding(receipt, "fasta")
        recorded = payload["files"]
        if set(recorded) != {ref.name, ref.name + ".fai"}:
            raise ValueError(f"incomplete FASTA binding: {receipt}")
    records = {}
    for role, path in (("fasta", ref), ("fai", Path(str(ref) + ".fai"))):
        records[role] = _reference_file_record(path, recorded.get(os.path.relpath(path, receipt.parent)))
    return {role: {k: v for k, v in rec.items() if k != "state"} for role, rec in records.items()}


def reference_component_files(reference: Path, kind: str, root: Path | None = None) -> list[Path]:
    if kind == "fasta":
        return [reference, Path(str(reference) + ".fai")]
    if kind == "rastair":
        directory = root or Path(str(reference) + ".bwa")
        return [directory / ("genome" + suffix) for suffix in BWA_INDEX_SUFFIXES]
    if kind == "biscuit":
        return [Path(str(reference) + suffix) for suffix in BISCUIT_INDEX_SUFFIXES]
    if kind in {"bismark-local", "bismark-combined"}:
        directory = root or Path(str(reference) + ".bismark")
        return [directory / "genome.fa", *backend_active_index_paths(
            reference, "bismark", bismark_root=directory,
            bismark_local_alignment=kind == "bismark-local")]
    raise ValueError(f"unsupported reference component: {kind}")


def validate_reference_component(reference: Path, kind: str, *, root: Path | None = None) -> dict:
    receipt = reference_binding_path(reference, kind, root=root)
    payload = _read_reference_binding(receipt, kind)
    expected = {os.path.relpath(p, receipt.parent) for p in reference_component_files(reference, kind, root)}
    if set(payload["files"]) != expected:
        raise ValueError(f"incomplete reference file family: {receipt}")
    if kind != "fasta" and payload["reference"] != reference_sequence_identity(reference):
        raise ValueError(f"reference/component generation mismatch: {receipt}")
    _check_bound_artifacts(receipt, payload)
    return payload


def validate_downstream_reference(reference: Path, artifact: Path, kind: str,
                                  annotation: Path | None = None) -> dict:
    receipt = Path(str(artifact) + ".reference.json")
    payload = _read_reference_binding(receipt, kind)
    if set(payload["files"]) != {artifact.name} or payload["reference"] != reference_sequence_identity(reference):
        raise ValueError(f"downstream/reference generation mismatch: {artifact}")
    _check_bound_artifacts(receipt, payload)
    if annotation is not None and annotation.is_file():
        if _reference_file_record(annotation, payload["annotation"])["sha256"] != payload["annotation"]["sha256"]:
            raise ValueError(f"annotation changed; rebuild downstream references: {artifact}")
    return payload


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
    if expected_species in MANAGED_SPECIES and any(
        dict(observed).get(name) != length for name, length in MANAGED_CONTIGS[expected_species]
    ):
        raise ValueError(
            f"{expected_species} reference contig/order/length contract mismatch: {resolved_ref}"
        )
    content = reference_sequence_identity(resolved_ref)
    ref_stat = resolved_ref.stat()
    return {
        "content": content,
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


# Doctor 维护：环境和参考构建

ENVS_SPEC_NAME = "envs.yaml"


def envs_section(envs_path: str | Path, role: str) -> dict[str, object]:
    """用 Conda 自带的 ruamel.yaml 读取角色 spec，支持未安装 PyYAML 的 Seed。"""
    from ruamel.yaml import YAML

    payload = YAML(typ="safe", pure=True).load(Path(envs_path).read_text(encoding="utf-8"))
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
