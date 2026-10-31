# Alopex v16 分析参数

> 正式源码以 Git tag `v16` 为准。
>
> 本文件维护用户可调参数、工具原生默认与 Pipeline 取值对照、内部资源/I/O取值、QC 计数口径及下游重算开关。`config.yaml` 的可调范围见下一节；完整参数表以 † 标记可调项，其余为固定值。Notebook/Processor 参数单列，不写入主 config。科学不变量、恢复边界与维护验收见 [实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/dna_pipeline.py)。
> 固定参数不得因为“提高 mapping rate”或“增加 CpG 数量”而临时修改；若需改变，应先做同一测试数据的 mapping / MAPQ / M-bias / CpG recovery 回归。

## 当前用户可调参数

YAML 同一映射内的键必须唯一，包括 `<<` 合并的来源映射；重复 section 或叶子键会报文件与行号，避免静默覆盖科学参数。合法锚点继承、合并序列的优先级及显式键覆盖继承值仍遵循 YAML 语义。

`00_config/config.yaml` 使用 schema 3，保留一份文件，由 `analysis.protocol` 选择同名协议模块。
全局设置包含 `species`、`references`、`analysis`、`runtime`、`slurm` 和 `retention`；
`common` 统一声明 demux、high_cph、bismark、biscuit，共用参数只写一次。
默认值唯一出自 `core/dna_pipeline.py::CONFIG_SECTION_DEFAULTS`；内置 common 默认 → 显式 common → 选定协议同名叶子覆盖，最后应用协议固定策略，列表整体替换。
未知字段/null section拒绝，合法null叶子保留；加载、深拷贝及未选模块的校验/消费边界见[实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/dna_pipeline.py)。
下表列出各协议可以读取和覆盖的公共分组；模板中的 Cabernet/SRD 只写可选 backend，TAPS/Droplet 无差异时为 `{}`。

| 模块 | backend | 可读取/覆盖的公共分组 |
| --- | --- | --- |
| `cabernet` | `methylation_backend: bismark`（默认）或 `biscuit` | demux、high_cph、bismark、biscuit |
| `srd` | `methylation_backend: bismark`（默认）或 `biscuit` | demux、high_cph、bismark、biscuit |
| `taps` | 自动 Rastair，无 backend 参数 | demux |
| `droplet` | 自动 Bismark，无 backend 参数 | demux（仅 min_matched_read_pairs）、high_cph、bismark（不含 library_type） |

| 设置类别 | 配置路径与默认值 | 用途与适用范围 |
| --- | --- | --- |
| 路线与参考 | `schema_version: 3`、`species`、`references.*`、`analysis.protocol` | 选择协议及物种/reference；仅 Cabernet/SRD 可选 backend |
| CpG 坐标 | `analysis.cpg_coordinate_base: 1`；可选整数 `0` | CpG、RawAdata 实际位点与 SingleCpG 共用；原生 BED/VCF/BAM 不变 |
| Demux | `common.demux.min_matched_read_pairs: 10`；`dna_w_spacer_len: 0` | 阈值可覆盖；Droplet spacer 固定，不接受该覆盖字段 |
| Conversion QC | `common.high_cph.excluded_contigs: [pUC19, lambda, chrM]` | Cabernet/SRD/Droplet 固定统计并剔除，无动作开关 |
| BISCUIT | `library_type: non_directional`；`high_cph_retention_threshold: 0.7`；`generate_snp: false` | 位于 `common.biscuit`，只由有效 BISCUIT 路线消费 |
| Bismark | `library_type: non_directional`；`local_alignment: false`；`non_conversion_percentage_cutoff: 70`；`non_conversion_minimum_count: 5` | 位于 `common.bismark`；Droplet library type 固定 |
| 本地运行 | `runtime.local_executor_cores: 32`；`latency_wait_seconds: 120` | 并发 CPU 预算与文件可见性等待，不改变科学参数 |
| Slurm | `jobs: 30`；`controller_cores: 2`；`qos: huge`；`node_tmpdir: /tmp`；`account/partition` | 位于 `slurm`；qos 为 null 时使用集群默认 |
| 最终 BAM | `retention.keep_final_bam: false` | true 保留当前 backend 的最终 BAM；TAPS/BISCUIT 有 BAI，Bismark 无 BAI；持锁成功发布时清理不再要求保留的 BAM，SRD RNA BAM 独立保留 |

> `common.biscuit.generate_snp` 是 BISCUIT-only SNP 输出开关，默认 `false`；选定 Cabernet/SRD 的有效 backend 为 Bismark 时必须关闭，可用该协议同名叶子覆盖公共值；TAPS/Droplet 不读取此分组。TAPS 的 SNP VCF 与 Non-CpG QC 默认生成，不受此开关控制，唯一科学契约见 [实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/dna_pipeline.py) §5。

以上公共参数均可在适用协议模块的同名路径覆盖，例如 `srd.demux.min_matched_read_pairs`；TAPS/Droplet 的 backend 及 Droplet 的 spacer/library_type 为内部固定策略，不接受这些 config 字段；common 中的 spacer/library_type 仅用于可配置的孔板协议。

### 运行时参数归属

| 配置 | 执行器用途 |
| --- | --- |
| `runtime.local_executor_cores` | `--cores`，受 affinity/Slurm CPU 预算限幅 |
| `runtime.latency_wait_seconds` | `--latency-wait` |
| `slurm.jobs/controller_cores` | 所有阶段共用的作业并发上限 / controller 本地核数；不改变 demux 块大小 |
| `slurm.partition/qos/account` | Slurm 提交参数 |
| `slurm.node_tmpdir` | 节点 scratch 根 |

这些设置不改变科学参数或输出 schema。资源取值见[固定策略与资源](#固定策略与资源)；重试、CPU 预算及 backend Python 选择见
[实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/dna_pipeline.py)、[实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/doctor.sh)；
Droplet配对暂存/calling/重放策略见 [实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/dna_pipeline.py)。
FASTQ 命名、manifest 选样与 refresh 操作见 [README §5](../README.md#5-导入与选样)，不增加分析参数。

## 完整参数表

表中的 Bismark `dedup` 与原生 `extract` 行适用于 Cabernet/SRD；Droplet 共用 alignment 与
non-CG filter，随后使用下文列出的位点 UMI caller，不执行这两项原生操作。

### CpG / 输出

| 路径/阶段 | 参数 | 软件原生默认 | Pipeline 默认 | 解读与约束 |
| --- | --- | --- | --- | --- |
| CpG / 输出 | `analysis.cpg_coordinate_base` † | 不适用 | 1 | 整数 1 为 1-based，0 为 0-based；CpG pos、RawAdata 实际位点与 SingleCpG 标识共用。sealed 记录坐标身份；下游按 base 适配分箱与 BED 区间，原生 BED/VCF/BAM 不变。 |

### 共用 / Cutadapt

| 路径/阶段 | 参数 | 软件原生默认 | Pipeline 默认 | 解读与约束 |
| --- | --- | --- | --- | --- |
| 共用 / Cutadapt | R1 adapter | 无 | 孔板：`CTGTCTCTTATACACATCT`；Droplet：`N{9}CTGTCTCTTATACACATCT` | R1 3′ adapter 序列；在 alignment 前物理剪除已知技术序列。 Droplet 读穿另一端时连同相邻 9 bp 一并剪除。 |
| 共用 / Cutadapt | R2 adapter | 无 | 孔板：`N{9}CTGTCTCTTATACACATCT`；Droplet：`N{9}CTRTCTCTTATACACATCT` | R2 3′ 技术结构；`N{9}` 与 adapter 一起识别并剪除。 Droplet 的 R2 5′ 端另由 demux 去除 9 bp；此处处理 3′ 读穿。 |
| 共用 / Cutadapt | R1 minimum overlap | 3 | 孔板：3；Droplet：12 | R1 adapter 至少匹配多少碱基才触发 trimming。 Droplet 两端均含 `N{9}`；孔板路线保留原有 overlap。 |
| 共用 / Cutadapt | R2 minimum overlap | 3 | 12 | R2 adapter 最小 overlap。 `N{9}` + 至少 3 个真实 adapter 匹配，因此固定 12。 |
| 共用 / Cutadapt | adapter max error rate `-e` | 0.1 | 0.25 | 错误数除以实际匹配 adapter 长度；`N` wildcard 不计入分母，允许 indel。共用入口显式传入 `--error-rate 0.25`；改变此值会影响短 overlap 的允许错配数，须验证科学结果。 |
| 共用 / Cutadapt | adapter indels | 允许 | 允许 | adapter 匹配时是否允许 insertion/deletion。 保持默认，避免不必要漏剪。 |
| 共用 / Cutadapt | adapter search times `-n` | 1 | 1 | 每条 read 最多搜索多少次 adapter。 当前只需要一次明确的 3′ adapter 识别。 |
| 共用 / Cutadapt | adapter wildcards | 启用 | 启用 | adapter 中 IUPAC/`N` 是否作为 wildcard。 R2 `N{9}` 必需。 |
| 共用 / Cutadapt | read wildcards | 不匹配 | 关闭 | read 中 wildcard 是否参与 adapter wildcard 匹配。 保持关闭，避免增加低质量模糊匹配。 |
| 共用 / Cutadapt | adapter action | trim | trim | 找到 adapter 后执行的动作。 固定 trim。 |
| 共用 / Cutadapt | quality cutoff `-q` | 不做 quality trimming | `20,20` | 按 baseQ 对 read 两端做质量剪切。 Q20 固定，避免低质量末端进入 alignment。 |
| 共用 / Cutadapt | R2 独立 quality cutoff `-Q` | 未设置，继承 `-q` | 继承 `-q` | 允许 R2 使用不同 quality cutoff。 当前 R1/R2 均采用 Q20，不单独配置。 |
| 共用 / Cutadapt | minimum length `-m` | 0 / 无下限 | 20 | trimming 后允许进入 alignment 的最短 read 长度。 固定 20 bp。 |
| 共用 / Cutadapt | cores `-j` | 1 | Pipeline 动态分配 | Cutadapt 并行线程数。 performance-only；资源策略内部管理。 |

### BISCUIT / align

| 路径/阶段 | 参数 | 软件原生默认 | Pipeline 默认 | 解读与约束 |
| --- | --- | --- | --- | --- |
| BISCUIT / align | library mode `-b` † | 0 non-directional | `non_directional → 0` | bisulfite library orientation；决定允许的 strand model。 设置错误会系统性影响 mapping。 |
| BISCUIT / align | strand search `-f` | 0 both | 0 | 进一步限制 BSW/BSC strand 搜索。 保持 both，不额外收窄。 |
| BISCUIT / align | minimum seed length `-k` | 19 | 19 | FM-index seeding 的最短 seed 长度。 变小更敏感但更慢。 |
| BISCUIT / align | band width `-w` | 100 | 100 | banded extension 的搜索宽度。 主要影响复杂 indel 搜索。 |
| BISCUIT / align | X-dropoff `-d` | 100 | 100 | extension 分数下降到何种程度后停止扩展。 |
| BISCUIT / align | internal seed ratio `-r` | 1.5 | 1.5 | 控制长 seed 内部重新 seeding。 底层敏感度参数，不暴露。 |
| BISCUIT / align | third-round seed occurrence `-y` | 20 | 20 | 第三轮 seed 搜索的 occurrence 门槛。 底层搜索策略。 |
| BISCUIT / align | adapter R1 `-J` | 无 | 未用 | BISCUIT 内部 adapter 处理。 已由 Cutadapt 完成，避免重复 trimming。 |
| BISCUIT / align | adapter R2 `-K` | 无 | 未用 | 同上。 已由 Cutadapt 完成。 |
| BISCUIT / align | end baseQ trim `-z` | 0 | 0 | aligner 内部按 baseQ 再剪 read 末端。 Cutadapt 已做 Q20，不重复。 |
| BISCUIT / align | extra 5′ clip `-5` | 0 | 0 | alignment 前强制剪掉 read 5′ 端固定碱基。 默认不做；改变有效序列坐标。 |
| BISCUIT / align | extra 3′ clip `-3` | 0 | 0 | alignment 前强制剪掉 read 3′ 端固定碱基。 默认不做。 |
| BISCUIT / align | max seed occurrence `-c` | 500 | 500 | 高重复 seed 的 occurrence 上限。 调小更快但损失重复区域敏感度。 |
| BISCUIT / align | chain drop fraction `-D` | 0.5 | 0.5 | 弱 alignment chain 的 pruning 阈值。 底层 mapper 参数。 |
| BISCUIT / align | minimum seeded bases `-W` | 0 | 0 | alignment chain 所需的最少 seeded bases。 |
| BISCUIT / align | mate rescue rounds `-m` | 50 | 50 | PE mate rescue 的最大尝试程度。 更高更慢。 |
| BISCUIT / align | skip mate rescue `-S` | OFF | OFF | 完全关闭 mate rescue。 关闭 rescue 会降低 paired mapping recovery。 |
| BISCUIT / align | skip pairing `-P` | OFF | OFF | 不进行 paired-end pairing。 PE 数据固定 OFF。 |
| BISCUIT / align | discard exact matches `-e` | OFF | OFF | 是否丢弃 full-length exact matches。 生产分析无理由开启。 |
| BISCUIT / align | barcode/UMI read-name extraction `-9` | OFF | OFF | 从 read name 解析 barcode/UMI。 当前 barcode 不通过该机制传递。 |
| BISCUIT / align | match score `-A` | 1 | 1 | alignment match reward。 改变会重塑整个评分系统。 |
| BISCUIT / align | mismatch penalty `-B` | 2 | 2 | 非 bisulfite conversion mismatch 的 penalty。 底层灵敏度参数，保持默认。 |
| BISCUIT / align | gap open `-O` | 6,6 | 6,6 | insertion/deletion 开启 penalty。 |
| BISCUIT / align | gap extension `-E` | 1,1 | 1,1 | gap 延伸 penalty。 |
| BISCUIT / align | clipping penalty `-L` | 10,10 | 10,10 | read-end clipping penalty。 |
| BISCUIT / align | unpaired penalty `-U` | 17 | 17 | PE 不能形成合理 pair 时的 penalty。 |
| BISCUIT / align | minimum output score `-T` | 30 | 30 | 低于该 alignment score 的记录不输出到 BAM。 保留软件默认。 |
| BISCUIT / align | max XA hits `-g` | 5,5 | 5,5 | alternative alignment tag 最多保留数量。 主要影响辅助 tag/文件体积。 |
| BISCUIT / align | output all SE/unpaired alignments `-a` | OFF | OFF | 是否输出所有 SE/unpaired alternative alignments。 关闭，避免扩张 ambiguous 记录。 |
| BISCUIT / align | insert-size model `-I` | 自动推断 | 自动推断 | paired-end insert-size distribution。 由数据自动推断。 |
| BISCUIT / align | threads `-@` | 1 | 动态，当前最多15 | BISCUIT align 线程。 performance-only。 |

### BISCUIT / Dupsifter

| 路径/阶段 | 参数 | 软件原生默认 | Pipeline 默认 | 解读与约束 |
| --- | --- | --- | --- | --- |
| BISCUIT / Dupsifter | duplicate removal | 工具默认不删除，只标记 | **物理删除 duplicate (`--remove-dups`)** | 识别 PCR duplicate；保留每个 duplicate family 的代表 molecule，并从输出 BAM 删除其余 duplicate records。 最终 CpG 与“保留 duplicate flag 后由 pileup 排除”使用同一类非重复 molecules；优点是 BAM 更小、后续 sort/bsconv/pileup I/O 更低。Duplicate 统计经 Dupsifter 的 MultiQC custom content 进入 QC。 |

### 共用 / high-CpH

| 路径/阶段 | 参数 | 软件原生默认 | Pipeline 默认 | 解读与约束 |
| --- | --- | --- | --- | --- |
| 共用 / high-CpH | 处理动作 | 无，Pipeline 自定义 | 固定「统计并剔除」 | 分类永远执行，CpG 恒消费过滤代 BAM；无 off/assess 用户开关。 Cabernet/SRD/Droplet 行为一致；被筛除 BAM 仅 SRD 恒定保留到 02_work/RNA/BAM。 |

### BISCUIT / high-CpH

| 路径/阶段 | 参数 | 软件原生默认 | Pipeline 默认 | 解读与约束 |
| --- | --- | --- | --- | --- |
| BISCUIT / high-CpH | retention threshold † | 无，Pipeline 自定义 | 0.7 | read-pair CpH retention **严格大于**该值时判为 high-CpH；恰好等于阈值不命中。 越低越严格；越高越宽松。 |
| BISCUIT / high-CpH | excluded contigs † | 无，Pipeline 自定义 | `pUC19, lambda, chrM` | 不参与 high-CpH/non-conversion 判断的 control/mitochondrial contigs。 避免 control biology 扰动 conversion QC。 |

### BISCUIT / pileup

| 路径/阶段 | 参数 | 软件原生默认 | Pipeline 默认 | 解读与约束 |
| --- | --- | --- | --- | --- |
| BISCUIT / pileup | threads `-@` | 3 | 1 | pileup 并行线程；资源策略各 tier 固定单线程。 performance-only；生产 pipeline 输出按单 worker 耗时的 walltime 分级预留。 |
| BISCUIT / pileup | minimum baseQ `-b` | 20 | 20（原生默认，未显式传入） | 单个 C/T observation 最低 base quality。 Q20；低质量 base 不贡献 methylation evidence。 |
| BISCUIT / pileup | minimum MAPQ `-m` | 40 | 40 | alignment 最低 mapping quality。 MAPQ<40 的 ambiguous/低置信 alignment 不参与 CpG calling。 |
| BISCUIT / pileup | minimum AS `-a` | 40 | 40 | 参与 pileup 的最低 alignment score。 与 MAPQ 共同形成高质量 evidence gate。 |
| BISCUIT / pileup | max C retention `-t` | 999999 | 999999 | pileup 内部的 C-retention 绝对数量过滤。 保持关闭状态；已有独立 high-CpH 比例过滤，避免重复 gate。 |
| BISCUIT / pileup | minimum read length `-l` | 10 | 20 | pileup 允许贡献 evidence 的最短 read。 与 Cutadapt minimum length 对齐为20。 |
| BISCUIT / pileup | minimum distance from 5′ `-5` | 3 | 3 | read 5′ 端前 N bp 不贡献 cytosine call。 BISCUIT 原生 read-end evidence contract；不修改 alignment。 |
| BISCUIT / pileup | minimum distance from 3′ `-3` | 3 | 3 | read 3′ 端后 N bp 不贡献 cytosine call。 同上。 |
| BISCUIT / pileup | redistribute ambiguous Y/R | ON | ON | 处理 SNP 场景下 ambiguous Y/R allele evidence。 |
| BISCUIT / pileup | filter secondary | ON | ON | 排除 secondary mappings。 固定 ON。 |
| BISCUIT / pileup | overlap double-count | OFF | OFF | 是否让 overlapping mates 对同一 C 重复计数。 固定 OFF，避免人为抬高 coverage。 |
| BISCUIT / pileup | filter duplicates | ON | ON | Dupsifter 已物理去重；pileup继续排除仍带 duplicate flag 的 reads。 |
| BISCUIT / pileup | filter improper pairs | ON | ON | 排除 improper paired alignments。 固定 ON。 |
| BISCUIT / pileup | max NM `-n` | 999999 | 999999 | 最大 non-conversion mismatch 数。 已有 MAPQ/AS，不叠加重复 mismatch gate。 |

### BISCUIT / pileup SNP

| 路径/阶段 | 参数 | 软件原生默认 | Pipeline 默认 | 解读与约束 |
| --- | --- | --- | --- | --- |
| BISCUIT / pileup SNP | error rate `-E` | 0.001 | 0.001 | SNP/genotype likelihood 的 sequencing error prior。 仅 SNP 模型。 |
| BISCUIT / pileup SNP | mutation rate `-M` | 0.001 | 0.001 | germline mutation prior。 仅 SNP 模型。 |
| BISCUIT / pileup SNP | somatic rate `-x` | 0.001 | 默认 | somatic mutation prior。 当前非 somatic workflow。 |
| BISCUIT / pileup SNP | contamination `-C` | 0.01 | 默认 | somatic contamination model。 当前不使用。 |
| BISCUIT / pileup SNP | heterozygous prior `-P` | 0.333 | 默认 | heterozygous genotype prior。 |
| BISCUIT / pileup SNP | homozygous prior `-Q` | 0.333 | 默认 | homozygous genotype prior。 |

### BISCUIT / vcf2bed

| 路径/阶段 | 参数 | 软件原生默认 | Pipeline 默认 | 解读与约束 |
| --- | --- | --- | --- | --- |
| BISCUIT / vcf2bed | context `-t` | CG | CG | 只导出 CpG context。 固定 CpG。 |
| BISCUIT / vcf2bed | min coverage `-k` | 1 | 1 | 最终 CpG 输出最低 coverage。 单细胞保留 singleton CpG。 |
| BISCUIT / vcf2bed | sample `-s` | FIRST | FIRST | 多样本 VCF 中选择 sample。 当前单 sample VCF。 |
| BISCUIT / vcf2bed | show context `-e` | OFF | OFF | 额外输出 context 信息。 SnapATAC2 canonical 四列表不需要 context，固定 OFF。 |
| BISCUIT / vcf2bed | Beta-M-U `-c` | OFF | ON | 输出 beta / methylated / unmethylated counts。 canonical CpG 表需要。 |

### BISCUIT / biscuit_qc

| 路径/阶段 | 参数 | 软件原生默认 | Pipeline 默认 | 解读与约束 |
| --- | --- | --- | --- | --- |
| BISCUIT / biscuit_qc | QC cytosine `vcf2bed -e -t c` | `-e` OFF，`-t cg`，`-k 1` | 显式 `-e -t c`，不传 `-k` | 导出带 context 的全胞嘧啶观测，服务 `totalBaseConversionRate`。 公开 CpG 使用 `-t cg -k 1 -c`，两条输出服务不同 consumer。 |

### Bismark / align

| 路径/阶段 | 参数 | 软件原生默认 | Pipeline 默认 | 解读与约束 |
| --- | --- | --- | --- | --- |
| Bismark / align | library type † | directional | non_directional | Bismark library orientation；non-directional 搜索完整 bisulfite strand space。 设置错误会系统性影响 mapping。 |
| Bismark / align | backend | Bowtie2 | Bowtie2 | 底层 short-read mapper。 当前 PE bisulfite short reads 固定 Bowtie2。 |
| Bismark / align | alignment/index mode † | end-to-end / per-strand | `local_alignment=false`：combined-index end-to-end；`true`：faithful CT/GA + local | 决定 soft-clipping 与实际使用的 index generation。 两种模式互斥；默认 end-to-end，runtime gate 只检查当前模式所需 index。 |
| Bismark / align | seed mismatches（Bowtie2 `-N`；Rust CLI `--seedmms`） | 0 | **0** | seed 允许的 mismatch 数量；固定 upstream 默认。 |
| Bismark / align | seed length（Bowtie2 `-L`；Rust CLI `--seedlen`） | 20 | **20** | Bowtie2 multiseed 长度；变短通常更敏感但更慢。 |
| Bismark / align | minimum insert `-I` | 0 | 0 | paired fragment 最小 insert。 |
| Bismark / align | maximum insert `-X` / `--maxins` | 500 | **1000** | paired fragment 最大 insert（bp）。 所有 Bismark 协议与 alignment mode 固定显式传入 `--maxins 1000`，见 [实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/dna_pipeline.py) §5。 |
| Bismark / align | no-mixed | ON | ON | pair 失败时是否尝试单 mate mapping。 固定 ON=禁止 mixed，维持 paired evidence contract。 |
| Bismark / align | no-discordant | ON | ON | 是否接受 discordant pairs。 固定 ON=禁止 discordant。 |
| Bismark / align | dovetail | ON | ON | 短 fragment 中 mates 交叉/重叠时仍可视为 concordant。 对短片段重要，固定 ON。 |
| Bismark / align | Bowtie2 effort `-D` | 15 | **15** | 连续 seed-extension 失败次数上限。 |
| Bismark / align | reseed `-R` | 2 | **2** | 重复 read 的 reseeding 次数。 |
| Bismark / align | score minimum | mode-dependent | local：`G,20,8`；end-to-end：`L,0,-0.2` | Bowtie2/Bismark 最低 alignment score function。 与 `local_alignment` 自动联动，禁止 `--local + L,...` 或 combined + local 的非法组合。 |
| Bismark / align | read gap penalty `--rdg` | 5,3 | 5,3 | read-side gap open/extend penalty。 |
| Bismark / align | reference gap penalty `--rfg` | 5,3 | 5,3 | reference-side gap open/extend penalty。 |
| Bismark / align | ignore qualities | ON | ON | alignment mismatch scoring 是否忽略 FASTQ qualities。 保持 Bismark/Bowtie2 bisulfite 默认语义。 |
| Bismark / align | temporary conversion gzip | OFF unless specified | OFF unless specified | Bismark 转换 FASTQ 临时文件是否 gzip。 按最终决策不显式启用；临时数据位于 node-local scratch。 |
| Bismark / align | temp dir | 输入/输出附近默认 | node-local scratch | combined spill / 临时文件位置。 性能/I/O 参数，固定 node-local scratch。 |
| Bismark / align | Bowtie2 threads `-p` | 1 | rule总CPU÷并发instance数 | 每个Bowtie2 instance的线程数。 faithful directional/PBAT有2个、non-directional有4个并发instance；combined模式同一时刻只有1个instance，禁止把rule总CPU原样传给faithful模式的每个instance。 |
| Bismark / align | Bismark `--multicore` | 1 | 1 | 同时启动多个 Bismark worker。 Rust 版本固定1，避免重复 index/RAM。 |

受管 Rust Bismark v3.1.0 的 seed 参数使用长选项 `--seedmms`/`--seedlen`，不依赖大写短选项支持；升级时核对实际 CLI 并做 mapping/MAPQ/CpG 回归。

### Bismark / dedup

| 路径/阶段 | 参数 | 软件原生默认 | Pipeline 默认 | 解读与约束 |
| --- | --- | --- | --- | --- |
| Bismark / dedup | SE/PE detection | auto | paired | 决定 duplicate 的 SE/PE 坐标定义。 孔板路线固定 paired；Droplet 使用逐胞嘧啶 UMI 共识（[实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/dna_pipeline.py) §7）。 |
| Bismark / dedup | output BAM | ON | ON | dedup 输出 BAM。 |
| Bismark / dedup | barcode/UMI `--barcode` | OFF | OFF | 将 UMI/barcode 纳入 duplicate identity。 孔板 cell barcode 不是 molecular UMI；Droplet 使用位点 UMI 共识，仍不开此 flag。 |
| Bismark / dedup | multiple files | OFF | OFF | 多个 BAM 是否视为同一 library 一起 dedup。 当前每个 cell 独立。 |
| Bismark / dedup | compression threads `--parallel` | 通常1 | `min(job threads,4)` | Rust dedup BGZF 压缩线程，属于内部资源策略。 |

### Bismark / non-CG filter

| 路径/阶段 | 参数 | 软件原生默认 | Pipeline 默认 | 解读与约束 |
| --- | --- | --- | --- | --- |
| Bismark / non-CG filter | absolute threshold | 3 | 未使用 | 按 methylated non-CG 绝对个数判断 conversion failure。 受 read 长度/CpH 数量影响，禁用；使用 percentage 模式。 |
| Bismark / non-CG filter | percentage cutoff † | 未启用 | 70% | 任一 mate 的 non-CG retention 达到该百分比且满足 minimum count 时移除整对。 每个 mate 独立判定，越低越严格。 |
| Bismark / non-CG filter | minimum count † | 5 | 5 | 单个 mate 的 informative non-CG C 数量下限；须与 percentage cutoff 同时满足才移除整对，不合并两端分母。 |
| Bismark / non-CG filter | consecutive | OFF | OFF | 是否要求异常 methylated non-CG 连续出现。 不叠加额外定义；固定 OFF。 |
| Bismark / non-CG filter | paired | auto | paired | 保留/剔除的单位为整对，异常阈值按 mate 判断。 PE 数据固定 paired。 |

### Bismark / extract

| 路径/阶段 | 参数 | 软件原生默认 | Pipeline 默认 | 解读与约束 |
| --- | --- | --- | --- | --- |
| Bismark / extract | input SE/PE | auto-detect | auto | extractor识别输入为SE/PE。 输入 BAM 已有 paired flags。 |
| Bismark / extract | no-overlap | ON for PE | ON | R1/R2 overlap 中同一 cytosine 只计一次。 固定 ON，避免 double-count。 |
| Bismark / extract | ignore R1 5′ `--ignore` | 0 | 0 | R1 5′ 前N bp不贡献 methylation call。 不是 Tn5 9-bp gap 修复；仅 M-bias 证明端部偏倚时才考虑修改。 |
| Bismark / extract | ignore R2 5′ `--ignore_r2` | 0 | 0 | R2 5′ 前N bp不贡献 methylation call。 默认0；技术序列已由 Cutadapt 处理。 |
| Bismark / extract | ignore R1 3′ | 0 | 0 | R1 3′ 后N bp不贡献 methylation call。 默认0；M-bias-driven。 |
| Bismark / extract | ignore R2 3′ | 0 | 0 | R2 3′ 后N bp不贡献 methylation call。 默认0；M-bias-driven。 |
| Bismark / extract | comprehensive | OFF | ON | 将四种 strand 的 context outputs 合并。 减少临时文件；不改变 pooled CpG 信息。 |
| Bismark / extract | merge non-CpG | OFF | OFF | CHG/CHH 保持分开供原生报告和 QC 使用；公开定量输出为 CpG。 |
| Bismark / extract | report | ON | ON | 生成 extraction report。 QC需要。 |
| Bismark / extract | gzip | split output默认OFF；coverage默认gz | ON | 压缩 extraction 输出。 降低磁盘/I/O。 |
| Bismark / extract | M-bias only | OFF | OFF | 只计算 M-bias 而不做正式 extraction。 production关闭。 |
| Bismark / extract | M-bias off | OFF | OFF | 关闭 M-bias report。 固定OFF，保留端部偏倚QC。 |
| Bismark / extract | `--parallel` | 1 | 1 | Rust extractor worker数；BAM decode和gzip output自带并行，不叠加workers。申请资源见[Backend资源档位](#backend-资源档位)。 |

### Bismark / extract bedGraph

| 路径/阶段 | 参数 | 软件原生默认 | Pipeline 默认 | 解读与约束 |
| --- | --- | --- | --- | --- |
| Bismark / extract bedGraph | bedGraph | OFF | ON | 生成 bedGraph/coverage。 孔板路线直接生成 `.bismark.cov.gz`，作为按配置 base 写出 canonical CpG 的计数来源；context 沿用 XM，保留 CIGAR deletion 边界形成的 observed CpG，不额外执行 reference-only CpG 筛选。 |
| Bismark / extract bedGraph | cutoff | 1 | 1 | coverage 输出最低 observation count。 单细胞保持1。 |
| Bismark / extract bedGraph | CX context | OFF | OFF | 是否输出所有 cytosine context。 只输出CpG，固定OFF。 |
| Bismark / extract bedGraph | zero-based extra output | OFF | OFF | 额外生成 zero-based coverage。 |
| Bismark / extract bedGraph | buffer size | 工具自动/default sort策略 | 默认 | bedGraph 排序缓冲。 performance-only。 |

### 共用 / processing

| 路径/阶段 | 参数 | 软件原生默认 | Pipeline 默认 | 解读与约束 |
| --- | --- | --- | --- | --- |
| 共用 / processing | 重复证据处理 | 无软件统一默认 | 按路线固定 | BISCUIT/Dupsifter 与 Bismark 孔板物理去重；TAPS 标记并保留 BAM 记录，调用时排除 duplicate；Droplet 保留 pairs、按胞嘧啶/UMI 汇集共识。计数单位见[跨backend QC](#cross-backend-qc-计数口径)。 |
| 共用 / processing | zero mapped policy | Pipeline 自定义 | empty | 全未比对或 trim 后全空的活跃 cell 按当前 backend 输出契约发布合法空产物、带表头 CpG 和 QC；合法零 primary reads 的 mapping 为 0%，无 CpH 观测为 NaN。 |

## 固定策略与资源

协议科学路线、恢复和调用边界由当前源码实现；参数取值与计数口径在本文维护，相关边界索引如下：

| 主题 | 维护契约 |
| --- | --- |
| 科学 helper 身份、实际消费者与增量重算 | [实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/dna_pipeline.py) |
| DNA/RNA manifest、外部关联与 Notebook 入口 | [实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/dna_pipeline.py)、[下游操作](../downstream/README.md#关联-rna_darlin-结果) |
| 锁、重试、CPU/内存预算、日志与执行器 | [实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/dna_pipeline.py) |
| inventory、稳定 stage、temp/BAM 与 demux 恢复 | [实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/dna_pipeline.py) |
| BISCUIT / Bismark / TAPS 科学路线与reference | [实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/dna_pipeline.py)；资源取值见[Backend资源档位](#backend-资源档位) |
| 环境与 Doctor fixture/调度/验证边界 | [实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/doctor.sh) |
| Droplet 结构、calling、白名单、位点 UMI、解复用与配套 RNA | [实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/dna_pipeline.py)；资源/I/O取值见[下表](#droplet-内部资源与-io) |
| 坐标、canonical CpG、single-CpG 矩阵 | [实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/downstream/downstream_qc.py) |
| MultiQC 展示、资源计算边界与下游缓存 | [实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/downstream/downstream_qc.py) |

以下资源计算值只供维护，不暴露为 config：

- **孔板 demux**：按 snapshot 中压缩 FASTQ 总字节分级；直接读取当前 snapshot 的 `size_bytes`，不回退读取实时路径。
  `≤5 GiB→S(4 GiB/120 min)`，`≤50 GiB→M(8 GiB/360 min)`，
  `>50 GiB→L/XL(12 GiB/720–1440 min)`；该画像不按字节继续外推内存。
- **controller 小任务**：`initial_cell_manifest` 为 512 MiB/10 min，
  `final_sample_manifest` 为 2048 MiB/30 min，`delivery_ready` 为 1024 MiB/15 min，
  唯一表为 `CONTROLLER_RESOURCE_REQUESTS`；final 最多 8 threads、每批最多 16 cells，其余单线程。
  重试与并发边界见[实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/dna_pipeline.py)。
- **MultiQC**：按已解析声明 input 数 N，首次 `mem=max(1024,ceil(1.5×(512＋0.075N)))` MiB，
  `runtime=max(60,ceil(N/600))` min；计算及重试边界见[实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/downstream/downstream_qc.py)。

### Backend 资源档位

下表按 S/M/L/XL 给出首次资源；重试仅增加内存，规则见[实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/dna_pipeline.py)。分组失败、部分重跑和 scratch 的恢复/验收约束见[实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/dna_pipeline.py)。

| Bismark阶段 | CPU | 内存GiB | runtime/适用范围 |
| --- | --- | --- | --- |
| combined alignment | 8/12/12/16 | 16/16/16/16 | 默认end-to-end模式 |
| faithful local alignment | 12/12/16/16 | 24/24/32/32 | 仅所选local模式 |
| SRD独立alignment（两模式） | 按所选模式 | 按所选模式 | alignment120/180/270/900min，Cutadapt30/45/60/90min，合计150/225/330/990min |
| 孔板/Droplet extraction | 8 | 16 | attempt2/3为24/32GiB |
| Cabernet/Droplet Bismark组 | 串行阶段最大值 | 串行阶段最大值 | 四规则完整600min（各150），部分重跑≤600 |

| TAPS阶段 | CPU/内存 | runtime或scratch |
| --- | --- | --- |
| alignment含Cutadapt | Slurm16核/16GiB；local偏好8/12/16/16核 | 150/225/300/330min |
| calling含Python后处理 | Slurm16核/4/8/12/16GiB；local偏好4核且只受总预算截断 | 60/120/180/240min；Rastair用实际分配线程 |
| taps_cell组 | 每cell两规则，串行最大16核/16GiB，重试24/32GiB | 完整210/345/480/570min；部分重跑算剩余成员 |
| 两阶段scratch | 独立于mem/runtime | 16/32/64/64GiB；空间余量及降低门禁见[实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/dna_pipeline.py) |

### Droplet 内部资源与 I/O

以下为当前内部取值，不进入analysis config。路线、整数守恒、FD与失败传播、提交/恢复边界见[实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/dna_pipeline.py)；实际CPU限幅和仅内存递增的重试规则见[实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/dna_pipeline.py)。调整须满足[实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/dna_pipeline.py)。

| 项目 | 当前取值与适用范围 |
| --- | --- |
| 作业组与阶段 | 单库droplet_run完整32CPU/首次至少16GiB/2890min；calling与demux各32CPU/1440min，demux首次16GiB；文库身份规则512MiB/10min |
| calling首次内存 | `max(16GiB, ceil(原R1计数请求/4) + 4GiB)`；原XL请求为`ceil((2 + 0.048×R1输入GiB) × 1.3 × 1024)` MiB，其余tier沿实现 |
| stage并行与分块 | rapidgzip0.16.0每端4解码线程、19计算workers（32CPU）；100k pairs/block，输入池为计算worker数+1且至少2，32CPU时20块；64个独立增长计数分区，每barcode保留34项质量计数 |
| replay准备 | 最多4blocks/估2GiB窗口，复用64个Mutex保护的spool读句柄；FD足够时启用4个确定性cell分区writer |
| 回写触发与并发 | 有效memory.max/high的四分之一，限64MiB–6GiB；同时以512MiB未结清写入触发，不等待内核限流。4个I/O workers每次处理≤1024个非空积压弱引用句柄，后台≤4个强引用占原FD预留 |
| FD缓存 | 预留普通64＋spool64FD，剩余每对输出2FD；不能全开时≤128对LRU并受可用FD限幅。其他协议128对/default BufWriter |
| 写缓冲与批量 | Droplet总缓冲为回写阈值一半，限512MiB–2GiB（未知512MiB），单流8–256KiB；paired批100k，count/分子/孔板批20k，gzip level1 |
| spool读缓存 | replay批次读完每8MiB回收已消费完整页；尾批包含不足8MiB的完整页 |
| 文件copy | Linux≤4workers，各复用1MiB缓冲，每8MiB尽力释放源/目标已消费缓存；无缓存提示的平台≤8native workers |
| scratch容量 | 估2×spool压缩量＋512MiB，并留10%余量；空间不足或文件系统不适用时按当前实现回退共享目录 |



### TAPS 固定参数与第三方依据

Rastair 2.2.0 的固定调用参数及字段如下；reference事务、文件提交和资源边界见[实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/dna_pipeline.py)。性能小样本不能证明完整字段等价。

| 阶段 | 参数、字段与科学定义 |
| --- | --- |
| BWA先验 | taps-bwa-insert-model取当前cell前10,000完整trimmed pairs跑原生FR；-I为mean/SD/max/min，保持mean/SD/min，max取max(1000,原生max)。仅有限正mean/SD及有效边界才传；无可靠FR/零SD用自动推断，命令失败传播。scratch采样/日志回收，不固定全库均值，不冒充比对真值 |
| Rastair调用 | 环境固定2.2.0并核原生version；`call --segment-max-length 100000 --guess-read-orientation -q 20 -Q 30 -f 3 -F 3852 --all --m-min-depth 1 --vcf-info-fields AS_SB,AS_SS_BQ,AS_SS_MQ,PIR,AF`写native BED/VCF；不用`--bed-include-empty`/`--keep-overlapping-reads`/`--unpaired`，ML默认0.50、pileup/methylation最大深度1000、`--v-min-depth 3`；缺证据/TG-CA平票用原生可复现分配，不当已测链方向 |
| Non-CpG定义 | marked BAM完整primary QNAME pair，MAPQ≥20/BQ≥30/include3/exclude3852，proper/nonduplicate合格mates；宿主chr*去chrM/chrMT，不混controls。高质量参考CpG相邻TG/CA证据唯一类型判OT/OB，两类都有或都无则ambiguous排除，不能用Rastair随机分配替代真实判链 |
| Non-CpG计数 | 正/反参考三碱基CHG/CHH，OT=T/(C+T)，OB=A/(G+A)，未知context不计；重叠坐标取高BQ，等质冲突不计，每pair/胞嘧啶一次，无位点深度上限。总modified/有效观测×100，零分母NaN；candidate/oriented无需有CHG/CHH观测。source=`taps_oriented_host_chg_chh_conversion`，保分项/分母/pair计数，不解释为真实mCH/转化效率/cDNA率 |
| pair funnel | BAM保留记录与调用集合分开。`backend_accepted_pairs`=至少一primary mate proper、双端mapped、非QC fail、MAPQ≥20（忽略duplicate）；只有duplicate满足则归`duplicate_pairs`；`postdedup_pairs/final_retained_pairs`至少一非dup mate可调用，不要求覆盖高质CpG。`both_primary_mapped_pairs`只计接受集合双mapped，`one_primary_mapped_pairs=0`；拒绝分`unmapped_pairs/discordant_pairs/quality_rejected_pairs`，不伪造ambiguous/no_primary |
| controls | lambda假阳性/pUC19转换均mod/(mod+unmod)，同报位点/观测，零覆盖未评估；实物修饰状态先核，不以子集原始转换率替代真实化学效率 |
| SNP解释 | PASS/`--v-min-depth 3`不保证FORMAT/DP≥3或高GC/QUAL，不额外按DP/GC/QUAL筛选；缺乏独立链方向/基因型真值，候选不宣称高置信，未调用非纯合参考。BISCUIT generate_snp开关不控制TAPS |

TAPS 第三方依据：[Rastair 2.2.0 源码](https://github.com/bsbludwig/rastair/tree/v2.2.0)、[链识别说明](https://docs.rastair.com/calling/strand-guesser.html)、[原生 BED](https://docs.rastair.com/formats/bed.html)、[Watchmaker 计数示例](https://github.com/watchmaker-genomics/taps-foundry/blob/main/rastair_summarize.py)、[samtools markdup](https://www.htslib.org/doc/samtools-markdup.html)、[samtools collate](https://www.htslib.org/doc/samtools-collate.html)、[BWA官方手册](https://github.com/lh3/bwa/blob/v0.7.19/bwa.1)、[Snakemake分组](https://snakemake.readthedocs.io/en/stable/executing/grouping.html)。版本行为以受管 2.2.0 的实际 help 与真值测试为准。

### Droplet 固定参数与配套 RNA 关联

执行路线固定为paired stage→calling/evidence→replay。分子证据只消费已提交spool，原始R1/R2由stage读取；replay按原始pair顺序输出。

| 阶段 | 参数、字段与科学定义 |
| --- | --- |
| R1结构 | CB17—UMI12—TSO13—17L17—ME19—gap9—insert；TSO=TTTCTTATATGGG，17L=CGTCCGTCGTTGCTCGT，ME=AGATGTGTATAAGAGACAG。49bp锚点独立评估C→T或G→A，不混模式；额外编辑≤2/indel≤1，零基起点28/29/30，同提取边界合并，不同边界最低编辑并列拒绝，最佳28缺完整CB/UMI拒绝；R1去实际锚点末端＋gap9，R2固定去前9 |
| 原始计数 | paired stage按R1结构计数，同时完整校验R2；观测ACGT拒N，含C保留，候选AGT17且count>100。counts.tsv四列barcode/count/below_q20/q30，后两列为实际CB起点17位逗号计数；多/缺列拒绝，质量单位原始结构通过reads，不是分子 |
| 谷底 | log10(count+1)，200等宽bins，23-bin零填充平滑，峰距≥20；按高度降序/bin顺序取最高两峰间最低谷，平台中心偶数向下、多等深取最前；阈值10**bin_center−1，count≥阈值保留，无有效双峰失败，无回退/expected-cell调参 |
| 错误条码筛除 | 谷底后不重算阈值；须唯一更高计数Hamming1父、父≥5×子，差异位点子Q30计数<谷底、Q<20比例≥20%且≥父同位5倍，并共享精确UMI12＋insert前32中≥3不同UMI及≥3不同insert；缺证据/多父/父也拟筛除保留，不传递合并 |
| 分子证据 | 有候选时Rust `--droplet-spool/--evidence-pairs/--evidence-output`两遍扫完整暂存，子保确定性hash最小≤4096签名、父全扫只查询目标；ACGT完整UMI/insert、精确88-bit编码（hash仅排序），重复不增共享计数，子饱和须记录且缺证据保留；不依赖alignment/Map |

位点 UMI caller 固定 MAPQ≥10、baseQ≥20、BAQ 关闭、深度不截断，保留默认 mate overlap。
`Duplicate_Pair_Rate%` 为不适用的 NaN，分子单位为胞嘧啶/UMI；`keep_final_bam` 保留过滤后
全部 pairs，包含 PCR 复制。参考 CpG、Hamming-1、平票、QC 与配套 RNA 关联见
[实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/dna_pipeline.py)。其他去重策略缺失 pair 去重率仍拒绝 QC 导出。

配置操作见 [README §4](../README.md#4-初始化与配置)；专用 Cutadapt adapter/overlap 见完整参数表。

### 生物学分母与零值解释

| 指标 | 当前计算与不可混淆的含义 |
| --- | --- |
| QC三参考Meth_CpG_Rate% | 分别对宿主/lambda/pUC19已覆盖的canonical胞嘧啶位点取m/(m+u)的算术均值×100，不是将reads计数合并后的比率；零位点为NaN |
| contig归属 | 位点统计与coverage分类共用原始名称：controls仅精确`lambda`/`pUC19`，mt按既有chrM/MT/chrMT规则，宿主为其余chr前缀；其他contig只进total分母。例如`chr_lambda_alt`属宿主，`lambda_decoy`只进total；不按子串或大小写别名归为controls。修正分类只失效CpG统计缓存 |
| SingleCpG X | 同一参考CpG区间内已覆盖C/G位点的fraction均值；例如1/4与9/12得到0.5，不是10/16。缺覆盖为稀疏缺项，有覆盖但全unmethyl为显式0；不把缺项补成观测0 |
| TSS | 同offset有覆盖bin的fraction均值；count为甲基化信号之和，不能当总coverage/read pairs；无覆盖不出行，有覆盖全unmethyl保fraction=0/count=0 |
| TAPS native controls | lambda假阳性/pUC19转换取总mod/总(mod+unmod)，与上述QC位点均值不是同一分母；必须同时报告观测数，不混用为化学效率 |

## Cross-backend QC 计数口径

最终 `sample_manifest.tsv` 的主 QC 采用 **read-pair/QNAME** 统一口径：

```text
cutadapt_input_pairs
  -> trimmed_pairs
  -> backend_accepted_pairs + backend_rejected_pairs
  -> postdedup_pairs = backend_accepted_pairs - duplicate_pairs
  -> final_retained_pairs = postdedup_pairs - high_cph_removed_pairs
```

其中 `both_primary_mapped_pairs + one_primary_mapped_pairs = backend_accepted_pairs`。各 backend 的接受定义以`mapping_policy`记录，不可观测字段见下方路线表，不填伪造的0。backend-native mapping 单独记录 policy/unit，统一 pair 单位不代表 uniqueness/concordance 判定相同。

统一比率按下表计算，不跨阶段替换分母：

| 指标 | 分子 | 分母 | 解释 |
| --- | --- | --- | --- |
| Trim retention | trimmed pairs | Cutadapt input pairs | trimming 后保留比例 |
| Duplicate pair rate | duplicate pairs | backend accepted pairs | 以pair为去重单位的路线适用，Droplet 为 NaN |
| PostDedup retention（`PostDedup_Retention_pct`） | postdedup pairs | backend accepted pairs | 分子/分母×100，不使用 trimmed 作为分母 |
| Final yield | final retained pairs | Cutadapt input pairs | 从 trimming 前到最终调用输入的 pair 产出 |
| High-CpH rate | flagged pairs | assessed pairs | 保留实际评估集合；TAPS 不适用 |

| 路线 | 接受 / 去重单位 | 不可误读的边界 |
| --- | --- | --- |
| Bismark 孔板 | unique concordant pairs；paired-end dedup；policy=`bismark_unique_concordant` | 原生 `ambiguous_pairs` 可观测，`no_primary_pairs` 留空 |
| BISCUIT / Dupsifter | any-primary QNAME pair funnel，允许 one-mate-primary；policy=`biscuit_any_primary_pair_funnel` | 原生互斥 `No_Mapped`/`No_Primary` 分别记为 `unmapped_pairs`/`no_primary_pairs`，两者之和为 `backend_rejected_pairs`；`ambiguous_pairs` 留空，不自行重判 BAM mapping |
| Droplet | alignment/filter 以 pairs 记账，甲基化按胞嘧啶/UMI 共识 | 不移除整对，`duplicate_pairs=0`、`postdedup=accepted`；pair 去重率为NaN，不填0%；位点共识计数单列 |
| TAPS / Rastair | 标记重复，允许进入调用的 pairs | 保留 BAM 含全部记录；final 不等同 BAM 总 pairs；high-CpH 为 NaN |

MultiQC 只对展示抽样，文库概要与科学结果仍覆盖全库；名单、持久 summary 重建和临时文件边界见[实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/downstream/downstream_qc.py)。

Cutadapt adapter QC 分开输出 `R1_Adapter_pct` 与 `R2_Adapter_pct`。`non_cpg_methylation_pct` 对 BISCUIT 来自 CpH retention table，对 Cabernet/SRD Bismark 来自原生 extraction splitting report 的 CHG+CHH 观察，对 Droplet 来自过滤后 CHG/CHH 的胞嘧啶/UMI 共识计数，对 TAPS 采用[判链宿主统计](#taps-固定参数与第三方依据)。均记录 `non_cpg_metric_source`，不作为跨 backend 的统一 hard threshold。

### 下游 QC 字段

`DNAme_QC_Information.csv`按下表分组顺序写20列；完整pair funnel、policy/unit及更细的诊断保留在sealed final manifest与MultiQC。
身份列按字符串读写，保留前导零及合法`NA/nan/None`；数值缺失保持NaN，CSV按round-trip精度复用。

| 列（按写出顺序） | 定义与单位 |
| --- | --- |
| `Sample_ID`、`CloneID`、`PlateID` | sealed细胞/克隆/孔板身份；不能用行号或静默交集替换 |
| `Native_Mapping%` | backend原生mapping百分比，保留其原生单位与分母，不当作跨backend同义指标 |
| `Sample_Unique_CpG_Sites`、`Sample_Meth_CpG_Rate%` | 宿主已覆盖canonical胞嘧啶位点数及逐位点甲基化率均值 |
| `Lambda_Unique_CpG_Sites`、`Lambda_Meth_CpG_Rate%` | 精确lambda contig的同类位点数及均值，非合并reads比率 |
| `pUC19_Unique_CpG_Sites`、`pUC19_Meth_CpG_Rate%` | 精确pUC19 contig的同类位点数及均值 |
| `Non_CpG_Methylation%`、`High_CpH_Flag_Rate%` | 各backend定义的non-CpG统计及flagged/assessed比例；TAPS前者为判链宿主CHG/CHH转换率，后者不适用为NaN |
| `Signal_Host_Rate%`、`Signal_Lambda_Rate%`、`Signal_pUC19_Rate%`、`Signal_mtDNA_Rate%` | 各类canonical CpG coverage(m+u)/全部contig总coverage×100；其他contig仍进入总分母 |
| `Gini_Index` | canonical CpG coverage在有覆盖固定宽bin的Gini；不默认排除chrM/chrY，只在空CpG细胞缺失 |
| `Duplicate_Pair_Rate%`、`Trim_Retention%`、`Final_Pair_Yield%` | 分母按上方pair表；Droplet不以分子数反推pair去重率，该列为NaN |

## Bismark index 构建与 readiness

构建命令、CT/GA 与 Combined 的独立提交、六件 index 文件族检查、CPU 预算及事务恢复，
统一见 [实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/doctor.sh)。
运行只依赖所选 alignment mode 的 index；线程数不构成 readiness 或重建 gate。
第三方依据：[Bowtie2 文件族选择](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)、
[Bismark genome preparation](https://docs.rs/bismark/latest/bismark/genome_prep/index.html)；以安装版本行为为准。

## 执行与缓存契约

执行器、科学身份、发布/恢复和缓存提交见[契约索引](#固定策略与资源)；数值口径见[生物学分母表](#生物学分母与零值解释)及[QC字段](#下游-qc-字段)。固定计算参数为：

| 计算 | value_type / summary_type |
| --- | --- |
| Gini分箱coverage | `total` / `sum`；流式float32累加与坐标约束见[实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/doctor.sh) |
| TSS fraction与count | 分别为`fraction` / `mean`、`target` / `sum`；两矩阵共用按sealed base适配的区间 |

Final QC 直方图和 CpG density 的展示参数、文件名及操作见[下游绘图说明](../downstream/README.md#指标与筛选)；样本集合与拟合边界见[实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/downstream/downstream_qc.py)。展示不改变QC计算或筛选。

## 下游重算与 TSS 抽检控制

实际重建 Raw 或统计前会撤销旧 QC 完成身份。失败后请重新运行 QC Processor，成功后再加载 Notebook QC 或生成 SingleCpG；这些入口拒绝缺失、过期或计算期间变化的完成身份。 Single 缓存命中且 QC/HQ/注释未变时，HDF5 与其完成标记均保持不变；比较时 NaN/null 均表示缺失，真实零值仍与缺失区分。

这些参数属于 Notebook / downstream Processor，不写入主 Pipeline config，不改变上游科学参数。RawAdata原生导入、CpG流式统计的完整校验及缓存提交边界见[实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/downstream/downstream_qc.py)，不在Raw重建前重复全表预扫。

| Notebook 参数 | Processor CLI | 默认行为 |
|---|---|---|
| `RECOMPUTE_RAW_ADATA` | `--recompute-raw-adata` | False；按当前 sealed 输入、reference、算法及矩阵身份复用 RawAdata 和空表集合，命中不解压/建链接；True 或身份失配时从 sealed CpG 重新导入，完整校验按上述契约执行 |
| `RECOMPUTE_CPG_STATS` | `--recompute-cpg-stats` | False；联合复用 Mapping_Stats.csv 与 CpG_Signal_Composition.csv，任一缺失/损坏或True时单次扫描重建两份统计，Raw/Gini/TSS独立复用 |
| `RECOMPUTE_GINI` | `--recompute-gini` | False；按身份和 bin size 复用 Gini |
| `RECOMPUTE_TSS` | `--recompute-tss` | False；按身份和 TSS BED 复用 TSS profile，SKIP_TSS 时不执行 |
| `TSS_SAMPLE_FRACTION` | `--tss-sample-fraction` | 0.10；有限 `(0,1]`，每个 CloneID 内按固定种子/Sample_ID 哈希排序抽检；少于 10 个全检，否则 `min(100,ceil(N*fraction))`；1.0 取消抽检、全部全量 |
| `TSS_SAMPLE_SEED` | `--tss-sample-seed` | 0；非负整数，行顺序不影响抽样 |

抽样保留原始CloneID身份（含前导零及合法NA/nan/None），包含空CpG，不按覆盖筛选。
实际方法、比例、种子、总数、各克隆数量及选中ID写入 `.cache/DNAme_QC_Input_Identity.json` 的 `tss_sampling`；CSV只保抽中且有覆盖的bins。
比例/种子/抽样算法变化仅失效TSS，其余QC与矩阵仍保全部细胞。子集复制、覆盖预筛与句柄清理见[实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/downstream/downstream_qc.py)，零值与缺覆盖见[分母表](#生物学分母与零值解释)。


### Notebook 末尾 single-CpG 矩阵

固定参数为`value_type="fraction"`、`summary_type="mean"`、`chunk_size=500`、`inplace=False`。
`RECOMPUTE_SINGLECPG=False`默认复用已验证矩阵，来源未知或身份变化时显式设为True；无Processor CLI选项。
Raw算法身份不匹配时须先重跑QC Processor，再按新generation重建SingleCpG，不能用此开关绕过。
坐标、原生reference顺序及缓存提交统一见[实现与边界](https://github.com/Felix-owo/Alopex/blob/v16/core/doctor.sh)，操作见[下游说明](../downstream/README.md#保存供下游复用的-single-cpg-矩阵)。
