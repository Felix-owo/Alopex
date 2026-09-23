# Alopex v13.13

Alopex 是面向 Cabernet、SRD、Droplet（DD-MET5）和 Cabernet–TAPS+ 单细胞 DNA 甲基化数据的 Snakemake Pipeline。它将原始 paired FASTQ 转为单细胞 CpG 结果、质量报告和可追溯的交付清单：

- paired-end FASTQ 样本识别与单细胞 demultiplex；
- Cabernet adapter trimming（R1/R2 Tn5 技术结构）；
- Cabernet/SRD 使用 BISCUIT 或 Bismark；Droplet 使用 Bismark non-directional 和 UMI-tools；TAPS 使用 BWA-MEM + Rastair；
- Cabernet/SRD 按位置去重；Droplet 按物理 R1 与 UMI 去重；TAPS 标记重复并在甲基化调用时排除；
- Cabernet/SRD/Droplet high-CpH / non-conversion 评估与过滤；TAPS 标为不适用；
- 直接生成 SnapATAC2 `pp.import_values` 所需的 0-based CpG 四列表，并生成兼容 `peak_file` 的 single-CpG BED3；
- BISCUIT 路线的可选 SNP 输出；
- cell QC、MultiQC 汇总与运行结果清单。

三个 backend 共用固定 Cutadapt `--error-rate 0.25`，R1/R2 overlap 为 3/12，
R2 同时剪除 9 bp 技术序列。各参数的含义与默认值见初始化项目生成的 `00_config/config.yaml` 模板行内注释；
参数依据与生物学表现验证所属的内部 benchmark 归档不随公开仓库分发。

## 关于本仓库

本仓库是 Alopex 的 GPL-3.0 公开发布副本，包含完整可运行的 Pipeline：核心 workflow、
Rust demultiplexer、下游 QC、Doctor 端到端回归 fixture（真实 Cabernet 抽样 reads，
数据集命名与来源身份已中性化）以及默认 Barcode Map 模板。以下内容不随公开副本分发：

- 维护规范（AGENTS.md）、参数表、流程图、测试套件与内部 benchmark 归档。

在本副本上 `core/doctor.sh` 的九个阶段均可完整执行。

普通用户只需要使用两个入口：

```text
core/doctor.sh        安装、构建 reference 并验证 Pipeline
core/run_pipeline.sh  创建项目、检查输入并运行分析
```

## 全流程

```text
从 GitHub 下载 Pipeline
  → 准备 4 个 FASTA + 2 个 GTF
  → 运行 Doctor，看到 FINAL STATUS: READY
  → 初始化项目
  → 检查 config.yaml（孔板协议另需 Barcode_Map.csv）
  → 将 paired FASTQ 放入 01_raw/
  → 刷新并检查 sample_manifest.tsv
  → dry-run
  → 正式运行
  → 检查 03_results/ 和最终 dry-run
```

## 1. 从 GitHub 下载 Pipeline

获取代码：

```bash
git clone --branch v13.13 https://github.com/Felix-owo/Alopex.git
cd Alopex
```

开始前需要：

- macOS 或 Linux；Linux/HPC 正式运行需要 Slurm；
- 一个可以从终端调用的 Conda，例如 Miniforge；
- `git`、`curl` 和 `gzip`；
- HPC 上可以使用 `sbatch` 和 `squeue`。

可先确认 Conda：

```bash
conda --version
```

不要跨机器复制已经构建的 `conda/`（含 Rust 构建状态与 demux 二进制）或 reference indexes；在目标机器上让 Doctor 构建。

## 2. 准备标准 reference

Pipeline 自动构建的标准 reference 只包括 `hg38` 和 `mm10`。Pipeline 不会自动下载原始 reference；首次运行 Doctor 前，需要准备以下六个文件：

```text
resources/reference_source/
├── hg38.fa
├── mm10.fa
├── lambda.fa
├── pUC19.fa
├── hg38.genes.gtf
└── mm10.genes.gtf
```

推荐来源：

| 固定文件名 | 来源 | 下载链接 |
| --- | --- | --- |
| `hg38.fa` | UCSC hg38 | [hg38.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) |
| `mm10.fa` | UCSC mm10 | [mm10.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz) |
| `lambda.fa` | NCBI `J02459.1` | [J02459.1 FASTA](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=J02459.1&rettype=fasta&retmode=text) |
| `pUC19.fa` | NCBI `L09137.2` | [L09137.2 FASTA](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=L09137.2&rettype=fasta&retmode=text) |
| `hg38.genes.gtf` | GENCODE v48 / GRCh38 | [GENCODE v48 GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz) |
| `mm10.genes.gtf` | GENCODE M25 / GRCm38 | [GENCODE M25 GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz) |

可以直接在 Pipeline 根目录执行：

```bash
mkdir -p resources/reference_source

curl -fL https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz \
  | gzip -dc > resources/reference_source/hg38.fa

curl -fL https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz \
  | gzip -dc > resources/reference_source/mm10.fa

curl -fL \
  'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=J02459.1&rettype=fasta&retmode=text' \
  > resources/reference_source/lambda.fa

curl -fL \
  'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=L09137.2&rettype=fasta&retmode=text' \
  > resources/reference_source/pUC19.fa

curl -fL \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz \
  | gzip -dc > resources/reference_source/hg38.genes.gtf

curl -fL \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz \
  | gzip -dc > resources/reference_source/mm10.genes.gtf
```

Doctor 会为两个物种构建主染色体加 lambda/pUC19 的 reference、BISCUIT/Bismark/BWA indexes，以及下游使用的 TSS 和 single-CpG reference。首次构建会根据可用 CPU/RAM 自动选择安全的并行度；高内存主机可同时构建两套 BISCUIT index。标准 Bismark CT/GA index 使用专用 `bismark_genome_preparation --parallel N`，使两个 `bowtie2-build` 各获得 N 个线程；16-core 预算时默认 N=8。所有 Bismark 构建显式绑定当前环境的 indexer。Combined 独立校验和提交，构建线程参数最多为 8；当前工具还会在候选目录生成 CT/GA，但修复 Combined 时保留已发布的 CT/GA。

## 3. 构建并验证 Pipeline

在 Pipeline 根目录运行：

```bash
bash core/doctor.sh
```

网络受限或希望把下载与构建分开时，可先在登录节点只预热下载缓存（Conda 包与 Rust 依赖，
不创建环境、不编译、不提交作业），随后各构建段离线复用：

```bash
bash core/doctor.sh download
```

Doctor 统一准备 Conda 环境（含 Rastair 2.2.0）、Rust demultiplexer、hg38/mm10 reference 与三后端索引（含普通 BWA），并执行三条隔离的端到端测试路线（A/B 并行，随后执行独立 TAPS 真值路线）。输出为连续 9 个阶段：平台、Conda、Rust、reference、原始输入完整性、A/B 项目配置、A/B 并行运行、A/B 结果验证、C 项目配置与运行及真值验证。每个阶段只有一组 START/OK，失败显示对应阶段 FAIL；检查后的修复以阶段内 REPAIR 提示。

共享文件系统上的 Conda 冷启动可能较慢。CLI 探针在 GNU timeout 可用时允许 300 秒；超时或命令失败保留原始退出码，只有 help 成功返回后才判断选项是否受支持。

标题、三条路线和 detail/summary 路径仅在启动时显示一次，概要实时写入 summary，结尾显示最终状态与沙箱去留。manifest 刷新、构建命令、退出码及原生诊断写入 detail；`logs/doctor/latest.detail.log` 和 `latest.summary.log` 始终指向最近一次运行。失败保留详细日志与沙箱，终端不重复播放整份日志。

Doctor 自身的诊断说明使用中文，保留 START/OK/FAIL、READY 和原生工具输出以便定位错误。

`logs/` 保存日志、实验与审计证据以及 Doctor 沙箱，不是生产源码目录；其中的历史实验脚本和源码副本不参与 Pipeline 调用或当前代码审计。新增的一次性维护脚本和备份使用系统临时目录。

HPC 上请从 login 节点运行构建。Doctor 会显示 Slurm 队列状态和日志尾部，并在作业结束后校验产物。固定的 hg38/Bismark 与 mm10/BISCUIT 测试各使用一个 2-CPU/20-min worker，内存分别为 16 GiB 与 24 GiB；若已在 allocation 内，前两条测试各用 2 cores 本地执行；TAPS 随后使用 2 CPU / 24 GiB / 20 min worker。`conda --check` 可在 allocation 内运行，环境构建与发布须在 login 节点进行。

`run_pipeline.sh --dry-run` 和正式运行均可在已有 Slurm allocation 内启动；launcher 加载 Doctor 公共函数时不会执行 Conda 构建节点检查。该限制仅在 `doctor.sh conda` 解析参数后、开始构建前生效。

首次构建 reference 需要较多时间和磁盘空间；在 Linux/HPC 上应从能够提交 Slurm 作业的 login/controller shell 启动。失败或按 Ctrl-C 中断时，构建候选会清理，失败测试沙箱与诊断日志保留。reference 提交中断时由持锁 builder 回滚；若回滚也失败，journal 与恢复所需备份保留到下次修复。原始 reference、已经验证的组件和可复用下载缓存会保留。

Doctor 复用未变化且可用的环境与 reference；每次调用都重新执行三条完整测试路线。READY 表示本次检查和回归通过。

| 路线 | 小型 FASTQ 输入 | 检验重点 |
| --- | --- | --- |
| hg38 / Cabernet / Bismark | hg38_cabernet 真实 raw 1,050 对 + 40 对合成 DNA | 真实人源比对、CpG、重复及 non-conversion 过滤 |
| mm10 / SRD / BISCUIT | mm10_srd 真实 Cabernet raw 1,050 对 + 40 对合成 DNA + 双管各 24 对合成 RNA | 真实小鼠比对、CpG、high-CpH 过滤、SRD 双管 RNA 保留 |
| hg38 / TAPS / Rastair | 独立合成 133 对 | 相反转化信号、双链/重叠/重复、对照真值及 BAM 保留 |

真实抽样只涉及 A1/B1 两个分析 barcode，另有低读数 C1、零输出 D1 和天然错配/未分配 reads；两个原始样本压缩后共 413,636 bytes（约 404 KiB），来自内部基准项目（数据集 hg38_cabernet/mm10_srd，来源身份已在 fixture manifest 中中性化）。原始 read 名称、序列和质量均未改写；抽样序号、源文件身份与 SHA256 固定在 `resources/doctor_benchmark/benchmark.json`，运行时无需连接 HPC。合成 reads 单独命名和记账，共三项目、五个分析 cell、2,361 对 reads。

A/C 同时核验外部 RNA 关联随发布保留且不触发本地 RNA 拆分；B 核验 SRD 本地双管输入。

成功报告另保留 `logs/doctor/doctor.<run>.<route>.json`，包含来源、环境/reference/source snapshot、逐 cell 漏斗、CpG 计数和交付 inventory；沙箱清理后仍可审计。SRD 的 RNA 布局和 TAPS 化学部分目前是合成验证，不能据此宣称真实 SRD/TAPS 文库已通过生产验证。混合输入的对照计数也不代表真实实验的转化效率。

最终显示以下状态表示本次工具链、真实 Cabernet 抽样及合成边界回归通过：

```text
FINAL STATUS: READY
```

如果 Doctor 失败，修正它报告的缺失文件或环境问题后，重新运行同一个命令即可。

## 4. 创建分析项目

Pipeline 安装目录与数据项目目录应分开。创建项目：

```bash
mkdir -p /path/to/Patient001
cd /path/to/Patient001
/path/to/Alopex/core/run_pipeline.sh --init-project
```

初始化后项目结构为：

```text
Patient001/
├── 00_config/
│   ├── config.yaml
│   ├── Barcode_Map.csv        # 仅孔板协议
│   └── sample_manifest.tsv
├── 01_raw/
├── 02_work/
├── 03_results/
├── 04_logs/
├── 05_tmp/
├── 06_downstream/              下游 QC 首次运行时按 delivery_id 创建
└── run_pipeline.sh
```

后续命令都从项目目录运行：

```bash
cd /path/to/Patient001
```

## 5. 检查项目配置

运行前至少检查：

1. `00_config/config.yaml`；
2. 孔板协议的 `00_config/Barcode_Map.csv`（Droplet 跳过）；
3. FASTQ 导入后生成的 `00_config/sample_manifest.tsv`。

### config.yaml

模板显式列出全部参数及其默认值（`high_cph` 在 Cabernet/SRD/Droplet 中固定为「统计并剔除」，TAPS 不适用，无开关；
模板列出其唯一参数 `excluded_contigs` 及默认对照列表），
通常只需要修改物种、实验协议和甲基化后端：

```yaml
schema_version: 2
species: hg38

analysis:
  protocol: cabernet           # cabernet / srd / taps / droplet
  methylation_backend: bismark # cabernet/srd: biscuit/bismark；droplet: bismark；taps: rastair
```

配置项说明：

| 配置 | 什么时候需要修改 |
| --- | --- |
| `species` | 选择当前项目物种。标准值为 `hg38` 或 `mm10`。 |
| `references` | hg38/mm10 通常保留默认值；custom species 需要加入物种名及最终 FASTA 路径。 |
| `analysis.protocol` | 按实验选择 `cabernet`、`srd`、`droplet` 或 `taps`。 |
| `analysis.methylation_backend` | Cabernet/SRD 选择 `biscuit` 或 `bismark`；Droplet 只用 `bismark`；TAPS 只用 `rastair`。 |
| `demux.min_matched_read_pairs` | 一个 barcode bucket 被保留所需的最少 matched read pairs；默认 10。 |
| `demux.dna_w_spacer_len` | DNA barcode 后 spacer 长度；只有建库设计不同时修改，默认 0。 |
| `high_cph.excluded_contigs` | 不参与 high-CpH/non-conversion 判定的 control/mitochondrial contig；默认 `pUC19/lambda/chrM`。 |
| `biscuit.library_type` | `directional` 或 `non_directional`；Pipeline 会转换成 BISCUIT `align -b` 参数。 |
| `biscuit.high_cph_retention_threshold` | BISCUIT high-CpH read 判定阈值；必须是 `[0,1]` 内的有限数值，默认 0.7。 |
| `biscuit.generate_snp` | BISCUIT-only SNP 输出功能开关；默认 `false`。 |
| `bismark.library_type` | `directional`、`non_directional` 或 `pbat`。 |
| `bismark.local_alignment` | `false`（默认）使用 combined-index end-to-end；`true` 使用 faithful CT/GA + local alignment。 |
| `bismark.non_conversion_percentage_cutoff` | Bismark non-conversion 百分比阈值；默认 70。 |
| `bismark.non_conversion_minimum_count` | 判定 percentage filter 所需的最少 informative non-CG count；默认 5。 |
| `runtime.local_executor_cores` | macOS local executor 的总 core 预算，默认 32；Linux/HPC 使用 Slurm executor。 |
| `runtime.latency_wait_seconds` | 输出在本地或共享文件系统变得可见的等待秒数；默认 120。 |
| `slurm.jobs` | Linux/HPC 上允许同时在途的最大作业数；默认 30。 |
| `slurm.controller_cores` | Linux login/controller 上供 Snakemake 本地控制任务使用的 core 数；默认 2。 |
| `slurm.account/partition` | 只有集群明确要求时填写，否则保留 `null`。 |
| `slurm.qos` | 提交作业的 QoS；默认 `huge`，设为 `null` 时沿用集群默认 QoS。 |
| `slurm.node_tmpdir` | 没有 `$SLURM_TMPDIR` 时使用的绝对临时目录；默认 `/tmp`，不可用时可设为 `null`。 |
| `retention.keep_final_bam` | 默认 `false`；`true` 保留 `02_work/` 最终 BAM。TAPS 保留全部记录的 marked BAM+BAI 供 SNP/SNV 重分析，其余后端保留过滤后 DNA BAM；磁盘代价为全部 cell 的最终 BAM。 |

> 参数解释与默认值见 `00_config/config.yaml` 模板行内注释；内部固定的参数不放入 `config.yaml`。

> Bismark alignment mode 由 `bismark.local_alignment` 选择，只有当前模式使用的 index 是 runtime 依赖；过滤后的 paired-end BAM 保持 read-pair 相邻直达 methylation extraction，避免重复排序 I/O。extraction 所有档位固定申请 8 CPU + 16 GiB，attempt 2/3 内存分别为 24/32 GiB。

重要组合限制：

- custom species 需自行准备 FASTA、FAI 与所选 backend/alignment mode 的 indexes，再填写 `references.<species>`。

### Cabernet–TAPS+ 配置

```yaml
analysis:
  protocol: taps
  methylation_backend: rastair
retention:
  keep_final_bam: true
```

此路线用于 Watchmaker **7BK0003-024 TAPS+**：普通 C 保留，5mC+5hmC 转为 T。barcode、接头和 gap filling 使用普通 C，保留 Cabernet 的 reads 布局。流程为 demux → Cutadapt → BWA-MEM → fixmate/排序/markdup → Rastair 2.2.0 → canonical CpG。

每个 cell 的 `02_work/rastair/<Sample_ID>.marked.bam` 与 `.bai` 保留测得序列、全部比对记录和重复标记，final manifest 记录路径与用途；它们不进入 `03_results` inventory。甲基化调用使用 MAPQ≥20、baseQ≥30，排除重复/不合格配对/非主比对，重叠 mate 只计一次。参考 CpG 任一端调用为变异时，两端均排除。后续 SNP/SNV 应使用支持 TAPS 化学的方法，不能把所有 C→T/G→A 当作变异。

TAPS 的 high-CpH/cDNA 剔除不适用，下游两项 CpH 指标为缺失值。MultiQC 显示 lambda CpG 假阳性率、pUC19 CpG 转化效率及观测分母；使用实物对照前应确认 lambda 未甲基化、pUC19 为 CpG 甲基化。首版没有真实 TAPS FASTQ 验证。[Watchmaker 原理](https://www.watchmakergenomics.com/taps.html)

### Barcode_Map.csv

孔板协议初始化时复制默认 barcode map；已配置 Droplet 的项目跳过此步骤。若实验使用不同 barcode、plate 名称或 cell 顺序，必须在运行前替换或修改 `00_config/Barcode_Map.csv`。

- Cabernet/TAPS 需要 `DNA_Barcode`、`PlateID`、`Cell_Order`；
- SRD 还需要 `RNA_Barcode`；
- barcode 只能包含 `A/C/G/T`，长度必须为 8 或 10 bp；
- barcode 之间必须能够进行唯一的单碱基纠错，即同长度 barcode 的 Hamming distance 必须大于 2；
- `PlateID` 和 `Cell_Order` 必须唯一。
- 表头去 BOM、外围空白并忽略大小写后不得重复或为空；每行列数必须与表头一致。`PlateID` 只接受字母、数字、点、下划线、连字符，拒绝单独的 `.`/`..`。

## 6. 准备 raw FASTQ 和 sample manifest

将 paired FASTQ 放入项目的 `01_raw/`。支持 `.fastq.gz` 和 `.fq.gz`，并支持以下命名：

```text
SampleA.R1.raw.fastq.gz       SampleA.R2.raw.fastq.gz
SampleA_L001_R1_001.fastq.gz SampleA_L001_R2_001.fastq.gz
SampleA_L02_R1.fastq.gz      SampleA_L02_R2.fastq.gz
SampleA_R1_001.fastq.gz      SampleA_R2_001.fastq.gz
SampleA_R1.fastq.gz          SampleA_R2.fastq.gz
SampleA_L1_1.fq.gz           SampleA_L1_2.fq.gz
SampleA_1.fq.gz              SampleA_2.fq.gz
SampleA.1.fq.gz              SampleA.2.fq.gz
SampleA.R1.fq.gz             SampleA.R2.fq.gz
样本_SampleA_L1_1.fq.gz      样本_SampleA_L1_2.fq.gz
```

同一样本可以包含多个 lane/chunk（lane 号 1-3 位且为 1-999，3 位 chunk 段可省略）。mate 可为 `R1/R2` 或 `1/2`，分隔符支持下划线、点和连字符；lane 按数值配对排序，`L1/L01/L001` 表示同一个 lane，缺省 chunk 等于 `001`。每个 R1 segment 必须有对应 R2，同一 segment 的重复 mate 仍会报错。

识别时保留原文件名，只去掉内部样本名开头的 `样本_`；例如 `样本_PT_Mon_1_FKDL123_L1_1.fq.gz` 对应 raw sample `PT_Mon_1_FKDL123`，`Mon_1` 和批次编号不会被删除。其余内部样本名只接受英文字母、数字、点、下划线和连字符，不能是单独的 `.` 或 `..`。有前缀与无前缀的两种基名若映射到同一 raw sample，即使 lane 不同也拒绝隐式合并。manifest、snapshot 和 workflow 使用同一识别结果，无需手工改名或增加 config 参数。

同一 raw sample 的单文件式命名对（`SampleA_R1/R2.fastq.gz`）如果分散在 `01_raw/` 下的多个子文件夹中，`--refresh-manifest` 会按文件夹路径排序自动改名为 chunk 式（`SampleA_R1_001/R2_001`、`_002/_002`…）并合并为同一个样本。若两份同 mate 文件字节级相同（疑似重复副本）、缺 R1/R2 配对，或该样本下还混有 lane/chunk 命名文件，则拒绝处理并提示手动清理。候选 manifest 在改名前完成校验，应用后验证或写出失败会回滚改名；内部 `dna_pipeline refresh-manifest --dry-run` 按相同计划预览，不写临时文件。

导入 FASTQ 后运行：

```bash
./run_pipeline.sh --refresh-manifest
./run_pipeline.sh --dry-run
```

只分析部分样本时，在刷新后编辑 `00_config/sample_manifest.tsv`，仅保留需要分析的行，
然后直接执行 `--dry-run` 或无参正式运行。匹配使用 `dna_raw_sample`；SRD 还读取本地
`rna_sample`，Cabernet/TAPS 的外部 RNA 关联不参与 DNA 选样。未登记 FASTQ 可继续放在
`01_raw/`，会提示并跳过，不进入 DAG、snapshot 或结果；它们的改动不会触发新 snapshot。
已登记输入缺失、缺配对或重名冲突仍报错。**选样后不要再执行 `--refresh-manifest`，
除非希望重新扫描并补登记所有有效输入。**

该命令会根据 `01_raw/` 更新 `00_config/sample_manifest.tsv`，然后检查配置、FASTQ、孔板协议的 barcode map、reference 和待执行任务，不会正式分析数据。若现有 manifest 仅有 CRLF/LF、末尾换行、行顺序、字段外围空白或 protocol 大小写差异，刷新会保持文件字节、mtime 和备份不变；只有样本映射、物种、协议或 notes 等行语义真实变化时才更新滚动 `.bak` 备份并写回，避免无意义刷新触发 demux 及下游全量重算。

`sample_manifest.tsv` 的 `rna_sample` 按协议解释：

| 协议 | `dna_raw_sample` | `rna_sample` |
|---|---|---|
| Cabernet / TAPS | 本项目 DNA FASTQ basename，必填 | 可选，外部 RNA_DARLIN 项目的 RNA FASTQ basename |
| SRD | 本项目 DNA tube basename | 本项目 RNA-enrichment tube basename；两个 tube 不能仅从文件名推断 |

Cabernet/TAPS 可将配对的 DNA/RNA 放在同一行，例如：

```text
sample_id  dna_raw_sample  rna_sample  species  protocol  notes
DNA_alias  RAW_DNA         RAW_RNA     mm10     taps
```

对应 RNA_DARLIN manifest 可使用 `sample_id=RNA_alias, rna_sample=RAW_RNA`。
两边项目样本名可以不同；`rna_sample` 填原始 RNA basename，不填 `RNA_alias` 或文件路径。
RNA FASTQ 放在 RNA 项目；DNA 项目只需要 `RAW_DNA` FASTQ，也无需为这项关联增加 RNA barcode。
一个外部 RNA raw 只能关联一个 DNA 样本；每个本地 FASTQ 也只能归属一个物理输入。
刷新 manifest 在本地 DNA FASTQ 仍存在时保留外部关联；DNA FASTQ 消失时会移除整行，并提示随行移除的 RNA 关联。疑似误放入 DNA 项目的 RNA FASTQ 不会自动注册为 DNA。

DNA 发布后，RNA 下游从当前已发布 manifest 读取关联，结合最终 `RNA.h5ad` 的 raw/library
身份和 `PlateID` 对齐细胞，输出 `Integrated_QC_Information.csv`。它校验两边参与细胞集合、
条码和身份，不静默丢掉未匹配细胞；DNA 20 列与 RNA 13 列 QC 表保持不变。
SRD 的 `rna_sample` 表示本地输入，跨项目整合使用共享的项目样本名与 `PlateID`。
关联修改后重新运行 Alopex 与 downstream，随后在 RNA notebook 填 `DNA_PROJECT_ROOT`。
修改 manifest 后再次检查：

```bash
./run_pipeline.sh --dry-run
```

## 7. 运行 Pipeline

macOS 和 Linux/HPC 使用相同的用户命令：

```bash
./run_pipeline.sh --dry-run
./run_pipeline.sh
```

macOS 使用 local executor。Linux/HPC 的当前 shell 承载 Snakemake controller，规则作业通过 Slurm 提交；已处于 compute allocation 时直接复用该节点。交互终端每 60 秒显示进度，输出重定向时每 300 秒显示；完整日志位于 `04_logs/controller/<run_id>.log`，阶段与 Slurm 日志位于 `04_logs/rules/` 和 `04_logs/slurm/`。

### 作业资源与时限

Cabernet/Bismark 和 Droplet/Bismark 的每个 cell 独立组成一个 Slurm 作业，依次执行比对去重（含 Cutadapt）、过滤、extraction 和 CpG 转换。全部 S/M/L/XL 档位均分组，SRD 与 BISCUIT 按规则独立调度。

| Alignment mode | S/M/L/XL CPU | 首次内存 | 第一次重试 | 第二次重试 |
| --- | --- | --- | --- | --- |
| combined-index end-to-end（默认） | 8/12/12/16 | 均为 16 GiB | 均为 24 GiB | 均为 32 GiB |
| faithful CT/GA local | 12/12/16/16 | 24/24/32/32 GiB | 36/36/48/48 GiB | 48/48/64/64 GiB |

**完整 cell 作业每次申请 10 小时，四个阶段共享这段时间。** 每阶段的 150 分钟仅用于 Snakemake 累加组资源，阶段没有独立 timeout；alignment 可以使用超过 150 分钟。部分阶段重跑申请不超过 10 小时。最多尝试三次，重试只增加内存；排队和多次尝试的累计耗时不计入单次时限。

### 失败恢复

修正问题后重新执行 `./run_pipeline.sh`。Snakemake 根据依赖与有效输出决定重算范围；失败作业的成员输出会清理，Cabernet/Bismark 组可能需要整组重跑。调用日志初始化 helper 的资源规则每次尝试均建立 `.attempt1/.attempt2/.attempt3` 对应日志硬链接，自动重试清理主日志后仍保留本轮各次证据；新一轮首次尝试会清理上一轮陈旧归档。共享文件系统的时钟偏差不作为作业失败依据，命令退出码、输出完整性与科学守恒仍会检查。

同一项目从 workflow 执行到发布均由 Snakemake 项目锁保护。遇到残留锁时，先确认没有运行中的 launcher，再删除项目的 `.snakemake/locks` 后重试。保留 `02_work/` 与 `.snakemake/` 可支持恢复；各规则会回收自身独占的临时 scratch。scratch 不允许 `.`/`..` 路径组件或符号链接父目录。

Snakemake 的 `onsuccess` 在释放项目锁前将结果以硬链接逐文件原子安装到 `03_results/`，完整验证后最后写入 `run_manifest.json`。`02_work/results_stage/` 与 `03_results/` 必须位于同一文件系统，公开结果的父目录必须是真实目录；链接失败不回退复制。发布中断时完成清单缺失，重新运行 launcher 可从稳定结果阶段完成发布。开启 `retention.keep_final_bam` 后，缺失的最终 BAM/索引会进入 DAG 恢复，不能由“交付已完成”提前退出。

发布清单只收录本次最终 sample manifest 声明的 CpG、SNP（含 index）和固定 QC 文件。移除样本或关闭 SNP 后，stage 中的旧文件保留供 Snakemake 复用，`03_results/` 中的对应旧文件会清理。
若全部 cell 都低于 demux 保留阈值，仍发布含 demux 统计的 MultiQC 和最终 manifest，便于检查筛除原因。

## 8. 数据产出与完成确认

主要结果位于：

```text
03_results/
├── CpG/                     每个细胞唯一一份、SnapATAC2-ready 的 0-based 四列 CpG 表
├── SNP/                     可选，BISCUIT 路线且 biscuit.generate_snp=true；BED.gz + tabix index
├── QC_Results/
│   ├── multiqc_report.html
│   └── sample_manifest.tsv  最终样本/细胞清单，内含下游需要的 QC 指标和结果路径
└── run_manifest.json        sealed delivery 的唯一完成/文件校验清单

02_work/RNA/                 仅 SRD；RNA 部分原始数据保留区（不进 03_results 发布）
├── FASTQ/<route>/           demux RNA FASTQ 的硬链接收集
└── BAM/                     被筛除 BAM（<cell>.rna.bam + index），恒定产出
```

CpG 文件固定为 SnapATAC2 `pp.import_values` 可直接读取的四列：`chrom pos methyl unmethyl`，其中 `pos` 为 0-based。coverage 在需要时由 `methyl + unmethyl` 即时计算。SRD 的 RNA 原始数据（FASTQ + 被筛除 BAM）保留在 `02_work/RNA/`：FASTQ 是 demux 输出的硬链接（发布后 demux scratch 回收仍存续），被筛除 BAM 来自 DNA tube 的 high-CpH/non-conversion 过滤，随 CpG 一起产出。DNA workflow 的 `02_work/` 与 `.snakemake/` 保留为 Snakemake 可恢复状态，不复制进 `03_results/`。

成功后的持久文件分成三类：

| 位置 | 是否必要 | 成功后的策略 |
|---|---|---|
| `03_results/` | 必要 | 保留唯一一代最终结果及 `run_manifest.json` |
| `04_logs/controller/<run_id>.log` | 必要 | 按 run id 保留 controller 审计日志 |
| `04_logs/provenance/` | 必要 | 保留内容寻址的 run snapshots 与共享小型输入副本；输入文件只读、目录对 owner 可写，可整体清理运行目录 |
| `06_downstream/<delivery_id>/QC_Results/` | 可选下游产物 | 与 sealed delivery 绑定；下游分析不得写回 `03_results/` |
| `02_work/`、`.snakemake/` | 增量恢复状态 | 成功/失败均保留，由 Snakemake 根据依赖决定复用或重算 |
| `05_tmp/` | rule scratch 根目录 | 临时 generation 由对应 rule 成功/失败清理；不得作为公开结果 |
| `04_logs/rules/`、`04_logs/slurm/` | 运行诊断 | 保留与当前运行有关的日志；不得因为其存在阻断下一次运行 |
| `00_config/sample_manifest.tsv.bak` | 可选恢复保护 | 最多保留一份滚动备份 |
| `conda/releases/<release_id>`、`logs/doctor.<id>`（失败沙箱） | 生成物缓存 | release 发布成功后自动回收更早代际（保留上一代作回退）；Doctor 失败沙箱保留 14 天后于下次启动清理。envs 未变化的宿主在 `conda/current` 复用同一环境，不会重复安装 |

正式运行成功时，controller log 会记录 `Published sealed delivery`。发布函数自动回收 `02_work/demux/`（先删 `state/` 提交边界，再删整棵 scratch）；state 删除失败则保留 FASTQ 并告警，下次启动可幂等补收。发布之前该目录始终完整保留以支持断点续跑。之后可执行 `./run_pipeline.sh --dry-run` 检查当前 DAG；需要重算哪些规则由 Snakemake 当前 dependency/rerun trigger 决定。

任何依赖 demux 输出的修改都会重新执行 demux checkpoint 并重新生成其 FASTQ scratch；该策略在发布后释放 demux 磁盘空间。重建时间戳本身不触发比对/CpG 重算；这些规则直接跟踪真实 raw 输入、孔板路线的 Barcode Map、样本声明、demux binary 和分析参数。不要手动只删除其中 FASTQ 而保留 completion。也不要单独删除 `02_work/results_stage/`：它与公开结果共享 hardlink 内容（不额外占用磁盘），同时也是 Snakemake 增量重算的稳定 biological output stage 与重发布来源。

MultiQC 的 high-CpH 展示 JSON 若已清理，报告规则会从保留的过滤摘要重建；这一展示恢复不会要求重新生成 BAM 或重新比对。

通过 demux 阈值的 cell 即使全未比对或 trim 后全空，也保留合法空 CpG 和 QC，不阻断整批交付；零 primary reads 的 mapping 为 0%，无 CpH 观测为缺失值。

最终至少应存在以下非空文件：

```text
03_results/QC_Results/multiqc_report.html
03_results/QC_Results/sample_manifest.tsv
03_results/run_manifest.json
```

交付清单的 `pipeline.name` 固定为 `Alopex`；`delivery-ready` 由内部 CLI 生成，
正式发布由 Snakemake 持锁的 `onsuccess` 调用 `publish_delivery` 完成。

## 常见问题

- Doctor 报 reference 缺失：确认六个 source 文件的固定路径和文件名，然后重新运行 `bash core/doctor.sh`。
- 找不到 paired FASTQ：检查文件名是否符合支持格式，并确认每个 R1 都有对应 R2。
- 配置校验失败：根据错误修改 `00_config/config.yaml`；只使用参数文档与 config schema 支持的字段。
- Bismark 配置失败：检查 `bismark.library_type` 与 reference/index。
- Slurm 作业长时间 PENDING：使用 `squeue -u "$USER"` 查看状态；资源或 partition 问题请联系集群管理员。
- Slurm 作业反复在已确认故障的节点启动失败：可从项目目录执行 `SBATCH_EXCLUDE=节点名 ./run_pipeline.sh` 排除该节点并断点续跑。launcher 保留此原生 Slurm 变量并显示排除列表；其它外部调度覆盖仍会清除。无需刷新 manifest 或删除已有结果。
- 分析失败：先查看 `04_logs/` 中对应 rule 的日志，修正问题后重新运行 `./run_pipeline.sh`。

## 9. 下游 QC 与可视化

使用 Doctor 管理的 notebook 环境运行 [下游分析](downstream/README.md)，结果写入项目的 `06_downstream/<delivery_id>/`。Notebook 首格接受 Pipeline、项目和 TSS reference 的绝对路径。第二格可按需强制重算；`RECOMPUTE_RAW_ADATA=True` 从 sealed CpG 重建 AnnData，默认按缓存身份自动复用。

Notebook 最后一格生成供 MethylTree 等分析复用的 `QC_Results/SingleCpG_Adata.h5ad`，按细胞 ID 保存完整 QC 和 `HQ` 注释，生成后关闭文件。重复运行复用已验证的矩阵，注释变化只更新 `obs`；需要替换过期或未验证的文件时使用 `RECOMPUTE_SINGLECPG=True`。具体参数和读取方式见[下游说明](downstream/README.md#保存供下游复用的-single-cpg-矩阵)。

`DNAme_QC_Information.csv` 固定输出 20 列，包含样本身份、甲基化、signal composition、Gini、mapping 与 read-pair retention 指标；完整计数与 backend 诊断保留在主 Pipeline 的最终 manifest 和 MultiQC 中。CpG Density 从 QC 表、sealed manifest 与 composition 缓存读取所需字段。

CpG 统计、Gini 与 TSS 分别缓存，缓存身份包含相关输入、参数、计算代码与库版本。修改绘图或说明文字可复用科学统计；修改 Gini bin size 或 TSS BED 只使对应统计失效。

TSS 无覆盖时保留仅表头 CSV 并清理旧 TSS 图；显式跳过 TSS 时清理该交付的旧 TSS CSV/PDF。

## 许可证

Alopex 以 [GPL-3.0](LICENSE) 发布。本仓库为公开发布副本，与开发仓库按发布版本同步。

## Droplet DD-MET5

在初始化前的 `00_config/config.yaml` 中设置 `analysis.protocol: droplet`、
`analysis.methylation_backend: bismark` 和 `bismark.library_type: non_directional`，然后使用同一
`--init-project`、`--refresh-manifest` 和无参数 launcher。Droplet 无需 Barcode Map，manifest 的
`rna_sample` 留空。每个原始文库独立进行结构计数与谷底调用，输出身份为项目样本加 17 bp barcode，
`plate_id` 为空；下游仍消费相同的四列 CpG、20 列 QC 和 sealed manifest。

R1 的 CB17、UMI12、TSO13、Linker17、ME19、gap9 在拆分时移除，R2 去除前 9 bp。
两种转换模式分别匹配，歧义拒绝；UMI 用于物理 R1 位置与方向上的 directional 去重，忽略 R2 端点。
谷底调用后，结合单碱基近邻、差异位点质量和多个共享 UMI＋insert 筛查错误衍生条码。先保留候选条码的有限分子特征，再完整扫描高丰度近邻查找共享分子；证据不足时保留。复核、筛除、证据不足与采样饱和数量进入 MultiQC。
筛查后再与固定的官方 ME5 设计白名单取交集；原谷底保持不变。设计外候选的条码、计数和纠错去向留在 calling 诊断中，不标成空液滴。拆分时精确匹配优先；合法但未 called 的设计条码保持未分配，其余观测只做唯一 Hamming-1，歧义拒绝，不按丰度强行分配。DNA 不使用 RNA 空液滴名单。白名单为 [SeekSoulMethyl 官方 ME5/U3CB_methylation](https://github.com/seekgene/SeekSoulMethyl/blob/nf_rna_methy/dependence/seeksoultools/utils/barcode/ME5/U3CB_methylation.txt.gz)（829,440 个 AGT17），随仓库分发于 `resources/`，Python 运行时与 Rust 编译均按固定 SHA256 校验。
Droplet 拆分会在文件句柄限额允许时保留全部细胞的输出 writer，以减少反复开关小文件；容量不足时自动使用有界缓存。
逐条判定及初始/最终数量进入 calling 日志，移除数进入 MultiQC。该筛查不改变原谷底阈值。
谷底调用无双峰时明确失败；全部 Droplet 科学参数为固定契约，不暴露为 config 参数。
官方对照采用 [单细胞 BAM 去重说明](https://github.com/seekgene/SeekSoulMethyl/blob/nf_rna_methy/docs/How_to_deduplicate_single_cell_bam.md)；谷底调用是本 Pipeline 的指定策略，不声称与官方完整流程等价。
真实数据验证与回归证据保留在开发仓库，代码按发布版本同步。


### Droplet 的配套 RNA 独立关联

配套 DD-MET5 RNA 在 Lemmus（RNA_DARLIN_Pipeline）的 `pipeline_mode: droplet` 路线独立前处理、
simpleaf/piscem/alevin-fry 定量、emptyDrops 调用与 RNA QC；DNA 的 called cells、阈值、UMI 去重和
QC 不受 RNA 名单约束，DNA droplet manifest 的 `rna_sample` 继续为空。

两模态由 RNA 项目 `00_config/dna_rna_library_map.tsv`（`dna_project_sample_id`/`rna_sample` 两列）
关联：RNA 下游读取当前 complete Alopex sealed manifest 的 `dna_barcode` 与同代 DNA QC，按明确对应
文库＋完整 17 bp barcode 精确匹配 finalized RNA.h5ad，不用 PlateID、Cell_Order 或样本名猜测；修改
该表无需重新发布 DNA。RNA 下游输出 `DNA_RNA_cell_links.csv`（双方细胞全集与归属）与
`Integrated_QC_Information.csv`（仅共同且双方 QC 可用细胞），重复身份、多对一映射与错误交付代际
直接报错。完整关联契约、RNA 侧 QC 阈值规则与 Lemmus 参数见 Lemmus 文档与开发仓库维护规范；
本轮真实 RNA／DNA 子集验证的范围和结果保留在开发仓库。
