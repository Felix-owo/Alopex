# Alopex v14.10

Alopex 将单细胞 DNA 甲基化 paired FASTQ 转为单细胞 CpG、MultiQC 报告和可追溯的完整交付。

| 协议 | 后端 | 主要处理 |
|---|---|---|
| Cabernet | BISCUIT / Bismark | 孔板解复用、trimming、位置去重、high-CpH/non-conversion 过滤 |
| SRD | BISCUIT / Bismark | DNA/RNA tube 解复用、DNA 分析，并保留 RNA FASTQ 与被筛除 BAM |
| Droplet DD-MET5 | Bismark non-directional | DNA 独立谷底 calling、官方设计白名单、逐胞嘧啶 UMI 共识 |
| Cabernet–TAPS+ | BWA-MEM + Rastair | 标记重复、调用 CpG 和 SNP，统计 Non-CpG；不执行 CpH/cDNA 筛除 |

用户入口只有 `core/doctor.sh`（构建与维护）和 `core/run_pipeline.sh`（项目运行）。
完整配置及参数默认值见初始化项目生成的 `00_config/config.yaml`；下游用法见
[下游 README](downstream/README.md)，验证范围和升级说明见本版 Release。

本仓库是 Alopex 的 GPL-3.0 公开发布副本，包含核心 workflow、Rust demultiplexer、下游 QC、
Doctor 回归 fixture（真实 Cabernet 抽样与合成对照）和默认 Barcode Map。
维护规范、内部文档、测试套件与 benchmark 归档不随公开副本分发；Doctor 可完整运行。

## 1. 安装

支持 macOS 和 Linux；Linux 正式运行使用 Slurm，需要 `sbatch`、`squeue`、`sacct`。
先准备可调用的 Conda、git、curl 和 gzip：

```bash
git clone https://github.com/Felix-owo/Alopex.git
cd Alopex
conda --version
```

也可从源码归档部署，无需 `.git`。Conda、Rust 和 reference indexes 必须在目标机器构建，
不要跨平台复制环境。版本以 Git tag 为准。

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

按表中链接下载，解压为上述固定文件名。例如：

```bash
mkdir -p resources/reference_source
curl -fL https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz \
  | gzip -dc > resources/reference_source/hg38.fa
```

## 3. 构建与检查

在 Pipeline 根目录执行：

```bash
bash core/doctor.sh
```

Doctor 检查或修复环境、Rust demultiplexer、hg38/mm10 reference 和三个 backend 的 indexes，
再运行 A（hg38/Cabernet/Bismark）、B（mm10/SRD/BISCUIT）、C（TAPS/Rastair）三条隔离路线。
HPC 登录节点提交 Slurm workers；macOS 或已有 allocation 内串行运行。

成功结尾为 `FINAL STATUS: READY`。日志在 `logs/doctor/`；成功清理沙箱，失败保留沙箱和诊断。
Doctor 每次重跑三路线，READY 仅证明本次工具链及小样本回归通过；完整验证边界见本版 Release。

可单独维护某一段，或先在登录节点填充下载缓存：

```bash
bash core/doctor.sh download
bash core/doctor.sh conda --check
bash core/doctor.sh reference --check
bash core/doctor.sh demux --check
```

去掉 `--check` 即检查并按需修复。Conda 构建在 macOS 或 HPC 登录节点执行；已有 allocation
允许只读 `conda --check`。环境与 reference 的角色及恢复机制由 Doctor 统一检查。

## 4. 初始化与配置

Pipeline 与数据项目分开存放：

```bash
mkdir -p /path/to/Patient001
cd /path/to/Patient001
/path/to/Alopex/core/run_pipeline.sh --init-project
```

也可使用 `--init-project /path/to/Patient001`。初始化得到：

```text
00_config/     config.yaml、sample_manifest.tsv、孔板协议的 Barcode_Map.csv
01_raw/        原始 paired FASTQ
02_work/       可恢复中间结果与 SRD RNA 原始数据
03_results/    已发布科学结果
04_logs/       controller、规则、Slurm 日志及 provenance
05_tmp/        临时文件
06_downstream/ 下游 QC_Results_年月日_时分秒 目录
run_pipeline.sh
```

编辑 `00_config/config.yaml`，确认物种、reference 和路线。例如：

```yaml
schema_version: 3
species: hg38
references:
  hg38: /path/to/Alopex/resources/hg38_reference/hg38.primary.lambda.puc19.fa
analysis:
  protocol: cabernet
  cpg_coordinate_base: 1
cabernet:
  methylation_backend: bismark
```

配置优先级为**内置默认 → common → 选定协议**，列表整体替换。共用 demux、high-CpH、
Bismark 和 BISCUIT 参数放在 `common`；协议模块只写 backend 与差异。例如
`common.bismark.local_alignment: true` 对 Bismark 路线生效，
`cabernet.bismark.local_alignment: false` 只覆盖 Cabernet。
物种/reference、`analysis`、`runtime`、`slurm`、`retention` 为全局设置。
当前只接受 schema 3 模块结构；完整参数、默认值与覆盖范围见初始化生成的 config 注释。

- **Cabernet/SRD**：默认 Bismark，可选择 BISCUIT。Bismark 所有协议/模式固定 `--maxins 1000`；
  `local_alignment` 默认 `false`（combined-index end-to-end），`true` 使用 CT/GA local。
- **Droplet**：`analysis.protocol: droplet`，自动采用 Bismark non-directional、spacer=0；
  无需 Barcode Map，`rna_sample` 留空。初始化前已有 Droplet config 时跳过 Barcode Map。
  non-conversion 过滤后按胞嘧啶/UMI 汇集共识；BAM 仍保留 PCR 复制，独立分子计数见 CpG。
  `Duplicate_Pair_Rate%` 不适用，保留 NaN；MAPQ≥10、baseQ≥20，未解决的平票不调用。
- **TAPS**：`analysis.protocol: taps`，自动采用 Rastair，生成 SNP VCF/索引和 Non-CpG QC。
  Non-CpG 是可判链宿主 CHG/CHH 的原始转换率，零观测为缺失；不作 high-CpH 筛除。
  SNP 未调用位置不等于纯合参考，原生 PASS 也不保证高置信基因型；须结合深度、质量与
  OT/OB 链支持评估。需要额外重分析时设置 `retention.keep_final_bam: true`。

TAPS/Droplet 的 backend，以及 Droplet 的 spacer/library_type 属于固定策略，不接受 config
声明；无差异时可写 `taps: {}`、`droplet: {}` 或省略模块。科学阈值见 config 注释；资源档位由 Pipeline 按输入规模自动选择。

孔板 Barcode Map 至少含 `DNA_Barcode,PlateID,Cell_Order`，SRD 还需 `RNA_Barcode`。
barcode 为 8 或 10 bp A/C/G/T，同长度序列间 Hamming distance 必须大于 2；PlateID 和
Cell_Order 各自唯一。表头忽略外围空白、BOM 和大小写后仍须非空且唯一，每行列数一致。
PlateID 和样本名只允许字母、数字、点、下划线、连字符，禁止单独的 `.`/`..`。

## 5. 导入与选样

把 `.fastq.gz` 或 `.fq.gz` paired FASTQ 放入 `01_raw/`。常见命名：

```text
SampleA_L001_R1_001.fastq.gz  SampleA_L001_R2_001.fastq.gz
SampleA_R1.fastq.gz          SampleA_R2.fastq.gz
SampleA.R1.raw.fastq.gz      SampleA.R2.raw.fastq.gz
样本_SampleA_L1_1.fq.gz      样本_SampleA_L1_2.fq.gz
```

mate 支持 `R1/R2` 或 `1/2`，分隔符支持下划线、点、连字符；lane 为 1–999，chunk 可省略。
内部样本名仅去掉一次 `样本_`，保留其余完整基名。缺 mate、重复 segment、归一化身份冲突会失败。

```bash
./run_pipeline.sh --refresh-manifest
```

刷新扫描全部输入。多个子目录中同一样本的 plain 命名对会按目录顺序改成 chunk 命名；
相同副本、缺 mate 或混合命名会拒绝，失败会回滚已应用改名。已有 manifest 仅在行语义变化时
更新并写滚动 `.bak`；换行、排序及外围空白差异不会触发无效重算。

`sample_manifest.tsv` 固定六列：`sample_id/dna_raw_sample/rna_sample/species/protocol/notes`。
正式运行只分析其中选定的样本。可删除不需要的行，直接 dry-run/运行；未登记 FASTQ 留在原地会
提示并跳过。**选样后再 refresh 会重新补登记有效输入。**

Cabernet/TAPS 的 `rna_sample` 可关联外部 RNA_DARLIN 项目的原始 RNA basename；SRD 该列是
本地 RNA-enrichment tube。DNA/RNA 对齐与 Droplet 的独立文库映射统一见
[下游关联说明](downstream/README.md#关联-rna_darlin-结果)。

## 6. 运行与恢复

```bash
./run_pipeline.sh --dry-run
./run_pipeline.sh
```

项目从当前目录向上定位。launcher 不接受额外 Snakemake flag 或 target；初始化、刷新和
运行分别执行。macOS 使用 local executor；Linux 当前 shell/allocation 承载 controller，
规则作业经 Slurm 提交。`runtime.*`/`slurm.*` 预算含义见 config 注释；失败任务最多执行三次，内存按首次请求的 1/1.5/2 倍重试。

完整日志：`04_logs/controller/<run_id>.log`。终端每 60 秒、重定向时每 300 秒展示进度；
逐规则证据在 `04_logs/rules/` 和 `04_logs/slurm/`，资源规则的
`.attempt1/.attempt2/.attempt3` 保留本轮失败与重试原因。

失败后先检查对应规则及 attempt 日志，修正原因，再运行同一 launcher。保留 `02_work/`
和 `.snakemake/` 支持恢复；Snakemake 决定最小重算范围。残留锁必须先确认无运行中的
launcher，再移除项目 `.snakemake/locks`。不要仅凭进度百分比判定完成。

发布与 demux 清理都在 Snakemake 持锁的成功回调内执行。结果从稳定 stage 硬链接原子落到
`03_results/`，最后提交 `run_manifest.json`；两处必须位于同一文件系统。发布中断可复用
stage 重试。已完成的快路只读退出；清理失败保留残留与告警，待后续持锁发布再回收。
`retention.keep_final_bam: true` 要求保留 BAM/索引齐全，缺失时会进入 DAG 恢复。
改回 `false` 后，成功发布会在同一个锁内清理当前 backend 活跃细胞的最终 BAM/索引；
残留清理失败时，下次启动会重试。SRD 的 RNA BAM 不受此开关影响。
下游读取期间若上游重新发布，旧交付 context 会拒绝继续读取或提交 QC，需要重新启动下游。

## 7. 结果与完成确认

```text
03_results/
├── CpG/<Sample_ID>.cpg.tsv.zst
├── QC_Results/multiqc_report.html
├── QC_Results/sample_manifest.tsv
├── SNP/                 TAPS 默认 VCF；BISCUIT generate_snp=true 时 BED
└── run_manifest.json
```

CpG 唯一格式为 `chrom, pos, methyl, unmethyl`，tab 分隔、pos 默认 1-based，计数为非负整数。`analysis.cpg_coordinate_base: 0` 可改为 0-based（仅接受整数 0/1，省略默认 1）。BED/VCF/BAM 仍按各自格式；下游根据 sealed 交付的坐标身份将输入归一化为 SnapATAC2 所需的 0-based。
坏行会报告来源和行号并阻止提交；合法空结果保留表头。SRD 的 RNA FASTQ 和被筛除 BAM
固定保存在 `02_work/RNA/{FASTQ,BAM}/`，不进入公开 inventory。

完成需同时满足：

1. controller 成功退出；若通过 Slurm 作业承载 controller，该作业须为 `COMPLETED|0:0`。
2. `03_results/run_manifest.json` 为 `status: complete`，声明的结果与最终 manifest 齐全。
3. 再次运行 launcher 显示 `Alopex delivery is complete and current`，不启动科学作业。

MultiQC 汇总原生 mapping、pair funnel、重复、过滤、对照和适用的 M-bias；不同 backend
的 mapping 分母不同。科学解释与空结果的判定见 config 注释及上文协议说明；已验证数据范围见本版 Release。
超过 1,000 个活跃 DNA 细胞时，MultiQC 固定抽取 1,000 个细胞，所有逐细胞模块使用同一集合；
较小项目保留全部输入。报告首页标明抽样范围，文库级解复用概要仍为全库计数。
抽样只用于报告展示，最终样本清单、全部细胞的科学结果及 BAM 保留策略不变。

## 8. 下游

使用 `conda/current/notebook` 环境打开 `downstream/DNA_QC_Visualization.ipynb`，填写
`DNA_PROJECT_ROOT`、`DNA_PIPELINE_ROOT`。输出固定在当前交付的 `06_downstream/QC_Results_<生成时间>/`。
Notebook 启动、20 列 QC、RawAdata/single-CpG、缓存重算和 RNA 关联均见
[下游 README](downstream/README.md)。

## 许可证

Alopex 以 [GPL-3.0](LICENSE) 发布。本仓库为公开发布副本，与开发仓库按发布版本同步。
