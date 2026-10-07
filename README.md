# Alopex v14.3

Alopex 将单细胞 DNA 甲基化 paired FASTQ 转为单细胞 CpG、MultiQC 报告和可追溯的完整交付。

| 协议 | 后端 | 主要处理 |
|---|---|---|
| Cabernet | BISCUIT / Bismark | 孔板解复用、trimming、位置去重、high-CpH/non-conversion 过滤 |
| SRD | BISCUIT / Bismark | DNA/RNA tube 解复用、DNA 分析，并保留 RNA FASTQ 与被筛除 BAM |
| Droplet DD-MET5 | Bismark non-directional | DNA 独立谷底 calling、官方设计白名单、逐胞嘧啶 UMI 共识 |
| Cabernet–TAPS+ | BWA-MEM + Rastair | 标记重复、调用 5mC+5hmC；不执行 CpH/cDNA 筛除 |

用户入口只有 `core/doctor.sh`（构建与维护）和 `core/run_pipeline.sh`（项目运行）。
各参数的含义与默认值见初始化项目生成的 `00_config/config.yaml` 模板行内注释。

## 关于本仓库

本仓库是 Alopex 的 GPL-3.0 公开发布副本，包含完整可运行的 Pipeline：核心 workflow、
Rust demultiplexer、下游 QC、Doctor 端到端回归 fixture（真实 Cabernet 抽样 reads，
数据集命名与来源身份已中性化）以及默认 Barcode Map 模板。维护规范、参数表、流程图、
测试套件与内部 benchmark 归档不随公开副本分发；Doctor 的九个阶段在本副本均可完整执行。

Bismark paired-end 比对在所有协议与比对模式下固定使用 `--maxins 1000`（等价 `-X 1000`），最大插入片段长度为 1000 bp。

## 1. 安装

支持 macOS 和 Linux；Linux 正式运行使用 Slurm，需要 `sbatch`、`squeue`、`sacct`。
先准备可调用的 Conda、git、curl 和 gzip：

```bash
git clone https://github.com/Felix-owo/Alopex.git
cd Alopex
conda --version
```

也可使用源码归档部署，无需 `.git`。在目标机器上构建 Conda、Rust 和 reference indexes，
不要跨平台复制已构建环境。版本以 Git tag 为准。
Doctor 使用 Seed Conda 自带的 `ruamel.yaml` 抽取环境配置，无需在 Seed 中额外安装 PyYAML。
科学规则的 Python 入口固定到所选 backend 环境，避免 HPC 已激活环境的 PATH 顺序影响依赖加载。

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

## 3. 构建与检查

在 Pipeline 根目录执行：

```bash
bash core/doctor.sh
```

Doctor 检查或修复环境、Rust demultiplexer、hg38/mm10 reference 与三个 backend 的 indexes，
然后运行 A（hg38/Cabernet/Bismark）、B（mm10/SRD/BISCUIT）和 C（TAPS/Rastair）三条隔离路线。
登录节点的 A/B 为独立 Slurm workers，可并行；macOS 或已有 allocation 内串行，避免重复使用
整份内存预算。每条 route 最多 2 CPU，Bismark 16 GiB、BISCUIT/Rastair 24 GiB、20 min。
C 在 A/B 成功后执行；真实 TAPS 子集已完成生命周期与下游回归；完整 TAPS/SRD 文库仍未完成科学验收。

成功结尾必须为 `FINAL STATUS: READY`。详细日志在 `logs/doctor/`，成功清理测试沙箱，
失败保留沙箱和原生诊断。Doctor 每次都重跑三条路线，READY 只代表本次检查和回归通过。

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
06_downstream/ 下游按 delivery_id 创建
run_pipeline.sh
```

编辑 `00_config/config.yaml`，确认物种、reference 和路线。例如：

```yaml
schema_version: 2
species: hg38
references:
  hg38: /path/to/Alopex/resources/hg38_reference/hg38.primary.lambda.puc19.fa
analysis:
  protocol: cabernet
  methylation_backend: bismark
```

完整模板、默认值、可调参数及计数单位见初始化生成的 config 行内注释。
Droplet 使用 `protocol: droplet`、Bismark `non_directional`，不需要 Barcode Map，
`rna_sample` 留空；若在初始化前写好 Droplet config，初始化会直接跳过 Barcode Map。
Droplet 在 non-conversion 过滤后按胞嘧啶位点/UMI 汇集分子共识，保留不同 read pairs 的联合覆盖。
固定 MAPQ≥10、baseQ≥20；未解决的共识平票不调用。保留 BAM 时其中仍有 PCR 复制，
独立分子计数以 CpG 输出为准；QC 的 pair 去重率为空，位点分子统计见 MultiQC 与 final manifest。
TAPS 使用 `protocol: taps`、`methylation_backend: rastair`。

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
规则作业经 Slurm 提交。`runtime.*`/`slurm.*` 预算及重试策略的含义见 `00_config/config.yaml` 模板行内注释。

完整日志：`04_logs/controller/<run_id>.log`。终端每 60 秒、重定向时每 300 秒展示进度；
进度按有界缓冲只读取新增日志。逐规则证据在 `04_logs/rules/` 和 `04_logs/slurm/`，资源规则的
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
├── SNP/                 仅 BISCUIT 且 generate_snp=true
└── run_manifest.json
```

CpG 唯一格式为 `chrom, pos, methyl, unmethyl`，tab 分隔、pos 为 0-based，计数为非负整数。
坏行会报告来源和行号并阻止提交；合法空结果保留表头。SRD 的 RNA FASTQ 和被筛除 BAM
固定保存在 `02_work/RNA/{FASTQ,BAM}/`，不进入公开 inventory。

完成需同时满足：

1. controller 成功退出；HPC controller 作业为 `COMPLETED|0:0`。
2. `03_results/run_manifest.json` 为 `status: complete`，声明的结果与最终 manifest 齐全。
3. 再次运行 launcher 显示 `Alopex delivery is complete and current`，不启动科学作业。

MultiQC 汇总原生 mapping、pair funnel、重复、过滤、对照和适用的 M-bias；不同 backend
的 mapping 分母不同。科学解释与空结果的判定以 `00_config/config.yaml` 模板注释为准；已验证数据范围见本版 release 说明。

## 8. 下游

使用 `conda/current/notebook` 环境打开 `downstream/DNA_QC_Visualization.ipynb`，填写
`DNA_PROJECT_ROOT`、`DNA_PIPELINE_ROOT`。输出固定在当前交付的 `06_downstream/<delivery_id>/`。
Notebook 启动、20 列 QC、RawAdata/single-CpG、缓存重算和 RNA 关联均见
[下游 README](downstream/README.md)。

## 许可证

Alopex 以 [GPL-3.0](LICENSE) 发布。本仓库为公开发布副本，与开发仓库按发布版本同步。
