# Alopex v15.3

正式源码以 Git tag `v15.3` 为准。Alopex 将单细胞 DNA 甲基化 paired FASTQ 转为单细胞 CpG、QC 报告和可追溯的完整交付。

| 协议 / `analysis.protocol` | 后端 | 主要处理 |
| --- | --- | --- |
| Cabernet / `cabernet` | 默认 Bismark，可选 BISCUIT | 孔板解复用、trimming、位置去重、high-CpH/non-conversion 过滤 |
| SRD / `srd` | 默认 Bismark，可选 BISCUIT | DNA/RNA tube 解复用、DNA 分析，保留 RNA FASTQ 与被筛除 BAM |
| Droplet DD-MET5 / `droplet` | 固定 Bismark non-directional | DNA 独立谷底 calling、官方设计白名单、逐胞嘧啶 UMI 共识 |
| Cabernet–TAPS+ / `taps` | 固定 BWA-MEM + Rastair | 标记重复，调用 CpG/SNP，统计 Non-CpG；不作 high-CpH 筛除 |

| 查阅内容 | 文档 |
| --- | --- |
| 配置默认值与可调参数 | 初始化生成的 `00_config/config.yaml` 与内置 schema |
| Notebook、矩阵、缓存重算与 RNA 关联 | [下游 README](downstream/README.md) |
| 升级说明与验证边界 | [本版 Release](https://github.com/Felix-owo/Alopex/releases/tag/v15.3) |

## 1. 安装

准备可调用的 Conda、git、curl 和 gzip；Pipeline 与数据项目分开存放。

```bash
git clone https://github.com/Felix-owo/Alopex.git
cd Alopex
conda --version
```

| 平台 / 部署方式 | 要求 |
| --- | --- |
| macOS | 本机执行 workflow |
| Linux / HPC | 使用 Slurm，需 `sbatch/squeue/sacct` 与作业提交权限；已有 allocation 承载 controller，规则作业仍经 Slurm 提交 |
| 源码归档 | 可无 `.git` 部署；环境、Rust 二进制和 reference indexes 必须在目标平台构建 |

## 2. 准备标准 reference

自动构建的标准 reference 为 `hg38` 和 `mm10`。首次运行 Doctor 前，将以下六个文件**下载并解压**到 `resources/reference_source/`；Pipeline 不自动下载原始 reference。

| 固定文件名 | 对应版本 / 来源 | 下载链接 |
| --- | --- | --- |
| `hg38.fa` | UCSC hg38 | [hg38.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) |
| `mm10.fa` | UCSC mm10 | [mm10.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz) |
| `lambda.fa` | NCBI J02459.1 | [FASTA](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=J02459.1&rettype=fasta&retmode=text) |
| `pUC19.fa` | NCBI L09137.2 | [FASTA](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=L09137.2&rettype=fasta&retmode=text) |
| `hg38.genes.gtf` | GENCODE v48 / GRCh38 | [GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz) |
| `mm10.genes.gtf` | GENCODE M25 / GRCm38 | [GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz) |

## 3. 构建与检查

在 Pipeline 根目录执行：

| 命令 | 用途 |
| --- | --- |
| `bash core/doctor.sh` | 检查或修复 Conda、Rust、reference/indexes，然后运行下表三路线 |
| `bash core/doctor.sh conda --check` | 只读检查环境 |
| `bash core/doctor.sh reference --check` | 只读检查 hg38/mm10 reference 及 indexes |
| `bash core/doctor.sh demux --check` | 只读核对 Rust 源码与二进制身份 |
| `bash core/doctor.sh download` | 预热 Conda/Rust 下载缓存，不编译、不发布、不提交作业 |

去掉 `--check` 可维护对应段；显式 `demux` 会测试、构建并原子替换生产二进制。完整 Doctor 先检查并复用就绪资源。Conda 构建在 macOS 或 HPC 登录节点执行，已有 allocation 可作只读检查。

| Doctor 路线 | 内容 |
| --- | --- |
| A | hg38 / Cabernet / Bismark：真实 Cabernet 抽样与合成对照 |
| B | mm10 / SRD / BISCUIT：真实 Cabernet DNA 抽样与合成 RNA 布局/对照 |
| C | hg38 / TAPS / Rastair：独立合成真值 |

HPC 登录节点提交路线 workers；macOS 或已有 allocation 内串行执行。成功结尾为 `FINAL STATUS: READY`，日志在 `logs/doctor/`；成功清理沙箱，失败保留诊断。READY 证明本次工具链与固定小样本通过，不代表真实全库科学等价或性能验证。

## 4. 初始化与配置

在独立数据项目中初始化：

```bash
mkdir -p /path/to/Patient001
cd /path/to/Patient001
/path/to/Alopex/core/run_pipeline.sh --init-project
```

也可用 `--init-project /path/to/Patient001`。之后使用项目根目录的 `./run_pipeline.sh`。

| 项目路径 | 用途 |
| --- | --- |
| `00_config/` | `config.yaml`、六列 `sample_manifest.tsv`、孔板协议的 `Barcode_Map.csv` |
| `01_raw/` | 原始 paired FASTQ |
| `02_work/` | 可恢复中间结果、SRD RNA 数据与稳定结果 stage |
| `03_results/` | 当前完整交付的科学结果和 QC |
| `04_logs/` | controller、规则、Slurm 日志与 provenance |
| `05_tmp/` | 作业临时文件 |
| `06_downstream/` | 按交付绑定的下游 QC 与矩阵 |

在生成的 `00_config/config.yaml` 上修改关键项，例如：

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

| 配置主题 | 使用规则 |
| --- | --- |
| 结构与覆盖 | 只接受 schema 3；内置默认 → `common` → 选定协议覆盖，列表整体替换，固定协议策略最终生效 |
| 共用 / 全局设置 | demux、high-CpH、Bismark、BISCUIT 参数放 `common`；物种/reference、`analysis`、`runtime`、`slurm`、`retention` 为全局设置 |
| Cabernet / SRD | 可选 `methylation_backend: bismark` 或 `biscuit`；Bismark 默认 combined-index end-to-end，`bismark.local_alignment: true` 改为 CT/GA local；Bismark 各协议均固定 `--maxins 1000` |
| Droplet | backend、spacer=0、non-directional 固定，不声明 `methylation_backend`、`demux.dna_w_spacer_len` 或 `bismark.library_type`；无需 Barcode Map，`rna_sample` 留空 |
| TAPS | backend 固定 Rastair，不声明 `methylation_backend`；输出 SNP VCF/索引与 Non-CpG QC |
| 无差异模块 | `taps` / `droplet` 可省略或为 `{}`；完整可调项见初始化生成的 config 模板和内置 schema |

| 孔板 Barcode Map | 要求 |
| --- | --- |
| 必需列 | `DNA_Barcode,PlateID,Cell_Order`；SRD 另需 `RNA_Barcode` |
| Barcode | 8 或 10 bp A/C/G/T；同长度序列间 Hamming distance 必须大于 2 |
| 身份 | PlateID、Cell_Order 各自唯一；Cell_Order 保留字符串，不丢前导零 |
| CSV 格式 | 表头去 BOM、外围空白并忽略大小写后非空且唯一；各行列数一致 |
| 名称 | PlateID 和样本名仅用 ASCII 字母、数字、点、下划线、连字符，禁止单独的 `.` / `..` |

## 5. 导入与选样

将 gzip 压缩的 `.fastq.gz` / `.fq.gz` 双端输入放入 `01_raw/`。

| 命名形式 | R1 示例 | R2 示例 |
| --- | --- | --- |
| lane / chunk | `SampleA_L001_R1_001.fastq.gz` | `SampleA_L001_R2_001.fastq.gz` |
| 普通配对 | `SampleA_R1.fastq.gz` | `SampleA_R2.fastq.gz` |
| raw 后缀 | `SampleA.R1.raw.fastq.gz` | `SampleA.R2.raw.fastq.gz` |
| 物理前缀 | `样本_SampleA_L1_1.fq.gz` | `样本_SampleA_L1_2.fq.gz` |

mate 支持 `R1/R2` 或 `1/2`，分隔符支持 `_ . -`；lane 为 1–999，chunk 可省略、显式时为三位数字。内部身份仅去掉一次 `样本_`，保留其余基名；缺 mate、重复 segment 或归一化冲突会拒绝。

```bash
./run_pipeline.sh --refresh-manifest
```

| 操作 / 清单 | 行为 |
| --- | --- |
| 刷新 | 全量发现输入，生成或更新 `00_config/sample_manifest.tsv`；多个目录中的同样本 plain 配对按目录顺序改为 chunk 命名 |
| 冲突与备份 | 相同副本、缺 mate、混合命名等冲突拒绝，已应用改名会回滚；清单内容未变不改写，有效旧清单变更时保留滚动 `.bak` |
| 清单六列 | `sample_id / dna_raw_sample / rna_sample / species / protocol / notes` |
| 选样 | 删除不需要的行后直接 dry-run / 运行；未登记 FASTQ 可留在原地，会提示并跳过 |
| 再次刷新 | **会补登记有效输入，可能重新加入手动删行的样本** |

| 协议 | `rna_sample` 的含义 |
| --- | --- |
| Cabernet / TAPS | 可选的外部 RNA_DARLIN **原始 RNA basename**，不是 RNA sample_id 或路径；不产生本地 RNA route，也不会被 refresh 自动当成本地 DNA 登记 |
| SRD | 本地 RNA-enrichment tube，需对应本地输入 |
| Droplet | 留空；配套 RNA 独立 calling，通过文库映射与完整 CB17 关联 |

DNA/RNA 关联与更新步骤见[下游关联说明](downstream/README.md#关联-rna_darlin-结果)。

## 6. 运行与恢复

```bash
./run_pipeline.sh --dry-run
./run_pipeline.sh
```

项目从当前目录向上定位；初始化、刷新和运行分别执行。launcher 不接受额外 Snakemake flag 或 target。

| 情况 | 处理方式 |
| --- | --- |
| 查看进度 / 失败 | 完整日志在 `04_logs/controller/<run_id>.log`，规则与 Slurm 日志在 `04_logs/{rules,slurm}/`；`.attempt1/.attempt2/.attempt3` 保留本轮重试证据 |
| 规则失败 | 先查对应日志并修正原因，再运行同一 launcher；保留 `02_work/` 和 `.snakemake/`，由 Snakemake 决定重算范围 |
| 残留锁 | 确认没有运行中的 launcher 后，才手动移除项目 `.snakemake/locks` |
| 发布中断 | 保留 `02_work/results_stage/` 后重跑；stage 与 `03_results/` 必须同一文件系统，公开 `run_manifest.json` 最后提交 |
| 最终 DNA BAM | `retention.keep_final_bam: true` 保留 BAM 及适用索引，缺失时进入 DAG 恢复；`false` 在成功发布后回收当前 backend 活跃细胞的最终 BAM/索引 |
| SRD RNA | FASTQ 和被筛除 BAM 固定保存在 `02_work/RNA/{FASTQ,BAM}/`，不受最终 DNA BAM 开关影响 |
| 上游重新发布 | 旧下游 context 会拒绝继续读取或提交 QC，需重新启动下游 |
| reference 构建中断 | 重跑对应 Doctor 段，恢复未完成事务并复验；未完成不报告 READY |

## 7. 结果与完成确认

| 结果路径 | 内容 |
| --- | --- |
| `03_results/CpG/<Sample_ID>.cpg.tsv.zst` | Tab 分隔四列：`chrom, pos, methyl, unmethyl` |
| `03_results/QC_Results/sample_manifest.tsv` | 最终细胞清单、身份、计数与 QC |
| `03_results/QC_Results/multiqc_report.html` | 汇总报告 |
| `03_results/SNP/` | TAPS 默认 VCF/索引；BISCUIT 开启 `generate_snp` 时为 BED/索引 |
| `03_results/run_manifest.json` | 当前交付的完整标记、配置身份及结果 inventory |

| 科学结果解读 | 边界 |
| --- | --- |
| CpG 坐标与计数 | `analysis.cpg_coordinate_base` 仅接受整数 0/1，默认 1；methyl/unmethyl 为非负整数，有效位点 coverage>0，合法空文件保留表头 |
| 坐标传递 | RawAdata / SingleCpG 沿用 sealed 交付的 base；BED、VCF、BAM 仍按各自格式解释 |
| Droplet 计数 | BAM 保留 PCR copies，CpG 按胞嘧啶/UMI 共识计独立分子；MAPQ≥10、baseQ≥20，未解决平票不调用；`Duplicate_Pair_Rate%` 为 NaN |
| TAPS 解读 | Non-CpG 为可判链宿主 CHG/CHH 原始转换率，零观测为缺失，不直接等同真实 mCH 或 cDNA 污染率；未调用 SNP 不等于纯合参考，PASS 仍需结合深度、质量与链支持 |
| MultiQC | backend 原生 mapping 分母不同；活跃 DNA>1000 时固定展示同一组1000细胞，文库概要仍为全库，全部细胞科学结果与 BAM 策略不变 |

| 完成判据 | 必须满足 |
| --- | --- |
| 执行结束 | controller 成功退出；若由 Slurm 作业承载，该作业为 `COMPLETED`、exit `0:0` |
| 交付完整 | `run_manifest.json` 为 `status: complete`，声明的结果与最终清单齐全 |
| 再次启动 | 显示 `Alopex delivery is complete and current`，不启动科学作业 |

## 8. 下游

| 操作 / 结果 | 使用方式 |
| --- | --- |
| 启动 | 用 `conda/current/notebook` 环境打开 `downstream/DNA_QC_Visualization.ipynb`，填写 `DNA_PROJECT_ROOT`、`DNA_PIPELINE_ROOT`；按需设置 `TSS_BED` / `SKIP_TSS` |
| 输出代际 | 写入 `06_downstream/QC_Results_<生成时间>/`；同交付复用目录/cache，新交付新建，不按 mtime 猜测输入，不写回公开结果 |
| TSS 抽检 | 默认按 CloneID 抽10%，不足10个全检、超过1000个取100；`TSS_SAMPLE_SEED=0` 可复现，`TSS_SAMPLE_FRACTION=1.0` 全检；抽样只影响 TSS，Raw/QC/Gini/Single 仍处理全部细胞 |
| 矩阵与零值 | 位点坐标跟随交付；有覆盖全 unmethyl 保留显式0，缺覆盖不补0；QC 位点均值、SingleCpG 区间均值与合并 reads 的比率不同 |
| 缓存重建 | CpG统计统一用`RECOMPUTE_CPG_STATS`/`--recompute-cpg-stats`同时重建位点统计与信号组成；过期缓存按 Processor 提示重算；Raw 过期时先完成 QC Processor 的必要重建，再用 `RECOMPUTE_SINGLECPG=True` 重建 SingleCpG |
| 绘图 | Final QC 按 CloneID 展示计数直方图，CpG sites 用 log 横轴；CpG density 展示 input pairs / unique CpG sites 散点及拟合，支持 PDF |
| RNA 关联 | 使用当前完整交付的 sealed manifest 与同代 DNA QC；更新关联后先重发布 DNA 和当前代际 downstream，再整合 RNA |

完整参数、矩阵解释与操作步骤见[下游 README](downstream/README.md)。

## 许可证

Alopex 以 [GPL-3.0](LICENSE) 发布。本仓库按发布版本同步核心 workflow、下游代码、Doctor fixture 和默认 Barcode Map。
