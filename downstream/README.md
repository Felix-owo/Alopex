# Alopex v16.0.1 下游 QC 与可视化

正式源码以 Git tag `v16.0.1` 为准。本文说明 Notebook、矩阵、QC 分母、缓存重算与 RNA 关联；安装与主流程见[主 README](../README.md)。

| 输入边界 | 要求 |
| --- | --- |
| 完整交付 | `03_results/run_manifest.json` 必须为当前 complete 交付；从其解析 final manifest、CpG inventory 与绑定的配置快照 |
| 快照与路径 | 交付绑定的 snapshot/config（通常在 `04_logs/provenance/`）及 reference/FAI 须在记录路径可读，项目/Pipeline 路径须匹配；只复制 `03_results/` 不能完成迁移 |
| 可清理中间件 | 下游不依赖 `02_work/` 中的 BAM、demux FASTQ 等中间结果 |
| 写出边界 | 只写 `06_downstream/`，不修改公开科学结果；上游重发布后须重新启动下游 |

## 运行

使用 Doctor 管理的 `conda/current/notebook` 环境（规格在 `core/envs.yaml`）：

```bash
PIPELINE=/path/to/Alopex
"$PIPELINE/conda/current/notebook/bin/python" \
  "$PIPELINE/downstream/downstream_qc.py" \
  --project-dir /path/to/Patient001
```

Processor 完成 QC 后，可打开 `downstream/DNA_QC_Visualization.ipynb` 绘图；Notebook 也可直接启动 Processor。交互式运行前先按下节注册并选择专用 kernel，子进程沿用 Notebook 的 Python。

| Notebook 第一格 | 设置与定位 |
| --- | --- |
| `DNA_PROJECT_ROOT` | 建议填写 `Path("/path/to/Patient001")`；显式值优先于同名环境变量，错误时不回退 |
| `DNA_PIPELINE_ROOT` | 建议填写 `Path("/path/to/Alopex")`；显式指定仍须通过交付与 Pipeline 身份验证 |
| 两根留空 | 项目从 Jupyter 当前目录及父目录查找 complete delivery，Pipeline 从 sealed `source_root` 定位；先检查同名环境变量 |
| `TSS_BED=None` | 自动使用当前物种的 TSS reference；可指定同一 genome build 的外部 BED |
| `SKIP_TSS=False` | 确实无需 TSS 时显式改为 `True`；不能同时指定外部 `TSS_BED` |
| CLI `--workers=0` | 自动选择，最多32；显式值同样受 affinity、Slurm allocation 与调用方预算限制，每个 worker 内 native 线程为1 |

### Jupyter / VS Code kernel

Notebook 使用名为 `alopex-qc`、显示为 `Alopex QC` 的 kernel。以下命令将启动路径固定到 `conda/current`，环境重建后无需重新注册；Pipeline 迁移目录后重新注册并选择 kernel。

```bash
PIPELINE=/path/to/Alopex
"$PIPELINE/conda/current/notebook/bin/python" - "$PIPELINE" <<'PY_KERNEL'
import json
import sys
from pathlib import Path
from jupyter_core.paths import jupyter_data_dir

pipeline = Path(sys.argv[1]).absolute()
kernel = Path(jupyter_data_dir()) / "kernels/alopex-qc"
kernel.mkdir(parents=True, exist_ok=True)
(kernel / "kernel.json").write_text(json.dumps({
    "argv": [str(pipeline / "conda/current/notebook/bin/python"),
             "-m", "ipykernel_launcher", "-f", "{connection_file}"],
    "display_name": "Alopex QC", "language": "python",
    "metadata": {"debugger": True},
}, indent=2))
PY_KERNEL
```

### 缓存与强制重算

默认按 sealed 输入、参数与算法身份复用缓存，以下开关均为 `False`。强制重算只影响下游；新 Raw generation 会使 Gini/TSS 重算，canonical CpG 统计仍可独立复用。

| Notebook 开关 | CLI 对应项 | 重建对象 |
| --- | --- | --- |
| `RECOMPUTE_RAW_ADATA` | `--recompute-raw-adata` | 从 sealed CpG 重建 `RawAdata.h5ad` |
| `RECOMPUTE_CPG_STATS` | `--recompute-cpg-stats` | 同时重建 CpG 位点/甲基化统计及信号组成；任一缓存缺失/损坏也会自动联合重建 |
| `RECOMPUTE_GINI` | `--recompute-gini` | Gini |
| `RECOMPUTE_TSS` | `--recompute-tss` | TSS profile |
| `RECOMPUTE_SINGLECPG` | 无；Notebook 最后一格使用 | `SingleCpG_Adata.h5ad`；完成后恢复 `False` |

| 情况 | 正确处理 |
| --- | --- |
| 正常重复运行 | 自动复用；输入身份复用 sealed inventory，不为 hash 再次遍历全部 CpG |
| 调整 TSS BED、比例或种子 | 仅 TSS 自动失效，无需强制重建 Raw/QC/Gini/Single |
| QC 完成身份缺失或计算期间变化 | 先重新运行 Processor 成功完成，再加载 QC 或生成 SingleCpG；Notebook 与 Single 使用同一已验证上下文 |
| Raw 输入或算法身份过期 | 先运行 Processor 完成 Raw 必要重建，再处理 SingleCpG；不能仅改 metadata 或用 Single 强制开关绕过 |
| SingleCpG 来源未知或身份变化 | 按提示设置 `RECOMPUTE_SINGLECPG=True`；身份绑定交付、Raw generation、BED、参数、算法与科学库版本 |
| 仅细胞注释变化 | 矩阵可复用，只更新 SingleCpG `obs`；原有 QC 列仍须与当前磁盘 QC 一致。NaN/null 仅在比较时视为同一缺失，真实零值独立；未变化不改写矩阵或完成记录 |
| 科学统计重算中断 | 旧完成身份撤销，下一次重新计算；有效 Raw 仍可复用。全部缓存命中时的导出失败不破坏有效缓存 |
| 计算中参考变化 | FAI、TSS BED或Single BED在身份读取到提交期间被替换/改写时拒绝提交完成身份；确认参考稳定后重新运行。仅复验本次文件状态，不增加哈希遍数；已替换但未完成的矩阵不能作为有效缓存复用 |
| 缓存计数自相矛盾 | 零位点/零 coverage 必须一致，宿主 coverage 须在 `sites` 至 `65535×sites`；异常从 canonical CpG 重算，不静默接受 |

| 科学输入校验 | 保留的含义 |
| --- | --- |
| CpG 四列 | `chrom, pos, methyl, unmethyl`；位置与计数为 ASCII 十进制整数，计数非负，coverage 为正且不超过65535 |
| 位点与 reference | 按 sealed base 和 FAI 校验真实边界；每条染色体位置严格递增，重复/乱序、未知 contig、缺表头、空行或解压失败均拒绝 |
| 合法空细胞 | 允许仅表头的空 CpG；不把格式失败当作 zero-mapped，缓存命中不重复解压 |
| 解压进程 | 数据保持流式读取，诊断暂存自动清理的临时文件，避免大量stderr堵塞读取。完整读取及空表必须检查解压退出码并报告原生错误；非空首行探测可提前结束，解析失败或早停均终止并回收本次解码器 |
| 原生导入 | Step 2 只探测表头、首条记录与空表，直接链接原 zst 调用 `snap.pp.import_values(chunk_size=200)`；私有链接目录优先使用 SLURM_TMPDIR/节点临时目录，不可写时回退输出目录，不复制 CpG 内容。不做 Python 全表预扫。完整行/边界/顺序/解压校验随原有 Step 3 统计完成，失败不提交下游完成身份；有效 Raw/统计缓存命中不解压 |
| 流式统计 | Step 3 在有界文本块中用编译循环完成同一套严格校验和逐行 Kahan 求和；释放 GIL 后采用最多32个工作线程，并按实际 CPU/调用方预算限幅，解码器单线程。跨块保持位点顺序和浮点累计顺序，任务全部退出后恢复 native 线程环境；失败关闭解压进程，不提交完成身份 |
| contig 分类 | 保留 FAI 全部原名，包括点号、`NA` 与前导零；controls 仅精确匹配 `lambda` / `pUC19`，含其子串的宿主 chr 不算 controls；mt 按 chrM/MT/chrMT 分类，其他 decoy 仅进 total |
| 身份与缺失 | `Sample_ID/CloneID/PlateID` 按字符串保留 `001`、`NA`、`nan`、`None`；数值列按数值/缺失解析，浮点缓存保持 round-trip 精度 |


Notebook 在本地与 HPC 均通过 `run_qc_processor(PATHS)` 启动 `PATHS.pipeline_root/downstream/downstream_qc.py`，
子进程优先使用同一 Pipeline 的 `conda/current/control/bin` 工具，避免继承 Jupyter 服务的其他 Conda 环境。
Raw 重建与 CpG 统计统一使用上述已验证方案，无需另设分块参数；第一格的 `Pipeline:` 应指向已更新的源码目录。
原 Notebook 已完成同资源全库对照；本地通过科学/并发回归，尚无本地2707细胞完整性能测量，不承诺所有平台固定耗时。
有效缓存仍优先复用；过期缓存按 Processor 提示处理。更新本版后请结束旧写任务并重启 Notebook kernel；SingleCpG 提示来源过期时需显式重建。


### 保存供下游复用的 single-CpG 矩阵

运行 Notebook 最后一格生成 `SingleCpG_Adata.h5ad`，CLI Processor 不生成它。矩阵保留全部细胞，按 `Sample_ID` 严格对齐完整 QC、`HQ` 与自定义注释，不按 HQ 筛选。

| 对象 | 使用与科学边界 |
| --- | --- |
| 计算 | `prepare_single_cpg_adata(PATHS, df_qc)`；固定 `value_type="fraction"`、`summary_type="mean"`、`chunk_size=500`、`inplace=False` |
| 完成记录 | `SingleCpG_Adata.input_schema.json` 与矩阵配套提交；未知/失配来源须显式重建 |
| 坐标 | 跟随 sealed `analysis.cpg_coordinate_base`；0-based `chr1:0-2` 为半开区间，1-based `chr1:1-2` 为闭区间，表示同一 CpG，数值和列顺序相同 |
| 坐标 metadata | `uns["cpg_coordinate_base"]`、`uns["cpg_interval_convention"]`、`uns["source_cpg_representation"]` 记录实际表示；传给 BED/0-based 工具前须转换 |
| Raw 索引 | 真实长度 `reference_bp_length=L`，索引长度 `reference_seq_length=L+base`，base1 时索引0留空；原生 `reference_sequences` 行序与索引长度决定染色体归属，不可按 FAI 重排 |
| 读取 | helper 返回路径并关闭自身句柄；用打印出的 `SINGLECPG_ADATA_PATH` 执行 `snap.read(..., backed="r")` |
| 后续注释 | 只读对象不能回填标签；使用独立可写副本或合适大小的内存子集。MethylTree 分箱后仍按细胞 ID 核对并对齐 `obs` |

### TSS profile

| reference / 参数 | 约定 |
| --- | --- |
| 默认 TSS BED | `resources/<genome>_reference/tss/<genome>_TSS_2000_2000_20.bed`；链方向一致的 ±2000 bp / 20 bp bins，第9列为 offset |
| 构建来源 | Doctor 从 hg38 GENCODE v48 或 mm10 GENCODE M25 的 protein-coding `gene` 记录生成；原始 FASTA/GTF 按主 README 准备，不自动下载 |
| single-CpG BED | `resources/<genome>_reference/cpg/<genome>.single_cpg.bed.gz`，BED3 / 0-based half-open `[C, G+1)`，不含 lambda/pUC19 |
| 资源分发 | 生成的 `resources/*_reference/` 留在目标平台，不进入源码包或强制提交 Git |
| 外部资源 | CLI `--tss-bed /path/to/genome_TSS_2000_2000_20.bed` 或 Notebook `TSS_BED`；须匹配 sealed reference build |
| BED 校验 | genomic interval 非负整数、非空且唯一，第9列 offset 为有限32-bit整数；先检查全部行，再限制原生reference范围；实际抽检子集另筛无覆盖区间，全检不额外预扫描Raw |
| 显式跳过 | `--skip-tss` 或 `SKIP_TSS=True`；清理该交付旧 TSS CSV/PDF，避免把旧曲线当作新结果 |

Notebook 第二格默认 `TSS_SAMPLE_FRACTION=0.10`、`TSS_SAMPLE_SEED=0`，对应 CLI `--tss-sample-fraction` / `--tss-sample-seed`。比例须有限且在 `(0,1]`，种子为非负整数。

| 每克隆细胞数 N | 默认10% | 调整比例后的规则 |
| --- | --- | --- |
| N<10 | 全部 | 仍全部 |
| 10≤N≤1000 | `ceil(N×0.10)` | `min(100, ceil(N×fraction))` |
| N>1000 | 100 | `min(100, ceil(N×fraction))` |
| 任意 N，fraction=1.0 | — | 全部，不受100上限限制 |

| 抽检与结果 | 解读 |
| --- | --- |
| 选样 | 每个原始 CloneID 内按 `SHA256(seed, Sample_ID)` 排序；行序不影响名单，含空 CpG 细胞，不按 coverage 选样 |
| 影响范围 | 仅 TSS 用选中细胞；Raw、完整 QC、Gini、SingleCpG 均保留全部细胞 |
| 抽检记录 | `.cache/DNAme_QC_Input_Identity.json` 的 `tss_sampling` 保存方法、比例、种子、逐克隆数量与选中 ID |
| profile 行数 | CSV 仅记录选中且有覆盖的 bins，细胞数可小于抽检数；全体无覆盖时仅表头 CSV 合法并可复用，旧 PDF 清理 |
| 零值与坐标 | 有 coverage 而甲基化为0仍保留；无覆盖不补0。区间自动适配 sealed base，标准 BED 不变，相对 TSS 距离不加减1 |

抽检副本只读取选中细胞的 CSR 行切片，保留原始整数记录、全局位点索引和原生 reference 顺序，再交给原生 TSS 矩阵计算；避免离散行 `subset` 先把全库解压到内存。复制内存随抽中细胞的实际观测数增长，完整库仍用已有 backed Raw；不需要新增资源参数或修改 Notebook。

## 输出

输出目录为 `<project>/06_downstream/QC_Results_YYYYMMDD_HHMMSS/`，时间取首次创建的本地时间，同秒冲突加微秒。`.delivery.json` 绑定交付：同交付复用目录/cache，新交付新建；多个匹配会拒绝，不按 mtime 或目录名猜测。

同一交付的 QC Processor 与 SingleCpG 更新只能有一个写入者；重复启动遇到 `Downstream QC is busy` 时，等当前任务结束后重试。目录初始化也受互斥保护，避免两个 Notebook 同时创建重复目录。`.qc-write.lock` 是持久锁文件，**不要删除或替换**；进程退出后系统会释放锁，文件仍存在不代表任务在运行。最终 QC CSV 替换前撤销旧完成身份，提交中断后重跑 Processor 恢复，不应继续消费未完成结果。更新此模块后应结束旧写任务并重启 Notebook kernel，确保所有写入者加载同一实现。

| 文件 / 目录 | 用途 |
| --- | --- |
| `DNAme_QC_Information.csv` | 下游与 adata 合并的唯一 QC 表，固定20列，见下表 |
| `TSS_Profile_Information.csv` | 未跳过 TSS 时生成，有覆盖的抽检 bins |
| `RawAdata.h5ad` / `RawAdata.input_schema.json` | 全细胞 canonical CpG 导入矩阵与复用记录 |
| `SingleCpG_Adata.h5ad` / `SingleCpG_Adata.input_schema.json` | Notebook 最后一格生成的全细胞 single-CpG 矩阵与完成记录 |
| `.cache/` | `Mapping_Stats.csv`、`CpG_Signal_Composition.csv`、`Gini_Index.csv`、`DNAme_QC_Input_Identity.json` 等内部缓存，不手动改写 |
| `plots/` | QC、Final QC、孔板、信号组成、TSS、CpG density 与 PCA 等 PDF |

| QC 分组 | 列名（实际列序由 `QC_INFO_COLUMNS` 定义） |
| --- | --- |
| 身份 | `Sample_ID`、`CloneID`、`PlateID` |
| mapping | `Native_Mapping%` |
| 宿主 | `Sample_Unique_CpG_Sites`、`Sample_Meth_CpG_Rate%` |
| lambda | `Lambda_Unique_CpG_Sites`、`Lambda_Meth_CpG_Rate%` |
| pUC19 | `pUC19_Unique_CpG_Sites`、`pUC19_Meth_CpG_Rate%` |
| CpH | `Non_CpG_Methylation%`、`High_CpH_Flag_Rate%` |
| 信号组成 | `Signal_Host_Rate%`、`Signal_Lambda_Rate%`、`Signal_pUC19_Rate%`、`Signal_mtDNA_Rate%` |
| 均匀度 | `Gini_Index` |
| pair 指标 | `Duplicate_Pair_Rate%`、`Trim_Retention%`、`Final_Pair_Yield%` |

完整 read-pair funnel 与 backend policy/unit 保留在主 Pipeline final manifest/MultiQC，不重复塞入20列表。Droplet 的 `Duplicate_Pair_Rate%` 为 NaN；Non-CpG 来自过滤后 CHG/CHH 胞嘧啶/UMI 共识，分子观测数不与 read-pair 数混算。

## 指标与筛选

| 指标 | 分母与零值含义 |
| --- | --- |
| QC 甲基化率 | 宿主、lambda、pUC19 分别对已覆盖位点的 `m/(m+u)` 取均值再×100；无观测为 NaN，不是合并 reads 的比例 |
| SingleCpG | 参考区间内已覆盖位点的 fraction 均值；如 `1/4` 与 `9/12` 的均值为0.5，不是 `10/16`；有覆盖全 unmethyl 保显式0，缺覆盖不补0 |
| TSS | 每个 offset 对有覆盖 bins 的 fraction 取均值；count 为甲基化信号之和，不是 read pairs 或总 coverage |
| Gini | 流式读取 backed CSR，以原生 float32 顺序累加有覆盖 bins 的 `m+u`，省去临时分箱矩阵；不纳入未覆盖 bins，有 coverage 的零甲基化不等同缺失。 |
| 信号组成 | host/lambda/pUC19/mtDNA coverage 各除以全部 contig 的总 coverage 再×100，其他 decoy 只进分母；不是位点比例 |
| 原生 mapping | 各 backend 单位/分母不同，不视为同义指标或默认跨 backend HQ 阈值 |
| TAPS controls | lambda 假阳性 / pUC19 转换率为 `sum(mod)/sum(mod+unmod)`，须同时报告观测数；与 QC 位点均值不同，不混用为化学效率 |
| TAPS Non-CpG | 可判链宿主分子的 CHG/CHH 原始转换率；零观测缺失，未经 SNP/测序错误校正，不能直接解释为真实 mCH 或 cDNA 污染率；High-CpH 为 NaN |



| 默认 HQ 条件（同时满足） | 阈值 |
| --- | --- |
| Lambda methylation | `<10%` |
| pUC19 methylation | `>90%` |
| 宿主 unique CpG sites | `>500000` |
| Gini | `<0.5` |
| mapping | 默认不参与；对应数据完成验证后才设置 `MAPPING_THRESHOLD` |
| pair 去重率 | 默认不参与，Droplet 为 NaN |

| 图 / Notebook 控制 | 解读与输出 |
| --- | --- |
| Final QC | `FINAL_QC_MODE=per_clone_samples` 按 CloneID 分面，`distribution` 叠加；`FINAL_QC_METRIC` 支持单指标/列表，`FINAL_QC_CLONE_ID` 选组，逐项展示并保存 PDF |
| Final QC 横纵轴 | 100 bins、组间共用边界；CpG 正值用 log 横轴，含非正值改 signed-log1p 并标明，其他指标线性；纵轴 Counts，保留有效 N/中位数，实际直方图加 sigma=1.5、半径6、对称边界的高斯平滑，在变换空间插值；平滑保总计数且不参与 QC 筛选，仅排缺失/非数值/无穷 |
| `CpG_Density_Unique_vs_Input_Pairs.pdf` | 横轴 trimming 前 input pairs，纵轴最终宿主 unique CpG sites，CloneID 配色；散点用原始计数，不作深度归一化或异常点剔除 |
| 拟合与 `CpG_Density.pdf` | 散点拟合线与斜率汇总共用联合有效样本、异常点筛选及 `min_slope_ratio_change`；过原点 `slope=sum(x*y)/sum(x*x)`，单位 sites/input pair，不能拟合则 N/A |
| density 输入 | QC 身份/sites、sealed final manifest pairs、composition cache host signals 的样本集合须完全一致，不静默取交集 |
| 其他图 | `WGBS_QC_Panel.pdf`、`Plate_View_Native_Mapping%.pdf`、`CpG_Signal_Composition.pdf`、可选 `TSS_Profile.pdf`、`Batch_Effect_PCA.pdf` |

## 关联 RNA_DARLIN 结果

| 孔板 / 协议 | 操作与校验 |
| --- | --- |
| Cabernet / TAPS | DNA manifest 的可选 `rna_sample` 填 RNA_DARLIN **原始 RNA basename**；两边 `sample_id` 可独立命名，Alopex 只处理本地 DNA |
| 执行顺序 | 完成 DNA 发布及同代 downstream，再在 RNA Notebook 填 `DNA_PROJECT_ROOT` |
| 对齐来源 | 当前 complete DNA 的 sealed manifest + 同代 DNA QC；RNA.h5ad 的 `rna_sample → canonical sample_id → PlateID` 一对一关联，不读取未发布 config 或挑最新 mtime |
| 失败边界 | 参与细胞集合须相等；多对一、缺失/多余细胞、barcode/Cell_Order 或身份冲突直接失败，不静默取交集 |
| SRD / 无外部关联 | 用共享项目 sample_id 与 PlateID；SRD 本地 `rna_sample` 不参与外部桥接 |
| RNA 输出 | `Integrated_QC_Information.csv` 整合 RNA 13列与 DNA 20列并保留 RNA 身份；实际映射和输入哈希写 `analysis_summary.json`，DNA QC 不增加关联列 |
| 修改关联 | 先重发布 DNA 和当前代际 downstream，再重新整合 RNA |

| Droplet 独立关联 | 操作与校验 |
| --- | --- |
| 文库映射 | DNA `rna_sample` 留空；RNA 项目 `00_config/dna_rna_library_map.tsv` 指定 DNA/RNA 文库，按完整 CB17 关联 |
| 细胞集合 | 两边独立 calling/QC；`DNA_RNA_cell_links.csv` 保留双方全集，整合 QC 只含共同且双方 QC 可用的细胞 |
| 身份冲突 | 文库冲突或错误交付代际失败；不套用孔板的细胞集合相等约束，未配对不补科学零值 |
