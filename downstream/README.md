# Alopex 下游 QC 与可视化

版本以当前 checkout 或正式发行包对应的 Git tag 为准。本文说明下游运行、结果读取、缓存恢复与 RNA 关联；安装和上游流程见[主 README](../README.md)。

## 从这里开始

1. **选对 kernel。** 打开 [DNA_QC_Visualization.ipynb](DNA_QC_Visualization.ipynb)，使用 `Alopex QC`；更新源码或更换 Pipeline 后先重启 kernel。
2. **确认第一格的两条路径。** `DNA_PROJECT_ROOT` 指向数据项目，`DNA_PIPELINE_ROOT` 指向 Pipeline；按[运行前检查](#运行)确认 complete 交付可读。
3. **按顺序运行。** 等 `qc.run_qc_processor()` 状态卡变为“已完成”后再加载 QC、绘图；矩阵与细胞注释可直接使用，见[SingleCpG 读取方式](#保存供下游复用的-single-cpg-矩阵)。

| 现在要做什么 | 直接进入 |
| --- | --- |
| 第一次注册 kernel，或只用终端计算 | [kernel 注册](#jupyter--vs-code-kernel) / [CLI](#cli) |
| 找 QC 表、SingleCpG | [输出](#输出) |
| 恢复失败计算、调整重算范围 | [缓存与强制重算](#缓存与强制重算) |
| 调整图、筛选或 TSS 抽检 | [指标与筛选](#指标与筛选) / [TSS](#tss-profile) |
| 连接 MethylTree 或 RNA_DARLIN | [MethylTree](#连接-methyltree) / [RNA 关联](#关联-rna_darlin-结果) |

**完成判据：** Processor 成功退出后才显示“已完成”；结果在项目 `06_downstream/`。Processor 同时生成 QC 与带细胞注释的 SingleCpG。

## 运行

运行前确认：

- 项目已有当前 complete 的 `03_results/run_manifest.json`；下游从交付记录读取 sealed manifest、CpG inventory 和配置快照。
- 交付绑定的 snapshot/config、reference/FAI 及记录的项目/Pipeline 路径可读。只复制 `03_results/` 不能完成迁移；下游不依赖 `02_work/` 的 BAM 或 demux FASTQ。
- 使用 Doctor 管理的 `conda/current/notebook` 环境。下游写入 `06_downstream/`；上游重发布后须重新启动下游。

### Notebook

打开 [DNA_QC_Visualization.ipynb](DNA_QC_Visualization.ipynb)，选择 `Alopex QC` kernel，确认第一格参数后按顺序运行。

Notebook 调用独立 Processor 完成计算，再加载 QC 绘图；Processor 子进程沿用 kernel 的 Python，并优先使用同一 Pipeline 的 `conda/current/control/bin` 工具。

| 第一格参数 | 填法与含义 |
| --- | --- |
| `DNA_PROJECT_ROOT` | 项目路径，如 `"/path/to/Patient001"`。使用 `Path` 时先导入 |
| `DNA_PIPELINE_ROOT` | 安装路径，如 `"/path/to/Alopex"`，须与交付记录吻合 |
| 两根留空 | 先查同名环境变量，再从当前目录向上找最近的项目，从交付的 `source_root` 取得 Pipeline |
| `TSS_BED=None` | 使用该物种的默认 TSS 参考。外部 BED 须与封存参考匹配 |
| `SKIP_TSS=False` | 跳过 TSS 时设为 `True` 并保留 `TSS_BED=None` |

遇到显式路径错误、项目身份不符、缺交付或读取失败时，按报错修正目标项目。已发现项目配置时先完成该项目的上游发布。更换 Pipeline、更新源码或上游重新发布后，重启 kernel，从第一格重跑。

若 BISCUIT 旧交付的 `native_mapping_unit` 仍是 `individual_primary_reads_postdedup`，请在数据项目目录用当前 Pipeline 的 `bash /path/to/Alopex/core/run_pipeline.sh` 重新交付，再重启 kernel。新的清单按去重后的非 secondary 比对记录计算原生比对率，含仅带 supplementary 标志的记录。

### QC 运行面板与详细日志

Notebook 的状态卡展示阶段、开始时间、耗时、CPU 预算及本机进程树采样。展开详情可查看缓存、TSS 抽样与告警。终端在原位更新，重定向输出只记阶段、告警、结束和每五分钟心跳。RSS 合计含共享页，峰值是采样值。子进程关闭句柄并成功退出后才显示完成。

失败先查 `04_logs/downstream/qc.<时间>.<PID>.log`。日志目录不可用时查看原生输出，修复目录后再运行。中断请使用 kernel 的中断操作并等待清理，关闭标签页会保留计算进程。交付 ID 变化时回到第一格。[JupyterLab 运行中的 kernel](https://jupyterlab.readthedocs.io/en/stable/user/running.html)

### CLI

可在终端运行 Processor，之后打开 Notebook 绘图。CLI 同时生成 `SingleCpG_Adata.h5ad`。

```bash
PIPELINE=/path/to/Alopex
PATH="$PIPELINE/conda/current/control/bin:$PATH" \
  PYTHONPATH="$PIPELINE/core" \
  "$PIPELINE/conda/current/notebook/bin/python" -m alopex.downstream.qc \
  --project-dir /path/to/Patient001
```

`--workers 0` 自动采用可用 CPU、affinity 和 Slurm allocation 允许的预算。第 3 阶段并行解码细胞、计算统计量和 SingleCpG 矩阵行，由一个线程按原序写盘。面板显示整个进程树的实测 CPU，用量会随参考准备、细胞覆盖深度和写盘速度变化。

### Jupyter / VS Code kernel

首次使用时注册名为 `alopex-qc`、显示为 `Alopex QC` 的 kernel。启动路径固定到 `conda/current`，环境重建后无需重新注册；Pipeline 迁移目录后须重新注册并选择 kernel。

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

## 输出

结果固定写入 `<project>/06_downstream/`。完整交付/输入身份匹配才复用缓存；新结果覆盖对应旧文件，不创建时间戳或哈希代际目录。完成身份最后提交，失败时不能把新旧混合文件当作完整结果；历史目录不自动搬迁或删除。

| 文件 / 目录 | 用途 |
| --- | --- |
| `DNAme_QC_Information.csv` | 下游与 adata 合并的唯一 QC 表，固定20列 |
| `DNA_Doublet_1Mb_Counts.npz` | Droplet 常染色体 1 Mb 窗口整数 mc/cov，hg38 为 chr1–22，mm10 为 chr1–19，供 Lemmus 读取 |
| `SingleCpG_Adata.h5ad` / `SingleCpG_Adata.input_schema.json` | Processor 单次解码生成的全细胞、主核染色体 single-CpG 矩阵与完成记录 |
| `Chromosome_QC_Information.csv` | 全细胞逐染色体整数计数、位点数与甲基化率 |
| `QC_Filter_Reasons.csv` | Processor 生成的逐细胞筛选结果与失败原因 |
| `TSS_Profile_Information.csv` | 原有 TSS profile：抽中细胞有覆盖 offset 的甲基化率与甲基化信号整数计数 |
| `plots/` | QC、Final QC、孔板、信号组成、TSS、CpG density 与 PCA 等 PDF |
| `.cache/` | 信号组成及 QC 完成身份等内部缓存，不手动改写 |
| `methyltree/allCpGs.h5ad` / `methyltree/sample_sheet.csv` | 显式运行 MethylTree 衔接入口时生成；独立的所选细胞矩阵及全部细胞注释 |

同一项目只能有一个下游写入者。遇到 `Downstream QC is busy`，等当前 Processor 或 SingleCpG 更新结束后重试。项目和 QC 目录内的 `.qc-write.lock` 是持久锁文件，**不要删除或替换**；进程退出后系统会释放锁，文件仍存在不代表任务在运行。

## 连接 MethylTree

1. 使用独立 MethylTree kernel（模板标记为 `MethylTree_3_10_lab`，按本机实际注册项选择）打开 [only_quality_checking_alopex.ipynb](only_quality_checking_alopex.ipynb)。从目标 Alopex 项目或子目录启动 kernel。此模板按工作目录发现最近的项目，从其完整交付读取 Pipeline 路径。Alopex 计算通过该 Pipeline 的受管 Python 子进程执行，无需向 MethylTree kernel 安装 Alopex 的整套依赖；无需先执行 Alopex 可视化 Notebook。
2. Processor 从当前 complete 交付构建或复用标准 QC/SingleCpG；`alopex.downstream.methyltree.prepare_methyltree_inputs` 复用同一 SingleCpG builder，导出至固定 `06_downstream/methyltree/`。
3. 默认不运行可选谱系推断；该可选分析需在自己的 MethylTree 环境单独验收。完成后确认全部细胞在 `sample_sheet.csv`，按导入门槛选中的细胞在 `allCpGs.h5ad`。H5AD 为标准 AnnData，X 保持 backed。Notebook 保留 MethylTree 原有绘图、筛选、注释、可选谱系和诊断流程，需安装 MethylTree；仅直接调用下方 Alopex 导出 API 时不要求该包。

在 macOS 上，如果开启谱系后内核退出，且崩溃栈位于 `libomp`，先检查该 MethylTree 环境是否同时加载了多份 OpenMP 运行库。本次混用 pip 与 Conda 的内核加载了三份。在启动该内核前设置 `NUMBA_THREADING_LAYER=workqueue` 后，同一数据完成了推断和保存。[Numba 线程层文档](https://numba.readthedocs.io/en/stable/user/threading-layer.html)

不修改公开 CpG、文件名、表头、坐标或整数计数，也不修改 MethylTree 包源码。sample/clone/plate/barcode/cell order 来自 sealed manifest；原20列QC保留，无观测保持NaN；signal比例是coverage比例，不是FASTQ reads比例。不要把原生mapping跨backend当成同一指标。

读图时，nCG 是有覆盖的宿主主核胞嘧啶位点数。分组图和汇总图均计入 nCG 为 0 的细胞，有零值时使用线性横轴。染色体散点图和分布图共用同一组 Pearson 系数。每对细胞仅使用共同有覆盖的常染色体，少于两条或任一侧在共同集合上无变异时，相关性记为缺失。分布图标题注明可计算配对数与总配对数。

矩阵沿用 SingleCpG 的坐标和双链 fraction/mean，原始逐位点整数计数保存在公开四列 CpG；导出按连续单行读取并直接复制 CSR 的 data/indices，保留显式零和缺失掩码。**稀疏未存储项代表缺失**，不能将 `toarray()` 后的零当成真实零甲基化。

每次导出重建并原子覆盖独立 H5AD。重新运行输入单元会关闭该 kernel 上一次打开的矩阵。打开矩阵、执行可选谱系分析及写入样本表或矩阵时，Notebook 在下游写锁内核对本次导出。看到 `export changed` 时从输入单元重新运行。看到 `Downstream QC is busy` 时等当前写者完成，再运行该单元。标准 SingleCpG/QC 仍按原有身份复用。临时H5AD验证和关闭前不替换旧目标；替换前撤销旧`allCpGs.input_schema.json`，两份输出成功后最后提交新完成身份。失败只清本次临时项；缺少完成身份时重跑入口，不接纳部分导出。后续MethylTree分析修改导出文件后不再符合该初始交付身份，衔接入口仍会重建。完整代码调用为：

```python
from alopex.downstream import qc, methyltree

paths = qc.resolve_notebook_paths(DNA_PROJECT_ROOT, DNA_PIPELINE_ROOT)
qc.run_qc_processor(paths, skip_tss=True)
samples, h5ad = methyltree.prepare_methyltree_inputs(
    paths, nCG_cutoff=100000, import_nCG_cutoff=500000,
)
```

单独读取CSV时为 `Sample_ID/CloneID/PlateID/sample/cell/clone_id/plate_id/plate_barcode/cell_order` 显式使用字符串converter，保留前导零和合法`NA/nan`；数值列仍按数值读取，不能统一关闭缺失解析后把缺失当字符串。

## 指标与筛选

标准分析范围为 hg38 的 chr1–22/X/Y、mm10 的 chr1–19/X/Y；辅助／decoy 参与比对，不进入宿主 nCG、甲基化率、Gini、TSS 和最终 Single-CpG 矩阵。TSS／single-CpG BED 必须携带与当前参考匹配的 `.reference.json`，不得把旧目录套到新交付。

读图时先区分位点均值、coverage 比例和 read-pair 指标。身份列按字符串保留前导零及合法 `NA/nan/None`，数值缺失保持 NaN；观测零值不能当作缺失。

| 指标 | 解读要点 |
| --- | --- |
| QC 甲基化率 | 宿主、lambda、pUC19 分别对已覆盖位点的 `m/(m+u)` 取均值再×100；无观测为 NaN |
| SingleCpG | 区间内已覆盖位点的 fraction 均值；`1/4` 与 `9/12` 的均值为0.5。有覆盖全 unmethyl 保显式0，缺覆盖不补0 |
| TSS | CSV 按细胞和 offset 报告有覆盖 gene/bin 的平均甲基化比例（float64）及甲基化信号整数和。图按克隆平均该 offset 有覆盖的细胞，包含实测零。阴影为固定种子 bootstrap 的 95% 区间。 |
| Gini | 仅主核染色体、按有 coverage 的固定宽 bins 计算，含 chrY、排 chrM 及辅助序列；无主核覆盖为缺失 |
| 信号组成 | 各类 canonical CpG coverage / 全部 contig 总 coverage；其他 decoy 仍进入分母。controls 仅精确匹配 `lambda` / `pUC19` |
| mapping / pair 指标 | mapping 的 backend 单位与分母不同；完整 pair funnel、policy/unit 保留在 sealed final manifest 与 MultiQC。DD-MET5 的 UMI 和 DD-MET3 的片段标签均按位点做共识，两者的 pair 去重率为 NaN |
| TAPS | controls 使用合并计数比率；Non-CpG 为可判链宿主 CHG/CHH 转换率，不能直接解释为真实 mCH 或 cDNA 污染率，High-CpH 为 NaN |

### HQ 标记

Notebook 默认同时要求 lambda methylation `<10%`、pUC19 methylation `>90%`、宿主 unique CpG sites `>500000`、Gini `<0.5`。mapping 默认不参与；对应数据完成验证后才设置 `MAPPING_THRESHOLD`。pair 去重率默认不参与。

`HQ` 是细胞注释；SingleCpG 保留全部细胞，并保存标准 QC 与 `HQ`。在 Processor 单元设置筛选门槛。

### 图与展示控制

染色体诊断图沿用前面 QC 单元的主题、字号和 CloneID 公共前后缀省略规则；Notebook 传入同一份 `clone_color_map`，使各图的组别颜色一致，完整身份仍保留在表格和矩阵中。热图分别用青绿表示覆盖、玫红表示甲基化；甲基化缺失为灰色，观测零仍显示为色阶起点。下方表格按 CloneID 和染色体展示数值列中位数，逐细胞原表保留在 `df_chromosome` 中。

染色体图支持独立调用：从当前完成 QC 读取数据后，可运行 `plotting.plot_chromosome_diagnostics(..., clone_color_map=clone_color_map)` 并导出图片；不传颜色表时使用现有默认配色规则。

#### 图 / Notebook 控制：解读与输出

**孔板图**

在绘图单元用 `PLATE_VIEW_METRIC` 选择数值列，默认是 `Native_Mapping%`。每个 CloneID 单独成图面，所有图面共用数值范围。孔内标签最多保留四位有效数字，较大的数值以 K、M 缩写。精确值在 `DNAme_QC_Information.csv`。

有观测的零显示为 `0`，缺失和无穷值的孔留空。当前指标没有可绘制的孔时，重新运行该单元会清理对应的旧 PDF。指标名拼错或列中含有非数值记录时，单元报错。

**Final QC**

- `FINAL_QC_MODE=per_clone_samples`：每个指标一张图，按 CloneID 分面，适合逐组检查；
- `distribution`：每个指标一张图，各 CloneID 按固定颜色叠加，适合比较组间分布，组数多时降低填充透明度以减少遮挡。
- 两种模式均严格使用 `FINAL_QC_METRIC` 的单指标/有序列表，不自动替换成预设指标；
- `FINAL_QC_CLONE_ID=None` 选全部，也支持单个 CloneID 或列表

**Final QC 指标与保存**

- 支持 QC 表全部17个数值字段（含 `Lambda_Meth_CpG_Rate%`、`pUC19_Meth_CpG_Rate%`）；未知/缺列指标明确报错。
- Notebook 统一调用 batch，每项保存 `Final_QC_Dashboard_<mode>_<metric>_<clone_selection>.pdf`（`%` 转为 `pct`）；两模式文件独立。
- 仅展示与保存绘图，无需重跑 Processor；更新绘图库后须重启 kernel，重新运行首格和 QC 加载格再绘图

**Final QC 横纵轴**

- 两模式共用配色、100 bins、组间共用边界；CpG 位点数正值用 log 横轴，含非正值改 signed-log1p 并标明；log/log1p 主刻度定位于原始位点数的十进数量级，显示 `1E1、1E2、1E3…`，零值单独标 `0`；log1p 近零区间与纯 log 的间距略有不同。
- 宽范围适当跳过数量级以免重叠，窄范围内无整数数量级时标实际中点的 E 记数值。关闭次刻度标签和轴外科学计数法偏移；不把 `9999` 伪标为 `1E4`。横轴在变换空间两端各留4%空白；
- 纵轴按全部已绘组的最高原始柱计数留峰顶20%、零基线下3%空白，完整容纳端点柱、曲线与线宽；留白不改变分箱或计数。甲基化率等其他指标线性；纵轴 Counts（受各组样本量影响，不是概率密度）。
- 分面标注各组有效 N/输入 N 与中位数，叠加标注全部所选组的有效 N/输入 N 与总体中位数（虚线）；无有限值明确标注 No finite data，不填0、不切换指标。实际直方图加 sigma=1.5、半径6、对称边界的高斯平滑，在变换空间插值；平滑保总计数且不参与 QC 筛选，仅排缺失/非数值/无穷

**`CpG_Density_Unique_vs_Input_Pairs.pdf`**

横轴 trimming 前 input pairs，纵轴最终宿主 unique CpG sites，CloneID 配色；散点用原始计数，不作深度归一化或异常点剔除

**拟合与 `CpG_Density.pdf`**

散点拟合线与斜率汇总共用联合有效样本、异常点筛选及 `min_slope_ratio_change`；过原点 `slope=sum(x*y)/sum(x*x)`，单位 sites/input pair，不能拟合则 N/A

**density 输入**

QC 身份/sites、sealed final manifest pairs、composition cache host signals 的样本集合须完全一致，不静默取交集

**其他图**

`WGBS_QC_Panel.pdf`、`Plate_View_Native_Mapping%.pdf`、`CpG_Signal_Composition.pdf`、可选 `TSS_Profile.pdf`、`Batch_Effect_PCA.pdf`

**PCA**

使用所选QC特征均有效的细胞，保留真实0值，缺失甲基化率不补成0%；图注列出纳入/输入细胞数及排除数。空输入、特征不足、少于2个完整细胞或无变异则跳过并清掉明确指定位置的旧PCA图；QC表及其他图的细胞集合不受影响

## 保存供下游复用的 single-CpG 矩阵

Processor 按 `Sample_ID` 对齐全细胞 QC 和筛选结果并保存到矩阵。Processor 成功结束后，使用该单元打印的 `SINGLECPG_ADATA_PATH` 读取。

```python
import anndata as ad

adata = ad.read_h5ad(SINGLECPG_ADATA_PATH, backed="r")
# 在此使用只读对象；结束后关闭句柄。
adata.file.close()
```

| 使用边界 | 要求 |
| --- | --- |
| 完成记录 | 矩阵与 `SingleCpG_Adata.input_schema.json` 配套使用；输入身份失配时由 Processor 从当前 sealed CpG 自动重建，见[缓存与强制重算](#缓存与强制重算) |
| 坐标 | 跟随 sealed `analysis.cpg_coordinate_base`：0-based `chr1:0-2` 为半开区间，1-based `chr1:1-2` 为闭区间，表示同一 CpG；传给 BED/0-based 工具前须按实际 metadata 转换 |
| 位点列序 | 按 SingleCpG BED 的原始行序排列。跨矩阵对齐时使用 `var_names` 中的区间和 `uns` 中的坐标基准。 |
| 后续注释 | 只读对象不能回填标签；使用独立可写副本或合适大小的内存子集。MethylTree 分箱后仍按细胞 ID 对齐 `obs` |

矩阵的 `uns` 记录实际坐标 base、区间约定与来源表示。

## TSS profile

默认 reference 为 `resources/<genome>_reference/tss/<genome>_TSS_2000_2000_20.bed`，使用链方向一致的 ±2000 bp / 20 bp bins，第9列为 offset。

外部 BED 须匹配 sealed genome build，并通过完整坐标、重复区间和 offset 校验。reference 由 Doctor 根据已准备的 FASTA/GTF 构建，准备步骤见[主 README](../README.md)。

自备 TSS BED 时，用制表符分列。第 7 列是基因 ID，第 9 列是有符号 32 位整数 offset。坐标和 offset 必须恰为整数。每个基因在同一 offset 只能有一个区间，不同基因可以共享坐标。遇到 NUL 字节、重复 offset 或裁切后重复的基因区间时，修正源 BED 后重新封存。

Notebook 第二格的 `TSS_SAMPLE_FRACTION`、`TSS_SAMPLE_SEED` 控制抽检，对应 CLI `--tss-sample-fraction` / `--tss-sample-seed`。默认每个原始 CloneID 抽10%：不足10个全检，其余向上取整且最多100个；比例设为1.0时全检。

选样按种子与 Sample_ID 确定，不受行序影响，也不按 coverage 选样。

抽检只影响 TSS，其他 QC 和矩阵仍保留全部细胞。CSV 仅记录选中且有覆盖的 bins，实际出现的细胞数可小于抽检数；有 coverage 的零甲基化仍保留，无覆盖不补0。全体无覆盖时仅表头 CSV 合法，旧 PDF 会清理。显式 `--skip-tss` / `SKIP_TSS=True` 会清理该交付旧 TSS CSV/PDF。

重新启用 TSS 时，先设 `SKIP_TSS=False` 并运行 QC Processor，再运行绘图单元。若本次应有的 CSV 丢失，也从 Processor 恢复。单独重跑绘图只使用本次完成记录绑定的曲线。

## 缓存与强制重算

默认核对封存输入、参数、算法和全部同代工件。强制重算使用 `RECOMPUTE_CPG_STATS=True` 或 `--recompute-cpg-stats`，一次解码生成 QC、SingleCpG 及适用的 1 Mb 计数。孔板项目重算或复用时清理旧 Droplet 计数文件。

| 遇到的情况 | 处理方式 |
| --- | --- |
| 调整 TSS BED、比例或种子 | Processor 自动重建整代产物，无需打开强制重算开关 |
| 只调整 HQ 门槛 | 重新运行 Processor；已完成矩阵复用，只更新 `obs`、筛选参数和筛选表 |
| QC 完成身份缺失、过期或重算中断 | 重新运行 Processor，完成后再加载 QC 或读取 SingleCpG |
| SingleCpG 来源未知或变化 | 重新运行 Processor，从 sealed CpG 自动恢复 |
| 输入或参考在计算期间变化 | 确认文件稳定后重跑；失败时不提交完成身份，未完成矩阵不能作为有效缓存使用 |
| 已发布 CpG 的 SHA-256 不符 | 按报错路径检查公开文件及其暂存原件。从完好的暂存重新发布，或先恢复有校验值的备份、重跑该产物的生成规则。交付恢复后从 Notebook 首格重跑 |
| CpG 格式错误或缓存计数矛盾 | 按错误定位 canonical CpG 输入并重算；不把格式失败当作合法空细胞，也不手改缓存绕过校验 |

升级后，算法身份变化会触发一次完整重建。完成身份缺失时重跑 Processor，确认提交成功后再消费结果。

## 关联 RNA_DARLIN 结果

Cabernet、TAPS、Droplet 均在 DNA manifest 的可选 `rna_sample` 填 RNA_DARLIN **原始 RNA basename**，不是 RNA `sample_id` 或路径；两边项目样本可独立命名，Alopex 只处理本地 DNA。

先完成 DNA 发布和同代 downstream QC，再在 RNA Notebook 填 `DNA_PROJECT_ROOT`。RNA 下游只读当前 complete DNA 的 sealed manifest 与同代 `DNAme_QC_Information.csv`，不读取可变 config 或按 mtime 挑目录。修改关联后须重发布 DNA、生成同代 downstream QC，再重新整合 RNA。

| 协议 | 文库与细胞关联 |
| --- | --- |
| Cabernet / TAPS | sealed `rna_sample` → finalized RNA.h5ad 原始 `rna_sample` → RNA canonical `sample_id`，再按 PlateID 一对一连接；未填外部关联时使用共享项目 `sample_id` 与 PlateID |
| Droplet | 同样从 sealed `rna_sample` 解析文库，再按协议连接：DD-MET5 精确匹配同一 CB17，DD-MET3 精确匹配官方 gex_cb 与 m_cb；留空保持未关联，不自动猜配 |

孔板协议要求参与 QC 的两侧细胞集合完全一致，barcode/Cell_Order 冲突、多对一或非身份列重名均失败，不静默取交集。Droplet 两侧独立 calling/QC，`DNA_RNA_cell_links.csv` 保留双方全集，整合 QC 只含共同且双方 QC 可用的细胞；未配对细胞不补科学零值。文库冲突或错误交付代际均拒绝。

`Integrated_QC_Information.csv` 保留 RNA 身份，RNA 13列与 DNA 20列的独立 QC schema 不变；RNA 下游记录实际关联和输入哈希。


## 诊断与矩阵注释

Notebook提供逐染色体覆盖/甲基化图与逐细胞筛选原因表；TSS QC沿用原有profile及可视化。默认HQ规则保持不变。

Processor 同时在 SingleCpG 的 obs 保存 HQ、筛选原因、sealed 条码与孔序，后续可直接使用adata.obs。全部细胞仍保留，选择HQ由下游显式执行。标准H5AD/CSR路径不需要SnapATAC2或Polars；大参考矩阵用磁盘CSR逐行处理，避免一次读取全体var标签。

TSS profile 有 `Sample_ID`、`position_coarse`、`meth_frac`、`count`、`CloneID` 五列。更新源码后重启 kernel，从第一格运行；缓存重建与阈值变化分别按[恢复表](#缓存与强制重算)处理。
