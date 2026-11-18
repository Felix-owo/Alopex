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

#### 第一格参数：填写方式

**`DNA_PROJECT_ROOT`**

建议填 `"/path/to/Patient001"`；也接受 `Path` 对象，但须先导入 `Path`

**`DNA_PIPELINE_ROOT`**

建议填 `"/path/to/Alopex"`；显式指定仍须通过交付与 Pipeline 身份验证

**两根留空**

先检查同名环境变量；项目从 Jupyter 当前目录向上查找 complete delivery，Pipeline 从 sealed `source_root` 定位。显式路径错误不会静默回退

**`TSS_BED=None`**

自动使用当前物种的 TSS reference；外部 BED 须与 sealed reference 使用同一 genome build

**`SKIP_TSS=False`**

确实无需 TSS 时设为 `True`；不能同时指定外部 `TSS_BED`

首格只负责参数、源码导入和路径确认；后续计算与绘图分别调用 `qc`、`plotting` 模块。显式项目缺少 manifest 时立即报错，不借用当前目录的其他项目定位 Pipeline。自动定位时，发现当前项目的 `00_config/config.yaml` 但缺少交付记录，会在该项目停止。先完成它的上游发布，再运行 Notebook。权限错误、I/O 错误、悬空链接或非普通 manifest 也会指出出错位置。更换 Pipeline 或更新源码后须重启 kernel，从首格重跑。

### QC 运行面板与详细日志

`qc.run_qc_processor()` 自动显示一个持续更新的中文状态卡，包含当前阶段、开始时间/耗时、实际worker预算、子进程及本机进程树CPU/RSS采样。阶段与缓存明细、TSS逐CloneID抽样和诊断默认折叠；完整原始输出在项目`04_logs/downstream/qc.<时间>.<PID>.log`。已有调用参数不变，无需设置输出模式。

如果运行 QC 前上游发布了新的交付 ID，Processor 会在写结果前停止。遇到交付变更提示时，从 Notebook 第一格重新运行以读取新交付。

终端运行时自动选择原位动态面板，非交互输出仅保留阶段变化、告警和五分钟心跳。资源缺失显示不可用，RSS合计和采样峰值不代表独占物理内存或精确峰值；界面显示“最终提交”时仍须等待完成，只有子进程（含句柄关闭）成功退出才显示“已完成”。

失败先看面板原因，再打开完整日志定位首条错误。Notebook中断会先等待本次QC清理，再回收其进程组内仍存活的后代；等待中断处理结束后再重跑Processor。

关闭 Notebook 标签页不会停止 kernel 中的计算。需要停止 QC 时，使用中断操作并等待清理结束；可在 JupyterLab 的 Running 面板确认 kernel 状态，见[官方说明](https://jupyterlab.readthedocs.io/en/stable/user/running.html)。
日志目录不可用时会提示并显示原生输出；中断仍会等待本次计算清理并回收其进程组。先修复提示中的日志路径，再重跑即可恢复状态卡和完整调用日志。

### CLI

可在终端运行 Processor，之后打开 Notebook 绘图。CLI 同时生成 `SingleCpG_Adata.h5ad`。

```bash
PIPELINE=/path/to/Alopex
PATH="$PIPELINE/conda/current/control/bin:$PATH" \
  PYTHONPATH="$PIPELINE/core" \
  "$PIPELINE/conda/current/notebook/bin/python" -m alopex.downstream.qc \
  --project-dir /path/to/Patient001
```

`--workers 0` 自动选择，最多32；显式值也受实际 CPU、affinity、Slurm allocation 与调用方预算限制，每个 worker 的 native 线程为1。完整选项可用上述命令的 `--help` 查看。

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
| `DNAme_QC_Information.csv` | 下游与 adata 合并的唯一 QC 表，固定20列；[字段、顺序与单位][qc-fields] |
| `DNA_Doublet_1Mb_Counts.npz` | Droplet 同代 chr1–chr22 窗口整数 mc/cov，供 Lemmus 读取 |
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

不修改公开 CpG、文件名、表头、坐标或整数计数，也不修改 MethylTree 包源码。sample/clone/plate/barcode/cell order 来自 sealed manifest；原20列QC保留，别名和单位见[指标与分母](#指标与分母)。无观测保持NaN；signal比例是coverage比例，不是FASTQ reads比例。不要把原生mapping跨backend当成同一指标。

读图时，nCG 是有覆盖的宿主主核胞嘧啶位点数。分组图和汇总图均计入 nCG 为 0 的细胞，有零值时使用线性横轴。染色体散点图和分布图共用同一组 Pearson 系数。每对细胞仅使用共同有覆盖的常染色体，少于两条或任一侧在共同集合上无变异时，相关性记为缺失。分布图标题注明可计算配对数与总配对数。

矩阵沿用 SingleCpG 的坐标和双链 fraction/mean，原始逐位点整数计数保存在公开四列 CpG；导出按连续单行读取并直接复制 CSR 的 data/indices，保留显式零和缺失掩码。**稀疏未存储项代表缺失**，不能将 `toarray()` 后的零当成真实零甲基化。

每次导出重建并原子覆盖独立 H5AD，不复用用户已写入的谱系结果；重跑前关闭同一 Notebook 打开的 backing file。标准 SingleCpG/QC 仍按原有身份复用。临时H5AD验证和关闭前不替换旧目标；替换前撤销旧`allCpGs.input_schema.json`，两份输出成功后最后提交新完成身份。失败只清本次临时项；缺少完成身份时重跑入口，不接纳部分导出。后续MethylTree分析修改导出文件后不再符合该初始交付身份，衔接入口仍会重建。完整代码调用为：

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

参考身份以当前封存回执为准。

读图时先区分位点均值、coverage 比例和 read-pair 指标。身份列按字符串保留前导零及合法 `NA/nan/None`，数值缺失保持 NaN；观测零值不能当作缺失。完整定义见[QC 字段][qc-fields]和[生物学分母][denominators]。

| 指标 | 解读要点 |
| --- | --- |
| QC 甲基化率 | 宿主、lambda、pUC19 分别对已覆盖位点的 `m/(m+u)` 取均值再×100；无观测为 NaN |
| SingleCpG | 区间内已覆盖位点的 fraction 均值；`1/4` 与 `9/12` 的均值为0.5。有覆盖全 unmethyl 保显式0，缺覆盖不补0 |
| TSS | 每个 offset 对有覆盖 bins 的 fraction 取均值；count 为甲基化信号之和 |
| Gini | 仅主核染色体、按有 coverage 的固定宽 bins 计算，含 chrY、排 chrM 及辅助序列；无主核覆盖为缺失 |
| 信号组成 | 各类 canonical CpG coverage / 全部 contig 总 coverage；其他 decoy 仍进入分母。controls 仅精确匹配 `lambda` / `pUC19` |
| mapping / pair 指标 | mapping 的 backend 单位与分母不同；完整 pair funnel、policy/unit 保留在 sealed final manifest 与 MultiQC。Droplet 的 pair 去重率为 NaN，不能用位点 UMI 数代替 pair 数 |
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
| reference 顺序 | 保留原生 `reference_sequences` 与索引顺序，不按 FAI 重排；坐标与原生索引约束见[矩阵契约][single-contract] |
| 后续注释 | 只读对象不能回填标签；使用独立可写副本或合适大小的内存子集。MethylTree 分箱后仍按细胞 ID 对齐 `obs` |

固定计算参数见[single-CpG 参数][single-params]；矩阵的 `uns` 记录实际坐标 base、区间约定与来源表示。

## TSS profile

默认 reference 为 `resources/<genome>_reference/tss/<genome>_TSS_2000_2000_20.bed`，使用链方向一致的 ±2000 bp / 20 bp bins，第9列为 offset。

外部 BED 须匹配 sealed genome build，并通过完整坐标、重复区间和 offset 校验。reference 由 Doctor 根据已准备的 FASTA/GTF 构建，准备步骤见[主 README](../README.md)。

自备 TSS BED 时，用制表符分列。第 7 列是原样读取的基因 ID，第 9 列是有符号 32 位整数 offset。坐标或 offset 带小数时，只有恰好等于整数才有效。同一基因的两个区间若在染色体末端裁成同一区间，需要修正源 BED 后重新封存。

Notebook 第二格的 `TSS_SAMPLE_FRACTION`、`TSS_SAMPLE_SEED` 控制抽检，对应 CLI `--tss-sample-fraction` / `--tss-sample-seed`。默认每个原始 CloneID 抽10%：不足10个全检，其余向上取整且最多100个；比例设为1.0时全检。

选样按种子与 Sample_ID 确定，不受行序影响，也不按 coverage 选样。完整范围与记录位置见[TSS 抽检控制][recompute]。

抽检只影响 TSS，其他 QC 和矩阵仍保留全部细胞。CSV 仅记录选中且有覆盖的 bins，实际出现的细胞数可小于抽检数；有 coverage 的零甲基化仍保留，无覆盖不补0。全体无覆盖时仅表头 CSV 合法，旧 PDF 会清理。显式 `--skip-tss` / `SKIP_TSS=True` 会清理该交付旧 TSS CSV/PDF。

重新启用 TSS 时，先设 `SKIP_TSS=False` 并运行 QC Processor，再运行绘图单元。若本次应有的 CSV 丢失，也从 Processor 恢复。单独重跑绘图只使用本次完成记录绑定的曲线。

## 缓存与强制重算

默认按 sealed 输入、参数与算法身份复用缓存，无需先删除输出目录。输入身份失配或强制重算时，逐细胞解码一次，重新生成全部同代产物。重算开关默认均为 `False`，只影响下游。

| Notebook 开关 | CLI 对应项 | 重建对象 |
| --- | --- | --- |
| `RECOMPUTE_CPG_STATS` | `--recompute-cpg-stats` | 同代 QC、SingleCpG 和适用的 1 Mb 计数 |
| `RECOMPUTE_GINI` | `--recompute-gini` | 同上 |
| `RECOMPUTE_TSS` | `--recompute-tss` | 同上，按当前 TSS 设置生成诊断 |

| 遇到的情况 | 处理方式 |
| --- | --- |
| 调整 TSS BED、比例或种子 | Processor 自动重建整代产物，无需打开强制重算开关 |
| 只调整 HQ 门槛 | 重新运行 Processor；已完成矩阵复用，只更新 `obs`、筛选参数和筛选表 |
| QC 完成身份缺失、过期或重算中断 | 重新运行 Processor，完成后再加载 QC 或读取 SingleCpG |
| SingleCpG 来源未知或变化 | 重新运行 Processor，从 sealed CpG 自动恢复 |
| 输入或参考在计算期间变化 | 确认文件稳定后重跑；失败时不提交完成身份，未完成矩阵不能作为有效缓存使用 |
| CpG 格式错误或缓存计数矛盾 | 按错误定位 canonical CpG 输入并重算；不把格式失败当作合法空细胞，也不手改缓存绕过校验；输入格式与守恒见[下游契约][qc-contract] |

模块更新后，结束旧写任务并重启 kernel。最终 QC 提交中断后须重跑 Processor 恢复，再继续消费结果。缓存与提交的完整规则见[下游契约][qc-contract]；实测与验证范围见[BENCHMARKS][benchmarks]。

## 关联 RNA_DARLIN 结果

Cabernet、TAPS、Droplet 均在 DNA manifest 的可选 `rna_sample` 填 RNA_DARLIN **原始 RNA basename**，不是 RNA `sample_id` 或路径；两边项目样本可独立命名，Alopex 只处理本地 DNA。

先完成 DNA 发布和同代 downstream QC，再在 RNA Notebook 填 `DNA_PROJECT_ROOT`。RNA 下游只读当前 complete DNA 的 sealed manifest 与同代 `DNAme_QC_Information.csv`，不读取可变 config 或按 mtime 挑目录。修改关联后须重发布 DNA、生成同代 downstream QC，再重新整合 RNA。

| 协议 | 文库与细胞关联 |
| --- | --- |
| Cabernet / TAPS | sealed `rna_sample` → finalized RNA.h5ad 原始 `rna_sample` → RNA canonical `sample_id`，再按 PlateID 一对一连接；未填外部关联时使用共享项目 `sample_id` 与 PlateID |
| Droplet | 同样从 sealed `rna_sample` 解析文库，再按对应文库内完整 CB17 精确连接；留空保持未关联，不自动猜配 |

孔板协议要求参与 QC 的两侧细胞集合完全一致，barcode/Cell_Order 冲突、多对一或非身份列重名均失败，不静默取交集。Droplet 两侧独立 calling/QC，`DNA_RNA_cell_links.csv` 保留双方全集，整合 QC 只含共同且双方 QC 可用的细胞；未配对细胞不补科学零值。文库冲突或错误交付代际均拒绝。

`Integrated_QC_Information.csv` 保留 RNA 身份，RNA 13列与 DNA 20列的独立 QC schema 不变；RNA 下游记录实际关联和输入哈希。完整关联边界见[Manifest 契约][manifest-contract]与[Droplet 关联契约][droplet-contract]。

[qc-fields]: #指标与分母
[denominators]: #指标与分母
[recompute]: #缓存与强制重算
[single-params]: #保存供下游复用的-single-cpg-矩阵
[single-contract]: #保存供下游复用的-single-cpg-矩阵
[qc-contract]: #指标与分母
[manifest-contract]: #关联-rna_darlin-结果
[droplet-contract]: #关联-rna_darlin-结果
[benchmarks]: #指标与分母

## 诊断与矩阵注释

Notebook提供逐染色体覆盖/甲基化图与逐细胞筛选原因表；TSS QC沿用原有profile及可视化。默认HQ规则保持不变。

末尾Single矩阵保存完整obs注释、HQ/筛选原因、sealed条码与孔序，后续可直接使用adata.obs。全部细胞仍保留，选择HQ由下游显式执行。标准H5AD/CSR路径不需要SnapATAC2或Polars；大参考矩阵用磁盘CSR逐行处理，避免一次读取全体var标签或obsm。

更新后重启kernel，从第一格运行。Processor按各自身份更新失效的下游缓存；本次调整只使TSS缓存失效，Raw、CpG、Gini、Single-CpG与染色体统计的算法身份不变。原TSS profile保留`Sample_ID/position_coarse/meth_frac/count/CloneID`五列。缓存与输入失败时保留原始日志，修复原因后运行Processor，再运行Notebook后续单元。

## 指标与分母

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
| `Gini_Index` | canonical CpG coverage在有覆盖固定宽bin的Gini；限主核染色体，含chrY、排除chrM及辅助序列；主核位点数为0时为NaN，即使线粒体、辅助序列或对照有覆盖；真实零甲基化且有主核coverage仍计算Gini |
| `Duplicate_Pair_Rate%`、`Trim_Retention%`、`Final_Pair_Yield%` | 分母按本节pair表；Droplet不以分子数反推pair去重率，该列为NaN |


| 指标 | 定义 |
| --- | --- |
| QC三参考Meth_CpG_Rate% | 分别对宿主/lambda/pUC19已覆盖的canonical胞嘧啶位点取m/(m+u)的算术均值×100，不是将reads计数合并后的比率；零位点为NaN |
| contig归属 | 位点统计与coverage分类共用原始名称：controls仅精确`lambda`/`pUC19`，mt按既有chrM/MT/chrMT规则，宿主仅主核染色体（标准人chr1–22/X/Y、小鼠chr1–19/X/Y）；其他contig只进total分母。例如`chr_lambda_alt`、`lambda_decoy`均不属宿主或controls；不按子串或大小写别名分类 |
| SingleCpG X | 同一参考CpG区间内已覆盖C/G位点的fraction均值；例如1/4与9/12得到0.5，不是10/16。缺覆盖为稀疏缺项，有覆盖但全unmethyl为显式0；不把缺项补成观测0 |
| TSS | 同offset有覆盖gene/bin的fraction均值；count为这些gene/bin甲基化信号之和，不能当总coverage/read pairs或全基因组唯一分子数。共享基因组区间仍保留每个gene关联/权重；无覆盖不出行，有覆盖全unmethyl保fraction=0/count=0 |
| TAPS native controls | lambda假阳性/pUC19转换取总mod/总(mod+unmod)，与上述QC位点均值不是同一分母；必须同时报告观测数，不混用为化学效率 |


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

MultiQC 的 native 模块仅去掉本 Pipeline 的末尾报告后缀，不使用默认中间子串截断或首尾字符裁剪；Cutadapt 固定按 `<Sample_ID>.cutadapt.json`、BISCUIT 按唯一细胞报告文件名取身份，重复清名也只去一次报告后缀；文件格式只判最后扩展名。

保留前导连字符、`NA`、前导零及 ID 内的点号/`.cutadapt`，以及与报告后缀同名的结尾。原始 native JSON/文本不改写。

MultiQC 只对展示抽样，文库概要与科学结果仍覆盖全库；完整名单与诊断保留在对应原始日志和已封存交付中。

Cutadapt adapter QC 分开输出 `R1_Adapter_pct` 与 `R2_Adapter_pct`。

`non_cpg_methylation_pct` 对 BISCUIT 来自 CpH retention table，对 Cabernet Bismark 来自原生 extraction splitting report 的 CHG+CHH 观察，对 Droplet 来自过滤后 CHG/CHH 的胞嘧啶/UMI 共识计数，对 TAPS 采用[判链宿主统计](#指标与分母)。

均记录 `non_cpg_metric_source`，不作为跨 backend 的统一 hard threshold。
