# Alopex v16.1 下游 QC 与可视化

正式源码以 Git tag `v16.1` 为准。本文说明下游运行、结果读取、缓存恢复与 RNA 关联；安装和上游流程见[主 README](../README.md)。

[运行](#运行) · [输出](#输出) · [指标与筛选](#指标与筛选) · [SingleCpG 矩阵](#保存供下游复用的-single-cpg-矩阵) · [TSS](#tss-profile) · [缓存与重算](#缓存与强制重算) · [RNA 关联](#关联-rna_darlin-结果)

## 运行

运行前确认：

- 项目已有当前 complete 的 `03_results/run_manifest.json`；下游从交付记录读取 sealed manifest、CpG inventory 和配置快照。
- 交付绑定的 snapshot/config、reference/FAI 及记录的项目/Pipeline 路径可读。只复制 `03_results/` 不能完成迁移；下游不依赖 `02_work/` 的 BAM 或 demux FASTQ。
- 使用 Doctor 管理的 `conda/current/notebook` 环境。下游写入 `06_downstream/`；上游重发布后须重新启动下游。

### Notebook

打开 [DNA_QC_Visualization.ipynb](DNA_QC_Visualization.ipynb)，选择 `Alopex QC` kernel，确认第一格参数后按顺序运行。Notebook 调用独立 Processor 完成计算，再加载 QC 绘图；Processor 子进程沿用 kernel 的 Python，并优先使用同一 Pipeline 的 `conda/current/control/bin` 工具。

| 第一格参数 | 填写方式 |
| --- | --- |
| `DNA_PROJECT_ROOT` | 建议填 `"/path/to/Patient001"`；也接受 `Path` 对象，但须先导入 `Path` |
| `DNA_PIPELINE_ROOT` | 建议填 `"/path/to/Alopex"`；显式指定仍须通过交付与 Pipeline 身份验证 |
| 两根留空 | 先检查同名环境变量；项目从 Jupyter 当前目录向上查找 complete delivery，Pipeline 从 sealed `source_root` 定位。显式路径错误不会静默回退 |
| `TSS_BED=None` | 自动使用当前物种的 TSS reference；外部 BED 须与 sealed reference 使用同一 genome build |
| `SKIP_TSS=False` | 确实无需 TSS 时设为 `True`；不能同时指定外部 `TSS_BED` |

### CLI

仅计算 QC 可直接运行 Processor，之后再打开 Notebook 绘图。CLI 不生成 `SingleCpG_Adata.h5ad`。

```bash
PIPELINE=/path/to/Alopex
PATH="$PIPELINE/conda/current/control/bin:$PATH" \
  "$PIPELINE/conda/current/notebook/bin/python" \
  "$PIPELINE/downstream/downstream_qc.py" \
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

结果写入 `<project>/06_downstream/QC_Results_YYYYMMDD_HHMMSS/`，时间取首次创建的本地时间，同秒冲突加微秒。`.delivery.json` 绑定交付：同交付复用目录与缓存，新交付新建目录；多个匹配会拒绝，不按 mtime 或目录名猜测。

| 文件 / 目录 | 用途 |
| --- | --- |
| `DNAme_QC_Information.csv` | 下游与 adata 合并的唯一 QC 表，固定20列；[字段、顺序与单位][qc-fields] |
| `RawAdata.h5ad` / `RawAdata.input_schema.json` | 全细胞 canonical CpG 导入矩阵与复用记录 |
| `SingleCpG_Adata.h5ad` / `SingleCpG_Adata.input_schema.json` | Notebook 最后一格生成的全细胞 single-CpG 矩阵与完成记录 |
| `TSS_Profile_Information.csv` | 未跳过 TSS 时生成，仅含抽中且有覆盖的 bins |
| `plots/` | QC、Final QC、孔板、信号组成、TSS、CpG density 与 PCA 等 PDF |
| `.cache/` | CpG 统计、信号组成、Gini 与 QC 完成身份等内部缓存，不手动改写 |

同一交付只能有一个下游写入者。遇到 `Downstream QC is busy`，等当前 Processor 或 SingleCpG 更新结束后重试。项目和 QC 目录内的 `.qc-write.lock` 是持久锁文件，**不要删除或替换**；进程退出后系统会释放锁，文件仍存在不代表任务在运行。

### QC 字段

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
| `Duplicate_Pair_Rate%`、`Trim_Retention%`、`Final_Pair_Yield%` | 去重率为duplicate/accepted pairs，trim retention为trimmed/input pairs，final yield为final/input pairs（均×100）；Droplet的pair去重率为NaN |

## 指标与筛选

读图时先区分位点均值、coverage 比例和 read-pair 指标。身份列按字符串保留前导零及合法 `NA/nan/None`，数值缺失保持 NaN；观测零值不能当作缺失。字段见[QC 字段][qc-fields]，分母要点见下表。

| 指标 | 解读要点 |
| --- | --- |
| QC 甲基化率 | 宿主、lambda、pUC19 分别对已覆盖位点的 `m/(m+u)` 取均值再×100；无观测为 NaN |
| SingleCpG | 区间内已覆盖位点的 fraction 均值；`1/4` 与 `9/12` 的均值为0.5。有覆盖全 unmethyl 保显式0，缺覆盖不补0 |
| TSS | 每个 offset 对有覆盖 bins 的 fraction 取均值；count 为甲基化信号之和 |
| Gini | 按有 coverage 的固定宽 bins 计算，不默认排除 chrM/chrY；仅空 CpG 细胞缺失 |
| 信号组成 | 各类 canonical CpG coverage / 全部 contig 总 coverage；其他 decoy 仍进入分母。controls 仅精确匹配 `lambda` / `pUC19` |
| mapping / pair 指标 | mapping 的 backend 单位与分母不同；完整 pair funnel、policy/unit 保留在 sealed final manifest 与 MultiQC。Droplet 的 pair 去重率为 NaN，不能用位点 UMI 数代替 pair 数 |
| TAPS | controls 使用合并计数比率；Non-CpG 为可判链宿主 CHG/CHH 转换率，不能直接解释为真实 mCH 或 cDNA 污染率，High-CpH 为 NaN |

### HQ 标记

Notebook 默认同时要求 lambda methylation `<10%`、pUC19 methylation `>90%`、宿主 unique CpG sites `>500000`、Gini `<0.5`。mapping 默认不参与；对应数据完成验证后才设置 `MAPPING_THRESHOLD`。pair 去重率默认不参与。

`HQ` 是细胞注释；Raw 与 SingleCpG 仍保留全部细胞，SingleCpG 同时保存 QC、`HQ` 和自定义注释。

### 图与展示控制

| 图 / Notebook 控制 | 解读与输出 |
| --- | --- |
| Final QC | `FINAL_QC_MODE=per_clone_samples` 按 CloneID 分面，`distribution` 叠加；`FINAL_QC_METRIC` 支持单指标/列表，`FINAL_QC_CLONE_ID` 选组，逐项展示并保存 PDF |
| Final QC 横纵轴 | 100 bins、组间共用边界；CpG 正值用 log 横轴，含非正值改 signed-log1p 并标明，其他指标线性；纵轴 Counts，保留有效 N/中位数，实际直方图加 sigma=1.5、半径6、对称边界的高斯平滑，在变换空间插值；平滑保总计数且不参与 QC 筛选，仅排缺失/非数值/无穷 |
| `CpG_Density_Unique_vs_Input_Pairs.pdf` | 横轴 trimming 前 input pairs，纵轴最终宿主 unique CpG sites，CloneID 配色；散点用原始计数，不作深度归一化或异常点剔除 |
| 拟合与 `CpG_Density.pdf` | 散点拟合线与斜率汇总共用联合有效样本、异常点筛选及 `min_slope_ratio_change`；过原点 `slope=sum(x*y)/sum(x*x)`，单位 sites/input pair，不能拟合则 N/A |
| density 输入 | QC 身份/sites、sealed final manifest pairs、composition cache host signals 的样本集合须完全一致，不静默取交集 |
| 其他图 | `WGBS_QC_Panel.pdf`、`Plate_View_Native_Mapping%.pdf`、`CpG_Signal_Composition.pdf`、可选 `TSS_Profile.pdf`、`Batch_Effect_PCA.pdf` |

## 保存供下游复用的 single-CpG 矩阵

运行 Notebook 最后一格的 `prepare_single_cpg_adata(PATHS, df_qc)`，按 `Sample_ID` 严格对齐全细胞 QC 与注释。完成后使用打印出的 `SINGLECPG_ADATA_PATH` 读取：

```python
import snapatac2 as snap

adata = snap.read(SINGLECPG_ADATA_PATH, backed="r")
# 在此使用只读对象；结束后关闭句柄。
adata.close()
```

| 使用边界 | 要求 |
| --- | --- |
| 完成记录 | 矩阵与 `SingleCpG_Adata.input_schema.json` 配套使用；来源未知或失配时按提示显式重建，见[缓存与强制重算](#缓存与强制重算) |
| 坐标 | 跟随 sealed `analysis.cpg_coordinate_base`：0-based `chr1:0-2` 为半开区间，1-based `chr1:1-2` 为闭区间，表示同一 CpG；传给 BED/0-based 工具前须按实际 metadata 转换 |
| reference 顺序 | 保留原生 `reference_sequences` 与索引顺序，不按 FAI 重排；真实长度为 `reference_bp_length=L`，索引长度为 `reference_seq_length=L+base`，base1 时索引0留空 |
| 后续注释 | 只读对象不能回填标签；使用独立可写副本或合适大小的内存子集。MethylTree 分箱后仍按细胞 ID 对齐 `obs` |

固定计算参数见[single-CpG 参数][single-params]；矩阵的 `uns["cpg_coordinate_base"]`、`uns["cpg_interval_convention"]`、`uns["source_cpg_representation"]` 记录实际坐标和来源。

## TSS profile

默认 reference 为 `resources/<genome>_reference/tss/<genome>_TSS_2000_2000_20.bed`，使用链方向一致的 ±2000 bp / 20 bp bins，第9列为 offset。外部 BED 须匹配 sealed genome build，并通过完整坐标、重复区间和 offset 校验。reference 由 Doctor 根据已准备的 FASTA/GTF 构建，准备步骤见[主 README](../README.md)。

Notebook 第二格的 `TSS_SAMPLE_FRACTION`、`TSS_SAMPLE_SEED` 控制抽检，对应 CLI `--tss-sample-fraction` / `--tss-sample-seed`。默认每个原始 CloneID 抽10%：不足10个全检，其余向上取整且最多100个；比例设为1.0时全检。选样按种子与 Sample_ID 确定，不受行序影响，也不按 coverage 选样。完整范围与记录位置见[TSS 抽检控制][recompute]。

抽检只影响 TSS，其他 QC 和矩阵仍保留全部细胞。CSV 仅记录选中且有覆盖的 bins，实际出现的细胞数可小于抽检数；有 coverage 的零甲基化仍保留，无覆盖不补0。全体无覆盖时仅表头 CSV 合法，旧 PDF 会清理。显式 `--skip-tss` / `SKIP_TSS=True` 会清理该交付旧 TSS CSV/PDF。

## 缓存与强制重算

默认按 sealed 输入、参数与算法身份复用缓存，无需先删除输出目录。新 Raw generation 会使 Gini/TSS 重算；canonical CpG 统计可独立复用。重算开关默认均为 `False`，只影响下游。

| Notebook 开关 | CLI 对应项 | 重建对象 |
| --- | --- | --- |
| `RECOMPUTE_RAW_ADATA` | `--recompute-raw-adata` | 从 sealed CpG 重建 Raw |
| `RECOMPUTE_CPG_STATS` | `--recompute-cpg-stats` | 联合重建 CpG 位点/甲基化统计与信号组成；任一缓存缺失或损坏也会自动联合重建 |
| `RECOMPUTE_GINI` | `--recompute-gini` | Gini |
| `RECOMPUTE_TSS` | `--recompute-tss` | TSS profile |
| `RECOMPUTE_SINGLECPG` | 无，仅 Notebook 最后一格 | 重建 SingleCpG；完成后恢复 `False` |

| 遇到的情况 | 处理方式 |
| --- | --- |
| 调整 TSS BED、比例或种子 | 仅 TSS 自动失效，无需强制重建其他 QC 或矩阵 |
| QC 完成身份缺失、过期或重算中断 | 重新运行 Processor，成功后再加载 QC 或生成 SingleCpG；有效 Raw 仍可复用 |
| Raw 输入或算法身份过期 | 先运行 Processor 完成必要重建，再处理 SingleCpG；不能仅修改 metadata 或用 Single 强制开关绕过 |
| SingleCpG 来源未知或变化 | 按提示设置 `RECOMPUTE_SINGLECPG=True`；仅注释变化时可复用矩阵并更新 `obs`，原 QC 列须与当前磁盘 QC 一致 |
| 输入或参考在计算期间变化 | 确认文件稳定后重跑；失败时不提交完成身份，未完成矩阵不能作为有效缓存使用 |
| CpG 格式错误或缓存计数矛盾 | 按错误定位 canonical CpG 输入并重算；不把格式失败当作合法空细胞，也不手改缓存绕过校验；输入格式与守恒见[下游契约][qc-contract] |

模块更新后，结束旧写任务并重启 kernel。最终 QC 提交中断后须重跑 Processor 恢复，再继续消费结果。缓存与提交的恢复方式见[缓存与强制重算][qc-contract]。

## 关联 RNA_DARLIN 结果

Cabernet、TAPS、Droplet 均在 DNA manifest 的可选 `rna_sample` 填 RNA_DARLIN **原始 RNA basename**，不是 RNA `sample_id` 或路径；两边项目样本可独立命名，Alopex 只处理本地 DNA。SRD 的 `rna_sample` 则是本地 RNA-enrichment tube，不作外部桥接。

先完成 DNA 发布和同代 downstream QC，再在 RNA Notebook 填 `DNA_PROJECT_ROOT`。RNA 下游只读当前 complete DNA 的 sealed manifest 与同代 `DNAme_QC_Information.csv`，不读取可变 config 或按 mtime 挑目录。修改关联后须重发布 DNA、生成同代 downstream QC，再重新整合 RNA。

| 协议 | 文库与细胞关联 |
| --- | --- |
| Cabernet / TAPS | sealed `rna_sample` → finalized RNA.h5ad 原始 `rna_sample` → RNA canonical `sample_id`，再按 PlateID 一对一连接；未填外部关联时使用共享项目 `sample_id` 与 PlateID |
| SRD | 使用共享项目 `sample_id` 与 PlateID |
| Droplet | 同样从 sealed `rna_sample` 解析文库，再按对应文库内完整 CB17 精确连接；留空保持未关联，不自动猜配 |

孔板协议要求参与 QC 的两侧细胞集合完全一致，barcode/Cell_Order 冲突、多对一或非身份列重名均失败，不静默取交集。Droplet 两侧独立 calling/QC，`DNA_RNA_cell_links.csv` 保留双方全集，整合 QC 只含共同且双方 QC 可用的细胞；未配对细胞不补科学零值。文库冲突或错误交付代际均拒绝。

`Integrated_QC_Information.csv` 保留 RNA 身份，RNA 13列与 DNA 20列的独立 QC schema 不变；RNA 下游记录实际关联和输入哈希。文库与细胞的关联边界按本节分别检查。

[qc-fields]: #qc-字段
[denominators]: #指标与筛选
[recompute]: #缓存与强制重算
[single-params]: #保存供下游复用的-single-cpg-矩阵
[single-contract]: #保存供下游复用的-single-cpg-矩阵
[qc-contract]: #缓存与强制重算
[manifest-contract]: #关联-rna_darlin-结果
[droplet-contract]: #关联-rna_darlin-结果
