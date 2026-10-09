# Alopex v14.6 下游 QC 与可视化

这套下游分析以当前 `Alopex` 的 sealed delivery 为输入，通过
`03_results/run_manifest.json` 解析 immutable config snapshot、reference FAI、
`03_results/QC_Results/sample_manifest.tsv` 与其声明的 `03_results/CpG/` 文件。
主 Pipeline 已把下游需要的 normalized QC 指标写入 final manifest，因此 Downstream
只依赖 sealed `03_results/` 与 reference；`02_work/` 是否保留不影响下游结果。

Processor 接受 `analysis.methylation_backend: biscuit`、`bismark` 与 `rastair`。三条链都从
final manifest 的 canonical read-pair funnel 进入相同计算；backend 原生 MAPQ、比对
报告和 M-bias 等诊断仍保留在主 Pipeline 的 backend-specific MultiQC section，分别标注其原生口径。

## 运行

Notebook/Downstream 使用 Doctor 管理的 `conda/current/notebook` 环境，规格唯一定义在 `core/envs.yaml`。可直接执行：

```bash
PIPELINE=/path/to/Alopex
"$PIPELINE/conda/current/notebook/bin/python" \
  "$PIPELINE/downstream/downstream_qc.py" \
  --project-dir /path/to/Patient001
```

`--workers=0` 自动选择 CpG worker 数（最多 32）；显式值同样受 affinity、Slurm allocation
和调用方预算限制，每个 worker 内 native 线程固定为 1。Notebook 启动器使用同一 CPU 预算规则。

完成后打开 `downstream/DNA_QC_Visualization.ipynb`。第一格包含可编辑参数。
Pipeline 与项目通常不在同一目录时，直接填写两个绝对路径；显式值优先于环境变量，
且不会静默回退到其他目录：

```python
DNA_PROJECT_ROOT = Path("/path/to/Patient001")
DNA_PIPELINE_ROOT = Path("/path/to/Alopex")
```

也可以改为导出同名环境变量。两者都留空时，Notebook 先检查 Jupyter 当前目录及父目录中的
completed delivery，再通过 sealed run manifest 解析 Pipeline。

### 缓存与强制重算

Notebook 第二格的 `RECOMPUTE_RAW_ADATA`、`RECOMPUTE_METHYL_STATS`、
`RECOMPUTE_COMPOSITION`、`RECOMPUTE_GINI`、`RECOMPUTE_TSS` 默认均为 `False`，
按输入和算法身份自动复用缓存，CSV 中的科学浮点值读回时保持原精度。RawAdata 与统计缓存的
CpG 输入身份复用 sealed inventory 中的 size/hash，不重复读取所有 CpG 计算哈希。需要强制重建 AnnData 时将 `RECOMPUTE_RAW_ADATA`
设为 `True`；它从 sealed CpG 重建 `RawAdata.h5ad`，结果仍只写入 `06_downstream/`。

RawAdata 的 schema、文件身份、样本集合及当前 sealed 输入全部匹配时，直接复用其空 CpG 集合，不逐细胞启动解压或创建临时链接。失配时重新检查表头和空表，确认 zero-mapped cells 为空后再导入；临时链接只在重建时创建。导入算法更新仍触发正常重建，不复用旧算法身份。
命令行对应 `--recompute-raw-adata`；其余开关对应 `--recompute-methyl-stats`、
`--recompute-composition`、`--recompute-gini`、`--recompute-tss`。

显式指定 Pipeline 路径也要求项目的交付状态为 `complete`、Pipeline 身份有效。
Processor 子进程沿用 Notebook 的 Python，运行前请使用上述专用 kernel。

### 保存供下游复用的 single-CpG 矩阵

运行 Notebook 最后一格会调用 `prepare_single_cpg_adata(PATHS, df_qc)`，自动使用当前物种的
single-CpG BED 和本次交付的 `RawAdata.h5ad`，生成同目录的 `SingleCpG_Adata.h5ad`。
固定参数为 `value_type="fraction"`、`summary_type="mean"`、`chunk_size=500`、`inplace=False`。
保留全部细胞，按 `Sample_ID` 严格对齐并保存完整 QC、`HQ` 及自行添加的注释，不按 HQ 筛选。
不要修改 `df_qc` 的原有 QC 列；它们会与磁盘上的当前 QC 表校验。

首次运行需要完成矩阵计算。再次运行按交付、RawAdata generation、BED 内容、参数与
SnapATAC2 版本复用矩阵，只有注释变化时仅更新 `obs`。已有文件未经验证或来源变化时会提示；
需要替换它时设置 `RECOMPUTE_SINGLECPG=True`，完成后恢复 `False`。新矩阵先写临时文件，
关闭后替换正式文件；`SingleCpG_Adata.input_schema.json` 是配套的完成及复用记录。
函数返回路径并关闭自身打开的全部句柄，不保留 `DNA_ADATA_CpG` 活跃写句柄。

下游把打印出的 `SINGLECPG_ADATA_PATH` 用于 `snap.read(..., backed="r")` 即可读取。
只读 SnapATAC2 对象不能回填标签；修改标签应使用独立可写副本或合适大小的内存子集。
MethylTree 分箱后仍须检查并按细胞 ID 对齐 obs，该预生成步骤不修改外部 MethylTree 的转换行为。

### Jupyter / VS Code kernel

QC Notebook 元数据要求名为 `alopex-qc` 的 kernel。交互式运行前把它注册到当前用户的
Jupyter kernel 目录；argv 经由 `conda/current` 符号链接指向 Doctor 管理的 notebook 环境，
后续环境重建不需要重新注册。Pipeline 迁移目录后应重新注册并在 Notebook 中选择 `Alopex QC`；确认新 kernel 可启动后移除指向旧目录的专用 kernel。以下命令自动选择当前系统的 Jupyter 数据目录：

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

命令行执行已有 delivery 时直接用 nbconvert：`jupyter nbconvert --to notebook --execute --inplace <notebook>`，
并把工作目录放在项目 `03_results/` 下以便第一格自动发现 project root。

### TSS profile

Processor 会按项目 `species` 在 `resources/hg38_reference/tss/` 等现有 reference
目录中寻找以下文件，并生成 TSS CSV/PDF：

```text
hg38_TSS_2000_2000_20.bed
```

`core/doctor.sh reference` 会从 `resources/reference_source/hg38.genes.gtf`（GENCODE
v48）或 `mm10.genes.gtf`（GENCODE M25）的 protein-coding `gene` 记录自动构建唯一的
TSS BED（链方向一致的 `±2000 bp / 20 bp` bins，第 9 列为 offset）。脚本还会生成
`resources/<genome>_reference/cpg/<genome>.single_cpg.bed.gz`，供 SnapATAC2
按需 `make_peak_matrix(..., peak_file=...)` 使用的 BED3：0-based、half-open，每条记录覆盖一个
CpG 二核苷酸 `[C, G+1)`，不含 lambda/pUC19。Notebook 最后一格从已导入 values 生成 single-CpG matrix；命令行 QC Processor 不执行这一步。Pipeline 的 `03_results/CpG/*.cpg.tsv.zst`
则直接采用 SnapATAC2 `pp.import_values` 四列输入：`chrom, pos, methyl, unmethyl`，其中
`pos` 为 0-based；Downstream 直接读取 CpG 文件。

`resources/*_reference/` 是体积很大的 HPC 本地 reference，受仓库 `.gitignore`
保护，不进入源代码包，也不应使用 `git add -f` 提交。普通 clone 在按根
README 准备 4 个 FASTA 与 2 个 GTF 后运行 Doctor，即会构建这些资源；本 Pipeline
不会自动下载 reference。

要覆盖其它 genome build 或外部资源，可显式提供与当前 reference 匹配的 TSS BED：

```bash
"$PIPELINE/conda/current/notebook/bin/python" "$PIPELINE/downstream/downstream_qc.py" \
  --project-dir /path/to/Patient001 \
  --tss-bed /path/to/genome_TSS_2000_2000_20.bed
```

Notebook 首格也提供 `TSS_BED` 和 `SKIP_TSS` 设置。默认
`TSS_BED=None`、`SKIP_TSS=False`，即自动使用当前 species 的 Pipeline
资源。TSS BED 必须与 sealed delivery 使用同一 reference build。其它 species 确实
没有 TSS 资源时，才显式设置 `SKIP_TSS=True` 或使用 CLI `--skip-tss`。

TSS profile 只记录有覆盖的 bins。有 CpG 的 cell 也可能没有 TSS 覆盖；全体无覆盖时生成并复用仅表头 CSV，其他 QC 正常完成，并清理旧 TSS PDF。显式跳过 TSS 时同时清理该交付的旧 TSS CSV/PDF，避免旧曲线被误认为本次结果。真实覆盖且甲基化为零的 bin 保留数值 `0`。
外部 TSS BED 的 genomic interval 必须唯一，第 9 列 offset 必须为有限的 32-bit 整数；重复 interval 或错误 offset 会被明确拒绝。

## 输出

输出写入与当前交付绑定的
`<project>/06_downstream/<delivery_id>/QC_Results/`，不会修改 sealed
`03_results/`：

- `DNAme_QC_Information.csv`：下游与 adata 合并的唯一 QC 表，固定 20 列且
  列序唯一出自 `downstream_qc.QC_INFO_COLUMNS`：`Sample_ID`、`CloneID`、`PlateID`、
  `Native_Mapping%`、`Sample_Unique_CpG_Sites`、`Sample_Meth_CpG_Rate%`、
  `Lambda_Unique_CpG_Sites`、`Lambda_Meth_CpG_Rate%`、`pUC19_Unique_CpG_Sites`、
  `pUC19_Meth_CpG_Rate%`、`Non_CpG_Methylation%`、`High_CpH_Flag_Rate%`、
  `Signal_Host_Rate%`、`Signal_Lambda_Rate%`、`Signal_pUC19_Rate%`、
  `Signal_mtDNA_Rate%`、`Gini_Index`、`Duplicate_Pair_Rate%`、`Trim_Retention%`、
  `Final_Pair_Yield%`。完整 read-pair funnel 与 backend policy/unit 诊断不复制进
  本表，保留在主 Pipeline 的 MultiQC 报告。Droplet 的 `Duplicate_Pair_Rate%` 为 NaN，
  `Non_CpG_Methylation%` 来自过滤后 CHG/CHH 的胞嘧啶/UMI 共识；分子观测数保留在
  final manifest 与 MultiQC，不与 read-pair 数混算。默认 HQ 判定不使用 pair 去重率。
- `TSS_Profile_Information.csv`（仅未跳过 TSS 时）
- `RawAdata.h5ad`
- `SingleCpG_Adata.h5ad` 与 `SingleCpG_Adata.input_schema.json`（运行 Notebook 最后一格后）
- `.cache/Mapping_Stats.csv`（内部缓存）
- `.cache/CpG_Signal_Composition.csv`（内部缓存）
- `.cache/Gini_Index.csv`（内部缓存）
- `plots/` 下的 PDF，包括 `WGBS_QC_Panel.pdf`、Final QC dashboards、
  `Plate_View_Native_Mapping%.pdf`、`CpG_Signal_Composition.pdf`、
  `TSS_Profile.pdf`（可选）、`CpG_Density.pdf` 与 `Batch_Effect_PCA.pdf`。

`CpG_Density.pdf` 的输入由 Notebook 经 `downstream_qc.load_cpg_density_frame`
合并三个来源：`Sample_ID`、`CloneID` 和 `Sample_Unique_CpG_Sites` 取自
20 列 QC 表，pair 计数取自 sealed `03_results/QC_Results/sample_manifest.tsv`
（`build_qc_records` 唯一加载链），`Sample_CpG_Signals` 取自
`.cache/CpG_Signal_Composition.csv`；三个来源的样本集合必须完全一致。
QC 表中的 `Sample_ID`、`CloneID`、`PlateID` 以及 TSS 的 `CloneID` 按字符串读取，保留 `001` 等前导零以及
`NA`、`nan`、`None` 等合法名称；数值列仍按数值和缺失值解析，并保留 round-trip 浮点精度。

RawAdata 的样本集合以 sealed final manifest 为准；临时 CpG 导入文件只调整格式后缀，样本名内部的 `.cpg.` 原样保留。

Notebook 的 `FINAL_QC_MODE` 可选 `per_clone_samples` 或 `distribution`；前者的 `FINAL_QC_METRIC` 可填单个指标或指标列表。


## 指标与筛选

20 列 QC 的分母、backend 原生 mapping 单位与 TAPS 控制含义的统一定义保留在开发仓库
参数表。不同 backend 的
`Native_Mapping%` 不视为同义指标；TAPS 的 `High_CpH_Flag_Rate%` 为缺失值；`Non_CpG_Methylation%` 默认显示可判链宿主分子的 CHG/CHH 原始转换率，零观测为缺失。该指标未经 SNP/测序错误校正，只代表可判链子集，不能直接解释为真实 mCH 或 cDNA 污染率。分项计数与判链分母见 final manifest/MultiQC。
Gini 使用有覆盖 bins 的 coverage，不纳入未覆盖 bins；TSS count 为甲基化信号之和。

Notebook 默认 HQ 使用 Lambda methylation `< 10`、pUC19 methylation `> 90`、unique CpG
sites `> 500000` 且 Gini `< 0.5`；mapping 默认不参与统一 HQ，只有完成对应数据回归后才设置
`MAPPING_THRESHOLD`。真实数据与合成对照的验证边界见本版 Release 说明。

## 关联 RNA_DARLIN 结果

Cabernet/TAPS 的 DNA 项目 manifest 可填写外部 `rna_sample`，值与 RNA_DARLIN 项目的
原始 RNA basename 一致；两边 `sample_id` 可独立命名。Alopex 只处理本地
`dna_raw_sample`，将关联保存在最终 `03_results/QC_Results/sample_manifest.tsv`。
先完成 DNA 发布及本 downstream，再在 RNA notebook 填写 `DNA_PROJECT_ROOT`。
RNA Processor 从当前 complete DNA delivery 读取 sealed manifest 与同代际本表，
用 RNA.h5ad 的 `rna_sample → sample_id` 和 PlateID 对齐，发布 RNA 侧
`Integrated_QC_Information.csv`。SRD 的本地 `rna_sample` 不参与外部桥接。
文库映射须一对一、参与细胞集合须相等，身份/条码/Cell_Order 冲突直接报错；
不会读取未发布的 config 修改或选择 mtime 较新的旧 QC。DNA 20 列 QC 不增加关联列，
实际文库映射及输入哈希保存在 RNA `analysis_summary.json`。更新关联后先重跑 DNA 发布
和当前代际 downstream，再重新整合。

Droplet 的配套 RNA 独立进行细胞调用与 QC。RNA 项目用
`00_config/dna_rna_library_map.tsv` 指定 DNA/RNA 文库，再按完整 CB17 关联；DNA manifest 的
`rna_sample` 留空。RNA 侧输出双方细胞全集 `DNA_RNA_cell_links.csv`，整合 QC 只含共同且双方
QC 可用的细胞；文库冲突或错误交付代际直接失败。该流程不同于孔板的完整细胞集合对应，
完整维护契约保留在开发仓库。
