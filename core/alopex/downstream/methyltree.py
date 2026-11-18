"""将当前 Alopex 交付接入 MethylTree；不转换公开 CpG，也不修改 MethylTree。"""
from __future__ import annotations

import importlib.metadata
import json
import os
from pathlib import Path
import sys
import uuid

_LOADED_MODULE_CODE = sys._getframe().f_code

from . import qc
from ..identity import capture_loaded_source_codes, assert_loaded_sources_current, PythonSourceGraph
from ..paths import CORE_ROOT

_LOADED_SOURCE_CODES = capture_loaded_source_codes()


def _sample_info(frame, cells, nCG_cutoff):
    """保留原 QC 和 sealed 身份；MethylTree 的别名只改变字段名或百分比单位。"""
    ids = [cell["sample_id"] for cell in cells]
    frame = qc._require_exact_sample_ids(frame, ids, "MethylTree QC").copy()
    frame["sample"] = frame["Sample_ID"]
    frame["cell"] = frame["Sample_ID"]
    frame["clone_id"] = frame["CloneID"]
    frame["plate_id"] = frame["PlateID"]
    frame["plate_barcode"] = [cell["dna_barcode"] for cell in cells]
    frame["cell_order"] = [cell["cell_order"] for cell in cells]
    for column, key in (("CloneID", "project_sample_id"), ("PlateID", "plate_id")):
        if frame[column].tolist() != [cell[key] for cell in cells]:
            raise ValueError(f"MethylTree {column} differs from the sealed manifest")
    frame["nCG"] = frame["Sample_Unique_CpG_Sites"].astype("int64")
    for target, source in (
        ("met_rate", "Sample_Meth_CpG_Rate%"),
        ("lambda_rate", "Lambda_Meth_CpG_Rate%"),
        ("pUC19_rate", "pUC19_Meth_CpG_Rate%"),
        ("lambda_reads_ratio", "Signal_Lambda_Rate%"),
        ("puc19_reads_ratio", "Signal_pUC19_Rate%"),
        ("chrM_reads_ratio", "Signal_mtDNA_Rate%"),
    ):
        frame[target] = frame[source] / 100.0
    frame["lambda_pUC19_reads_ratio"] = frame["lambda_reads_ratio"] + frame["puc19_reads_ratio"]
    frame["lambda_nCG"] = frame["Lambda_Unique_CpG_Sites"].astype("int64")
    frame["pUC19_nCG"] = frame["pUC19_Unique_CpG_Sites"].astype("int64")
    frame["Mapping_efficiency"] = frame["Native_Mapping%"]
    # 此处 reads 是拆分后本细胞、trimming 前的 R1 数，即 input pairs；不是文库原始总量。
    frame["reads"] = [int(cell["cutadapt_input_pairs"]) for cell in cells]
    for column in ("native_mapping_policy", "native_mapping_unit", "methylation_backend"):
        frame[column] = [cell[column] for cell in cells]
    frame["source"] = "met"
    frame["pass_metQC"] = frame["nCG"] > nCG_cutoff
    frame["HQ"] = frame["pass_metQC"]
    return frame


def _write_csr_subset(source, destination, obs, provenance, before_commit):
    """有界复制CSR及HDF5注释；不将全参考var标签或整个X载入内存。"""
    import h5py
    from anndata._io.specs import write_elem
    from . import matrix

    reader = matrix.CSRReader(source)
    writer = None
    try:
        rows = reader.obs_names.get_indexer(obs.index)
        if (rows < 0).any():
            raise ValueError("MethylTree selection contains cells absent from SingleCpG")
        writer = matrix.CSRWriter(destination, obs, n_vars=reader.n_vars, uns=reader.uns)
        for row in rows:
            writer.append(reader.read_rows(int(row), int(row) + 1))
        writer.close()
        writer = None
        before_commit()
        with h5py.File(source, "r") as src, h5py.File(destination, "r+") as dst:
            del dst["var"]
            src.copy("var", dst)
            del dst["uns"]
            src.copy("uns", dst)
            write_elem(dst["uns"], "alopex_methyltree_input", provenance)
        before_commit()
    finally:
        reader.close()
        if writer is not None:
            writer.close()


def prepare_methyltree_inputs(
    paths: qc.NotebookPaths,
    *,
    nCG_cutoff: int = 100000,
    import_nCG_cutoff: int = 500000,
):
    """从已完成 QC 导出全部细胞注释及所选细胞的独立、标准 backed H5AD。

    先运行 qc.run_qc_processor。公开输出只读；复用标准 SingleCpG 的科学校验，
    导出始终覆盖 06_downstream/methyltree/allCpGs.h5ad，不复用用户分析后的文件。
    双链按已覆盖胞嘧啶 fraction 的均值合并，不改成合并计数之比；原始整数保存在已发布的四列 CpG。
    """
    import pandas as pd

    for name, value in (("nCG_cutoff", nCG_cutoff), ("import_nCG_cutoff", import_nCG_cutoff)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f"{name} must be a nonnegative integer")
    assert_loaded_sources_current(_LOADED_SOURCE_CODES)
    with qc._qc_write_lock(paths.qc_results_dir):
        sealed = qc.load_sealed_qc_context(paths.project_dir, expected_pipeline_root=paths.pipeline_root)
        if sealed.delivery_id != paths.delivery_id:
            raise ValueError("Delivery changed; rerun the Notebook from its first cell")
        run_paths = qc.build_paths(paths.project_dir)
        if run_paths.qc_results_dir != paths.qc_results_dir:
            raise ValueError("Notebook QC directory differs from the current delivery")
        completed = qc._require_completed_qc(run_paths, sealed)
        cells = qc.read_active_cells(paths.project_dir, sealed_context=sealed)
        frame = pd.read_csv(paths.qc_info, converters=dict.fromkeys(qc.QC_INFO_COLUMNS[:3], str),
                            float_precision="round_trip")
        frame = qc._require_exact_sample_ids(frame, list(completed.sample_ids), "MethylTree QC")
        samples = _sample_info(frame, cells, nCG_cutoff)
        selected = samples.loc[samples["nCG"] >= import_nCG_cutoff].copy()
        if selected.empty:
            raise ValueError("No cells meet import_nCG_cutoff; no MethylTree outputs were replaced")
        single = run_paths.single_cpg
        schema = qc._require_current_single(run_paths, completed.identity)

        def assert_current():
            assert_loaded_sources_current(_LOADED_SOURCE_CODES)
            completed.assert_current()

        code = PythonSourceGraph((CORE_ROOT,), loaded_sources=_LOADED_SOURCE_CODES).dependencies(
            ("alopex.downstream.methyltree._sample_info", "alopex.downstream.methyltree._write_csr_subset"))
        provenance = {
            "delivery_id": sealed.delivery_id,
            "single_input": json.dumps(schema["input"], sort_keys=True),
            "source_code": json.dumps(code, sort_keys=True),
            "anndata": importlib.metadata.version("anndata"),
            "h5py": importlib.metadata.version("h5py"),
            "nCG_cutoff": nCG_cutoff, "import_nCG_cutoff": import_nCG_cutoff,
            "missing_value": "absent CSR entry; observed zero is explicitly stored",
            "strand_merge": "mean of covered cytosine fractions, not pooled counts",
        }
        output_dir = qc._project_directory(paths.project_dir, Path("06_downstream/methyltree"), create=True)
        output = output_dir / "allCpGs.h5ad"
        sample_sheet = output_dir / "sample_sheet.csv"
        completion = output_dir / "allCpGs.input_schema.json"
        temporary = qc.project_temporary_parent(output) / f".methyltree.{uuid.uuid4().hex}.h5ad"
        try:
            selected.index = pd.Index(selected["sample"], name=None)
            _write_csr_subset(single, temporary, selected, provenance, assert_current)
            qc._fsync_file(temporary)
            assert_current()
            qc._invalidate_completion(completion)
            os.replace(temporary, output)
            qc._fsync_directory(output_dir)
            qc._atomic_write_csv(samples, sample_sheet, before_commit=assert_current)
            qc._atomic_write_json({"input": provenance, "h5ad": qc._published_file_binding(output),
                                   "sample_sheet": qc._published_file_binding(sample_sheet)},
                                  completion, before_commit=assert_current)
            assert_current()
            return samples, output
        finally:
            temporary.unlink(missing_ok=True)
