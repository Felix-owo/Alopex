"""逐染色体、TSS profile及逐细胞筛选解释；不以缺失补零或改变细胞集合。"""
from __future__ import annotations

import sys
_LOADED_MODULE_CODE = sys._getframe().f_code

import numpy as np
import pandas as pd
from scipy import sparse
from . import matrix

CHROMOSOME_INFO = "Chromosome_QC_Information.csv"
FILTER_INFO = "QC_Filter_Reasons.csv"


def tss_intervals(path, reference, host_contigs):
    """完整校验原始 gene/bin 成员，再按实际主核参考构建数字索引。"""
    columns = {0: "chrom", 1: "start", 2: "end", 6: "gene", 8: "offset"}
    try:
        frame = pd.concat(matrix.bed_blocks(path, columns, {"start": np.int64, "end": np.int64,
                                                            "offset": np.int32}), ignore_index=True)
    except pd.errors.ParserError as exc:
        raise ValueError(f"TSS BED must have at least 9 columns: {path}") from exc
    if frame.empty:
        raise ValueError(f"TSS BED contains no intervals: {path}")
    if (frame.start < 0).any() or (frame.end <= frame.start).any():
        raise ValueError(f"TSS BED must contain nonempty nonnegative intervals: {path}")
    if frame.gene.str.strip().eq("").any():
        raise ValueError(f"TSS BED contains an empty gene_id: {path}")
    lengths = pd.Series(reference.reference_bp_length.to_numpy(), index=reference.reference_seq_name,
                        dtype="Int64")
    limits = frame.chrom.map(lengths)
    clipped = (frame.start < limits) & (frame.end > limits)
    frame.loc[clipped, "end"] = limits[clipped].to_numpy(dtype=np.int64)
    if frame.duplicated(["chrom", "start", "end", "gene"]).any():
        raise ValueError(f"TSS BED contains duplicate gene/bin membership: {path}")
    if frame.duplicated(["gene", "offset"]).any():
        raise ValueError(f"TSS BED assigns multiple intervals to the same gene/offset: {path}")
    keys = ["chrom", "start", "end"]
    frame = frame.loc[frame.chrom.isin(host_contigs) & (frame.start < limits), [*keys, "offset"]]
    # 首次出现顺序即区间行号；同一区间的每个 gene/bin 成员仍各贡献一次。
    bins = frame.groupby(keys, sort=False).ngroup().to_numpy(dtype=np.int64)
    regions = frame.loc[~pd.Index(bins).duplicated(), keys].reset_index(drop=True)
    codes, offsets = pd.factorize(frame.offset, sort=True)
    membership = sparse.csr_matrix((np.ones(len(frame), dtype=np.int64), (bins, codes)),
                                   shape=(len(regions), len(offsets)))
    index = {
        chrom: matrix.interval_lookup(rows.start.to_numpy(), rows.end.to_numpy(), rows.index.to_numpy())
        for chrom, rows in regions.groupby("chrom", sort=False)
    }
    return index, membership, np.asarray(offsets, dtype=np.int64)


def tss_profile_row(reader, row, sid, clone, tss_index):
    index, membership, offsets = tss_index
    cols, frac, methyl = matrix.aggregate_row(reader, row, index)
    if not len(cols):
        return None
    weights = membership[cols]
    count = np.asarray(weights.sum(axis=0)).ravel()
    present = count > 0
    means = np.asarray(frac @ weights).ravel()[present] / count[present]
    # 多个基因可重复使用同一位点。仅当上界可能溢出时才用任意精度核算各 offset。
    if int(methyl.max()) * int(count.max()) <= np.iinfo(np.int64).max:
        methyl_sum = np.asarray(methyl @ weights).ravel()[present]
    else:
        entries = weights.tocoo()
        exact = np.zeros(len(offsets), dtype=object)
        np.add.at(exact, entries.col, methyl[entries.row].astype(object) * entries.data)
        if max(exact) > np.iinfo(np.int64).max:
            raise ValueError(f"{sid}: TSS methyl count exceeds int64")
        methyl_sum = exact[present].astype(np.int64)
    return pd.DataFrame({"Sample_ID": sid, "position_coarse": offsets[present], "meth_frac": means,
                         "count": methyl_sum, "CloneID": clone})


def tss_profile_frame(parts):
    columns = ["Sample_ID", "position_coarse", "meth_frac", "count", "CloneID"]
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=columns)


def chromosome_row(reference, offsets, row, sid, classes):
    """按实际 contig 核计数，未覆盖主染色体及 controls 保留 NaN 率。"""
    boundaries = np.searchsorted(row.indices, offsets)
    for j, (chrom, group) in enumerate(zip(reference.reference_seq_name, classes, strict=True)):
        lo, hi = boundaries[j], boundaries[j + 1]
        if lo == hi and group == "auxiliary":
            continue
        methyl, coverage = row.methyl[lo:hi], row.coverage[lo:hi]
        methyl_total, depth = int(methyl.sum()), int(coverage.sum())
        yield {"Sample_ID": sid, "chrom": str(chrom),
               "chromosome_class": group,
               "observed_sites": hi - lo, "methyl": methyl_total, "unmethyl": depth - methyl_total,
               "coverage": depth, "meth_frac": float(np.mean(methyl / coverage)) if lo < hi else np.nan}


def filter_reasons(frame, *, mapping_threshold=None, lambda_max=10.0, puc19_min=90.0,
                   cpg_min=500000, gini_max=0.5):
    """逐项解释原HQ门槛；未知值明确归为missing，观测零按实际阈值判定。"""
    for value in (lambda_max, puc19_min, cpg_min, gini_max):
        if not np.isfinite(value):
            raise ValueError("QC thresholds must be finite")
    if mapping_threshold is not None and not np.isfinite(mapping_threshold):
        raise ValueError("Mapping threshold must be finite or None")
    result = frame[["Sample_ID", "CloneID", "PlateID"]].copy()
    reasons = [[] for _ in range(len(frame))]
    criteria = [("Lambda", "Lambda_Meth_CpG_Rate%", lambda_max, "less"),
                ("pUC19", "pUC19_Meth_CpG_Rate%", puc19_min, "greater"),
                ("nCG", "Sample_Unique_CpG_Sites", cpg_min, "greater"),
                ("Gini", "Gini_Index", gini_max, "less")]
    if mapping_threshold is not None:
        criteria.append(("Mapping", "Native_Mapping%", mapping_threshold, "greater"))
    for label, column, threshold, direction in criteria:
        values = pd.to_numeric(frame[column], errors="raise").to_numpy(dtype=float)
        known = np.isfinite(values)
        passed = known & ((values < threshold) if direction == "less" else (values > threshold))
        result[f"QC_Pass_{label}"] = passed
        for i in np.flatnonzero(~passed):
            reasons[i].append(f"{label}_missing" if not known[i] else f"{label}_threshold")
    result["QC_Mapping_Evaluated"] = mapping_threshold is not None
    result["QC_Failed_Criteria"] = [len(x) for x in reasons]
    result["QC_Filter_Reasons"] = [";".join(x) if x else "pass" for x in reasons]
    result["HQ"] = result.QC_Failed_Criteria.eq(0)
    return result
