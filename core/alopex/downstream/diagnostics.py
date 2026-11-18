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


def tss_intervals(path, reference, is_primary, covered=None, base=0):
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
    frame = frame.loc[frame.chrom.map(is_primary) & (frame.start < limits)].copy()
    if covered is not None:
        offset = np.r_[0, np.cumsum(reference.reference_seq_length.to_numpy(dtype=np.int64))]
        reference_offsets = dict(zip(reference.reference_seq_name, offset[:-1]))
        starts = frame.chrom.map(reference_offsets).to_numpy(dtype=np.int64) + frame.start.to_numpy(dtype=np.int64) + base
        ends = frame.chrom.map(reference_offsets).to_numpy(dtype=np.int64) + frame.end.to_numpy(dtype=np.int64) + base
        if covered.dtype == np.uint8:
            keep = matrix.covered_intervals(covered, starts, ends)
        else:
            locations = np.searchsorted(covered, starts)
            keep = locations < len(covered)
            keep[keep] &= covered[locations[keep]] < ends[keep]
        frame = frame.loc[keep].copy()
    regions = frame[["chrom", "start", "end"]].drop_duplicates().reset_index(drop=True)
    regions["bin"] = np.arange(len(regions))
    members = frame.merge(regions, on=["chrom", "start", "end"], how="left", validate="many_to_one")
    codes, offsets = pd.factorize(members.offset, sort=True)
    membership = sparse.csr_matrix((np.ones(len(members), dtype=np.int64), (members.bin, codes)),
                                   shape=(len(regions), len(offsets)))
    return matrix.IntervalIndex(regions), membership, np.asarray(offsets, dtype=np.int64)


def tss_profile_row(reader, row, sid, clone, tss_index):
    index, membership, offsets = tss_index
    cols, frac, methyl, _, _ = matrix.aggregate_row(reader, row, index)
    if not len(cols):
        return None
    weights = membership[cols]
    count = np.asarray(weights.sum(axis=0)).ravel()
    present = count > 0
    means = np.asarray(frac.astype(np.float64) @ weights).ravel()[present] / count[present]
    methyl_sum = np.asarray(methyl @ weights).ravel()[present]
    return pd.DataFrame({"Sample_ID": sid, "position_coarse": offsets[present], "meth_frac": means,
                         "count": methyl_sum, "CloneID": clone})


def tss_profile_frame(parts):
    columns = ["Sample_ID", "position_coarse", "meth_frac", "count", "CloneID"]
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=columns)


def chromosome_row(reference, offsets, row, sid, is_primary, classify):
    """按实际 contig 核计数，未覆盖主染色体及 controls 保留 NaN 率。"""
    for j, chrom in enumerate(reference.reference_seq_name):
        lo, hi = np.searchsorted(row.indices, offsets[j:j + 2])
        covered = row.data[lo:hi]
        is_lambda, is_puc, is_mt, _ = classify(str(chrom))
        primary = is_primary(str(chrom))
        if not len(covered) and not (primary or is_lambda or is_puc or is_mt):
            continue
        methyl = (covered >> 16).astype(np.uint64)
        coverage = (covered & 65535).astype(np.uint64)
        methyl_total, depth = int(methyl.sum()), int(coverage.sum())
        yield {"Sample_ID": sid, "chrom": str(chrom),
               "chromosome_class": ("host" if primary else "lambda" if is_lambda else
                                    "pUC19" if is_puc else "mtDNA" if is_mt else "auxiliary"),
               "observed_sites": len(covered), "methyl": methyl_total, "unmethyl": depth - methyl_total,
               "coverage": depth, "meth_frac": float(np.mean(methyl / coverage)) if len(covered) else np.nan}


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
