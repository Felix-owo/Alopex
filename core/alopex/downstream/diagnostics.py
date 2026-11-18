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
    """一次读取gene/bin成员，先完整校验，再按实际主核参考构建数字索引。"""
    try:
        frame = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=[0, 1, 2, 6, 8],
                            converters={0: str, 6: str}, keep_default_na=False)
    except (ValueError, pd.errors.ParserError) as exc:
        raise ValueError(f"TSS BED must have at least 9 columns: {path}") from exc
    frame.columns = ["chrom", "start", "end", "gene", "offset"]
    for col in ("start", "end", "offset"):
        values = pd.to_numeric(frame[col], errors="coerce").to_numpy(dtype=float)
        if (not np.isfinite(values).all() or not np.equal(values, np.floor(values)).all()
                or (np.abs(values) >= (np.iinfo(np.int32).max if col == "offset" else np.iinfo(np.int64).max)).any()):
            raise ValueError(f"TSS BED {col} must contain finite integers: {path}")
        frame[col] = values.astype(np.int64)
    if (frame.start < 0).any() or (frame.end <= frame.start).any():
        raise ValueError(f"TSS BED must contain nonempty nonnegative intervals: {path}")
    if frame.gene.str.strip().eq("").any():
        raise ValueError(f"TSS BED contains an empty gene_id: {path}")
    if frame.duplicated(["chrom", "start", "end", "gene"]).any():
        raise ValueError(f"TSS BED contains duplicate gene/bin membership: {path}")
    lengths = dict(zip(reference.reference_seq_name, reference.reference_bp_length))
    frame = frame.loc[frame.chrom.map(is_primary) & frame.chrom.isin(lengths)].copy()
    frame = frame.loc[frame.start < frame.chrom.map(lengths)].copy()
    frame["end"] = np.minimum(frame.end, frame.chrom.map(lengths))
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


def tss_profiles(reader, selected, bed, is_primary):
    """直接对抽中行聚合，重叠gene/bin权重保留；不复制Raw或生成临时TSS矩阵。"""
    rows, parts = reader.obs_names.get_indexer(selected.Sample_ID), []
    if (rows < 0).any():
        raise ValueError("TSS Sample_ID set differs from RawAdata")
    covered = matrix.covered_positions(reader, rows) if len(rows) < reader.n_obs else None
    index, membership, offsets = tss_intervals(bed, reader.reference, is_primary, covered,
                                                        int(reader.uns["cpg_coordinate_base"]))
    for row, (sid, clone) in zip(rows, selected[["Sample_ID", "CloneID"]].itertuples(index=False, name=None)):
        cols, frac, methyl, _, _ = matrix.aggregate_row(reader, reader.read_rows(int(row), int(row) + 1), index)
        if not len(cols):
            continue
        weights = membership[cols]
        count = np.asarray(weights.sum(axis=0)).ravel()
        present = count > 0
        means = np.asarray(frac.astype(np.float64) @ weights).ravel()[present] / count[present]
        methyl_sum = np.asarray(methyl @ weights).ravel()[present]
        parts.append(pd.DataFrame({"Sample_ID": sid, "position_coarse": offsets[present], "meth_frac": means,
            "count": methyl_sum, "CloneID": clone}))
    columns = ["Sample_ID", "position_coarse", "meth_frac", "count", "CloneID"]
    profile = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=columns)
    return profile



def chromosome_summary(reader, is_primary, classify):
    """每个细胞的主染色体/controls均保留；未观测为零计数、NaN甲基化率。"""
    reference = reader.reference
    offsets = np.r_[0, np.cumsum(reference.reference_seq_length.to_numpy(dtype=np.int64))]
    parts = []
    for block, start, stop in reader.iter_rows(4):
        for row in range(stop - start):
            lo, hi = block.indptr[row:row + 2]
            indices, values = block.indices[lo:hi], block.data[lo:hi]
            for j, chrom in enumerate(reference.reference_seq_name):
                left, right = np.searchsorted(indices, offsets[j:j + 2])
                covered = values[left:right]
                labels = classify(str(chrom))
                if not len(covered) and not (is_primary(str(chrom)) or any(labels[:3])):
                    continue
                methyl = (covered >> 16).astype(np.uint64)
                coverage = (covered & 65535).astype(np.uint64)
                parts.append({"Sample_ID": str(reader.obs_names[start + row]), "chrom": str(chrom),
                    "chromosome_class": "host" if is_primary(str(chrom)) else "lambda" if labels[0] else "pUC19" if labels[1] else "mtDNA" if labels[2] else "auxiliary",
                    "observed_sites": len(covered), "methyl": int(methyl.sum()),
                    "unmethyl": int(coverage.sum() - methyl.sum()), "coverage": int(coverage.sum()),
                    "meth_frac": float(np.mean(methyl / coverage)) if len(covered) else np.nan})
    return pd.DataFrame(parts)


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
