"""逐染色体、TSS覆盖及逐细胞筛选解释；不以缺失补零或改变细胞集合。"""
from __future__ import annotations

import sys
_LOADED_MODULE_CODE = sys._getframe().f_code

import numpy as np
import pandas as pd
from scipy import sparse
from . import matrix

CHROMOSOME_INFO = "Chromosome_QC_Information.csv"
TSS_SCORE_INFO = "TSS_Score_Information.csv"
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
    eligible = {"eligible_bins": len(frame), "eligible_offsets": int(frame.offset.nunique())}
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
    return matrix.IntervalIndex(regions), membership, np.asarray(offsets, dtype=np.int64), eligible


def tss_profiles(reader, selected, bed, is_primary):
    """直接对抽中行聚合，重叠gene/bin权重保留；不复制Raw或生成临时TSS矩阵。"""
    rows, parts = reader.obs_names.get_indexer(selected.Sample_ID), []
    if (rows < 0).any():
        raise ValueError("TSS Sample_ID set differs from RawAdata")
    covered = matrix.covered_positions(reader, rows) if len(rows) < reader.n_obs else None
    index, membership, offsets, eligible = tss_intervals(bed, reader.reference, is_primary, covered,
                                                        int(reader.uns["cpg_coordinate_base"]))
    for row, (sid, clone) in zip(rows, selected[["Sample_ID", "CloneID"]].itertuples(index=False, name=None)):
        cols, frac, methyl, coverage, sites = matrix.aggregate_row(reader, reader.read_rows(int(row), int(row) + 1), index)
        if not len(cols):
            continue
        weights = membership[cols]
        count = np.asarray(weights.sum(axis=0)).ravel()
        present = count > 0
        means = np.asarray(frac.astype(np.float64) @ weights).ravel()[present] / count[present]
        methyl_sum = np.asarray(methyl @ weights).ravel()[present]
        coverage_sum = np.asarray(coverage @ weights).ravel()[present]
        sites_sum = np.asarray(sites @ weights).ravel()[present]
        parts.append(pd.DataFrame({"Sample_ID": sid, "position_coarse": offsets[present], "meth_frac": means,
            "count": methyl_sum, "CloneID": clone, "coverage": coverage_sum,
            "covered_bins": count[present], "observed_sites": sites_sum}))
    columns = ["Sample_ID", "position_coarse", "meth_frac", "count", "CloneID", "coverage", "covered_bins", "observed_sites"]
    profile = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=columns)
    return profile, eligible


def tss_scores(profile, full_qc, selected_ids, eligible, *, skipped=False):
    """TSS score为距离绝对值与甲基化率的Pearson相关；覆盖信息独立报告，不作分子加权。"""
    groups = dict(tuple(profile.groupby("Sample_ID", sort=False))) if len(profile) else {}
    rows = []
    for sid, clone in full_qc[["Sample_ID", "CloneID"]].itertuples(index=False, name=None):
        sampled = not skipped and sid in selected_ids
        row = {"Sample_ID": sid, "CloneID": clone, "TSS_sampled": sampled,
               "TSS_score": np.nan, "TSS_score_status": "skipped" if skipped else "not_sampled",
               "TSS_covered_offsets": np.nan, "TSS_covered_bins": np.nan,
               "TSS_eligible_offsets": np.nan, "TSS_eligible_bins": np.nan,
               "TSS_coverage_fraction": np.nan, "TSS_coverage": np.nan, "TSS_observed_sites": np.nan}
        if sampled:
            data = groups.get(sid)
            row.update(TSS_covered_offsets=0, TSS_covered_bins=0,
                       TSS_eligible_offsets=eligible["eligible_offsets"], TSS_eligible_bins=eligible["eligible_bins"],
                       TSS_coverage_fraction=0.0 if eligible["eligible_bins"] else np.nan,
                       TSS_coverage=0, TSS_observed_sites=0, TSS_score_status="no_coverage")
            if data is not None and len(data):
                x, y = np.abs(data.position_coarse.to_numpy(dtype=float)), data.meth_frac.to_numpy(dtype=float)
                status = "valid" if len(x) >= 2 and np.ptp(x) > 0 and np.ptp(y) > 0 else "undefined_correlation"
                score = float(np.corrcoef(x, y)[0, 1]) if status == "valid" else np.nan
                bins = int(data.covered_bins.sum())
                row.update(TSS_score=score, TSS_score_status=status, TSS_covered_offsets=len(data), TSS_covered_bins=bins,
                           TSS_coverage_fraction=bins / eligible["eligible_bins"] if eligible["eligible_bins"] else np.nan,
                           TSS_coverage=int(data.coverage.sum()), TSS_observed_sites=int(data.observed_sites.sum()))
        rows.append(row)
    return pd.DataFrame(rows)


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


def validate_tss_scores(scores, full_qc, selected_ids):
    """拒绝伪零score、错误抽样集合和无分母覆盖比例。"""
    try:
        if (len(scores) != len(full_qc) or scores.Sample_ID.duplicated().any()
                or set(scores.Sample_ID) != set(full_qc.Sample_ID)):
            return None
        scores = scores.set_index('Sample_ID').loc[full_qc.Sample_ID].reset_index()
        if not scores.CloneID.tolist() == full_qc.CloneID.tolist():
            return None
        flags = scores.TSS_sampled
        if not flags.isin([True, False]).all() or set(scores.loc[flags, 'Sample_ID']) != selected_ids:
            return None
        numeric = ['TSS_score', 'TSS_covered_offsets', 'TSS_covered_bins', 'TSS_eligible_offsets',
                   'TSS_eligible_bins', 'TSS_coverage_fraction', 'TSS_coverage', 'TSS_observed_sites']
        if not scores.loc[~flags, numeric].isna().all(axis=None):
            return None
        measured = scores.loc[flags]
        for name in ['TSS_covered_offsets', 'TSS_covered_bins', 'TSS_eligible_offsets', 'TSS_eligible_bins', 'TSS_coverage', 'TSS_observed_sites']:
            x = measured[name].to_numpy(dtype=float)
            if not np.isfinite(x).all() or (x < 0).any() or not np.equal(x, np.floor(x)).all():
                return None
        if (measured.TSS_covered_bins > measured.TSS_eligible_bins).any():
            return None
        valid = scores.TSS_score_status.eq('valid')
        if not np.isfinite(scores.loc[valid, 'TSS_score']).all() or scores.loc[valid, 'TSS_score'].abs().gt(1).any():
            return None
        if not scores.loc[~valid, 'TSS_score'].isna().all():
            return None
        expected = measured.TSS_covered_bins / measured.TSS_eligible_bins.replace(0, np.nan)
        if not np.allclose(measured.TSS_coverage_fraction, expected, rtol=0, atol=1e-14, equal_nan=True):
            return None
        return scores
    except (ValueError, KeyError, TypeError, AttributeError):
        return None
