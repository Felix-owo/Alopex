"""可移植 H5AD/CSR：有界读取、原始整数计数和观测零/缺失语义。"""
from __future__ import annotations

import sys
_LOADED_MODULE_CODE = sys._getframe().f_code

import csv
from contextlib import closing, ExitStack
from decimal import Decimal, InvalidOperation
from itertools import islice
from pathlib import Path

import anndata as ad
from anndata._io.specs import read_elem, write_elem
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from numba import njit


class CSRReader:
    """直接切片磁盘CSR，不由AnnData将整个obsm或var标签载入内存。"""

    def __init__(self, path, group="X"):
        self.path = Path(path)
        self.handle = h5py.File(path, "r")
        try:
            self.raw_counts = group == "obsm/__values__"
            self.group = self.handle[group]
            if self.group.attrs.get("encoding-type") != "csr_matrix":
                raise ValueError("Expected a CSR H5AD matrix")
            self.n_obs, self.n_vars = map(int, self.group.attrs["shape"])
            self.obs = read_elem(self.handle["obs"])
            self.obs_names = self.obs.index
            if len(self.obs) != self.n_obs or not self.obs.index.is_unique:
                raise ValueError("CSR row identities are not unique or have a wrong length")
            self.indptr = self.group["indptr"][:]
            if (len(self.indptr) != self.n_obs + 1 or self.indptr[0] != 0
                    or np.any(np.diff(self.indptr) < 0)
                    or self.indptr[-1] != len(self.group["data"])
                    or len(self.group["data"]) != len(self.group["indices"])):
                raise ValueError("Invalid CSR pointers")
            self.uns = {key: read_elem(self.handle["uns"][key]) for key in
                        ("cpg_coordinate_base", "source_cpg_representation", "cpg_interval_convention", "reference_sequences", "qc_filter_parameters")
                        if key in self.handle["uns"]}
            self.reference = self.uns.get("reference_sequences")
        except BaseException:
            self.close()
            raise

    def read_rows(self, start, stop):
        if not 0 <= start <= stop <= self.n_obs:
            raise IndexError("CSR row range is outside the matrix")
        lo, hi = map(int, self.indptr[[start, stop]])
        indices, values = self.group["indices"][lo:hi], self.group["data"][lo:hi]
        if (indices < 0).any() or (indices >= self.n_vars).any():
            raise ValueError("CSR indices are outside the reference")
        value = sparse.csr_matrix((values, indices, self.indptr[start:stop + 1] - lo), shape=(stop - start, self.n_vars))
        if not value.has_canonical_format:
            raise ValueError("CSR indices must be sorted and unique")
        if self.raw_counts:
            if values.dtype != np.dtype('uint32'):
                raise ValueError("Raw CSR must store uint32 packed counts")
            coverage, methyl = values & 65535, values >> 16
            if (coverage == 0).any() or (methyl > coverage).any():
                raise ValueError("Raw CSR has invalid methyl/coverage counts")
        return value


    def iter_rows(self, size=4):
        for start in range(0, self.n_obs, size):
            stop = min(start + size, self.n_obs)
            yield self.read_rows(start, stop), start, stop

    def feature_names(self, start=0, stop=None):
        group = self.handle["var"]
        return group[group.attrs["_index"]].asstr()[start:stop]

    def close(self):
        self.handle.close()


class CSRWriter:
    """增量写一个由调用方负责事务提交的标准H5AD；显式零不会被消除。"""

    def __init__(self, path, obs, *, group="X", n_vars=0, dtype=np.float32, uns=None):
        ad.AnnData(obs=obs, uns=uns or {}).write_h5ad(path)
        self.handle = h5py.File(path, "r+")
        self.group = self.handle.create_group(group)
        self.group.attrs.update({"encoding-type": "csr_matrix", "encoding-version": "0.1.0",
                                 "shape": np.asarray([len(obs), n_vars], dtype=np.int64)})
        self.rows, self.nnz, self.n_vars = 0, 0, n_vars
        for name, data_type in (("data", dtype), ("indices", np.int64)):
            self.group.create_dataset(name, (0,), maxshape=(None,), chunks=(65536,), dtype=data_type,
                                      compression="gzip", compression_opts=1)
        self.group.create_dataset("indptr", data=np.zeros(len(obs) + 1, dtype=np.int64))
        self.var = None
        if group == "X":
            var = self.handle["var"]
            key = var.attrs["_index"]
            del var[key]
            self.var = var.create_dataset(key, (0,), maxshape=(None,), chunks=(65536,), dtype=h5py.string_dtype())
            self.var.attrs.update({"encoding-type": "string-array", "encoding-version": "0.2.0"})

    def append_features(self, names):
        end = self.n_vars + len(names)
        self.var.resize((end,))
        self.var[self.n_vars:end] = names
        self.n_vars = end
        shape = self.group.attrs["shape"]
        shape[1] = end
        self.group.attrs.modify("shape", shape)

    def append(self, value):
        if value.shape[1] != self.n_vars or self.rows + value.shape[0] >= len(self.group["indptr"]):
            raise ValueError("CSR block shape differs from the output")
        if not value.has_canonical_format:
            raise ValueError("CSR indices must be sorted and unique")
        stop = self.nnz + value.nnz
        for name, values in (("data", value.data), ("indices", value.indices)):
            target = self.group[name]
            target.resize((stop,))
            target[self.nnz:stop] = values
        # SciPy 可用 int32 表示单块，磁盘行指针累计的是所有细胞的观测项。
        pointers = np.add(value.indptr[1:], self.nnz, dtype=np.int64)
        self.group["indptr"][self.rows + 1:self.rows + value.shape[0] + 1] = pointers
        self.nnz, self.rows = stop, self.rows + value.shape[0]

    def close(self):
        self.handle.close()


@njit(cache=True, nogil=True)
def _parse_cpg_block(data, names, name_offsets, lengths, chromosome_offsets, previous_positions, base, previous_chrom):
    # 与公共四列契约相同的ASCII整数、u16 coverage及逐染色体唯一顺序校验。
    indices = np.empty(len(data) // 7 + 1, dtype=np.int64)
    values = np.empty(len(indices), dtype=np.uint32)
    cursor, used = 0, 0
    while cursor < len(data):
        start = cursor
        while cursor < len(data) and data[cursor] not in (9, 10):
            cursor += 1
        if cursor == start or cursor >= len(data) or data[cursor] != 9:
            return indices[:used], values[:used], previous_chrom, 1
        chrom = previous_chrom
        if chrom >= 0:
            if name_offsets[chrom + 1] - name_offsets[chrom] != cursor - start:
                chrom = -1
            else:
                for j in range(cursor - start):
                    if data[start + j] != names[name_offsets[chrom] + j]:
                        chrom = -1
                        break
        if chrom < 0:
            for ref in range(len(lengths)):
                if name_offsets[ref + 1] - name_offsets[ref] != cursor - start:
                    continue
                same = True
                for j in range(cursor - start):
                    if data[start + j] != names[name_offsets[ref] + j]:
                        same = False
                        break
                if same:
                    chrom = ref
                    break
        if chrom < 0:
            return indices[:used], values[:used], previous_chrom, 2
        numbers = np.zeros(3, dtype=np.int64)
        for column in range(3):
            cursor += 1
            digit_start, number = cursor, 0
            limit = lengths[chrom] - 1 + base if column == 0 else 65535
            while cursor < len(data) and data[cursor] not in (9, 10):
                digit = int(data[cursor]) - 48
                if digit < 0 or digit > 9:
                    return indices[:used], values[:used], previous_chrom, 1
                if number > (limit - digit) // 10:
                    return indices[:used], values[:used], previous_chrom, 2
                number = number * 10 + digit
                cursor += 1
            if cursor == digit_start or cursor >= len(data) or data[cursor] != (9 if column < 2 else 10):
                return indices[:used], values[:used], previous_chrom, 1
            numbers[column] = number
        position, methyl, unmethyl = numbers
        coverage = methyl + unmethyl
        if position < base or coverage <= 0 or coverage > 65535:
            return indices[:used], values[:used], previous_chrom, 2
        if position <= previous_positions[chrom]:
            return indices[:used], values[:used], previous_chrom, 3
        previous_positions[chrom], previous_chrom = position, chrom
        indices[used] = chromosome_offsets[chrom] + position
        values[used] = (methyl << 16) + coverage
        used += 1
        cursor += 1
    return indices[:used], values[:used], previous_chrom, 0


def read_cpg_counts(handle, path, reference, base):
    """严格解析单细胞原始四列流，保留整数与染色体内顺序，不生成稠密数组。"""
    if handle.readline().rstrip("\n") != "chrom\tpos\tmethyl\tunmethyl":
        raise ValueError(f"CpG input has an invalid header: {path}")
    encoded = [str(name).encode() for name in reference.reference_seq_name]
    names = np.frombuffer(b"".join(encoded), dtype=np.uint8)
    offsets = np.asarray([0, *np.cumsum([len(x) for x in encoded])], dtype=np.int64)
    lengths = reference.reference_bp_length.to_numpy(dtype=np.int64)
    chromosome_offsets = np.r_[0, np.cumsum(lengths + base)].astype(np.int64)
    previous_positions = np.full(len(lengths), -1, dtype=np.int64)
    parts, values, pending, previous_chrom, line = [], [], "", -1, 2
    while True:
        text = handle.read(4 * 1024 * 1024)
        if text:
            pending += text
            end = pending.rfind("\n")
            if end < 0:
                continue
            block, pending = pending[:end + 1], pending[end + 1:]
        elif pending:
            block, pending = pending + "\n", ""
        else:
            break
        indices, packed, previous_chrom, error = _parse_cpg_block(np.frombuffer(block.encode(), dtype=np.uint8),
            names, offsets, lengths, chromosome_offsets, previous_positions, base, previous_chrom)
        if error:
            detail = "duplicate or out-of-order CpG position" if error == 3 else "invalid CpG fields or coordinate/counts"
            raise ValueError(f"{path}:{line + len(indices)}: {detail}")
        parts.append(indices)
        values.append(packed)
        line += len(indices)
    indices = np.concatenate(parts) if parts else np.empty(0, dtype=np.int64)
    packed = np.concatenate(values) if values else np.empty(0, dtype=np.uint32)
    if len(indices) and np.any(np.diff(indices) <= 0):
        order = np.argsort(indices, kind="stable")
        indices, packed = indices[order], packed[order]
    return sparse.csr_matrix((packed, indices, np.asarray([0, len(indices)], dtype=np.int64)),
                             shape=(1, int(chromosome_offsets[-1])))


def cpg_reference(chrom_sizes, base):
    """Canonical 全参考列布局，保留所有 contig 和实际坐标基准。"""
    names = sorted(chrom_sizes)
    lengths = np.asarray([chrom_sizes[name] for name in names], dtype=np.uint64)
    return pd.DataFrame({"reference_seq_name": names, "reference_seq_length": lengths + base,
                         "reference_bp_length": lengths})


class IntervalIndex:
    """紧凑的数字区间索引，保留输入列序以及重叠区间的独立观测。"""

    def __init__(self, frame):
        self.size = len(frame)
        self.lookup = {}
        for chrom, positions in frame.groupby("chrom", sort=False).indices.items():
            order = np.argsort(frame.iloc[positions].start.to_numpy(), kind="stable")
            columns = np.asarray(positions, dtype=np.int64)[order]
            starts = frame.start.to_numpy(dtype=np.int64)[columns]
            ends = frame.end.to_numpy(dtype=np.int64)[columns]
            self.lookup[str(chrom)] = (starts, ends, np.maximum.accumulate(ends), columns)


@njit(cache=True, nogil=True)
def _region_observations(positions, packed, starts, ends, max_end, columns):
    count = 0
    for i in range(len(positions)):
        pos = positions[i]
        j = np.searchsorted(starts, pos, side="right") - 1
        while j >= 0 and max_end[j] > pos:
            if ends[j] > pos:
                count += 1
            j -= 1
    cols = np.empty(count, dtype=np.int64)
    fractions = np.empty(count, dtype=np.float64)
    methyls = np.empty(count, dtype=np.int64)
    coverages = np.empty(count, dtype=np.int64)
    used = 0
    for i in range(len(positions)):
        pos = positions[i]
        j = np.searchsorted(starts, pos, side="right") - 1
        while j >= 0 and max_end[j] > pos:
            if ends[j] > pos:
                coverage = int(packed[i] & 65535)
                methyl = int(packed[i] >> 16)
                cols[used], fractions[used] = columns[j], methyl / coverage
                methyls[used], coverages[used] = methyl, coverage
                used += 1
            j -= 1
    return cols, fractions, methyls, coverages



def aggregate_row(reader, row, index):
    """返回有观测区间的fraction/mean、methyl/cov整数和及位点数；保留真实零。"""
    reference = reader.reference
    lengths = reference.reference_seq_length.to_numpy(dtype=np.int64)
    offsets = np.r_[0, np.cumsum(lengths)]
    base = int(reader.uns["cpg_coordinate_base"])
    parts = []
    for k, chrom in enumerate(reference.reference_seq_name):
        lookup = index.lookup.get(str(chrom))
        if lookup is None:
            continue
        lo, hi = np.searchsorted(row.indices, offsets[k:k + 2])
        if lo < hi:
            parts.append(_region_observations(row.indices[lo:hi] - offsets[k] - base, row.data[lo:hi], *lookup))
    if not parts or not any(len(p[0]) for p in parts):
        return (np.empty(0, dtype=np.int64), np.empty(0, dtype=np.float32),
                np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64))
    columns, fractions, methyl, coverage = [np.concatenate([p[i] for p in parts]) for i in range(4)]
    order = np.argsort(columns, kind="stable")
    columns, fractions, methyl, coverage = [v[order] for v in (columns, fractions, methyl, coverage)]
    starts = np.r_[0, np.flatnonzero(np.diff(columns)) + 1]
    counts = np.diff(np.r_[starts, len(columns)])
    return (columns[starts], (np.add.reduceat(fractions, starts) / counts).astype(np.float32),
            np.add.reduceat(methyl, starts), np.add.reduceat(coverage, starts), counts)


def exact_integers(values, dtype):
    """从原文或已精确解析的整数列读取，保留有符号边界和整值十进制写法。"""
    bounds = np.iinfo(dtype)
    if values.dtype.kind not in "iu":
        if values.dtype.kind == "f":
            raise ValueError(f"{values.name} requires original integer text")

        def parse(text):
            if not isinstance(text, str) or not text.isascii() or "_" in text:
                raise ValueError(f"{values.name} must contain finite integers")
            try:
                number = Decimal(text)
            except InvalidOperation as exc:
                raise ValueError(f"{values.name} must contain finite integers") from exc
            if (not number.is_finite() or not bounds.min <= number <= bounds.max
                    or number != number.to_integral_value()):
                raise ValueError(f"{values.name} must contain finite integers in the {np.dtype(dtype)} range")
            return int(number)

        return values.map(parse).astype(np.int64)
    if len(values) and (int(values.min()) < bounds.min or int(values.max()) > bounds.max):
        raise ValueError(f"{values.name} is outside the {np.dtype(dtype)} range")
    return values.astype(np.int64, copy=False)


def bed_blocks(path, columns, integer_types):
    """按字面制表符读 BED。只有遇到非整数推断时才重读原文，避免浮点舍入。"""
    options = dict(sep="\t", header=None, names=range(max(columns) + 1), usecols=list(columns),
                   index_col=False, quoting=csv.QUOTE_NONE, keep_default_na=False, chunksize=100000)

    text_types = {i: str for i, name in columns.items() if name not in integer_types}
    with ExitStack() as handles:
        blocks = handles.enter_context(pd.read_csv(path, dtype=text_types, **options))
        raw_blocks, previous_raw = None, -1
        for index, block in enumerate(blocks):
            if any(block[i].dtype.kind not in "iu" for i, name in columns.items() if name in integer_types):
                # 第二个 reader 只取有歧义的块。两边使用相同的空行与分块规则。
                if raw_blocks is None:
                    raw_blocks = handles.enter_context(pd.read_csv(path, dtype=str, **options))
                block = next(islice(raw_blocks, index - previous_raw - 1, None))
                previous_raw = index
            block = block.rename(columns=columns, copy=False)
            chromosomes = block.chrom.unique()
            if any(not name.strip() for name in chromosomes):
                raise ValueError(f"BED chromosome name is empty: {path}")
            comments = [name for name in chromosomes if name.startswith("#")]
            if comments:
                block = block.loc[~block.chrom.isin(comments)].copy()
            for name, dtype in integer_types.items():
                try:
                    block[name] = exact_integers(block[name], dtype)
                except ValueError as exc:
                    raise ValueError(f"Invalid BED {name} in {path}: {exc}") from exc
            yield block


def write_single_cpg(reader, bed, output, obs):
    """逐块写参考标签、逐细胞写fraction CSR；不物化全库X或全参考字符串列表。"""
    base = int(reader.uns["cpg_coordinate_base"])
    writer = CSRWriter(output, obs.loc[list(reader.obs_names)], uns={
        "cpg_coordinate_base": base, "cpg_interval_convention": "closed" if base else "half-open",
        "source_cpg_representation": reader.uns.get("source_cpg_representation", "")})
    pieces, offset = {}, 0
    try:
        with closing(bed_blocks(bed, {0: "chrom", 1: "start", 2: "end"},
                                {"start": np.int64, "end": np.int64})) as chunks:
            for block in chunks:
                if (block.start < 0).any() or (block.end <= block.start).any():
                    raise ValueError("Invalid single-CpG BED interval")
                labels = np.asarray([f"{c}:{s + base}-{e}" for c, s, e in block.itertuples(index=False, name=None)], dtype=object)
                writer.append_features(labels)
                block.index = np.arange(offset, offset + len(block))
                for chrom, rows in block.groupby("chrom", sort=False):
                    pieces.setdefault(chrom, []).append((rows.start.to_numpy(), rows.end.to_numpy(), rows.index.to_numpy()))
                offset += len(block)
        if not offset:
            raise ValueError("Single-CpG BED is empty")
        index = IntervalIndex.__new__(IntervalIndex)
        index.size, index.lookup = offset, {}
        for chrom in list(pieces):
            starts, ends, columns = [np.concatenate([p[i] for p in pieces[chrom]]) for i in range(3)]
            del pieces[chrom]
            order = np.argsort(starts, kind="stable")
            starts, ends, columns = starts[order], ends[order], columns[order]
            index.lookup[str(chrom)] = (starts, ends, np.maximum.accumulate(ends), columns)
        for i in range(reader.n_obs):
            row = reader.read_rows(i, i + 1)
            columns, fractions, _, _, _ = aggregate_row(reader, row, index)
            writer.append(sparse.csr_matrix((fractions, columns, np.asarray([0, len(columns)], dtype=np.int64)), shape=(1, offset)))
    finally:
        writer.close()


def update_obs(path, obs, uns=None):
    """只更新标准obs，不读取或重写X。调用方持锁并先撤完成身份。"""
    with h5py.File(path, "r+") as handle:
        del handle["obs"]
        write_elem(handle, "obs", obs)
        for key, value in (uns or {}).items():
            if key in handle["uns"]:
                del handle["uns"][key]
            write_elem(handle["uns"], key, value)


@njit(cache=True, nogil=True)
def _mark_coverage(bits, indices):
    for index in indices:
        bits[index >> 3] |= np.uint8(1 << (index & 7))


@njit(cache=True, nogil=True)
def covered_intervals(bits, starts, ends):
    result = np.zeros(len(starts), dtype=np.bool_)
    for row in range(len(starts)):
        for position in range(starts[row], ends[row]):
            if bits[position >> 3] & (1 << (position & 7)):
                result[row] = True
                break
    return result


def covered_positions(reader, rows):
    """抽检覆盖并集：小输入一次排序，大输入使用有界位图；不扫描未选细胞。"""
    count = sum(int(reader.indptr[row + 1] - reader.indptr[row]) for row in rows)
    if count <= 2_000_000:
        parts = [reader.read_rows(int(row), int(row) + 1).indices for row in rows]
        return np.unique(np.concatenate(parts)) if parts else np.empty(0, dtype=np.int64)
    bits = np.zeros((reader.n_vars + 7) // 8, dtype=np.uint8)
    for row in rows:
        _mark_coverage(bits, reader.read_rows(int(row), int(row) + 1).indices)
    return bits
