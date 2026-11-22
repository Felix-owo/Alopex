"""可移植 H5AD/CSR：有界读取、原始整数计数和观测零/缺失语义。"""
from __future__ import annotations

import sys
_LOADED_MODULE_CODE = sys._getframe().f_code

import csv
from contextlib import closing, contextmanager, ExitStack
from decimal import Decimal, InvalidOperation
from io import TextIOBase
from itertools import islice
from pathlib import Path
from typing import NamedTuple

import anndata as ad
from anndata.experimental import read_elem, write_elem
import h5py
import numpy as np
import pandas as pd
from pandas.io.common import get_handle
from scipy import sparse
from numba import njit


class CSRReader:
    """直接切片磁盘CSR，不由AnnData将整个obsm或var标签载入内存。"""

    def __init__(self, path):
        self.path = Path(path)
        self.handle = h5py.File(path, "r")
        try:
            self.group = self.handle["X"]
            if (self.group.attrs.get("encoding-type") != "csr_matrix"
                    or self.group.attrs.get("encoding-version") != "0.1.0"):
                raise ValueError("Expected the AnnData CSR matrix encoding version 0.1.0")
            shape = np.asarray(self.group.attrs["shape"])
            if (shape.shape != (2,) or shape.dtype.kind not in "iu"
                    or int(shape.min()) < 0 or int(shape.max()) > np.iinfo(np.int64).max):
                raise ValueError("Invalid CSR shape")
            self.n_obs, self.n_vars = map(int, shape)
            var = self.handle["var"]
            features = var[var.attrs["_index"]]
            if not isinstance(features, h5py.Dataset) or features.shape != (self.n_vars,):
                raise ValueError("CSR feature index length differs from the matrix width")
            self.obs = read_elem(self.handle["obs"])
            self.obs_names = self.obs.index
            if len(self.obs) != self.n_obs or not self.obs.index.is_unique:
                raise ValueError("CSR row identities are not unique or have a wrong length")
            pointers, indices, data = (self.group[key] for key in ("indptr", "indices", "data"))
            if (pointers.shape != (self.n_obs + 1,) or pointers.dtype.kind not in "iu"
                    or indices.ndim != 1 or indices.dtype.kind not in "iu" or data.shape != indices.shape):
                raise ValueError("Invalid CSR index or data arrays")
            self.indptr = pointers[:]
            # 相邻比较避免无符号行指针回绕；先验整数类型，SciPy 才不会截断错误列号。
            if (self.indptr[0] != 0 or np.any(self.indptr[1:] < self.indptr[:-1])
                    or int(self.indptr[-1]) != len(data)):
                raise ValueError("Invalid CSR pointers")
            self.uns = {key: read_elem(self.handle["uns"][key]) for key in
                        ("cpg_coordinate_base", "source_cpg_representation", "cpg_interval_convention", "qc_filter_parameters")
                        if key in self.handle["uns"]}
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

    def __init__(self, path, obs, *, n_vars=0, uns=None):
        ad.AnnData(obs=obs, uns=uns or {}).write_h5ad(path)
        self.handle = h5py.File(path, "r+")
        try:
            self.group = self.handle.create_group("X")
            self.group.attrs.update({"encoding-type": "csr_matrix", "encoding-version": "0.1.0",
                                     "shape": np.asarray([len(obs), n_vars], dtype=np.int64)})
            self.rows, self.nnz, self.n_vars = 0, 0, n_vars
            for name, data_type in (("data", np.float32), ("indices", np.int64)):
                self.group.create_dataset(name, (0,), maxshape=(None,), chunks=(65536,), dtype=data_type,
                                          compression="gzip", compression_opts=1, shuffle=True)
            # 保持 dataset 打开，让 HDF5 缓存跨细胞的压缩块，关闭文件时统一刷盘。
            self.data, self.indices = self.group["data"], self.group["indices"]
            self.indptr = self.group.create_dataset("indptr", data=np.zeros(len(obs) + 1, dtype=np.int64))
            var = self.handle["var"]
            key = var.attrs["_index"]
            del var[key]
            self.var = var.create_dataset(key, (0,), maxshape=(None,), chunks=(65536,), dtype=h5py.string_dtype())
            self.var.attrs.update({"encoding-type": "string-array", "encoding-version": "0.2.0"})
        except BaseException:
            self.close()
            raise

    def append_features(self, names):
        end = self.n_vars + len(names)
        self.var.resize((end,))
        self.var[self.n_vars:end] = names
        self.n_vars = end
        shape = self.group.attrs["shape"]
        shape[1] = end
        self.group.attrs.modify("shape", shape)

    def append(self, value):
        if value.shape[1] != self.n_vars or self.rows + value.shape[0] >= len(self.indptr):
            raise ValueError("CSR block shape differs from the output")
        if not value.has_canonical_format:
            raise ValueError("CSR indices must be sorted and unique")
        stop = self.nnz + value.nnz
        for target, values in ((self.data, value.data), (self.indices, value.indices)):
            target.resize((stop,))
            target[self.nnz:stop] = values
        # SciPy 可用 int32 表示单块，磁盘行指针累计的是所有细胞的观测项。
        pointers = np.add(value.indptr[1:], self.nnz, dtype=np.int64)
        self.indptr[self.rows + 1:self.rows + value.shape[0] + 1] = pointers
        self.nnz, self.rows = stop, self.rows + value.shape[0]

    def close(self):
        self.handle.close()


class CpGCounts(NamedTuple):
    """一个细胞的有序观测，共用全参考索引，保留原始整数计数。"""

    indices: np.ndarray
    methyl: np.ndarray
    coverage: np.ndarray


@njit(cache=True, nogil=True)
def _parse_cpg_block(data, names, name_offsets, lengths, chromosome_offsets, previous_positions,
                     base, previous_chrom, remaining, observations):
    numbers = np.empty(3, dtype=np.int64)
    cursor, used = 0, 0
    while cursor < len(data):
        start = cursor
        while cursor < len(data) and data[cursor] not in (9, 10):
            cursor += 1
        if cursor == start or cursor >= len(data) or data[cursor] != 9:
            return used, previous_chrom, remaining, 1
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
            return used, previous_chrom, remaining, 2
        for column in range(3):
            cursor += 1
            digit_start, number = cursor, 0
            limit = lengths[chrom] - 1 + base if column == 0 else (1 << 53) - 1
            while cursor < len(data) and data[cursor] not in (9, 10):
                digit = int(data[cursor]) - 48
                if digit < 0 or digit > 9:
                    return used, previous_chrom, remaining, 1
                if number > (limit - digit) // 10:
                    return used, previous_chrom, remaining, 2
                number = number * 10 + digit
                cursor += 1
            if cursor == digit_start or cursor >= len(data) or data[cursor] != (9 if column < 2 else 10):
                return used, previous_chrom, remaining, 1
            numbers[column] = number
        position, methyl, unmethyl = numbers
        coverage = methyl + unmethyl
        if position < base or coverage <= 0:
            return used, previous_chrom, remaining, 2
        if position <= previous_positions[chrom]:
            return used, previous_chrom, remaining, 3
        # 覆盖均非负。整行可表示时，染色体、QC 分组、Gini 和 1 Mb 窗口的各个子和均可表示。
        if coverage > remaining:
            return used, previous_chrom, remaining, 4
        remaining -= coverage
        previous_positions[chrom], previous_chrom = position, chrom
        observations[0, used] = chromosome_offsets[chrom] + position
        observations[1, used] = methyl
        observations[2, used] = coverage
        used += 1
        cursor += 1
    return used, previous_chrom, remaining, 0


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
    parts, pending, previous_chrom, line = [], "", -1, 2
    remaining = np.iinfo(np.int64).max
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
        observations = np.empty((3, block.count("\n")), dtype=np.int64)
        used, previous_chrom, remaining, error = _parse_cpg_block(
            np.frombuffer(block.encode(), dtype=np.uint8), names, offsets, lengths, chromosome_offsets,
            previous_positions, base, previous_chrom, remaining, observations)
        if error:
            detail = {3: "duplicate or out-of-order CpG position", 4: "cell CpG coverage exceeds int64"}.get(
                error, "invalid CpG fields or coordinate/counts")
            raise ValueError(f"{path}:{line + used}: {detail}")
        parts.append(observations)
        line += used
    if not parts:
        return CpGCounts(*np.empty((3, 0), dtype=np.int64))
    observations = parts[0] if len(parts) == 1 else np.concatenate(parts, axis=1)
    indices = observations[0]
    if np.any(indices[1:] < indices[:-1]):
        observations = observations[:, np.argsort(indices, kind="stable")]
    return CpGCounts(*observations)


def cpg_reference(chrom_sizes, base):
    """Canonical 全参考列布局，保留所有 contig 和实际坐标基准。"""
    names = sorted(chrom_sizes)
    lengths = np.asarray([chrom_sizes[name] for name in names], dtype=np.uint64)
    return pd.DataFrame({"reference_seq_name": names, "reference_seq_length": lengths + base,
                         "reference_bp_length": lengths})


def interval_lookup(starts, ends, columns):
    """保留参考列身份；仅对尚未按起点排列的区间作稳定排序。"""
    if np.any(starts[1:] < starts[:-1]):
        order = np.argsort(starts, kind="stable")
        starts, ends, columns = starts[order], ends[order], columns[order]
    max_end = np.maximum.accumulate(ends) if np.any(ends[1:] < ends[:-1]) else ends
    return starts, ends, max_end, columns


@njit(cache=True, nogil=True)
def _region_observations(positions, methyl, coverage, starts, ends, max_end, columns):
    last = np.searchsorted(starts, positions, side="right") - 1
    count = 0
    for i in range(len(positions)):
        pos = positions[i]
        j = last[i]
        while j >= 0 and max_end[j] > pos:
            if ends[j] > pos:
                count += 1
            j -= 1
    cols = np.empty(count, dtype=np.int64)
    fractions = np.empty(count, dtype=np.float64)
    methyls = np.empty(count, dtype=np.int64)
    used = 0
    for i in range(len(positions)):
        pos = positions[i]
        j = last[i]
        while j >= 0 and max_end[j] > pos:
            if ends[j] > pos:
                cols[used], fractions[used] = columns[j], methyl[i] / coverage[i]
                methyls[used] = methyl[i]
                used += 1
            j -= 1
    return cols, fractions, methyls


def aggregate_row(reader, row, index):
    """返回有观测区间的列号、逐胞嘧啶平均 fraction 和 methyl 整数和，保留真实零。"""
    reference = reader.reference
    lengths = reference.reference_seq_length.to_numpy(dtype=np.int64)
    offsets = np.r_[0, np.cumsum(lengths)]
    boundaries = np.searchsorted(row.indices, offsets)
    base = int(reader.uns["cpg_coordinate_base"])
    parts = []
    for k, chrom in enumerate(reference.reference_seq_name):
        lookup = index.get(str(chrom))
        if lookup is None:
            continue
        lo, hi = boundaries[k], boundaries[k + 1]
        if lo < hi:
            matched = _region_observations(row.indices[lo:hi] - offsets[k] - base,
                                           row.methyl[lo:hi], row.coverage[lo:hi], *lookup)
            if len(matched[0]):
                parts.append(matched)
    if not parts:
        return (np.empty(0, dtype=np.int64), np.empty(0, dtype=np.float64),
                np.empty(0, dtype=np.int64))
    # 大多数 BED 的列按染色体成段。交错或嵌套的列仍按原列号稳定归并。
    parts.sort(key=lambda part: part[0][0])
    columns, fractions, methyl = (np.concatenate(items) for items in zip(*parts, strict=True))
    if np.any(columns[1:] < columns[:-1]):
        order = np.argsort(columns, kind="stable")
        columns, fractions, methyl = columns[order], fractions[order], methyl[order]
    starts = np.r_[0, np.flatnonzero(columns[1:] != columns[:-1]) + 1]
    counts = np.diff(np.r_[starts, len(columns)])
    return (columns[starts], np.add.reduceat(fractions, starts) / counts,
            np.add.reduceat(methyl, starts))


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


class _BedText(TextIOBase):
    """在 pandas 的 C 解析器截短 NUL 字段前检查解压后的原文。"""

    def __init__(self, source, path):
        self.source = source
        self.path = path

    def read(self, size=-1):
        text = self.source.read(size)
        if "\x00" in text:
            raise ValueError(f"BED contains a NUL byte: {self.path}")
        return text


def bed_blocks(path, columns, integer_types):
    """按字面制表符读 BED。只有遇到非整数推断时才重读原文，避免浮点舍入。"""
    options = dict(sep="\t", header=None, names=range(max(columns) + 1), usecols=list(columns),
                   index_col=False, quoting=csv.QUOTE_NONE, keep_default_na=False, chunksize=100000)

    text_types = {i: str for i, name in columns.items() if name not in integer_types}
    with ExitStack() as handles:
        def open_blocks(dtype):
            source = handles.enter_context(get_handle(path, "r", encoding="utf-8", compression="infer"))
            return handles.enter_context(pd.read_csv(_BedText(source.handle, path), dtype=dtype, **options))

        blocks = open_blocks(text_types)
        raw_blocks, previous_raw = None, -1
        for index, block in enumerate(blocks):
            if any(block[i].dtype.kind not in "iu" for i, name in columns.items() if name in integer_types):
                # 第二个 reader 只取有歧义的块。两边使用相同的空行与分块规则。
                if raw_blocks is None:
                    raw_blocks = open_blocks(str)
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


@contextmanager
def single_cpg_writer(reader, bed, output, host_lengths):
    """按 BED 原列序写出参考标签，提供只读区间索引，由调用线程顺序提交稀疏行。"""
    base = int(reader.uns["cpg_coordinate_base"])
    writer = CSRWriter(output, pd.DataFrame(index=reader.obs_names), uns={
        "cpg_coordinate_base": base, "cpg_interval_convention": "closed" if base else "half-open",
        "source_cpg_representation": reader.uns["source_cpg_representation"]})
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
                    if chrom not in host_lengths or rows.end.max() > host_lengths[chrom]:
                        raise ValueError(f"Single-CpG BED interval outside the host reference: {chrom} in {bed}")
                    pieces.setdefault(chrom, []).append((rows.start.to_numpy(), rows.end.to_numpy(), rows.index.to_numpy()))
                offset += len(block)
        if not offset:
            raise ValueError("Single-CpG BED is empty")
        index = {}
        for chrom in list(pieces):
            starts, ends, columns = [np.concatenate([p[i] for p in pieces[chrom]]) for i in range(3)]
            del pieces[chrom]
            lookup = interval_lookup(starts, ends, columns)
            starts, ends = lookup[:2]
            ties = starts[1:] == starts[:-1]
            if ties.any():
                tied = np.r_[False, ties] | np.r_[ties, False]
                if pd.MultiIndex.from_arrays((starts[tied], ends[tied])).has_duplicates:
                    raise ValueError(f"Duplicate Single-CpG BED interval on {chrom}: {bed}")
            index[chrom] = lookup
        yield writer, index
    finally:
        writer.close()


def single_cpg_row(reader, row, index, n_vars):
    """将单细胞的区间均值编码为可按原序提交的 float32 稀疏行。"""
    columns, fractions, _ = aggregate_row(reader, row, index)
    return sparse.csr_matrix((fractions, columns, np.asarray([0, len(columns)], dtype=np.int64)),
                             shape=(1, n_vars), dtype=np.float32)


def update_obs(path, obs, uns=None):
    """只更新标准obs，不读取或重写X。调用方持锁并先撤完成身份。"""
    with h5py.File(path, "r+") as handle:
        del handle["obs"]
        write_elem(handle, "obs", obs)
        for key, value in (uns or {}).items():
            if key in handle["uns"]:
                del handle["uns"][key]
            write_elem(handle["uns"], key, value)
