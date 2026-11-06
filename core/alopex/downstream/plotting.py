"""准备好数据后的QC绘图、样式与标签；不读取sealed交付或提交科学缓存。"""
from __future__ import annotations

import logging
import math
import os
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


matplotlib = mcolors = plt = ticker = sns = Line2D = PCA = StandardScaler = None


_PLOT_IMPORT_ERROR = None


def _load_plot_dependencies() -> None:
    global matplotlib, mcolors, plt, ticker, sns, Line2D, PCA, StandardScaler
    global _PLOT_IMPORT_ERROR

    if plt is not None:
        return
    if _PLOT_IMPORT_ERROR is not None:
        raise RuntimeError("Downstream plotting dependencies are unavailable") from _PLOT_IMPORT_ERROR

    logging.getLogger("matplotlib.font_manager").setLevel(logging.WARNING)
    logging.getLogger("fontTools").setLevel(logging.ERROR)
    try:
        import matplotlib as _matplotlib
        import matplotlib.colors as _mcolors
        import matplotlib.pyplot as _plt
        import matplotlib.ticker as _ticker
        import seaborn as _sns
        from matplotlib.lines import Line2D as _Line2D
        from sklearn.decomposition import PCA as _PCA
        from sklearn.preprocessing import StandardScaler as _StandardScaler
    except ImportError as exc:
        _PLOT_IMPORT_ERROR = exc
        raise RuntimeError("Downstream plotting dependencies are unavailable") from exc

    matplotlib = _matplotlib
    mcolors = _mcolors
    plt = _plt
    ticker = _ticker
    sns = _sns
    Line2D = _Line2D
    PCA = _PCA
    StandardScaler = _StandardScaler


class Colors:
    TITLE = "#414956"
    SUBTITLE = "#525B69"
    LABEL = "#646D7A"
    TICK = "#7A8390"
    ANNOT = "#8B93A0"
    SPINE = "#B8C0CB"
    REFLINE = "#D5DBE3"
    MEDIAN = "#CF718C"


NATURE_BASE_PALETTE = [
    "#4FA29C", "#CF718C", "#6F90BC", "#C9A45A",
    "#7EA271", "#9A83BE", "#D08B55", "#8B9097",
    "#498A86", "#C15473", "#567BAD", "#B99244",
    "#6D8E61", "#8368AD", "#C1773D", "#787D84",
    "#74B3AF", "#D69BAC", "#96ABC8", "#CFB787",
    "#9CB494", "#B5A7CC", "#D4A784", "#A7AAAF",
]


COMPOSITION_PALETTE = [
    "#6F8197", "#C49772", "#8B9A84", "#A88998", "#D1A2A6"
]


class FontSizes:
    TICK = 13
    LABEL = 14.5
    TITLE = 15.5
    LEGEND = 11.5
    LEGEND_TITLE = 12.5
    PANEL_LETTER = 18


class PlotParams:
    DPI = 250
    SAVE_DPI = 300
    SPINE_WIDTH = 0.95
    SCATTER_SIZE = 40
    SCATTER_ALPHA = 0.7
    LINE_WIDTH = 2.0


def setup_nature_style() -> None:
    """应用全部 DNAme QC 图共用的绘图主题。"""
    _load_plot_dependencies()
    sns.set_theme(
        style="white",
        context="talk",
        font="sans-serif",
        font_scale=1.14,
        rc={
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans", "sans-serif"],
            "figure.dpi": PlotParams.DPI,
            "savefig.dpi": PlotParams.SAVE_DPI,
            "axes.linewidth": PlotParams.SPINE_WIDTH,
            "xtick.major.width": 0.9,
            "ytick.major.width": 0.9,
            "xtick.major.size": 4.0,
            "ytick.major.size": 4.0,
            "xtick.labelsize": FontSizes.TICK,
            "ytick.labelsize": FontSizes.TICK,
            "axes.labelsize": FontSizes.LABEL,
            "axes.titlesize": FontSizes.TITLE,
            "legend.fontsize": FontSizes.LEGEND,
            "legend.title_fontsize": FontSizes.LEGEND_TITLE,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        },
    )


def format_ax(
    ax: matplotlib.axes.Axes,
    title: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
) -> None:
    if title:
        ax.set_title(title, fontsize=FontSizes.TITLE, pad=15, weight="bold", color=Colors.SUBTITLE)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=FontSizes.LABEL, labelpad=9, color=Colors.LABEL)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=FontSizes.LABEL, labelpad=9, color=Colors.LABEL)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(Colors.SPINE)
    ax.spines["bottom"].set_color(Colors.SPINE)
    ax.spines["left"].set_linewidth(PlotParams.SPINE_WIDTH)
    ax.spines["bottom"].set_linewidth(PlotParams.SPINE_WIDTH)
    ax.grid(False)
    ax.tick_params(axis="both", which="major", labelsize=FontSizes.TICK, colors=Colors.TICK)


def add_panel_letter(ax: matplotlib.axes.Axes, letter: str) -> None:
    ax.text(
        -0.13, 1.065, letter,
        transform=ax.transAxes,
        fontsize=FontSizes.PANEL_LETTER,
        fontweight="bold",
        ha="left",
        va="top",
        color=Colors.SUBTITLE,
        alpha=0.8,
    )


def panel_label(index: int) -> str:
    """生成不限数量的 Excel 式面板标签：A..Z、AA..AZ、BA…"""
    if index < 0:
        raise ValueError("panel index must be non-negative")
    value = index + 1
    chars = []
    while value:
        value, remainder = divmod(value - 1, 26)
        chars.append(chr(ord("A") + remainder))
    return "".join(reversed(chars))


def add_right_legend(
    fig: plt.Figure,
    group_levels: list[str],
    color_map: dict[str, str],
    title: str = "CloneID",
    marker: str = "o",
    linewidth: float = 0.6,
) -> None:
    legend_handles = []
    for grp in group_levels:
        if marker == "o":
            handle = Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                label=str(grp),
                markerfacecolor=color_map.get(grp, "#8B9097"),
                markeredgecolor="white",
                markeredgewidth=linewidth,
                markersize=8.8,
                alpha=0.8,
            )
        else:
            handle = Line2D(
                [0],
                [0],
                color=color_map.get(grp, "#8B9097"),
                linewidth=linewidth,
                alpha=0.8,
                label=str(grp),
            )
        legend_handles.append(handle)

    legend = fig.legend(
        handles=legend_handles,
        labels=[str(g) for g in group_levels],
        title=title,
        loc="center left",
        bbox_to_anchor=(0.925, 0.5),
        frameon=False,
        fontsize=FontSizes.LEGEND,
        title_fontsize=FontSizes.LEGEND_TITLE,
    )
    legend.get_title().set_color(Colors.SUBTITLE)
    for text in legend.get_texts():
        text.set_color(Colors.TICK)


def _normalize_label_value(value: Any) -> str:
    if value is None:
        return ""
    return str(value)


def _common_suffix(labels: list[str]) -> str:
    if not labels:
        return ""
    reversed_labels = [label[::-1] for label in labels]
    reversed_suffix = os.path.commonprefix(reversed_labels)
    return reversed_suffix[::-1]


def get_diff_only_labels(
    labels: list[Any],
    ellipsis: str = "…",
    protect_single_label: bool = True,
) -> tuple[list[str], dict[str, str]]:
    clean_labels = [_normalize_label_value(label) for label in labels]

    if protect_single_label and len(clean_labels) <= 1:
        return clean_labels, {"prefix": "", "suffix": ""}

    common_prefix = os.path.commonprefix(clean_labels)
    common_suffix = _common_suffix(clean_labels)

    min_len = min((len(label) for label in clean_labels), default=0)
    if len(common_prefix) + len(common_suffix) >= min_len:
        common_suffix = ""
        if len(common_prefix) >= min_len:
            common_prefix = ""

    display_labels = []
    used = set()

    for original in clean_labels:
        start = len(common_prefix)
        end = len(original) - len(common_suffix) if common_suffix else len(original)
        diff = original[start:end]
        display = diff if diff else ellipsis


        if display in used:
            display = f"{display} ({original})"

        display_labels.append(display)
        used.add(display)

    return display_labels, {"prefix": common_prefix, "suffix": common_suffix}


def get_diff_only_label_map(
    labels: list[Any],
    ellipsis: str = "…",
) -> tuple[dict[str, str], dict[str, str]]:
    unique_labels = []
    for label in labels:
        clean_label = _normalize_label_value(label)
        if clean_label not in unique_labels:
            unique_labels.append(clean_label)

    display_labels, omitted = get_diff_only_labels(unique_labels, ellipsis=ellipsis)
    return dict(zip(unique_labels, display_labels, strict=True)), omitted


UNKNOWN_COLOR = "#C7CDD6"


def normalize_clone_id_value(value: Any) -> Any:
    if pd.isna(value):
        return np.nan
    clean_value = str(value).strip()
    if clean_value == "":
        return np.nan
    return clean_value


def make_cloneid_diff_labels(
    values: list[Any],
    ellipsis: str = "…",
) -> tuple[dict[str, str], dict[str, str]]:
    clean_values = []
    for value in values:
        clean_value = normalize_clone_id_value(value)
        if isinstance(clean_value, str) and clean_value not in clean_values:
            clean_values.append(clean_value)
    return get_diff_only_label_map(clean_values, ellipsis=ellipsis)


def cloneid_omission_note(omitted: dict[str, str]) -> str:
    parts = []
    if omitted.get("prefix"):
        parts.append(f"prefix '{omitted['prefix']}'")
    if omitted.get("suffix"):
        parts.append(f"suffix '{omitted['suffix']}'")
    return " and ".join(parts)


def add_cloneid_omission_note(
    fig: plt.Figure,
    omitted: dict[str, str],
    *,
    y: float,
) -> None:
    """以统一样式标注 CloneID 显示名中被省略的公共部分。"""

    note = cloneid_omission_note(omitted)
    if note:
        fig.text(
            0.5,
            y,
            f'Note: common CloneID {note} omitted',
            ha='center',
            fontsize=10,
            color=Colors.ANNOT,
            style='italic',
        )


def _panel_omission_note(
    sample_labels_shortened: bool,
    clone_note: str,
    *,
    include_clone_note: bool = True,
) -> str:
    if sample_labels_shortened and clone_note and include_clone_note:
        return (
            f"Note: common CloneID {clone_note} omitted; "
            "Sample_ID labels shortened within each CloneID panel by omitting shared text"
        )
    if sample_labels_shortened:
        return "Note: Sample_ID labels shortened within each CloneID panel by omitting shared text"
    if clone_note and include_clone_note:
        return f"Note: common CloneID {clone_note} omitted"
    return ""


def generate_clone_color_map(unique_clones: list[str]) -> dict[str, str]:
    n_clones = len(unique_clones)
    if n_clones <= 2:
        global_palette = ["#E78AAE", "#63C2C0"]
    else:
        _load_plot_dependencies()
        global_palette = (
            NATURE_BASE_PALETTE[:n_clones]
            if n_clones <= len(NATURE_BASE_PALETTE)
            else sns.color_palette("Set2", n_clones)
        )
    return {unique_clones[i]: global_palette[i] for i in range(n_clones)}


def plot_wgbs_qc_panel(
    df: pd.DataFrame,
    batch_col: str = 'CloneID',
    clone_color_map: dict = None,
    out_path: Path = None,
    figsize: tuple = (16.5, 15.5),
):
    """绘制 WGBS 对照与样本甲基化核心 QC 面板。"""
    if df.empty:
        return
    _load_plot_dependencies()

    plot_df = df.copy()
    if batch_col not in plot_df.columns:
        plot_df[batch_col] = 'All'
    if batch_col == 'CloneID':
        plot_df[batch_col] = plot_df[batch_col].map(normalize_clone_id_value)

    col_map = {
        'Lambda_Meth_CpG_Rate%': 'lambda_rate',
        'pUC19_Meth_CpG_Rate%': 'pUC19_rate',
        'Sample_Meth_CpG_Rate%': 'sample_rate',
        'Lambda_Unique_CpG_Sites': 'lambda_nCG',
        'pUC19_Unique_CpG_Sites': 'pUC19_nCG',
        'Sample_Unique_CpG_Sites': 'sample_nCG',
    }
    rename_dict = {key: value for key, value in col_map.items() if key in plot_df.columns}
    plot_df = plot_df.rename(columns=rename_dict)
    plot_df[batch_col] = plot_df[batch_col].astype('category')
    panel_cfg = [
        ('sample_nCG', 'sample_rate', 'Sample Methylation vs CpG Coverage'),
        ('lambda_rate', 'pUC19_rate', 'Bisulfite Conversion Controls'),
        ('lambda_rate', 'sample_rate', 'Sample Methylation vs Lambda Control'),
        ('pUC19_rate', 'sample_rate', 'Sample Methylation vs pUC19 Control'),
        ('lambda_nCG', 'lambda_rate', 'Lambda Coverage vs Methylation'),
        ('pUC19_nCG', 'pUC19_rate', 'pUC19 Coverage vs Methylation'),
    ]
    available_cols = set(plot_df.columns)
    panel_cfg = [
        (x_col, y_col, title)
        for x_col, y_col, title in panel_cfg
        if x_col in available_cols and y_col in available_cols
    ]
    if not panel_cfg:
        print('Not enough columns to plot WGBS QC panel.')
        return

    n_panels = len(panel_cfg)
    nrows = math.ceil(n_panels / 2)
    group_levels = list(plot_df[batch_col].cat.categories)
    if batch_col == 'CloneID':
        label_map, omitted = make_cloneid_diff_labels(group_levels)
        legend_levels = [label_map.get(g, g) for g in group_levels]
        color_map = (
            {grp: clone_color_map.get(grp, UNKNOWN_COLOR) for grp in group_levels}
            if clone_color_map
            else generate_clone_color_map(group_levels)
        )
        legend_color_map = {label_map.get(g, g): color_map[g] for g in group_levels}
    else:
        label_map, omitted = {}, {'prefix': '', 'suffix': ''}
        legend_levels = group_levels
        color_map = (
            {grp: clone_color_map.get(grp, UNKNOWN_COLOR) for grp in group_levels}
            if clone_color_map
            else generate_clone_color_map(group_levels)
        )
        legend_color_map = color_map

    fig, axes = plt.subplots(nrows, 2, figsize=(figsize[0], 5 * nrows), dpi=PlotParams.DPI)
    if nrows == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    def pretty_label(value):
        labels = {
            'sample_nCG': 'Sample nCpG',
            'sample_rate': 'Sample Rate',
            'lambda_rate': 'Lambda Rate',
            'pUC19_rate': 'pUC19 Rate',
            'lambda_nCG': 'Lambda nCpG',
            'pUC19_nCG': 'pUC19 nCpG',
        }
        return labels.get(value, value.replace('_', ' ').title())

    for i, (x_col, y_col, title) in enumerate(panel_cfg):
        ax = axes[i]
        for group in group_levels:
            subset = plot_df[plot_df[batch_col] == group]
            ax.scatter(
                subset[x_col],
                subset[y_col],
                s=PlotParams.SCATTER_SIZE,
                alpha=PlotParams.SCATTER_ALPHA,
                color=color_map[group],
                edgecolor='white',
                linewidth=0.55,
                zorder=3,
            )
        x_med = plot_df[x_col].median()
        y_med = plot_df[y_col].median()
        ax.axvline(x_med, ls='--', lw=0.95, color=Colors.REFLINE, alpha=0.8, zorder=0)
        ax.axhline(y_med, ls='--', lw=0.95, color=Colors.REFLINE, alpha=0.8, zorder=0)
        ax.text(
            0.02,
            0.975,
            f'{pretty_label(x_col)} Median: {x_med:.3f}\n'
            f'{pretty_label(y_col)} Median: {y_med:.3f}',
            transform=ax.transAxes,
            ha='left',
            va='top',
            fontsize=10.2,
            color=Colors.ANNOT,
        )
        format_ax(
            ax,
            title=title,
            xlabel=pretty_label(x_col),
            ylabel=pretty_label(y_col),
        )
        add_panel_letter(ax, panel_label(i))
        if x_col in ['sample_nCG', 'lambda_nCG'] and (plot_df[x_col] > 0).all():
            if plot_df[x_col].max() / max(plot_df[x_col].min(), 1) > 20:
                ax.set_xscale('log')
        if 'rate' in y_col:
            finite_y = pd.to_numeric(plot_df[y_col], errors='coerce')
            finite_y = finite_y[np.isfinite(finite_y)]
            if not finite_y.empty:
                ymin = finite_y.min()
                ymax = finite_y.max()
                pad = max((ymax - ymin) * 0.1, 0.006)
                ax.set_ylim(ymin - pad, ymax + pad)

    for index in range(len(panel_cfg), len(axes)):
        axes[index].axis('off')
    add_right_legend(fig, legend_levels, legend_color_map, title=batch_col)
    add_cloneid_omission_note(fig, omitted, y=0.015)
    plt.tight_layout(rect=[0, 0.03, 0.885, 0.958], w_pad=2.2, h_pad=2.4)
    if out_path:
        fig.savefig(out_path, bbox_inches='tight')
    plt.show()
    plt.close(fig)


def normalize_clone_id_list(clone_ids: Any) -> list[str]:
    if clone_ids is None:
        return []

    if isinstance(clone_ids, str):
        candidate_values = [clone_ids]
    else:
        try:
            candidate_values = list(clone_ids)
        except TypeError:
            candidate_values = [clone_ids]

    clean_values = []
    for value in candidate_values:
        clean_value = normalize_clone_id_value(value)
        if pd.isna(clean_value):
            continue
        clean_value = str(clean_value)
        if clean_value not in clean_values:
            clean_values.append(clean_value)
    return clean_values


def make_clone_selection_suffix(clone_ids: Any) -> str:
    clone_list = normalize_clone_id_list(clone_ids)
    if not clone_list:
        return 'all'
    if len(clone_list) == 1:
        return clone_list[0].replace('/', '_')
    return f'{len(clone_list)}clones'


def normalize_metric_list(metrics: Any) -> list[str]:
    if metrics is None:
        return []
    if isinstance(metrics, str):
        return [metrics]
    try:
        return [str(metric) for metric in metrics]
    except TypeError:
        return [str(metrics)]


def _format_metric_stat(value) -> str:
    if pd.isna(value):
        return 'NA'
    if value >= 1000000.0:
        return f'{value / 1000000.0:.2f}M'
    if value >= 1000.0:
        return f'{value / 1000.0:.2f}K'
    if value >= 1:
        return f'{value:.2f}'
    return f'{value:.3f}'


def _final_qc_histogram_bins(values, use_log):
    scale = 'log' if use_log and np.all(values > 0) else 'log1p' if use_log else 'linear'
    transformed = (np.log10(values) if scale == 'log' else
                   np.sign(values) * np.log1p(np.abs(values)) if scale == 'log1p' else values)
    edges = np.histogram_bin_edges(transformed, bins=100)
    edges = (np.power(10.0, edges) if scale == 'log' else
             np.sign(edges) * np.expm1(np.abs(edges)) if scale == 'log1p' else edges)
    if values.min() != values.max():
        edges[0], edges[-1] = values.min(), values.max()
    return edges, scale


def _draw_final_qc_histogram(ax, counts, edges, scale, color, fill_alpha=0.32):
    transformed = (np.log10(edges) if scale == 'log' else
                   np.sign(edges) * np.log1p(np.abs(edges)) if scale == 'log1p' else edges)
    offsets = np.arange(-6, 7, dtype=float)
    kernel = np.exp(-0.5 * (offsets / 1.5) ** 2)
    kernel /= kernel.sum()
    smoothed = np.convolve(np.pad(np.asarray(counts, dtype=float), 6, mode='symmetric'), kernel, mode='valid')
    positions = np.linspace(0, len(counts), len(counts) * 8 + 1)
    curve_x = np.interp(positions, np.arange(len(edges)), transformed)
    curve_x = (np.power(10.0, curve_x) if scale == 'log' else
               np.sign(curve_x) * np.expm1(np.abs(curve_x)) if scale == 'log1p' else curve_x)
    curve_y = np.interp(positions, np.arange(len(counts)) + 0.5, smoothed)
    ax.stairs(counts, edges, fill=True, color=color, alpha=fill_alpha / 2, linewidth=0)
    ax.fill_between(curve_x, curve_y, color=color, alpha=fill_alpha, linewidth=0)
    ax.plot(curve_x, curve_y, color=color, linewidth=1.6)
    if scale == 'log':
        ax.set_xscale('log')
    elif scale == 'log1p':
        ax.set_xscale('function', functions=(
            lambda x: np.sign(x) * np.log1p(np.abs(x)),
            lambda x: np.sign(x) * np.expm1(np.abs(x)),
        ))
    if scale == 'log':
        ax.xaxis.set_major_locator(ticker.LogLocator(base=10, subs=(1, 2, 5), numticks=7))
        ax.xaxis.set_minor_formatter(ticker.NullFormatter())
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(
            lambda x, _: f'{x / 1e6:g}M' if x >= 1e6 else f'{x / 1e3:g}K' if x >= 1e3 else f'{x:g}'
        ))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))


def plot_final_qc_dashboard(
    df: pd.DataFrame,
    group_col: str = 'CloneID',
    clone_id: Any = None,
    mode: str = 'distribution',
    metric: Any = 'Native_Mapping%',
    clone_color_map: dict = None,
    out_path: Path = None,
    figsize: tuple = (16.5, 12),
):
    """绘制所选指标：distribution 叠加各组，per_clone_samples 按组分面。

    distribution 支持单指标或有序列表；per_clone_samples 的多指标由 batch 逐项绘制。
    位点数使用对数横轴，真实零值保留并使用 signed-log1p；不填补缺失值。
    """
    if mode not in ('distribution', 'per_clone_samples'):
        raise ValueError("Unsupported mode. Use 'distribution' or 'per_clone_samples'.")
    if df.empty or group_col not in df.columns:
        print(f"Error: column '{group_col}' not found.")
        return
    _load_plot_dependencies()

    plot_df = df.copy()
    if group_col == 'CloneID':
        plot_df[group_col] = plot_df[group_col].map(normalize_clone_id_value)

    plot_df = plot_df.dropna(subset=[group_col]).copy()
    plot_df[group_col] = plot_df[group_col].astype(str)

    all_group_ids = sorted(plot_df[group_col].dropna().astype(str).unique().tolist())

    selected_clone_ids = normalize_clone_id_list(clone_id)

    if selected_clone_ids:
        missing = [cid for cid in selected_clone_ids if cid not in all_group_ids]
        if missing:
            print("Error: some clone_id values were not found:")
            for value in missing:
                print(f"  - {value}")
            print("Available CloneID values:")
            for value in all_group_ids:
                print(f"  - {value}")
            return

        plot_df = plot_df.loc[plot_df[group_col].isin(selected_clone_ids)].copy()
        group_ids = [cid for cid in all_group_ids if cid in selected_clone_ids]
    else:
        group_ids = all_group_ids

    if plot_df.empty:
        print("Error: no data available after CloneID filtering.")
        return

    metric_configs = {
        'Sample_Unique_CpG_Sites': ('Unique CpG Sites', True),
        'Gini_Index': ('Library Uniformity (Gini Index)', False),
        'Non_CpG_Methylation%': ('Non-CpG Methylation (%)', False),
        'Native_Mapping%': ('Native Mapping (%)', False),
        'Final_Pair_Yield%': ('Final Pair Yield (%)', False),
        'Duplicate_Pair_Rate%': ('Duplicate Pair Rate (%)', False),
        'Sample_Meth_CpG_Rate%': ('Host CpG Methylation (%)', False),
        'Lambda_Unique_CpG_Sites': ('Lambda Unique CpG Sites', True),
        'Lambda_Meth_CpG_Rate%': ('Lambda CpG Methylation (%)', False),
        'pUC19_Unique_CpG_Sites': ('pUC19 Unique CpG Sites', True),
        'pUC19_Meth_CpG_Rate%': ('pUC19 CpG Methylation (%)', False),
        'High_CpH_Flag_Rate%': ('High-CpH Flag Rate (%)', False),
        'Trim_Retention%': ('Trim Retention (%)', False),
        'Signal_Host_Rate%': ('Host CpG Signal Share (%)', False),
        'Signal_Lambda_Rate%': ('Lambda CpG Signal Share (%)', False),
        'Signal_pUC19_Rate%': ('pUC19 CpG Signal Share (%)', False),
        'Signal_mtDNA_Rate%': ('mtDNA CpG Signal Share (%)', False),
    }
    valid_metrics = normalize_metric_list(metric)
    if not valid_metrics:
        raise ValueError('No QC metric selected; specify a metric or a nonempty list.')
    unsupported = [name for name in valid_metrics if name not in metric_configs]
    if unsupported:
        raise ValueError(f'Unsupported QC metrics: {unsupported}. Supported: {list(metric_configs)}')
    missing_metrics = [name for name in valid_metrics if name not in plot_df.columns]
    if missing_metrics:
        raise ValueError(f'QC metrics not found in dataframe: {missing_metrics}')
    if mode == 'per_clone_samples' and len(valid_metrics) != 1:
        raise ValueError('Use plot_final_qc_dashboard_batch for multiple per_clone_samples metrics.')
    metric = valid_metrics[0]

    if mode == 'distribution':
        label_map, omitted = make_cloneid_diff_labels(group_ids)
        plot_df[f'{group_col}_Display'] = plot_df[group_col].map(label_map).fillna(plot_df[group_col])
        plot_df[f'{group_col}_Display'] = plot_df[f'{group_col}_Display'].astype('category')

        display_levels = [label_map.get(g, g) for g in group_ids]
        if clone_color_map:
            color_map = {
                label_map.get(orig, orig): clone_color_map.get(orig, UNKNOWN_COLOR)
                for orig in group_ids
            }
        else:
            color_map = generate_clone_color_map(display_levels)

        n_cols = min(2, len(valid_metrics))
        n_rows = math.ceil(len(valid_metrics) / n_cols)
        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=(figsize[0] if n_cols == 2 else figsize[0] * 0.68, 4.9 * n_rows),
            dpi=PlotParams.DPI,
        )
        axes = np.array(axes).reshape(-1)

        for (i, metric_name) in enumerate(valid_metrics):
            ax = axes[i]
            display_name, use_log = metric_configs[metric_name]

            xvals = pd.to_numeric(plot_df[metric_name], errors='coerce')
            panel_df = plot_df.loc[np.isfinite(xvals), [f'{group_col}_Display', metric_name]].copy()
            panel_df[metric_name] = pd.to_numeric(panel_df[metric_name], errors='coerce')
            panel_df = panel_df.dropna(subset=[metric_name])

            if panel_df.empty:
                format_ax(ax, title=display_name, xlabel=display_name, ylabel='Counts')
                ax.text(
                    0.5, 0.5, f'No finite data (n = 0 / {len(plot_df)})', transform=ax.transAxes,
                    ha='center', va='center', fontsize=12, color=Colors.ANNOT
                )
                continue

            bin_edges, scale = _final_qc_histogram_bins(panel_df[metric_name].to_numpy(dtype=float), use_log)
            for display_group in display_levels:
                sub = panel_df[panel_df[f'{group_col}_Display'] == display_group]
                if sub.empty:
                    continue
                vals = sub[metric_name].dropna()
                if len(vals) == 0:
                    continue
                counts, _ = np.histogram(vals, bins=bin_edges)
                # 组数多时降低填充遮挡；分箱、计数和平滑算法与分面模式完全相同。
                _draw_final_qc_histogram(
                    ax, counts, bin_edges, scale, color_map.get(display_group, UNKNOWN_COLOR),
                    fill_alpha=0.32 * min(1.0, math.sqrt(2 / len(group_ids))),
                )

            median_val = panel_df[metric_name].median()
            ax.axvline(median_val, ls='--', lw=1.0, color=Colors.MEDIAN, alpha=0.9, zorder=5)
            ax.text(
                0.02, 0.975,
                f'All groups: n = {len(panel_df)} / {len(plot_df)} | Median: {_format_metric_stat(median_val)}',
                transform=ax.transAxes, ha='left', va='top',
                fontsize=10.2, color=Colors.ANNOT,
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.9, pad=2), zorder=6,
            )

            format_ax(
                ax,
                title=display_name,
                xlabel=f"{display_name}{' (log scale)' if scale == 'log' else ' (log1p scale)' if scale == 'log1p' else ''}",
                ylabel='Counts',
            )
            add_panel_letter(ax, panel_label(i))
            ax.set_xlim(bin_edges[0], bin_edges[-1])

            if not use_log:
                ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
                ax.ticklabel_format(style='plain', axis='x')

        for j in range(len(valid_metrics), len(axes)):
            axes[j].axis('off')

        legend_title = group_col if not selected_clone_ids else f'{group_col} selected'
        add_right_legend(fig, display_levels, color_map, title=legend_title, marker='line', linewidth=2.2)

        if selected_clone_ids:
            display_summary = ', '.join(display_levels[:3]) + (' ...' if len(display_levels) > 3 else '')
            fig.suptitle(
                f'Final QC Dashboard - {display_summary}',
                fontsize=FontSizes.TITLE + 3,
                weight='bold',
                color=Colors.TITLE,
                y=0.995,
            )

        if not selected_clone_ids:
            add_cloneid_omission_note(fig, omitted, y=0.015)

        plt.tight_layout(rect=[0, 0.03, 0.885, 0.97], w_pad=2.2, h_pad=2.4)

        if out_path:
            fig.savefig(out_path, bbox_inches='tight')
        plt.show()
        plt.close(fig)
        return

    display_name, use_log = metric_configs[metric]
    group_totals = plot_df.groupby(group_col, observed=True).size()
    plot_df[metric] = pd.to_numeric(plot_df[metric], errors='coerce')
    plot_df = plot_df.loc[np.isfinite(plot_df[metric])].copy()
    bin_edges, scale = (_final_qc_histogram_bins(plot_df[metric].to_numpy(dtype=float), use_log)
                        if not plot_df.empty else (None, 'linear'))
    clones = group_ids
    clone_display_map, omitted = make_cloneid_diff_labels(clones)

    if clone_color_map is None:
        clone_color_map = generate_clone_color_map(clones)

    n_cols = min(2, len(clones))
    n_rows = max(1, math.ceil(len(clones) / n_cols))
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(figsize[0] if n_cols == 2 else figsize[0] * 0.68, 4.9 * n_rows),
        dpi=PlotParams.DPI,
    )
    axes = np.array(axes).reshape(-1)
    for i, clone in enumerate(clones):
        ax = axes[i]
        values = plot_df.loc[plot_df[group_col] == clone, metric]

        if values.empty:
            format_ax(
                ax, title=f"{clone_display_map.get(clone, clone)} - {display_name}",
                xlabel=display_name, ylabel='Counts',
            )
            ax.text(0.5, 0.5, f'No finite data (n = 0 / {group_totals[clone]})',
                    transform=ax.transAxes, ha='center', va='center', color=Colors.ANNOT)
            continue

        counts, _ = np.histogram(values.to_numpy(dtype=float), bins=bin_edges)
        _draw_final_qc_histogram(ax, counts, bin_edges, scale, clone_color_map.get(clone, UNKNOWN_COLOR))

        median_val = values.median()
        ax.axvline(median_val, ls='--', lw=1.0, color=Colors.MEDIAN, alpha=0.9, zorder=4)
        ax.text(
            0.02, 0.975,
            f'n = {len(values)} / {group_totals[clone]} | Median: {_format_metric_stat(median_val)}',
            transform=ax.transAxes, ha='left', va='top',
            fontsize=10.2, color=Colors.ANNOT,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.9, pad=2), zorder=6,
        )

        format_ax(
            ax,
            title=f"{clone_display_map.get(clone, clone)} - {display_name}",
            xlabel=f"{display_name}{' (log scale)' if scale == 'log' else ' (log1p scale)' if scale == 'log1p' else ''}",
            ylabel='Counts',
        )
        ax.set_xlim(bin_edges[0], bin_edges[-1])
        ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        add_panel_letter(ax, panel_label(i))
        if not use_log:
            ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
            ax.ticklabel_format(style='plain', axis='x')

    for j in range(len(clones), len(axes)):
        axes[j].axis('off')

    note_text = cloneid_omission_note(omitted) if not selected_clone_ids else ''

    if note_text:
        fig.text(0.5, 0.01, note_text, ha='center', fontsize=10, color=Colors.ANNOT, style='italic')

    plt.tight_layout(rect=[0, 0.03, 1, 0.98])
    if out_path:
        fig.savefig(out_path, bbox_inches='tight')
    plt.show()
    plt.close(fig)


def plot_plate_metric_map(
    df: pd.DataFrame,
    metric: str = 'Sample_Unique_CpG_Sites',
    clone_color_map: dict = None,
    out_path: Path = None,
):
    """按 CloneID 分面绘制一个 QC 指标的孔板位置热图。"""
    if df.empty or 'PlateID' not in df.columns or 'CloneID' not in df.columns:
        return
    _load_plot_dependencies()

    row_order = list('ABCDEFGH')
    row_map = {letter: index for index, letter in enumerate(row_order)}
    plot_df = df.copy()
    plot_df['CloneID'] = plot_df['CloneID'].map(normalize_clone_id_value)
    unique_for_labels = sorted(plot_df['CloneID'].dropna().unique().tolist())
    clone_display_map, omitted = make_cloneid_diff_labels(unique_for_labels)
    plate_id = plot_df['PlateID'].astype('string').str.strip().str.upper()
    if plate_id.fillna('').eq('').all():
        print('当前细胞无 PlateID，跳过孔板图。')
        return
    valid_plate = plate_id.str.fullmatch(r'[A-H](?:[1-9]|1[0-2])', na=False)
    if (~valid_plate).any():
        invalid = plot_df.loc[~valid_plate]
        examples = []
        for idx, row in invalid.head(10).iterrows():
            sample = str(row.get('Sample_ID', idx))
            value = row.get('PlateID')
            value_label = '<empty>' if pd.isna(value) or not str(value).strip() else str(value)
            examples.append(f'{sample}={value_label}')
        remainder = len(invalid) - len(examples)
        suffix = f'; 另有 {remainder} 行' if remainder > 0 else ''
        print(
            f"⚠️ 忽略 {len(invalid)} 行无效 PlateID（应为 A1-H12）："
            f"{'; '.join(examples)}{suffix}"
        )
    plot_df = plot_df.loc[valid_plate].copy()
    if plot_df.empty:
        print('⚠️ 没有可用于 Plate View 的有效 PlateID')
        return
    plot_df['PlateID'] = plate_id.loc[valid_plate]
    plot_df['row_idx'] = plot_df['PlateID'].str[0].map(row_map)
    plot_df['col_idx'] = plot_df['PlateID'].str[1:].astype(int)
    unique_clones = sorted(plot_df['CloneID'].dropna().unique().tolist())
    if not unique_clones:
        print('⚠️ 没有同时包含有效 PlateID 和 CloneID 的记录；请重新生成 DNAme QC 表')
        return

    num_clones = len(unique_clones)
    n_cols = 2
    n_rows = math.ceil(num_clones / n_cols)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(n_cols * 9, n_rows * 6),
        dpi=PlotParams.DPI,
    )
    if num_clones == 1:
        axes = np.array([axes])
    axes_flat = axes.flatten()

    for index, clone in enumerate(unique_clones):
        ax = axes_flat[index]
        subset = plot_df[plot_df['CloneID'] == clone]
        clone_title = clone_display_map.get(clone, clone)
        ax.set_facecolor('white')
        ax.grid(
            True,
            which='major',
            color='#EAEAEA',
            linestyle='--',
            linewidth=0.8,
            zorder=0,
        )
        if clone_color_map and clone in clone_color_map:
            base_color = clone_color_map[clone]
            cmap = mcolors.LinearSegmentedColormap.from_list(
                f'custom_{clone}', ['#F5F5F5', base_color]
            )
        else:
            cmap = 'YlOrBr'
        ax.scatter(
            subset['col_idx'],
            subset['row_idx'],
            s=850,
            c=subset[metric] if metric in subset.columns else np.zeros(len(subset)),
            cmap=cmap,
            edgecolors='white',
            linewidth=0.8,
            alpha=0.85,
            zorder=2,
        )
        if metric in subset.columns:
            for _, row in subset.iterrows():
                val = row[metric]
                if pd.isna(val):
                    continue
                if val >= 1000000.0:
                    label = f'{int(val / 1000000.0)}M'
                elif val >= 1000.0:
                    label = f'{int(val / 1000.0)}K'
                else:
                    label = f'{int(val)}'
                ax.text(
                    row['col_idx'],
                    row['row_idx'],
                    label,
                    ha='center',
                    va='center',
                    fontsize=9,
                    color=Colors.TITLE,
                    fontweight='normal',
                    zorder=3,
                )
        ax.set_xticks(range(1, 13))
        ax.set_yticks(range(8))
        ax.set_yticklabels(row_order)
        ax.invert_yaxis()
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color(Colors.SPINE)
            spine.set_linewidth(1.5)
        ax.tick_params(
            axis='both',
            which='major',
            color=Colors.TICK,
            labelcolor=Colors.LABEL,
            labelsize=11,
            pad=6,
        )
        metric_title = metric.replace('_', ' ').title()
        ax.set_title(
            f'{clone_title} - {metric_title}',
            fontsize=14,
            pad=15,
            weight='bold',
            color=Colors.SUBTITLE,
        )
    for index in range(num_clones, len(axes_flat)):
        axes_flat[index].axis('off')
    add_cloneid_omission_note(fig, omitted, y=0.01)
    plt.tight_layout(rect=[0, 0.03, 1, 1])
    if out_path:
        fig.savefig(out_path, bbox_inches='tight')
    plt.show()
    plt.close(fig)


def plot_cpg_signal_composition(
    df: pd.DataFrame,
    out_path: Path = None,
    simplify_sample_ids: bool = True,
):
    """绘制 host、lambda、pUC19 与 mtDNA 的 CpG signal 比例。"""
    if df.empty or 'Sample_ID' not in df.columns or 'CloneID' not in df.columns:
        return

    cols = ['Signal_Host_Rate%', 'Signal_Lambda_Rate%', 'Signal_pUC19_Rate%', 'Signal_mtDNA_Rate%']
    valid_cols = [c for c in cols if c in df.columns]
    if not valid_cols:
        return
    _load_plot_dependencies()

    plot_df = df.copy()
    plot_df['CloneID'] = plot_df['CloneID'].map(normalize_clone_id_value)
    plot_df = plot_df.dropna(subset=['CloneID', 'Sample_ID'])

    colors = COMPOSITION_PALETTE[:len(valid_cols)]
    clones = sorted(plot_df['CloneID'].dropna().unique().tolist())
    if len(clones) == 0:
        return

    clone_display_map, omitted = make_cloneid_diff_labels(clones)
    n_cols = 2
    n_rows = math.ceil(len(clones) / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 4.8 * n_rows), dpi=PlotParams.DPI)
    axes = np.array(axes).reshape(-1)

    sample_label_note_needed = False

    for (i, clone) in enumerate(clones):
        ax = axes[i]
        sub = plot_df.loc[plot_df['CloneID'] == clone, ['Sample_ID', *valid_cols]].copy()

        if sub.empty:
            ax.axis('off')
            continue

        sub['Sample_ID'] = sub['Sample_ID'].astype(str)

        if simplify_sample_ids:
            sample_label_map, sample_omitted = get_diff_only_label_map(sub['Sample_ID'].tolist())
            sub['Sample_ID_Display'] = sub['Sample_ID'].map(sample_label_map)
            if sample_omitted.get('prefix') or sample_omitted.get('suffix'):
                sample_label_note_needed = True
        else:
            sub['Sample_ID_Display'] = sub['Sample_ID']

        sub = sub.set_index('Sample_ID_Display')[valid_cols]
        sub.plot(
            kind='bar',
            stacked=True,
            ax=ax,
            width=0.85,
            color=colors,
            alpha=0.85,
            edgecolor='white',
            linewidth=0.4,
            legend=False,
        )

        ax.set_title(
            f'{clone_display_map.get(clone, clone)} - CpG Signal Composition',
            weight='bold',
            color=Colors.SUBTITLE,
            pad=12,
        )
        ax.set_xlabel('')
        ax.set_ylabel('Percentage of CpG Signals (%)', color=Colors.LABEL, labelpad=10)

        n_samples = sub.shape[0]
        if n_samples > 60:
            ax.set_xticklabels([])
        else:
            ax.tick_params(axis='x', labelrotation=90, labelsize=5.8)
            for tick_label in ax.get_xticklabels():
                tick_label.set_color(Colors.TICK)

        ax.tick_params(axis='x', length=0, colors=Colors.TICK, labelcolor=Colors.TICK)
        ax.tick_params(axis='y', colors=Colors.TICK, labelcolor=Colors.TICK)

        sns.despine(ax=ax)
        ax.spines['left'].set_color(Colors.SPINE)
        ax.spines['bottom'].set_color(Colors.SPINE)
        ax.spines['left'].set_linewidth(1.2)
        ax.spines['bottom'].set_linewidth(1.2)

    for j in range(len(clones), len(axes)):
        axes[j].axis('off')

    legend_labels = [c.replace('Reads_', '').replace('_Rate', '').replace('_', ' ') for c in valid_cols]
    handles = [
        plt.Rectangle((0, 0), 1, 1, fc=colors[k], ec='white', lw=0.4, alpha=0.85)
        for k in range(len(valid_cols))
    ]
    legend = fig.legend(
        handles=handles,
        labels=legend_labels,
        loc='center left',
        bbox_to_anchor=(0.92, 0.5),
        frameon=False,
        fontsize=FontSizes.LEGEND,
    )
    for text in legend.get_texts():
        text.set_color(Colors.TICK)

    note_text = _panel_omission_note(
        sample_label_note_needed, cloneid_omission_note(omitted)
    )

    if note_text:
        fig.text(
            0.5,
            0.01,
            note_text,
            ha='center',
            fontsize=10,
            color=Colors.ANNOT,
            style='italic',
        )

    plt.tight_layout(rect=[0, 0.03, 0.9, 0.98])
    if out_path:
        fig.savefig(out_path, bbox_inches='tight')
    plt.show()
    plt.close(fig)


def plot_tss_profile(
    df_tss: pd.DataFrame | str | Path,
    clone_color_map: dict = None,
    out_path: Path = None,
):
    """绘制 TSS 甲基化比例与信号数，保留 CSV 身份和浮点精度；无可绘数据时清理旧图。"""
    if isinstance(df_tss, (str, Path)):
        df_tss = pd.read_csv(
            df_tss, converters=dict.fromkeys(("Sample_ID", "CloneID"), str),
            float_precision="round_trip",
        )
    required = {'position_coarse', 'meth_frac', 'CloneID'}
    if df_tss.empty or not required.issubset(df_tss.columns):
        if out_path is not None:
            Path(out_path).unlink(missing_ok=True)
        return
    _load_plot_dependencies()

    plot_df = df_tss.copy()
    plot_df['CloneID'] = plot_df['CloneID'].map(normalize_clone_id_value)
    group_levels = sorted(plot_df['CloneID'].dropna().astype(str).unique().tolist())
    if not group_levels:
        if out_path is not None:
            Path(out_path).unlink(missing_ok=True)
        print('No valid CloneID values: TSS profile skipped.')
        return
    clone_display_map, omitted = make_cloneid_diff_labels(group_levels)
    if clone_color_map is None:
        clone_color_map = generate_clone_color_map(group_levels)

    num_clones = len(group_levels)
    n_rows = num_clones
    fig, axes = plt.subplots(
        n_rows,
        2,
        figsize=(16.5, max(5.5, n_rows * 5.2)),
        dpi=PlotParams.DPI,
    )
    if num_clones == 1:
        axes = np.array([axes])
    for index, clone in enumerate(group_levels):
        subset = plot_df[plot_df['CloneID'] == clone]
        clone_title = clone_display_map.get(clone, clone)
        ax1, ax2 = axes[index]
        sns.lineplot(
            data=subset,
            x='position_coarse',
            y='meth_frac',
            color=clone_color_map.get(clone, None),
            ax=ax1,
            linewidth=PlotParams.LINE_WIDTH,
            legend=False,
        )
        format_ax(
            ax1,
            title=f'{clone_title} - TSS Methylation Profile',
            xlabel='Distance to TSS (bp)',
            ylabel='mCpG Fraction',
        )
        if 'count' in subset.columns:
            sns.lineplot(
                data=subset,
                x='position_coarse',
                y='count',
                color=clone_color_map.get(clone, None),
                ax=ax2,
                linewidth=PlotParams.LINE_WIDTH,
                legend=False,
            )
            format_ax(
                ax2,
                title=f'{clone_title} - TSS mCpG Signal Profile',
                xlabel='Distance to TSS (bp)',
                ylabel='mCpG Signal Sum',
            )
        else:
            ax2.axis('off')
    add_cloneid_omission_note(fig, omitted, y=0.01)
    plt.tight_layout(rect=[0, 0.03, 1, 1])
    if out_path:
        plt.savefig(out_path, bbox_inches='tight')
    plt.show()
    plt.close(fig)


CPG_SLOPE_RATIO_MIN_CHANGE = 0.05


def _robust_inlier_mask(
    values: np.ndarray,
    z_thresh: float = 3.5,
    min_n: int = 6,
) -> np.ndarray:

    values = np.asarray(values, dtype=float)
    mask = np.isfinite(values)
    if mask.sum() < min_n:
        return mask

    finite_values = values[mask]
    median = float(np.median(finite_values))
    deviations = np.abs(finite_values - median)
    mad = float(np.median(deviations))
    if not np.isfinite(mad):
        return mask

    if mad < 1e-12:


        keep_finite = np.isclose(
            finite_values,
            median,
            rtol=1e-9,
            atol=1e-12,
        )
    else:
        modified_z = 0.6745 * deviations / mad
        keep_finite = modified_z <= z_thresh

    mask[np.flatnonzero(mask)] = keep_finite
    return mask


def _leave_one_out_slope_ratio_inlier_mask(
    x: np.ndarray,
    y_unique: np.ndarray,
    y_signal: np.ndarray,
    z_thresh: float = 3.5,
    min_n: int = 6,
    min_relative_change: float = CPG_SLOPE_RATIO_MIN_CHANGE,
) -> np.ndarray:

    n = len(x)
    if n < min_n:
        return np.ones(n, dtype=bool)
    if min_relative_change < 0:
        raise ValueError('Minimum slope-ratio change must be non-negative')

    xx = x * x
    xy_unique = x * y_unique
    xy_signal = x * y_signal
    total_xx = np.sum(xx)
    total_unique = np.sum(xy_unique)
    total_signal = np.sum(xy_signal)
    remaining_xx = total_xx - xx
    remaining_unique = total_unique - xy_unique
    remaining_signal = total_signal - xy_signal

    with np.errstate(divide='ignore', invalid='ignore'):
        full_slope_unique = total_unique / total_xx
        full_slope_signal = total_signal / total_xx
        full_slope_ratio = full_slope_unique / full_slope_signal
        slope_unique = remaining_unique / remaining_xx
        slope_signal = remaining_signal / remaining_xx
        leave_one_out_ratio = slope_unique / slope_signal
        log_slope_ratio = np.log(leave_one_out_ratio)
        relative_change = np.abs(leave_one_out_ratio / full_slope_ratio - 1.0)


    robust_inlier = _robust_inlier_mask(
        log_slope_ratio,
        z_thresh=z_thresh,
        min_n=min_n,
    )
    material_change = np.isfinite(relative_change) & (
        relative_change >= min_relative_change
    )
    return robust_inlier | ~material_change


def _fit_cpg_slope_pair(
    x: np.ndarray,
    y_unique: np.ndarray,
    y_signal: np.ndarray,
    z_thresh: float = 3.5,
    min_n: int = 6,
    min_ratio_change: float = CPG_SLOPE_RATIO_MIN_CHANGE,
) -> tuple[float, float, np.ndarray]:

    x = np.asarray(x, dtype=float)
    y_unique = np.asarray(y_unique, dtype=float)
    y_signal = np.asarray(y_signal, dtype=float)
    if not (x.shape == y_unique.shape == y_signal.shape):
        raise ValueError('CpG slope arrays must have identical shapes')
    if min_ratio_change < 0:
        raise ValueError('Minimum slope-ratio change must be non-negative')

    valid = (
        np.isfinite(x)
        & np.isfinite(y_unique)
        & np.isfinite(y_signal)
        & (x > 0)
        & (y_unique > 0)
        & (y_signal > 0)
    )
    inlier_mask = np.zeros(len(x), dtype=bool)
    if valid.sum() < 2:
        return np.nan, np.nan, inlier_mask

    valid_positions = np.flatnonzero(valid)
    x_valid = x[valid]
    unique_valid = y_unique[valid]
    signal_valid = y_signal[valid]
    keep = np.ones(len(x_valid), dtype=bool)

    if len(x_valid) >= min_n:
        with np.errstate(divide='ignore', invalid='ignore'):
            log_point_ratio = np.log(unique_valid / signal_valid)

        ratio_outlier = ~_robust_inlier_mask(
            log_point_ratio, z_thresh, min_n
        )
        ratio_influential = ~_leave_one_out_slope_ratio_inlier_mask(
            x_valid,
            unique_valid,
            signal_valid,
            z_thresh=z_thresh,
            min_n=min_n,
            min_relative_change=min_ratio_change,
        )


        confirmed_outlier = ratio_outlier & ratio_influential
        keep &= ~confirmed_outlier


        if keep.sum() < 2:
            keep = np.ones(len(x_valid), dtype=bool)

    inlier_mask[valid_positions] = keep
    x_fit = x[inlier_mask]
    unique_fit = y_unique[inlier_mask]
    signal_fit = y_signal[inlier_mask]
    denominator = float(np.dot(x_fit, x_fit))
    if not np.isfinite(denominator) or denominator <= 0:
        return np.nan, np.nan, inlier_mask

    slope_unique = float(np.dot(x_fit, unique_fit) / denominator)
    slope_signal = float(np.dot(x_fit, signal_fit) / denominator)
    return slope_unique, slope_signal, inlier_mask


def plot_cpg_density(
    df: pd.DataFrame,
    clone_color_map: dict = None,
    out_path: Path = None,
    min_slope_ratio_change: float = CPG_SLOPE_RATIO_MIN_CHANGE,
):
    """展示逐样本 input pairs–unique CpG 散点及密度汇总，分别保存 PDF。

    散点保留有限的正 input pairs 与非负 unique sites；汇总仅在异常点
    实质影响 slope ratio 时共同剔除，不改变散点中的原始计数。
    每个 CloneID 的过原点拟合线与斜率标注复用 input-pair 汇总的同一次拟合。
    """
    if min_slope_ratio_change < 0:
        raise ValueError('Minimum slope-ratio change must be non-negative')

    metric_a = 'Sample_Unique_CpG_Sites'
    metric_b = 'Sample_CpG_Signals'

    denominator_groups = [
        ('Input_Read_Pairs', 'Input Pair', 'input pair', 'Input-Pair'),
        ('Final_Retained_Pairs', 'Final Retained Pair', 'final retained pair', 'Final-Pair'),
    ]

    if df.empty or 'CloneID' not in df.columns:
        return
    missing_denominators = [
        x_col for x_col, *_ in denominator_groups if x_col not in df.columns
    ]
    if missing_denominators:
        print(f"Missing required denominator for CpG slopes: {', '.join(missing_denominators)}")
        return
    _load_plot_dependencies()

    plot_df = df.copy()
    plot_df['CloneID'] = plot_df['CloneID'].map(normalize_clone_id_value)
    plot_df = plot_df.dropna(subset=['CloneID'])
    clones = sorted(plot_df['CloneID'].dropna().unique().tolist())
    if not clones:
        return
    clone_display_map, omitted = make_cloneid_diff_labels(clones)
    if clone_color_map is None:
        clone_color_map = generate_clone_color_map(clones)
    clone_color_map = {
        normalize_clone_id_value(key): value
        for key, value in clone_color_map.items()
        if isinstance(normalize_clone_id_value(key), str)
    }

    def calc_slopes_by_clone(x_col: str) -> tuple[pd.Series, pd.Series]:
        required = [x_col, metric_a, metric_b]
        if any(column not in plot_df.columns for column in required):
            return pd.Series(dtype=float), pd.Series(dtype=float)

        unique_rows = []
        signal_rows = []
        for clone in clones:
            subset_columns = required + (
                ['Sample_ID'] if 'Sample_ID' in plot_df.columns else []
            )
            subset = plot_df.loc[
                plot_df['CloneID'] == clone,
                subset_columns,
            ].copy()
            x = pd.to_numeric(subset[x_col], errors='coerce').to_numpy(dtype=float)
            y_unique = pd.to_numeric(
                subset[metric_a], errors='coerce'
            ).to_numpy(dtype=float)
            y_signal = pd.to_numeric(
                subset[metric_b], errors='coerce'
            ).to_numpy(dtype=float)
            slope_unique, slope_signal, inlier_mask = _fit_cpg_slope_pair(
                x,
                y_unique,
                y_signal,
                min_ratio_change=min_slope_ratio_change,
            )
            if not (np.isfinite(slope_unique) and np.isfinite(slope_signal)):
                continue

            unique_rows.append((clone, slope_unique))
            signal_rows.append((clone, slope_signal))

            valid = (
                np.isfinite(x)
                & np.isfinite(y_unique)
                & np.isfinite(y_signal)
                & (x > 0)
                & (y_unique > 0)
                & (y_signal > 0)
            )
            excluded = valid & ~inlier_mask
            if excluded.any():
                excluded_rows = subset.iloc[np.flatnonzero(excluded)]
                if 'Sample_ID' in excluded_rows.columns:
                    labels = excluded_rows['Sample_ID'].astype(str).tolist()
                else:
                    labels = [f'row {index}' for index in excluded_rows.index]
                print(
                    'Excluded CpG slope-ratio influential outlier(s) '
                    f'[{x_col}, CloneID={clone}, '
                    f'min_change={min_slope_ratio_change:.1%}]: '
                    f'{", ".join(labels)}'
                )

        return pd.Series(dict(unique_rows)), pd.Series(dict(signal_rows))

    def _format_slope_value(v: float) -> str:
        if not np.isfinite(v):
            return ''
        absolute = abs(v)
        if absolute != 0 and (absolute >= 10000.0 or absolute < 0.01):
            return f'{v:.2e}'
        return f'{v:.3f}'.rstrip('0').rstrip('.')

    def _annotate_bars_h(ax, bars, values: np.ndarray) -> None:
        if len(values) == 0:
            return
        xmax = float(np.nanmax(values))
        if np.isfinite(xmax) and xmax > 0:
            ax.set_xlim(0, xmax * 1.24)
            offset = xmax * 0.02
        else:
            offset = 0.0
        for bar, value in zip(bars, values, strict=True):
            label = _format_slope_value(float(value))
            if not label:
                continue
            x_position = bar.get_width() + offset
            y_position = bar.get_y() + bar.get_height() / 2.0
            ax.text(
                x_position,
                y_position,
                label,
                ha='left',
                va='center',
                fontsize=9,
                clip_on=False,
            )

    def _plot_slope_panel(
        ax,
        slope: pd.Series,
        title: str,
        xlabel: str,
        letter: str,
    ) -> None:
        values = slope.reindex(clones).dropna()
        if values.empty:
            ax.text(0.5, 0.5, 'Data Missing', transform=ax.transAxes, ha='center', va='center')
            ax.axis('off')
            return
        y = np.arange(len(values))
        bars = ax.barh(
            y,
            values.values,
            color=[clone_color_map.get(clone, UNKNOWN_COLOR) for clone in values.index],
            edgecolor='white',
            linewidth=0.6,
            alpha=0.9,
        )
        ax.set_yticks(y)
        ax.set_yticklabels(
            [clone_display_map.get(clone, clone) for clone in values.index],
            fontsize=8,
        )
        format_ax(ax, title=title, xlabel=xlabel, ylabel='CloneID')
        add_panel_letter(ax, letter)
        _annotate_bars_h(ax, bars, values.values)

    def _plot_ratio_panel(
        ax,
        slope_a: pd.Series,
        slope_b: pd.Series,
        title: str,
        letter: str,
    ) -> None:
        comp_df = pd.DataFrame(
            {
                'Unique': slope_a.reindex(clones),
                'Theoretical': slope_b.reindex(clones),
            }
        ).dropna(subset=['Unique', 'Theoretical'], how='any')
        if not comp_df.empty:
            comp_df = comp_df[
                np.isfinite(comp_df['Unique'])
                & np.isfinite(comp_df['Theoretical'])
                & (comp_df['Theoretical'] != 0)
            ]
        if comp_df.empty:
            ax.text(0.5, 0.5, 'Data Missing', transform=ax.transAxes, ha='center', va='center')
            ax.axis('off')
            return
        ratio = (comp_df['Unique'] / comp_df['Theoretical']).astype(float)
        clone_order = ratio.index.tolist()
        y = np.arange(len(ratio))
        bars = ax.barh(
            y,
            ratio.to_numpy(dtype=float),
            color=[clone_color_map.get(clone, UNKNOWN_COLOR) for clone in clone_order],
            edgecolor='white',
            linewidth=0.6,
            alpha=0.9,
        )
        ax.set_yticks(y)
        ax.set_yticklabels(
            [clone_display_map.get(clone, clone) for clone in clone_order],
            fontsize=8,
        )
        ax.axvline(
            1.0,
            color=Colors.REFLINE,
            linestyle='--',
            linewidth=1.2,
            alpha=0.9,
            zorder=0,
        )
        format_ax(ax, title=title, xlabel='Ratio', ylabel='CloneID')
        add_panel_letter(ax, letter)
        _annotate_bars_h(ax, bars, ratio.to_numpy(dtype=float))

    legend_levels = [clone_display_map.get(clone, clone) for clone in clones]
    legend_color_map = {
        clone_display_map.get(clone, clone): clone_color_map.get(clone, UNKNOWN_COLOR)
        for clone in clones
    }
    input_slopes = calc_slopes_by_clone('Input_Read_Pairs')
    scatter_legend_levels = []
    scatter_legend_colors = {}
    fig_scatter, ax_scatter = plt.subplots(figsize=(10, 7), dpi=PlotParams.DPI)
    x = pd.to_numeric(plot_df['Input_Read_Pairs'], errors='coerce')
    y = pd.to_numeric(plot_df[metric_a], errors='coerce')
    valid = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y >= 0)
    for clone in clones:
        selected = valid & (plot_df['CloneID'] == clone)
        ax_scatter.scatter(
            x[selected], y[selected],
            color=clone_color_map.get(clone, UNKNOWN_COLOR),
            s=24, alpha=0.75, edgecolors='none',
        )
        slope = input_slopes[0].get(clone, np.nan)
        if np.isfinite(slope) and selected.any():
            line_x = np.array([0.0, float(x[selected].max())])
            ax_scatter.plot(
                line_x, slope * line_x,
                color=clone_color_map.get(clone, UNKNOWN_COLOR),
                linewidth=1.5, alpha=0.9,
            )
        slope_label = _format_slope_value(slope) if np.isfinite(slope) else 'N/A'
        label = f'{clone_display_map.get(clone, clone)} (slope={slope_label})'
        scatter_legend_levels.append(label)
        scatter_legend_colors[label] = clone_color_map.get(clone, UNKNOWN_COLOR)
    format_ax(
        ax_scatter, title='Unique CpG Sites vs Input Pairs',
        xlabel='Input read pairs', ylabel='Unique CpG sites',
    )
    ax_scatter.set_xlim(left=0)
    ax_scatter.set_ylim(bottom=0)
    add_right_legend(
        fig_scatter, scatter_legend_levels, scatter_legend_colors, title='CloneID'
    )
    add_cloneid_omission_note(fig_scatter, omitted, y=0.015)
    fig_scatter.tight_layout(rect=[0, 0.03, 0.80, 0.985])
    if out_path:
        scatter_path = Path(out_path)
        fig_scatter.savefig(
            scatter_path.with_name(
                f'{scatter_path.stem}_Unique_vs_Input_Pairs{scatter_path.suffix}'
            ),
            bbox_inches='tight',
        )
    plt.show()
    plt.close(fig_scatter)

    fig_h = max(28.0, 0.56 * len(clones) + 12.0)
    fig = plt.figure(figsize=(22, fig_h), dpi=PlotParams.DPI)
    gs = fig.add_gridspec(4, 2, height_ratios=[1.0, 1.2, 1.0, 1.2])
    panel_letters = iter('ABCDEF')
    for group_idx, (
        x_col,
        title_denominator,
        xlabel_denominator,
        ratio_prefix,
    ) in enumerate(denominator_groups):
        slope_a, slope_b = input_slopes if group_idx == 0 else calc_slopes_by_clone(x_col)
        row = group_idx * 2
        ax_a = fig.add_subplot(gs[row, 0])
        ax_b = fig.add_subplot(gs[row, 1])
        ax_ratio = fig.add_subplot(gs[row + 1, :])
        _plot_slope_panel(
            ax_a,
            slope_a,
            f'Unique CpG Sites per {title_denominator}',
            f'Unique CpG sites / {xlabel_denominator}',
            next(panel_letters),
        )
        _plot_slope_panel(
            ax_b,
            slope_b,
            f'Host CpG Signals per {title_denominator}',
            f'Host CpG coverage signals / {xlabel_denominator}',
            next(panel_letters),
        )
        _plot_ratio_panel(
            ax_ratio,
            slope_a,
            slope_b,
            f'{ratio_prefix} Normalization: Unique / Coverage-Signal Slope Ratio',
            next(panel_letters),
        )

    add_right_legend(fig, legend_levels, legend_color_map, title='CloneID')
    add_cloneid_omission_note(fig, omitted, y=0.015)
    plt.tight_layout(rect=[0, 0, 0.885, 0.985], w_pad=2.6, h_pad=2.8)
    if out_path:
        plt.savefig(out_path, bbox_inches='tight')
    plt.show()
    plt.close(fig)


def plot_pca_batch_effect(
    df: pd.DataFrame,
    clone_color_map: dict = None,
    out_path: Path = None,
):
    """仅用指标完整的细胞绘制 QC PCA；保留观测零，注明缺失/非法值排除数。"""
    def skip(message: str) -> None:
        if out_path is not None:
            Path(out_path).unlink(missing_ok=True)
        print(message)

    if df.empty:
        skip('No cells: PCA skipped.')
        return
    features = [
        'Native_Mapping%',
        'Sample_Unique_CpG_Sites',
        'Sample_Meth_CpG_Rate%',
        'Final_Pair_Yield%',
    ]
    valid_f = [f for f in features if f in df.columns]
    if len(valid_f) >= 3 and len(df) >= 2:
        _load_plot_dependencies()
        numeric = df[valid_f].apply(pd.to_numeric, errors='coerce')
        complete = np.isfinite(numeric.to_numpy(dtype=float, na_value=np.nan)).all(axis=1)
        excluded = int((~complete).sum())
        plot_df = df.loc[complete].copy()
        if len(plot_df) < 2:
            skip(f'PCA skipped: fewer than 2 complete cells ({len(plot_df)}/{len(df)}; excluded {excluded}).')
            return
        plot_df['CloneID'] = plot_df['CloneID'].map(normalize_clone_id_value)
        # 未定义的宿主甲基化率不等于真实0%；不能通过填零改变投影输入。
        values = StandardScaler().fit_transform(numeric.loc[complete])
        if np.all(np.var(values, axis=0) <= np.finfo(float).eps):
            skip('No variation: PCA skipped.')
            return
        pca = PCA(n_components=2)
        coords = pca.fit_transform(values)
        if not np.all(np.isfinite(pca.explained_variance_ratio_)):
            skip('No variation: PCA skipped.')
            return
        group_levels = sorted(
            plot_df['CloneID'].dropna().astype(str).unique().tolist()
        )
        label_map, omitted = make_cloneid_diff_labels(group_levels)
        plot_df['CloneID_Display'] = plot_df['CloneID'].map(label_map)
        if clone_color_map is None:
            clone_color_map = generate_clone_color_map(group_levels)
        display_levels = [label_map.get(group, group) for group in group_levels]
        display_color_map = {
            label_map.get(group, group): clone_color_map.get(group, UNKNOWN_COLOR)
            for group in group_levels
        }
        fig, ax = plt.subplots(figsize=(8, 6), dpi=PlotParams.DPI)
        fig.text(
            0.5, 0.98,
            f'PCA: {len(plot_df)} / {len(df)} cells; excluded {excluded} with missing or invalid QC',
            ha='center', va='top', fontsize=FontSizes.LEGEND,
        )
        sns.scatterplot(
            x=coords[:, 0],
            y=coords[:, 1],
            hue=plot_df.get('CloneID_Display'),
            palette=display_color_map,
            ax=ax,
            s=PlotParams.SCATTER_SIZE,
            alpha=PlotParams.SCATTER_ALPHA,
            edgecolor='white',
            linewidth=0.55,
            legend=False,
        )
        format_ax(
            ax,
            title='PCA: Batch Effect Check',
            xlabel=f'PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)',
            ylabel=f'PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)',
        )
        add_right_legend(fig, display_levels, display_color_map, title='CloneID')
        add_cloneid_omission_note(fig, omitted, y=0.015)
        plt.tight_layout(rect=[0, 0.03, 0.885, 0.958], w_pad=2.2, h_pad=2.4)
        if out_path:
            plt.savefig(out_path, bbox_inches='tight')
        plt.show()
        plt.close(fig)
    else:
        skip('⚠️ 核心特征不足，无法进行 PCA。')


def plot_final_qc_dashboard_batch(
    df: pd.DataFrame,
    metrics: Any,
    group_col: str = 'CloneID',
    clone_id: Any = None,
    clone_color_map: dict = None,
    out_dir: Path = None,
    filename_prefix: str = 'Final_QC_Dashboard',
    mode: str = 'per_clone_samples',
    figsize: tuple = (16.5, 12),
):
    """为指标列表生成一个或多个最终 QC 面板。"""
    metric_list = normalize_metric_list(metrics)
    if not metric_list:
        raise ValueError('No QC metric selected; specify a metric or a nonempty list.')

    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

    clone_suffix = make_clone_selection_suffix(clone_id)

    for metric_name in metric_list:
        safe_metric = (
            str(metric_name)
            .replace('%', 'pct')
            .replace('/', '_')
            .replace('\\', '_')
            .replace(' ', '_')
            .replace(':', '_')
        )
        out_path = None
        if out_dir is not None:
            out_path = out_dir / f'{filename_prefix}_{mode}_{safe_metric}_{clone_suffix}.pdf'

        plot_final_qc_dashboard(
            df=df,
            group_col=group_col,
            clone_id=clone_id,
            mode=mode,
            metric=metric_name,
            clone_color_map=clone_color_map,
            out_path=out_path,
            figsize=figsize,
        )
