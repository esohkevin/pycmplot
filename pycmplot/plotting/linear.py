"""
pycmplot.plotting.linear
========================
Multi-track linear Manhattan plot.
"""

from __future__ import annotations

import logging
from typing import Optional

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from matplotlib.patches import FancyArrowPatch
from natsort import natsort_keygen

from pycmplot.constants import CHROM_ORDER
from pycmplot.stats import get_highlight_snps
from pycmplot.io import get_output_paths

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Annotation helpers (cluster-aware label spreading)
# ---------------------------------------------------------------------------

def _cluster_annotations_by_chr(
    annot_df,
    chr_col: str = "CHR",
    x_col: str = "x",
    window_size: float = 50e6,
) -> list[list]:
    """Cluster annotations within each chromosome by genomic proximity."""
    clusters: list[list] = []
    for _chr_name, df_chr in annot_df.groupby(chr_col):
        df_chr = df_chr.sort_values(x_col)
        current_cluster = [df_chr.index[0]]
        last_x = df_chr.iloc[0][x_col]

        for idx, row in df_chr.iloc[1:].iterrows():
            x = row[x_col]
            if x - last_x <= window_size:
                current_cluster.append(idx)
            else:
                clusters.append(current_cluster)
                current_cluster = [idx]
            last_x = x

        clusters.append(current_cluster)
    return clusters


def _draw_annotation_arrows(
    ax,
    annot_df,
    chr_col: str,
    label_col: str,
    offsets: dict,
    chr_max: dict,
    spread_width: float = 60e6,
    y_tip: float = 0.0,
    y_text: float = 0.55,
) -> None:
    """Draw angled FancyArrowPatch arrows from text labels to signal positions."""
    annot_df = annot_df.sort_values(by=[chr_col, "x"], key=natsort_keygen())
    last_xtext = 0 - spread_width

    for chr_name, df_chr in annot_df.groupby(chr_col, sort=False):
        df_chr = df_chr.sort_values("x")
        chr_start = offsets[chr_name]
        chr_end = offsets[chr_name] + chr_max[chr_name]

        x_signals = df_chr["x"].values
        labels = df_chr[label_col].values
        n = len(df_chr)

        # Adaptive spread
        chr_range = chr_end - chr_start
        sw = spread_width
        pad = sw / int(str(sw)[:2]) / 2
        while sw > chr_range:
            sw -= pad

        sig_start = df_chr["x"].iloc[0]
        xmin = sig_start - sw
        xmax = xmin + n * sw
        x_texts = np.arange(xmin, xmax, sw)

        first_xtext = x_texts[0]
        while first_xtext <= last_xtext:
            x_texts = [xv + sw for xv in x_texts]
            first_xtext = x_texts[0]

        for x_sig, x_txt, label in zip(x_signals, x_texts, labels):
            dx = x_txt - x_sig
            rad = 0.15 * np.sign(dx)

            arrow = FancyArrowPatch(
                (x_txt, y_text),
                (x_sig, y_tip - 0.05),
                arrowstyle="-|>",
                mutation_scale=12,
                lw=0.6,
                color="grey",
                alpha=0.5,
                connectionstyle=f"arc3,rad={rad}",
            )
            ax.add_patch(arrow)

            ax.text(
                x_txt,
                y_text + 0.02,
                str(label),
                rotation=45,
                ha="left",
                va="bottom",
                fontsize=10,
                clip_on=False,
                color="black",
                fontstyle="italic",
                fontweight="regular",
            )

        last_xtext = x_texts[-1]


# ---------------------------------------------------------------------------
# Public function
# ---------------------------------------------------------------------------

def plot_linearm(
    tracks: list,
    track_labels: Optional[list[str]] = None,
    annot_df=None,
    highlight: bool = False,
    highlight_thresh: float = 5e-8,
    trim_pval: Optional[float] = None,
    logp: bool = True,
    label_col: Optional[str] = "label",
    chr_order: Optional[list[str]] = None,
    chr_spacing: float = 9e6,
    track_heights: Optional[list[float]] = None,
    track_spacing: float = 0.10,
    point_size: float = 5,
    colors: Optional[list[str]] = ['steelblue','orange'],
    sig_lines: Optional[list[dict]] = None,
    plt_name: Optional[str] = None,
    no_track_labels: bool = False,
    fig_format: Optional[str] = None,
    dpi: int = 300,
    figsize: tuple = (15, 9),
):

    chr_col = "CHR"
    pos_col = "POS"
    p_col = "P"

    if chr_order is None:
        chr_order = CHROM_ORDER

    chr_to_idx = {c: i for i, c in enumerate(chr_order)}

    # ------------------------------------------------------------------
    # Prep DataFrames
    # ------------------------------------------------------------------
    def _prep(df):
        df = df.copy()
        """ ALREADY HANDLED IN DATA LOADER FUNCTION
        if trim_pval:
            df = df[df[p_col] < trim_pval]
        if logp:
            df["logP"] = -np.log10(df[p_col])
        """

        df[chr_col] = (
            df[chr_col]
            .astype(str)
            .str.replace("chr", "", regex=False)
            .str.upper()
            .replace({"23": "X", "24": "Y", "M": "MT", "MTDNA": "MT"})
        )

        if highlight:
            df, _ = get_highlight_snps(
                df=df,
                window=500_000,
                highlight_thresh=highlight_thresh,
                logp=logp,
            )

        df = df[df[chr_col].isin(chr_order)]
        df["chr_idx"] = df[chr_col].map(chr_to_idx)
        return df.sort_values(["chr_idx", pos_col])

    tracks = [_prep(df) for df in tracks]
    if annot_df is not None:
        annot_df = _prep(annot_df)

    # ------------------------------------------------------------------
    # Cumulative x-axis positions
    # ------------------------------------------------------------------
    chr_max: dict[str, float] = {}
    offsets: dict[str, float] = {}
    offset = 0.0

    for c in chr_order:
        max_pos = max(
            [df[df[chr_col] == c][pos_col].max() for df in tracks if c in df[chr_col].values]
            + [0]
        )
        chr_max[c] = max_pos
        offsets[c] = offset
        offset += max_pos + chr_spacing

    def _add_cum(df):
        df = df.copy()
        df["x"] = df.apply(lambda r: r[pos_col] + offsets[r[chr_col]], axis=1)
        return df

    tracks = [_add_cum(df) for df in tracks]
    if annot_df is not None:
        annot_df = _add_cum(annot_df)

    # ------------------------------------------------------------------
    # Figure layout
    # ------------------------------------------------------------------
    n_tracks = len(tracks)

    if track_heights is None:
        track_heights = [1] + [3] * n_tracks

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(
        n_tracks + 1, 1,
        height_ratios=track_heights,
        hspace=track_spacing,
    )

    ax_annot = fig.add_subplot(gs[0, 0])
    axes = [ax_annot]
    for i in range(n_tracks):
        axes.append(fig.add_subplot(gs[i + 1, 0], sharex=ax_annot))

    if colors is None:
        colors = ["gray", "steelblue"]

    # Per-track highlight colours from tab20 colormap
    cmap = plt.get_cmap("tab20")
    hex_colors = [mcolors.to_hex(cmap(i / n_tracks)) for i in range(n_tracks)]

    # ------------------------------------------------------------------
    # Plot data tracks
    # ------------------------------------------------------------------
    t_labels = track_labels or [f"Track {i+1}" for i in range(n_tracks)]

    for i, (ax, df, t_label, h_color) in enumerate(
        zip(axes[1:], tracks, t_labels, hex_colors)
    ):
        color_cycle = [colors[j % len(colors)] for j in df["chr_idx"]]
        df = df[df[p_col] >= 0]

        y_vals = df["logP"] if logp else df[p_col]
        ax.scatter(df["x"], y_vals, c=color_cycle, s=point_size)

        if highlight:
            sig = df[df["in_locus"]]
            if not sig.empty:
                sig_y = sig["logP"] if logp else sig[p_col]
                ax.scatter(sig["x"].to_numpy(), sig_y.to_numpy(), s=point_size,
                           marker="o", color="brown")

        if no_track_labels:
            pass
        else:
            ax.set_ylabel(t_label, color="black")

        if sig_lines is not None and i < len(sig_lines):
            sl = sig_lines[i]
            if "genome" in sl:
                ax.axhline(y=sl["genome"], color="red", linestyle="--", linewidth=0.6)
            if "suggestive" in sl:
                ax.axhline(y=sl["suggestive"], color="grey", linestyle="--", linewidth=0.5)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        left_pad = chr_spacing * 0.2
        xmax = max(offsets[c] + chr_max[c] for c in chr_order)
        ax.set_xlim(-left_pad, xmax)

    # ------------------------------------------------------------------
    # Annotation track
    # ------------------------------------------------------------------
    if annot_df is not None:
        # Vertical lines across all data tracks
        for x in annot_df["x"].values:
            for ax in axes[1:]:
                ax.axvline(x, color="grey", alpha=0.45, linewidth=0.7,
                           linestyle="--", zorder=0)

        _draw_annotation_arrows(
            ax_annot,
            annot_df,
            chr_col=chr_col,
            label_col=label_col,
            offsets=offsets,
            chr_max=chr_max,
            spread_width=60e6,
        )

    ax_annot.set_ylim(0, 1)
    ax_annot.axis("off")

    # ------------------------------------------------------------------
    # Chromosome labels on x-axis
    # ------------------------------------------------------------------
    xticks, xlabels = [], []
    for c in chr_order:
        if chr_max[c] == 0:
            continue
        start = offsets[c]
        end = offsets[c] + chr_max[c]
        mid = (start + end) / 2
        xticks.append(mid)
        xlabels.append(c)
        for ax in axes:
            ax.axvline(end, color="lightgray", linewidth=0.1, alpha=0.05)

    axes[-1].set_xticks(xticks)
    axes[-1].set_xticklabels(xlabels)
    axes[-1].set_xlabel("Chromosome", fontsize=12)

    for ax in axes[:-1]:
        ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
        ax.spines["bottom"].set_visible(False)

    plt.subplots_adjust(hspace=track_spacing, left=0.08)
    plt.tight_layout()

    fig.text(
        0, 0.5,
        "-log\u2081\u2080(p-value)" if logp else p_col,
        va="center",
        rotation="vertical",
        fontsize=12,
    )

    if plt_name:
        fmt = fig_format or Path(plt_name).suffix.lstrip(".") or "png"
        plt.savefig(plt_name.lower(), format=fmt, dpi=dpi)
        logger.info("Saved linear Manhattan plot: %s", plt_name.lower())

    return fig, axes


def plot_linear(
    sumstats_loaded: list[str],
    trim_pval: Optional[float] = None,
    track_heights: list[float] = None,
    logp: bool = False,
    point_size: Optional[float] = None,
    highlight: bool = False,
    highlight_thresh: float = 5e-8,
    hits_table: Optional[pd.DataFrame] = None,
    label_col: Optional[str] = None,
    chr_spacing: Optional[float] = None,
    track_spacing: Optional[float] = None,
    colors: list[str] = None,
    signif_lines: Optional[dict] = None,
    plot_title: Optional[str] = None,
    no_track_labels: bool = False,
    dpi: Optional[int] = None,
    output_format: Optional[str] = None,
    output_dir: Optional[str] = '.',
    figsize: Optional[tuple] = None,
):
    """Generate a multi-track linear Manhattan plot.

    Parameters
    ----------
    sumstats_loaded:
        List of DataFrames, one per GWAS trait.  Each must have columns
    highlight:
        Whether to highlight significant signals.
    hits_table:
        Optional DataFrame of lead SNPs to annotate (must contain *chr_col*,
        *pos_col*, *label_col*).
    label_col:
        Column to use in the annot_df e.g. column containing gene names.
    highlight:
        Highlight loci within ``500 kb`` of a lead SNP.
    chr_spacing:
        Gap (bp) inserted between chromosomes on the x-axis.
    signif_lines:
        List of ``{"genome": float, "suggestive": float}`` dicts, one per track.
    plot_title:
        Output file path (extension determines format when *fig_format* is ``None``).
    output_format:
        Override output format (e.g. ``'png'``, ``'pdf'``).
    figsize:
        Figure size e.g. (15, 8)
    dpi:
        FIgure resolution (default: 300)


    Returns
    -------
    (fig, axes)
    """

    dfs      = [v[0] for v in sumstats_loaded.values()]
    t_labels = list(sumstats_loaded.keys())

    if not track_heights:
        t_heights = None
    else:
        t_heights = [float(x) for x in track_heights]

    # plot name
    (
        plt_name, 
        table_out
    ) = get_output_paths(
        labels = t_labels,
        mode='lm', 
        logp=logp, 
        output_dir=output_dir, 
        plot_title=plot_title, 
        output_format=output_format
    )

    fig, axes = plot_linearm(
        tracks=dfs,
        track_labels=t_labels,
        trim_pval=trim_pval,
        logp=True if logp else False,
        point_size=point_size,
        highlight=highlight,
        highlight_thresh=highlight_thresh,
        annot_df=hits_table if not hits_table.empty else None,
        label_col=label_col,
        chr_spacing=chr_spacing,
        track_heights=t_heights,
        track_spacing=track_spacing,
        colors=colors,
        sig_lines=signif_lines,
        plt_name=plt_name,
        no_track_labels = no_track_labels,
        dpi=dpi,
        fig_format=output_format,
        figsize=figsize,
    )

    return fig, axes