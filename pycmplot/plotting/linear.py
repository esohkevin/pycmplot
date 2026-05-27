"""
pycmplot.plotting.linear
========================

Multi-track stacked linear Manhattan plot.

The module exposes two public functions:

* :func:`plot_linear` — the user-facing entry point.  Accepts the
  ``sumstats_loaded`` dict produced by
  :func:`~pycmplot.io.get_sumstats_and_merged_sector_list`, resolves
  output paths, and delegates rendering to :func:`plot_linearm`.
* :func:`plot_linearm` — the core rendering engine.  Accepts a list of
  DataFrames and a fully resolved set of plotting parameters, builds the
  matplotlib figure, draws all tracks, and saves the file.

Internal helpers:

* :func:`_cluster_annotations_by_chr` — groups annotation points that
  are close together on the same chromosome so label-spreading can be
  applied per cluster rather than globally.
* :func:`_draw_annotation_arrows` — places angled
  :class:`~matplotlib.patches.FancyArrowPatch` arrows from spread gene
  labels down to their corresponding signal positions.
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
from pycmplot.annotation import get_annotation_column

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
    """Cluster annotation points within each chromosome by genomic proximity.

    Groups rows of *annot_df* on the same chromosome into clusters such that
    consecutive points within *window_size* base-pairs are placed in the same
    cluster.  Clusters are used by :func:`_draw_annotation_arrows` to
    determine independent label-spreading regions.

    Parameters
    ----------
    annot_df : pandas.DataFrame
        Annotation DataFrame containing at least *chr_col* and *x_col*
        (cumulative x-axis position).
    chr_col : str, optional
        Name of the chromosome column.  Default ``'CHR'``.
    x_col : str, optional
        Name of the cumulative x-axis position column.  Default ``'x'``.
    window_size : float, optional
        Maximum gap in cumulative x-axis units between two points that should
        be considered part of the same cluster.  Default ``50e6`` (50 Mb).

    Returns
    -------
    list of list
        Each inner list contains the integer index values of *annot_df* rows
        that form one cluster.
    """

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
    y_text: float = 0.35,
) -> None:
    """Draw angled arrow annotations from gene-label text to signal positions.

    For each significant locus in *annot_df*, places the gene/SNP label text
    above the track and connects it to the corresponding scatter-plot point
    with a curved :class:`~matplotlib.patches.FancyArrowPatch` arrow.  Labels
    within the same chromosome are spread horizontally to avoid overlap; the
    arrow curvature direction (clockwise/counter-clockwise arc) is determined
    automatically from the sign of the horizontal displacement.

    The function operates on *ax* (the annotation sub-panel at the top of the
    figure) and returns ``None`` — it modifies the axes in place.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The annotation axes (top sub-panel, invisible background).
    annot_df : pandas.DataFrame
        Lead-SNP annotation DataFrame with columns *chr_col*, ``'x'``
        (cumulative x-axis position), and *label_col*.
    chr_col : str
        Name of the chromosome column.
    label_col : str
        Name of the column containing the text labels (gene symbols or rsIDs).
    offsets : dict
        Mapping of ``chromosome → cumulative_x_start`` as computed from the
        main track data.
    chr_max : dict
        Mapping of ``chromosome → max_position`` for each chromosome.
    spread_width : float, optional
        Horizontal spacing between adjacent labels within a cluster.
        Default ``60e6`` (60 Mb).
    y_tip : float, optional
        y-coordinate (in axes data units, 0–1) of the arrowhead tip.
        Default ``0.0``.
    y_text : float, optional
        y-coordinate of the label text anchor.  Default ``0.55``.
    """

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
                #connectionstyle="arc,angleA=5,angleB=90,armA=0,armB=100,rad=5",

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


def _layout_annotation_labels(
    annot_df,
    chr_col: str,
    label_col: str,
    offsets: dict,
    chr_max: dict,
    fsize: float = 8,
    char_width_factor: float = 2.5e5,
):
    """Compute label x positions and rail assignments without drawing.

    Runs the same two-phase placement algorithm as
    :func:`_draw_annotation_arrows_3` (horizontal width-aware spread,
    then greedy multi-rail stacking) but returns the layout instead of
    rendering it.  Used by :func:`plot_linearm` to size the annotation
    panel's gridspec row *before* the figure is created — without it,
    the panel always gets a fixed ratio and stacked rails end up
    physically squashed regardless of how generous ``y_stack_step`` is
    in axes-relative units.

    Returns
    -------
    x_texts : numpy.ndarray
        Final anchor x positions for each label.
    rail_ids : numpy.ndarray of int
        Rail index per label (0 = bottom).
    max_rail : int
        Maximum rail index used (rail count = ``max_rail + 1``).
    rail_start, rail_end : float
        Effective horizontal bounds used for placement.
    """
    annot_df = annot_df.sort_values(by=[chr_col, "x"], key=natsort_keygen())
    x_signals = annot_df["x"].values
    labels    = annot_df[label_col].astype(str).values
    n = len(x_signals)
    if n == 0:
        return np.array([]), np.array([], dtype=int), 0, 0.0, 0.0

    genome_start = min(offsets.values())
    genome_end = max(offsets[c] + chr_max[c] for c in chr_max)
    genome_width = float(genome_end - genome_start) or 1.0
    left_margin  = 0.06 * genome_width
    right_margin = 0.02 * genome_width
    rail_start = genome_start + left_margin
    rail_end   = genome_end   - right_margin
    rail_width = max(rail_end - rail_start, 1.0)

    # initial cumulative-scaled positions
    x_sorted = np.sort(x_signals)
    deltas = np.diff(x_sorted)
    cumdist = np.concatenate([[0], np.cumsum(deltas)])
    scaled = (cumdist / cumdist[-1]) if cumdist[-1] else np.zeros_like(cumdist)
    x_texts = rail_start + scaled * rail_width

    # width-aware horizontal spread + proportional compress
    char_width = fsize * char_width_factor
    label_widths = np.array([len(lbl) * char_width for lbl in labels])
    left_ext  = 0.72 * label_widths
    right_ext = 0.05 * label_widths
    for i in range(1, n):
        required_gap = (left_ext[i] + right_ext[i - 1]) * 0.9
        if x_texts[i] < x_texts[i - 1] + required_gap:
            x_texts[i] = x_texts[i - 1] + required_gap
    if n >= 2 and x_texts[-1] > rail_end:
        used = x_texts[-1] - x_texts[0]
        avail = rail_end - x_texts[0]
        if used > avail and used > 0:
            x_texts = x_texts[0] + (x_texts - x_texts[0]) * (avail / used)

    # greedy lowest-free-rail stacking
    x_lo = x_texts - left_ext
    x_hi = x_texts + right_ext
    rail_ids = np.zeros(n, dtype=int)
    rail_extents: list[list[tuple[float, float]]] = []
    for i in range(n):
        for r in range(len(rail_extents) + 1):
            if r == len(rail_extents):
                rail_extents.append([])
            if not any(
                not (x_hi[i] < lo or x_lo[i] > hi)
                for lo, hi in rail_extents[r]
            ):
                rail_ids[i] = r
                rail_extents[r].append((x_lo[i], x_hi[i]))
                break

    x_texts = np.clip(x_texts, rail_start, rail_end)
    return x_texts, rail_ids, int(rail_ids.max()), rail_start, rail_end


# Using cumulative distance for anntations and separating clusters
def _draw_annotation_arrows_3(
    ax,
    annot_df,
    chr_col: str,
    label_col: str,
    offsets: dict,
    chr_max: dict,
    spread_width: float = 60e6,
    y_text_base: float = 0.25,
    max_rad: float = 0.35,
    y_tip: float = 0.0,
    fsize: float = 8,
) -> None:
    """Draw angled arrow annotations from gene-label text to signal positions.

    For each significant locus in *annot_df*, places the gene/SNP label text
    above the track and connects it to the corresponding scatter-plot point
    with a curved :class:`~matplotlib.patches.FancyArrowPatch` arrow.  Labels
    within the same chromosome are spread horizontally to avoid overlap; the
    arrow curvature direction (clockwise/counter-clockwise arc) is determined
    automatically from the sign of the horizontal displacement.

    The function operates on *ax* (the annotation sub-panel at the top of the
    figure) and returns ``None`` — it modifies the axes in place.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The annotation axes (top sub-panel, invisible background).
    annot_df : pandas.DataFrame
        Lead-SNP annotation DataFrame with columns *chr_col*, ``'x'``
        (cumulative x-axis position), and *label_col*.
    chr_col : str
        Name of the chromosome column.
    label_col : str
        Name of the column containing the text labels (gene symbols or rsIDs).
    offsets : dict
        Mapping of ``chromosome → cumulative_x_start`` as computed from the
        main track data.
    chr_max : dict
        Mapping of ``chromosome → max_position`` for each chromosome.
    spread_width : float, optional
        Horizontal spacing between adjacent labels within a cluster.
        Default ``60e6`` (60 Mb).
    y_tip : float, optional
        y-coordinate (in axes data units, 0–1) of the arrowhead tip.
        Default ``0.0``.
    y_text : float, optional
        y-coordinate of the label text anchor.  Default ``0.55``.
    """
    annot_df = annot_df.sort_values(by=[chr_col, "x"], key=natsort_keygen())

    x_signals = annot_df["x"].values
    labels    = annot_df[label_col].astype(str).values
    n    = len(x_signals)

    # ------------------------------------------------------------------
    # Genome span and rail bounds
    # ------------------------------------------------------------------
    genome_start = min(offsets.values())
    genome_end = max(
        offsets[c] + chr_max[c]
        for c in chr_max
    )
    genome_width = float(genome_end - genome_start) or 1.0

    # Reserve a wider margin on the LEFT than on the RIGHT.  With
    # ``rotation=315, ha="right"`` the visible glyphs extend up-AND-LEFT
    # of the anchor by roughly ``0.7 * text_width`` (the cos(45°)
    # projection of a 45°-rotated bbox).
    left_margin  = 0.06 * genome_width
    right_margin = 0.02 * genome_width
    rail_start = genome_start + left_margin
    rail_end   = genome_end   - right_margin
    rail_width = max(rail_end - rail_start, 1.0)

    # ------------------------------------------------------------------
    # Density-aware initial spacing — start each label near its signal x
    # by linearly mapping the sorted cumulative distribution onto the
    # available rail.  This preserves the "label is near its signal"
    # invariant before any collision-resolution work.
    # ------------------------------------------------------------------
    x_sorted = np.sort(x_signals)
    deltas = np.diff(x_sorted)
    cumdist = np.concatenate([[0], np.cumsum(deltas)])
    if cumdist[-1] == 0:
        scaled = np.zeros_like(cumdist)
    else:
        scaled = cumdist / cumdist[-1]
    x_texts = rail_start + scaled * rail_width

    # ------------------------------------------------------------------
    # Per-label horizontal footprint estimate (data units).  With
    # ``rotation=315, ha="right"`` the visible glyphs extend ~0.72 ×
    # text_width to the LEFT of the anchor and a small overshoot to
    # the right.  ``char_width`` is an empirical conversion from
    # characters at the requested font size to data-x units; the
    # default works well at fsize=8 in a 15-inch figure.
    # ------------------------------------------------------------------
    char_width = fsize * 2.5e5
    label_widths = np.array([len(lbl) * char_width for lbl in labels])
    left_ext  = 0.72 * label_widths
    right_ext = 0.05 * label_widths

    # ------------------------------------------------------------------
    # Phase 1 — Horizontal collision spread (along ONE rail).
    #
    # Walk labels in x-sorted order and push the next label rightward
    # whenever it would overlap the previous one's visible bbox.  This
    # is the classic "min_sep" collision strategy, but the required
    # gap is **width-aware**: long labels get a wider gap than short
    # rsIDs.  If the cumulative push runs past ``rail_end``, scale the
    # whole spread back proportionally to fit — at which point we
    # accept that some pairs will still overlap, and let phase 2 sort
    # them out by moving them to higher rails.
    # ------------------------------------------------------------------
    for i in range(1, n):
        required_gap = (left_ext[i] + right_ext[i - 1]) * 0.9
        min_next = x_texts[i - 1] + required_gap
        if x_texts[i] < min_next:
            x_texts[i] = min_next

    if n >= 2 and x_texts[-1] > rail_end:
        used = x_texts[-1] - x_texts[0]
        avail = rail_end - x_texts[0]
        if used > avail and used > 0:
            scale = avail / used
            x_texts = x_texts[0] + (x_texts - x_texts[0]) * scale

    # Recompute bbox extents after the spread / compress.
    x_lo = x_texts - left_ext
    x_hi = x_texts + right_ext

    # ------------------------------------------------------------------
    # Phase 2 — Vertical stacking for residual horizontal overlap.
    #
    # Greedy lowest-free-rail assignment: each label goes on the
    # lowest-numbered rail whose existing labels don't horizontally
    # overlap this one.  After phase 1 most labels will land on rail 0;
    # only the residual overlap from compressed-dense clusters fans
    # upward into rails 1, 2, ….  Worst case is O(n²) but ``n`` is
    # the number of lead SNPs (typically <200), so the cost is
    # negligible.
    # ------------------------------------------------------------------
    rail_ids = np.zeros(n, dtype=int)
    rail_extents: list[list[tuple[float, float]]] = []
    for i in range(n):
        for r in range(len(rail_extents) + 1):
            if r == len(rail_extents):
                rail_extents.append([])
            collides = any(
                not (x_hi[i] < lo or x_lo[i] > hi)
                for lo, hi in rail_extents[r]
            )
            if not collides:
                rail_ids[i] = r
                rail_extents[r].append((x_lo[i], x_hi[i]))
                break

    # Hard safety clamp — anchors must still sit inside the rail even
    # after spread + stacking (label widths near the rail edge could
    # otherwise push the anchor itself outside).
    x_texts = np.clip(x_texts, rail_start, rail_end)

    # ------------------------------------------------------------------
    # Compute label y positions from rail assignment.  Rails span the
    # vertical space between ``y_text_base`` and ``y_text_ceiling`` (just
    # below the annotation panel's ylim of ~1.0).  ``y_stack_step`` is
    # picked adaptively so all rails fit inside this band even when there
    # are many labels — preferred default is 0.06 per rail, but it
    # shrinks when needed.
    # ------------------------------------------------------------------
    max_rail = int(rail_ids.max()) if n else 0
    y_text_ceiling = 0.95
    available = max(y_text_ceiling - y_text_base, 0.05)
    y_stack_step = (
        min(0.15, available / max_rail) if max_rail > 0 else 0.06
    )
    y_texts = y_text_base + rail_ids * y_stack_step

    # ------------------------------------------------------------------
    # Draw arrows and labels
    # ------------------------------------------------------------------
    for x_sig, x_txt, y_txt, label in zip(x_signals, x_texts, y_texts, labels):
        dx  = x_txt - x_sig
        rad = np.clip(dx / (spread_width * 2), -max_rad, max_rad)

        arrow = FancyArrowPatch(
            (x_txt, y_txt),
            (x_sig, y_tip),
            arrowstyle="-|>",
            mutation_scale=12,
            lw=0.6,
            color="grey",
            alpha=0.5,
            # ``arc3,rad=<small>`` draws a single gentle curve from
            # label anchor to signal, rather than the previous L-shaped
            # "arc,angleA=...,armB=90" path that took a long horizontal
            # segment before turning down.  The horizontal segments from
            # different arrows overlapped chaotically in dense regions,
            # which is what made arrows look like they were crossing.
            #connectionstyle=f"arc3,rad={rad:.3f}",
            connectionstyle="arc,angleA=0,angleB=90,armA=5,armB=90,rad=2",
            transform=ax.transData,
        )
        ax.add_patch(arrow)

        ax.text(
            x_txt,
            y_txt + 0.001,
            str(label),
            #rotation=45,
            #ha="left",
            rotation=315,
            ha="right",
            va="bottom",
            fontsize=fsize,
            clip_on=False,
            color="black",
            fontstyle="italic",
            fontweight="regular",
        )


def _draw_annotation_arrows_4(
    ax,
    annot_df,
    chr_col: str,
    label_col: str,
    offsets: dict,
    chr_max: dict,
    spread_width: float = 60e6,
    y_text_base: float = 0.25,
    y_stack_step: float = 0.05,
    max_rad: float = 0.35,
    y_tip: float = 0.0,
    fsize: float = 8,
    rail_frac: float = 0.95,
    min_sep: float = 60e6,
    char_width_factor: float = 2.5e5,
) -> None:
    """
    Dense annotation renderer with:
        - density-aware cumulative spacing
        - width-aware relaxation
        - multi-rail stacking
        - curved arrows
        - adaptive ylim
    """

    # -------------------------------------------------------------
    # Sort annotations
    # -------------------------------------------------------------
    annot_df = (
        annot_df
        .sort_values(by=[chr_col, "x"], key=natsort_keygen())
        .reset_index(drop=True)
    )

    x_signals = annot_df["x"].to_numpy(dtype=float)
    labels = annot_df[label_col].astype(str).to_numpy()

    n = len(x_signals)
    if n == 0:
        return

    # -------------------------------------------------------------
    # Genome span
    # -------------------------------------------------------------
    genome_start = min(offsets.values())
    genome_end = max(
        offsets[c] + chr_max[c]
        for c in chr_max
    )
    genome_width = genome_end - genome_start

    # -------------------------------------------------------------
    # Density-aware cumulative spacing
    # -------------------------------------------------------------
    x_sorted = np.sort(x_signals)
    deltas = np.diff(x_sorted)
    cumdist = np.concatenate([
        [0],
        np.cumsum(deltas)
    ])

    if cumdist[-1] == 0:
        scaled = np.zeros_like(cumdist)
    else:
        scaled = cumdist / cumdist[-1]

    rail_width = genome_width * rail_frac
    rail_start = (
        genome_start
        + (genome_width - rail_width) / 2
    )
    x_texts = rail_start + scaled * rail_width

    # -------------------------------------------------------------
    # Width-aware label sizing
    # -------------------------------------------------------------
    char_width = fsize * char_width_factor

    label_widths = np.array([
        len(lbl) * char_width
        for lbl in labels
    ])

    # -------------------------------------------------------------
    # Relax labels using label widths
    # -------------------------------------------------------------
    for i in range(1, n):
        required_sep = max(
            (
                label_widths[i - 1]
                + label_widths[i]
            ) / 2,
            min_sep,
        )
        min_allowed = x_texts[i - 1] + required_sep

        if x_texts[i] < min_allowed:
            x_texts[i] = min_allowed

    # -------------------------------------------------------------
    # Multi-rail stacking
    # -------------------------------------------------------------
    # Density-aware rail assignment
    # -------------------------------------------------------------
    rail_ids = np.zeros(n, dtype=int)
    density_window = spread_width * 2

    for i in range(n):
        current_x = x_signals[i]
        nearby = np.where(
            np.abs(x_signals[:i] - current_x)
            < density_window
        )[0]

        if len(nearby) == 0:
            rail_ids[i] = 0
        else:
            used_rails = set(rail_ids[nearby])
            rail = 0
            while rail in used_rails:
                rail += 1
            rail_ids[i] = rail

    # -------------------------------------------------------------
    # Convert rails to y positions
    # -------------------------------------------------------------
    y_texts = np.array([
        y_text_base + r * y_stack_step
        for r in rail_ids
    ])

    # -------------------------------------------------------------
    # Arrow tip jitter
    # -------------------------------------------------------------
    jitter = np.linspace(
        -spread_width * 0.03,
        spread_width * 0.03,
        n,
    )

    # -------------------------------------------------------------
    # Draw arrows + labels
    # -------------------------------------------------------------
    for i, (
        x_sig,
        x_txt,
        y_txt,
        label,
    ) in enumerate(
        zip(
            x_signals,
            x_texts,
            y_texts,
            labels,
        )
    ):

        x_tip = x_sig + jitter[i]
        dx = x_txt - x_sig
        rad = np.clip(
            dx / (genome_width * 0.15),
            -max_rad,
            max_rad,
        )

        arrow = FancyArrowPatch(
            (x_txt, y_txt),
            (x_tip, y_tip),
            arrowstyle="-|>",
            mutation_scale=12,
            lw=0.6,
            color="grey",
            alpha=0.5,
            #connectionstyle=f"arc3,rad={rad}",
            connectionstyle="arc,angleA=5,angleB=90,armA=5,armB=90,rad=5",
            transform=ax.transData,
        )

        ax.add_patch(arrow)

        ax.text(
            #x_txt,
            x_txt + rail_ids[i] * spread_width * 0.12,
            y_txt + 0.001,
            str(label),
            rotation=290,
            ha="right",
            va="bottom",
            fontsize=fsize,
            clip_on=False,
            color="black",
            fontstyle="italic",
            fontweight="regular",
        )

    # -------------------------------------------------------------
    # Adaptive ylim
    # -------------------------------------------------------------
    max_rail = max(rail_ids)

    ax.set_ylim(
        y_tip - 0.05,
        y_text_base + (max_rail + 2) * y_stack_step
    )

# ---------------------------------------------------------------------------
# Public function
# ---------------------------------------------------------------------------

def plot_linearm(
    tracks: list,
    track_labels: Optional[list[str]] = None,
    annot_df: pd.DataFrame = None,
    annotate: bool = False,
    highlight: bool = False,
    highlight_thresh: float = 5e-8,
    highlight_color: str = 'brown',
    highlight_line: bool = False,
    highlight_line_color: str = 'grey',
    trim_pval: Optional[float] = None,
    logp: bool = True,
    label_col: Optional[str] = 'SNP',
    chr_order: Optional[list[str]] = None,
    chr_spacing: float = 9e6,
    track_heights: Optional[list[float]] = None,
    linear_track_spacing: float = 0.10,
    point_size: float = 8,
    colors: Optional[list[str]] = ['steelblue','silver'],
    sig_lines: Optional[list[dict]] = None,
    plt_name: Optional[str] = None,
    no_track_labels: bool = False,
    ylabel: Optional[str] = None,
    fig_format: Optional[str] = None,
    dpi: int = 300,
    figsize: Optional[list[float]] = [10, 4],
):
    """Core rendering engine for the multi-track stacked linear Manhattan plot.

    Builds a :class:`~matplotlib.figure.Figure` with one annotation sub-panel
    at the top (for gene/SNP labels and connecting arrows) and *n* data tracks
    below it, one per element of *tracks*.  All tracks share the same
    cumulative x-axis.

    This function is called by the higher-level :func:`plot_linear` wrapper and
    is not normally invoked directly.

    Parameters
    ----------
    tracks : list of pandas.DataFrame
        One DataFrame per GWAS trait.  Each must have columns ``CHR``,
        ``POS``, ``P``, and ``logP`` (when *logp* is ``True``).  The
        DataFrames are pre-processed internally (chromosome normalisation,
        highlighting, cumulative x-axis computation).
    track_labels : list of str, optional
        Y-axis label for each track, in the same order as *tracks*.
        Defaults to ``['Track 1', 'Track 2', …]``.
    annot_df : pandas.DataFrame, optional
        Annotation DataFrame (typically the hits summary table).  Must have
        columns ``CHR``, ``POS``, and *label_col*.  When provided, a dashed
        vertical guide line is drawn through every annotated position across
        all tracks, and gene/SNP label arrows are drawn in the top sub-panel.
    highlight : bool, optional
        If ``True``, variants within 500 kb of a lead SNP (as determined by
        :func:`~pycmplot.stats.get_highlight_snps`) are rendered in brown.
        Default ``False``.
    highlight_thresh : float, optional
        P-value threshold used for locus highlighting.  Default ``5e-8``.
    trim_pval : float, optional
        Reserved for future use; trimming is currently handled upstream in
        :func:`~pycmplot.io.get_sumstats_and_merged_sector_list`.
    logp : bool, optional
        If ``True``, plot –log₁₀(p) on the y-axis; requires a ``logP``
        column in each DataFrame.  Default ``True``.
    label_col : str, optional
        Column in *annot_df* to use as annotation text labels (e.g.
        ``'top_gene'``).  Default ``'label'``.
    chr_order : list of str, optional
        Chromosome display order.  Defaults to
        :data:`~pycmplot.constants.CHROM_ORDER` (autosomes 1–22, X, Y, MT).
    chr_spacing : float, optional
        Gap in base-pairs inserted between consecutive chromosomes on the
        x-axis.  Default ``9e6``.
    track_heights : list of float, optional
        Relative height ratios for the gridspec rows.  The first element
        controls the annotation sub-panel; subsequent elements control the
        data tracks.  When ``None``, the annotation panel is given a weight
        of 1 and each data track a weight of 3.
    linear_track_spacing : float, optional
        Vertical ``hspace`` between tracks as a fraction of average track
        height.  Default ``0.10``.
    point_size : float, optional
        Scatter-plot point size passed to :func:`matplotlib.axes.Axes.scatter`.
        Default ``5``.
    colors : list of str, optional
        Two alternating matplotlib colour strings for even/odd chromosomes.
        Default ``['steelblue', 'orange']``.
    sig_lines : list of dict, optional
        One ``{'genome': float, 'suggestive': float}`` dict per track.
        ``'genome'`` values are drawn as red dashed lines; ``'suggestive'``
        values as grey dashed lines.
    plt_name : str, optional
        Full output file path (including extension).  When provided the figure
        is saved to disk.
    no_track_labels : bool, optional
        If ``True``, per-track labels are suppressed.  Default ``False``.
    ylabel : str, optional
        Shared y-axis label rendered on the left of the figure.  Use this
        to set a sensible label for non-p-value statistics such as iHS,
        FST or XP-EHH (e.g. ``ylabel="iHS"``).  When ``None`` (the
        default), the label is derived automatically:
        ``"-log₁₀(p-value)"`` when *logp* is ``True``, otherwise the
        p-value column name (``"P"``).
    fig_format : str, optional
        Output image format (``'png'``, ``'pdf'``, ``'svg'``).  Inferred
        from *plt_name*'s extension when ``None``.
    dpi : int, optional
        Output resolution in dots per inch.  Default ``300``.
    figsize : tuple of (float, float), optional
        Figure dimensions ``(width, height)`` in inches.  Default ``(15, 9)``.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The completed figure.
    axes : list of matplotlib.axes.Axes
        All axes in the figure: ``axes[0]`` is the annotation sub-panel;
        ``axes[1:]`` are the data-track axes in the same order as *tracks*.

    See Also
    --------
    plot_linear :
        User-facing wrapper that extracts DataFrames from the
        ``sumstats_loaded`` dict and resolves output paths before calling
        this function.
    _draw_annotation_arrows :
        Called internally to render gene/SNP label arrows in ``axes[0]``.
    pycmplot.stats.get_highlight_snps :
        Called internally when *highlight* is ``True``.
    """

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

        # Fast path: io.py's loader stores CHR as a Categorical with
        # ``CHROM_ORDER`` as the canonical category set.  When that's the
        # case we can derive ``chr_idx`` from ``cat.codes`` directly,
        # skipping ~0.7 s of redundant string normalisation on a 500K-row
        # frame (and proportionally more on larger frames).
        chr_series = df[chr_col]
        is_canonical_cat = (
            isinstance(chr_series.dtype, pd.CategoricalDtype)
            and list(chr_series.cat.categories) == list(chr_order)
        )
        if is_canonical_cat:
            df["chr_idx"] = chr_series.cat.codes.to_numpy()
        else:
            # Defensive normalisation for callers that bypass io.py (e.g.
            # users feeding hand-built DataFrames straight to plot_linearm).
            df[chr_col] = (
                chr_series
                .astype(str)
                .str.replace("chr", "", regex=False)
                .str.upper()
                .replace({"23": "X", "24": "Y", "M": "MT", "MTDNA": "MT"})
            )
            df = df[df[chr_col].isin(chr_order)]
            df["chr_idx"] = df[chr_col].map(chr_to_idx)

        return df.sort_values(["chr_idx", pos_col])

    tracks = [_prep(df) for df in tracks]

    chr_maxes = []

    for df in tracks:
        # ``observed=True`` skips categorical levels with no rows in this
        # particular track; ``s.get(c, 0)`` below handles those cases.
        # Without this, pandas returns ``NA`` for the empty groups, which
        # later propagates as ``TypeError: boolean value of NA is
        # ambiguous`` in the ``xmax = max(...)`` reduction.
        chr_maxes.append(
            df.groupby(chr_col, observed=True)[pos_col].max()
        )

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
            [s.get(c, 0) for s in chr_maxes]
        )
        chr_max[c] = max_pos
        offsets[c] = offset
        offset += max_pos + chr_spacing

    def _add_cum(df):
        df = df.copy()
        df["x"] = (
            df[pos_col].to_numpy()
            + df[chr_col].map(offsets).to_numpy()
        )
        return df

    tracks = [_add_cum(df) for df in tracks]
    if annot_df is not None:
        annot_df = _add_cum(annot_df)

    # ------------------------------------------------------------------
    # Figure layout
    # ------------------------------------------------------------------
    n_tracks = len(tracks)

    if track_heights is None:
        if annotate and annot_df is not None and len(annot_df.index) > 0:
            # Pre-compute how many rails the annotation labels will need so
            # we can size the annotation panel proportionally.  Without
            # this step the panel always has the same gridspec ratio (1)
            # and stacked rails end up physically squashed even when
            # ``y_stack_step`` is generous in axes-relative units.
            try:
                _, _, _max_rail, _, _ = _layout_annotation_labels(
                    annot_df,
                    chr_col=chr_col,
                    label_col=label_col,
                    offsets=offsets,
                    chr_max=chr_max,
                    fsize=8,
                )
            except Exception:
                _max_rail = 0
            # Annotation panel sizing is driven by BOTH the rail count
            # (how many vertical levels we'll stack onto) AND the total
            # label count (so dense single-rail clusters still get more
            # physical height; otherwise the rotated text overlaps
            # visually even when bbox-wise it doesn't, since the rotation
            # mixes axes units).  Cap at 5× a data track so the
            # annotation panel doesn't completely drown out the data
            # tracks underneath.
            n_labels = len(annot_df.index)
            annot_ratio = min(
                5.0,
                1.0 + 0.5 * max(_max_rail, 0) + 0.03 * n_labels,
            )
            track_heights = [annot_ratio] + [3] * n_tracks
        elif annotate:
            track_heights = [1] + [3] * n_tracks
        else:
            track_heights = [3] * n_tracks

    fig = plt.figure(figsize=figsize)
    gs_tracks = n_tracks + 1 if annotate else n_tracks
    gs = fig.add_gridspec(
        gs_tracks, 1,
        height_ratios=track_heights,
        hspace=linear_track_spacing,
    )

    if annotate:
        ax_annot = fig.add_subplot(gs[0, 0])
        axes = [ax_annot]
        for i in range(n_tracks):
            axes.append(fig.add_subplot(gs[i + 1, 0], sharex=ax_annot))        
    else:
        first_ax = fig.add_subplot(gs[0, 0])
        axes = [first_ax]        
        for i in range(n_tracks-1):
            axes.append(fig.add_subplot(gs[i+1, 0], sharex=first_ax))

    # Per-track highlight colours from tab20 colormap
    cmap = plt.get_cmap("tab20")
    hex_colors = [mcolors.to_hex(cmap(i / n_tracks)) for i in range(n_tracks)]

    # ------------------------------------------------------------------
    # Plot data tracks
    # ------------------------------------------------------------------
    t_labels = track_labels or [f"Track {i+1}" for i in range(n_tracks)]

    if annotate:
        loop_axes = axes[1:]
    else:
        loop_axes = axes

    left_pad = chr_spacing * 0.2
    xmax = max(offsets[c] + chr_max[c] for c in chr_order)

    for i, (ax, df, t_label, h_color) in enumerate(
        zip(loop_axes, tracks, t_labels, hex_colors)
    ):
        logger.info("Plotting: %s ...", t_label)

        # Sanity filter.  When plotting -log10(p) we assume the P column is a
        # p-value and drop non-positive entries (a common product of upstream
        # NaN -> 0 coercion or mistyped headers).  When plotting raw values
        # we keep every ro(w, since selection statistics such as iHS, XP-EHH
        # or Fay & Wu's H can legitimately be negative.  color_cycle is
        # computed *after* filtering so the color array always aligns with
        # the scatter x/y arrays.
        if logp:
            df = df[df[p_col] >= 0]

        # Render markers via ``ax.plot`` (one Line2D per chromosome) rather
        # than ``ax.scatter`` (which builds a PathCollection: one Path per
        # point, with a per-point ``should_simplify`` check and a 4×N RGBA
        # array).  At 10 M variants this cuts the matplotlib draw step from
        # several minutes to a few seconds while producing visually
        # identical output (the rasterised PNG/PDF backend treats the
        # vertices the same way either way).  Marker size is converted from
        # the scatter ``s`` unit (points²) to the plot ``markersize`` unit
        # (points) via ``sqrt`` so existing point-size values keep meaning.
        y_col = "logP" if logp else p_col
        x_all = df["x"].to_numpy()
        y_all = df[y_col].to_numpy()
        chr_idx_all = df["chr_idx"].to_numpy()

        # ``point_size`` is in scatter units (points²); convert to plot's
        # markersize (points) so existing user-supplied values still mean
        # the same visual size.  Guard against ``None`` propagated from
        # callers that don't set the parameter explicitly.
        _ps = 4.0 if point_size is None else float(point_size)
        markersize = float(np.sqrt(_ps))
        for i_chr in np.unique(chr_idx_all):
            mask = chr_idx_all == i_chr
            ax.plot(
                x_all[mask], y_all[mask],
                marker="o",
                linestyle="none",
                color=colors[int(i_chr) % len(colors)],
                markersize=markersize,
                markeredgecolor="none",
                rasterized=True,
            )
        if y_all.size:
            ax.set_ylim(y_all.min(), y_all.max()+2)

        # Track labels — rendered vertically in the right margin (just past
        # the last chromosome), orthogonal to the track-stacking direction.
        # This keeps the label out of the data region, survives arbitrarily
        # tight packing (``linear_track_spacing=0`` onward), and never
        # collides with the shared y-axis label on the left.
        if not no_track_labels:
            ax.text(
                1.005, 0.5, t_label,
                transform=ax.transAxes,
                ha="left", va="center",
                rotation=-90,
                fontsize=10,
            )

        if highlight:
            sig = df[df["in_locus"]]
            if not sig.empty:
                sig_y = sig["logP"] if logp else sig[p_col]
                # Highlight points need an explicit ``zorder`` higher than
                # the background scatter.  The background is drawn with
                # ``ax.plot`` (Line2D, default zorder=2), while
                # ``ax.scatter`` returns a PathCollection whose default
                # zorder is **1**.  Without the explicit bump below the
                # highlight ends up *underneath* the background despite
                # being called later in code.  ``zorder=3`` puts the
                # highlight cleanly on top of both background scatter and
                # any axhline / axvline guides.
                ax.scatter(
                    sig["x"].to_numpy(),
                    sig_y.to_numpy(),
                    s=point_size,
                    marker="o",
                    color=highlight_color,
                    edgecolors="none",
                    zorder=3,
                    rasterized=True,
                )
                # Vertical lines across all data tracks at highlight positions
                if highlight_line:
                    for x in sig["x"].values:
                        for _ax in loop_axes:
                            _ax.axvline(x, color=highlight_line_color, alpha=0.1, linewidth=0.2,
                                    linestyle="--", zorder=0)

        if sig_lines is not None and i < len(sig_lines):
            sl = sig_lines[i]
            if "genome" in sl:
                ax.axhline(y=sl["genome"], color="red", linestyle="--", linewidth=0.5)
            if "suggestive" in sl:
                ax.axhline(y=sl["suggestive"], color="blue", linestyle="--", linewidth=0.4)

        ax.spines[["top", "right"]].set_visible(False)

        ax.set_xlim(-left_pad, xmax)

    # ------------------------------------------------------------------
    # Annotation track
    # ------------------------------------------------------------------
    if annotate and annot_df is not None:
        less_than_spread_width = []
        for chr, df in annot_df.groupby(chr_col):
            df_chr = df[df[chr_col]==chr]
            differences = np.diff(df_chr['POS']).tolist()
            less_than_spread_width.append(list(filter(lambda x: x < 60e6, differences)))

        if len(less_than_spread_width) < 5:
            _draw_annotation_arrows(
                ax_annot,
                annot_df,
                chr_col=chr_col,
                label_col=label_col,
                offsets=offsets,
                chr_max=chr_max,
                spread_width=60e6,
            )
        else:
            _draw_annotation_arrows_3(
                ax=ax_annot,
                annot_df=annot_df,
                chr_col=chr_col,
                label_col=label_col,
                offsets=offsets,
                chr_max=chr_max,
                spread_width=60e6,
                y_text_base=0.25,
                max_rad=0.35,
                y_tip=0.0,
            )
            """
            _draw_annotation_arrows_4(
                ax=ax_annot,
                annot_df=annot_df,
                chr_col=chr_col,
                label_col=label_col,
                offsets=offsets,
                chr_max=chr_max,
                spread_width=60e6,
                y_text_base=0.55,
                max_rad=0.35,
                y_tip=0.0,                
            )
            """
        
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
        ax.spines[["top", "right", "bottom"]].set_visible(False)

    # Explicit margins keep the shared y-axis label (added below as
    # ``fig.text``) in its own reserved strip on the left of the figure.
    # We avoid ``tight_layout`` here because the annotation panel shares
    # ``sharex`` with the data tracks inside a custom ``GridSpec``, which
    # triggers a matplotlib "not compatible with tight_layout" warning
    # without actually helping the layout.
    # Reserve margins explicitly: the left strip holds the shared y-axis
    # label, the right strip holds the per-track labels placed just past
    # the last chromosome.  ``hspace`` is forwarded verbatim so users can
    # set ``linear_track_spacing=0`` to stack tracks flush against each
    # other — the right-side track labels remain readable because they
    # are orthogonal to the stacking direction.
    fig.subplots_adjust(
        left=0.09,
        right=0.94,
        top=0.90,
        bottom=0.15,
        hspace=linear_track_spacing,
    )

    if ylabel is None:
        ylabel_text = "-log\u2081\u2080(P)" if logp else p_col
    else:
        ylabel_text = ylabel

    fig.text(
        0.025, 0.5,
        ylabel_text,
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
    point_size: Optional[float] = 8,
    highlight: bool = False,
    highlight_thresh: float = 5e-8,
    highlight_color: str = 'brown',
    highlight_line: bool = False,
    highlight_line_color: str = 'grey',    
    hits_table: Optional[pd.DataFrame] = None,
    annotate: str = None,
    label_col: Optional[str] = None,
    chr_spacing: Optional[float] = 9e6,
    linear_track_spacing: Optional[float] = None,
    colors: list[str] = ['steelblue','silver'],
    signif_lines: Optional[dict] = None,
    plot_title: Optional[str] = None,
    no_track_labels: bool = False,
    ylabel: Optional[str] = None,
    dpi: Optional[int] = None,
    output_format: Optional[str] = 'png',
    output_dir: Optional[str] = '.',
    figsize: Optional[list[float]] = [10, 4],
):
    """Generate a multi-track stacked linear Manhattan plot.

    This is the primary user-facing entry point for linear Manhattan plots.
    It extracts DataFrames and labels from *sumstats_loaded*, resolves output
    file paths, then delegates all rendering to :func:`plot_linearm`.

    Parameters
    ----------
    sumstats_loaded : dict
        Mapping of ``label → [DataFrame, n_chroms]`` as returned by
        :func:`~pycmplot.io.get_sumstats_and_merged_sector_list`.  One
        track is created per key; the DataFrames must have canonical columns
        ``CHR``, ``POS``, ``P``, and ``logP`` (when *logp* is ``True``).
    trim_pval : float, optional
        Reserved; trimming is applied upstream.  Default ``None``.
    track_heights : list of float, optional
        Comma-parsed relative track heights (e.g. ``[2.0, 2.0, 1.5]``).
        Passed directly to :func:`plot_linearm` as the gridspec height
        ratios for the data tracks.  When ``None``, all tracks are equal.
    logp : bool, optional
        Plot –log₁₀(p) on the y-axis.  Default ``False``.
    point_size : float, optional
        Scatter-plot point size.  Default uses :func:`plot_linearm`'s
        default (``5``).
    highlight : bool, optional
        Render significant-locus variants in brown.  Default ``False``.
    highlight_thresh : float, optional
        P-value threshold for locus highlighting.  Default ``5e-8``.
    hits_table : pandas.DataFrame, optional
        Annotation DataFrame (hits summary table from
        :func:`~pycmplot.annotation.get_hits_summary_table`).  Must contain
        columns ``CHR``, ``POS``, and *label_col*.  Passed to
        DataFrame suppresses annotations.
    label_col : str, optional
        Column in *hits_table* to use as annotation text (e.g.
        ``'top_gene'``).  Default ``None`` (falls back to :func:`plot_linearm`
        default ``'label'``).
    chr_spacing : float, optional
        Horizontal gap between chromosomes in base-pairs.  Default ``9e6``.
    linear_track_spacing : float, optional
        Vertical space between tracks as a fraction of average track height.
        Default ``0.10``.
    colors : list of str, optional
        Two alternating chromosome colours.  Default ``['steelblue', 'silver']``.
    signif_lines : list of dict, optional
        One ``{'genome': float, 'suggestive': float}`` dict per track, in
        the same order as *sumstats_loaded*.  Produced by
        :func:`~pycmplot.io.get_sumstats_and_merged_sector_list`.
    plot_title : str, optional
        Human-readable title used as the output file-name stem.  Passed to
        :func:`~pycmplot.io.get_output_paths`.
    no_track_labels : bool, optional
        Suppress per-track labels.  Default ``False``.
    ylabel : str, optional
        Override the shared y-axis label (left margin).  Useful for
        non-p-value statistics such as iHS, F_ST or XP-EHH (e.g.
        ``ylabel="iHS"``).  When ``None`` (the default), the label is
        ``"-log₁₀(p-value)"`` if *logp* is ``True`` and ``"P"`` otherwise.
    dpi : int, optional
        Output resolution in dots per inch.  Default ``300``.
    output_format : str, optional
        Image format (``'png'``, ``'pdf'``, ``'svg'``, ``'jpg'``).
        Default ``'png'``.
    output_dir : str or pathlib.Path, optional
        Directory in which to save the output files.  Default ``'.'``.
    figsize : tuple of (float, float), optional
        Figure dimensions ``(width, height)`` in inches.  Default ``(15, 9)``.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The completed figure.
    axes : list of matplotlib.axes.Axes
        All axes: ``axes[0]`` is the annotation sub-panel; ``axes[1:]``
        are the per-track data axes.

    See Also
    --------
    plot_linearm :
        The underlying rendering engine called by this function.
    pycmplot.io.get_sumstats_and_merged_sector_list :
        Produces *sumstats_loaded* and *signif_lines*.

    Examples
    --------
    >>> from pycmplot.plotting.linear import plot_linear
    >>> fig, axes = plot_linear(
    ...     sumstats_loaded=loaded,
    ...     logp=True,
    ...     highlight=True,
    ...     hits_table=hits,
    ...     label_col="top_gene",
    ...     signif_lines=sig_lines,
    ...     plot_title="RBC_Traits",
    ...     output_dir="./results",
    ... )
    """

    dfs      = [v[0] for v in sumstats_loaded.values()]
    t_labels = list(sumstats_loaded.keys())

    if not track_heights:
        t_heights = None
    else:
        t_heights = [float(x) for x in track_heights]

    label = 'SNP'
    if annotate:
        label = get_annotation_column(
            annotate=annotate, 
            hits_table=hits_table,
            label_col=label_col
        )

    print(label)

    # plot name
    (
        plt_name, 
        table_out,
        plt_base,
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
        highlight_color = highlight_color,
        highlight_line = highlight_line,
        highlight_line_color = highlight_line_color,
        annotate=annotate,        
        annot_df=hits_table if hits_table is not None and not hits_table.empty else None,
        label_col=label,
        chr_spacing=chr_spacing,
        track_heights=t_heights,
        linear_track_spacing=linear_track_spacing,
        colors=colors,
        sig_lines=signif_lines,
        plt_name=plt_name,
        no_track_labels = no_track_labels,
        ylabel=ylabel,
        dpi=dpi,
        fig_format=output_format,
        figsize=figsize,
    )

    return axes
