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

* :func:`_draw_annotation_arrows` — places angled
  :class:`~matplotlib.patches.FancyArrowPatch` arrows from spread gene
  labels down to their corresponding signal positions.
* :func:`_draw_annotation_arrows_multirail` — places angled
  :class:`~matplotlib.patches.FancyArrowPatch` arrows from spread gene
  labels down to their corresponding signal positions with multirail capability 
  and single sort + rank-reassignment to avoid arrow crossing.
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
import matplotlib.transforms as mtransforms
from natsort import natsort_keygen

from pycmplot.constants import CHROM_ORDER
from pycmplot.io import get_output_paths
from pycmplot.annotation import get_annotation_column

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Annotation helpers
# ---------------------------------------------------------------------------
# Using cumulative distance for anntations and separating clusters
def _draw_annotation_arrows(
    ax,
    annot_df,
    chr_col: str,
    label_col: str,
    offsets: dict,
    chr_max: dict,
    spread_width: float = 60e6,
    isolation_threshold: float = 80e6,
    stack_threshold: float = 10e6,
    y_text_base: float = 0.55,
    y_stack_step: float = 0.02,
    y_tip: float = 0.0,
    asize: float = 8,
) -> None:
    """
    Single-rail annotation renderer with tiered label placement,
    chromosome-boundary spreading, cumulative-distance stacking, and
    curved arrows.

    Layout pipeline
    ---------------
    **Chromosome boundary detection**
        For each pair of adjacent chromosomes, the genomic gap between
        the last annotation on chromosome *k* and the first annotation
        on chromosome *k+1* is computed.  If the gap is narrower than
        ``spread_width``, both boundary annotations receive an
        ``x_bound`` value (negative for the left boundary label,
        positive for the right) so that the spreading step knows to
        push them apart.

    **Tiered label x-placement** (per chromosome)
        Each label is assigned an x position according to a two-tier
        rule based on its proximity to neighbours and chromosome
        boundaries:

        *Tier 1 — isolated labels*
            If the nearest neighbour and the nearest chromosome boundary
            are both at least ``isolation_threshold`` bp away, the label
            is placed directly above its signal (``x_text = x_signal``).

        *Tier 2 — dense / boundary labels*
            Labels that are too close to a neighbour or a chromosome
            boundary are spread.  Boundary labels are shifted left or
            right by ``spread_width / 3`` depending on which side of
            the boundary they fall.  All remaining dense labels are
            placed at regular ``spread_width`` intervals starting from
            the leftmost dense signal, with the interval reduced
            iteratively if the chromosome is narrower than
            ``spread_width``.

    **Cumulative-distance stacking** (per chromosome)
        After x positions are assigned, y positions are computed by
        walking labels left-to-right.  If two adjacent labels are
        closer than ``stack_threshold`` in x, the right label is
        stacked above the left by ``y_stack_step`` scaled by how much
        closer than ``stack_threshold`` the pair is.  Labels that are
        far enough apart reset to ``y_text_base``.

    **Arrows**
        A ``FancyArrowPatch`` with ``arc3`` connectionstyle connects
        each label's text position to its signal tip at ``y_tip``.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes on which annotations are drawn.
    annot_df : pd.DataFrame
        Annotation table.  Must contain ``chr_col``, ``"x"``
        (cumulative genomic position in bp), and ``label_col``.
        Rows are processed in chromosome then position order.
    chr_col : str
        Column name for chromosome identifiers.
    label_col : str
        Column name for annotation labels.
    offsets : dict
        Mapping of chromosome name → cumulative start offset (bp).
    chr_max : dict
        Mapping of chromosome name → chromosome length (bp).
    spread_width : float, optional
        Genomic interval (bp) used as the base spacing between dense
        labels and as the threshold for chromosome-boundary detection.
        Default 60e6.
    isolation_threshold : float, optional
        Minimum distance (bp) to both the nearest neighbour and the
        nearest chromosome boundary for a label to be considered
        isolated (Tier 1).  Default 80e6.
    stack_threshold : float, optional
        Maximum x-distance (bp) between adjacent labels below which
        the right label is stacked upward rather than reset to
        ``y_text_base``.  Default 10e6.
    y_text_base : float, optional
        Axes-fraction y-coordinate for the baseline label row.
        Default 0.55.
    y_stack_step : float, optional
        Base axes-fraction increment applied when stacking a label
        above its left neighbour.  The actual increment is scaled by
        ``1 - cum_dist / stack_threshold`` so tighter pairs stack
        higher.  Default 0.02.
    y_tip : float, optional
        Axes-fraction y-coordinate for arrow tips (signal level).
        Default 0.0.
    asize : float, optional
        Font size (pt) for rendered label text.  Default 8.
    """

    from natsort import natsorted

    annot_df = annot_df.sort_values(by=[chr_col, "x"], key=natsort_keygen())
    last_xtext = 0 - spread_width

    # ------------------------------------------------------------------
    # Chromosome boundary detection
    # ------------------------------------------------------------------
    # For each adjacent chromosome pair, mark the boundary annotations
    # with x_bound if the inter-chromosome gap is narrower than
    # spread_width.  x_bound is negative for the left boundary label
    # and positive for the right, encoding both the direction and the
    # magnitude of the gap for use in the spreading step below.
    chr_order = natsorted(set(annot_df[chr_col]))
    annot_df["x_bound"] = np.nan

    for chr1, chr2 in zip(chr_order[:-1], chr_order[1:]):
        df1 = annot_df[annot_df[chr_col] == chr1]
        df2 = annot_df[annot_df[chr_col] == chr2]

        if df1.empty or df2.empty:
            continue

        idx1 = df1["x"].idxmax()   # last annotation in chr1
        idx2 = df2["x"].idxmin()   # first annotation in chr2

        gap = (
            annot_df.loc[idx2, "x"]
            - annot_df.loc[idx1, "x"]
        )

        if gap <= spread_width:
            annot_df.loc[idx1, "x_bound"] = -gap
            annot_df.loc[idx2, "x_bound"] =  gap

    # ------------------------------------------------------------------
    # Per-chromosome label placement and rendering
    # ------------------------------------------------------------------
    for chr_name, df_chr in annot_df.groupby(chr_col, sort=False):
        df_chr    = df_chr.sort_values("x")
        chr_start = offsets[chr_name]
        chr_end   = offsets[chr_name] + chr_max[chr_name]
        chr_range = chr_end - chr_start

        x_signals = df_chr["x"].values
        x_bounds  = df_chr["x_bound"].values
        labels    = df_chr[label_col].values
        n         = len(x_signals)

        # --------------------------------------------------------------
        # Tiered x-position assignment
        # --------------------------------------------------------------
        # Tier 1 (isolated): label sits directly above its signal.
        # Tier 2 (dense / boundary): label is spread away from its
        #   signal.  Boundary labels are shifted by spread_width/3 in
        #   the direction away from the adjacent chromosome.  All other
        #   dense labels receive evenly-spaced positions at spread_width
        #   intervals, shifted further right if they would overlap the
        #   previous chromosome's last label.
        x_texts = []
        sw  = spread_width
        pad = sw / int(str(sw)[:2]) / 2

        for k, (x_sig, x_bnd) in enumerate(zip(x_signals, x_bounds)):
            neighbours         = np.delete(x_signals, k)
            min_neighbor_dist  = (
                np.min(np.abs(neighbours - x_sig))
                if len(neighbours) else np.inf
            )
            min_dist = min(min_neighbor_dist, np.abs(x_bnd))

            if min_dist >= isolation_threshold:
                # Tier 1: isolated — place directly above signal
                x_texts.append(x_sig)
            else:
                if not np.isnan(x_bnd):
                    # Chromosome boundary: shift away from the boundary
                    x_sig = x_sig - sw / 3 if x_bnd < 0 else x_sig + sw / 3
                    x_texts.append(x_sig)
                else:
                    # Tier 2: dense interior — defer to spread pass
                    x_texts.append(None)

        # Spread deferred (Tier 2) labels at regular intervals
        spread_indices = [k for k, v in enumerate(x_texts) if v is None]
        if spread_indices:
            # Reduce spread interval until it fits within the chromosome
            while sw > chr_range and sw > pad:
                sw -= pad

            sig_start = x_signals[spread_indices[0]]
            xmin      = sig_start - sw
            positions = np.arange(
                xmin,
                xmin + len(spread_indices) * sw,
                sw,
            )

            # Shift right if positions overlap the previous chromosome
            while positions[0] <= last_xtext:
                positions = positions + sw

            for j, k in enumerate(spread_indices):
                x_texts[k] = positions[j]

        # --------------------------------------------------------------
        # Cumulative-distance y stacking
        # --------------------------------------------------------------
        # Walk labels left-to-right.  If the x gap to the previous label
        # is below stack_threshold, stack the current label upward by
        # y_stack_step scaled by how tight the pair is.  Otherwise reset
        # to y_text_base.
        y_texts = [y_text_base] * n

        for k in range(1, n):
            cum_dist = abs(x_texts[k] - x_texts[k - 1])
            if cum_dist <= stack_threshold:
                y_texts[k] = (
                    y_texts[k - 1]
                    + y_stack_step
                    + y_stack_step * (1 - cum_dist / stack_threshold)
                )
            else:
                y_texts[k] = y_text_base

        # --------------------------------------------------------------
        # Draw arrows and labels
        # --------------------------------------------------------------
        for x_sig, x_txt, y_txt, label in zip(
            x_signals, x_texts, y_texts, labels
        ):
            # Straight arrows: this function handles sparse annotations
            # only; curvature is intentionally zero.  Dense annotations
            # are handled by _draw_annotation_arrows_multirail which uses
            # arc3 connectionstyle with non-zero rad.
            rad = 0

            arrow = FancyArrowPatch(
                (x_txt, y_txt),
                (x_sig, y_tip - 0.05),
                arrowstyle="-|>",
                mutation_scale=8,
                lw=0.6,
                color="grey",
                alpha=0.5,
                connectionstyle=f"arc3,rad={rad}",
                transform=ax.transData,
            )
            ax.add_patch(arrow)

            ax.text(
                x_txt,
                y_txt + 0.02,
                str(label),
                rotation=90,
                ha="center",
                va="bottom",
                fontsize=asize,
                clip_on=False,
                color="black",
                fontstyle="italic",
                fontweight="regular",
            )

        last_xtext = max(x_texts)

#----- multi-rail annotations
def _draw_annotation_arrows_multirail(
    ax,
    annot_df,
    chr_col: str,
    label_col: str,
    offsets: dict,
    chr_max: dict,
    spread_width: float = 60e6,
    y_text_base: float = 0.25,
    y_stack_step: float = 0.1,
    max_rad: float = 0.35,
    y_tip: float = 0.0,
    fsize: float = 8,
    rail_frac: float = 0.95,
    min_sep: float = 6e6,
    asize: float = 8,
) -> None:
    """
    Dense annotation renderer with relaxation-driven multi-rail
    stacking, linspace rank-reassignment, curved arrows, and adaptive
    ylim.

    Layout pipeline
    ---------------
    1. **Relaxation pass**:
       All labels are sorted by ``x_signal`` and a bidirectional
       relaxation loop enforces ``min_sep`` between every adjacent pair.
       Labels are pushed apart until no two are closer than ``min_sep``.
       The relaxed positions are stored as ``x_relaxed``.

    2. **Rail assignment**:
       Each label's rail is determined by how far its relaxed position
       drifted from its signal::

           drift    = x_relaxed - x_signal
           rail_id  = floor(drift / (rail_width / max_rails))

       Labels that drifted further are assigned to higher rails,
       distributing stacking proportionally to local density.  No
       per-rail queue processing or max_drift threshold is needed.

    3. **linspace rank-reassignment**:
       The layout is sorted by ``x_signal``.  ``x_text`` values are
       replaced with ``np.linspace(rail_start, rail_end, n)`` — evenly
       spaced slots spanning the full rail.  This guarantees:

       - ``x_text`` rank == ``x_signal`` rank → no arrow crossings by
         construction.
       - Full rail coverage regardless of ``rail_frac`` or signal
         clustering.
       - Even slot spacing = ``rail_width / (n − 1)``.

    4. **Rendering pass**:
       ``rail_id`` is read here for the first time to compute
       ``y = y_text_base + rail_id * y_stack_step``.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    annot_df : pd.DataFrame
        Annotation table. Must contain ``chr_col``, ``"x"`` (cumulative
        genomic position in bp), and ``label_col``.
    chr_col : str
        Column name for chromosome identifiers.
    label_col : str
        Column name for annotation labels.
    offsets : dict
        Mapping of chromosome name → cumulative start offset (bp).
    chr_max : dict
        Mapping of chromosome name → chromosome length (bp).
    spread_width : float, optional
        Genomic window (bp) used for arrow tip jitter. Default 60e6.
    y_text_base : float, optional
        Axes-fraction y-coordinate for rail 0 labels. Default 0.25.
    y_stack_step : float, optional
        Axes-fraction increment per rail. Default 0.1.
    max_rad : float, optional
        Maximum arc curvature for ``FancyArrowPatch``. Default 0.35.
    y_tip : float, optional
        Axes-fraction y-coordinate for arrow tips. Default 0.0.
    fsize : float, optional
        Font size (pt) used to estimate label widths. Default 8.
    rail_frac : float, optional
        Fraction of genome width occupied by the label rail. Default 0.95.
    min_sep : float, optional
        Minimum genomic separation (bp) between any two adjacent label
        centres. Default 6e6.
    asize : float, optional
        Font size (pt) for rendered label text. Default 8.
    """

    # ------------------------------------------------------------------
    # Deduplication
    # ------------------------------------------------------------------
    annot_df = annot_df.drop_duplicates(subset=[chr_col, "x", label_col])

    # ------------------------------------------------------------------
    # Sort annotations
    # ------------------------------------------------------------------
    annot_df = (
        annot_df
        .sort_values(by=[chr_col, "x"], key=natsort_keygen())
        .reset_index(drop=True)
    )

    x_signals = annot_df["x"].to_numpy(dtype=float)
    labels    = annot_df[label_col].astype(str).to_numpy()
    n         = len(x_signals)

    if n == 0:
        return

    # ------------------------------------------------------------------
    # Genome span and rail bounds
    # ------------------------------------------------------------------
    genome_start = min(offsets.values())
    genome_end   = max(offsets[c] + chr_max[c] for c in chr_max)
    genome_width = genome_end - genome_start
    rail_width   = genome_width * rail_frac
    rail_start   = genome_start + (genome_width - rail_width) / 2
    rail_end     = rail_start + rail_width

    # ------------------------------------------------------------------
    # Auto char_width from axes geometry
    # ------------------------------------------------------------------
    # For vertical text (rotation=90°) the horizontal footprint of every
    # label is one character wide regardless of string length.
    try:
        fig         = ax.get_figure()
        renderer    = fig.canvas.get_renderer()
        ax_bbox     = ax.get_window_extent(renderer=renderer)
        xmin, xmax  = ax.get_xlim()
        px_per_bp   = ax_bbox.width / (xmax - xmin)
        char_width  = 0.6 * fsize * (fig.dpi / 72.0) / px_per_bp
    except Exception:
        char_width  = genome_width * 0.01 / fsize

    # ------------------------------------------------------------------
    # 1. Relaxation pass
    # ------------------------------------------------------------------
    # Start from x_signal positions and enforce min_sep between every
    # adjacent pair.  Bidirectional passes (rightward then leftward)
    # distribute pressure symmetrically so labels spread around their
    # signals rather than cascading in one direction.
    x_relaxed      = x_signals.copy()
    max_relax_iter = 100

    for _ in range(max_relax_iter):
        moved = False

        # Rightward pass
        for i in range(1, n):
            gap = x_relaxed[i] - x_relaxed[i - 1]
            if gap < min_sep:
                x_relaxed[i] = x_relaxed[i - 1] + min_sep
                moved = True

        # Leftward pass
        for i in range(n - 2, -1, -1):
            gap = x_relaxed[i + 1] - x_relaxed[i]
            if gap < min_sep:
                x_relaxed[i] = x_relaxed[i + 1] - min_sep
                moved = True

        if not moved:
            break

    # ------------------------------------------------------------------
    # 2. Rail assignment from relaxation drift
    # ------------------------------------------------------------------
    # How far each label drifted from its signal during relaxation is a
    # direct measure of local density: labels in dense regions drift more
    # and should be stacked higher.  We bin the drift into rails using
    # a stride of (rail_width / max_rails) so rails fill proportionally.
    max_rails  = 10
    rail_stride = rail_width / max_rails
    drift       = np.abs(x_relaxed - x_signals)
    rail_ids    = np.clip(
        (drift / rail_stride).astype(int),
        0,
        max_rails - 1,
    )

    # ------------------------------------------------------------------
    # 3. linspace rank-reassignment
    # ------------------------------------------------------------------
    # Sort by x_signal and assign evenly-spaced x_text slots across the
    # full [rail_start, rail_end] range.  This guarantees:
    #   - x_text rank == x_signal rank → no arrow crossings
    #   - Full, even rail coverage regardless of rail_frac
    sig_order = np.argsort(x_signals)
    x_texts   = np.empty(n)
    slots     = (
        np.linspace(rail_start, rail_end, n)
        if n > 1
        else np.array([(rail_start + rail_end) / 2])
    )
    x_texts[sig_order] = slots

    # ------------------------------------------------------------------
    # Build layout table
    # ------------------------------------------------------------------
    layout = pd.DataFrame({
        "label"    : labels,
        "x_signal" : x_signals,
        "x_text"   : x_texts,
        "rail_id"  : rail_ids,
    }).sort_values("x_signal").reset_index(drop=True)

    # ------------------------------------------------------------------
    # 4. Rendering pass — rail_id first used here
    # ------------------------------------------------------------------
    jitter = np.linspace(
        -spread_width * 0.03,
        spread_width * 0.03,
        len(layout),
    )

    for i, row in enumerate(layout.itertuples(index=False)):
        x_sig = row.x_signal
        x_txt = row.x_text
        r_idx = int(row.rail_id)
        label = row.label
        x_tip = x_sig + jitter[i]
        y_txt = y_text_base + r_idx * y_stack_step

        dx  = x_txt - x_sig
        rad = np.clip(
            dx / (genome_width * 0.15),
            -max_rad,
            max_rad,
        )

        if r_idx == 0:
            anglea = 0
            arma   = 0
            armb   = 30
        else:
            anglea = -90
            arma   = 90 * r_idx
            armb   = 30

        arrow = FancyArrowPatch(
            (x_txt, y_txt),
            (x_tip, y_tip),
            arrowstyle="-|>",
            mutation_scale=6,
            lw=0.4,
            color="grey",
            alpha=0.5,
            connectionstyle=(
                f"arc,"
                f"angleA={anglea},"
                f"armA={arma},"
                f"angleB=90,"
                f"armB={armb},"
                f"rad={rad}"
            ),
            transform=ax.transData,
        )
        ax.add_patch(arrow)

        ax.text(
            x_txt + r_idx * spread_width * 0.12,
            y_txt + 0.001,
            str(label),
            rotation=90,
            ha="center",
            va="bottom",
            fontsize=asize,
            clip_on=False,
            color="black",
            fontstyle="italic",
            fontweight="regular",
        )

    # ------------------------------------------------------------------
    # Adaptive ylim
    # ------------------------------------------------------------------
    max_rail = int(layout["rail_id"].max())
    ax.set_ylim(
        y_tip - 0.05,
        y_text_base + (max_rail + 2) * y_stack_step,
    )


def _draw_annotation_arrows_multirail1(
    ax,
    annot_df,
    chr_col: str,
    label_col: str,
    offsets: dict,
    chr_max: dict,
    spread_width: float = 60e6,
    y_text_base: float = 0.25,
    y_stack_step: float = 0.1,
    max_rad: float = 0.35,
    y_tip: float = 0.0,
    fsize: float = 8,
    rail_frac: float = 0.95,
    min_sep: float = 6e6,
    asize: float = 8,
) -> None:
    """
    Dense annotation renderer with relaxation-driven multi-rail
    stacking, linspace rank-reassignment, alternating rail stagger,
    curved arrows, and adaptive ylim.

    Layout pipeline
    ---------------
    1. **Relaxation pass**:
       All labels are sorted by ``x_signal`` and a bidirectional
       relaxation loop enforces ``min_sep`` between every adjacent pair.
       Labels are pushed apart until no two are closer than ``min_sep``.
       The relaxed positions are stored as ``x_relaxed``.

    2. **Rail assignment from relaxation drift**:
       Each label's rail is determined by how far its relaxed position
       drifted from its signal::

           drift    = |x_relaxed − x_signal|
           rail_id  = clip(floor(drift / rail_stride), 0, max_rails − 1)

       Labels in dense regions drift more and receive higher rail
       indices proportionally.  ``rail_stride = rail_width / max_rails``
       so rail assignment scales correctly with ``rail_frac``.

    3. **linspace rank-reassignment**:
       Labels are sorted by ``x_signal`` and assigned evenly-spaced
       ``x_text`` slots via ``np.linspace(rail_start, rail_end, n)``.
       This guarantees:

       - ``x_text`` rank == ``x_signal`` rank → no arrow crossings by
         construction.
       - Full rail coverage regardless of ``rail_frac`` or signal
         clustering.
       - Even slot spacing = ``rail_width / (n − 1)``.

    4. **Cyclic rail stagger (k_min-computed)**:
       The minimum number of rails required to avoid horizontal slot
       overlap is computed from the actual slot interval and
       ``char_width``, with a 1.5× safety factor to account for
       underestimation of rendered glyph widths by the ``0.6 * fsize``
       approximation::

           slot_interval  = rail_width / (n − 1)
           k_min          = ceil(1.5 * char_width / slot_interval)

       Labels are then assigned a stagger offset of ``rank % k_min``
       (cycling through 0, 1, …, k_min−1 in left-to-right
       ``x_signal`` order) which is added to their drift-assigned
       ``rail_id``.  This guarantees same-rail neighbours are at least
       ``k_min * slot_interval >= 1.5 * char_width`` apart.  At lower
       ``rail_frac`` or larger ``fsize``, ``slot_interval`` shrinks and
       ``k_min`` grows automatically, creating more rails as needed.

    5. **Rendering pass**:
       ``rail_id`` is read here for the first time to compute
       ``y = y_text_base + rail_id * y_stack_step``.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    annot_df : pd.DataFrame
        Annotation table. Must contain ``chr_col``, ``"x"`` (cumulative
        genomic position in bp), and ``label_col``.
    chr_col : str
        Column name for chromosome identifiers.
    label_col : str
        Column name for annotation labels.
    offsets : dict
        Mapping of chromosome name → cumulative start offset (bp).
    chr_max : dict
        Mapping of chromosome name → chromosome length (bp).
    spread_width : float, optional
        Genomic window (bp) used for arrow tip jitter. Default 60e6.
    y_text_base : float, optional
        Axes-fraction y-coordinate for rail 0 labels. Default 0.25.
    y_stack_step : float, optional
        Axes-fraction increment per rail. Default 0.1.
    max_rad : float, optional
        Maximum arc curvature for ``FancyArrowPatch``. Default 0.35.
    y_tip : float, optional
        Axes-fraction y-coordinate for arrow tips. Default 0.0.
    fsize : float, optional
        Font size (pt) used to estimate label widths. Default 8.
    rail_frac : float, optional
        Fraction of genome width occupied by the label rail. Default 0.95.
    min_sep : float, optional
        Minimum genomic separation (bp) between any two adjacent label
        centres. Default 6e6.
    asize : float, optional
        Font size (pt) for rendered label text. Default 8.
    """

    # ------------------------------------------------------------------
    # Deduplication
    # ------------------------------------------------------------------
    annot_df = annot_df.drop_duplicates(subset=[chr_col, "x", label_col])

    # ------------------------------------------------------------------
    # Sort annotations
    # ------------------------------------------------------------------
    annot_df = (
        annot_df
        .sort_values(by=[chr_col, "x"], key=natsort_keygen())
        .reset_index(drop=True)
    )

    x_signals = annot_df["x"].to_numpy(dtype=float)
    labels    = annot_df[label_col].astype(str).to_numpy()
    n         = len(x_signals)

    if n == 0:
        return

    # ------------------------------------------------------------------
    # Genome span and rail bounds
    # ------------------------------------------------------------------
    genome_start = min(offsets.values())
    genome_end   = max(offsets[c] + chr_max[c] for c in chr_max)
    genome_width = genome_end - genome_start
    rail_width   = genome_width * rail_frac
    rail_start   = genome_start + (genome_width - rail_width) / 2
    rail_end     = rail_start + rail_width

    # ------------------------------------------------------------------
    # Auto char_width from axes geometry
    # ------------------------------------------------------------------
    # For vertical text (rotation=90°) the horizontal footprint of every
    # label is one character wide regardless of string length.
    try:
        fig         = ax.get_figure()
        renderer    = fig.canvas.get_renderer()
        ax_bbox     = ax.get_window_extent(renderer=renderer)
        xmin, xmax  = ax.get_xlim()
        px_per_bp   = ax_bbox.width / (xmax - xmin)
        char_width  = 0.6 * fsize * (fig.dpi / 72.0) / px_per_bp
    except Exception:
        char_width  = genome_width * 0.01 / fsize

    # ------------------------------------------------------------------
    # 1. Relaxation pass
    # ------------------------------------------------------------------
    # Start from x_signal positions and enforce min_sep between every
    # adjacent pair.  Bidirectional passes (rightward then leftward)
    # distribute pressure symmetrically so labels spread around their
    # signals rather than cascading in one direction.
    x_relaxed      = x_signals.copy()
    max_relax_iter = 50

    for _ in range(max_relax_iter):
        moved = False

        # Rightward pass
        for i in range(1, n):
            gap = x_relaxed[i] - x_relaxed[i - 1]
            if gap < min_sep:
                x_relaxed[i] = x_relaxed[i - 1] + min_sep
                moved = True

        # Leftward pass
        for i in range(n - 2, -1, -1):
            gap = x_relaxed[i + 1] - x_relaxed[i]
            if gap < min_sep:
                x_relaxed[i] = x_relaxed[i + 1] - min_sep
                moved = True

        if not moved:
            break

    # ------------------------------------------------------------------
    # 2. Rail assignment from relaxation drift
    # ------------------------------------------------------------------
    # How far each label drifted from its signal during relaxation is a
    # direct measure of local density: labels in dense regions drift more
    # and should be stacked higher.  We bin the drift into rails using
    # a stride of (rail_width / max_rails) so rails fill proportionally.
    max_rails   = 10
    rail_stride = rail_width / max_rails
    drift       = np.abs(x_relaxed - x_signals)
    rail_ids    = np.clip(
        (drift / rail_stride).astype(int),
        0,
        max_rails - 1,
    )

    # ------------------------------------------------------------------
    # 3. linspace rank-reassignment
    # ------------------------------------------------------------------
    # Sort by x_signal and assign evenly-spaced x_text slots across the
    # full [rail_start, rail_end] range.  This guarantees:
    #   - x_text rank == x_signal rank → no arrow crossings
    #   - Full, even rail coverage regardless of rail_frac
    sig_order = np.argsort(x_signals)
    x_texts   = np.empty(n)
    slots     = (
        np.linspace(rail_start, rail_end, n)
        if n > 1
        else np.array([(rail_start + rail_end) / 2])
    )
    x_texts[sig_order] = slots

    # ------------------------------------------------------------------
    # 4. Compute minimum rails required (k_min) and apply cyclic stagger
    # ------------------------------------------------------------------
    # With n labels evenly spaced across rail_width, the slot interval is:
    #   slot_interval = rail_width / (n - 1)
    #
    # With a cyclic stagger of k rails, same-rail neighbours are k slots
    # apart, giving a same-rail interval of k * slot_interval.
    # The no-overlap condition requires:
    #   k * slot_interval >= char_width
    #   → k_min = ceil(char_width / slot_interval)
    #
    # A safety factor of 1.5 is applied to char_width to account for
    # the fact that 0.6 * fsize underestimates the true rendered glyph
    # width, which varies by font and renderer.  Without this correction
    # k_min is systematically too small at larger font sizes and lower
    # rail_frac values, causing labels to still visually overlap even
    # after staggering.
    #
    # The stagger cycles through k_min values (0, 1, …, k_min-1) in
    # left-to-right x_signal order.  The drift-based rail_id is used as
    # a base elevation and the stagger offset is added on top, so dense
    # regions are still elevated relative to sparse ones while adjacent
    # labels are guaranteed to be on different rails.
    slot_interval      = rail_width / max(n - 1, 1)
    char_width_safe    = char_width * 3
    k_min              = max(1, int(np.ceil(char_width_safe / slot_interval)))

    for rank, idx in enumerate(sig_order):
        stagger         = rank % k_min
        rail_ids[idx]   = min(rail_ids[idx] + stagger, max_rails - 1)

    # ------------------------------------------------------------------
    # Build layout table
    # ------------------------------------------------------------------
    layout = pd.DataFrame({
        "label"    : labels,
        "x_signal" : x_signals,
        "x_text"   : x_texts,
        "rail_id"  : rail_ids,
    }).sort_values("x_signal").reset_index(drop=True)

    # ------------------------------------------------------------------
    # 4. Rendering pass — rail_id first used here
    # ------------------------------------------------------------------
    jitter = np.linspace(
        -spread_width * 0.03,
        spread_width * 0.03,
        len(layout),
    )

    for i, row in enumerate(layout.itertuples(index=False)):
        x_sig = row.x_signal
        x_txt = row.x_text
        r_idx = int(row.rail_id)
        label = row.label
        x_tip = x_sig + jitter[i]
        y_txt = y_text_base + r_idx * y_stack_step

        dx  = x_txt - x_sig
        rad = np.clip(
            dx / (genome_width * 0.15),
            -max_rad,
            max_rad,
        )

        if r_idx == 0:
            anglea = 0
            arma   = 0
            armb   = 30
        else:
            anglea = -90
            arma   = 90 * r_idx
            armb   = 30

        arrow = FancyArrowPatch(
            (x_txt, y_txt),
            (x_tip, y_tip),
            arrowstyle="-|>",
            mutation_scale=6,
            lw=0.4,
            color="grey",
            alpha=0.5,
            connectionstyle=(
                f"arc,"
                f"angleA={anglea},"
                f"armA={arma},"
                f"angleB=90,"
                f"armB={armb},"
                f"rad={rad}"
            ),
            transform=ax.transData,
        )
        ax.add_patch(arrow)

        ax.text(
            x_txt + r_idx * spread_width * 0.12,
            y_txt + 0.001,
            str(label),
            rotation=90,
            ha="center",
            va="bottom",
            fontsize=asize,
            clip_on=False,
            color="black",
            fontstyle="italic",
            fontweight="regular",
        )

    # ------------------------------------------------------------------
    # Adaptive ylim
    # ------------------------------------------------------------------
    max_rail = int(layout["rail_id"].max())
    ax.set_ylim(
        y_tip - 0.05,
        y_text_base + (max_rail + 2) * y_stack_step,
    )



# ---------------------------------------------------------------------------
# Public function
# ---------------------------------------------------------------------------
def plot_linearm(
    tracks: list,
    track_labels: Optional[list[str]] = None,
    annot_df: pd.DataFrame = None,
    annotate: bool = False,
    annotation_size: float = 8,
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
    annot_rail_frac: float = 0.95,
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
    annot_rail_frac : float, optional
        Fraction of horizontal space covering the center of the annotation track within
        which to place annotation texts. Default ``0.98`` (annotation texts will cover
        98% of annotation track horizontally)
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

    annot_df = annot_df.drop_duplicates(subset=[chr_col, pos_col, label_col])

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
        #annot_df = _add_cum(annot_df)
        x_lookup = pd.concat([
            df[[chr_col, pos_col, 'LABEL', "x"]]
            for (label, df) in zip(track_labels, tracks)
        ])

        annot_df = annot_df.merge(
            x_lookup,
            on=[chr_col, pos_col, 'LABEL'],
            how="left"
        )

    # ------------------------------------------------------------------
    # Figure layout
    # ------------------------------------------------------------------
    n_tracks = len(tracks)

    # ------------------------------------------------------------------
    # Track heights sanity check
    # ------------------------------------------------------------------
    # When annotating, track_heights must have n_tracks + 1 values:
    # one for the annotation track (first/top) and one per data track.
    # When not annotating, track_heights must have exactly n_tracks values.
    expected_n = n_tracks + 1 if annotate else n_tracks

    if track_heights is None:
        track_heights = ([1] + [3] * n_tracks) if annotate else ([3] * n_tracks)
    else:
        try:
            if len(track_heights) != expected_n:
                raise ValueError(
                    f"track_heights has {len(track_heights)} values but "
                    f"{expected_n} are required"
                    + (" (one extra for the annotation track)" if annotate else "")
                )
        except TypeError:
            raise TypeError(
                "track_heights must be a sized iterable (e.g. list or tuple), "
                f"got {type(track_heights).__name__}"
            )

    # ------------------------------------------------------------------
    # y-label position
    # ------------------------------------------------------------------
    # Tracks are laid out top-to-bottom, so the annotation track occupies
    # the top portion of the axes and data tracks fill the remainder below.
    # The y-label (-log10(P)) should sit at the vertical midpoint of the
    # data track region in axes-fraction coordinates (0 = bottom, 1 = top).
    #
    # Data tracks span [0, data_total/total_height] from the bottom,
    # so their midpoint in axes-fraction is data_total / (2 * total_height).
    total_height = sum(track_heights)
    data_heights = track_heights[1:] if annotate else track_heights
    data_total   = sum(data_heights)
    y_lab_pos    = data_total / (2 * total_height)


    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(
        expected_n, 1,
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
        s_width = 20e6
        for chr, df in annot_df.groupby(chr_col):
            df_chr = df[df[chr_col]==chr]
            differences = np.diff(df_chr['POS']).tolist()
            less_than_spread_width.append(list(filter(lambda x: x < s_width, differences)))
            less_than_spread_width = [l for l in less_than_spread_width if not len(l) == 0]
        print(len(less_than_spread_width))
        if len(less_than_spread_width) < 5:
            _draw_annotation_arrows(
                ax_annot,
                annot_df,
                chr_col=chr_col,
                label_col=label_col,
                offsets=offsets,
                chr_max=chr_max,
                isolation_threshold=80e6,
                stack_threshold=10e6,
                y_tip=0.0,
                y_text_base=0.3, 
                y_stack_step=0.17,
                spread_width=60e6,
                asize=annotation_size,
            )         
        else:        
            _draw_annotation_arrows_multirail(
                ax_annot,
                annot_df,
                chr_col=chr_col,
                label_col=label_col,
                offsets=offsets,
                chr_max=chr_max,
                spread_width=s_width,
                asize=annotation_size,
                rail_frac=annot_rail_frac,
                y_tip=0.0,
                y_text_base=0.3, 
                y_stack_step=0.17, 
                min_sep=6e6,
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
        0.025, y_lab_pos,
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
    annotation_size: float = 8,
    label_col: Optional[str] = None,
    chr_spacing: Optional[float] = 9e6,
    linear_track_spacing: Optional[float] = None,
    annot_rail_frac: Optional[float] = 0.98,
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
    annot_rail_frac : float, optional
        Fraction of horizontal space covering the center of the annotation track within
        which to place annotation texts. Default ``0.98`` (annotation texts will cover
        98% of annotation track horizontally)        
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

    if track_heights is not None:        
        t_heights = [float(x) for x in track_heights]

    label = 'SNP'
    if annotate:
        label = get_annotation_column(
            annotate=annotate, 
            hits_table=hits_table,
            label_col=label_col
        )

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
        annotation_size=annotation_size,      
        annot_df=hits_table if hits_table is not None and not hits_table.empty else None,
        label_col=label,
        chr_spacing=chr_spacing,
        track_heights=t_heights,
        linear_track_spacing=linear_track_spacing,
        annot_rail_frac=annot_rail_frac,
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
