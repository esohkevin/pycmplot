"""
pycmplot.plotting.circular
==========================
Per-chromosome circular (Circos-style) Manhattan track plotter and
track-radius calculator.
"""

from __future__ import annotations

import logging
import math
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Track radius calculator
# ---------------------------------------------------------------------------

def compute_track_radii_dict(
    n_tracks: int,
    r_min: float = 0,
    r_max: float = 100,
    pad: float = 1,
    annotate: bool = False,
) -> dict[str, tuple[float, float]]:
    """Compute (r_start, r_end) tuples for *n_tracks* evenly-spaced tracks.

    Parameters
    ----------
    n_tracks:
        Number of data tracks.
    r_min, r_max:
        Inner and outer radius of the full plotting area.
    pad:
        Spacing between consecutive tracks.
    annotate:
        If ``True``, add one extra track slot for an annotation ring.

    Returns
    -------
    dict
        ``{"track_1": (start, end), "track_2": (start, end), …}``
    """
    if annotate:
        n_tracks += 1

    total_space = r_max - r_min
    usable_space = total_space - pad * (n_tracks - 1)

    if usable_space <= 0:
        raise ValueError(
            f"Padding ({pad}) is too large for {n_tracks} tracks in "
            f"radius range [{r_min}, {r_max}]."
        )

    track_height = usable_space / n_tracks
    radii: dict[str, tuple[float, float]] = {}
    current = float(r_min)

    for i in range(n_tracks):
        radii[f"track_{i + 1}"] = (current, current + track_height)
        current += track_height + pad

    return radii


# ---------------------------------------------------------------------------
# Per-chromosome circular Manhattan track
# ---------------------------------------------------------------------------

def plot_circosm(
    sector=None,
    sector_radius=None,
    annotation_r=None,
    assoc: Optional[pd.DataFrame] = None,
    sector_sizes: Optional[dict] = None,
    chrom_label_loc=None,
    chrom_label_size: float = 6,
    track_label_size: float = 6,
    track_label_orientation: str = "vertical",
    track_index: int = 0,
    assoc_label: Optional[str] = None,
    logp: bool = True,
    signif_line: Optional[float] = None,
    signif_threshold: Optional[float] = None,
    suggest_line: Optional[float] = None,
    suggest_threshold: Optional[float] = None,
    highlight: bool = False,
    highlight_thresh: Optional[float] = None,
    colors: Optional[list[str]] = None,
) -> None:
    """Plot a single chromosome's data onto a pycirclize sector track.

    This function is called once per (sector, sumstat) pair in the main
    circular Manhattan loop.  It mutates *sector* in-place and returns
    ``None``.  Lead-SNP collection is handled in the calling code.

    Parameters
    ----------
    sector:
        A :class:`pycirclize.Sector` object.
    sector_radius:
        ``(r_start, r_end)`` tuple for this track on the sector.
    assoc:
        Summary statistics DataFrame for **all** chromosomes (filtered to the
        current sector chromosome inside the function).  Must have columns
        ``CHR``, ``POS``, ``P`` (and ``logP`` if *logp* is ``True``).
    sector_sizes:
        Ordered dict of ``{chrom: [min_pos, max_pos]}`` for all sectors,
        used to place labels on the first/last sector.
    track_index:
        0-based index of the current sumstat track (used for chromosome labels).
    colors:
        Two alternating colours for even/odd chromosomes.
    """
    if colors is None:
        colors = ["steelblue", "silver"]

    logger.info("Processing sector: %s", sector.name)

    assoc = assoc.copy()
    assoc["POS"] = assoc["POS"].fillna(0).astype(int)

    genome_wide_sig = signif_threshold
    suggestive = suggest_threshold

    assoc_uniq_chroms = list(assoc["CHR"].unique())

    v_min = float(math.floor(min(assoc["logP"]))) if logp else float(math.floor(min(assoc["P"])))
    v_max = float(math.ceil(max(assoc["logP"]))) if logp else float(math.ceil(max(assoc["P"])))
    if logp:
        v_max += 2
        

    if pd.isna(v_max):
        v_max = 0.0

    sector_keys = list(sector_sizes.keys())

    # ------------------------------------------------------------------
    # Track label on the last (spacer) sector
    # ------------------------------------------------------------------
    if sector.name == sector_keys[-1]:
        lbl_track = sector.add_track(sector_radius)
        lbl_track.axis(fc="white", alpha=0)

        lbl_track.text(
            assoc_label,
            x=(sector.end - sector.start) / 6,
            adjust_rotation=True,
            orientation=track_label_orientation,
            size=float(track_label_size),
            color="black",
            fontstyle="normal",
            fontweight="regular",
            multialignment="left",
        )

    if sector.name not in assoc_uniq_chroms:
        return

    # ------------------------------------------------------------------
    # Chromosome label (first track only, or chrX)
    # ------------------------------------------------------------------
    if track_index == 0 or sector.name == "X":
        sector.text(
            sector.name.replace("23", "X"),
            r=chrom_label_loc,
            size=chrom_label_size,
        )

    sector.axis(fc="none", lw=0, ec="none", alpha=0.5)

    # ------------------------------------------------------------------
    # Y-axis ticks on the first chromosome
    # ------------------------------------------------------------------
    if sector.name == sector_keys[0]:
        yax_track = sector.add_track(sector_radius)
        yax_track.axis(fc="white", alpha=0.08)

        if logp:
            tick_step = 1
            yticks = []
            while len(yticks) < 2 or len(yticks) > 5:
                yticks = np.arange(v_min, v_max, tick_step)
                tick_step += 1
        else:
            yticks = np.arange(v_min, v_max)

        yax_track.yticks(
            yticks,
            labels=[str(int(t)) for t in yticks],
            side="left",
            vmin=v_min,
            vmax=v_max,
            label_size=5,
        )

    # ------------------------------------------------------------------
    # Data track
    # ------------------------------------------------------------------
    assoc_chr = assoc.loc[assoc["CHR"] == sector.name]
    track = sector.add_track(sector_radius, r_pad_ratio=0.05)
    track.axis(fc="lightgrey", alpha=0.08)

    chrom_num = sector.name.replace("X", "23").replace("Y", "24")
    color = colors[0] if int(chrom_num) % 2 == 0 else colors[1]

    y_col = "logP" if logp else "P"

    if highlight:
        sig = assoc_chr[assoc_chr["in_locus"]]
        bg = assoc_chr[~assoc_chr["in_locus"]]

        track.scatter(
            data=bg,
            x=list(bg["POS"].astype(float)),
            y=list(bg[y_col].astype(float)),
            vmin=v_min, vmax=v_max,
            marker="o", s=6, color=color, alpha=1,
        )

        if not sig.empty:
            track.scatter(
                sig["POS"].to_numpy(),
                sig[y_col].to_numpy(),
                vmin=v_min, vmax=v_max,
                s=6, marker="o", color="brown",
            )
    else:
        track.scatter(
            data=assoc_chr,
            x=list(assoc_chr["POS"].astype(float)),
            y=list(assoc_chr[y_col].astype(float)),
            vmin=v_min, vmax=v_max,
            marker="o", s=6, color=color, alpha=1,
        )

    # ------------------------------------------------------------------
    # Significance lines
    # ------------------------------------------------------------------
    if signif_line:
        track.line(
            x=[sector.start, sector.end],
            y=[genome_wide_sig, genome_wide_sig],
            vmin=v_min, vmax=v_max,
            color="orangered", linestyle="--",
        )

    if suggest_line:
        track.line(
            x=[sector.start, sector.end],
            y=[suggestive, suggestive],
            vmin=v_min, vmax=v_max,
            color="lightblue", linestyle="--",
        )
