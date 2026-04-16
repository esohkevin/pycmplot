from __future__ import annotations

CIRCULAR_MODULE = """
pycmplot.plotting.circular
===========================

Circos-style multi-track circular Manhattan plot.

The module exposes two public functions and one internal per-sector helper:

* :func:`plot_circular` — user-facing entry point.  Configures the
  :class:`pycirclize.Circos` canvas, computes track radii, iterates over
  sectors and tracks, renders gene/SNP annotations, and saves the figure.
* :func:`compute_track_radii_dict` — divides the radial space between
  *r_min* and *r_max* into *n_tracks* evenly-spaced, padded bands and
  returns their ``(r_start, r_end)`` limits.
* :func:`plot_circosm` — internal per-sector renderer called once per
  ``(sector, sumstat)`` pair inside the main loop of :func:`plot_circular`.
  Mutates the :class:`pycirclize.Sector` object in place and returns
  ``None``.
"""

import logging
import math
from typing import Optional

import numpy as np
import pandas as pd

from pycmplot.io import get_output_paths
from pycmplot.stats import get_highlight_snps

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Track radius calculator
# ---------------------------------------------------------------------------

def compute_track_radii_dict(
    n_tracks: int,
    r_min: float = 20,
    r_max: float = 100,
    pad: float = 1,
    annotate: bool = False,
) -> dict[str, tuple[float, float]]:
    COMPUTE_RADII = """Compute ``(r_start, r_end)`` tuples for *n_tracks* evenly-spaced radial bands.

    Divides the usable radial space between *r_min* and *r_max* into
    *n_tracks* bands of equal height, separated by gaps of *pad* units.  The
    tracks are ordered from innermost (``'track_1'``) to outermost
    (``'track_n'``).

    Parameters
    ----------
    n_tracks : int
        Number of data tracks to accommodate.
    r_min : float, optional
        Inner boundary of the full plotting area (as a percentage of the
        figure radius).  Default ``20``.
    r_max : float, optional
        Outer boundary of the full plotting area.  Default ``100``.
    pad : float, optional
        Gap in the same radius units between consecutive tracks.  Default ``1``.
    annotate : bool, optional
        If ``True``, an extra slot is reserved for the annotation ring by
        incrementing *n_tracks* before computing heights.  The extra slot
        is always placed at the outermost position.  Default ``False``.

    Returns
    -------
    dict
        Mapping of ``'track_i' → (r_start, r_end)`` for
        ``i`` in ``1 … n_tracks`` (plus one extra entry when *annotate* is
        ``True``).

    Raises
    ------
    ValueError
        If the total padding ``pad × (n_tracks − 1)`` exceeds the available
        radial space ``r_max − r_min``.

    Examples
    --------
    >>> from pycmplot.plotting.circular import compute_track_radii_dict
    >>> radii = compute_track_radii_dict(n_tracks=3, r_min=20, r_max=100, pad=2)
    >>> list(radii.items())
    [('track_1', (20.0, 45.33...)),
    ('track_2', (47.33..., 72.66...)),
    ('track_3', (74.66..., 100.0))]
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
    chrom_label_loc: Optional[float] = -3,
    chrom_label_size: float = 6,
    track_label_size: float = 6,
    track_label_orientation: Optional[str] = "vertical",
    track_index: int = 0,
    assoc_label: Optional[str] = None,
    logp: bool = True,
    signif_line: Optional[float] = 5e-8,
    signif_threshold: Optional[float] = 5e-8,
    suggest_line: Optional[float] = 1e-5,
    suggest_threshold: Optional[float] = 1e-5,
    highlight: bool = False,
    highlight_thresh: Optional[float] = 5e-8,
    highlight_color: str = 'brown',
    colors: Optional[list[str]] = ['steelblue','orange'],
    no_track_labels: bool = False
) -> None:
    PLOT_CIRCOSM = """Plot one track of summary statistics onto a single pycirclize sector.

    This is a low-level internal function called once for every
    ``(sector, sumstat)`` combination in the :func:`plot_circular` main loop.
    It adds a scatter track to *sector* in-place and optionally draws
    significance lines, y-axis ticks (on the first chromosome only), and
    chromosome labels.  Returns ``None``.

    Parameters
    ----------
    sector : pycirclize.Sector
        The pycirclize Sector object representing one chromosome arc.
    sector_radius : tuple of (float, float)
        ``(r_start, r_end)`` radial limits for this track within *sector*,
        as returned by :func:`compute_track_radii_dict`.
    annotation_r : tuple of (float, float) or None
        Radial limits reserved for the annotation ring.  Passed for context
        but not used directly inside this function.
    assoc : pandas.DataFrame, optional
        Full summary statistics DataFrame (all chromosomes).  Filtered to the
        current sector's chromosome internally.  Must have columns ``CHR``,
        ``POS``, ``P``, and ``logP`` (when *logp* is ``True``).
    sector_sizes : dict, optional
        Ordered mapping of ``chrom → [min_pos, max_pos]`` as returned by
        :func:`~pycmplot.io.get_sumstats_and_merged_sector_list`.  Used to
        identify the first and last sectors for y-axis ticks and track labels.
    chrom_label_loc : float or None
        Radial position at which to draw the chromosome label.  Computed in
        :func:`plot_circular` from *chrom_label_side*, *r_min*, and *r_max*.
    chrom_label_size : float, optional
        Font size for chromosome labels.  Default ``6``.
    track_label_size : float, optional
        Font size for the track (sumstat) label written on the spacer sector.
        Default ``6``.
    track_label_orientation : {'vertical', 'horizontal'}, optional
        Orientation of the track label text.  Default ``'vertical'``.
    track_index : int, optional
        0-based index of the current sumstat track.  Chromosome labels are
        only drawn on ``track_index == 0`` (or for chromosome X).
        Default ``0``.
    assoc_label : str, optional
        Track label text (sumstat name) rendered on the spacer sector.
    logp : bool, optional
        If ``True``, use the ``logP`` column for y-values and threshold
        comparisons.  Default ``True``.
    signif_line : float, optional
        Y-value at which to draw the genome-wide significance dashed line
        (orange-red).  Default ``5e-8``.
    signif_threshold : float, optional
        Significance threshold used for y-axis scaling.  Default ``5e-8``.
    suggest_line : float or bool, optional
        Y-value for the suggestive significance dashed line (light blue).
        Pass ``False`` or ``None`` to suppress.  Default ``1e-5``.
    suggest_threshold : float, optional
        Suggestive threshold value used for y-axis scaling.  Default ``1e-5``.
    highlight : bool, optional
        If ``True``, variants within significant loci (``in_locus == True``
        after :func:`~pycmplot.stats.get_highlight_snps`) are rendered in
        ``highlight_color`` (see below).  Default ``False``.
    highlight_thresh : float, optional
        P-value threshold passed to
        :func:`~pycmplot.stats.get_highlight_snps` when *highlight* is
        ``True``.  Default ``5e-8``.
    highlight_color : str, optional
        Color of highlighted positions when *highlight* is ``True``.
        Default ``brown``. 
    colors : list of str, optional
        Two alternating colours for even/odd chromosome numbers.
        Default ``['steelblue', 'orange']``.
    no_track_labels : bool, optional
        Suppress the track label on the spacer sector.  Default ``False``.
    """

    if colors is None:
        colors = ["steelblue", "orange"]

    if highlight:
        assoc, _ = get_highlight_snps(
            df=assoc,
            window=500_000,
            highlight_thresh=highlight_thresh,
            logp=logp,
        )

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

        if no_track_labels:
            pass
        else:
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
        if chrom_label_loc > 100:
            chr_label = str("chr") + str(sector.name.replace("23", "X"))
        else: 
            chr_label = sector.name.replace("23", "X")
        sector.text(
            chr_label,
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
                s=6, marker="o", color=highlight_color,
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


def plot_circular(
    sumstats_loaded: dict,
    sector_sizes: dict = None,
    signif_lines: dict = None,
    logp: bool = False,
    pad: float = 1,
    r_min: float = 20,
    r_max: float = 100,
    annotate: str = None,
    label_col: str = None,
    chrom_label_side: str = 'inside',
    signif_line: float = 5e-8,
    highlight: bool = False,
    highlight_thresh: float = 5e-8,
    highlight_color: str = 'brown',
    highlight_line: bool = False,
    highlight_line_color: str = 'grey',
    colors: list[str] = ['steelblue', 'grey'],
    chrom_label_size: float = 6,
    track_label_size: float = 6,
    track_label_orientation: str = 'vertical',
    hits_table: pd.DataFrame = None,
    annotation_size: float = 6,
    plot_title: Optional[str] = None,
    plot_title_size: float = 12,
    dpi: Optional[int] = None,
    output_format: Optional[str] = None,
    output_dir: Optional[str] = '.',
    no_track_labels: bool = False
):
    PLOT_CIRCULAR = """Generate a multi-track Circos-style circular Manhattan plot.

    Sets up a :class:`pycirclize.Circos` canvas with one arc sector per
    chromosome, computes radial track extents, and calls :func:`plot_circosm`
    once per ``(sector, sumstat)`` pair to populate each track with scatter
    data and significance lines.  After all tracks are rendered, gene or SNP
    annotations from *hits_table* are added to a dedicated annotation ring,
    and a shared y-axis label is placed on the spacer sector.

    Parameters
    ----------
    sumstats_loaded : dict
        Mapping of ``label → [DataFrame, n_chroms]`` as returned by
        :func:`~pycmplot.io.get_sumstats_and_merged_sector_list`.  One
        radial track is created per key.  The outermost track corresponds to
        the first key after reversal of the radii dict.
    sector_sizes : dict, optional
        Ordered mapping of ``chrom → [min_pos, max_pos]`` defining the arc
        length of each chromosome sector.  The last key is expected to be
        ``'Spacer1'`` (automatically added by
        :func:`~pycmplot.io.get_sumstats_and_merged_sector_list`).
    signif_lines : list of dict, optional
        One ``{'genome': float, 'suggestive': float}`` dict per track in the
        same order as *sumstats_loaded*, as returned by
        :func:`~pycmplot.io.get_sumstats_and_merged_sector_list`.
    logp : bool, optional
        Plot –log₁₀(p) radially.  Default ``False``.
    pad : float, optional
        Gap in radius units between consecutive tracks.  Default ``1``.
    r_min : float, optional
        Inner radius of the innermost track (as a percentage of the figure
        radius).  Default ``0``.
    r_max : float, optional
        Outer radius of the outermost track.  Default ``100``.
    annotate : {'SNP', 'GENE'} or falsy, optional
        Annotation content for significant loci.  ``'GENE'`` uses
        ``nearest_upstream_gene`` for genic hits and ``top_gene`` for
        intergenic hits (italic text); ``'SNP'`` uses the ``SNP`` column
        (regular text).  Pass ``None`` or ``False`` to disable annotations.
        Default ``'SNP'``.
    chrom_label_side : {'inside', 'outside'}, optional
        Radial position of chromosome labels.  ``'inside'`` places them just
        inside the innermost track; ``'outside'`` places them beyond the
        outermost track.  Default ``'inside'``.
    signif_line : float, optional
        Genome-wide significance threshold value for the orange-red dashed
        line.  Default ``5e-8``.
    highlight : bool, optional
        Render significant-locus variants in brown.  Default ``False``.
    highlight_thresh : float, optional
        P-value threshold for locus highlighting.  Default ``5e-8``.
    highlight_color : str, optional
        Color of highlighted positions when *highlight* is ``True``.
        Default ``brown``.         
    colors : list of str, optional
        Two alternating chromosome colours.  Default ``['steelblue', 'grey']``.
    chrom_label_size : float, optional
        Chromosome label font size.  Default ``6``.
    track_label_size : float, optional
        Track (sumstat) label font size.  Default ``6``.
    track_label_orientation : {'vertical', 'horizontal'}, optional
        Track label text orientation.  Default ``'vertical'``.
    hits_table : pandas.DataFrame, optional
        Hits summary table from
        :func:`~pycmplot.annotation.get_hits_summary_table`.  Required for
        annotations (``annotate`` truthy and ``hits_table`` non-empty).
    annotation_size : float, optional
        Font size for annotation labels.  Default ``6``.
    highlight_line : bool, optional
        Draw a dashed radial line from the innermost track to the annotation
        ring for each annotated position.  Default ``False``.
    highlight_line_color : str, optional
        Color of highlight line when *highlight_line* is ``True``.
    plot_title : str, optional
        Text placed in the centre of the circle and used as the output
        file-name stem.
    plot_title_size : float, optional
        Font size for the centre title.  Default ``12``.
    dpi : int, optional
        Output resolution in dots per inch.  Default ``300``.
    output_format : str, optional
        Image format (``'png'``, ``'pdf'``, ``'svg'``, ``'jpg'``).
        Default ``'png'``.
    output_dir : str or pathlib.Path, optional
        Output directory.  Default ``'.'``.
    no_track_labels : bool, optional
        Suppress track labels on the spacer sector.  Default ``False``.

    Returns
    -------
    matplotlib.figure.Figure
        The completed circular Manhattan figure (also saved to *output_dir*).

    See Also
    --------
    plot_circular is the circular counterpart to
    :func:`pycmplot.plotting.linear.plot_linear`.

    compute_track_radii_dict :
        Computes the ``(r_start, r_end)`` limits for each track.
    plot_circosm :
        Per-sector rendering function called inside the main loop.
    pycmplot.io.get_sumstats_and_merged_sector_list :
        Produces *sumstats_loaded*, *sector_sizes*, and *signif_lines*.

    Examples
    --------
    >>> from pycmplot.plotting.circular import plot_circular
    >>> fig = plot_circular(
    ...     sumstats_loaded=loaded,
    ...     sector_sizes=sectors,
    ...     signif_lines=sig_lines,
    ...     logp=True,
    ...     highlight=True,
    ...     annotate="GENE",
    ...     hits_table=hits,
    ...     plot_title="RBC_Traits",
    ...     output_dir="./results",
    ... )
    """

    from pycirclize import Circos

    # plot name
    labels = list(sumstats_loaded.keys())
    (
        plt_name, 
        table_out
    ) = get_output_paths(
        labels,
        mode='cm', 
        logp=logp, 
        output_dir=output_dir, 
        plot_title=plot_title, 
        output_format=output_format
    )

    circos = Circos(sector_sizes, space=0.8)

    if plot_title:
        circos.text(text=plot_title, size=plot_title_size, weight="normal")

    n_studies = len(sumstats_loaded)

    radii = compute_track_radii_dict(
        n_tracks=n_studies,
        pad=pad,
        r_min=r_min,
        r_max=r_max,
        annotate=bool(annotate),
    )

    annotation_track_key    = next(reversed(radii))
    annotation_track_radius = radii[annotation_track_key]

    # Reverse so outermost track is plotted first
    radii_reversed = dict(reversed(list(radii.items())))

    inside_loc  = r_min - 3
    outside_loc = 105
    chrom_label_loc = outside_loc if chrom_label_side == "outside" else inside_loc

    if annotate:
        annot_key = next(iter(radii_reversed))
        annot_r   = radii_reversed.pop(annot_key)
        radii_reversed["annot_track_r"] = annot_r

    for index, (sector_radius, sumstats_key, sumstats_value, signif_dict) in enumerate(
        zip(
            radii_reversed.values(),
            sumstats_loaded.keys(),
            sumstats_loaded.values(),
            signif_lines,
        )
    ):
        assoc = sumstats_value[0].copy()
        assoc["P"]   = assoc["P"].dropna()
        assoc["CHR"] = assoc["CHR"].replace("23", "X").replace("24", "Y")
        sumstat_name = sumstats_key

        sig_thresh = signif_dict["genome"]
        sug_thresh = signif_dict["suggestive"]

        logger.info(f"Plotting : {sumstat_name}")
        #logger.info(f"SUGGESTIVE THRESHOLD: {sug_thresh}")

        for sector in circos.sectors:
            plot_circosm(
                sector=sector,
                sector_radius=sector_radius,
                annotation_r=annotation_track_radius if annotate else None,
                sector_sizes=sector_sizes,
                track_index=index,
                chrom_label_loc=chrom_label_loc,
                chrom_label_size=chrom_label_size,
                track_label_size=track_label_size,
                track_label_orientation=track_label_orientation,
                assoc=assoc,
                assoc_label=sumstat_name,
                logp=logp,
                signif_line=sig_thresh,
                signif_threshold=sig_thresh,
                suggest_line=True if signif_line else False,
                suggest_threshold=sug_thresh,
                highlight=highlight,
                highlight_thresh=highlight_thresh,
                highlight_color=highlight_color,
                colors=colors,
                no_track_labels=no_track_labels
            )

    # ------------------------------------------------------------------
    # Circular: gene/SNP annotations
    # ------------------------------------------------------------------
    if annotate and not hits_table.empty:
        for i, (_, row) in enumerate(hits_table.iterrows()):
            label = row['SNP']
            fstyle = "normal" 
            if label_col:
                label_col = str(label_col)
                try:
                    if label_col == "GENE":
                        if row["genic"]:
                            label = row["nearest_upstream_gene"]
                        else:
                            label = row["top_gene"]
                            fstyle = "italic"                       
                    elif label_col != "SNP":
                        label = row[label_col]
                        fstyle = "italic"
                except Exception:
                    logger.info("'SNP' column is used for annotation since '%s' column could not be resolved in hits table.", label_col)
                    pass                
            
            for sector in circos.sectors:
                if str(row["CHR"]) == sector.name:
                    a_track = sector.add_track(annotation_track_radius)
                    a_track.axis(fc="none", lw=0, ec="none", alpha=0)

                    r_low  = annotation_track_radius[0]
                    r_high = annotation_track_radius[1]
                    r_pos  = r_low if i % 2 == 0 else r_high
                    pos    = row["POS"]

                    a_track.annotate(
                        x=pos,
                        label=str(label),
                        min_r=r_low,
                        max_r=r_low + 3,
                        label_size=annotation_size,
                        text_kws={
                            "size": "large",
                            "color": "black",
                            "alpha": 1,
                            "fontstyle": fstyle,
                            "fontweight": "normal",
                            "multialignment": "left",
                        },
                    )

                    if highlight_line:
                        if not highlight_line_color:
                            highlight_line_color = 'grey'
                        sector_rlim = [t.r_lim for t in sector.tracks]
                        sector_min_r = min(sector_rlim)[0]
                        sector.line(
                            r=[sector_min_r, r_low],
                            start=pos,
                            end=pos,
                            color=highlight_line_color,
                            lw=0.5,
                            ls="--",
                        )

    # ------------------------------------------------------------------
    # Circular: single y-axis label on last sector
    # ------------------------------------------------------------------
    for sector in circos.sectors:
        if sector.name == list(sector_sizes.keys())[-1]:
            if logp:
                SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
                y_label = "-log10(p-value)".translate(SUB)
            else:
                y_label = "p-value"

            sector_rlim  = [t.r_lim for t in sector.tracks]
            sector_min_r = min(sector_rlim)[0]
            sector_max_r = max(sector_rlim)[1]

            sector.text(
                y_label,
                x=sector.end - (sector.end - sector.start) / 5,
                r=(sector_min_r + sector_max_r) / 2
                    + (sector_min_r + sector_max_r) / 12,
                adjust_rotation=False,
                ignore_range_error=True,
                size=float(track_label_size),
                color="black",
                fontstyle="italic",
                fontweight="regular",
                rotation=92,
                rotation_mode="default",
                va="top",
                ha="right",
            )

    fig = circos.plotfig()

    if plt_name:
        fig.savefig(fname=plt_name.lower(), dpi=dpi)
        logger.info("Saved circular Manhattan plot: %s", plt_name)

    return fig