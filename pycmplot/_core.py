"""
pycmplot._core
==============

Main entry point that orchestrates CLI argument parsing, data loading, and
plot dispatch.  This module is intentionally thin: it delegates all heavy
work to :mod:`pycmplot.io`, :mod:`pycmplot.plotting.linear`,
:mod:`pycmplot.plotting.circular`, and :mod:`pycmplot.plotting.qq`.

All imports are deferred inside :func:`main` so that
``import pycmplot`` remains fast regardless of the size of the dependency
tree.
"""

from __future__ import annotations

import logging
import warnings
import sys

# Suppress noisy font-manager warnings before any matplotlib import
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
warnings.filterwarnings("ignore")

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s", stream=sys.stdout)
logger = logging.getLogger(__name__)


def main() -> None:
    """Orchestrate the full pycmplot pipeline from the command line.

    This function is registered as the ``pycmplot`` console-script entry point
    in ``pyproject.toml`` / ``setup.cfg``.  It performs the following steps in
    order:

    1. **Parse CLI arguments** via :func:`~pycmplot.cli.get_arguments`.
    2. **Parse comma-separated inputs** (files, labels, colours, track heights,
       builds) into Python lists via
       :func:`~pycmplot.io.strip_comma_separated_input_streams`.
    3. **Construct output paths** (plot image and locus summary table TSV) via
       :func:`~pycmplot.io.get_output_paths`.
    4. **Resolve column names** for every input file via
       :func:`~pycmplot.io.prep_pycmplot_input_info`.
    5. **Load data** — reads summary statistics, normalises chromosome names,
       runs hg19 → hg38 liftover if needed, extracts lead SNPs, generates the
       hits summary table, and computes merged Circos sector sizes via
       :func:`~pycmplot.io.get_sumstats_and_merged_sector_list`.
    6. **Dispatch Manhattan plot** — calls
       :func:`~pycmplot.plotting.circular.plot_circular` when ``--mode cm``,
       or :func:`~pycmplot.plotting.linear.plot_linear` otherwise.
    7. **Optional QQ plot** — when ``--qq_plot`` is set, dispatches to one of
       :func:`~pycmplot.plotting.qq.plot_qq_combined` (default),
       :func:`~pycmplot.plotting.qq.plot_qq_separate` (``--qq_separate``), or
       :func:`~pycmplot.plotting.qq.plot_qq_overlay` (``--qq_overlay``).

    Returns
    -------
    None
        Saves the plot image(s) and locus summary table to the directory
        specified by ``--output_dir``.

    Raises
    ------
    SystemExit
        If required arguments are missing, column names cannot be resolved,
        or the number of summary stats files and labels do not match.

    See Also
    --------
    pycmplot.cli.get_arguments :
        Defines and parses all command-line arguments consumed here.
    pycmplot.io.get_sumstats_and_merged_sector_list :
        The primary data-loading function called in step 5.
    pycmplot.plotting.linear.plot_linear :
        Linear Manhattan plotter called for ``--mode lm`` (default).
    pycmplot.plotting.circular.plot_circular :
        Circular Manhattan plotter called for ``--mode cm``.
    """

    # ------------------------------------------------------------------
    # Deferred imports so ``import pycmplot`` remains fast
    # ------------------------------------------------------------------
    from pycmplot.cli import get_arguments, DESCMSG
    from pycmplot.io import (
        get_sumstats_and_merged_sector_list,
        prep_pycmplot_input_info,
        get_output_paths,
        strip_comma_separated_input_streams,
        #detect_delimiter,
        #resolve_delimiter,
        #get_file_header,
    )
    from pycmplot.plotting.linear import plot_linear
    from pycmplot.plotting.circular import plot_circular
    from pycmplot.plotting.qq import plot_qq_combined, plot_qq_separate, plot_qq_overlay
    from pycmplot.resources import ResourceConfig
    from pycmplot.annotation import get_annotation_column

    # ------------------------------------------------------------------
    # Parse CLI
    # ------------------------------------------------------------------
    args = get_arguments(DESCMSG)
    print(DESCMSG)

    mode             = args.mode
    sum_stats_raw    = args.sum_stats
    chrom_arg        = args.chrom_column
    pos_arg          = args.pos_column
    snp_arg          = args.snp_column
    build_arg        = args.build
    buildc_arg       = args.build_column
    labels_raw       = args.labels
    pcol_arg         = args.pval_column
    logp             = args.logp
    qq               = args.qq_plot
    qq_separate      = args.qq_separate
    qq_ncols         = args.qq_ncols
    qq_thin          = args.qq_thin
    thin_below       = args.thin_below
    qq_max_points    = args.qq_max_points
    qq_overlay       = args.qq_overlay    
    chrom_label_size = args.chrom_label_size
    chrom_label_side = args.chrom_label_side
    track_label_size = args.track_label_size
    track_label_orientation = args.track_label_orientation
    sort_track       = args.sort_track
    trim_pval        = args.trim_pval
    signif_threshold = args.signif_threshold
    signif_line      = args.signif_line
    suggest_threshold= args.suggest_threshold
    annotate         = args.annotate
    annotation_size  = args.annotation_size
    point_size       = args.point_size
    highlight        = args.highlight
    highlight_thresh = args.highlight_thresh
    highlight_color   = args.highlight_color
    highlight_line   = args.highlight_line
    highlight_line_color = args.highlight_line_color
    colors_raw       = args.colors
    r_min            = args.min_radius
    r_max            = args.max_radius
    pad              = args.circular_track_spacing
    output_format    = args.output_format
    output_dir       = args.output_dir
    dpi              = args.dpi
    plot_title       = args.plot_title
    plot_title_size  = args.plot_title_size
    track_heights    = args.track_heights
    linear_track_spacing    = args.linear_track_spacing
    no_track_labels  = args.no_track_labels
    ylabel           = args.ylabel
    chr_spacing      = args.chr_spacing
    figure_size      = args.figure_size


    # ------------------------------------------------------------------
    # Sumstat, labels, colours, track heights [build] str to list
    # ------------------------------------------------------------------
    (
        sum_stats, 
        labels, 
        colors, 
        t_heights,
        builds
    ) = strip_comma_separated_input_streams(
        sum_stats = sum_stats_raw,
        labels = labels_raw,
        colors_raw = colors_raw,
        track_heights = track_heights,
        builds = build_arg if build_arg else None,
    )

    # ------------------------------------------------------------------
    # Output paths
    # ------------------------------------------------------------------
    (
        plt_name, 
        table_out,
        plt_base, 
    ) = get_output_paths(
        labels,
        mode = mode,
        logp = logp,
        output_dir = output_dir,
        plot_title = plot_title,
        output_format = output_format
    )

    # ------------------------------------------------------------------
    # Resolve column names
    # ------------------------------------------------------------------
    sumstats_hdr_dic = prep_pycmplot_input_info(
        sum_stats = sum_stats,
        labels = labels,
        delim = args.delim,
        chrom = chrom_arg,
        pos = pos_arg,
        snp = snp_arg,
        pcol = pcol_arg,
        build_column = buildc_arg,
        build_list = builds
    )

    # ------------------------------------------------------------------
    # ResourceConfig — picks up environment variables automatically
    # ------------------------------------------------------------------
    resources = ResourceConfig()

    # ------------------------------------------------------------------
    # Load data, compute sectors, get hits table
    # ------------------------------------------------------------------
    pycmplot_dict = get_sumstats_and_merged_sector_list(
        sum_stats=sum_stats,
        labels=labels,
        trim_pval=trim_pval,
        logp=logp,
        file_info=sumstats_hdr_dic,
        sort_tracks=sort_track,
        table_out=plt_base,
        signif_threshold=signif_threshold,
        signif_line=signif_line,
        suggest_threshold=suggest_threshold,
        highlight=highlight,
        highlight_thresh=highlight_thresh,
        resources=resources,
        # Only materialise the per-track p-value arrays used for QQ
        # plotting when a QQ render was actually requested.  At 10 M
        # variants per file this avoids an ~80 MB copy + dropna for
        # Manhattan-only / circular-only invocations.
        compute_pvals=bool(qq),
        # Default-on, gwaslab-style density-aware sub-sampling for the
        # Manhattan / circular scatter.  Skipped when the user supplied
        # --no_auto_thin, and a no-op when the dataset is already small
        # enough or when --trim_pval already removed the background.
        auto_thin=not getattr(args, "no_auto_thin", False),
        auto_thin_threshold=getattr(args, "auto_thin_threshold", 2.0),
        auto_thin_max_below=getattr(args, "auto_thin_max_below", 200_000),
    )

    merged_assoc_sector_sizes = pycmplot_dict["sectors"]
    sumstats_loaded = pycmplot_dict["dfs"]
    hits_table = pycmplot_dict["annot"]
    signif_lines = pycmplot_dict["lines"]
    pval_dict = pycmplot_dict["pvals"]

    # ------------------------------------------------------------------
    # CIRCULAR MANHATTAN
    # ------------------------------------------------------------------
    if mode.upper() == "CM":
        logger.info("Generating CIRCULAR MANHATTAN Plot ...")
        plot_circular(
            sumstats_loaded = sumstats_loaded,
            logp = logp,
            signif_line = signif_line,
            signif_lines = signif_lines,
            highlight = highlight,
            highlight_thresh = highlight_thresh,
            highlight_color = highlight_color,
            highlight_line = highlight_line,
            highlight_line_color = highlight_line_color,
            colors = colors,
            point_size=point_size,
            chrom_label_side = chrom_label_side,
            chrom_label_size = chrom_label_size,
            track_label_size = track_label_size,
            track_label_orientation = track_label_orientation,
            annotate = annotate,
            annotation_size = annotation_size,
            hits_table = hits_table,
            sector_sizes = merged_assoc_sector_sizes,
            pad = pad,
            r_min = r_min,
            r_max = r_max,
            plot_title = plot_title,
            plot_title_size = plot_title_size,
            no_track_labels = no_track_labels,
            dpi = dpi,
            output_format=output_format,
            output_dir=output_dir
        )

    # ------------------------------------------------------------------
    # LINEAR MANHATTAN
    # ------------------------------------------------------------------
    else:
        logger.info("Generating LINEAR MANHATTAN Plot ...")
        fsize = figure_size.strip(" ").split(",")
        fsize = [float(v) for v in fsize]
        logger.info(f"FIGURE SIZE: {fsize}")
        plot_linear(
            sumstats_loaded=sumstats_loaded,
            track_heights=t_heights,
            trim_pval=trim_pval,
            logp=True if logp else False,
            point_size=point_size,
            highlight=highlight,
            highlight_thresh=highlight_thresh,
            highlight_color=highlight_color,
            highlight_line=highlight_line,
            highlight_line_color=highlight_line_color,
            annotate=annotate,
            hits_table=hits_table if not hits_table.empty else None,
            chr_spacing=chr_spacing,
            linear_track_spacing=linear_track_spacing,
            colors=colors,
            signif_lines=signif_lines,
            plot_title=plot_title,
            no_track_labels=no_track_labels,
            ylabel=ylabel,
            dpi=dpi,
            output_format=output_format,
            output_dir=output_dir,
            figsize=fsize
        )

    # ------------------------------------------------------------------
    # QQ PLOT
    # ------------------------------------------------------------------
    if qq and sumstats_loaded:
        logger.info("Generating QQ Plot(s) ...")
        qq_stem = f"{plt_base}_qq"
 
        if qq_separate:
            plot_qq_separate(
                pval_dict=pval_dict,
                base_name=plot_title,
                thin=qq_thin,
                thin_below=thin_below,
                max_points=qq_max_points,                
                output_path=qq_stem,
                colors=colors,
                signif_threshold=signif_threshold or 5e-8,
                dpi=dpi,
                fig_format=output_format,
            )
        elif qq_overlay:
            plot_qq_overlay(
                pval_dict=pval_dict,
                thin=qq_thin,
                thin_below=thin_below,
                max_points=qq_max_points,                
                colors=colors,
                signif_threshold=signif_threshold or 5e-8,
                dpi=dpi,
                title=plot_title,
                output_path=f"{qq_stem}_overlay",
                fig_format=output_format,
            )
        else:
            plot_qq_combined(
                pval_dict=pval_dict,
                thin=qq_thin,
                thin_below=thin_below,
                max_points=qq_max_points,
                colors=colors,
                ncols=qq_ncols,
                signif_threshold=signif_threshold or 5e-8,
                dpi=dpi,
                title=plot_title,
                output_path=f"{qq_stem}_combined",
                fig_format=output_format,
            )

if __name__ == "__main__":
    main()