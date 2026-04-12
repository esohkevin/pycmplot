CORE_MODULE = '''"""
pycmplot._core
==============

Main entry point that orchestrates CLI argument parsing, data loading, and
plot dispatch.  This module is intentionally thin: it delegates all heavy
work to :mod:`pycmplot.io`, :mod:`pycmplot.plotting.linear`, and
:mod:`pycmplot.plotting.circular`.

All imports are deferred inside :func:`main` so that
``import pycmplot`` remains fast regardless of the size of the dependency
tree.
"""'''

from __future__ import annotations

import logging
import warnings

# Suppress noisy font-manager warnings before any matplotlib import
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
warnings.filterwarnings("ignore")

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def main() -> None:
    MAIN = '''"""Orchestrate the full pycmplot pipeline from the command line.

    This function is registered as the ``pycmplot`` console-script entry point
    in ``pyproject.toml`` / ``setup.cfg``.  It performs the following steps in
    order:

    1. **Parse CLI arguments** via :func:`~pycmplot.cli.get_arguments`.
    2. **Parse comma-separated inputs** (files, labels, colours, track heights)
    into Python lists via
    :func:`~pycmplot.io.strip_comma_separated_input_streams`.
    3. **Construct output paths** (plot image and locus summary table TSV) via
    :func:`~pycmplot.io.get_output_paths`.
    4. **Resolve column names** for every input file via
    :func:`~pycmplot.io.prep_pycmplot_input_info`.
    5. **Load data** — reads summary statistics, normalises chromosome names,
    runs hg19 → hg38 liftover if needed, extracts lead SNPs, generates the
    hits summary table, and computes merged Circos sector sizes via
    :func:`~pycmplot.io.get_sumstats_and_merged_sector_list`.
    6. **Dispatch plotting** — calls
    :func:`~pycmplot.plotting.circular.plot_circular` when ``--mode cm``,
    or :func:`~pycmplot.plotting.linear.plot_linear` otherwise.

    Parameters
    ----------
    None
        All input is taken from ``sys.argv`` via :mod:`argparse`.

    Returns
    -------
    None
        Saves the plot image and locus summary table to the directory
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
    """'''

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
    from pycmplot.resources import ResourceConfig

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
    build_arg        = args.build_column
    labels_raw       = args.labels
    pcol_arg         = args.pval_column
    logp             = args.logp
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
    highlight_line   = args.highlight_line
    colors_raw       = args.colors
    r_min            = args.r_min
    r_max            = args.r_max
    pad              = args.pad
    output_format    = args.output_format
    output_dir       = args.output_dir
    dpi              = args.dpi
    plot_title       = args.plot_title
    plot_title_size  = args.plot_title_size
    track_heights    = args.track_heights
    track_spacing    = args.track_spacing
    no_track_labels  = args.no_track_labels
    chr_spacing      = args.chr_spacing


    # ------------------------------------------------------------------
    # Sumstat, labels, colours, track heights str to list
    # ------------------------------------------------------------------
    (
        sum_stats, 
        labels, 
        colors, 
        t_heights
    ) = strip_comma_separated_input_streams(
        sum_stats = sum_stats_raw,
        labels = labels_raw,
        colors_raw = colors_raw,
        track_heights = track_heights,
    )

    # ------------------------------------------------------------------
    # Output paths
    # ------------------------------------------------------------------
    (
        plt_name, 
        table_out 
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
        build = build_arg
    )

    # ------------------------------------------------------------------
    # ResourceConfig — picks up environment variables automatically
    # ------------------------------------------------------------------
    resources = ResourceConfig()

    # ------------------------------------------------------------------
    # Load data, compute sectors, get hits table
    # ------------------------------------------------------------------
    (
        merged_assoc_sector_sizes,
        sumstats_loaded,
        hits_table,
        signif_lines,
    ) = get_sumstats_and_merged_sector_list(
        sum_stats=sum_stats,
        labels=labels,
        trim_pval=trim_pval,
        logp=logp,
        file_info=sumstats_hdr_dic,
        sort_tracks=sort_track,
        table_out=table_out,
        signif_threshold=signif_threshold,
        signif_line=signif_line,
        suggest_threshold=suggest_threshold,
        resources=resources,
    )

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
            highlight_line = highlight_line,
            colors = colors,
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
        plot_linear(
            sumstats_loaded = sumstats_loaded,
            track_heights = t_heights,
            trim_pval=trim_pval,
            logp=True if logp else False,
            point_size=point_size,
            highlight=highlight,
            highlight_thresh=highlight_thresh,
            annot_df=hits_table if not hits_table.empty else None,
            label_col="top_gene",
            chr_spacing=chr_spacing,
            track_spacing=track_spacing,
            colors=colors,
            signif_lines=signif_lines,
            plot_title=plot_title,
            no_track_labels = no_track_labels,
            dpi=dpi,
            output_format=output_format,
            output_dir=output_dir,
            figsize=(15, 9)
        )


if __name__ == "__main__":
    main()
