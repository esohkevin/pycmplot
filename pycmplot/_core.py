"""
pycmplot._core
==============
Main entry point — orchestrates CLI parsing, data loading, and plotting.
"""

from __future__ import annotations

import logging
import re
import sys
import warnings
from pathlib import Path

import numpy as np

# Suppress noisy font-manager warnings before any matplotlib import
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
warnings.filterwarnings("ignore")

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def main() -> None:
    """CLI entry point — ``pycmplot`` console script."""

    # ------------------------------------------------------------------
    # Deferred imports so ``import pycmplot`` remains fast
    # ------------------------------------------------------------------
    from pycirclize import Circos

    from pycmplot.cli import get_arguments, DESCMSG
    from pycmplot.io import (
        get_sumstats_and_merged_sector_list,
        detect_delimiter,
        resolve_delimiter,
        get_file_header,
    )
    from pycmplot.plotting.linear import multi_track_linear_manhattan
    from pycmplot.plotting.circular import plot_circosm, compute_track_radii_dict
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

    # ------------------------------------------------------------------
    # Resolve delimiter
    # ------------------------------------------------------------------
    if args.delim:
        sep = resolve_delimiter(args.delim)
    else:
        sep = None  # autodetect per file

    # ------------------------------------------------------------------
    # Output paths
    # ------------------------------------------------------------------
    out_path = Path(output_dir).resolve()
    out_path.mkdir(parents=True, exist_ok=True)

    pltitle    = re.sub(r"[^a-zA-Z0-9\s]", "", plot_title).replace(" ", "_")
    labels     = [lbl.strip() for lbl in labels_raw.strip().split(",")]
    sum_stats  = [s.strip() for s in sum_stats_raw.strip().split(",")]

    if len(sum_stats) != len(labels):
        sys.exit(
            "Error: number of summary stats files and labels must match.\n"
            f"  Files:  {sum_stats}\n"
            f"  Labels: {labels}"
        )

    plt_base   = str(out_path / f"{pltitle}_{'_'.join(labels)}_{mode.lower()}")
    suffix     = "_logp" if logp else "_pval"
    plt_name   = f"{plt_base}{suffix}.{output_format.lower()}"
    table_out  = f"{plt_base}{suffix}_locus_summary_table.tsv"

    # ------------------------------------------------------------------
    # Column-name candidate lists for auto-resolution
    # ------------------------------------------------------------------
    chr_candidates = [chrom_arg, "CHR", "CHROM", "Chromosome", "#CHROM", "#CHR",
                      "Chrom", "chrom", "chr", "chromosome", "#chr", "#chrom"]
    pos_candidates = [pos_arg,   "BP", "POS", "bp", "pos", "Basepair"]
    snp_candidates = [snp_arg,   "SNP", "RSID", "rsID", "MarkerName", "MarkerID",
                      "Predictor", "Marker", "SNPID", "ID"]
    pvl_candidates = [pcol_arg,  "P", "P-value", "Wald_P", "pvalue", "p_val", "pval"]
    bld_candidates = [build_arg, "BUILD", "Genome", "Genome_Build", "Genome-build"]

    # Remove None entries
    chr_candidates = [c for c in chr_candidates if c]
    pos_candidates = [c for c in pos_candidates if c]
    snp_candidates = [c for c in snp_candidates if c]
    pvl_candidates = [c for c in pvl_candidates if c]
    bld_candidates = [c for c in bld_candidates if c]

    # ------------------------------------------------------------------
    # Resolve column names per file
    # ------------------------------------------------------------------
    sumstats_hdr_dic: dict = {}

    for name, fpath in zip(labels, sum_stats):
        if sep:
            file_sep, dialect = sep, None
        else:
            file_sep, dialect = detect_delimiter(fpath, sample_size=5_000)

        hdr = get_file_header(fpath, delim=file_sep, dialect=dialect)

        try:
            chrom_col = next(c for c in hdr if c in set(chr_candidates))
            pos_col   = next(c for c in hdr if c in set(pos_candidates))
            snp_col   = next(c for c in hdr if c in set(snp_candidates))
            pcol      = next(c for c in hdr if c in set(pvl_candidates))
            bcol      = next(c for c in hdr if c in set(bld_candidates))
        except StopIteration as exc:
            sys.exit(
                f"Error: could not find a required column in {fpath}.\n"
                f"  Header: {hdr}\n"
                f"  Details: {exc}"
            )

        old_cols = [chrom_col, pos_col, snp_col, pcol, bcol]
        new_cols = {
            chrom_col: "CHR",
            pos_col:   "POS",
            snp_col:   "SNP",
            pcol:      "P",
            bcol:      "BUILD",
        }
        col_dtypes = {
            chrom_col: str,
            pos_col:   object,
            snp_col:   str,
            pcol:      float,
            bcol:      str,
        }

        sumstats_hdr_dic[name] = [old_cols, col_dtypes, new_cols, file_sep]

    # ------------------------------------------------------------------
    # Colours
    # ------------------------------------------------------------------
    colors = [c.strip() for c in colors_raw.strip().split(",")]

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
        highlight=highlight,
        highlight_thresh=highlight_thresh,
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
        circos = Circos(merged_assoc_sector_sizes, space=0.8)

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
        outside_loc = 101
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
            logger.info(f"SIGNIFICANCE THRESHOLD: {sig_thresh}")
            logger.info(f"SUGGESTIVE THRESHOLD: {sug_thresh}")

            if logp:
                assoc["logP"] = assoc["logP"].dropna()

            for sector in circos.sectors:
                plot_circosm(
                    sector=sector,
                    sector_radius=sector_radius,
                    annotation_r=annotation_track_radius if annotate else None,
                    sector_sizes=merged_assoc_sector_sizes,
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
                    colors=colors,
                )

        # ------------------------------------------------------------------
        # Circular: gene/SNP annotations
        # ------------------------------------------------------------------
        if annotate and not hits_table.empty:
            for i, (_, row) in enumerate(hits_table.iterrows()):
                if str(annotate).upper() != "GENE":
                    label  = row["SNP"]
                    fstyle = "normal"
                else:
                    if row["genic"]:
                        label = row["nearest_upstream_gene"]
                    else:
                        label = row["top_gene"]
                    fstyle = "italic"

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
                            sector_rlim = [t.r_lim for t in sector.tracks]
                            sector_min_r = min(sector_rlim)[0]
                            sector.line(
                                r=[sector_min_r, r_low],
                                start=pos,
                                end=pos,
                                color="lightgrey",
                                lw=0.4,
                                ls="--",
                            )

        # ------------------------------------------------------------------
        # Circular: single y-axis label on last sector
        # ------------------------------------------------------------------
        for sector in circos.sectors:
            if sector.name == list(merged_assoc_sector_sizes.keys())[-1]:
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

        circos.plotfig().savefig(fname=plt_name.lower(), dpi=dpi)
        logger.info("Saved circular Manhattan plot: %s", plt_name)

    # ------------------------------------------------------------------
    # LINEAR MANHATTAN
    # ------------------------------------------------------------------
    else:
        logger.info("Generating LINEAR MANHATTAN Plot ...")
        dfs      = [v[0] for v in sumstats_loaded.values()]
        t_labels = list(sumstats_loaded.keys())

        if not track_heights:
            t_heights = None
        else:
            t_heights = [float(x) for x in track_heights.strip().split(",")]

        multi_track_linear_manhattan(
            tracks=dfs,
            track_labels=t_labels,
            chr_col="CHR",
            pos_col="POS",
            value_col="P",
            trim_p=trim_pval,
            logp=True if logp else False,
            point_size=6,
            highlight=highlight,
            annot_df=hits_table if not hits_table.empty else None,
            label_col="top_gene",
            chr_spacing=9e6,
            track_heights=t_heights,
            track_spacing=0.10,
            colors=colors,
            sig_lines=signif_lines,
            title=plt_name,
            dpi=dpi,
            fig_format=output_format,
            figsize=(15, 9),
        )
        logger.info("Saved linear Manhattan plot: %s", plt_name)


if __name__ == "__main__":
    main()
