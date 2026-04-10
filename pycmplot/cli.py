"""
pycmplot.cli
============
Command-line argument definitions.
"""

from __future__ import annotations

import argparse
from pathlib import Path

DESCMSG = """
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        |  PACKAGE FOR CIRCULAR AND LINEAR MANHATTAN PLOTTING  |
        |                    Kevin Esoh, 2026                  |
        |                    kesohku1@jh.edu                   |
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
"""


def get_arguments(descmsg: str = DESCMSG) -> argparse.Namespace:
    """Parse and return command-line arguments."""

    parser = argparse.ArgumentParser(
        prog="pycmplot",
        description=descmsg,
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=False,
    )

    req = parser.add_argument_group("Required")
    opt = parser.add_argument_group("Optional")
    cio = parser.add_argument_group("Circular Only")
    lio = parser.add_argument_group("Linear Only")

    # ------------------------------------------------------------------
    # Required
    # ------------------------------------------------------------------
    req.add_argument(
        "-s", "--sum_stats",
        help="Comma-separated list of GWAS summary stats files (e.g. file1.txt.gz,file2.tsv).",
        required=True, type=str, metavar="str",
    )
    req.add_argument(
        "-l", "--labels",
        help=(
            "Comma-separated track labels, same order as --sum_stats.\n"
            "E.g. HbF,MCV,MCH"
        ),
        required=True, type=str, metavar="str",
    )
    req.add_argument(
        "-b",   "--build_column",  required=True, type=str, metavar="str",
        help="Genome build column name (containing hg18/hg19/hg38)."
    )

    # ------------------------------------------------------------------
    # Optional
    # ------------------------------------------------------------------
    opt.add_argument(
        "-m", "--mode",
        help="Plot mode: lm (linear Manhattan) or cm (circular Manhattan). Default: lm.",
        choices=["lm", "cm"], default="lm", type=str,
    )
    opt.add_argument(
        "-chr", "--chrom_column",  type=str, metavar="str",
        help="Chromosome column name in sumstats (e.g. CHR)."
    )
    opt.add_argument(
        "-pos", "--pos_column",    type=str, metavar="str",
        help="Position column name (e.g. BP)."
    )
    opt.add_argument(
        "-snp", "--snp_column",    type=str, metavar="str",
        help="SNP ID column name (e.g. ID)."
    )
    opt.add_argument(
        "-p",   "--pval_column",   type=str, metavar="str",
        help="P-value column name (e.g. P)."
    )
    opt.add_argument(
        "-d",   "--delim",
        choices=["space", "tab", "comma", "colon", "semi-colon"],
        type=str, metavar="str",
        help="File delimiter (autodetected if omitted)."
    )
    opt.add_argument(
        "--logp", action="store_true",
        help="Plot −log₁₀(p) instead of raw p-values."
    )
    opt.add_argument(
        "-qq", "--qq_plot", action="store_true",
        help="Also generate a QQ-plot."
    )
    opt.add_argument(
        "-tp", "--trim_pval", type=float, metavar="float",
        help="Trim variants with p > this value before plotting."
    )
    opt.add_argument(
        "-sig", "--signif_threshold",
        default=None, const=5e-8, nargs="?", type=float, metavar="float",
        help="Genome-wide significance threshold (default: 5e-8)."
    )
    opt.add_argument(
        "-sigl", "--signif_line",
        default=None, const=5e-8, nargs="?", type=float, metavar="float",
        help="Value for genome-wide significance line if different from `-sig` (default: 5e-8)."
    )
    opt.add_argument(
        "-sug", "--suggest_threshold",
        default=None, const=1e-5, nargs="?", type=float, metavar="float",
        help="Suggestive significance threshold (default: 1e-5)."
    )
    opt.add_argument(
        "-a", "--annotate",
        choices=["SNP", "GENE"], nargs="?",
        default="SNP", const="SNP", type=str, #metavar="str",
        help="Annotate significant loci by SNP ID or nearest gene."
    )
    opt.add_argument(
        "-p_size", "--point_size", default=6, type=float, metavar="float",
        help="Size of each point of scatter plot (default: 6)."
    )
    opt.add_argument(
        "-a_size", "--annotation_size", default=6, type=float, metavar="float",
        help="Annotation label font size (default: 6)."
    )
    opt.add_argument(
        "-hl",  "--highlight", action="store_true",
        help="Highlight significant loci."
    )
    opt.add_argument(
        "-ht", "--highlight_thresh", default=5e-8, type=float, metavar="float",
        help="P-value threshold for highlighting (default: 5e-8)."
    )
    opt.add_argument(
        "-hl_line", "--highlight_line", action="store_true",
        help="Draw vertical lines through highlighted positions."
    )
    opt.add_argument(
        "--colors", default="steelblue,silver", type=str, metavar="str",
        help="Two comma-separated alternating chromosome colours (default: steelblue,silver)."
    )
    opt.add_argument(
        "-st", "--sort_track",
        choices=["chrom_len", "label"], nargs="?",
        const="chrom_len", default=None, type=str, #metavar="str",
        help="Sort tracks by chromosome count or label."
    )
    opt.add_argument(
        "-plt", "--plot_title", default="MyCMplot", type=str, metavar="str",
        help="Plot plot_title / output file stem."
    )
    opt.add_argument(
        "-pts", "--plot_title_size", default=8, type=float, metavar="float",
        help="Plot plot_title font size (default: 8)."
    )
    opt.add_argument(
        "-od", "--output_dir", default=".", type=Path, metavar="path",
        help="Output directory (default: current directory)."
    )
    opt.add_argument(
        "-of", "--output_format",
        choices=["png", "pdf", "svg", "jpg"],
        default="png", type=str, metavar="str",
        help="Output image format (default: png)."
    )
    opt.add_argument(
        "--dpi", default=300, type=int, metavar="int",
        help="Output resolution in DPI (default: 300)."
    )
    opt.add_argument(
        "-f", "--force", action="store_true",
        help="Overwrite existing output files."
    )

    # circular only
    cio.add_argument(
        "--pad", default=1, type=int, metavar="int",
        help="Space between circular tracks (default: 1)."
    )
    cio.add_argument(
        "-cl_size", "--chrom_label_size",  default=6, type=float, metavar="float",
        help="Chromosome label font size (default: 6)."
    )
    cio.add_argument(
        "-cl_side", "--chrom_label_side", choices=["inside", "outside"],
        nargs="?", default="inside", const="inside", type=str,
        help="Chromosome label placement (default: inside)."
    )
    cio.add_argument(
        "-tl_size", "--track_label_size", default=6, type=float, metavar="float",
        help="Track label font size (default: 6)."
    )
    cio.add_argument(
        "-tl_orient", "--track_label_orientation",
        choices=["vertical", "horizontal"], nargs="?",
        default="vertical", const="vertical", type=str,
        help="Track label orientation (default: vertical)."
    )
    cio.add_argument(
        "--r_min", default=20, type=int, metavar="int",
        help="Inner radius proportion (circular mode, default: 20)."
    )
    cio.add_argument(
        "--r_max", default=100, type=int, metavar="int",
        help="Outer radius (circular mode, default: 100)."
    )

    # linear only
    lio.add_argument(
        "-th", "--track_heights", type=str, metavar="str",
        help="Comma-separated relative track heights (e.g. 2,2,1.5)."
    )
    lio.add_argument(
        "-cs","--chr_spacing", default=9e6, type=float, metavar="float",
        help="Spacing between chromosomes. Useful to reduce chromosome overlap (default: 9e6 or 9000000)."
    )
    lio.add_argument(
        "-t_space", "--track_spacing", default=0.10, type=float, metavar="float",
        help="Space between linear tracks (default: 0.10)."
    )

    opt.add_argument(
        "-h", "--help", action="help",
        help="Show this help message and exit."
    )

    return parser.parse_args()
