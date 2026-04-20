"""
pycmplot.cli
============

Command-line argument definitions for the ``pycmplot`` console script.

All argument parsing is handled by :func:`get_arguments`, which returns
an :class:`argparse.Namespace` consumed by :func:`~pycmplot._core.main`.
Arguments are organised into four groups:

* **Required** — inputs that must always be supplied.
* **Optional** — flags that control data loading, thresholds, annotation,
  colours, and output format (apply to both plot modes).
* **Circular Only** — arguments specific to ``--mode cm``.
* **Linear Only** — arguments specific to ``--mode lm`` (default).
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
    """Parse and return command-line arguments for the pycmplot entry point.

    Parameters
    ----------
    descmsg : str, optional
        Banner string printed in the ``--help`` output and echoed to stdout
        at the start of every run.  Defaults to :data:`DESCMSG`.

    Returns
    -------
    argparse.Namespace
        Parsed argument namespace containing the attributes listed below.

    Notes
    -----
    **Required arguments**

    ``sum_stats`` : str
        Comma-separated list of summary statistics file paths.
    ``labels`` : str
        Comma-separated track labels, same order as ``sum_stats``.

    **Input / column resolution**

    ``mode`` : {'lm', 'cm'}
        Plot mode: ``'lm'`` (linear, default) or ``'cm'`` (circular).
    ``chrom_column``, ``pos_column``, ``snp_column``, ``pval_column`` : str or None
        Column names in the summary statistics files.  Auto-detected when
        ``None``.
    ``delim`` : {'tab', 'space', 'comma', 'colon', 'semi-colon'} or None
        File delimiter name; auto-detected when ``None``.
    ``build_column`` : str or None
        Column name containing per-variant genome-build values
        (``hg19`` / ``hg38``).
    ``build`` : str or None
        Comma-separated list of genome builds per summary statistics file,
        in the same order as ``sum_stats``.  Alternative to ``build_column``.

    **Data filtering**

    ``logp`` : bool
        Plot −log₁₀(p) on the y-axis when ``True``.
    ``trim_pval`` : float or None
        Drop variants with p > this value before plotting.
        Strongly recommended for large files (e.g. ``0.01``).

    **QQ plot**

    ``qq_plot`` : bool
        Generate a QQ-plot alongside the Manhattan plot.
    ``qq_separate`` : bool
        Save one QQ-plot file per sumstat instead of a combined figure.
    ``qq_overlay`` : bool
        Overlay all sumstats on a single QQ axis.
    ``qq_ncols`` : int
        Number of columns in the combined QQ grid.  Default ``3``.
    ``qq_thin`` : bool
        Enable log-uniform p-value thinning for fast QQ plotting.
    ``thin_below`` : float
        P-value threshold below which all points are retained; points
        above are downsampled.  Default ``0.01``.
    ``qq_max_points`` : int
        Maximum points plotted per QQ track after thinning.
        Default ``50_000``.

    **Significance thresholds**

    ``signif_threshold`` : float or None
        Genome-wide significance threshold for lead-SNP extraction.
        Defaults to ``5e-8`` when the flag is passed without a value.
    ``signif_line`` : float or None
        Value of the significance line drawn on the plot.
        Defaults to ``5e-8`` when the flag is passed without a value.
    ``suggest_threshold`` : float or None
        Suggestive significance threshold for a second dashed line.
        Defaults to ``1e-5`` when the flag is passed without a value.

    **Annotation and highlighting**

    ``annotate`` : str
        Annotation column name in the hits table: ``'snp'`` (rsID),
        ``'top_gene'``, ``'nearest_upstream_gene'``,
        ``'nearest_downstream_gene'``, or ``'gene'`` (let the package
        decide one of ``top_gene`` / ``nearest_upstream_gene``).
    ``annotation_size`` : float
        Font size for annotation labels.  Default ``6``.
    ``point_size`` : float
        Scatter-plot point size.  Default ``6``.
    ``highlight`` : bool
        Colour all variants in significant loci distinctly.
    ``highlight_thresh`` : float
        P-value threshold for locus highlighting.  Default ``5e-8``.
    ``highlight_color`` : str
        Colour for highlighted points.  Default ``'brown'``.
    ``highlight_line`` : bool
        Draw vertical lines through highlighted locus positions.
    ``highlight_line_color`` : str
        Colour of highlight lines.  Default ``'grey'``.

    **Appearance and output**

    ``colors`` : str
        Two comma-separated alternating chromosome colours.
        Default ``'steelblue,silver'``.
    ``sort_track`` : {'chrom_len', 'label'} or None
        Track sort order.
    ``no_track_labels`` : bool
        Suppress track label rendering when ``True``.
    ``plot_title`` : str
        Plot title and output file stem.  Default ``'MyCMplot'``.
    ``plot_title_size`` : float
        Plot title font size.  Default ``8``.
    ``output_dir`` : pathlib.Path
        Output directory.  Default ``Path('.')``.
    ``output_format`` : {'png', 'pdf', 'svg', 'jpg'}
        Image format.  Default ``'png'``.
    ``dpi`` : int
        Output resolution in dots per inch.  Default ``300``.
    ``force`` : bool
        Overwrite existing output files when ``True``.

    **Circular-only arguments** (``--mode cm``)

    ``circular_track_spacing`` : int
        Gap between circular tracks.  Default ``1``.
    ``chrom_label_size`` : float
        Chromosome label font size.  Default ``6``.
    ``chrom_label_side`` : {'inside', 'outside'}
        Chromosome label placement.  Default ``'inside'``.
    ``track_label_size`` : float
        Track label font size.  Default ``6``.
    ``track_label_orientation`` : {'vertical', 'horizontal'}
        Track label orientation.  Default ``'vertical'``.
    ``min_radius`` : int
        Inner radius of the innermost track.  Default ``20``.
    ``max_radius`` : int
        Outer radius of the outermost track.  Default ``100``.

    **Linear-only arguments** (``--mode lm``)

    ``track_heights`` : str or None
        Comma-separated relative track heights (e.g. ``'2,2,1.5'``).
    ``chr_spacing`` : float
        Horizontal gap between chromosomes in base-pairs.  Default ``9e6``.
    ``linear_track_spacing`` : float
        Vertical gap between tracks as a fraction of track height.
        Default ``0.10``.

    See Also
    --------
    pycmplot._core.main :
        Consumes the :class:`~argparse.Namespace` returned by this function.
    """

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
        "-bc",   "--build_column", default=None, required=False, type=str, metavar="str",
                     help=("Name of column containing genome build (hg18/hg19/hg38)." 
                         "Or use ``--build`` below to supply genome builds per summary stat file."
                    ))
    opt.add_argument(
        "-b","--build", default=None, required=False, type=str, metavar='str',
        help=
        """Comma-sperated list of genome build of summary stats file(s) listed 
        in the same order as sumstats files. e.g. hg19,hg38,hg38,hg19 means:
        file1.txt.gz --> hg19
        file2.txt.gz --> hg38
        file3.tsv --> hg38 ... etc
        """
    )
    opt.add_argument(
        "--logp", action="store_true",
        help="Plot −log₁₀(p) instead of raw p-values."
    )
    opt.add_argument("-qq", "--qq_plot", action="store_true",
                     help="Generate QQ-plot(s) alongside the Manhattan plot.")
    opt.add_argument("-qq_sep", "--qq_separate", action="store_true",
                     help=(
                         "Save one QQ-plot file per sumstat instead of a "
                         "combined multi-panel figure. Only used when -qq is set."
                     ))
    opt.add_argument("-qq_cols", "--qq_ncols", default=3, type=int, metavar="int",
                     help="Number of columns in the combined QQ-plot grid (default: 3).")
    opt.add_argument("-qq_thin", "--qq_thin", action="store_true", default=False,
                     help=(
                         "Thin null-like p-values before QQ plotting for speed (default: off)."
                         "Include this flag to turn on for speed."
                    ))
    opt.add_argument("-thin_below", "--thin_below", type=float, metavar="float", default=0.01,
                     help=(
                         "P-value threshold below which all points are always kept."
                         "Points above this threshold are downsampled (default: 0.01)."
                     ))
    opt.add_argument("-qq_max_pts", "--qq_max_points", default=50000, type=int, metavar="int",
                     help="Max points to plot per QQ track after thinning (default: 50000).")
    opt.add_argument("-qq_ov", "--qq_overlay", action="store_true",
                     help=(
                         "Plot all sumstats on a single overlaid QQ-plot, "
                         "each coloured by label with lambda in the legend. "
                         "Only used when -qq is set."
                     ))
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

    # CLASS TO HANDLE ANNOTATION VALUES NOT IN CHOICE LIST
    class AllowAll(list):
        def __contains__(self, item):
            return True

    opt.add_argument(
        "-a", "--annotate",
        choices=AllowAll(["snp", "gene", "top_gene", "nearest_upstream_gene", "nearest_downstream_gene"]), nargs="?",
        default=None, type=str, metavar="{snp,gene,top_gene,nearest_upstream_gene,nearest_downstream_gene,...}", const="SNP",
        help="Annotate loci by column name in hits table (defaults to 'snp' if provided and no value set)."
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
        "-hc", "--highlight_color", default="brown", type=str, metavar="str",
        help="Color of highlighted positions (default: brown)."
    )     
    opt.add_argument(
        "-hll", "--highlight_line", action="store_true",
        help="Draw vertical dashed lines through highlighted positions."
    )     
    opt.add_argument(
        "-hlc", "--highlight_line_color", default="grey", type=str, metavar="str",
        help="Color of highlight line (default: grey)."
    )    
    opt.add_argument(
        "-col", "--colors", default="steelblue,silver", type=str, metavar="str",
        help="Two comma-separated alternating chromosome colours (default: steelblue,silver)."
    )
    opt.add_argument(
        "-st", "--sort_track",
        choices=["chrom_len", "label"], nargs="?",
        const="chrom_len", default=None, type=str, 
        help="Sort tracks by chromosome count or label."
    )
    opt.add_argument(
        "-ntl", "--no_track_labels",
        help=(
            "Exclude track labels from plot. (default: False)"
        ),
        action="store_true"
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
        "-pad", "--circular_track_spacing", default=1, type=int, metavar="int",
        help="Space between circular tracks (default: 1)."
    )
    cio.add_argument(
        "-cl_size", "--chrom_label_size",  default=6, type=float, metavar="float",
        help="Chromosome label font size (default: 6)."
    )
    cio.add_argument(
        "-cl_side", "--chrom_label_side", choices=["inside", "outside"],
        nargs="?", default='inside', const="inside", type=str,
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
        "-r_min", "--min_radius", default=20, type=int, metavar="int",
        help="Inner radius proportion (circular mode, default: 20)."
    )
    cio.add_argument(
        "-r_max", "--max_radius", default=100, type=int, metavar="int",
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
        "-t_space", "--linear_track_spacing", default=0.10, type=float, metavar="float",
        help="Space between linear tracks (default: 0.10)."
    )

    opt.add_argument(
        "-h", "--help", action="help",
        help="Show this help message and exit."
    )

    return parser.parse_args()
