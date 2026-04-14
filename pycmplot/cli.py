from __future__ import annotations

CLI_MODULE = '''"""
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
"""'''

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
    GET_ARGUMENTS = '''"""Parse and return command-line arguments for the pycmplot entry point.

    Parameters
    ----------
    descmsg : str, optional
        Banner string printed in the ``--help`` output and echoed to stdout
        at the start of every run.  Defaults to :data:`DESCMSG`.

    Returns
    -------
    argparse.Namespace
        Parsed argument namespace.  Attributes are grouped below.

        **Required**

        .. list-table::
        :widths: 28 12 60
        :header-rows: 1

        * - Attribute
            - Type
            - Description
        * - ``sum_stats``
            - str
            - Comma-separated list of summary statistics file paths.
        * - ``labels``
            - str
            - Comma-separated track labels, same order as ``sum_stats``.
        * - ``build_column``
            - str
            - Column name containing genome-build values (``hg19`` /
            ``hg38``).

        **Optional — input / column resolution**

        .. list-table::
        :widths: 28 12 60
        :header-rows: 1

        * - Attribute
            - Type
            - Description
        * - ``mode``
            - str
            - Plot mode: ``'lm'`` (linear, default) or ``'cm'`` (circular).
        * - ``chrom_column``
            - str or None
            - Chromosome column name; auto-detected when ``None``.
        * - ``pos_column``
            - str or None
            - Base-pair position column name; auto-detected when ``None``.
        * - ``snp_column``
            - str or None
            - Variant / marker ID column name; auto-detected when ``None``.
        * - ``pval_column``
            - str or None
            - P-value column name; auto-detected when ``None``.
        * - ``delim``
            - str or None
            - Delimiter name (``'tab'``, ``'space'``, ``'comma'``,
            ``'colon'``, ``'semi-colon'``); auto-detected when ``None``.

        **Optional — data filtering**

        .. list-table::
        :widths: 28 12 60
        :header-rows: 1

        * - Attribute
            - Type
            - Description
        * - ``logp``
            - bool
            - Plot –log₁₀(p) on the y-axis when ``True``.
        * - ``qq_plot``
            - bool
            - Generate a QQ-plot alongside the Manhattan plot
            *(not yet implemented)*.
        * - ``trim_pval``
            - float or None
            - Drop variants with p > this value before plotting.
            Strongly recommended for large files (e.g. ``0.01``).

        **Optional — significance thresholds**

        .. list-table::
        :widths: 28 12 60
        :header-rows: 1

        * - Attribute
            - Type
            - Description
        * - ``signif_threshold``
            - float or None
            - Genome-wide significance threshold for lead-SNP extraction.
            Defaults to ``5e-8`` when the flag is passed without a value.
        * - ``signif_line``
            - float or None
            - Explicit value for the significance line drawn on the plot.
            Defaults to ``5e-8`` when the flag is passed without a value.
        * - ``suggest_threshold``
            - float or None
            - Suggestive significance threshold for a second dashed line.
            Defaults to ``1e-5`` when the flag is passed without a value.

        **Optional — annotation and highlighting**

        .. list-table::
        :widths: 28 12 60
        :header-rows: 1

        * - Attribute
            - Type
            - Description
        * - ``annotate``
            - str
            - Annotation content: ``'SNP'`` (rsID) or ``'GENE'`` (nearest
            gene symbol).  Default ``'SNP'``.
        * - ``annotation_size``
            - float
            - Font size for annotation labels.  Default ``6``.
        * - ``point_size``
            - float
            - Scatter-plot point size.  Default ``6``.
        * - ``highlight``
            - bool
            - Colour all variants in significant loci distinctly.
        * - ``highlight_thresh``
            - float
            - P-value threshold for locus highlighting.  Default ``5e-8``.
        * - ``highlight_line``
            - bool
            - Draw vertical lines through highlighted locus positions.

        **Optional — appearance and output**

        .. list-table::
        :widths: 28 12 60
        :header-rows: 1

        * - Attribute
            - Type
            - Description
        * - ``colors``
            - str
            - Two comma-separated alternating chromosome colours.
            Default ``'steelblue,silver'``.
        * - ``sort_track``
            - str or None
            - Track sort order: ``'chrom_len'`` or ``'label'``.
        * - ``no_track_labels``
            - bool
            - Suppress track label rendering when ``True``.
        * - ``plot_title``
            - str
            - Plot title and output file stem.  Default ``'MyCMplot'``.
        * - ``plot_title_size``
            - float
            - Plot title font size.  Default ``8``.
        * - ``output_dir``
            - pathlib.Path
            - Output directory.  Default ``Path('.')``.
        * - ``output_format``
            - str
            - Image format: ``'png'``, ``'pdf'``, ``'svg'``, or ``'jpg'``.
            Default ``'png'``.
        * - ``dpi``
            - int
            - Output resolution in dots per inch.  Default ``300``.
        * - ``force``
            - bool
            - Overwrite existing output files when ``True``.

        **Circular-only arguments** (``--mode cm``)

        .. list-table::
        :widths: 28 12 60
        :header-rows: 1

        * - Attribute
            - Type
            - Description
        * - ``pad``
            - int
            - Gap between circular tracks.  Default ``1``.
        * - ``chrom_label_size``
            - float
            - Chromosome label font size.  Default ``6``.
        * - ``chrom_label_side``
            - str
            - Chromosome label placement: ``'inside'`` or ``'outside'``.
            Default ``'inside'``.
        * - ``track_label_size``
            - float
            - Track label font size.  Default ``6``.
        * - ``track_label_orientation``
            - str
            - Track label orientation: ``'vertical'`` or ``'horizontal'``.
            Default ``'vertical'``.
        * - ``r_min``
            - int
            - Inner radius proportion for the innermost track.
            Default ``20``.
        * - ``r_max``
            - int
            - Outer radius proportion for the outermost track.
            Default ``100``.

        **Linear-only arguments** (``--mode lm``)

        .. list-table::
        :widths: 28 12 60
        :header-rows: 1

        * - Attribute
            - Type
            - Description
        * - ``track_heights``
            - str or None
            - Comma-separated relative track heights (e.g. ``'2,2,1.5'``).
        * - ``chr_spacing``
            - float
            - Horizontal gap between chromosomes in base-pairs.
            Default ``9e6``.
        * - ``track_spacing``
            - float
            - Vertical gap between tracks as a fraction of track height.
            Default ``0.10``.

    See Also
    --------
    pycmplot._core.main :
        Consumes the :class:`~argparse.Namespace` returned by this function.
    """'''

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
        "-hc", "--highight_color", default="brown", type=str, metavar="str",
        help="Color of highlighted positions (default: brown)."
    )     
    opt.add_argument(
        "-hll", "--highlight_line", action="store_true",
        help="Draw vertical dashed lines through highlighted positions."
    )     
    opt.add_argument(
        "-hlc", "--highight_line_color", default="grey", type=str, metavar="str",
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
