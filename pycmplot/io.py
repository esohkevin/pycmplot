from __future__ import annotations

MODULE_DOCSTRING = """
pycmplot.io
===========

Functions for loading, validating, and pre-processing GWAS summary statistics
files.  Handles delimiter auto-detection (whitespace, tab, comma), gzip
decompression, and resolution of column-name variants to the canonical set
used throughout the package.

The primary entry point for the plotting pipeline is
:func:`get_sumstats_and_merged_sector_list`, which loads all tracks, runs
coordinate liftover when needed, extracts lead SNPs, generates the hits
summary table, and computes the merged Circos sector-size dictionary — all
in a single call.

Notes
-----
This module is called automatically by the command-line entry point and by
:func:`pycmplot._core.main`; most users will not need to import it directly.
It is documented here for users who wish to load and pre-process summary
statistics programmatically before passing them to the plotting functions.
"""

import csv
import gzip
import sys
import re
import logging
from collections import defaultdict
from pathlib import Path
from typing import Optional

import natsort
import numpy as np
import pandas as pd

from pycmplot.stats import get_lead_snps, get_highlight_snps
from pycmplot.annotation import get_hits_summary_table
from pycmplot.resources import ResourceConfig, default_resources
from pycmplot.constants import hg38_chr_lengths

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# File utilities
# ---------------------------------------------------------------------------

def smart_open(file_path: str):
    SMART_OPEN = """Open a plain-text or gzip-compressed file transparently.

    Detects gzip compression from the ``.gz`` file suffix; all other paths are
    opened as plain text.

    Parameters
    ----------
    file_path : str or pathlib.Path
        Path to the file to open.

    Returns
    -------
    io.TextIOWrapper or gzip.GzipFile
        An open, readable text-mode file object.  Must be used as a context
        manager (``with smart_open(...) as f: ...``).

    Examples
    --------
    >>> from pycmplot.io import smart_open
    >>> with smart_open("HbF.tsv.gz") as f:
    ...     header = f.readline()
    """

    path = Path(file_path)
    if path.suffix == ".gz":
        return gzip.open(file_path, "rt")
    return open(file_path, "r")


def resolve_delimiter(delim: str) -> str:
    RESOLVE_DELIMITER = """Map a human-readable delimiter name to its single-character representation.

    Parameters
    ----------
    delim : str
        A delimiter name — one of ``'space'``, ``'tab'``, ``'comma'``,
        ``'colon'``, ``'semi-colon'``, ``'semicolon'`` — or a single bare
        character (e.g. ``'|'``).  Matching is case-insensitive.

    Returns
    -------
    str
        The corresponding single-character separator string.

    Raises
    ------
    TypeError
        If *delim* is not a string.
    ValueError
        If *delim* is neither a recognised name nor a single character.

    Examples
    --------
    >>> from pycmplot.io import resolve_delimiter
    >>> resolve_delimiter("tab")
    '\\t'
    >>> resolve_delimiter(",")
    ','
    """

    if not isinstance(delim, str):
        raise TypeError("Delimiter must be a string.")

    mapping = {
        "space":      " ",
        "tab":        "\t",
        "comma":      ",",
        "colon":      ":",
        "semi-colon": ";",
        "semicolon":  ";",
    }
    key = delim.strip().lower()
    if key in mapping:
        return mapping[key]
    if len(key) == 1:
        return key  # allow bare characters like '\t'
    raise ValueError(
        f"Invalid delimiter '{delim}'. "
        "Choose from: space, tab, comma, colon, semi-colon."
    )


def detect_delimiter(file_path: str, sample_size: int = 5_000):
    DETECT_DELIMITER = """Infer the field delimiter of a summary statistics file automatically.

    Reads the first *sample_size* bytes of *file_path* and passes the content
    to :class:`csv.Sniffer`.  Falls back to a character-frequency heuristic
    (testing ``','``, ``'\\t'``, ``' '``, ``';'``, ``'|'``) if
    :class:`csv.Sniffer` raises :class:`csv.Error`.

    Parameters
    ----------
    file_path : str or pathlib.Path
        Path to the summary statistics file.  Gzip-compressed files (``.gz``)
        are supported transparently via :func:`smart_open`.
    sample_size : int, optional
        Number of bytes to read for delimiter detection.  Default is ``5000``.

    Returns
    -------
    delimiter : str
        The inferred single-character field separator (e.g. ``'\\t'``,
        ``','``, ``' '``).
    dialect : csv.Dialect or None
        The :class:`csv.Dialect` object returned by :class:`csv.Sniffer`, or
        ``None`` when the fallback heuristic was used.

    Examples
    --------
    >>> from pycmplot.io import detect_delimiter
    >>> delim, dialect = detect_delimiter("HbF.tsv.gz")
    >>> delim
    '\\t'
    """

    with smart_open(file_path) as f:
        sample = f.read(sample_size)

    try:
        dialect = csv.Sniffer().sniff(sample)
        return dialect.delimiter, dialect
    except csv.Error:
        return _fallback_delimiter(sample), None


def _fallback_delimiter(sample: str) -> str:
    candidates = [",", "\t", " ", ";", "|"]
    counts = {d: sample.count(d) for d in candidates}
    best = max(counts, key=counts.get)
    if counts[best] == 0:
        raise ValueError("Unable to detect delimiter automatically.")
    return best


def get_file_header(
    file_path: str,
    delim: Optional[str] = None,
    dialect=None,
) -> list[str]:
    GET_FILE_HEADER = """Read and return the column names from the header line of a file.

    Opens *file_path*, reads the first row using :class:`csv.DictReader`
    configured with the supplied delimiter or dialect, and returns the field
    names as an ordered list of strings.

    Parameters
    ----------
    file_path : str or pathlib.Path
        Path to the summary statistics file (plain text or ``.gz``).
    delim : str, optional
        Field separator character (e.g. ``'\\t'``).  Takes priority over
        *dialect* when both are provided.
    dialect : csv.Dialect, optional
        A :class:`csv.Dialect` object (e.g. as returned by
        :func:`detect_delimiter`).  Used only when *delim* is ``None``.

    Returns
    -------
    list of str
        Ordered list of column names exactly as they appear in the file header.
        Returns an empty list and logs a warning if the header cannot be
        determined.

    Examples
    --------
    >>> from pycmplot.io import detect_delimiter, get_file_header
    >>> delim, dialect = detect_delimiter("HbF.tsv.gz")
    >>> header = get_file_header("HbF.tsv.gz", delim=delim)
    >>> header[:4]
    ['CHR', 'POS', 'SNP', 'P']
    """

    with smart_open(file_path) as f:
        try:
            if delim:
                reader = csv.DictReader(f)
                hdr = f"{delim}".join(reader.fieldnames or []).split(delim)
            elif dialect:
                reader = csv.DictReader(f, dialect=dialect)
                hdr = reader.fieldnames or []
            else:
                reader = csv.DictReader(f)
                hdr = reader.fieldnames or []
        except csv.Error:
            logger.warning("Header could not be determined for %s", file_path)
            hdr = []
    return list(hdr)



def strip_comma_separated_input_streams(
    sum_stats,
    labels,
    colors_raw = 'steelblue,grey',
    track_heights = None,
    builds = None
):
    STRIP_COMMA = """Parse comma-separated CLI strings into Python lists.

    Converts the raw string arguments produced by ``argparse`` (e.g.
    ``"HbF.tsv.gz,MCV.txt.gz,MCH.tsv.gz"``) into the lists expected by the
    rest of the API.  Validates that *sum_stats* and *labels* have the same
    number of elements.

    Parameters
    ----------
    sum_stats : str
        Comma-separated list of summary statistics file paths.
    labels : str
        Comma-separated list of track labels.  Must contain the same number
        of elements as *sum_stats*.
    colors_raw : str, optional
        Comma-separated list of matplotlib colour strings, one per track.
        Default is ``'steelblue,grey'``.
    track_heights : str, optional
        Comma-separated list of relative track heights (floats), one per track.

    Returns
    -------
    sum_stats : list of str
        Parsed file paths, whitespace-stripped.
    labels : list of str
        Parsed track labels, whitespace-stripped.
    colors : list of str
        Parsed colour strings, whitespace-stripped.
    t_heights : list of float
        Parsed track heights converted to ``float``.

    Raises
    ------
    SystemExit
        If *sum_stats* and *labels* have different numbers of elements.
    """

    # ------------------------------------------------------------------
    # Sumstat, labels str to list
    # ------------------------------------------------------------------
    labels     = [lbl.strip() for lbl in labels.strip().split(",")]
    
    sum_stats  = [s.strip() for s in sum_stats.strip().split(",")]

    if builds:
        builds  = [s.strip() for s in builds.strip().split(",")]
        if len(sum_stats) == len(labels) == len(builds):
            pass
        else:
            sys.exit(
                "Error: number of summary stats files, labels, and builds must match.\n"
                f"  Files:  {sum_stats}\n"
                f"  Labels: {labels}"
                f"  Builds: {builds}"
            )

    if len(sum_stats) != len(labels):
        sys.exit(
            "Error: number of summary stats files and labels must match.\n"
            f"  Files:  {sum_stats}\n"
            f"  Labels: {labels}"
        )

    # ------------------------------------------------------------------
    # Colours str to list
    # ------------------------------------------------------------------
    colors = [c.strip() for c in colors_raw.strip().split(",")]

    # ------------------------------------------------------------------
    # Linear track heights str to list
    # ------------------------------------------------------------------
    if track_heights:
        t_heights = [float(x) for x in track_heights.strip().split(",")]
    else:
        t_heights = None

    return sum_stats, labels, colors, t_heights, builds


# ------------------------------------------------------------------
# Random string for output paths
# ------------------------------------------------------------------
def generate_random_string(length):
    GENERATE_RANDOM = """Generate a random alphanumeric string.

    Used internally to create a unique output file-name component when no
    ``--plot_title`` is provided.

    Parameters
    ----------
    length : int
        Number of characters in the output string.

    Returns
    -------
    str
        Random string drawn from ASCII letters (upper- and lower-case) and
        digits (``[A-Za-z0-9]``).

    Examples
    --------
    >>> from pycmplot.io import generate_random_string
    >>> s = generate_random_string(10)
    >>> len(s)
    10
    """

    import random
    import string
    # Combine uppercase, lowercase, and digits
    characters = string.ascii_letters + string.digits
    # random.choices picks multiple characters with replacement
    return ''.join(random.choices(characters, k=length))


# ------------------------------------------------------------------
# Output paths
# ------------------------------------------------------------------
def get_output_paths(
    labels,
    mode: Optional[str] = 'lm',
    logp: bool = False,
    output_dir: Optional[str] = '.',
    plot_title: Optional[str] = None,
    output_format: Optional[str] = 'png'
):
    GET_OUTPUT_PATHS = """Construct output file paths for the plot image and locus summary table.

    Creates *output_dir* (including any missing parent directories) and derives
    deterministic, human-readable file names from the plot title, track labels,
    plot mode, and y-axis scale.

    Parameters
    ----------
    labels : list of str
        Track labels joined with underscores in the output file name.
    mode : {'lm', 'cm'}, optional
        Plot mode: ``'lm'`` for linear Manhattan, ``'cm'`` for circular.
        Default is ``'lm'``.
    logp : bool, optional
        When ``True`` the string ``'_logp'`` is appended to the base name;
        otherwise ``'_pval'`` is appended.  Default is ``False``.
    output_dir : str or pathlib.Path, optional
        Directory in which output files will be written.  Created with
        ``mkdir(parents=True, exist_ok=True)`` if it does not already exist.
        Default is ``'.'``.
    plot_title : str, optional
        Human-readable plot title.  Non-alphanumeric characters are stripped and
        spaces replaced with underscores for safe use in file names.  When
        ``None`` a 10-character random alphanumeric string is used instead.
    output_format : str, optional
        Image file extension without the leading dot (e.g. ``'png'``, ``'pdf'``,
        ``'svg'``).  Default is ``'png'``.

    Returns
    -------
    plt_name : str
        Absolute path to the output plot image file.
    table_out : str
        Absolute path to the output locus summary table TSV file.

    Examples
    --------
    >>> from pycmplot.io import get_output_paths
    >>> plt_name, table_out = get_output_paths(
    ...     labels=["HbF", "MCV"],
    ...     mode="lm",
    ...     logp=True,
    ...     output_dir="./results",
    ...     plot_title="RBC Traits",
    ... )
    >>> plt_name
    '.../results/RBC_Traits_HbF_MCV_lm_logp.png'
    """

    out_path = Path(output_dir).resolve()

    out_path.mkdir(parents=True, exist_ok=True)

    if plot_title:
        pltitle = re.sub(r"[^a-zA-Z0-9\s]", "", plot_title).replace(" ", "_")
    else:
        pltitle = generate_random_string(10)

    labels = [re.sub(r"[^a-zA-Z0-9\s]", "", x).replace(" ", "_") for x in labels]

    suffix     = "_logp" if logp else "_pval"

    plt_base = str(out_path / f"{pltitle}_{'_'.join(labels)}_{mode.lower()}{suffix}")

    plt_name   = f"{plt_base}.{output_format.lower()}"
    
    table_out  = f"{plt_base}_locus_summary_table.tsv"


    return plt_name, table_out, plt_base



# ---------------------------------------------------------------------------
# input formatter
# ---------------------------------------------------------------------------
def prep_pycmplot_input_info(
    sum_stats: list[str],
    labels: list[str],
    buildc: Optional[str] = None,
    build: list[str] = None,
    delim: Optional[str] = None,
    chrom: Optional[str] = None,
    pos: Optional[str] = None,
    snp: Optional[str] = None,
    pcol: Optional[str] = None,
):
    PREP_INPUT_INFO = """Resolve column names and delimiters for each summary statistics file.

    Iterates over every file in *sum_stats*, auto-detects (or uses the supplied)
    delimiter, reads the file header, and maps each required column
    (chromosome, position, SNP ID, p-value, genome build) to the first matching
    entry in an ordered candidate-name list.  Returns a per-label mapping that
    tells :func:`get_sumstats_and_merged_sector_list` exactly which columns to
    read and how to rename them.

    Parameters
    ----------
    sum_stats : list of str
        Paths to one or more summary statistics files (gzip supported).
    labels : list of str
        Track labels in the same order as *sum_stats*.
    build : str
        Genome-build column name (candidates: ``'BUILD'``, ``'Genome'``,
        ``'Genome_Build'``, ``'Genome-build'``, …).
        Or list of genome builds supplied via ``--build``.
    delim : str, optional
        Field delimiter shared by all files.  Accepts human-readable names
        (``'tab'``, ``'space'``, ``'comma'``) or single characters.  When
        ``None`` the delimiter is auto-detected independently for each file
        using :func:`detect_delimiter`.
    chrom : str, optional
        Chromosome column name.  When ``None``, the first header field that
        matches any built-in candidate (``'CHR'``, ``'CHROM'``, ``'#CHROM'``,
        ``'chrom'``, ``'chr'``, …) is used.
    pos : str, optional
        Base-pair position column name (candidates: ``'BP'``, ``'POS'``,
        ``'bp'``, ``'pos'``, ``'Basepair'``).
    snp : str, optional
        Variant / marker ID column name (candidates: ``'SNP'``, ``'RSID'``,
        ``'rsID'``, ``'MarkerName'``, ``'MarkerID'``, ``'SNPID'``, ``'ID'``,
        …).
    pcol : str, optional
        P-value column name (candidates: ``'P'``, ``'P-value'``,
        ``'pvalue'``, ``'p_val'``, ``'pval'``, ``'Wald_P'``).

    Returns
    -------
    dict
        Mapping of ``label → [old_columns, col_dtypes, rename_map, sep]``:

        * **old_columns** – list of the five original column names as found
        in the file header.
        * **col_dtypes** – ``{column_name: dtype}`` passed to
        :func:`pandas.read_csv`.
        * **rename_map** – ``{old_name: canonical_name}`` for
        ``CHR``, ``POS``, ``SNP``, ``P``, ``BUILD``.
        * **sep** – the resolved delimiter character for this file.

    Raises
    ------
    SystemExit
        If any required column (chromosome, position, SNP ID, p-value, or
        build) cannot be resolved from the file header.

    See Also
    --------
    get_sumstats_and_merged_sector_list :
        The main loading function that consumes the mapping returned here.
    detect_delimiter :
        Auto-detects the file delimiter when *delim* is ``None``.
    """

    # ------------------------------------------------------------------
    # Resolve delimiter
    # ------------------------------------------------------------------
    if delim:
        sep = resolve_delimiter(delim)
    else:
        sep = None  # autodetect per file

    # ------------------------------------------------------------------
    # Column-name candidate lists for auto-resolution
    # ------------------------------------------------------------------
    chr_candidates = ["CHR", "CHROM", "Chromosome", "#CHROM", "#CHR",
                    "Chrom", "chrom", "chr", "chromosome", "#chr", "#chrom"]
    chr_candidates_l = [x.lower() for x in chr_candidates]
    chr_candidates_u = [x.upper() for x in chr_candidates]
    chr_candidates = [chrom] + chr_candidates + chr_candidates_l + chr_candidates_u
                
    pos_candidates = ["BP", "POS", "bp", "pos", "Basepair"]
    pos_candidates_l = [x.lower() for x in pos_candidates]
    pos_candidates_u = [x.upper() for x in pos_candidates]
    pos_candidates = [pos] + pos_candidates + pos_candidates_l + pos_candidates_u

    snp_candidates = ["SNP", "RSID", "rsID", "MarkerName", "MarkerID",
                    "Predictor", "Marker", "SNPID", "ID"]
    snp_candidates_l = [x.lower() for x in snp_candidates]
    snp_candidates_u = [x.upper() for x in snp_candidates]
    snp_candidates = [snp] + snp_candidates + snp_candidates_l + snp_candidates_u

    pvl_candidates = ["P", "P-value", "Wald_P", "pvalue", "p_val", "pval"]
    pvl_candidates_l = [x.lower() for x in pvl_candidates]
    pvl_candidates_u = [x.upper() for x in pvl_candidates]
    pvl_candidates = [pcol] + pvl_candidates + pvl_candidates_l + pvl_candidates_u

    # Remove None entries
    chr_candidates = [c for c in chr_candidates if c]
    pos_candidates = [c for c in pos_candidates if c]
    snp_candidates = [c for c in snp_candidates if c]
    pvl_candidates = [c for c in pvl_candidates if c]

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
            pcol           = next(c for c in hdr if c in set(pvl_candidates))
        except StopIteration as exc:
            sys.exit(
                f"Error: could not find a required column in {fpath}.\n"
                f"  Header: {hdr}\n"
                f"  Details: {exc}"
            )

        old_cols = [chrom_col, pos_col, snp_col, pcol]
        new_cols = {
            chrom_col: "CHR",
            pos_col:   "POS",
            snp_col:   "SNP",
            pcol:      "P",
        }
        col_dtypes = {
            chrom_col: str,
            pos_col:   int,
            snp_col:   str,
            pcol:      float,
        }

    # CHECK BUILD COLUMN
    if buildc:
        bld_candidates = buildc
    else:
        bld_candidates = ["BUILD", "Genome", "Genome_Build", "Genome-build"]
        bld_candidates_l = [x.lower() for x in bld_candidates]
        bld_candidates_u = [x.upper() for x in bld_candidates]
        bld_candidates = [buildc] + bld_candidates + bld_candidates_l + bld_candidates_u
        bld_candidates = [c for c in bld_candidates if c]

    try:
        bcol = next(c for c in hdr if c in set(bld_candidates))
    except Exception:
        bcol = False


    if not bcol and isinstance(build, list):
        for name, fpath, build in zip(labels, sum_stats, build):
            sumstats_hdr_dic[name] = [old_cols, col_dtypes, new_cols, file_sep, build]

    elif bcol:
        for name, fpath in zip(labels, sum_stats):
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
    else:
        logger.warning("""No build column or builds detected. Summary stats will be
        plotted in their respective coordinate systems.
        If your data are in different coordinate systems, putting them in one plot
        is not advisable, especially if ``--annotate`` is set!""")
        for name, fpath in zip(labels, sum_stats):
            sumstats_hdr_dic[name] = [old_cols, col_dtypes, new_cols, file_sep]

    return sumstats_hdr_dic


# ---------------------------------------------------------------------------
# Sector-size helpers
# ---------------------------------------------------------------------------

def _merge_min_max_lists(dicts: list[dict]) -> dict:
    """Merge per-chromosome [min, max] lists across multiple sumstats."""
    temp: dict = defaultdict(list)
    for d in dicts:
        for key, values in d.items():
            temp[key].extend(values)
    return {k: [min(v), max(v)] for k, v in temp.items()}


# ---------------------------------------------------------------------------
# Main loader
# ---------------------------------------------------------------------------

def get_sumstats_and_merged_sector_list(
    sum_stats: list[str],
    labels: list[str],
    logp: bool = False,
    trim_pval: Optional[float] = None,
    file_info: Optional[dict] = None,
    sort_tracks: Optional[str] = "chrom_len",
    table_out: Optional[str] = None,
    signif_threshold: Optional[float] = None,
    signif_line: Optional[float] = None,
    suggest_threshold: Optional[float] = None,
    resources: Optional[ResourceConfig] = None,
    hg38_chr_lengths = hg38_chr_lengths,
):
    GET_SUMSTATS = """Load summary statistics, run liftover, extract lead SNPs, and compute
    merged Circos sector sizes.

    This is the primary data-loading function for the plotting pipeline.  For
    each track it: reads the file using the column mapping from *file_info*,
    optionally filters by *trim_pval*, normalises chromosome names
    (``'chr'`` prefix stripped; ``'23'``→``'X'``, ``'24'``→``'Y'``,
    ``'M'``/``'MTDNA'``→``'MT'``), lifts over hg19 coordinates when a build
    column is present, and extracts lead SNPs.  After all tracks are loaded it
    builds the hits summary table, derives significance thresholds, optionally
    sorts tracks, and computes the merged sector-size dict consumed by both
    plotters.

    Parameters
    ----------
    sum_stats : list of str
        Paths to summary statistics files (gzip supported).
    labels : list of str
        Track labels in the same order as *sum_stats*.
    logp : bool, optional
        If ``True``, a ``logP`` column (–log₁₀(P)) is added to every loaded
        DataFrame and used for lead-SNP ranking and threshold-line computation.
        Default is ``False``.
    trim_pval : float, optional
        Drop variants with ``P > trim_pval`` before any further processing.
        Strongly recommended for large files (e.g. ``0.01``).  Default is
        ``None`` (no trimming; variants with ``P > 1`` are still removed).
    file_info : dict, optional
        Column-resolution mapping as returned by
        :func:`prep_pycmplot_input_info`.  Must be supplied for data to be
        loaded.
    sort_tracks : {'label', 'chrom_len', None}, optional
        Track ordering after loading.  ``'label'`` sorts alphabetically;
        ``'chrom_len'`` sorts by the number of distinct chromosomes (most
        chromosomes first).  ``None`` preserves input order.
        Default is ``'chrom_len'``.
    table_out : str, optional
        File path at which to write the locus summary table TSV.  Passed
        through to :func:`~pycmplot.annotation.get_hits_summary_table`.
    signif_threshold : float, optional
        Genome-wide significance threshold for lead-SNP extraction and the
        significance line.  When ``None``, computed as
        ``max(0.05 / N, 5e-8)`` where *N* is the variant count in the last
        loaded track; falls back to ``5e-8`` when *trim_pval* is set.
    signif_line : float, optional
        Explicit significance-line value drawn on the plot.  When ``None``,
        *signif_threshold* is used.  If *logp* is ``True`` and the value is
        < 1, it is converted to –log₁₀ scale automatically.
    suggest_threshold : float, optional
        Suggestive significance threshold for a second dashed line.  Defaults
        to ``1e-5``.
    resources : ResourceConfig, optional
        :class:`~pycmplot.resources.ResourceConfig` instance supplying paths to
        the liftover chain file and gene-info reference files.  Falls back to
        :data:`~pycmplot.resources.default_resources`.

    Returns
    -------
    merged_sector_sizes : dict
        Mapping of ``chromosome → [min_pos, max_pos]`` across all tracks, in
        natural chromosome order (``'1'``, ``'2'``, …, ``'X'``, ``'Y'``),
        with a ``'Spacer1'`` entry appended for y-axis labelling.
    sumstats_loaded : dict
        Mapping of ``label → [DataFrame, n_chroms]``.  Each DataFrame
        contains canonical columns ``CHR``, ``POS``, ``SNP``, ``P``,
        ``LABEL``, and optionally ``logP``, ``BUILD``, ``OLD_POS``,
        ``OLD_BUILD`` (when build and liftover were applied).
    hits_table : pandas.DataFrame
        Clumped locus summary table with nearest-gene annotations.  Empty
        DataFrame when no variants pass the significance threshold.
    signif_lines : list of dict
        One ``{'genome': float, 'suggestive': float}`` dict per track, in
        the final sorted order.

    See Also
    --------
    prep_pycmplot_input_info :
        Resolves column names and delimiters; its output is passed as
        *file_info*.
    pycmplot.annotation.get_hits_summary_table :
        Gene annotation and distance-based clumping of the locus table.
    pycmplot.liftover.liftover_position :
        hg19 → hg38 coordinate conversion applied row-wise.

    Examples
    --------
    >>> from pycmplot.io import prep_pycmplot_input_info
    >>> from pycmplot.io import get_sumstats_and_merged_sector_list
    >>> files  = ["HbF.tsv.gz", "MCV.txt.gz"]
    >>> labels = ["HbF", "MCV"]
    >>> file_info = prep_pycmplot_input_info(files, labels)
    >>> sectors, loaded, hits, sig_lines = get_sumstats_and_merged_sector_list(
    ...     sum_stats=files,
    ...     labels=labels,
    ...     logp=True,
    ...     trim_pval=0.01,
    ...     file_info=file_info,
    ...     signif_threshold=5e-8,
    ... )
    >>> list(sectors.keys())[:4]
    ['1', '2', '3', '4']
    >>> hits.shape
    (18, 12)
    """

    if resources is None:
        resources = default_resources

    from pycmplot.liftover import liftover_position

    # Build a label → file path mapping
    sumstats: dict[str, list] = {
        name: [path] for name, path in zip(labels, sum_stats)
    }

    sumstats_loaded: dict[str, list] = {}
    pval_dict: dict[str, np.ndarray | pd.Series] = {}
    all_lead_snps: list[pd.DataFrame] = []

    for label in sumstats.keys() & (file_info or {}).keys():
        sumstat_cols   = file_info[label][0]
        sumstat_dtypes = file_info[label][1]
        sumstat_newcols= file_info[label][2]
        sep            = file_info[label][3]

        build = None
        try:
            build      = file_info[label][4]
        except Exception:
            pass

        logger.info("Loading %s from %s …", label, sumstats[label][0])
        df = pd.read_csv(
            sumstats[label][0],
            sep=sep,
            header=0,
            usecols=sumstat_cols,
            dtype=sumstat_dtypes,
        ).rename(columns=sumstat_newcols)

        # Get dict of p-values for qq-plotting before applying trim_pval
        logger.info("Extracting raw p-values for QQ-plotting ...")
        pval_dict[label] = df["P"].dropna().astype(float).values


        # Add build column if not exist and build supplied
        if build:
            df['BUILD'] = build

        # Trim insignificant variants for faster plotting
        if trim_pval:
            logger.info("Excluding variants with p-value less than %s to speed up Manhattan plotting ...", trim_pval)
            df = df[df["P"].astype(float) <= float(trim_pval)]
        else:
            df = df[df["P"].astype(float) <= 1]

        if logp:
            logger.info("Adding a 'logP' column ...")
            df["logP"] = -np.log10(df["P"])

        df["LABEL"] = label

        # Normalise chromosome names
        logger.info('Normalizing chromosome names {"23": "X", "24": "Y", "M": "MT", "MTDNA": "MT"} ...')
        df["CHR"] = (
            df["CHR"]
            .str.replace("chr", "", regex=False)
            .dropna()
            .str.upper()
            .replace({"23": "X", "24": "Y", "M": "MT", "MTDNA": "MT"})
        )

        # Number of distinct chromosomes (for track sorting)
        n_chroms = len(df["CHR"].unique()) - 1
        sumstats_loaded[label] = [df, n_chroms]

        # Liftover hg19 data if needed
        if "BUILD" in df.columns and "hg19" in df["BUILD"].unique():
            logger.info("Converting hg19 coordinates to hg38 ...")
            sumstats_loaded[label][0] = liftover_position(df, resources=resources)
            liftover = True

        # Lead SNPs
        logger.info("Extracting variants to highlight ...")
        leads = get_lead_snps(
            df=sumstats_loaded[label][0],
            signif_threshold=signif_threshold or 5e-8,
            logp=True,
        )

        all_lead_snps.append(leads)

    # Combine lead SNPs and filter to significance threshold
    all_lead_snps_df = (
        pd.concat(all_lead_snps, ignore_index=True).drop_duplicates()
        if all_lead_snps
        else pd.DataFrame()
    )
    if not all_lead_snps_df.empty and signif_threshold:
        all_lead_snps_df = all_lead_snps_df[
            all_lead_snps_df["P"] <= signif_threshold
        ]

    hits_table = (
        get_hits_summary_table(
            leads_df=all_lead_snps_df,
            table_out=table_out,
            window_kb=2_000,
            resources=resources,
        )
        if not all_lead_snps_df.empty
        else pd.DataFrame()
    )

    # Derive significance/suggestive thresholds
    if not signif_threshold:
        if trim_pval:
            signif_threshold = 5e-8
        elif sumstats_loaded:
            last_label = list(sumstats_loaded)[-1]
            n = len(sumstats_loaded[last_label][0]["P"])
            signif_threshold = max(0.05 / n, 5e-8)
        else:
            signif_threshold = 5e-8

    if not suggest_threshold:
        suggest_threshold = 1e-5

    suggest_line = suggest_threshold
    if logp:
        suggest_line = -np.log10(suggest_threshold)

    if signif_line is None:
        signif_line = signif_threshold
        if logp:
            signif_line = -np.log10(signif_threshold)
    else:
        if logp and signif_line < 1:
            signif_line = -np.log10(signif_line)

    signif_lines = [
        {"genome": signif_line, "suggestive": suggest_line}
        for _ in sumstats
    ]


    # sort dicts by user-supplied order
    sumstats_loaded = {key: sumstats_loaded[key] for key in labels if key in sumstats_loaded}
    pval_dict = {key: pval_dict[key] for key in labels if key in pval_dict}
    

    # or sort by user option
    if sort_tracks is not None:
        if sort_tracks.lower() == "label":
            sumstats_loaded = dict(sorted(sumstats_loaded.items()))
        else:  # chrom_len
            sumstats_loaded = dict(
                sorted(
                    sumstats_loaded.items(),
                    key=lambda item: (item[0], natsort.natsort_keygen()(item[1][1])),
                )
            )
       

    # Compute per-sumstat sector sizes (chrom → [min_pos, max_pos])
    assoc_sector_sizes_list: list[dict] = []
    min_dic_val = None

    for df, _n in sumstats_loaded.values():
        assoc = df[~(df["CHR"].str.len() > 2)].copy()
        assoc["POS"] = assoc["POS"].fillna(0).astype(int)

        assoc_dic: dict[str, list] = {}
        for chrom in assoc["CHR"].unique():
            sub = assoc[assoc["CHR"] == chrom]
            lo_val = max(sub["POS"].min() - 1_000_000, 0)
            hi_val = sub["POS"].max()
            chrom_max = hi_val

            # Ensure sector sizes are within chrom ranges if liftover
            if liftover:
                hg38_chr_lengths = {k.replace("chr",""): v for k, v in hg38_chr_lengths.items()}
                chrom_max = hg38_chr_lengths[chrom]

            hi_val = min(hi_val, chrom_max)
            assoc_dic[str(chrom)] = [lo_val, hi_val]

        min_dic_val = min(assoc_dic.values())
        assoc_sector_sizes_list.append(assoc_dic)
        # Ensure POS values out of the chrom range are excluded
        assoc = assoc[~assoc["POS"] > hi_val]

    merged = _merge_min_max_lists(assoc_sector_sizes_list)
    merged = dict(natsort.natsorted(merged.items(), key=lambda item: item[0]))

    if "23" in merged:
        merged["X"] = merged.pop("23")

    # Add spacer sector for y-axis labelling
    if min_dic_val is not None:
        #if len(labels) <= 5:
        #    merged["Spacer1"] = [x + x / 2 for x in min_dic_val]
        #else:
        merged["Spacer1"] = [x * 2 for x in min_dic_val]

    return merged, sumstats_loaded, hits_table, signif_lines, pval_dict
