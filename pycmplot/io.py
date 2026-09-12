"""
pycmplot.io
===========

Functions for loading, validating, and pre-processing GWAS summary statistics
files.  Handles delimiter auto-detection (whitespace, tab, comma), gzip
decompression, and resolution of column-name variants to the canonical set
used throughout the package.

The primary entry point for the plotting pipeline is
:func:`get_sumstats_and_merged_sector_list`, which loads all tracks, runs
coordinate liftover when needed, applies density-aware auto-thinning,
extracts lead SNPs, generates the hits summary table, and computes the
merged Circos sector-size dictionary — all in a single call.

Density-aware sub-sampling
--------------------------
:func:`auto_thin_for_manhattan` keeps every variant whose "interestingness"
signal is at or above a configurable threshold (default ``-log10(P) >= 2``,
i.e. ``P <= 0.01``) and uniformly down-samples the dense null background
to at most ``max_below`` rows per track.  It also works for non-p-value
selection scans (iHS, XP-EHH, F_ST, Fay & Wu's H, Tajima's D) via a
``logp=False`` mode that compares ``|statistic|`` to the threshold.  On a
10 M-variant scan this typically cuts the plotted point count from 10 M
to ~200 K + a few hundred peaks — visually indistinguishable above the
suggestive band, ~one to two orders of magnitude faster to render.

Per-track Stage-1 cache
-----------------------
When :func:`get_sumstats_and_merged_sector_list` is invoked with
``cache=True``, each track's post-load / post-liftover / post-thinning
DataFrame is written to ``<cache_dir>/tracks/<label>.<short_key>.parquet``
alongside its lead-SNP table and (when ``compute_pvals=True``) its full
raw p-value array.  A JSON metadata file records a per-track ``cache_key``
computed from the raw file's SHA-256, the pycmplot version, and every
Stage-1 parameter that affects the cached data.  Subsequent invocations
skip Stage 1 for tracks whose cache_key matches — a 2-10x speedup that
scales with input size.  See :mod:`pycmplot.cache` for the on-disk layout,
cache-key semantics, the user-editable hits overlay
(``<cache_dir>/annotations/hits.tsv``), and multi-panel safety notes.

The cache is track-content-aware, not label-aware: two loader calls that
share a ``cache_dir`` and reuse a track label (say ``"Hb"`` in panel A
and ``"Hb"`` in panel B, pointing at different sumstats) get separate
cache entries that coexist safely — the ``<short_key>`` suffix in the
per-track filename disambiguates them and metadata is keyed on the full
cache_key rather than the label alone.

Public functions
----------------
* :func:`prep_pycmplot_input_info` — resolve delimiters and column-name
  mappings for each input file; must be called before the loader.
* :func:`get_sumstats_and_merged_sector_list` — the main loader.  Loads
  all tracks (with optional cache/resume), applies liftover / auto-thin
  / trim, extracts leads, builds the hits summary, and computes the
  merged Circos sector-size dict.  Returns a bundle consumed directly by
  :func:`pycmplot.plotting.linear.plot_linear`,
  :func:`pycmplot.plotting.circular.plot_circular`, and the QQ plotters.
* :func:`auto_thin_for_manhattan` — the density-aware sub-sampling
  helper described above; exported for advanced pipelines.
* :func:`get_output_paths` — resolve the deterministic output file
  names used by the plotters when *plot_title* is supplied.

Notes
-----
This module is called automatically by the command-line entry point and by
:func:`pycmplot._core.main`; most users will not need to import it directly.
It is documented here for users who wish to load and pre-process summary
statistics programmatically before passing them to the plotting functions.

See Also
--------
pycmplot.cache : per-track Stage-1 cache implementation.
pycmplot.annotation : lead-SNP annotation with nearest-gene lookup.
pycmplot.liftover : hg18/hg19 -> hg38 coordinate conversion.
"""

from __future__ import annotations

import csv
import gzip
import os
import sys
import re
import logging
from collections import defaultdict
from pathlib import Path
from typing import Optional

import natsort
import numpy as np
import pandas as pd

from pycmplot.constants import CHROM_ORDER
from pycmplot.stats import get_lead_snps, get_highlight_snps
from pycmplot.annotation import get_hits_summary_table
from pycmplot.resources import ResourceConfig, default_resources

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s", stream=sys.stdout)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# File utilities
# ---------------------------------------------------------------------------

def smart_open(file_path: str):
    """Open a plain-text or gzip-compressed file transparently.

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
    """Map a human-readable delimiter name to its single-character representation.

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
    """Infer the field delimiter of a summary statistics file automatically.

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
    """Read and return the column names from the header line of a file.

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
                reader = csv.DictReader(f, delimiter=delim)
                hdr = reader.fieldnames or []
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
    """Parse comma-separated CLI strings into Python lists.

    Converts the raw string arguments produced by ``argparse`` (e.g.
    ``"HbF.tsv.gz,MCV.txt.gz,MCH.tsv.gz"``) into the lists expected by the
    rest of the API.  Validates that *sum_stats*, *labels* and *builds*
    (when supplied) have the same number of elements.

    Parameters
    ----------
    sum_stats : str
        Comma-separated list of summary statistics file paths.
    labels : str
        Comma-separated list of track labels.  Must contain the same number
        of elements as *sum_stats*.
    colors_raw : str, optional
        Comma-separated list of matplotlib colour strings.  Default is
        ``'steelblue,grey'``.
    track_heights : str, optional
        Comma-separated list of relative track heights (floats), one per
        track.
    builds : str, optional
        Comma-separated list of genome builds (e.g.
        ``'hg19,hg38,hg38,hg19'``), one per summary statistics file.

    Returns
    -------
    sum_stats : list of str
        Parsed file paths, whitespace-stripped.
    labels : list of str
        Parsed track labels, whitespace-stripped.
    colors : list of str
        Parsed colour strings, whitespace-stripped.
    t_heights : list of float or None
        Parsed track heights converted to ``float``.  ``None`` when
        *track_heights* was not supplied.
    builds : list of str or None
        Parsed build strings, whitespace-stripped.  ``None`` when
        *builds* was not supplied.

    Raises
    ------
    SystemExit
        If *sum_stats*, *labels* and *builds* have mismatched lengths.
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
    """Generate a random alphanumeric string.

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
    output_dir: Optional[str] = None,
    plot_title: Optional[str] = None,
    output_format: Optional[str] = 'png'
):
    """Construct output file paths for the plot image and locus summary table.

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
    plt_base : str
        Absolute path base (no extension) used to derive the QQ-plot output
        stems.

    Examples
    --------
    >>> from pycmplot.io import get_output_paths
    >>> plt_name, table_out, plt_base = get_output_paths(
    ...     labels=["HbF", "MCV"],
    ...     mode="lm",
    ...     logp=True,
    ...     output_dir="./results",
    ...     plot_title="RBC Traits",
    ... )
    >>> plt_name
    '.../results/RBC_Traits_HbF_MCV_lm_logp.png'
    """

    if output_dir is None:
        output_dir = '.'
    out_path = Path(output_dir).resolve()

    out_path.mkdir(parents=True, exist_ok=True)

    if plot_title:
        #pltitle = re.sub(r"[^a-zA-Z0-9\s]", "", plot_title).replace(" ", "_")
        pltitle = [ re.sub(r"[^a-zA-Z0-9\s]", "", p) for p in plot_title.split('_') ]
        pltitle = '_'.join(pltitle).replace(" ", "_")
    else:
        pltitle = generate_random_string(10)

    #labels = [re.sub(r"[^a-zA-Z0-9\s]", "", x).replace(" ", "_") for x in labels]
    labels = [ x.replace(" ", "_") for x in labels ]

    suffix     = "_logp" if logp else "_pval"

    plt_base = str(out_path / f"{pltitle}_{'_'.join(labels)}_{mode.lower()}{suffix}")

    plt_name   = f"{plt_base}.{output_format.lower()}"
    
    table_out  = f"{plt_base}_locus_summary_table.tsv"


    return plt_name, table_out, plt_base



# ---------------------------------------------------------------------------
# input formatter
# ---------------------------------------------------------------------------
def prep(
    sum_stats: list[str],
    labels: list[str],
    build_column: Optional[str] = None,
    build_list: list[str] = None,
    delim: Optional[str] = None,
    chrom: Optional[str] = None,
    pos: Optional[str] = None,
    snp: Optional[str] = None,
    pcol: Optional[str] = None,
):
    """Resolve column names and delimiters for each summary statistics file.

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
    build_column : str, optional
        Genome-build column name (candidates: ``'BUILD'``, ``'Genome'``,
        ``'Genome_Build'``, ``'Genome-build'``, …).
        Or list of genome builds supplied via ``--build``.
    build_list : list, optional
        List of genome builds in same order as sumstats and labels
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
        Mapping of ``label -> [old_columns, col_dtypes, rename_map, sep]``:

        * **old_columns** -- list of the five original column names as
          found in the file header.
        * **col_dtypes** -- ``{column_name: dtype}`` passed to
          :func:`pandas.read_csv`.
        * **rename_map** -- ``{old_name: canonical_name}`` for ``CHR``,
          ``POS``, ``SNP``, ``P``, ``BUILD``.
        * **sep** -- the resolved delimiter character for this file.

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
                
    pos_candidates = ["BP", "POS", "Basepair", "position"]
    pos_candidates_l = [x.lower() for x in pos_candidates]
    pos_candidates_u = [x.upper() for x in pos_candidates]
    pos_candidates = [pos] + pos_candidates + pos_candidates_l + pos_candidates_u

    snp_candidates = ["SNP", "RSID", "rsID", "MarkerName", "MarkerID",
                    "Predictor", "Marker", "SNPID", "ID"]
    snp_candidates_l = [x.lower() for x in snp_candidates]
    snp_candidates_u = [x.upper() for x in snp_candidates]
    snp_candidates = [snp] + snp_candidates + snp_candidates_l + snp_candidates_u

    pvl_candidates = ["P", "P-value", "Wald_P", "pvalue", "p_val", "pval",
                    "IHS", "RSB", "LOGP", "LOGPVALUE"]
    pvl_candidates_l = [x.lower() for x in pvl_candidates]
    pvl_candidates_u = [x.upper() for x in pvl_candidates]
    pvl_candidates = [pcol] + pvl_candidates + pvl_candidates_l + pvl_candidates_u

    # Remove None entries
    chr_candidates = [c for c in chr_candidates if c]
    pos_candidates = [c for c in pos_candidates if c]
    snp_candidates = [c for c in snp_candidates if c]
    pvl_candidates = [c for c in pvl_candidates if c]

    # ------------------------------------------------------------------
    # Per-file build tokens
    # ------------------------------------------------------------------
    # Each entry of ``build_list`` can now take one of two shapes:
    #
    #   • Literal build name — ``"hg19"``, ``"hg38"``, ``"GRCh37"``, …
    #     (all values recognised by the ``BUILD_MAP`` further down).
    #     Treated identically to the pre-0.4.x behaviour: the file
    #     has no per-row build column, so we stamp every row with the
    #     given build.
    #
    #   • Column-reference token — ``"col:<colname>"`` (case-insensitive
    #     prefix; ``"C:<colname>"`` is also accepted as the shorthand
    #     Kevin proposed).  Says "for this file, source per-row builds
    #     from the column named ``<colname>``".  Useful when *most*
    #     files declare their build via ``--build`` but one file happens
    #     to carry a per-row build column with a non-standard name
    #     (e.g. ``my_build``) that the auto-detector can't find.
    #
    # The parser returns a tuple ``(kind, value)`` where ``kind`` is
    # either ``"literal"`` or ``"column"``.  Anything unparseable falls
    # through as ``("literal", token)`` and is validated later by the
    # BUILD_MAP lookup — that keeps the error message pointing at the
    # actual offending token instead of at this parser.
    def _parse_build_token(token: str) -> tuple[str, Optional[str]]:
        """Classify a ``--build`` token as either a literal or a column ref.

        Returns ``("literal", "hg19")`` for a plain build name and
        ``("column", "my_build")`` for an explicit column reference.

        Ergonomic fallback: a bare marker with no name — ``col:``,
        ``col``, ``column``, ``C:``, ``C``, etc. — resolves to
        ``("column", None)``, meaning "this file has a build column
        somewhere; auto-detect its name using the standard candidate
        list (``BUILD`` / ``Genome`` / …).  We only complain if that
        auto-detection also fails.  This matches what a user
        naturally types when they know a file has *some* build column
        and don't want to type its exact name.
        """
        if token is None:
            return ("literal", None)
        s = str(token).strip()
        # Prefixed form.  Empty name after the prefix is legal and
        # means "auto-detect using standard build-column candidates".
        for prefix in ("col:", "COL:", "column:", "Column:", "C:", "c:"):
            if s.startswith(prefix):
                colname = s[len(prefix):].strip()
                return ("column", colname if colname else None)
        # Bare short-form marker: `c`, `col`, `column` (any case).
        # Same "auto-detect" intent as the prefixed empty form.
        if s.lower() in ("c", "col", "column"):
            return ("column", None)
        return ("literal", s)

    # Pre-parse every token so downstream code inspects a normalised
    # tuple instead of re-running the prefix check per iteration.
    if build_list is None:
        _build_tokens: list[tuple[str, str] | None] = []
    else:
        _build_tokens = [_parse_build_token(b) for b in build_list]

    # ------------------------------------------------------------------
    # Build-column candidate list (shared across all files)
    # ------------------------------------------------------------------
    if build_column:
        # User supplied a specific build-column name: look only for that name
        bld_candidates = [build_column]
    else:
        bld_candidates = ["BUILD", "Genome", "Genome_Build", "Genome-build"]
        bld_candidates_l = [x.lower() for x in bld_candidates]
        bld_candidates_u = [x.upper() for x in bld_candidates]
        bld_candidates = bld_candidates + bld_candidates_l + bld_candidates_u
        bld_candidates = [c for c in bld_candidates if c]

    # ------------------------------------------------------------------
    # Resolve column names per file
    # ------------------------------------------------------------------
    sumstats_hdr_dic: dict = {}
    user_pcol = pcol  # preserve user-supplied p-column hint across iterations

    for idx, (name, fpath) in enumerate(zip(labels, sum_stats)):
        if sep:
            file_sep, dialect = sep, None
        else:
            file_sep, dialect = detect_delimiter(fpath, sample_size=5_000)

        hdr = get_file_header(fpath, delim=file_sep, dialect=dialect)

        # Reuse the p-value candidate list built once above.  The
        # previous code duplicated a hard-coded base list here, which
        # silently discarded any additions made to ``pvl_candidates``
        # at function scope (e.g. ``IHS``, ``RSB``, ``LOGP``).  A
        # per-iteration rebuild isn't necessary — ``user_pcol`` is not
        # mutated inside the loop, so the outer list is stable.
        pvl_cands = pvl_candidates

        # Two-pass column resolution:
        #   Pass 1: user-supplied hint wins.  Case-insensitive match
        #           against the file header — if the user typed
        #           ``pcol="P"`` and the header contains ``P`` (in any
        #           case), that column is used regardless of what
        #           other candidates (e.g. ``IHS``, ``LOGPVALUE``)
        #           happen to sit earlier in the header.  Fixes the
        #           bug where a header like
        #           ``SNP CHR POSITION IHS LOGPVALUE P BH_adj_P``
        #           resolved p-value to ``IHS`` because the leftmost
        #           header column matching the candidate *set* wins
        #           and ``IHS`` had been added to the p-value
        #           candidates.
        #   Pass 2: leftmost header column matching any built-in
        #           candidate (previous behaviour).
        def _resolve(hint: Optional[str], candidates: list[str],
                     kind: str) -> Optional[str]:
            if hint:
                _hint_l = str(hint).strip().lower()
                for c in hdr:
                    if c.lower() == _hint_l:
                        return c
                # Hint typed a name that isn't in the header — fail
                # fast rather than silently picking a different column.
                sys.exit(
                    f"Error: --{kind} entry for {fpath} references "
                    f"column {hint!r}, but the file header does not "
                    f"contain that column.\n  Header: {hdr}"
                )
            for c in hdr:
                if c in set(candidates):
                    return c
            return None

        chrom_col = _resolve(chrom, chr_candidates, "chrom")
        pos_col   = _resolve(pos,   pos_candidates, "pos")
        snp_col   = _resolve(snp,   snp_candidates, "snp")
        pcol_col  = _resolve(user_pcol, pvl_cands,  "pcol")

        _missing = [
            k for k, v in (
                ("chrom", chrom_col), ("pos", pos_col),
                ("snp",  snp_col),    ("pcol", pcol_col),
            ) if v is None
        ]
        if _missing:
            sys.exit(
                f"Error: could not find required column(s) "
                f"{_missing} in {fpath}.\n  Header: {hdr}"
            )

        # Detect build column in this file's header.
        #
        # Precedence (highest wins):
        #   1. Per-file ``col:<name>`` token from ``build_list``
        #      (user says "for this file, look at THIS column").
        #   2. Auto-detected column matching one of ``bld_candidates``.
        # If (1) resolves to a name that doesn't exist in the file's
        # header we fail fast — the user asked for something specific
        # and shouldn't get silent fallback.
        bcol = None
        _tok = _build_tokens[idx] if idx < len(_build_tokens) else None
        if _tok is not None and _tok[0] == "column":
            _requested = _tok[1]
            if _requested is None:
                # Bare `col` / `col:` / `C` / etc. — user is saying
                # "use this file's build column" but didn't name it.
                # Fall back to the standard candidate list; only
                # error out if nothing matches.
                for c in hdr:
                    if c in set(bld_candidates):
                        bcol = c
                        break
                if bcol is None:
                    sys.exit(
                        f"Error: --build entry for {fpath} is a bare "
                        "column marker (no explicit column name), and "
                        "none of the standard build-column names "
                        f"({sorted(set(bld_candidates))[:8]}...) appears "
                        f"in the file's header.\n"
                        f"  Header: {hdr}\n"
                        "  Use `col:<colname>` to name the column "
                        "explicitly, or a literal build "
                        "(hg18/hg19/hg38) to stamp every row."
                    )
            else:
                # Match case-insensitively against the header — file
                # authors are often inconsistent about column casing.
                _hits = [c for c in hdr if c.lower() == _requested.lower()]
                if not _hits:
                    sys.exit(
                        f"Error: --build entry for {fpath} references "
                        f"column {_requested!r}, but the file's header "
                        f"does not contain that column.\n  Header: {hdr}"
                    )
                bcol = _hits[0]
        else:
            for c in hdr:
                if c in set(bld_candidates):
                    bcol = c
                    break

        if bcol is not None:
            # File has an explicit build column — use it
            old_cols = [chrom_col, pos_col, snp_col, pcol_col, bcol]
            new_cols = {
                chrom_col: "CHR",
                pos_col:   "POS",
                snp_col:   "SNP",
                pcol_col:  "P",
                bcol:      "BUILD",
            }
            col_dtypes = {
                chrom_col: 'category',
                pos_col:   object,
                snp_col:   str,
                pcol_col:  float,
                bcol:      'category',
            }
            sumstats_hdr_dic[name] = [old_cols, col_dtypes, new_cols, file_sep]

        elif (idx < len(_build_tokens)
              and _build_tokens[idx] is not None
              and _build_tokens[idx][0] == "literal"
              and _build_tokens[idx][1] is not None):
            # No build column, but a literal per-file build was
            # supplied via --build (e.g. "hg19").  The ``col:`` case
            # is handled by the branch above via the ``bcol`` path;
            # only literal build tokens flow through here.
            old_cols = [chrom_col, pos_col, snp_col, pcol_col]
            new_cols = {
                chrom_col: "CHR",
                pos_col:   "POS",
                snp_col:   "SNP",
                pcol_col:  "P",
            }
            col_dtypes = {
                chrom_col: 'category',
                pos_col:   object,
                snp_col:   str,
                pcol_col:  float,
            }
            sumstats_hdr_dic[name] = [
                old_cols, col_dtypes, new_cols, file_sep, _build_tokens[idx][1]
            ]

        else:
            # No build info at all
            old_cols = [chrom_col, pos_col, snp_col, pcol_col]
            new_cols = {
                chrom_col: "CHR",
                pos_col:   "POS",
                snp_col:   "SNP",
                pcol_col:  "P",
            }
            col_dtypes = {
                chrom_col: 'category',
                pos_col:   object,
                snp_col:   str,
                pcol_col:  float,
            }
            sumstats_hdr_dic[name] = [old_cols, col_dtypes, new_cols, file_sep]

    def _has_build_info(info: list) -> bool:
        """A file has build info when either (a) its header had a build
        column (which is stored as a fifth entry in ``old_cols``), or
        (b) a per-file build was supplied via ``--build`` (stored as a
        fifth entry in the top-level list)."""
        old_cols = info[0]
        return len(old_cols) == 5 or len(info) == 5

    if not any(_has_build_info(info) for info in sumstats_hdr_dic.values()):
        # Neither build column nor --build was available for any file
        logger.warning(
            "No build column or --build values detected. Summary stats will "
            "be plotted in their native coordinate systems. If your data "
            "are in different coordinate systems, combining them in one plot "
            "is not advisable, especially if ``--annotate`` is set!"
        )

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
# Memory usage
# ---------------------------------------------------------------------------
def _get_memory_usage(mem_df: int):
    if mem_df > 1e6:
        df_mem = mem_df / 1e9
        unit = 'GB'
    else:
        df_mem = mem_df / 1e6
        unit = 'MB'
    if df_mem >= 0.001 and df_mem < 1:
        df_mem = df_mem * 1000
        unit = 'MB'
    if df_mem < 0.001:
        df_mem = df_mem * 1000
        unit = 'KB'

    return f"{df_mem:.3g} {unit}"


# ---------------------------------------------------------------------------
# Density-aware "auto" thinning for Manhattan / circular plotting
# ---------------------------------------------------------------------------

def auto_thin_for_manhattan(
    df: "pd.DataFrame",
    keep_threshold: float = 2.0,
    max_below: int = 200_000,
    logp: bool = True,
    logp_col: str = "logP",
    p_col: str = "P",
    seed: int = 42,
) -> "pd.DataFrame":
    """Density-aware sub-sampling for Manhattan-style scatter plots.

    Inspired by ``gwaslab``'s default behaviour, this helper preserves *every*
    variant whose "interestingness" signal is at or above ``keep_threshold``
    (so peaks, suggestive hits, genome-wide-significant hits, and extreme
    selection-scan values are kept verbatim) and uniformly sub-samples the
    dense bulk below the threshold down to at most ``max_below`` rows in
    total.  For a 10 M-variant scan with the defaults below, this typically
    cuts the plotted point count from 10 M to ~200 K + a few hundred
    peaks — visually indistinguishable above the suggestive band, but two
    orders of magnitude faster to render.

    Two modes, switched by *logp*:

    * **P-value mode** (*logp=True*, the default).  ``signal = -log10(P)``.
      ``keep_threshold`` is in ``-log10(P)`` units (default ``2.0``,
      i.e. ``P <= 0.01``).  Variants with ``-log10(P) >= keep_threshold``
      are all retained.
    * **Raw-statistic mode** (*logp=False*).  ``signal = |value|`` of
      *p_col* — the column carrying the test statistic.  This is the
      right mode for non-p-value scans such as iHS, XP-EHH, F_ST,
      Fay & Wu's H, Tajima's D, etc., where "interesting" means large
      magnitude (positive or negative).  ``keep_threshold`` is then in
      the units of the underlying statistic (default still ``2.0``,
      which is a sensible cutoff for standardised selection scans;
      override with e.g. ``0.05`` for F_ST).

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame.  In p-value mode, must contain either *logp_col*
        (preferred) or *p_col*.  In raw-statistic mode, must contain
        *p_col*.  When the relevant column is absent, *df* is returned
        unchanged.
    keep_threshold : float, optional
        Threshold above which all variants are retained.  Interpreted in
        ``-log10(P)`` units when *logp=True* (default ``2.0``), or in the
        natural units of the underlying statistic when *logp=False*.
    max_below : int, optional
        Maximum number of below-threshold rows to retain, sampled
        uniformly at random.  Default ``200_000``.
    logp : bool, optional
        When ``True`` (default), interpret the data as p-values and use
        ``-log10(P)`` as the signal.  When ``False``, treat *p_col* as a
        raw statistic and use ``|value|`` as the signal.
    logp_col : str, optional
        Name of the precomputed ``-log10(P)`` column for p-value mode.
        Default ``'logP'``.
    p_col : str, optional
        Name of the raw p-value column (p-value mode) or test-statistic
        column (raw-statistic mode).  Default ``'P'``.
    seed : int, optional
        Seed for the RNG used to sub-sample the bulk.  Default ``42``.

    Returns
    -------
    pandas.DataFrame
        Sub-sampled view of *df* preserving the original index ordering.
        When the below-threshold count is already <= *max_below*, the
        input is returned unchanged.

    Examples
    --------
    GWAS p-values (default):

    >>> thinned = auto_thin_for_manhattan(df, keep_threshold=2.0)

    iHS / XP-EHH (signed selection statistics, ``|value|`` >= 2):

    >>> thinned = auto_thin_for_manhattan(
    ...     df, logp=False, keep_threshold=2.0, p_col="iHS",
    ... )

    F_ST (unsigned, 0–1, outlier cutoff e.g. 0.05):

    >>> thinned = auto_thin_for_manhattan(
    ...     df, logp=False, keep_threshold=0.05, p_col="FST",
    ... )
    """
    if df is None or len(df.index) == 0:
        return df

    if logp:
        # p-value mode: use precomputed logP if present, else derive it.
        if logp_col in df.columns:
            signal = df[logp_col].to_numpy()
        elif p_col in df.columns:
            with np.errstate(divide="ignore", invalid="ignore"):
                signal = -np.log10(df[p_col].to_numpy())
        else:
            return df
    else:
        # Raw-statistic mode: |value| of the test-statistic column.  Works
        # for signed statistics (iHS, XP-EHH, Fay & Wu's H, Tajima's D) as
        # well as unsigned ones (F_ST).
        if p_col not in df.columns:
            return df
        signal = np.abs(df[p_col].to_numpy(dtype=float))

    above = signal >= keep_threshold
    below_idx = np.flatnonzero(~above & np.isfinite(signal))

    if below_idx.size <= max_below:
        return df

    rng = np.random.default_rng(seed)
    keep_below = rng.choice(below_idx, size=max_below, replace=False)

    keep_mask = above.copy()
    keep_mask[keep_below] = True

    # Preserve the input DataFrame's positional ordering so chromosomes
    # remain sorted as the caller left them.
    return df.iloc[np.flatnonzero(keep_mask)].copy()


# Resolve significance line
def process_signif_line(
    signif_line, #: Union[float, bool, None], 
    calculated_threshold: float
) -> Optional[float]:
    """
    Resolves the target y-value for the significance line.

    Parameters
    ----------
    signif_line : float, bool, or None
        Option passed from CLI/API (False, None, True, or custom float).
    calculated_threshold : float
        The auto-calculated threshold derived from the dataset.

    Returns
    -------
    float or None
        The exact numerical y-value to draw, or None if suppressed.
    """
    # 1. Suppressed state (user passed False or None)
    if signif_line is False or signif_line is None:
        return None

    # 2. Flag-only state (user passed --signif-line without a value -> True)
    if signif_line is True:
        return calculated_threshold

    # 3. Explicit numerical value state (user passed e.g. 5e-8)
    if isinstance(signif_line, (int, float)):
        return float(signif_line)

    return None

# ---------------------------------------------------------------------------
# Main loader
# ---------------------------------------------------------------------------

def load(
    sum_stats: list[str],
    labels: list[str],
    logp: bool = False,
    trim_pval: Optional[float] = None,
    file_info: Optional[dict] = None,
    sort_tracks: Optional[str] = None,
    table_out: Optional[str] = None,
    signif_threshold: Optional[float] = None,
    signif_line: Optional[bool | float] = None,
    suggest_threshold: Optional[float] = 1e-5,
    highlight: Optional[bool] = False,
    highlight_thresh: Optional[float] = None,
    resources: Optional[ResourceConfig] = None,
    compute_pvals: bool = False,
    auto_thin: bool = True,
    auto_thin_threshold: float = 2.0,
    auto_thin_max_below: int = 200_000,
    cache: bool = False,
    cache_dir: str | os.PathLike = ".pycmplot_cache",
    resume: bool = True,
):
    """Load summary statistics, run liftover, extract lead SNPs, and compute merged Circos sector sizes.

    This is the primary data-loading function for the plotting pipeline.
    For each track it reads the file using the column mapping from
    ``file_info``, optionally filters by ``trim_pval``, normalises
    chromosome names (``chr`` prefix stripped; ``23`` to ``X``, ``24`` to
    ``Y``, ``M`` / ``MTDNA`` to ``MT``), lifts over hg19 coordinates when a
    build column is present, and extracts lead SNPs. After all tracks are
    loaded it builds the hits summary table, derives significance
    thresholds, optionally sorts tracks, and computes the merged
    sector-size dict consumed by both plotters.

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
    compute_pvals : bool, optional
        When ``True``, the full untrimmed p-value array is materialised
        per track and returned in ``bundle['pvals']`` for QQ plotting.
        When ``False`` (the default as of 0.4.0), ``bundle['pvals']``
        contains ``None`` per label and the ~80 MB-at-10M copy is
        skipped entirely.  **Python-API callers that plan to draw a QQ
        figure must set this to** ``True`` **explicitly** — otherwise
        the QQ plotters will raise ``ValueError`` with a message
        pointing back to this argument.  The CLI entry point sets this
        automatically from ``bool(args.qq_plot)``, so ``pycmplot --qq``
        just works without user intervention.

        .. versionchanged:: 0.4.0
           Default changed from ``True`` to ``False``.  Callers that
           relied on the old default and did not set ``compute_pvals``
           explicitly will now see ``None`` per label; pass
           ``compute_pvals=True`` to restore the previous behaviour.
    auto_thin : bool, optional
        Enable density-aware sub-sampling of the null background before
        rendering.  Default ``True``.  See :func:`auto_thin_for_manhattan`
        for the algorithm.
    auto_thin_threshold : float, optional
        ``-log10(P)`` (when *logp*) or ``|statistic|`` (otherwise) at or
        above which every variant is retained verbatim.  Default ``2.0``.
    auto_thin_max_below : int, optional
        Cap on the number of below-threshold "background" variants
        retained per track during auto-thinning.  Default ``200_000``.
    cache : bool, optional
        Enable the per-track Stage-1 cache.  When ``True``, each track's
        post-load / post-liftover / post-thinning DataFrame is written to
        ``<cache_dir>/tracks/<label>.parquet`` alongside its lead-SNP
        table and (when *compute_pvals* is ``True``) its full raw
        p-value array.  Subsequent invocations skip Stage 1 entirely for
        tracks whose SHA-256 + parameter fingerprint matches the cached
        key.  A parameter change or file rewrite invalidates the affected
        entries automatically.  Default ``False``.
    cache_dir : str or os.PathLike, optional
        Directory used to persist the cache.  Created lazily on first
        write.  Default ``".pycmplot_cache"`` (relative to the current
        working directory when the loader is invoked).
    resume : bool, optional
        Reserved for future use.  Currently a no-op because completed
        tracks are atomically committed to the cache one-at-a-time, so
        a run that crashes on track *N* automatically resumes from track
        *N* on the next invocation of the loader with ``cache=True`` —
        the "resume" semantics are inherent in the on-disk layout and
        do not need to be opted into.  Default ``True``; setting to
        ``False`` currently has no effect but is reserved so the API
        can grow a stricter "no resume, always regenerate" mode later.

    Returns
    -------
    dict
        A dictionary with the following keys:

        * ``'sectors'`` — ``dict`` mapping ``chromosome → [min_pos, max_pos]``
          across all tracks, in natural chromosome order (``'1'``, ``'2'``,
          …, ``'X'``, ``'Y'``), with a ``'Spacer1'`` entry appended for
          y-axis labelling.
        * ``'dfs'`` — ``dict`` mapping ``label → [DataFrame, n_chroms]``.
          Each DataFrame contains canonical columns ``CHR``, ``POS``,
          ``SNP``, ``P``, ``LABEL`` and optionally ``logP``, ``BUILD``,
          ``OLD_POS``, ``OLD_BUILD`` (when a build column and liftover
          were applied).
        * ``'annot'`` — :class:`pandas.DataFrame` containing the clumped
          locus summary with nearest-gene annotations.  Empty when no
          variants pass the significance threshold.
        * ``'lines'`` — ``list`` of ``{'genome': float, 'suggestive': float}``
          dicts, one per track, in the final sorted order.
        * ``'pvals'`` — ``dict`` mapping ``label → numpy.ndarray`` of raw
          (un-trimmed) p-values for QQ plotting.

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
    >>> result = get_sumstats_and_merged_sector_list(
    ...     sum_stats=files,
    ...     labels=labels,
    ...     logp=True,
    ...     trim_pval=0.01,
    ...     file_info=file_info,
    ...     signif_threshold=5e-8,
    ... )
    >>> sorted(result.keys())
    ['annot', 'dfs', 'lines', 'pvals', 'sectors']
    >>> list(result["sectors"].keys())[:4]
    ['1', '2', '3', '4']
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
    snp_counts: dict[str, np.ndarray | pd.Series] = {}
    signif_lines: list[dict[str, float]] = []
    all_lead_snps: list[pd.DataFrame] = []

    # ------------------------------------------------------------------
    # Optional per-track cache (Stage-1 checkpointing).  When cache=True
    # we hash each raw file + the Stage-1 parameters and skip the CSV
    # read / dtype conversion / chromosome normalisation / liftover /
    # auto-thinning / lead extraction for any label whose cached
    # DataFrame matches the current key.  Data cached across runs but
    # invalidated automatically when any parameter changes (see
    # :func:`pycmplot.cache.compute_cache_key`).
    # ------------------------------------------------------------------
    _track_cache = None
    if cache:
        from pycmplot import __version__ as _pcm_version
        from pycmplot.cache import TrackCache, compute_cache_key, sha256_file
        _track_cache = TrackCache(cache_dir, _pcm_version)

    def _resolve_signif_lines_entry(label: str) -> dict[str, float]:
        """Recompute the per-track signif_lines dict from current params.

        Kept outside the cache payload so ``--signif_line`` /
        ``--suggest_threshold`` / ``--logp`` remain live-tunable
        between runs without invalidating the cache.
        """
        n_local = int(snp_counts[label]) or 1
        _suggest = 1e-5 if suggest_threshold is None else suggest_threshold
        if logp:
            _suggest = -np.log10(_suggest)
        _resolved = max(0.05 / n_local, 5e-8)
        if signif_line not in (False, None):
            if isinstance(signif_line, (int, float)) and not isinstance(signif_line, bool):
                _resolved = float(signif_line)
        if logp and _resolved < 1:
            _resolved = -np.log10(_resolved)
        return {"genome": _resolved, "suggestive": _suggest}

    # Iterate the user-supplied ``labels`` list rather than the set
    # intersection of the two dicts.  ``dict.keys() & other.keys()``
    # returns a plain ``set``, whose iteration order is *not*
    # deterministic across runs — that broke cache hits for the second
    # and later tracks because ``signif_threshold`` is auto-computed
    # from the first-loaded track's SNP count and then fed into every
    # subsequent track's cache_key.  Fixing the iteration order keeps
    # cache_keys stable across warm re-runs.
    _ordered_labels = [
        label for label in labels
        if label in sumstats and label in (file_info or {})
    ]
    # Collect each track's full cache_key so we can derive a
    # hits-overlay group_key later on (see
    # :func:`pycmplot.cache.compute_hits_group_key` and the multi-panel
    # hits-overlay block near the end of this function).  Keying the
    # group on the *cache_key* (which includes every Stage-1 param)
    # means each ``(files, parameters)`` combination gets its own
    # hits overlay — this matches the per-track cache-invalidation
    # semantics and the reproducibility mental model.  A parameter
    # change therefore spawns a fresh overlay; user edits from an
    # earlier setting linger on disk but are not consulted (users who
    # want to carry a hand-edit forward can copy rows across manually).
    _track_cache_keys: list[str] = []

    # ------------------------------------------------------------------
    # Cross-file build detection (group-level liftover trigger)
    # ------------------------------------------------------------------
    # Before pycmplot 0.4.x the liftover check inside the per-track
    # loop only inspected each file's *own* BUILD column and fired
    # liftover when that single file had mixed hg19/hg38 rows.  When
    # two files were each *internally* single-build but had *different*
    # declared builds — e.g. one hg19 and one hg38 — no liftover fired,
    # every track was plotted in its native coordinate system, and
    # highlight lines drawn across tracks no longer traced the same
    # locus.  Fix: scan the whole loader group up front, and if the
    # group contains more than one distinct build, unify to hg38 by
    # lifting every hg18/hg19 track (an hg38-only or single-build
    # group is unchanged).
    #
    # ``_declared_build_for`` looks in two places, in order:
    #   1. The ``build_list=`` value stored in file_info[label][4]
    #      (populated from the API kwarg or the CLI ``--build``).
    #   2. A light-touch peek at the file's BUILD column (up to 2000
    #      rows), for files that declare builds *inside* the sumstats
    #      TSV instead of via the loader argument.
    # Unknown / unpeekable files return ``None`` and are ignored in
    # the mixed-build test — we only trigger when we can *prove* the
    # group is mixed.  The peek uses ``pd.read_csv(nrows=2000)`` which
    # is negligible even against 100 M-row files.
    _BUILD_MAP = {
        "hg38": "hg38", "grch38": "hg38", "b38": "hg38", "38": "hg38",
        "hg19": "hg19", "grch37": "hg19", "b37": "hg19", "19": "hg19",
        "hg18": "hg18", "ncbi36": "hg18", "18": "hg18",
    }
    _BUILD_COL_CANDIDATES = ("BUILD", "Genome", "Genome_Build", "Genome-build")

    def _declared_build_for(label: str) -> Optional[str]:
        # (1) explicit ``build_list=`` entry — but only when it's a
        # literal build.  A ``col:<name>`` token means "look at that
        # column in the file", so it routes to the peek path with
        # the user-supplied column name overriding the default
        # candidate list.
        _b = None
        _forced_col = None
        try:
            _b = file_info[label][4]
        except (IndexError, KeyError, TypeError):
            _b = None
        if _b:
            _s = str(_b).strip()
            _prefixed = False
            for _p in ("col:", "COL:", "column:", "Column:", "C:", "c:"):
                if _s.startswith(_p):
                    _forced_col = _s[len(_p):].strip() or None
                    _prefixed = True
                    break
            # Bare-marker forms (``c`` / ``col`` / ``column``) also
            # mean "look at the file's build column" without pinning
            # a specific name.  Both prefixed-empty and bare forms
            # fall through to auto-detection below.
            if not _prefixed and _s.lower() in ("c", "col", "column"):
                _forced_col = None
                _prefixed = True
            if not _prefixed:
                return _BUILD_MAP.get(_s.lower())
        # (2) peek at the file's BUILD column
        try:
            _path = sumstats[label][0]
            _sep = file_info[label][3]
        except (IndexError, KeyError, TypeError):
            return None
        try:
            _peek = pd.read_csv(_path, sep=_sep, nrows=2000)
        except Exception:
            return None
        if _forced_col is not None:
            # User told us exactly which column to read; look case-
            # insensitively and fall through to None if the column
            # isn't there (``prep_pycmplot_input_info`` will already
            # have raised a clean error before we reach this point).
            _match = next(
                (c for c in _peek.columns if c.lower() == _forced_col.lower()),
                None,
            )
            _col = _match
        else:
            _col = next(
                (c for c in _BUILD_COL_CANDIDATES if c in _peek.columns), None,
            )
        if _col is None:
            return None
        _values = (
            _peek[_col].dropna().astype(str).str.strip().str.lower().unique()
        )
        _normed = {_BUILD_MAP[v] for v in _values if v in _BUILD_MAP}
        if not _normed:
            return None
        # A file with per-row mixed builds already trips the existing
        # per-file liftover check downstream — nothing extra needed
        # from us.  When only one distinct build lives in the column,
        # advertise that single build so the group-level detector
        # can see whether the tracks disagree.
        return next(iter(_normed)) if len(_normed) == 1 else "mixed"

    _declared_builds = {
        label: _declared_build_for(label) for label in _ordered_labels
    }
    _group_build_set = {
        b for b in _declared_builds.values()
        if b in ("hg18", "hg19", "hg38")
    }
    # More than one distinct declared build across the group means the
    # tracks live in different coordinate systems and must be unified
    # to hg38 for cross-track highlight lines to line up.
    _cross_file_liftover = len(_group_build_set) > 1
    if _cross_file_liftover:
        logger.info(
            "Cross-file mixed builds detected across the loader group "
            "(%s); hg18/hg19 tracks will be lifted over to hg38 so that "
            "highlight lines and shared coordinates align.",
            "/".join(sorted(_group_build_set)),
        )

    for label in _ordered_labels:
        # ---- Cache lookup fast path ----------------------------------
        _cache_key = None
        if _track_cache is not None:
            _raw_path = sumstats[label][0]
            try:
                _raw_hash = sha256_file(_raw_path)
            except OSError as _exc:
                logger.warning("Cache: cannot hash %s (%s); bypassing cache for this track.",
                               _raw_path, _exc)
                _raw_hash = None

            if _raw_hash is not None:
                # NOTE: compute_pvals is deliberately EXCLUDED from the
                # cache_key.  The cached DataFrame is unaffected by
                # whether pvals are also materialised — pvals are a
                # sidecar that either exists or doesn't.  Keeping the
                # key stable across compute_pvals toggles means a user
                # can flip QQ on/off between runs without invalidating
                # the (expensive) main track cache.  See the fallback
                # below for the case where compute_pvals=True but the
                # sidecar is absent.
                _cache_params = dict(
                    raw_sha256=_raw_hash,
                    logp=bool(logp),
                    trim_pval=None if trim_pval is None else float(trim_pval),
                    auto_thin=bool(auto_thin),
                    auto_thin_threshold=float(auto_thin_threshold),
                    auto_thin_max_below=int(auto_thin_max_below),
                    highlight=bool(highlight),
                    highlight_thresh=None if highlight_thresh is None else float(highlight_thresh),
                    signif_threshold=None if signif_threshold is None else float(signif_threshold),
                    build=file_info[label][4] if len(file_info[label]) > 4 else None,
                )
                _cache_key = compute_cache_key(**_cache_params)
                _track_cache_keys.append(_cache_key)
                _hit = _track_cache.get(label, _cache_key)
                if _hit is not None:
                    # If pvals are needed but weren't cached in a
                    # previous run (e.g. that run had compute_pvals=False),
                    # treat this as a miss so Stage 1 re-runs and
                    # populates the pvals sidecar for next time.
                    if compute_pvals and _hit.pvals is None:
                        logger.info(
                            "Cache MISS %s (compute_pvals=True but no "
                            "pvals sidecar; re-parsing to populate it).",
                            label,
                        )
                    else:
                        sumstats_loaded[label] = [_hit.df, int(_hit.extras.get("n_chroms", 0))]
                        pval_dict[label] = _hit.pvals if compute_pvals else None
                        snp_counts[label] = int(_hit.extras.get("snp_count", 0))
                        all_lead_snps.append(
                            _hit.leads if _hit.leads is not None else pd.DataFrame()
                        )
                        # Mirror the cold-path side-effect on
                        # ``signif_threshold`` (which is auto-computed
                        # from the first track's SNP count and then
                        # feeds the *next* track's cache_key).  Without
                        # this, subsequent tracks on a warm re-run would
                        # compute their cache_keys against a stale
                        # ``signif_threshold=None`` and MISS.
                        if signif_threshold is None:
                            _n_hit = int(snp_counts[label]) or 1
                            signif_threshold = max(0.05 / _n_hit, 5e-8)
                        signif_lines.append(_resolve_signif_lines_entry(label))
                        continue

        # ---- Regular Stage-1 processing --------------------------------
        # Any exception raised inside the loop body will propagate up and
        # abort the run.  Because tracks that already completed have been
        # written to the cache at the end of their own iteration, a
        # subsequent invocation of this loader with cache=True will hit
        # those cached tracks and re-run Stage 1 only from the point of
        # failure onward — that's the "resume" semantics from to-do.md.
        sumstat_cols   = file_info[label][0]
        sumstat_dtypes = file_info[label][1]
        sumstat_newcols= file_info[label][2]
        sep            = file_info[label][3]

        build = None
        try:
            build      = file_info[label][4]
        except Exception:
            pass

        logger.info("Loading %s [%s] ...", label, sumstats[label][0])
        # Prefer the ``pyarrow`` CSV engine when available — it's typically
        # 1.5–3× faster than the default C engine on large GWAS summary
        # statistics, especially for files with many numeric columns.  The
        # ``pyarrow`` engine ignores the ``dtype=`` argument for category
        # casts, so we cast the chromosome column to ``Categorical`` after
        # the read.  Any pyarrow-side failure (missing package, unsupported
        # option, etc.) falls back to the default C engine.
        read_kwargs = dict(
            filepath_or_buffer=sumstats[label][0],
            sep=sep,
            header=0,
            usecols=sumstat_cols,
        )
        try:
            df = pd.read_csv(
                **read_kwargs,
                engine="pyarrow",
                dtype_backend="numpy_nullable",
            ).rename(columns=sumstat_newcols)
        except (ImportError, ValueError, TypeError):
            df = pd.read_csv(
                **read_kwargs,
                dtype=sumstat_dtypes,
            ).rename(columns=sumstat_newcols)


        # Normalise chromosome names — done once here and stored as a
        # ``Categorical`` with ``CHROM_ORDER`` as the canonical category
        # set.  Downstream plotting code can recognise this dtype and skip
        # repeating the (string-heavy) normalisation, and any aliasing /
        # filtering on chromosome name becomes integer-code work rather
        # than per-element Python string ops.
        #
        # Critically, when CHR comes in as a ``Categorical`` (the dtype we
        # request in ``prep_pycmplot_input_info`` for any non-build file)
        # the actual normalisation is applied to the **categories**, not
        # to the underlying N-row code array.  That turns a 500K (or 10M)
        # per-element ``str.replace + str.upper + replace`` chain into the
        # equivalent work on ~25 distinct chromosome labels.
        logger.info('Normalizing chromosome names {"23": "X", "24": "Y", "M": "MT", "MTDNA": "MT"} ...')
        chr = df["CHR"]

        if pd.api.types.is_numeric_dtype(chr):
            chr = chr.astype("Int64").astype(str)
            chr = chr.replace({
                "23": "X",
                "24": "Y",
            })
        else:
            chr = (
                chr.astype(str)
                .str.upper()
                .str.replace("CHR", "", regex=False)
                .replace({
                    "23": "X",
                    "24": "Y",
                    "M": "MT",
                    "MTDNA": "MT",
                })
            )

        df["CHR"] = pd.Categorical(
            chr,
            categories=CHROM_ORDER,
            ordered=True,
        )

        # Coerce POS to numeric, drop rows that fail to parse, then store
        # as plain int64 (not nullable ``Int64``) so downstream arithmetic
        # / ``max()`` / categorical groupby reductions cannot leak
        # ``pd.NA`` into mixed-type expressions.
        df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
        df = df.dropna(subset=["POS"]).copy()
        df["POS"] = df["POS"].astype("int64")
        pre_trim_mem = _get_memory_usage(df.memory_usage(deep=True).sum())
        pre_trim_vars = len(df.index)
        logger.info("Loaded %s variants from summary stat file, using %s of memory", pre_trim_vars, pre_trim_mem)

        # Get dict of p-values for QQ-plotting before applying trim_pval.
        # Computing this is only meaningful when a QQ plot will actually be
        # rendered downstream; for Manhattan-only or circular-only runs we
        # skip the ~80 MB copy at 10 M variants entirely.
        if compute_pvals:
            logger.info("Extracting raw p-values for QQ-plotting ...")
            pval_dict[label] = df["P"].dropna().astype(float).values
        else:
            pval_dict[label] = None

        # Get SNP counts for significance threshold calculation
        snp_counts[label] = len(df["P"].dropna().astype(float).values)

        # ------------------------------------------------------------------
        # Signed-statistic detection (iHS / XP-EHH / Fay & Wu's H / etc.)
        # ------------------------------------------------------------------
        # When ``logp=False`` and the score column carries negative values,
        # switch into signed mode: add a companion |value| column
        # (``P_UNSIGNED``) and use *that* for lead-SNP extraction and
        # highlight-window selection.  The signed original stays in ``P``
        # and drives the y-axis, so peaks and troughs both stay visible.
        # ``auto_thin_for_manhattan`` already operates on ``|value|`` so
        # the two paths agree.
        _signed = (not logp) and bool(
            np.any(df["P"].to_numpy(dtype=float, na_value=0.0) < 0)
        )
        if _signed:
            df["P_UNSIGNED"] = df["P"].abs()
            _score_col = "P_UNSIGNED"
            _score_asc = False
            if signif_threshold is None:
                raise ValueError(
                    "signif_threshold is required for signed selection "
                    "statistics (iHS/XP-EHH/etc.). The p-value fallback "
                    "``max(0.05/N, 5e-8)`` is meaningless on |value| — "
                    "pass e.g. ``signif_threshold=4`` for iHS."
                )
        else:
            _score_col = None
            _score_asc = None

        # Derive significance/suggestive thresholds
        n = snp_counts[label]
        if signif_threshold is None:
            #last_label = list(sumstats_loaded)[-1]
            signif_threshold = max(0.05 / n, 5e-8)

        if suggest_threshold is not None:
            suggest_line = suggest_threshold
        else:
            suggest_line = 1e-5
        if logp:
            suggest_line = -np.log10(suggest_line)

        # Initialise from the effective significance threshold — this
        # is the p-value auto-fill for unsigned data (already resolved
        # above by the ``if signif_threshold is None`` block), or the
        # user-supplied |value| cutoff for signed statistics.  The old
        # code re-derived from ``max(0.05/n, 5e-8)`` which silently
        # discarded ``signif_threshold=4`` on signed loads, leaving the
        # drawn reference line at ~5e-8 (a p-value) instead of at 4
        # (the |iHS| cutoff the user asked for).
        resolved_signif_line = float(signif_threshold)

        # Check if signif_line was requested (i.e. not False and not None)
        if signif_line not in (False, None):
            if signif_line is True:
                pass
            elif isinstance(signif_line, (int, float)):
                # Use user-supplied custom float
                resolved_signif_line = float(signif_line)

        if logp and resolved_signif_line < 1:
            resolved_signif_line = -np.log10(resolved_signif_line)

        # In signed mode ``resolved_signif_line`` is a threshold on
        # |value| — record the mirrored negative sibling so plotters can
        # draw ±threshold bands.  Same for the suggestive line.
        #
        # Clamp each reference-line y-value to the observed data range
        # so a threshold that exceeds the data extremes still draws at
        # the plot edge rather than floating off-screen.  Positive-tail
        # lines are capped at ``max(P)``; negative-tail lines are
        # floored at ``min(P)``.
        if _signed:
            _vals = df["P"].to_numpy(dtype=float)
            _finite = _vals[np.isfinite(_vals)]
            if _finite.size:
                _data_max = float(_finite.max())
                _data_min = float(_finite.min())
                _pos_line = min(float(resolved_signif_line), _data_max)
                _neg_line = max(-float(resolved_signif_line), _data_min)
                _sug_pos = (
                    None if suggest_line is None
                    else min(float(suggest_line), _data_max)
                )
                _sug_neg = (
                    None if suggest_line is None
                    else max(-float(suggest_line), _data_min)
                )
            else:
                _pos_line = float(resolved_signif_line)
                _neg_line = -float(resolved_signif_line)
                _sug_pos = suggest_line
                _sug_neg = None if suggest_line is None else -float(suggest_line)
            _line_dict = {"genome": _pos_line, "suggestive": _sug_pos}
            _line_dict["genome_neg"] = _neg_line
            if _sug_neg is not None:
                _line_dict["suggestive_neg"] = _sug_neg
        else:
            _line_dict = {"genome": resolved_signif_line, "suggestive": suggest_line}
        signif_lines.append(_line_dict)

        # Density-aware auto-thinning for Manhattan / circular rendering.
        # Applied after lead-SNP extraction so the leads come from the full
        # dataset, and after liftover so coordinates are final.  Variants
        # at or above ``auto_thin_threshold`` (default ``-log10(P) >= 2``)
        # are kept verbatim, so all suggestive / significant hits and their
        # surrounding LD bumps survive untouched — only the dense null
        # background is sub-sampled.  Skipped automatically when the
        # below-threshold count is already small.
        if auto_thin:
            n_before = len(df.index)
            df = auto_thin_for_manhattan(
                df,
                keep_threshold=auto_thin_threshold,
                max_below=auto_thin_max_below,
                logp=logp,
            )
            n_after = len(df.index)
            if n_after < n_before:
                signal_desc = (
                    "-log10(P)" if logp else "|value| of test statistic"
                )
                logger.info(
                    "Auto-thinning: %s -> %s variants (kept all %s >= %s; "
                    "down-sampled below-threshold background to <=%s).",
                    n_before, n_after, signal_desc,
                    auto_thin_threshold, auto_thin_max_below,
                )

        # Add build column if not exist and build supplied
        BUILD_MAP = {
            "hg38": "hg38", "grch38": "hg38", "38": "hg38",
            "hg19": "hg19", "grch37": "hg19", "b37": "hg19", "19": "hg19",
            "hg18": "hg18", "ncbi36": "hg18", "18": "hg18",
        }
        if build is not None:
            build_key = str(build).strip().lower()
            if build_key not in BUILD_MAP:
                raise ValueError(
                    f"Unsupported genome build '{build}'. "
                    f"Supported options: {', '.join(sorted(set(BUILD_MAP.values())))}."
                )
            df['BUILD'] = BUILD_MAP[build_key]
            df['BUILD'] = df['BUILD'].astype('category')

        # Trim insignificant variants for faster plotting
        if trim_pval:
            logger.info("Excluding variants with p-value less than %s to speed up Manhattan plotting ...", trim_pval)
            df = df[df["P"].astype(float) <= float(trim_pval)]
            post_trim_mem = _get_memory_usage(df.memory_usage(deep=True).sum())
            post_trim_vars = len(df.index)
            logger.info("%s variants remain after trimming, using %s of memory", post_trim_vars, post_trim_mem)
        # this breaks with statistics that have both negative and positive values
        # such as iHS
        #else:
        #    df = df[df["P"].astype(float) <= 1] 

        if logp:
            logger.info("Adding a 'logP' column ...")
            df["logP"] = -np.log10(df["P"])

        df["LABEL"] = label

        # Liftover hg18/hg19 data if needed.
        #
        # ``sumstats_loaded[label]`` is not populated until the very end of
        # this iteration (line ``sumstats_loaded[label] = [df, n_chroms]``),
        # so the result must be assigned to the *local* ``df`` — writing
        # into ``sumstats_loaded[label][0]`` here raised ``KeyError`` (e.g.
        # ``KeyError: 'MCV'``) the first time the liftover branch fired on
        # a given track.
        # Guard the BUILD lookup: not every input file has a BUILD column
        # (single-build workflows don't need one).  Without this check
        # ``df["BUILD"]`` raises ``KeyError`` before we reach the
        # ``if "BUILD" in df.columns`` test.
        if "BUILD" in df.columns:
            builds = df["BUILD"].unique()
        else:
            builds = []
        # Per-file liftover (unchanged): fires when this file's own
        # BUILD column carries hg18 rows, or a mix of hg19 and hg38
        # rows within the same file.
        _needs_lift_file = "BUILD" in df.columns and (
            "hg18" in builds or ("hg19" in builds and "hg38" in builds)
        )
        # Group-level liftover (0.4.x): fires when the whole loader
        # group carries more than one distinct declared build AND this
        # particular file is hg18/hg19.  Without this branch, a
        # ``[hg19, hg38]`` group left every track in its native
        # coordinate system, so cross-track highlight lines drew at
        # different genomic positions (bug reported 2026-09-05).
        _this_declared = _declared_builds.get(label)
        _needs_lift_group = (
            _cross_file_liftover
            and _this_declared in ("hg18", "hg19")
            and "BUILD" in df.columns
            and any(b in builds for b in ("hg18", "hg19"))
        )
        if _needs_lift_file or _needs_lift_group:
            builds_present = sorted(
                b for b in builds if b in {"hg18", "hg19"}
            )
            _reason = (
                "same-file mixed builds"
                if _needs_lift_file
                else "cross-file mixed builds in the loader group"
            )
            logger.info(
                "Converting %s coordinates to hg38 (%s) ...",
                "/".join(builds_present), _reason,
            )
            df = liftover_position(df, resources=resources)

        # get highlight SNPs
        if highlight:
            logger.info("Extracting lead variants and variants to highlight ...")
        else:
            logger.info("Extracting lead variants ...")

        df, leads = get_highlight_snps(
            df=df,
            window=500_000,
            highlight=highlight,
            highlight_thresh=highlight_thresh if highlight_thresh is not None else signif_threshold,
            logp=logp,
            score_col=_score_col,
            ascending=_score_asc,
        )

        ## Lead SNPs
        #logger.info("Extracting lead variants ...")
        #leads = get_lead_snps(
        #    df=sumstats_loaded[label][0],
        #    signif_threshold=signif_threshold or 5e-8,
        #    logp=logp,
        #)

        if not leads.empty:
            if _signed:
                # Signed statistics: keep both tails (|value| >= threshold).
                leads = leads[leads["P_UNSIGNED"] >= float(signif_threshold)]
            else:
                leads = leads[leads["P"] <= signif_threshold]

        all_lead_snps.append(leads)

        # Number of distinct chromosomes (for track sorting)
        n_chroms = len(df["CHR"].unique()) - 1
        sumstats_loaded[label] = [df, n_chroms]

        # ---- Persist to cache (Stage-1 checkpoint) --------------------
        # Runs only when cache=True and we computed a valid cache_key
        # above.  Best-effort: a failure to write the cache logs a
        # warning but does not abort the run.
        if _track_cache is not None and _cache_key is not None:
            try:
                _track_cache.put(
                    label,
                    _cache_key,
                    df,
                    leads=leads if isinstance(leads, pd.DataFrame) and not leads.empty else None,
                    pvals=pval_dict.get(label) if compute_pvals else None,
                    extras={
                        "n_chroms": int(n_chroms),
                        "snp_count": int(snp_counts.get(label, 0)),
                    },
                    raw_source=str(sumstats[label][0]),
                )
            except Exception as _exc:
                logger.warning("Cache write failed for %r (%s); continuing.",
                               label, _exc)

    # Combine lead SNPs and filter to significance threshold
    all_lead_snps_df = (
        pd.concat(all_lead_snps, ignore_index=True).drop_duplicates()
        if all_lead_snps
        else pd.DataFrame()
    )
  
    # ------------------------------------------------------------------
    # Hits table with optional user-editable overlay.
    #
    # When ``cache=True``, the auto-generated table is written to
    # ``<cache_dir>/annotations/hits.tsv`` with a ``source`` column
    # ('auto' | 'user').  Regeneration is skipped when the leads +
    # resource + window key matches the previously cached value.  User-
    # authored rows (``source="user"``) survive every regeneration so
    # analysts can hand-add / hand-edit annotations without touching
    # Python.  See ``pycmplot.cache`` for the file layout.
    # ------------------------------------------------------------------
    if all_lead_snps_df.empty:
        hits_table = pd.DataFrame()
    elif cache:
        from pycmplot.cache import (
            hits_auto_key, read_hits_overlay, resources_fingerprint,
            write_hits_overlay, compute_hits_group_key,
            AUTO_TAG, USER_TAG, SOURCE_COL,
        )
        _res_sig = resources_fingerprint(resources)
        _auto_key = hits_auto_key(
            all_lead_snps_df, resources_signature=_res_sig, window_kb=2_000,
        )
        # ``_group_key`` scopes the hits overlay to THIS loader call's set
        # of tracks.  A second call with the same ``cache_dir`` but a
        # different set of sumstats (a common multi-panel pattern) gets
        # a different group_key and therefore a separate
        # ``hits.<group>.tsv`` file — they don't clobber each other, and
        # each group keeps its own user-authored rows across regenerations.
        # ``_track_cache_keys`` was populated per-track in the main loop.
        _group_key = compute_hits_group_key(_track_cache_keys)

        _existing_overlay, _existing_meta = read_hits_overlay(cache_dir, _group_key)
        _cached_key = (_existing_meta or {}).get("auto_key")

        if _existing_overlay is not None and _cached_key == _auto_key:
            _hits_file = os.path.join(
                str(cache_dir), "annotations", f"hits.{_group_key}.tsv",
            )
            logger.info(
                "Hits overlay: cache HIT (%s rows; user-editable at %s)",
                len(_existing_overlay), _hits_file,
            )
            hits_table = _existing_overlay
        else:
            _auto_hits = get_hits_summary_table(
                leads_df=all_lead_snps_df,
                table_out=table_out,
                window_kb=2_000,
                resources=resources,
            )
            _written = write_hits_overlay(
                cache_dir, _auto_hits,
                auto_key=_auto_key,
                group_key=_group_key,
                preserve_user_from=_existing_overlay,
            )
            # Re-read to pick up any preserved user rows.
            hits_table, _ = read_hits_overlay(cache_dir, _group_key)
            if hits_table is None:
                hits_table = _auto_hits
            logger.info(
                "Hits overlay: regenerated (%s auto rows; edit %s to add "
                "custom loci)",
                int((_auto_hits[SOURCE_COL] == AUTO_TAG).sum())
                if SOURCE_COL in _auto_hits.columns else len(_auto_hits),
                _written,
            )
    else:
        hits_table = get_hits_summary_table(
            leads_df=all_lead_snps_df,
            table_out=table_out,
            window_kb=2_000,
            resources=resources,
        )

    # sort dicts by user-supplied order
    sumstats_loaded = {key: sumstats_loaded[key] for key in labels if key in sumstats_loaded}
    pval_dict = {key: pval_dict[key] for key in labels if key in pval_dict}
    

    # or sort by user option
    if sort_tracks is not None:
        if sort_tracks.lower() == "label":
            sumstats_loaded = dict(sorted(sumstats_loaded.items()))
        else:  # chrom_len: most chromosomes first (descending n_chroms)
            sumstats_loaded = dict(
                sorted(
                    sumstats_loaded.items(),
                    key=lambda item: -int(item[1][1]),
                )
            )
       

    # Compute per-sumstat sector sizes (chrom → [min_pos, max_pos])
    assoc_sector_sizes_list: list[dict] = []
    min_dic_val = None

    logger.info("Computing per-sumstat sector sizes (chrom → [min_pos, max_pos])")
    for df, _n in sumstats_loaded.values():
        assoc = df[~(df["CHR"].str.len() > 2)].copy()
        assoc["POS"] = assoc["POS"].fillna(0).astype(int)

        assoc_dic: dict[str, list] = {}
        for chrom in assoc["CHR"].unique():
            sub = assoc[assoc["CHR"] == chrom]
            lo_val = max(sub["POS"].min() - 1_000_000, 0)
            hi_val = sub["POS"].max()

            # Ensure sector sizes are within chrom ranges if liftover
            #chrom_max = hi_val
            #if liftover:
            #    hg38_chr_lengths = {k.replace("chr",""): v for k, v in hg38_chr_lengths.items()}
            #    chrom_max = hg38_chr_lengths[chrom]
            #hi_val = min(hi_val, chrom_max)

            assoc_dic[str(chrom)] = [lo_val, hi_val]

        min_dic_val = min(assoc_dic.values())
        assoc_sector_sizes_list.append(assoc_dic)

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
        
    logger.info("All processes completed successfully!")
    return {"sectors": merged, "dfs": sumstats_loaded, "annot": hits_table, "lines": signif_lines, "pvals": pval_dict}


# ---------------------------------------------------------------------------
# Backwards-compatible aliases (deprecated in 0.4.x)
# ---------------------------------------------------------------------------
# The two entry points above were renamed from ``prep_pycmplot_input_info``
# → ``prep`` and ``get_sumstats_and_merged_sector_list`` → ``load`` so that
# ``pycmplot as pcm; pcm.prep(...); pcm.load(...)`` reads naturally.  The
# old names remain importable for one release cycle; each call emits a
# single DeprecationWarning and delegates to the new function.  Internal
# pycmplot code was updated to use the new names on the same commit, so
# only user scripts that pinned the old names will see the warning.
from pycmplot._deprecation import _deprecated_alias as _da

prep_pycmplot_input_info = _da(
    prep, old_name="prep_pycmplot_input_info", new_name="prep",
)
get_sumstats_and_merged_sector_list = _da(
    load,
    old_name="get_sumstats_and_merged_sector_list",
    new_name="load",
)
