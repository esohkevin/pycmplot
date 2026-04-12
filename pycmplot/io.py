"""
pycmplot.io
===========
Summary statistics loading, delimiter detection, and sector-size computation.
"""

from __future__ import annotations

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

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# File utilities
# ---------------------------------------------------------------------------

def smart_open(file_path: str):
    """Open a regular or gzip-compressed file transparently."""
    path = Path(file_path)
    if path.suffix == ".gz":
        return gzip.open(file_path, "rt")
    return open(file_path, "r")


def resolve_delimiter(delim: str) -> str:
    """Map a human-readable delimiter name to the actual separator character."""
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
    """Automatically detect the delimiter using :mod:`csv.Sniffer`.

    Returns
    -------
    (delimiter_str, dialect_or_None)
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
    """Return the column names from the first line of *file_path*."""
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
):

    if len(sum_stats) != len(labels):
        sys.exit(
            "Error: number of summary stats files and labels must match.\n"
            f"  Files:  {sum_stats}\n"
            f"  Labels: {labels}"
        )

    # ------------------------------------------------------------------
    # Sumstat, labels str to list
    # ------------------------------------------------------------------
    labels     = [lbl.strip() for lbl in labels.strip().split(",")]
    
    sum_stats  = [s.strip() for s in sum_stats.strip().split(",")]

    # ------------------------------------------------------------------
    # Colours str to list
    # ------------------------------------------------------------------
    colors = [c.strip() for c in colors_raw.strip().split(",")]

    # ------------------------------------------------------------------
    # Linear track heights str to list
    # ------------------------------------------------------------------
    t_heights = [float(x) for x in track_heights.strip().split(",")]

    return sum_stats, labels, colors, t_heights


# ------------------------------------------------------------------
# Random string for output paths
# ------------------------------------------------------------------
def generate_random_string(length):
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

    out_path = Path(output_dir).resolve()

    out_path.mkdir(parents=True, exist_ok=True)

    if plot_title:
        pltitle = re.sub(r"[^a-zA-Z0-9\s]", "", plot_title).replace(" ", "_")
    else:
        pltitle = generate_random_string(10)

    plt_base = str(out_path / f"{pltitle}_{'_'.join(labels)}_{mode.lower()}")

    suffix     = "_logp" if logp else "_pval"

    plt_name   = f"{plt_base}{suffix}.{output_format.lower()}"
    
    table_out  = f"{plt_base}{suffix}_locus_summary_table.tsv"

    return plt_name, table_out



# ---------------------------------------------------------------------------
# input formatter
# ---------------------------------------------------------------------------
def prep_pycmplot_input_info(
    sum_stats: list[str],
    labels: list[str],
    delim: Optional[str] = None,
    chrom: Optional[str] = None,
    pos: Optional[str] = None,
    snp: Optional[str] = None,
    pcol: Optional[str] = None,
    build: Optional[str] = None
):
    """Resolve column names and delimiter

    Parameters
    ----------
    sum_stats:
        List of file paths to GWAS summary statistics (possibly gzip-compressed).
    labels:
        Track labels in the same order as *sum_stats*.
    delim:
        File delimiter (autodetected if omitted)
    chrom:
        Chromosome column
    pos:
        Position column
    snp:
        SNP or Marker ID column
    pcol:
        P-value column
    build:
        Build version column

    Returns
    -------
    {old_columns, column_dtypes, new_columns, delim}

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
    chr_candidates = [chrom, "CHR", "CHROM", "Chromosome", "#CHROM", "#CHR",
                    "Chrom", "chrom", "chr", "chromosome", "#chr", "#chrom"]
    pos_candidates = [pos,   "BP", "POS", "bp", "pos", "Basepair"]
    snp_candidates = [snp,   "SNP", "RSID", "rsID", "MarkerName", "MarkerID",
                    "Predictor", "Marker", "SNPID", "ID"]
    pvl_candidates = [pcol,  "P", "P-value", "Wald_P", "pvalue", "p_val", "pval"]
    bld_candidates = [build, "BUILD", "Genome", "Genome_Build", "Genome-build"]

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
):
    """Load summary statistics and compute merged Circos sector sizes.

    Parameters
    ----------
    sum_stats:
        List of file paths to GWAS summary statistics (possibly gzip-compressed).
    labels:
        Track labels in the same order as *sum_stats*.
    file_info:
        Dict keyed by label; each value is a list
        ``[col_names, col_dtypes, rename_map, sep]``.
    sort_tracks:
        ``'label'`` — sort tracks alphabetically by label.
        ``'chrom_len'`` — sort by number of chromosomes (default).
        ``None`` — preserve input order.
    signif_threshold:
        Threshold of significance to create hits table.
    resources:
        :class:`~pycmplot.resources.ResourceConfig` instance.

    Returns
    -------
    (merged_sector_sizes, sumstats_loaded, hits_table, signif_lines)
    """
    if resources is None:
        resources = default_resources

    from pycmplot.liftover import liftover_position

    # Build a label → file path mapping
    sumstats: dict[str, list] = {
        name: [path] for name, path in zip(labels, sum_stats)
    }

    sumstats_loaded: dict[str, list] = {}
    all_lead_snps: list[pd.DataFrame] = []

    for label in sumstats.keys() & (file_info or {}).keys():
        sumstat_cols   = file_info[label][0]
        sumstat_dtypes = file_info[label][1]
        sumstat_newcols= file_info[label][2]
        sep            = file_info[label][3]

        logger.info("Loading %s from %s …", label, sumstats[label][0])
        df = pd.read_csv(
            sumstats[label][0],
            sep=sep,
            header=0,
            usecols=sumstat_cols,
            dtype=sumstat_dtypes,
        ).rename(columns=sumstat_newcols)

        # Trim insignificant variants for faster plotting
        if trim_pval:
            logger.info("Excluding variants with p-value less than %s ...", trim_pval)
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

    # Optionally sort tracks
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
            hi_val = sub["POS"].max() + 1_000_000
            assoc_dic[str(chrom)] = [lo_val, hi_val]

        min_dic_val = min(assoc_dic.values())
        assoc_sector_sizes_list.append(assoc_dic)

    merged = _merge_min_max_lists(assoc_sector_sizes_list)
    merged = dict(natsort.natsorted(merged.items(), key=lambda item: item[0]))

    if "23" in merged:
        merged["X"] = merged.pop("23")

    # Add spacer sector for y-axis labelling
    if min_dic_val is not None:
        if len(labels) <= 5:
            merged["Spacer1"] = [x + x / 2 for x in min_dic_val]
        else:
            merged["Spacer1"] = [x * 2 for x in min_dic_val]

    return merged, sumstats_loaded, hits_table, signif_lines
