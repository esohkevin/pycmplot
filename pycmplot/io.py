"""
pycmplot.io
===========
Summary statistics loading, delimiter detection, and sector-size computation.
"""

from __future__ import annotations

import csv
import gzip
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
    chrom_col: Optional[str] = None,
    pos_col: Optional[str] = None,
    pval_col: Optional[str] = None,
    logp: bool = False,
    trim_pval: Optional[float] = None,
    snp_col: Optional[str] = None,
    delim: Optional[str] = None,
    file_info: Optional[dict] = None,
    sort_tracks: Optional[str] = "chrom_len",
    table_out: Optional[str] = None,
    highlight: bool = False,
    highlight_thresh: float = 5e-8,
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
    highlight:
        Whether to flag loci for highlighting.
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

        # Lead SNPs / highlight SNPs
        if highlight:
            logger.info("Extracting variants to highlight ...")
            sumstats_loaded[label][0], leads = get_highlight_snps(
                df=sumstats_loaded[label][0],
                window=2_000_000,
                highlight_thresh=highlight_thresh,
                logp=True,
            )
        else:
            leads = get_lead_snps(
                df=sumstats_loaded[label][0],
                highlight_thresh=signif_threshold or 5e-8,
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
