"""
pycmplot.stats
==============
Statistical helper functions for identifying lead SNPs and loci to highlight.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def get_lead_snps(
    df: pd.DataFrame,
    highlight_thresh: float = 5e-8,
    logp: bool = False,
    window: int = 500_000,
) -> pd.DataFrame:
    """Identify independent lead SNPs by greedy distance clumping.

    Starting from the most significant SNP, each subsequent SNP is kept only
    if it is > *window* bp away from all previously kept leads on the same
    chromosome.

    Parameters
    ----------
    df:
        Summary statistics DataFrame containing columns ``CHR``, ``POS``,
        ``P`` (and ``logP`` when *logp* is ``True``).
    highlight_thresh:
        P-value (or −log₁₀(p) when *logp* is ``True``) significance cutoff.
    logp:
        If ``True``, filter and rank by the ``logP`` column instead of ``P``.
    window:
        Clumping window in base-pairs (default 500 kb).

    Returns
    -------
    pd.DataFrame
        Subset of *df* containing only the lead SNPs.
    """
    if logp:
        thresh = -np.log10(float(highlight_thresh))
        sig = df[df["logP"] >= thresh].copy()
        p_col = "logP"
        ascending = False
    else:
        sig = df[df["P"] <= highlight_thresh].copy()
        p_col = "P"
        ascending = True

    sig = sig.sort_values(p_col, ascending=ascending)
    leads: list[pd.Series] = []

    while not sig.empty:
        top = sig.iloc[0]
        leads.append(top)
        sig = sig[
            ~(
                (sig["CHR"] == top["CHR"])
                & (abs(sig["POS"] - top["POS"]) <= window)
            )
        ]

    return pd.DataFrame(leads)


def get_highlight_snps(
    df: pd.DataFrame,
    highlight_thresh: float = 5e-8,
    logp: bool = False,
    window: int = 500_000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Mark all SNPs within *window* bp of a lead SNP.

    Adds an ``in_locus`` boolean column to *df* and returns the annotated
    DataFrame together with the lead SNP DataFrame.

    Parameters
    ----------
    df, highlight_thresh, logp, window:
        See :func:`get_lead_snps`.

    Returns
    -------
    (df_annotated, leads_df)
    """
    df = df.copy()
    df["in_locus"] = False

    leads_df = get_lead_snps(
        df=df,
        highlight_thresh=highlight_thresh,
        logp=False,
        window=window,
    )

    for _, row in leads_df.iterrows():
        min_pos = row["POS"] - window
        max_pos = row["POS"] + window
        chrom = row["CHR"]

        mask = (df["CHR"] == chrom) & (df["POS"] >= min_pos) & (df["POS"] <= max_pos)
        df.loc[mask, "in_locus"] = True

    return df, leads_df
