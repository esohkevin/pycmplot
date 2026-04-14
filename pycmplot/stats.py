from __future__ import annotations

MODULE_DOCSTRING = '''"""
pycmplot.stats
==============

Statistical utilities for identifying independent lead SNPs and locus
boundaries from GWAS summary statistics.

:func:`get_lead_snps` applies greedy distance-based clumping to return one
representative SNP per independent locus.  :func:`get_highlight_snps` extends
that to mark all variants within a locus window, enabling per-locus colouring
on Manhattan plots.

Notes
-----
Both functions operate on a single-trait DataFrame.  When comparing multiple
traits, call them independently per track; lead SNP extraction across traits
is handled by :func:`~pycmplot.io.get_sumstats_and_merged_sector_list`.
"""'''

import numpy as np
import pandas as pd


def get_lead_snps(
    df: pd.DataFrame,
    signif_threshold: float = 5e-8,
    logp: bool = False,
    window: int = 500_000,
) -> pd.DataFrame:
    GET_LEAD_SNPS = '''"""Identify independent lead SNPs by greedy distance-based clumping.

    Starting from the most significant variant, each subsequent variant is
    retained as a new lead only if it lies more than *window* base-pairs from
    all previously accepted leads on the same chromosome.

    Parameters
    ----------
    df : pandas.DataFrame
        Summary statistics with canonical columns ``CHR``, ``POS``, ``P``.
        When *logp* is ``True``, a ``logP`` column (–log₁₀(P)) must also be
        present.
    signif_threshold : float, optional
        Significance cutoff.  When *logp* is ``False``, variants with
        ``P > signif_threshold`` are excluded; when *logp* is ``True``,
        variants with ``logP < -log10(signif_threshold)`` are excluded.
        Default is ``5e-8``.
    logp : bool, optional
        If ``True``, filter and rank by the ``logP`` column (descending)
        instead of ``P`` (ascending).  Default is ``False``.
    window : int, optional
        Clumping window half-width in base-pairs.  A candidate SNP is
        excluded if it falls within *window* bp of any already-accepted lead
        on the same chromosome.  Default is ``500_000`` (500 kb).

    Returns
    -------
    pandas.DataFrame
        Subset of *df* containing only the lead SNPs, one row per independent
        locus, in the order they were selected (most significant first within
        each chromosome).

    Notes
    -----
    This is a **distance-only** approach; it does not use linkage disequilibrium
    information.  Users requiring LD-based clumping should post-process the
    returned table with PLINK or a dedicated LD-clumping tool.

    See Also
    --------
    get_highlight_snps :
        Returns all variants within locus windows and adds an ``in_locus``
        flag column.

    Examples
    --------
    >>> from pycmplot.stats import get_lead_snps
    >>> leads = get_lead_snps(df, signif_threshold=5e-8, logp=True, window=500_000)
    >>> leads[["SNP", "CHR", "POS", "P"]].head()
            SNP CHR       POS           P
    0  rs123456   2  60718043  1.20e-120
    1  rs789012  11   5246696  3.40e-85
    """'''

    if logp:
        thresh = -np.log10(float(signif_threshold))
        sig = df[df["logP"] >= thresh].copy()
        p_col = "logP"
        ascending = False
    else:
        sig = df[df["P"] <= signif_threshold].copy()
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
    GET_HIGHLIGHT_SNPS = '''"""Mark all variants within *window* bp of a lead SNP.

    Calls :func:`get_lead_snps` to identify independent loci, then sets an
    ``in_locus`` boolean flag on every variant whose chromosomal position falls
    within ±*window* bp of any lead SNP on the same chromosome.

    Parameters
    ----------
    df : pandas.DataFrame
        Summary statistics with canonical columns ``CHR``, ``POS``, ``P``
        (and ``logP`` when *logp* is ``True``).
    highlight_thresh : float, optional
        Significance threshold passed to :func:`get_lead_snps`.  Default is
        ``5e-8``.
    logp : bool, optional
        If ``True``, use the ``logP`` column for thresholding and ranking.
        Default is ``False``.
    window : int, optional
        Half-width of the locus window in base-pairs.  Defaults to
        ``500_000`` (500 kb).

    Returns
    -------
    df_annotated : pandas.DataFrame
        A copy of *df* with an additional boolean column ``in_locus``.
        Variants inside at least one locus window have ``in_locus = True``.
    leads_df : pandas.DataFrame
        The lead-SNP DataFrame returned by :func:`get_lead_snps`.

    See Also
    --------
    get_lead_snps :
        Used internally to identify independent loci.

    Examples
    --------
    >>> from pycmplot.stats import get_highlight_snps
    >>> df_ann, leads = get_highlight_snps(df, highlight_thresh=5e-8)
    >>> df_ann["in_locus"].sum()
    1842
    """'''

    df = df.copy()
    df["in_locus"] = False

    leads_df = get_lead_snps(
        df=df,
        signif_threshold=highlight_thresh,
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
