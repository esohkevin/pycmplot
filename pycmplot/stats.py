"""
pycmplot.stats
==============

Statistical utilities for identifying independent lead SNPs and locus
boundaries from GWAS summary statistics.

:func:`get_lead_snps` applies greedy distance-based clumping to return one
representative SNP per independent locus.  :func:`get_highlight_snps` extends
that to mark all variants within a locus window, enabling per-locus colouring
on Manhattan plots.  :func:`clump` applies LD-based clumping using an
external PLINK ``.ld`` reference panel — the LD-aware equivalent of
:func:`get_lead_snps`.

Notes
-----
Both functions operate on a single-trait DataFrame.  When comparing multiple
traits, call them independently per track; lead SNP extraction across traits
is handled by :func:`~pycmplot.io.get_sumstats_and_merged_sector_list`.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import Optional, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def get_lead_snps(
    df: pd.DataFrame,
    signif_threshold: float = 5e-8,
    logp: bool = False,
    window: int = 250_000,
    score_col: Optional[str] = None,
    ascending: Optional[bool] = None,
    tiebreak_score_col: Optional[str] = None,
) -> pd.DataFrame:
    """Identify independent lead SNPs by greedy distance-based clumping.

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
    """

    # Resolve column + comparison direction.  Callers can override the
    # default (P-value semantics) to point at any column — e.g.
    # ``score_col="P_UNSIGNED", ascending=False`` for signed selection
    # statistics like iHS or XP-EHH, where "more significant" means larger
    # |value| rather than smaller p-value.
    if score_col is None and ascending is None:
        # Legacy calling convention — preserved verbatim.
        if logp:
            thresh = -np.log10(float(signif_threshold))
            sig = df[df["logP"] >= thresh].copy()
            score_col = "logP"
            ascending = False
        else:
            sig = df[df["P"] <= float(signif_threshold)].copy()
            score_col = "P"
            ascending = True
    else:
        if score_col is None:
            score_col = "logP" if logp else "P"
        if ascending is None:
            ascending = not logp
        thresh = float(signif_threshold)
        if ascending:
            sig = df[df[score_col] <= thresh].copy()
        else:
            sig = df[df[score_col] >= thresh].copy()

    # Sort so that (a) the most-significant variant is first (primary
    # key = score_col), and (b) among tied-P variants the highest-
    # tiebreak variant is first (secondary key = tiebreak_score_col
    # descending).  When the caller supplies a per-variant
    # informativeness column (e.g. GH / eQTL / cCRE presence),
    # tied-P LD-block members get greedy-picked by function rather
    # than by POS accident.
    if tiebreak_score_col is not None and tiebreak_score_col in sig.columns:
        sig = sig.sort_values(
            [score_col, tiebreak_score_col],
            ascending=[ascending, False],
            kind="mergesort",
        )
    else:
        sig = sig.sort_values(score_col, ascending=ascending, kind="mergesort")
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


_SIG_KEEP_COLS = ("CHR", "POS", "SNP", "P", "logP",
                  "OLD_POS", "OLD_BUILD", "BUILD", "LABEL")


def get_signif_snps(
    df: pd.DataFrame,
    signif_threshold: float = 5e-8,
    logp: bool = False,
    score_col: Optional[str] = None,
    ascending: Optional[bool] = None,
    keep_cols: Optional[list[str]] = None,
) -> pd.DataFrame:
    """Extract every variant crossing *signif_threshold*.

    A thin subset of :func:`get_lead_snps`'s pre-clumping step: filter
    the rows and return a minimal DataFrame with only the columns
    downstream annotation / clumping steps need.  No greedy clumping,
    no LD lookup, no priority scoring — just the sig subset.

    Motivating the split: :func:`pycmplot.io.load` used to run
    ``get_highlight_snps`` (which includes clumping) per track and
    then concatenate.  On multi-track runs with an LD reference the
    reference was reloaded per track (~30 s for a 31 M-pair PLINK
    file), which dominated wall time.  Now the loader calls
    :func:`get_signif_snps` per track (essentially free) and defers
    clumping / prioritisation to a single post-loop pass on the
    concatenated set — LD graph + functional tracks are loaded once
    across all tracks.

    Parameters
    ----------
    df : pandas.DataFrame
        Summary statistics with canonical columns ``CHR``, ``POS``,
        ``P`` (and ``logP`` when *logp* is True).
    signif_threshold : float, optional
        P (or ``|value|`` for signed statistics) cutoff.  Default 5e-8.
    logp, score_col, ascending :
        Same semantics as :func:`get_lead_snps` — used to decide the
        comparison direction (ascending for p-values, descending for
        logP or |value|).
    keep_cols : list[str], optional
        Column names to retain in the output.  Defaults to the
        pycmplot pipeline's canonical set: ``CHR``, ``POS``, ``SNP``,
        ``P``, ``logP``, ``OLD_POS``, ``OLD_BUILD``, ``BUILD``,
        ``LABEL``.  Missing columns are silently skipped.

    Returns
    -------
    pandas.DataFrame
        Sig subset with only *keep_cols* present.  Order preserved
        from the input.
    """
    if score_col is None and ascending is None:
        # Legacy calling convention preserved.
        if logp:
            thresh = -np.log10(float(signif_threshold))
            sig = df[df["logP"] >= thresh]
        else:
            sig = df[df["P"] <= float(signif_threshold)]
    else:
        if score_col is None:
            score_col = "logP" if logp else "P"
        if ascending is None:
            ascending = not logp
        thresh = float(signif_threshold)
        if ascending:
            sig = df[df[score_col] <= thresh]
        else:
            sig = df[df[score_col] >= thresh]

    #if keep_cols is None:
    #    keep_cols = list(_SIG_KEEP_COLS)
    #present = [c for c in keep_cols if c in sig.columns]
    #return sig[present].copy()
    return sig

def get_highlight_snps(
    df: pd.DataFrame,
    highlight: bool = False,
    highlight_thresh: float = 5e-8,
    logp: bool = False,
    window: int = 250_000,
    score_col: Optional[str] = None,
    ascending: Optional[bool] = None,
    tiebreak_score_col: Optional[str] = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Mark all variants within *window* bp of a lead SNP.

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
    """

    df = df.copy()

    sig_df = get_signif_snps(
        df=df,
        signif_threshold=highlight_thresh,
        logp=logp,
        score_col=score_col,
        ascending=ascending,
    )

    #leads_df = get_lead_snps(
    #    df=df,
    #    signif_threshold=highlight_thresh,
    #    logp=logp,
    #    window=window,
    #    score_col=score_col,
    #    ascending=ascending,
    #    tiebreak_score_col=tiebreak_score_col,
    #)

    if highlight:
        df["in_locus"] = False
        for _, row in sig_df.iterrows():
            min_pos = row["POS"] - window
            max_pos = row["POS"] + window
            chrom = row["CHR"]

            mask = (df["CHR"] == chrom) & (df["POS"] >= min_pos) & (df["POS"] <= max_pos)
            df.loc[mask, "in_locus"] = True

    return df, sig_df


# ---------------------------------------------------------------------------
# LD-based clumping (PLINK-style)
# ---------------------------------------------------------------------------

def _parse_plink_ld(path: str) -> dict:
    """Parse a PLINK ``--r2`` (``.ld``) file into a symmetric r² dict.

    Expected columns (space or tab separated, PLINK's default): at least
    ``SNP_A``, ``SNP_B``, ``R2``.  ``CHR_A`` / ``BP_A`` / ``CHR_B`` /
    ``BP_B`` are ignored (positions come from the summary-stats
    DataFrame).

    Returns ``{snp_a: {snp_b: r2, ...}, snp_b: {snp_a: r2, ...}}``.
    Symmetry is stored explicitly for O(1) lookup in either
    direction.  Pairs are stored only if ``R2`` parses as a float
    in ``[0, 1]``.
    """
    logger.info("Parsing PLINK LD reference: %s", path)
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    required = {"SNP_A", "SNP_B", "R2"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"PLINK .ld file at {path} is missing required column(s) "
            f"{sorted(missing)}.  Header seen: {list(df.columns)}"
        )
    ld: dict[str, dict[str, float]] = defaultdict(dict)
    snp_a = df["SNP_A"].astype(str).to_numpy()
    snp_b = df["SNP_B"].astype(str).to_numpy()
    r2 = df["R2"].astype(float).to_numpy()
    for a, b, r in zip(snp_a, snp_b, r2):
        if not (0.0 <= r <= 1.0):
            continue
        ld[a][b] = r
        ld[b][a] = r
    logger.info(
        "  parsed %s SNP pairs; %s unique SNPs in reference panel",
        f"{len(df.index):,}", f"{len(ld):,}",
    )
    return dict(ld)


def clump(
    df: pd.DataFrame,
    ld_reference: Union[str, dict, "LDGraph", None] = None,
    r2: float = 0.1,
    kb: int = 500,
    p_threshold: float = 5e-8,
    logp: bool = False,
    snp_col: str = "SNP",
    p_col: str = "P",
    logp_col: str = "logP",
    chr_col: str = "CHR",
    pos_col: str = "POS",
    harmonize_variants: bool = False,
    variant_matcher: Optional["VariantMatcher"] = None,
) -> pd.DataFrame:
    """LD-based greedy clumping (PLINK ``--clump`` semantics).

    Starting from the most-significant variant, each subsequent variant
    is retained as an independent lead only if it is *both* more than
    ``kb`` kilobases from every already-accepted lead on the same
    chromosome AND has ``r² < r2`` with every accepted lead in the
    reference panel.

    Parameters
    ----------
    df : pandas.DataFrame
        Summary statistics with columns ``CHR``, ``POS``, ``SNP``, and
        ``P`` (or ``logP`` when ``logp=True``).
    ld_reference : str, dict, or None
        Path to a PLINK ``.ld`` file (produced by ``plink --r2``), or a
        pre-parsed ``{snp_a: {snp_b: r2}}`` dict.  ``None`` falls back
        to distance-only clumping equivalent to :func:`get_lead_snps`.
    r2 : float, optional
        LD threshold in [0, 1].  Variants with ``r² >= r2`` to any
        already-accepted lead are clumped away.  Default 0.1 (PLINK
        ``--clump-r2`` default is 0.1 in modern GWAS pipelines).
    kb : int, optional
        Distance window in kilobases.  A candidate variant more than
        ``kb`` bp from every accepted lead is retained regardless of
        r² (no LD lookup needed).  Default 500 kb.
    p_threshold : float, optional
        Only variants with ``P <= p_threshold`` (or ``logP >= -log10(p_threshold)``
        when ``logp=True``) enter the clumping pool.  Default ``5e-8``.
    logp : bool, optional
        Rank by ``logP`` (descending) instead of ``P`` (ascending).
        Default ``False``.
    snp_col, p_col, logp_col, chr_col, pos_col : str
        Column-name overrides for non-standard schemas.
    harmonize_variants : bool, optional
        Resolve naming mismatches between the sumstats and an
        :class:`~pycmplot.ld.LDGraph` reference (``rs``/``CHR:POS``/
        ``chr``-prefix/allele-order/separator differences) before
        clumping.  Only the significant hits are mapped onto the
        graph's existing node names via
        :func:`pycmplot.variant_matcher.ld_id_map`; the graph itself is
        not modified and the returned ``SNP`` column is unchanged.
        Ignored for dict / ``None`` references.  Default ``False``.
    variant_matcher : VariantMatcher, optional
        Deprecated.  Any non-``None`` value is treated as
        ``harmonize_variants=True``; the matcher's contents are not
        used.

    Returns
    -------
    pandas.DataFrame
        Subset of *df* containing only the independent lead variants,
        in the order they were selected (most-significant first within
        each chromosome).

    Notes
    -----
    * When a lead or candidate is absent from the LD reference, r²
      cannot be computed and the candidate is clumped by distance
      alone (folded into the lead if within ``kb``) rather than being
      promoted to an independent lead.  The fraction of pair checks
      that used this fallback is logged.
    * The distance window is a fast short-circuit: it lets the
      algorithm skip the LD dict lookup for variants that are
      obviously unlinked by distance alone.

    Examples
    --------
    >>> from pycmplot.stats import clump
    >>> leads = clump(sumstats, ld_reference="ref.ld",
    ...               r2=0.1, kb=500, p_threshold=5e-8, logp=True)
    """
    # ------------------------------------------------------------------
    # 1. Resolve LD reference
    # ------------------------------------------------------------------
    # Accepted forms (checked in order):
    #   * ``LDGraph`` instance          → used directly
    #   * ``.npz`` path                 → LDGraph.load
    #   * ``.ld`` / other path          → LDGraph.from_plink_ld
    #     (the sparse CSR representation is 5-10× smaller than the
    #      dict-of-dict fallback and supports ``any_linked`` batched
    #      lookup — preferred for reference panels >10 K variants)
    #   * ``dict``                      → legacy dict-of-dict
    #   * ``None``                      → distance-only clumping
    from pycmplot.ld import LDGraph  # local import — avoid module cycle

    ld: Union[dict, LDGraph]
    if ld_reference is None:
        ld = {}
    elif isinstance(ld_reference, LDGraph):
        ld = ld_reference
    elif isinstance(ld_reference, dict):
        ld = ld_reference
    elif isinstance(ld_reference, (str, bytes)) or hasattr(ld_reference, "__fspath__"):
        _p = str(ld_reference)
        if _p.endswith(".npz"):
            ld = LDGraph.load(_p)
        else:
            ld = LDGraph.from_plink_ld(_p, r2_threshold=float(r2))
    else:
        raise TypeError(
            f"ld_reference must be an LDGraph, a path, a dict, or "
            f"None; got {type(ld_reference).__name__}."
        )

    if variant_matcher is not None:
        import warnings
        warnings.warn(
            "clump(variant_matcher=...) is deprecated; pass "
            "harmonize_variants=True instead.  The matcher is rebuilt "
            "internally from the significant hits.",
            DeprecationWarning, stacklevel=2,
        )
        harmonize_variants = True

    # ------------------------------------------------------------------
    # 2. Filter to sig variants and sort by significance
    # ------------------------------------------------------------------
    if logp:
        thresh = -np.log10(float(p_threshold))
        sig = df[df[logp_col] >= thresh].copy()
        sig = sig.sort_values(logp_col, ascending=False, kind="mergesort")
    else:
        sig = df[df[p_col] <= float(p_threshold)].copy()
        sig = sig.sort_values(p_col, ascending=True, kind="mergesort")

    if sig.empty:
        return sig

    # Optional variant-ID harmonisation.  Map only the significant
    # hits onto the LD graph's existing node names; r² lookups then
    # use the mapped name while the ``SNP`` column is left as-is.  The
    # graph is never relabelled, so a shared / cached LDGraph is safe
    # to reuse and two tracks naming the same variant differently both
    # resolve to the same node.
    ld_ids = sig[snp_col].astype(str)
    if harmonize_variants and isinstance(ld, LDGraph):
        try:
            from pycmplot.variant_matcher import ld_id_map
            _m = ld_id_map(sig, ld.snp_names, ld.snp_index,
                           snp_col=snp_col, chr_col=chr_col,
                           pos_col=pos_col,
                           ld_chrom=ld.snp_chrom, ld_pos=ld.snp_pos)
            n_renamed = sum(1 for k, v in _m.items() if k != v)
            ld_ids = ld_ids.map(lambda x: _m.get(x, x))
            logger.info(
                "clump: harmonize_variants resolved %s of %s significant "
                "hits to the LD reference (%s via ID harmonisation)",
                f"{len(_m):,}", f"{sig[snp_col].nunique():,}", f"{n_renamed:,}",
            )
        except Exception as exc:  # pragma: no cover  (defensive)
            logger.warning(
                "clump: variant harmonisation failed (%s); continuing "
                "with raw IDs.", exc,
            )

    # ------------------------------------------------------------------
    # 3. Greedy clumping
    # ------------------------------------------------------------------
    # When either the lead or candidate is MISSING from the LD
    # reference, we can't compute r² — treat as clumped-by-distance
    # instead of assuming independence.  This is the biology-safe
    # default: an unreferenced variant *might* be in LD with the
    # lead, so it's conservative to fold it in rather than promote
    # it to a new independent lead.  Missing decisions are counted
    # and logged at the end so the user knows how much of the
    # clumping was actually LD-driven vs. distance-driven fallback.
    window = int(kb) * 1_000
    r2_thresh = float(r2)
    is_graph = isinstance(ld, LDGraph)
    have_ld_ref = bool(ld)
    leads: list[pd.Series] = []
    lead_ld_ids: list[str] = []
    n_pair_checks = 0
    n_fallback_missing = 0

    def _in_ref(snp: str) -> bool:
        if not have_ld_ref:
            return False
        if is_graph:
            return snp in ld
        return snp in ld

    # Positional pairing (not index lookup) so duplicate DataFrame
    # indices from a caller's concat can't cross-wire LD names.
    for (_, row), this_snp in zip(sig.iterrows(), ld_ids.tolist()):
        this_chr = row[chr_col]
        this_pos = int(row[pos_col])
        clumped = False
        for lead, lead_snp in zip(leads, lead_ld_ids):
            if lead[chr_col] != this_chr:
                continue
            if abs(int(lead[pos_col]) - this_pos) > window:
                continue
            # Same chromosome, within distance window.
            n_pair_checks += 1
            if lead_snp == this_snp:
                # Same variant (e.g. the same hit in two tracks, or two
                # names harmonised onto one LD node).  The graph stores
                # no self-edges, so r2() would return 0 and keep it as
                # a second lead; a variant is in perfect LD with itself.
                clumped = True
                break
            if have_ld_ref and _in_ref(lead_snp) and _in_ref(this_snp):
                if is_graph:
                    r = ld.r2(lead_snp, this_snp)
                else:
                    r = ld.get(lead_snp, {}).get(this_snp, 0.0)
                if r >= r2_thresh:
                    clumped = True
                    break
            else:
                # LD unknown — biology-safe fallback: within kb
                # window means "assume possibly linked" → clumped.
                n_fallback_missing += 1
                clumped = True
                break
        if not clumped:
            leads.append(row)
            lead_ld_ids.append(this_snp)

    if have_ld_ref and n_pair_checks:
        _frac_fallback = 100.0 * n_fallback_missing / n_pair_checks
        logger.info(
            "clump: %s pair checks, %s (%.1f%%) fell back to "
            "distance-only because at least one variant was absent "
            "from the LD reference",
            f"{n_pair_checks:,}", f"{n_fallback_missing:,}",
            _frac_fallback,
        )
        if _frac_fallback > 20.0:
            logger.warning(
                "clump: LD-reference overlap is low — %.1f%% of pair "
                "checks used distance-only fallback.  Consider a "
                "reference panel whose SNP naming matches your "
                "summary stats, or upstream harmonisation.",
                _frac_fallback,
            )

    return pd.DataFrame(leads)

