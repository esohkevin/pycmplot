"""
pycmplot.annotation
====================

Nearest-gene annotation for GWAS lead SNPs and generation of the structured
locus summary table.

The main public function, :func:`get_hits_summary_table`, accepts the lead
SNP DataFrame produced by :func:`~pycmplot.stats.get_lead_snps`, annotates
each lead with the nearest (and most biologically plausible) gene using a
two-pass strategy — strand-aware boundary distance followed by a composite
priority score — and writes a tab-delimited locus summary file alongside the
plot.

Gene reference files
--------------------
Annotation relies on a bundled Ensembl gene-info TSV (hg38 or hg19).  The
file is resolved through :class:`~pycmplot.resources.ResourceConfig`; custom
paths can be supplied via the ``PYCMPLOT_GENEINFO_HG38`` /
``PYCMPLOT_GENEINFO_HG19`` environment variables.
"""

from __future__ import annotations

import bisect
import logging
import math
from typing import Optional

import natsort
import numpy as np
import pandas as pd

from pycmplot.constants import BIOTYPE_WEIGHTS
from pycmplot.resources import ResourceConfig, default_resources

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal: gene dictionary builder
# ---------------------------------------------------------------------------

def _build_genes_dict(genes_df: pd.DataFrame) -> dict:
    """Build a chromosome-keyed interval dictionary with sorted start positions.

    Pre-processes the gene reference DataFrame into a structure that supports
    efficient O(log N) binary-search lookup of genes near a query position.

    Parameters
    ----------
    genes_df : pandas.DataFrame
        Gene reference with columns ``CHR``, ``START``, ``END``,
        ``STRAND``, ``GENE``.  The ``START`` and ``END`` columns must be
        coercible to ``int``.

    Returns
    -------
    dict
        Mapping of ``chromosome_string → {'intervals': [...], 'starts': [...]}``.

        * **intervals** – list of ``(start, end, strand, gene_symbol)`` tuples
        sorted by ``start`` position.
        * **starts** – flat list of ``start`` values, used as the sorted key
        sequence for :func:`bisect.bisect_left`.

    Notes
    -----
    This function is called once per :func:`get_hits_summary_table` invocation;
    the result is passed to :func:`_annotate_variant` for each lead SNP.
    """

    genes_df = genes_df.sort_values(["CHR", "START"])
    genes_dict: dict = {}
    # Include BIOTYPE in the tuple so the merged
    # :func:`_annotate_variant` (which now covers both nearest-gene
    # detection *and* biotype-weighted prioritisation in one pass)
    # doesn't have to touch the source DataFrame again.  Falls back to
    # ``None`` when the reference lacks a BIOTYPE column so older gene
    # references still work.
    has_biotype = "BIOTYPE" in genes_df.columns
    for chrom, group in genes_df.groupby("CHR"):
        if has_biotype:
            intervals = list(
                zip(
                    group["START"].astype(int),
                    group["END"].astype(int),
                    group["STRAND"],
                    group["GENE"],
                    group["BIOTYPE"],
                )
            )
        else:
            intervals = list(
                zip(
                    group["START"].astype(int),
                    group["END"].astype(int),
                    group["STRAND"],
                    group["GENE"],
                    [None] * len(group),
                )
            )
        starts = [g[0] for g in intervals]
        genes_dict[str(chrom)] = {"intervals": intervals, "starts": starts}

    return genes_dict


# ---------------------------------------------------------------------------
# Internal: single-pass variant annotation
# ---------------------------------------------------------------------------

def _annotate_variant(
    chrom: str,
    pos: int,
    genes_dict: dict,
    window: int = 500_000,
    promoter_window: int = 2_000,
    biotype_weights: Optional[dict] = None,
) -> dict:
    """Single-pass positional + biotype-weighted gene annotation for a lead SNP.

    Walks every candidate gene inside ``[pos - window, pos + window]``
    on *chrom* exactly once, and returns the union of the fields
    previously produced by two separate passes
    (``_annotate_variant`` and ``_annotate_and_prioritize_variant``,
    since merged here).  Consumers of the hits summary table therefore
    see the same columns as before, but the loop body runs half as many
    times and one source of truth governs the distance / genic /
    promoter semantics.

    The merger is safe because the two former passes computed
    everything from the same candidate set, only differing in what
    they *reported*.  Now they report jointly.

    Parameters
    ----------
    chrom : str
        Chromosome (without ``'chr'`` prefix, e.g. ``'11'``, ``'X'``).
    pos : int
        Variant position in base-pairs (1-based coordinates).
    genes_dict : dict
        Pre-built chromosome-keyed interval dictionary from
        :func:`_build_genes_dict`.  Tuples are now
        ``(START, END, STRAND, GENE, BIOTYPE)``.
    window : int, optional
        Search radius in base-pairs.  Default is ``500_000`` (500 kb).
    promoter_window : int, optional
        Distance upstream of the TSS considered a promoter region.
        Default is ``2_000`` (2 kb).
    biotype_weights : dict, optional
        Mapping of Ensembl biotype → numeric weight.  Defaults to
        :data:`~pycmplot.constants.BIOTYPE_WEIGHTS`.

    Returns
    -------
    dict
        Positional fields (from the pre-merge ``_annotate_variant``):

        * ``genic`` (bool)
        * ``nearest_gene`` / ``nearest_gene_distance``
        * ``nearest_upstream_gene`` / ``upstream_distance`` (positional
          left flanker: gene body ends at a lower coordinate than *pos*)
        * ``nearest_downstream_gene`` / ``downstream_distance``
          (positional right flanker)
        * ``promoter_upstream_flag`` (strand-aware — retained because
          "promoter" is genuinely a biological concept)
        * ``gene_density``

        Prioritisation fields (from the pre-merge
        ``_annotate_and_prioritize_variant``):

        * ``top_gene`` — highest-priority gene when genic, otherwise
          the ``"LEFT_GENE-RIGHT_GENE"`` positional flanker pair in
          genomic order (or a single symbol when only one side has a
          gene in the window).  Flanker selection uses raw
          base-pair distance, not priority score, so the joined label
          always brackets the variant.
        * ``biotype`` — Ensembl biotype of ``top_gene`` (``'intergenic'``
          when no genic overlap).
        * ``priority_score`` — composite score ``2·genic + 1·promoter +
          2·biotype_w · dist_score`` (genic hits only).
        * ``distance``, ``promoter_flag``, ``distance_score``,
          ``biotype_weight``, ``promoter_bonus`` — components of the
          winning candidate's score.
    """
    if biotype_weights is None:
        biotype_weights = BIOTYPE_WEIGHTS

    _empty = {
        "genic": False,
        "nearest_gene": None,
        "nearest_gene_distance": None,
        "nearest_upstream_gene": None,
        "upstream_distance": None,
        "nearest_downstream_gene": None,
        "downstream_distance": None,
        "promoter_upstream_flag": False,
        "bidirectional_promoter_flag": False,
        "gene_density": 0,
        "top_gene": None,
        "biotype": None,
        "priority_score": None,
        "distance": None,
        "promoter_flag": None,
        "distance_score": None,
        "biotype_weight": None,
        "promoter_bonus": None,
    }
    if chrom not in genes_dict:
        return _empty

    chrom_data = genes_dict[chrom]
    genes = chrom_data["intervals"]
    starts = chrom_data["starts"]
    left_bound = pos - window
    right_bound = pos + window

    # Advance to the first gene whose START could still overlap the
    # search window; walk forward until STARTs exceed the right bound.
    i = bisect.bisect_left(starts, left_bound)

    gene_density = 0
    containing_gene: Optional[str] = None
    containing_biotype: Optional[str] = None

    # Positional trackers — closest gene entirely to the left of *pos*
    # (by distance to END), closest entirely to the right (by distance
    # to START), and the true closest gene overall (any side, by
    # distance to gene body).
    nearest_left: Optional[str] = None
    nearest_left_dist = float("inf")
    nearest_right: Optional[str] = None
    nearest_right_dist = float("inf")
    nearest_any: Optional[str] = None
    nearest_any_dist = float("inf")

    promoter_upstream_flag = False

    # Priority-score tracker.  Keeps the winner's *full* component
    # breakdown so the returned record has the same shape as the
    # pre-merge output — this is what preserves numerical identity
    # against the two-pass version.
    top_score = -float("inf")
    top: Optional[dict] = None

    while i < len(genes):
        _entry = genes[i]
        # Genes_dict tuples are now (start, end, strand, gene, biotype).
        # Accept a 4-tuple too for backwards compat with cached
        # references built by an older pycmplot version.
        if len(_entry) >= 5:
            start, end, strand, gene, biotype = _entry[:5]
        else:
            start, end, strand, gene = _entry[:4]
            biotype = None

        if start > right_bound:
            break
        if end < left_bound:
            i += 1
            continue

        gene_density += 1
        is_genic = start <= pos <= end

        # Distance to gene body — 0 when inside, else min gap to either
        # edge.  Also update the positional flanker slots + the
        # "true nearest" slot.
        if is_genic:
            if containing_gene is None:
                # If multiple genes contain pos (nested / overlapping
                # bodies), the first one wins the ``nearest_gene``
                # slot; the higher-priority one still wins ``top_gene``
                # via the score tracker below.
                containing_gene = gene
                containing_biotype = biotype
            distance = 0
            if nearest_any_dist > 0:
                nearest_any_dist = 0
                nearest_any = gene
        else:
            if end < pos:
                distance = pos - end
                if distance < nearest_left_dist:
                    nearest_left_dist = distance
                    nearest_left = gene
            else:
                distance = start - pos
                if distance < nearest_right_dist:
                    nearest_right_dist = distance
                    nearest_right = gene
            if distance < nearest_any_dist:
                nearest_any_dist = distance
                nearest_any = gene

        # Strand-aware promoter check.  TSS = START for '+' genes,
        # END for '-'.  This is the one place strand-awareness is
        # preserved because "promoter" is genuinely a biological
        # concept — every other field is positional.
        tss = start if strand == "+" else end
        if strand == "+":
            in_promoter = (tss - promoter_window) <= pos < tss
        else:
            in_promoter = tss < pos <= (tss + promoter_window)
        if in_promoter:
            promoter_upstream_flag = True

        # Composite priority score — identical formula to the pre-merge
        # ``_annotate_and_prioritize_variant``.  ``distance_score`` uses
        # ``log10(distance + 10)`` so it stays finite at distance 0.
        distance_score = 1.0 / math.log10(distance + 10)
        biotype_weight = biotype_weights.get(biotype, 0) if biotype else 0
        promoter_bonus = 0.5 * (1 if in_promoter else 0)
        priority_score = (
            (2 if is_genic else 0)
            + (1 if in_promoter else 0)
            + biotype_weight * 2 * distance_score
        )
        if priority_score > top_score:
            top_score = priority_score
            top = {
                "gene": gene,
                "biotype": biotype,
                "distance": distance,
                "promoter_flag": in_promoter,
                "distance_score": distance_score,
                "biotype_weight": biotype_weight,
                "promoter_bonus": promoter_bonus,
                "priority_score": priority_score,
                "is_genic": is_genic,
            }

        i += 1

    is_genic_any = containing_gene is not None

    # ---- Assemble positional fields ------------------------------------
    result = {
        "genic": is_genic_any,
        "nearest_gene": containing_gene if is_genic_any else nearest_any,
        "nearest_gene_distance": (
            0 if is_genic_any else (
                int(nearest_any_dist) if nearest_any is not None else None
            )
        ),
        "nearest_upstream_gene": nearest_left,
        "upstream_distance": (
            int(nearest_left_dist) if nearest_left is not None else None
        ),
        "nearest_downstream_gene": nearest_right,
        "downstream_distance": (
            int(nearest_right_dist) if nearest_right is not None else None
        ),
        "promoter_upstream_flag": promoter_upstream_flag,
        "gene_density": gene_density,
    }

    # ---- Assemble prioritisation fields --------------------------------
    if top is None:
        # No candidates at all — populate with None.
        result.update({
            "top_gene": None, "biotype": None, "priority_score": None,
            "distance": None, "promoter_flag": None,
            "distance_score": None, "biotype_weight": None,
            "promoter_bonus": None,
        })
    elif is_genic_any:
        # Genic branch: report the highest-priority gene overall.
        # (Almost always this is the containing gene because the
        # ``genic·2`` term dominates, but the pre-merge code took the
        # global argmax so we preserve that exactly.)
        result.update({
            "top_gene": top["gene"],
            "biotype": top["biotype"],
            "priority_score": top["priority_score"],
            "distance": top["distance"],
            "promoter_flag": top["promoter_flag"],
            "distance_score": top["distance_score"],
            "biotype_weight": top["biotype_weight"],
            "promoter_bonus": top["promoter_bonus"],
        })
    else:
        # Intergenic branch: ``top_gene`` is the ``"LEFT-RIGHT"``
        # positional flanker pair by base-pair distance (not priority
        # score).  This is the guarantee introduced by the 2026-09-05
        # fix — the joined label always brackets the variant.
        parts: list[str] = []
        dist_parts: list[str] = []
        if nearest_left is not None:
            parts.append(str(nearest_left))
            dist_parts.append(str(int(nearest_left_dist)))
        if nearest_right is not None:
            parts.append(str(nearest_right))
            dist_parts.append(str(int(nearest_right_dist)))
        if not parts:
            # Everything in the window somehow failed both flanker
            # checks (shouldn't happen because is_genic_any is False
            # and we saw at least one candidate); fall back to the
            # top-priority pick so the field is never empty.
            parts = [str(top["gene"])]
            dist_parts = [str(int(top["distance"]))]
        result.update({
            "top_gene": "-".join(parts),
            "biotype": "intergenic",
            "priority_score": None,
            "distance": "-".join(dist_parts),
            "promoter_flag": None,
            "distance_score": None,
            "biotype_weight": None,
            "promoter_bonus": None,
        })

    return result




# ---------------------------------------------------------------------------
# Internal: clumping
# ---------------------------------------------------------------------------

def _clump_by_distance(df: pd.DataFrame, window_kb: int = 500) -> pd.DataFrame:
    """Reduce a lead-SNP table to one representative SNP per locus.

    Applies greedy distance-based clumping within each chromosome group,
    starting from the most significant SNP (lowest ``P`` or highest ``logP``).
    Candidate SNPs within *window_kb* kilobases of an already-accepted lead are
    discarded.

    Parameters
    ----------
    df : pandas.DataFrame
        Lead-SNP DataFrame with columns ``CHR``, ``POS``, and either ``P`` or
        ``logP``.
    window_kb : int, optional
        Clumping window half-width in kilobases.  Default is ``500``.

    Returns
    -------
    pandas.DataFrame
        Deduplicated locus representatives sorted by chromosome and position
        (natural sort order).
    """

    window = window_kb * 1000
    clumped: list[pd.Series] = []

    for _chrom, group in df.groupby("CHR"):
        if "logP" in df.columns:
            group = group.sort_values("logP", ascending=False)
        else:
            group = group.sort_values("P", ascending=True)

        kept_positions: list[int] = []
        for _, row in group.iterrows():
            if all(abs(row["POS"] - p) > window for p in kept_positions):
                clumped.append(row)
                kept_positions.append(row["POS"])

    return pd.DataFrame(clumped).sort_values(
        ["CHR", "POS"], key=natsort.natsort_keygen()
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def get_hits_summary_table(
    leads_df: pd.DataFrame,
    window_kb: int = 500,
    table_out: Optional[str] = None,
    resources: Optional[ResourceConfig] = None,
) -> pd.DataFrame:
    """Annotate lead SNPs with nearest genes and write the locus summary table.

    For each lead SNP in *leads_df*, runs two complementary annotation passes:

    1. **Strand-aware boundary search** (:func:`_annotate_variant`) — identifies
    the nearest upstream and downstream genes and detects genic / promoter
    overlap.
    2. **Priority scoring** (:func:`_annotate_and_prioritize_variant`) — ranks
    all candidate genes within *window_kb* by a composite score that weights
    biotype, promoter proximity, and distance, then selects the single
    top-ranked gene (or the two flanking genes for intergenic hits).

    After annotation, the table is deduplicated with distance-based clumping
    (:func:`_clump_by_distance`) and optionally written to *table_out*.

    Parameters
    ----------
    leads_df : pandas.DataFrame
        DataFrame of lead SNPs as returned by
        :func:`~pycmplot.stats.get_lead_snps`.  Must contain columns
        ``CHR``, ``POS``, ``P``, ``BUILD``.
    window_kb : int, optional
        Search radius in kilobases around each lead SNP.  Default is ``500``.
    table_out : str or None, optional
        File path at which to write the annotated locus summary table as a
        tab-delimited TSV.  Set to ``None`` to suppress file output.
    resources : ResourceConfig, optional
        :class:`~pycmplot.resources.ResourceConfig` instance providing paths to
        the Ensembl gene-info TSV (hg38 or hg19).  Defaults to
        :data:`~pycmplot.resources.default_resources`.

    Returns
    -------
    pandas.DataFrame
        Clumped locus summary table.  Contains all columns from *leads_df*
        plus annotation fields from both passes, including:

        - ``genic`` — ``True`` when the lead SNP overlaps a gene body.
        - ``nearest_gene`` — closest gene by base-pair distance
          regardless of side (the containing gene when genic, else the
          nearest flanking gene).
        - ``nearest_gene_distance`` — bp distance to ``nearest_gene``;
          0 when genic.
        - ``nearest_upstream_gene`` — closest gene positionally to the
          left of the SNP (lower genomic coordinate).  As of pycmplot
          0.4.x this is a positional definition, not strand-aware; see
          :func:`_annotate_variant` for the migration note.
        - ``upstream_distance`` — distance to ``nearest_upstream_gene``
          in bp.
        - ``nearest_downstream_gene`` — closest gene positionally to
          the right of the SNP.
        - ``downstream_distance`` — distance to ``nearest_downstream_gene``
          in bp.
        - ``promoter_upstream_flag`` — ``True`` when the SNP is within
          2 kb upstream of a TSS.  This is the only field that remains
          strand-aware.
        - ``gene_density`` — number of genes within the search window.
        - ``top_gene`` — top-priority gene from the scoring pass.  For
          intergenic hits, ``"LEFT-RIGHT"`` where LEFT is the nearest
          gene positionally to the left of the SNP and RIGHT is the
          nearest gene positionally to the right (order is always
          genomic).  When only one side has a gene in the window, a
          single symbol is returned.
        - ``biotype`` — Ensembl biotype of ``top_gene`` (``'intergenic'`` when
          no genic overlap).
        - ``priority_score`` — composite priority score (genic hits only).

    Notes
    -----
    The gene reference (hg38 or hg19) is selected automatically based on the
    ``BUILD`` column in *leads_df*.  hg19 builds are matched to the GRCh37
    gene-info file; all others use the GRCh38 file.

    See Also
    --------
    pycmplot.stats.get_lead_snps :
        Provides the *leads_df* input to this function.
    pycmplot.resources.ResourceConfig :
        Controls the paths to the gene-info reference files.

    Examples
    --------
    >>> from pycmplot.annotation import get_hits_summary_table
    >>> hits = get_hits_summary_table(
    ...     leads_df=leads,
    ...     window_kb=500,
    ...     table_out="./results/HbF_locus_summary.tsv",
    ... )
    >>> hits[["SNP", "CHR", "POS", "top_gene", "biotype"]].head()
            SNP CHR       POS  top_gene           biotype
    0  rs123456   2  60718043    BCL11A    protein_coding
    1  rs789012  11   5246696       HBB    protein_coding
    """

    if resources is None:
        resources = default_resources

    # Choose gene info file based on build
    if 'BUILD' in leads_df.columns:
        if "OLD_POS" not in leads_df.columns and list(set(leads_df["BUILD"])) == ["hg19"]:
            geneinfo_path = resources.require("geneinfo_hg19")
        else:
            geneinfo_path = resources.require("geneinfo_hg38")

        logger.info("Loading gene info from: %s", geneinfo_path)
        geneinfo = pd.read_csv(geneinfo_path, header=0, sep="\t")
        genes_dict = _build_genes_dict(geneinfo)

        window = window_kb * 1_000
        records: list[dict] = []


        logger.info("Annotating lead variants and generating hits summary table ...")
        for _, row in leads_df.iterrows():
            # ``_annotate_variant`` now emits both positional
            # nearest-gene fields *and* biotype-weighted prioritisation
            # fields in a single window walk (the two pre-merge passes
            # were consolidated in 0.4.x for a ~2× speedup on the
            # annotation step and one source of truth for the shared
            # distance / genic / promoter semantics).
            annotation = _annotate_variant(
                chrom=row["CHR"],
                pos=row["POS"],
                genes_dict=genes_dict,
                window=window,
            )

            record = {
                **(row.to_dict()),
                **(annotation if annotation is not None else {}),
            }
            records.append(record)

        locus_table = pd.DataFrame(records).sort_values(
            ["CHR", "POS"], key=natsort.natsort_keygen()
        )
    else:
        locus_table = leads_df

    if table_out is not None:
        outpath = table_out.replace(" ", "_").lower() + '.tsv'
        locus_table.to_csv(outpath, index=False, sep="\t", na_rep="None")
        logger.info("Locus summary written to: %s", outpath)

    return _clump_by_distance(locus_table, window_kb=window_kb)


def get_annotation_column(
    annotate: str = None, 
    hits_table: pd.DataFrame = None,
    label_col: str = None,
):
    label_clm = 'SNP'
    if annotate is not None and not hits_table.empty:
        if label_col is not None and label_col in hits_table.columns:
            label_clm = label_col
        elif annotate in hits_table.columns:
            label_clm = annotate
        else:
            if str(annotate).upper() == "GENE":
                for i, (_, row) in enumerate(hits_table.iterrows()):
                    try:
                        if row["genic"]:
                            # Prefer ``nearest_gene`` (the containing gene
                            # for genic variants).  Older cached hits
                            # tables predate this column, so we fall back
                            # to ``nearest_upstream_gene`` when it's
                            # missing.
                            label_clm = ("nearest_gene"
                                         if "nearest_gene" in hits_table.columns
                                         else "nearest_upstream_gene")
                            label_msg = (
                                f"Signal {row['SNP']} at {row['POS']} is genic "
                                f"[{row.get(label_clm)}]"
                            )
                        else:
                            label_clm = "top_gene"
                            label_msg = "'POS' is not genic"
                        #logger.info("%s", label_msg)
                    except Exception:
                        logger.warning(
                            "Annotation columns '%s' and '%s' not found in hits table: %s; "
                            "falling back to 'SNP'.", annotate, label_col, hits_table.columns.values,
                        )
                        #label_clm = 'SNP'
              
    logger.info("Annotating by: %s", label_clm)

    return label_clm


# ---------------------------------------------------------------------------
# Per-locus highlight color resolution
# ---------------------------------------------------------------------------

HIGHLIGHT_COLOR_COL = "highlight_color"
HIGHLIGHT_COLOR_AUTO = "auto"

CATEGORY_COL = "category"
CATEGORY_DEFAULT = "significant"


def ensure_highlight_color_column(
    hits_table: "pd.DataFrame",
    default: str = HIGHLIGHT_COLOR_AUTO,
) -> "pd.DataFrame":
    """Make sure the hits table has a ``highlight_color`` column.

    Backward-compat helper.  When cached hits overlays predate the
    per-locus color feature, the column may be missing entirely; older
    Python-API callers may build their own hits tables without it.  We
    inject the column with the sentinel *default* so the downstream
    color-lookup path can uniformly assume the column exists.  Existing
    values (including user-supplied hex codes / names / ``"auto"``) are
    preserved untouched.
    """
    import pandas as pd
    if hits_table is None or not isinstance(hits_table, pd.DataFrame):
        return hits_table
    if HIGHLIGHT_COLOR_COL not in hits_table.columns:
        hits_table = hits_table.copy()
        hits_table[HIGHLIGHT_COLOR_COL] = default
    return hits_table


def ensure_category_column(
    hits_table: "pd.DataFrame",
    default: str = CATEGORY_DEFAULT,
) -> "pd.DataFrame":
    """Make sure the hits table has a ``category`` column.

    Companion to :func:`ensure_highlight_color_column`; injects the
    sentinel default so downstream code can assume the column exists.
    """
    import pandas as pd
    if hits_table is None or not isinstance(hits_table, pd.DataFrame):
        return hits_table
    if CATEGORY_COL not in hits_table.columns:
        hits_table = hits_table.copy()
        hits_table[CATEGORY_COL] = default
    return hits_table


def resolve_highlight_colors(
    sig_df: "pd.DataFrame",
    hits_table: "pd.DataFrame",
    default_color: str,
    *,
    chr_col: str = "CHR",
    pos_col: str = "POS",
    window_kb: int = 500,
) -> list[str]:
    """Return a per-row highlight color for the highlighted variants.

    For every row in *sig_df* (the ``in_locus`` subset of a track), we
    find the closest lead in *hits_table* on the same chromosome and,
    within *window_kb*, use its ``highlight_color`` value.  The sentinel
    ``"auto"``, empty strings, NaNs, and colors matplotlib can't parse
    all fall back to *default_color*.  Invalid colors emit a
    ``logger.warning`` naming the offending value so users can fix a
    typo in their overlay TSV.

    Parameters
    ----------
    sig_df
        Subset of a track's DataFrame with ``in_locus == True``.  Must
        have *chr_col* and *pos_col*.
    hits_table
        Locus summary table (per-lead).  Must have *chr_col* and
        *pos_col*.  If it doesn't have ``highlight_color``, every row
        falls back to *default_color* (backward compat).
    default_color
        Global highlight color passed to the plotter.  Used when the
        cached row's color is ``"auto"``, missing, or invalid.
    window_kb
        Search radius in kb for matching a highlighted variant to its
        lead.  Defaults to 500 kb, which mirrors the highlight
        window in :func:`~pycmplot.stats.get_highlight_snps`.

    Returns
    -------
    list of str
        One matplotlib-parseable color per row of *sig_df*, in row
        order.
    """
    import numpy as np
    import pandas as pd
    from matplotlib.colors import is_color_like

    n = len(sig_df.index)
    if n == 0:
        return []

    # Fast path: no hits table (or missing color column) → uniform default.
    if hits_table is None or hits_table.empty \
            or HIGHLIGHT_COLOR_COL not in hits_table.columns:
        return [default_color] * n

    window_bp = int(window_kb) * 1000

    # Group leads by chromosome once so per-row lookup is cheap.
    hits_by_chr: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for chrom, sub in hits_table.groupby(chr_col, observed=True):
        try:
            pos_arr = sub[pos_col].astype("int64").to_numpy()
        except (TypeError, ValueError):
            continue
        col_arr = sub[HIGHLIGHT_COLOR_COL].astype(str).to_numpy()
        hits_by_chr[str(chrom)] = (pos_arr, col_arr)

    _warned: set = set()

    def _fallback_for(reason: str, offending: str) -> str:
        key = (reason, offending)
        if key not in _warned:
            _warned.add(key)
            logger.warning(
                "highlight_color %r rejected (%s); using default %r.",
                offending, reason, default_color,
            )
        return default_color

    colors_out: list[str] = []
    for chrom, pos in zip(
        sig_df[chr_col].astype(str).to_numpy(),
        sig_df[pos_col].astype("int64").to_numpy(),
    ):
        entry = hits_by_chr.get(str(chrom))
        if entry is None:
            colors_out.append(default_color)
            continue
        lead_pos, lead_col = entry
        d = np.abs(lead_pos - int(pos))
        j = int(d.argmin())
        if d[j] > window_bp:
            colors_out.append(default_color)
            continue
        raw = lead_col[j].strip() if isinstance(lead_col[j], str) else ""
        if not raw or raw.lower() == HIGHLIGHT_COLOR_AUTO \
                or raw.lower() in {"nan", "none", "na"}:
            colors_out.append(default_color)
            continue
        if not is_color_like(raw):
            colors_out.append(_fallback_for("not a matplotlib color", raw))
            continue
        colors_out.append(raw)
    return colors_out


# ---------------------------------------------------------------------------
# Category resolution + legend entry construction
# ---------------------------------------------------------------------------

def resolve_highlight_categories(
    sig_df: "pd.DataFrame",
    hits_table: "pd.DataFrame",
    default_category: str = CATEGORY_DEFAULT,
    *,
    chr_col: str = "CHR",
    pos_col: str = "POS",
    window_kb: int = 500,
) -> list:
    """Return a per-row category label for the highlighted variants.

    Same lookup discipline as :func:`resolve_highlight_colors`: each
    ``sig_df`` row is matched to its nearest lead on the same
    chromosome within *window_kb*, and the lead's ``category`` value
    is returned.  Missing / blank / NaN values fall back to
    *default_category* (typically ``"significant"``).
    """
    import numpy as np
    import pandas as pd

    n = len(sig_df.index)
    if n == 0:
        return []
    if hits_table is None or hits_table.empty \
            or CATEGORY_COL not in hits_table.columns:
        return [default_category] * n

    window_bp = int(window_kb) * 1000
    hits_by_chr = {}
    for chrom, sub in hits_table.groupby(chr_col, observed=True):
        try:
            pos_arr = sub[pos_col].astype("int64").to_numpy()
        except (TypeError, ValueError):
            continue
        cat_arr = sub[CATEGORY_COL].astype(str).to_numpy()
        hits_by_chr[str(chrom)] = (pos_arr, cat_arr)

    cats_out = []
    for chrom, pos in zip(
        sig_df[chr_col].astype(str).to_numpy(),
        sig_df[pos_col].astype("int64").to_numpy(),
    ):
        entry = hits_by_chr.get(str(chrom))
        if entry is None:
            cats_out.append(default_category)
            continue
        lead_pos, lead_cat = entry
        d = np.abs(lead_pos - int(pos))
        j = int(d.argmin())
        if d[j] > window_bp:
            cats_out.append(default_category)
            continue
        raw = lead_cat[j].strip() if isinstance(lead_cat[j], str) else ""
        if not raw or raw.lower() in {"nan", "none", "na"}:
            cats_out.append(default_category)
            continue
        cats_out.append(raw)
    return cats_out


def build_highlight_legend_entries(
    hits_table: "pd.DataFrame",
    default_color: str,
    default_category: str = CATEGORY_DEFAULT,
) -> list:
    """Return a list of ``(category, color)`` pairs for a custom legend.

    Only returns entries when the user has customised the overlay in a
    way worth surfacing:

    * Any row has a category other than *default_category*, OR
    * Any row has a highlight_color other than the sentinel
      ``"auto"`` (i.e. an explicit user-chosen colour).

    When nothing has been customised, returns ``[]`` -- the plotter
    then skips the legend entirely, matching the pre-feature layout.

    Entries are deduplicated by ``category`` in first-appearance order
    (so users can control legend order by reordering rows in their
    TSV).  If the same category appears with multiple colours, the
    first colour wins and a warning is logged.
    """
    from matplotlib.colors import is_color_like

    if hits_table is None or hits_table.empty:
        return []

    has_category = CATEGORY_COL in hits_table.columns
    has_color = HIGHLIGHT_COLOR_COL in hits_table.columns
    if not has_category and not has_color:
        return []

    def _norm(x, default):
        if x is None:
            return default
        s = str(x).strip()
        if not s or s.lower() in {"nan", "none", "na", HIGHLIGHT_COLOR_AUTO}:
            return default
        return s

    any_custom_cat = False
    any_custom_color = False
    if has_category:
        any_custom_cat = any(
            _norm(v, default_category) != default_category
            for v in hits_table[CATEGORY_COL]
        )
    if has_color:
        any_custom_color = any(
            _norm(v, "auto") != "auto"
            for v in hits_table[HIGHLIGHT_COLOR_COL]
        )
    if not any_custom_cat and not any_custom_color:
        return []

    seen = {}
    _warned = set()
    for _, row in hits_table.iterrows():
        cat = _norm(row.get(CATEGORY_COL) if has_category else None,
                    default_category)
        col = _norm(row.get(HIGHLIGHT_COLOR_COL) if has_color else None,
                    default_color)
        if not is_color_like(col):
            col = default_color
        if cat not in seen:
            seen[cat] = col
        elif seen[cat] != col and (cat, col) not in _warned:
            _warned.add((cat, col))
            logger.warning(
                "Category %r appears with multiple colours; using first "
                "(%r) for the legend and ignoring %r.",
                cat, seen[cat], col,
            )
    return list(seen.items())


# ---------------------------------------------------------------------------
# Legend placement helper (shared by the linear + circular plotters)
# ---------------------------------------------------------------------------

# Matplotlib's standard 9-way loc grid, exposed as a docstring / help
# hint so users don't have to look it up.  Anything outside this set is
# still accepted by matplotlib (it will raise a clear error) but these
# are the recommended values.
HIGHLIGHT_LEGEND_LOCS = (
    "upper center", "upper left", "upper right",
    "center", "center left", "center right",
    "lower center", "lower left", "lower right",
    "best",
)


def resolve_highlight_legend_placement(
    loc,
    default: str = "upper center",
) -> dict:
    """Translate a user-facing ``highlight_legend_loc`` value into legend kwargs.

    Accepts either:

    * A matplotlib ``loc`` string (``"upper center"``, ``"upper right"``,
      ``"lower center"``, ``"best"``, …).  Returned as ``{"loc": loc}``.
    * A 2-tuple ``(x, y)`` in axes coordinates (``0.0-1.0``, but values
      outside are accepted for anchoring outside the axes — useful for
      the circular plotter where placing the legend *below* the polar
      frame at e.g. ``(0.5, -0.05)`` avoids overlapping the sectors).
      Returned as ``{"loc": "center", "bbox_to_anchor": (x, y)}``.
    * ``None``, in which case *default* (``"upper center"``) is used.

    Anything else is coerced to a string and passed through as ``loc`` —
    that way matplotlib emits its own clear error when the user typos
    something like ``"upper centre"``.
    """
    if loc is None:
        return {"loc": default}
    # 2-tuple (or list) → bbox anchoring at that point.
    if isinstance(loc, (tuple, list)) and len(loc) == 2:
        try:
            x = float(loc[0]); y = float(loc[1])
        except (TypeError, ValueError):
            pass
        else:
            return {"loc": "center", "bbox_to_anchor": (x, y)}
    return {"loc": str(loc)}
