MODULE_DOCSTRING = '''"""
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
"""'''

from __future__ import annotations

import bisect
import logging
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
    BUILD_GENES_DICT = '''"""Build a chromosome-keyed interval dictionary with sorted start positions.

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
    """'''

    genes_df = genes_df.sort_values(["CHR", "START"])
    genes_dict: dict = {}

    for chrom, group in genes_df.groupby("CHR"):
        intervals = list(
            zip(
                group["START"].astype(int),
                group["END"].astype(int),
                group["STRAND"],
                group["GENE"],
            )
        )
        starts = [g[0] for g in intervals]
        genes_dict[str(chrom)] = {"intervals": intervals, "starts": starts}

    return genes_dict


# ---------------------------------------------------------------------------
# Internal: strand-aware variant annotation
# ---------------------------------------------------------------------------

def _annotate_variant(
    chrom: str,
    pos: int,
    genes_dict: dict,
    window: int = 500_000,
    promoter_window: int = 2_000,
) -> dict:
    ANNOTATE_VARIANT = '''"""Return strand-aware nearest-gene annotation for a single variant.

    Searches the pre-built *genes_dict* within *window* bp of *pos* on
    *chrom*.  Reports the nearest upstream and downstream genes (relative to
    each gene's TSS and strand), whether the variant falls inside a gene body,
    and whether it is within a promoter window upstream of any gene.

    Parameters
    ----------
    chrom : str
        Chromosome (without ``'chr'`` prefix, e.g. ``'11'``, ``'X'``).
    pos : int
        Variant position in base-pairs (hg38 1-based coordinates).
    genes_dict : dict
        Pre-built chromosome-keyed interval dictionary from
        :func:`_build_genes_dict`.
    window : int, optional
        Search radius in base-pairs.  Default is ``500_000`` (500 kb).
    promoter_window : int, optional
        Distance upstream of the TSS considered a promoter region.
        Default is ``2_000`` (2 kb).

    Returns
    -------
    dict
        Keys:

        * ``genic`` (bool) – ``True`` if the variant overlaps a gene body.
        * ``nearest_upstream_gene`` (str or None) – symbol of the closest
        gene whose TSS is downstream of the variant (relative to strand).
        * ``upstream_distance`` (int or None) – bp distance to
        ``nearest_upstream_gene``; ``0`` when genic.
        * ``nearest_downstream_gene`` (str or None) – symbol of the closest
        gene whose TSS is upstream of the variant (relative to strand).
        * ``downstream_distance`` (int or None) – bp distance to
        ``nearest_downstream_gene``.
        * ``promoter_upstream_flag`` (bool) – ``True`` when the variant is
        within *promoter_window* bp upstream of any TSS.
        * ``gene_density`` (int) – number of genes with any overlap in the
        search window.
    """'''

    _empty = {
        "genic": False,
        "nearest_upstream_gene": None,
        "upstream_distance": None,
        "nearest_downstream_gene": None,
        "downstream_distance": None,
        "promoter_upstream_flag": False,
        "bidirectional_promoter_flag": False,
        "gene_density": 0,
    }

    if chrom not in genes_dict:
        return _empty

    chrom_data = genes_dict[chrom]
    genes = chrom_data["intervals"]
    starts = chrom_data["starts"]

    left_bound = pos - window
    right_bound = pos + window

    i = bisect.bisect_left(starts, left_bound)

    gene_density = 0
    nearest_upstream: Optional[str] = None
    nearest_downstream: Optional[str] = None
    min_up_dist = float("inf")
    min_down_dist = float("inf")
    promoter_upstream_flag = False

    while i < len(genes):
        start, end, strand, gene = genes[i]

        if start > right_bound:
            break

        if end >= left_bound:
            gene_density += 1

            if start <= pos <= end:
                return {
                    "genic": True,
                    "nearest_upstream_gene": gene,
                    "upstream_distance": 0,
                    "nearest_downstream_gene": None,
                    "downstream_distance": None,
                    "promoter_upstream_flag": False,
                    "gene_density": gene_density,
                }

            tss = start if strand == "+" else end
            distance = abs(pos - tss)

            if distance <= window:
                if strand == "+":
                    is_upstream = pos < tss
                    in_promoter = (tss - promoter_window) <= pos < tss
                else:
                    is_upstream = pos > tss
                    in_promoter = tss < pos <= (tss + promoter_window)

                if is_upstream:
                    if distance < min_up_dist:
                        min_up_dist = distance
                        nearest_upstream = gene
                    if in_promoter:
                        promoter_upstream_flag = True
                else:
                    if distance < min_down_dist:
                        min_down_dist = distance
                        nearest_downstream = gene

        i += 1

    return {
        "genic": False,
        "nearest_upstream_gene": nearest_upstream,
        "upstream_distance": min_up_dist if nearest_upstream else None,
        "nearest_downstream_gene": nearest_downstream,
        "downstream_distance": min_down_dist if nearest_downstream else None,
        "promoter_upstream_flag": promoter_upstream_flag,
        "gene_density": gene_density,
    }


# ---------------------------------------------------------------------------
# Internal: prioritisation scorer
# ---------------------------------------------------------------------------

def _annotate_and_prioritize_variant(
    chrom: str,
    pos: int,
    genes_df: pd.DataFrame,
    lead_snps_df: pd.DataFrame,
    window: int = 500_000,
    promoter_window: int = 2_000,
    biotype_weights: Optional[dict] = None,
) -> Optional[dict]:
    ANNOTATE_PRIORITIZE = '''"""Score and rank candidate genes for a single variant using a composite
    priority metric.

    Builds a candidate gene set within *window* bp of *pos* on *chrom*, then
    scores each candidate on four additive components:

    * **Genic overlap** (weight 2) – variant falls inside the gene body.
    * **Promoter proximity** (weight 1) – variant is within *promoter_window*
    bp upstream of the TSS (strand-aware).
    * **Biotype weight** (from *biotype_weights*, scaled by distance score) –
    penalises pseudogenes and non-coding features relative to protein-coding
    genes.
    * **Distance score** = 1 / log₁₀(distance + 10) – continuously rewards
    closeness.

    The top-ranked gene (or the two closest intergenic flanking genes joined by
    ``'-'``) is returned.

    Parameters
    ----------
    chrom : str
        Chromosome string (no ``'chr'`` prefix).
    pos : int
        Variant position in base-pairs.
    genes_df : pandas.DataFrame
        Full gene reference DataFrame with columns ``CHR``, ``START``,
        ``END``, ``STRAND``, ``GENE``, ``BIOTYPE``.
    lead_snps_df : pandas.DataFrame
        Lead-SNP DataFrame (currently passed for context; reserved for future
        co-localisation scoring).
    window : int, optional
        Search radius in base-pairs.  Default is ``500_000``.
    promoter_window : int, optional
        Promoter window in bp upstream of TSS.  Default is ``2_000``.
    biotype_weights : dict, optional
        Mapping of Ensembl biotype → numeric weight.  Defaults to
        :data:`~pycmplot.constants.BIOTYPE_WEIGHTS`.

    Returns
    -------
    dict or None
        Keys: ``top_gene``, ``biotype``, ``priority_score``, ``distance``,
        ``promoter_flag``, ``distance_score``, ``biotype_weight``,
        ``promoter_bonus``, ``gene_density``.  Returns ``None`` if *chrom* has
        no gene entries in *genes_df*.

        For intergenic variants, ``top_gene`` contains the two nearest flanking
        gene symbols joined by ``'-'`` (e.g. ``'HBB-HBD'``) and ``biotype``
        is set to ``'intergenic'``.
    """'''

    if biotype_weights is None:
        biotype_weights = BIOTYPE_WEIGHTS

    genes_df = genes_df.copy()
    genes_df["TSS"] = np.where(
        genes_df["STRAND"] == "+",
        genes_df["START"],
        genes_df["END"],
    )

    chr_genes = genes_df[genes_df["CHR"] == chrom]
    if chr_genes.empty:
        return None

    candidates = chr_genes[
        (chr_genes["START"] <= pos + window) & (chr_genes["END"] >= pos - window)
    ].copy()

    if candidates.empty:
        return None

    gene_density = len(candidates)

    candidates["distance"] = np.where(
        (pos >= candidates["START"]) & (pos <= candidates["END"]),
        0,
        np.minimum(
            abs(pos - candidates["START"]),
            abs(pos - candidates["END"]),
        ),
    )

    candidates["genic"] = (pos >= candidates["START"]) & (pos <= candidates["END"])

    candidates["promoter_flag"] = (
        (candidates["STRAND"] == "+")
        & (pos >= candidates["TSS"] - promoter_window)
        & (pos <= candidates["TSS"])
    ) | (
        (candidates["STRAND"] == "-")
        & (pos <= candidates["TSS"] + promoter_window)
        & (pos >= candidates["TSS"])
    )

    candidates["distance_score"] = 1 / np.log10(candidates["distance"] + 10)
    candidates["biotype_weight"] = candidates["BIOTYPE"].map(
        lambda x: biotype_weights.get(x, 0)
    )
    candidates["promoter_bonus"] = candidates["promoter_flag"].astype(int) * 0.5
    candidates["priority_score"] = (
        candidates["genic"].astype(int) * 2
        + candidates["promoter_flag"].astype(int) * 1
        + candidates["biotype_weight"] * 2 * candidates["distance_score"]
    )

    candidates = candidates.sort_values("priority_score", ascending=False)

    if candidates.empty:
        return {
            "top_gene": None, "biotype": None, "priority_score": None,
            "distance": None, "promoter_flag": None, "distance_score": None,
            "biotype_weight": None, "promoter_bonus": None, "gene_density": None,
        }

    if candidates["genic"].any():
        top = candidates.iloc[0]
        return {
            "top_gene": top["GENE"],
            "biotype": top["BIOTYPE"],
            "priority_score": top["priority_score"],
            "distance": top["distance"],
            "promoter_flag": top["promoter_flag"],
            "distance_score": top["distance_score"],
            "biotype_weight": top["biotype_weight"],
            "promoter_bonus": top["promoter_bonus"],
            "gene_density": gene_density,
        }
    else:
        top2 = candidates.head(2)
        return {
            "top_gene": "-".join(top2["GENE"]),
            "biotype": "intergenic",
            "priority_score": None,
            "distance": "-".join(map(str, top2["distance"])),
            "promoter_flag": None,
            "distance_score": None,
            "biotype_weight": None,
            "promoter_bonus": None,
            "gene_density": None,
        }


# ---------------------------------------------------------------------------
# Internal: clumping
# ---------------------------------------------------------------------------

def _clump_by_distance(df: pd.DataFrame, window_kb: int = 500) -> pd.DataFrame:
    CLUMP_BY_DISTANCE = '''"""Reduce a lead-SNP table to one representative SNP per locus.

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
    """'''

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
    GET_HITS_SUMMARY_TABLE = '''"""Annotate lead SNPs with nearest genes and write the locus summary table.

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

        .. list-table::
        :widths: 30 70
        :header-rows: 1

        * - Column
            - Description
        * - ``genic``
            - ``True`` when the lead SNP overlaps a gene body
        * - ``nearest_upstream_gene``
            - Nearest upstream gene symbol (strand-aware)
        * - ``upstream_distance``
            - Distance to ``nearest_upstream_gene`` in bp
        * - ``nearest_downstream_gene``
            - Nearest downstream gene symbol (strand-aware)
        * - ``downstream_distance``
            - Distance to ``nearest_downstream_gene`` in bp
        * - ``promoter_upstream_flag``
            - ``True`` when the SNP is within 2 kb upstream of a TSS
        * - ``gene_density``
            - Number of genes within the search window
        * - ``top_gene``
            - Top-priority gene from the scoring pass
        * - ``biotype``
            - Ensembl biotype of ``top_gene`` (``'intergenic'`` when no
            genic overlap)
        * - ``priority_score``
            - Composite priority score (genic hits only)

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
    """'''

    if resources is None:
        resources = default_resources

    # Choose gene info file based on build
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
        annotation = _annotate_variant(
            chrom=row["CHR"],
            pos=row["POS"],
            genes_dict=genes_dict,
            window=window,
        )
        prioritized = _annotate_and_prioritize_variant(
            chrom=row["CHR"],
            pos=row["POS"],
            genes_df=geneinfo,
            lead_snps_df=leads_df,
            window=window,
        )

        record = {
            **(row.to_dict()),
            **(annotation if annotation is not None else {}),
            **(prioritized if prioritized is not None else {}),
        }
        records.append(record)

    locus_table = pd.DataFrame(records).sort_values(
        ["CHR", "POS"], key=natsort.natsort_keygen()
    )

    if table_out is not None:
        locus_table.to_csv(table_out, index=False, sep="\t", na_rep="None")
        logger.info("Locus summary written to: %s", table_out)

    return _clump_by_distance(locus_table, window_kb=window_kb)
