"""
pycmplot.annotation
===================
Nearest-gene annotation and locus summary table generation.
"""

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
    """Build a chromosome-keyed interval dict with sorted start positions.

    Parameters
    ----------
    genes_df:
        DataFrame with columns ``CHR``, ``START``, ``END``, ``STRAND``, ``GENE``.

    Returns
    -------
    dict keyed by chromosome string; each value is
    ``{"intervals": [...], "starts": [...]}``.
    """
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
    """Return strand-aware nearest-gene annotation for a single variant.

    Returns a dict with keys:
    ``genic``, ``nearest_upstream_gene``, ``upstream_distance``,
    ``nearest_downstream_gene``, ``downstream_distance``,
    ``promoter_upstream_flag``, ``gene_density``.
    """
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
    """Annotate lead SNPs with nearest genes and write a summary table.

    Parameters
    ----------
    leads_df:
        DataFrame of lead SNPs (output of :func:`~pycmplot.stats.get_lead_snps`).
        Must contain columns ``CHR``, ``POS``, ``P``, ``BUILD``.
    window_kb:
        Window in kb around each lead SNP to search for genes (default 500 kb).
    table_out:
        If provided, write the clumped table to this TSV file path.
    resources:
        :class:`~pycmplot.resources.ResourceConfig` instance.

    Returns
    -------
    pd.DataFrame
        Clumped locus summary table with gene annotations.
    """
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
