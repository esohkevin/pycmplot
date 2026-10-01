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


# Biotypes treated as high-confidence coding for the ``top_gene``
# preference rule.  When any of these appear among the candidates
# (per :func:`_annotate_variant`), ``top_gene`` prefers them over
# non-coding winners even if a non-coding candidate scored higher —
# encoding the biological intuition that a divergent lncRNA
# containing the SNP is almost always really about the neighboring
# PC gene.  Kept as a module-level constant so tests and future
# tuning have a single source of truth.
PC_LIKE_BIOTYPES = {
    "protein_coding",
    "protein_coding_LoF",
    "protein_coding_CDS_not_defined",
    "nonsense_mediated_decay",
    # Immunoglobulin & T-cell receptor coding genes — legitimate
    # protein-coding at the biological level even though Ensembl
    # gives them their own biotypes.
    "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene",
    "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene",
}


# ---------------------------------------------------------------------------
# GeneHancer regulatory-element loader + per-position lookup
# ---------------------------------------------------------------------------

# Divisor applied to the raw GeneHancer connection score to derive a
# priority bonus.  score/5 with no cap and no distance decay.  The
# variant must fall *inside* a GH element to trigger the bonus at
# all — biologically, a variant regulates a target only if it
# disrupts a specific enhancer, and GeneHancer's confidence score
# for that element-gene edge tells us the strength of the
# hypothesis.  No proximity or "nearby" softening; if the SNP is
# outside every element, GH contributes nothing and the algorithm
# falls back to geometry + strand + biotype.
#
# ``score=20`` -> bonus 4.0 (rare, dominates all other components);
# ``score=10`` -> bonus 2.0 (comparable to genic + upstream);
# ``score=1``  -> bonus 0.2 (nudge only).
_GENEHANCER_DIVISOR = 5.0


_LOADER_CACHE: dict = {}


def _file_key(path):
    import os
    try:
        st = os.stat(path)
        return (str(path), st.st_size, st.st_mtime_ns)
    except (OSError, TypeError):
        return None


def load_genehancer(path: Optional[str] = None) -> Optional[dict]:
    """Parse the packaged GeneHancer TSV into a chrom-keyed lookup.

    Returns a dict mapping ``chromosome_string -> {"starts": list[int],
    "ends": list[int], "rows": list[(gene, score)]}`` sorted by START.
    ``rows[i]`` gives the (gene, GH connection score) tuple for the
    element starting at ``starts[i]``.  Multiple elements at the same
    position (one GH element -> N connected genes) each get their own
    row so a per-position lookup can return every ``(gene, score)``
    the element predicts.

    Passing ``None`` (or a non-existent path) returns ``None`` — the
    caller then skips GH-based scoring entirely.
    """
    import os

    if path is None or not os.path.exists(path):
        return None
    _k = ("genehancer", _file_key(path))
    if _k in _LOADER_CACHE:          # shared by the tie-break and scoring
        return _LOADER_CACHE[_k]

    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return None
    df["CHR"] = df["CHR"].astype(str)
    df["START"] = df["START"].astype(int)
    df["END"] = df["END"].astype(int)
    df["GENE"] = df["GENE"].astype(str)
    df["SCORE"] = df["SCORE"].astype(float)

    out: dict = {}
    for chrom, group in df.groupby("CHR"):
        group = group.sort_values("START")
        out[str(chrom)] = {
            "starts": group["START"].to_list(),
            "ends":   group["END"].to_list(),
            "rows":   list(zip(group["GENE"], group["SCORE"])),
        }
    _LOADER_CACHE[_k] = out
    return out


def _genehancer_bonuses(chrom: str, pos: int,
                        gh_dict: Optional[dict]) -> dict:
    """Return ``{gene: bonus}`` for every gene the GH element(s)
    containing *pos* link to.

    Strict-overlap scoring: only elements where ``start <= pos <= end``
    contribute.  Within each containing element, every connected
    gene gets ``bonus = score / 5`` (uncapped).  Multiple containing
    elements linking the same gene collapse to the ``max`` bonus.

    Returns an empty dict when *pos* falls outside every element on
    *chrom*, or when *gh_dict* is ``None``.
    """
    if gh_dict is None:
        return {}
    entry = gh_dict.get(str(chrom))
    if entry is None:
        return {}
    starts = entry["starts"]
    ends = entry["ends"]
    rows = entry["rows"]

    # Elements are sorted by START.  ``bisect_right(starts, pos)`` gives
    # the count of elements with START <= pos; walk from 0 up to that
    # index and skip any whose END < pos.  The condition
    # ``start <= pos <= end`` then identifies exactly the elements
    # containing *pos*.
    idx = bisect.bisect_right(starts, pos)
    bonuses: dict = {}
    for i in range(idx):
        if ends[i] < pos:
            continue
        gene, score = rows[i]
        bonus = score / _GENEHANCER_DIVISOR
        if bonus > bonuses.get(gene, 0.0):
            bonuses[gene] = bonus
    return bonuses


# ---------------------------------------------------------------------------
# UCSC functional-annotation loaders + per-position lookups
# ---------------------------------------------------------------------------

# eQTL bonus scaling.  A CAVIAR CPP of 1.0 (100% causal posterior)
# gives a bonus of 3.0 — larger than any GH bonus at score 15,
# because a fine-mapped eQTL is the *strongest* possible piece of
# evidence that a variant regulates a specific gene.
_EQTL_SCALE = 3.0


def _load_interval_track(path: Optional[str],
                         cols: list[str]) -> Optional[dict]:
    """Generic chrom-keyed interval loader for cCRE / CpG / DNase /
    TFBS tracks.  Returns a dict mapping ``chrom -> {starts, ends, rows}``
    with ``starts`` sorted for :func:`bisect` lookup.  Passing ``None``
    or a missing path returns ``None`` so downstream code can skip
    the track silently.
    """
    import os

    if path is None or not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return None
    df["CHR"] = df["CHR"].astype(str)
    for c in ("START", "END"):
        df[c] = df[c].astype(int)
    out: dict = {}
    for chrom, group in df.groupby("CHR"):
        group = group.sort_values("START")
        out[str(chrom)] = {
            "starts": group["START"].to_list(),
            "ends":   group["END"].to_list(),
            "rows":   [tuple(row) for row in group[cols].itertuples(index=False, name=None)],
        }
    return out


def load_ccre(path: Optional[str]) -> Optional[dict]:
    """ENCODE cCRE lookup.  ``rows[i]`` = ``(label, score)``."""
    return _load_interval_track(path, cols=["LABEL", "SCORE"])


def load_cpg(path: Optional[str]) -> Optional[dict]:
    """CpG island lookup.  ``rows[i]`` = ``()`` (position-only track)."""
    import os
    if path is None or not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return None
    df["CHR"] = df["CHR"].astype(str)
    df["START"] = df["START"].astype(int)
    df["END"] = df["END"].astype(int)
    out: dict = {}
    for chrom, group in df.groupby("CHR"):
        group = group.sort_values("START")
        out[str(chrom)] = {
            "starts": group["START"].to_list(),
            "ends":   group["END"].to_list(),
            "rows":   [() for _ in range(len(group))],
        }
    return out


def load_dnase(path: Optional[str]) -> Optional[dict]:
    """ENCODE DNase HS cluster lookup.  ``rows[i]`` = ``(score,)``."""
    return _load_interval_track(path, cols=["SCORE"])


def load_tfbs(path: Optional[str]) -> Optional[dict]:
    """ENCODE TFBS cluster lookup.  ``rows[i]`` = ``(tf, score)``."""
    return _load_interval_track(path, cols=["TF", "SCORE"])


def load_eqtl(path: Optional[str]) -> Optional[dict]:
    """GTEx CAVIAR eQTL lookup keyed by (chrom, pos).  Returns dict
    mapping ``chrom -> {pos -> [(gene, cpp), ...]}``.
    """
    import os
    if path is None or not os.path.exists(path):
        return None
    _k = ("eqtl", _file_key(path))
    if _k in _LOADER_CACHE:
        return _LOADER_CACHE[_k]
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return None
    df["CHR"] = df["CHR"].astype(str)
    df["POS"] = df["POS"].astype(int)
    out: dict[str, dict[int, list[tuple[str, float]]]] = {}
    for chrom, group in df.groupby("CHR"):
        chrom_d: dict[int, list[tuple[str, float]]] = {}
        for pos, gene, cpp in zip(group["POS"], group["GENE"], group["CPP"]):
            chrom_d.setdefault(int(pos), []).append((str(gene), float(cpp)))
        out[str(chrom)] = chrom_d
    _LOADER_CACHE[_k] = out
    return out


def _interval_hit(chrom: str, pos: int,
                  track: Optional[dict]) -> Optional[tuple]:
    """Return the ``rows[i]`` payload of the first interval on *chrom*
    that contains *pos*, or ``None`` when no overlap exists.  Uses
    ``bisect`` on the pre-sorted starts list.
    """
    if track is None:
        return None
    entry = track.get(str(chrom))
    if entry is None:
        return None
    starts = entry["starts"]
    ends = entry["ends"]
    idx = bisect.bisect_right(starts, pos)
    for i in range(idx):
        if ends[i] >= pos:
            return entry["rows"][i]
    return None


def _eqtl_bonuses(chrom: str, pos: int,
                  eqtl_dict: Optional[dict]) -> dict:
    """Return ``{gene: eqtl_bonus}`` for a fine-mapped eQTL at *pos*.

    ``bonus = max(CPP) * _EQTL_SCALE`` per gene.  Empty dict when
    the SNP isn't a CAVIAR eQTL on this chromosome.
    """
    if eqtl_dict is None:
        return {}
    chrom_d = eqtl_dict.get(str(chrom))
    if chrom_d is None:
        return {}
    hits = chrom_d.get(int(pos))
    if not hits:
        return {}
    out: dict = {}
    for gene, cpp in hits:
        bonus = cpp * _EQTL_SCALE
        if bonus > out.get(gene, 0.0):
            out[gene] = bonus
    return out


def position_informativeness_score(chrom: str, pos: int, *,
                                   pool_genes: Optional[set] = None,
                                   genehancer: Optional[dict] = None,
                                   ccre: Optional[dict] = None,
                                   cpg: Optional[dict] = None,
                                   eqtl: Optional[dict] = None,
                                   dnase: Optional[dict] = None,
                                   tfbs: Optional[dict] = None,
                                   functional_tracks: Optional[dict] = None) -> float:
    """Compute a functional-evidence score for a single genomic
    position, used to pick the best-annotated variant within an LD
    block.  Higher = more likely to be a functional variant.

    Position-agnostic across candidate genes; the score reflects
    "how much regulatory machinery does this base pair sit on".
    The per-gene score (eqtl_bonus, genehancer_bonus) drives which
    candidate wins the ``top_gene`` slot, but this function decides
    which SNP position within a tied LD block gets annotated at all.

    *functional_tracks* is the output of :func:`load_functional_tracks`;
    when given, the cCRE / CpG / DNase / TFBS checks use it (with each
    file's coordinate convention) in place of the *ccre* / *cpg* /
    *dnase* / *tfbs* interval dicts.
    """
    score = 0.0
    # Per-gene evidence tracks (GH / eQTL): the SNP counts if it
    # links to any gene, and doubly so if that gene is a PC in the
    # candidate pool.
    if eqtl is not None:
        eb = _eqtl_bonuses(chrom, pos, eqtl)
        if eb:
            if pool_genes and (set(eb) & pool_genes):
                score += 3.0
            else:
                score += 1.5
    if genehancer is not None:
        gb = _genehancer_bonuses(chrom, pos, genehancer)
        if gb:
            if pool_genes and (set(gb) & pool_genes):
                score += 2.0
            else:
                score += 1.0
    # Position-only tracks contribute a smaller nudge each.
    if functional_tracks is not None:
        for name, weight in (("ccre", 1.5), ("cpg", 1.0), ("dnase", 1.0), ("tfbs", 0.5)):
            kind = _FUNCTIONAL_TRACKS[name][2]
            if _overlaps(functional_tracks.get(name), chrom, int(pos), kind) is not None:
                score += weight
        return score
    if ccre is not None and _interval_hit(chrom, pos, ccre) is not None:
        score += 1.5
    if cpg is not None and _interval_hit(chrom, pos, cpg) is not None:
        score += 1.0
    if dnase is not None and _interval_hit(chrom, pos, dnase) is not None:
        score += 1.0
    if tfbs is not None and _interval_hit(chrom, pos, tfbs) is not None:
        score += 0.5
    return score


# ---------------------------------------------------------------------------
# Functional-annotation columns for the hits table (reporting only)
# ---------------------------------------------------------------------------

#: Columns added to the hits table by :func:`annotate_functional`.
FUNCTIONAL_COLUMNS = (
    "gh_id", "gh_feature", "eqtl", "tfbs",
    "ccre_id", "ccre", "ccre_class", "dnase_score", "dnase_sources",
    "cpg_island",
)

# Which bundled file each reported track comes from, the columns to
# read, and how to test overlap for a 1-based variant position:
# GeneHancer is 1-based inclusive (START <= pos <= END); UCSC tracks are
# BED (0-based start, START < pos <= END); eQTLs match POS exactly.
_FUNCTIONAL_TRACKS = {
    "genehancer": ("genehancer_hg38", ["CHR", "START", "END", "FEATURE", "GH_ID"], "closed"),
    "eqtl":       ("eqtl_tissues_hg38", ["CHR", "POS", "GENE", "TISSUE", "CPP"], "point"),
    "tfbs":       ("tfbs_hg38", ["CHR", "START", "END", "TF", "SCORE"], "bed"),
    "ccre":       ("ccre_hg38", ["CHR", "START", "END", "CCRE_ID", "DESCRIPTION", "CCRE_CLASS"], "bed"),
    "dnase":      ("dnase_hg38", ["CHR", "START", "END", "SCORE", "SOURCES"], "bed"),
    "cpg":        ("cpg_hg38", ["CHR", "START", "END", "NAME"], "bed"),
}


_FUNCTIONAL_CACHE: dict = {}


def load_functional_tracks(resources: Optional[ResourceConfig] = None) -> dict:
    """Load the functional tracks used for hits-table reporting.

    Returns ``{track: {chrom: DataFrame}}`` for every track whose file
    exists and carries the reporting columns (files built by
    ``scripts/prep_functional_tracks.py``).  Tracks that are missing, or
    are older builds without those columns, are skipped with a log
    message; their hits-table columns are then left empty.
    """
    import os
    r = resources or default_resources
    # Cached per process, keyed on each file's path, size and mtime, so
    # the LD-block tie-break and the hits-table reporting share one load
    # (and a swapped file is re-read).
    _key = []
    for _attr, _c, _k in _FUNCTIONAL_TRACKS.values():
        _p = getattr(r, _attr, None)
        try:
            _st = os.stat(_p) if _p else None
            _key.append((_p, _st.st_size, _st.st_mtime_ns) if _st else (_p,))
        except OSError:
            _key.append((_p,))
    _key = tuple(_key)
    if _key in _FUNCTIONAL_CACHE:
        return _FUNCTIONAL_CACHE[_key]
    out: dict = {}
    for name, (attr, cols, _kind) in _FUNCTIONAL_TRACKS.items():
        path = getattr(r, attr, None)
        if path is None or not os.path.exists(path):
            logger.info("Functional track %r not found (%s); its columns stay empty.", name, path)
            continue
        header = pd.read_csv(path, sep="\t", nrows=0).columns
        missing = [c for c in cols if c not in header]
        if missing:
            logger.warning(
                "Functional track %r (%s) lacks columns %s; rebuild it with "
                "scripts/prep_functional_tracks.py. Its columns stay empty.",
                name, path, missing,
            )
            continue
        df = pd.read_csv(path, sep="\t", usecols=cols, dtype={"CHR": str})
        out[name] = {c: g.reset_index(drop=True) for c, g in df.groupby("CHR", sort=False)}
    _FUNCTIONAL_CACHE.clear()          # keep only the latest set in memory
    _FUNCTIONAL_CACHE[_key] = out
    return out


def _overlaps(track: Optional[dict], chrom: str, pos: int, kind: str) -> Optional[pd.DataFrame]:
    if track is None:
        return None
    g = track.get(str(chrom))
    if g is None:
        return None
    if kind == "point":
        hit = g[g["POS"].to_numpy() == pos]
    elif kind == "closed":
        hit = g[(g["START"].to_numpy() <= pos) & (g["END"].to_numpy() >= pos)]
    else:  # BED
        hit = g[(g["START"].to_numpy() < pos) & (g["END"].to_numpy() >= pos)]
    return hit if len(hit.index) else None


def _join(values) -> Optional[str]:
    vals = [str(v) for v in values if v is not None and str(v) != "nan"]
    return ",".join(dict.fromkeys(vals)) if vals else None


def annotate_functional(chrom, pos, tracks: dict) -> dict:
    """Functional-annotation fields for one variant (1-based *pos*).

    Multi-valued fields are comma-separated, strongest first:

    * ``gh_id`` / ``gh_feature`` — GeneHancer element(s) containing the
      variant and their type (e.g. ``GH01F011652``, ``Enhancer``).
    * ``eqtl`` — GTEx CAVIAR fine-mapped eQTLs at this variant as
      ``GENE_TISSUE`` (e.g. ``HBG2_Whole_Blood``), highest CPP first.
    * ``tfbs`` — ENCODE TF clusters (score >= 500) as ``TF_SCORE``.
    * ``ccre_id`` / ``ccre`` / ``ccre_class`` — ENCODE cCRE accession,
      signature description (e.g. ``proximal enhancer-like signature``)
      and class (e.g. ``pELS,CTCF-bound``).
    * ``dnase_score`` / ``dnase_sources`` — DNase cluster score
      (>= 250) and number of contributing experiments.
    * ``cpg_island`` — CpG island name (e.g. ``CpG:361``).

    Fields are ``None`` when the variant overlaps nothing in that track.
    """
    out = dict.fromkeys(FUNCTIONAL_COLUMNS)
    pos = int(pos)
    kinds = {k: v[2] for k, v in _FUNCTIONAL_TRACKS.items()}

    gh = _overlaps(tracks.get("genehancer"), chrom, pos, kinds["genehancer"])
    if gh is not None:
        gh = gh.drop_duplicates("GH_ID")
        out["gh_id"] = _join(gh["GH_ID"])
        out["gh_feature"] = _join(gh["FEATURE"])

    eq = _overlaps(tracks.get("eqtl"), chrom, pos, kinds["eqtl"])
    if eq is not None:
        # Highest CPP first; exact ties in gene, then tissue, order.
        eq = eq.sort_values(["CPP", "GENE", "TISSUE"], ascending=[False, True, True], kind="mergesort")
        out["eqtl"] = _join(eq["GENE"].astype(str) + "_" + eq["TISSUE"].astype(str))

    tf = _overlaps(tracks.get("tfbs"), chrom, pos, kinds["tfbs"])
    if tf is not None:
        tf = tf.sort_values(["SCORE", "TF"], ascending=[False, True], kind="mergesort")
        out["tfbs"] = _join(tf["TF"].astype(str) + "_" + tf["SCORE"].astype(int).astype(str))

    cc = _overlaps(tracks.get("ccre"), chrom, pos, kinds["ccre"])
    if cc is not None:
        out["ccre_id"] = _join(cc["CCRE_ID"])
        out["ccre"] = _join(cc["DESCRIPTION"])
        out["ccre_class"] = _join(cc["CCRE_CLASS"])

    dn = _overlaps(tracks.get("dnase"), chrom, pos, kinds["dnase"])
    if dn is not None:
        best = dn.sort_values("SCORE", ascending=False).iloc[0]
        out["dnase_score"] = int(best["SCORE"])
        out["dnase_sources"] = int(best["SOURCES"])

    cg = _overlaps(tracks.get("cpg"), chrom, pos, kinds["cpg"])
    if cg is not None:
        out["cpg_island"] = _join(cg["NAME"])
    return out


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


def load_genes_dict(path: str) -> dict:
    """Read a gene-info TSV into the :func:`_build_genes_dict` structure,
    cached per process (keyed on path, size and mtime)."""
    _k = ("genes", _file_key(path))
    if _k not in _LOADER_CACHE:
        _LOADER_CACHE[_k] = _build_genes_dict(pd.read_csv(path, header=0, sep="\t"))
    return _LOADER_CACHE[_k]


def candidate_gene_pool(chrom, pos: int, genes_dict: Optional[dict],
                        window: int = 500_000) -> set:
    """Protein-coding genes whose body lies within *window* bp of *pos*.

    The same window and protein-coding set the gene assignment
    (:func:`_annotate_variant`) searches, so an eQTL / GeneHancer link
    to a gene in this pool is a link to a plausible target of the locus.
    Genes without a biotype are included.
    """
    if not genes_dict:
        return set()
    entry = genes_dict.get(str(chrom))
    if entry is None:
        return set()
    lo, hi = int(pos) - int(window), int(pos) + int(window)
    pool = set()
    # Tuples are (START, END, STRAND, GENE, BIOTYPE), sorted by START.
    for start, end, _strand, gene, biotype in entry["intervals"][: bisect.bisect_right(entry["starts"], hi)]:
        if end >= lo and (biotype is None or biotype in PC_LIKE_BIOTYPES):
            pool.add(str(gene))
    return pool


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
    genehancer: Optional[dict] = None,
    eqtl: Optional[dict] = None,
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
        "promoter_proximal": None,
        "dist_to_tss": None,
        "upstream_of_gene": None,
        "snp_position": None,
        "distance_score": None,
        "biotype_weight": None,
        "promoter_bonus": None,
        "upstream_bonus": None,
        "genic_bonus": None,
        "genehancer_bonus": None,
        "gh_score": None,
        "strand": None,
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

    # 0.5.x: collect *every* in-window candidate's full priority-score
    # breakdown so the PC-preference selector below can inspect the
    # whole pool (rather than just the running argmax).  This is what
    # lets a nearby PC beat a containing lncRNA even when the lncRNA's
    # score is higher on the raw formula, and what enables the
    # ``LEFT-RIGHT`` PC-flanker join when the SNP sits inside a
    # non-PC gene.
    all_candidates: list[dict] = []

    # ------------------------------------------------------------------
    # Biologically-meaningful priority score (0.5.x rework)
    # ------------------------------------------------------------------
    # For each candidate gene inside [pos-window, pos+window]:
    #
    #   • distance_score = exp(-|distance_to_gene_body| / TAU)
    #     with TAU = 100 kb; nearby genes dominate sharply.
    #   • genic_bonus    = 2.0 if pos is inside gene body else 0.0
    #   • upstream_bonus = distance_score if SNP is *strand-aware*
    #     upstream of gene (5' side; TSS side) and outside the gene
    #     body, else 0.0.  Encodes the biological rule that regulatory
    #     variants tend to sit 5' of their target.  It decays with
    #     distance like distance_score (max 1.0 when adjacent): a flat
    #     +1.0 let any gene whose 5' end faced the SNP anywhere in the
    #     window outrank an immediately adjacent gene (e.g. PAPOLG at
    #     306 kb, 1.047, beat BCL11A at 520 bp, 0.995).
    #   • promoter_bonus = 2.0 if |pos - TSS| <= promoter_window,
    #     TSS = START for '+' strand, END for '-' strand.  Fires
    #     regardless of whether the SNP is inside the gene body.
    #
    #   priority = biotype_weight * (distance_score + genic_bonus
    #                                + upstream_bonus + promoter_bonus)
    #
    # Multiplying by biotype_weight (rather than adding it) is what
    # lets a distant protein-coding gene beat a nearby pseudogene
    # the SNP happens to sit inside — the pseudogene's 0.2 weight
    # scales its whole 3.0 genic-bonus stack down to 0.6, while a
    # PC gene 20 kb upstream gets 1.0 * (0.82 + 0.82) = 1.64.
    _TAU = 100_000.0
    # Per-SNP GeneHancer + eQTL lookups — computed once, applied
    # per candidate below.  Each ``{gene: bonus}`` dict is empty by
    # default when no connection exists.  The bonuses are already
    # scaled and can be added to the score sum directly.
    gh_bonuses = _genehancer_bonuses(chrom, pos, genehancer)
    eqtl_bonuses = _eqtl_bonuses(chrom, pos, eqtl)
    while i < len(genes):
        _entry = genes[i]
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

        # ---- Positional flanker + true-nearest tracking --------------
        if is_genic:
            if containing_gene is None:
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

        # ---- Strand-aware upstream/downstream + promoter check -------
        # Upstream of gene = 5' side.  + strand: 5' at START, so
        # upstream = pos < start.  - strand: 5' at END, so
        # upstream = pos > end.  "Genic" positions are neither upstream
        # nor downstream — the flag applies only when SNP is *outside*
        # the gene body.
        if is_genic:
            snp_position = "genic"
            upstream_of_gene = False
        elif strand == "+":
            snp_position = "upstream" if pos < start else "downstream"
            upstream_of_gene = pos < start
        else:  # '-' strand
            snp_position = "upstream" if pos > end else "downstream"
            upstream_of_gene = pos > end

        tss = start if strand == "+" else end
        dist_to_tss = abs(pos - tss)
        promoter_proximal = dist_to_tss <= promoter_window
        if promoter_proximal and upstream_of_gene:
            promoter_upstream_flag = True

        # ---- Priority-score components -------------------------------
        distance_score = math.exp(-distance / _TAU)
        biotype_weight = biotype_weights.get(biotype, 0.0) if biotype else 0.0
        genic_bonus     = 2.0 if is_genic else 0.0
        upstream_bonus  = distance_score if upstream_of_gene else 0.0
        promoter_bonus  = 2.0 if promoter_proximal else 0.0
        # Multiplying by biotype_weight is the load-bearing choice —
        # it's what lets a distant PC gene beat a nearby pseudogene.
        genehancer_bonus = gh_bonuses.get(gene, 0.0)
        eqtl_bonus = eqtl_bonuses.get(gene, 0.0)
        priority_score = biotype_weight * (
            distance_score + genic_bonus + upstream_bonus
            + promoter_bonus + genehancer_bonus + eqtl_bonus
        )
        all_candidates.append({
            "gene": gene,
            "biotype": biotype,
            "start": int(start),
            "end": int(end),
            "distance": distance,
            "promoter_flag": promoter_proximal,
            "promoter_proximal": promoter_proximal,
            "dist_to_tss": int(dist_to_tss) if promoter_proximal else None,
            "distance_score": distance_score,
            "biotype_weight": biotype_weight,
            "promoter_bonus": promoter_bonus,
            "upstream_bonus": upstream_bonus,
            "genic_bonus": genic_bonus,
            "genehancer_bonus": genehancer_bonus,
            "eqtl_bonus": eqtl_bonus,
            "snp_position": snp_position,
            "priority_score": priority_score,
            "is_genic": is_genic,
            "strand": strand,
        })

        i += 1

    # ---------------------------------------------------------------
    # PC-preference selection
    # ---------------------------------------------------------------
    # Rules, applied in order:
    #   1. No candidates at all  -> top = None
    #   2. No PC in the window   -> top = overall highest-priority (may
    #      be a lncRNA / pseudogene / etc.; nothing better is available)
    #   3. A PC contains the SNP -> top = highest-priority PC that is
    #      genic (matches the "SNP inside a PC gene" case; a nearby
    #      non-PC never wins here)
    #   4. Otherwise             -> take the highest-priority PC that
    #      is positional-left of the SNP AND the highest-priority PC
    #      that is positional-right of the SNP.  If both exist ->
    #      report as "LEFT-RIGHT".  If only one side has a PC -> that
    #      single PC is top_gene.
    top: Optional[dict] = None
    joined_pcs: Optional[tuple[dict, dict]] = None
    if all_candidates:
        _pcs = [c for c in all_candidates if c["biotype"] in PC_LIKE_BIOTYPES]
        if not _pcs:
            top = max(all_candidates, key=lambda c: c["priority_score"])
        else:
            _pcs_genic = [c for c in _pcs if c["is_genic"]]
            if _pcs_genic:
                top = max(_pcs_genic, key=lambda c: c["priority_score"])
            else:
                _pcs_left  = [c for c in _pcs if c["end"] < pos]
                _pcs_right = [c for c in _pcs if c["start"] > pos]
                _left_top  = (
                    max(_pcs_left,  key=lambda c: c["priority_score"])
                    if _pcs_left  else None
                )
                _right_top = (
                    max(_pcs_right, key=lambda c: c["priority_score"])
                    if _pcs_right else None
                )
                if _left_top and _right_top:
                    joined_pcs = (_left_top, _right_top)
                    top = _left_top  # base record; joined label built below
                elif _left_top:
                    top = _left_top
                elif _right_top:
                    top = _right_top
                else:
                    # Every PC was either at pos or partially straddled
                    # it but not flagged genic (shouldn't happen given
                    # is_genic covers start<=pos<=end, but be defensive).
                    top = max(_pcs, key=lambda c: c["priority_score"])

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
    # 0.5.x with PC-preference:
    #   * If PCs bracket the SNP but none contain it, ``top_gene`` is
    #     reported as ``"LEFT-RIGHT"`` (the highest-priority PC on each
    #     side).  Kept the historical joined-label semantic but scoped
    #     to PCs only.
    #   * Otherwise ``top_gene`` is a single gene chosen by the
    #     PC-preference cascade above (containing PC > single-side PC
    #     > non-PC fall-back).
    if top is None:
        result.update({
            "top_gene": None, "biotype": None, "priority_score": None,
            "distance": None, "promoter_flag": None,
            "promoter_proximal": None, "dist_to_tss": None,
            "upstream_of_gene": None, "snp_position": None,
            "distance_score": None, "biotype_weight": None,
            "promoter_bonus": None, "upstream_bonus": None,
            "genic_bonus": None, "genehancer_bonus": None,
            "eqtl_bonus": None, "gh_score": None, "strand": None,
        })
    elif joined_pcs is not None:
        _l, _r = joined_pcs
        result.update({
            "top_gene": f"{_l['gene']}-{_r['gene']}",
            "biotype": "intergenic_PC",
            "priority_score": (_l["priority_score"] + _r["priority_score"]) / 2,
            "distance": f"{_l['distance']}-{_r['distance']}",
            "promoter_flag": bool(_l["promoter_flag"] or _r["promoter_flag"]),
            "promoter_proximal": bool(_l["promoter_proximal"] or _r["promoter_proximal"]),
            "dist_to_tss": None,
            "upstream_of_gene": None,   # ambiguous — 5' of the right
                                        # gene, 3' of the left gene
            "snp_position": "intergenic_PC",
            "distance_score": None,
            "biotype_weight": None,
            "promoter_bonus": None,
            "upstream_bonus": None,
            "genic_bonus": None,
            "genehancer_bonus": max(
                _l.get("genehancer_bonus", 0.0) or 0.0,
                _r.get("genehancer_bonus", 0.0) or 0.0,
            ),
            "eqtl_bonus": max(
                _l.get("eqtl_bonus", 0.0) or 0.0,
                _r.get("eqtl_bonus", 0.0) or 0.0,
            ),
            "gh_score": None,
            "strand": None,
        })
    else:
        _gh = top.get("genehancer_bonus", 0.0) or 0.0
        result.update({
            "top_gene": top["gene"],
            "biotype": top["biotype"],
            "priority_score": top["priority_score"],
            "distance": top["distance"],
            "promoter_flag": top["promoter_flag"],
            "promoter_proximal": top["promoter_proximal"],
            "dist_to_tss": top["dist_to_tss"],
            "upstream_of_gene": top["snp_position"] == "upstream",
            "snp_position": top["snp_position"],
            "distance_score": top["distance_score"],
            "biotype_weight": top["biotype_weight"],
            "promoter_bonus": top["promoter_bonus"],
            "upstream_bonus": top["upstream_bonus"],
            "genic_bonus": top["genic_bonus"],
            "genehancer_bonus": _gh,
            "eqtl_bonus": top.get("eqtl_bonus", 0.0) or 0.0,
            # Report the reconstructed raw GH connection score
            # (bonus * divisor) for hits-table readability — makes
            # it obvious whether a strong or weak GH connection
            # decided the pick.
            "gh_score": _gh * _GENEHANCER_DIVISOR if _gh else None,
            "strand": top["strand"],
        })

    return result




# ---------------------------------------------------------------------------
# Internal: clumping
# ---------------------------------------------------------------------------

def _clump_by_distance(df: pd.DataFrame, window_kb: int = 250) -> pd.DataFrame:
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
            group = group.sort_values("logP", ascending=False, kind="mergesort")
        else:
            group = group.sort_values("P", ascending=True, kind="mergesort")

        kept_positions: list[int] = []
        for _, row in group.iterrows():
            if all(abs(row["POS"] - p) > window for p in kept_positions):
                clumped.append(row)
                kept_positions.append(row["POS"])

    return pd.DataFrame(clumped).sort_values(
        ["CHR", "POS"], key=natsort.natsort_keygen()
    )


def _gene_label(row) -> Optional[str]:
    """The gene label the plotters show for a hits row.

    Mirrors :func:`get_annotation_column` with ``annotate="GENE"``:
    the containing gene (``nearest_gene``) for genic hits, otherwise
    ``top_gene`` (a single gene or a ``LEFT-RIGHT`` flanker pair).
    """
    col = "nearest_gene" if bool(row.get("genic", False)) else "top_gene"
    val = row.get(col)
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return None
    val = str(val).strip()
    return val if val and val.lower() not in ("none", "nan") else None


GENE_LABEL_COL = "gene_label"
ANNOT_WINDOW_COL = "annot_window_kb"
LABEL_SHOWN_COL = "label_shown"
# Label columns whose values are gene symbols; repeated values are
# collapsed per locus when a plot is labelled by one of these.
GENE_LABEL_COLUMNS = frozenset({
    GENE_LABEL_COL, "nearest_gene", "top_gene",
    "nearest_upstream_gene", "nearest_downstream_gene",
})


def _significance_order(df: pd.DataFrame) -> pd.Index:
    """Row index of *df*, most significant first (stable)."""
    if "logP" in df.columns:
        return df["logP"].astype(float).sort_values(ascending=False, kind="mergesort").index
    if "P_UNSIGNED" in df.columns:
        return df["P_UNSIGNED"].astype(float).sort_values(ascending=False, kind="mergesort").index
    if "P" in df.columns:
        return df["P"].astype(float).sort_values(ascending=True, kind="mergesort").index
    return df.index


def _collapse_repeated_labels(df: pd.DataFrame, keys: pd.Series,
                              window_bp: int) -> pd.Series:
    """Mark one row per (chromosome, key, locus) as the label to show.

    Returns a boolean Series aligned to *df*: ``True`` for rows whose
    label should be drawn.  Rows are duplicates only when they share a
    non-null *key*, sit on the same chromosome, and lie within
    *window_bp* of each other; the most significant one is shown.
    Rows with a null key are always shown, and same-key rows farther
    apart than *window_bp* are shown as separate signals.
    """
    shown = pd.Series(True, index=df.index)
    if df.empty or "CHR" not in df.columns or "POS" not in df.columns:
        return shown

    kept_pos: dict = {}  # (chr, key) -> [positions]
    for idx in _significance_order(df):
        key = keys.get(idx)
        if key is None:
            continue
        ck = (str(df.at[idx, "CHR"]), key)
        pos = int(df.at[idx, "POS"])
        if any(abs(pos - p) <= window_bp for p in kept_pos.get(ck, ())):
            shown.at[idx] = False
            continue
        kept_pos.setdefault(ck, []).append(pos)

    n_hidden = int((~shown).sum())
    if n_hidden:
        logger.info(
            "Gene labels: %s hit(s) share a gene label with a more significant "
            "hit within %s kb and are drawn without a repeated label "
            "(points, guide lines and colours unchanged).",
            n_hidden, window_bp // 1000,
        )
    return shown


def prepare_annotation_labels(
    annotate: Optional[str],
    hits_table: Optional[pd.DataFrame],
    label_col: Optional[str] = None,
    window_kb: Optional[int] = None,
) -> tuple:
    """Resolve the label column and the rows to label for a plot.

    Returns ``(hits, label_column)``.  ``hits`` is a copy of the table
    with a boolean ``label_shown`` column.  When the labels are gene
    symbols, a hit repeating the same gene label within *window_kb* of a
    more significant hit gets ``label_shown=False``: several
    LD-independent signals in one gene read as one gene label, while
    each signal keeps its guide line and highlight colour.  When the
    labels are rsIDs or any other non-gene column, every hit is shown.
    The input table is not modified.

    *window_kb* defaults to the window the hits table was annotated
    with (``annot_window_kb`` column), else 500 kb.
    """
    if hits_table is None:
        return hits_table, "SNP"
    if hits_table.empty or not annotate:
        return hits_table, get_annotation_column(annotate, hits_table, label_col)

    table = hits_table.copy()
    table[LABEL_SHOWN_COL] = True
    wants_gene = str(annotate).upper() == "GENE" and label_col is None
    if wants_gene and GENE_LABEL_COL not in table.columns and (
        "nearest_gene" in table.columns or "top_gene" in table.columns
    ):
        # Older cached overlays predate ``gene_label``.
        table[GENE_LABEL_COL] = [
            _gene_label(r) or (str(r["SNP"]) if "SNP" in r else None)
            for _, r in table.iterrows()
        ]

    label = get_annotation_column(annotate, table, label_col)
    if label not in GENE_LABEL_COLUMNS:
        return table, label

    if window_kb is None:
        if ANNOT_WINDOW_COL in table.columns and table[ANNOT_WINDOW_COL].notna().any():
            window_kb = int(pd.to_numeric(table[ANNOT_WINDOW_COL], errors="coerce").max())
        else:
            window_kb = 500

    if label == GENE_LABEL_COL:
        # Rows whose gene_label fell back to the rsID are not genes.
        keys = pd.Series({i: _gene_label(r) for i, r in table.iterrows()}, dtype=object)
    else:
        def _k(v):
            if v is None or (isinstance(v, float) and math.isnan(v)):
                return None
            v = str(v).strip()
            return v if v and v.lower() not in ("none", "nan") else None
        keys = table[label].map(_k)
    table[LABEL_SHOWN_COL] = _collapse_repeated_labels(table, keys, int(window_kb) * 1_000)
    return table, label


def labels_to_draw(hits: pd.DataFrame) -> pd.DataFrame:
    """Rows of *hits* whose label should be drawn (all rows if unmarked)."""
    if hits is None or LABEL_SHOWN_COL not in hits.columns:
        return hits
    return hits[hits[LABEL_SHOWN_COL].astype(bool)]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def get_hits_summary_table(
    leads_df: pd.DataFrame,
    window_kb: int = 500,
    clump_window_kb: int = 250,
    table_out: Optional[str] = None,
    resources: Optional[ResourceConfig] = None,
    use_genehancer: bool = True,
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
        genes_dict = load_genes_dict(geneinfo_path)

        # Functional-annotation tracks (loaded once per hits-table
        # generation).  All are optional; a missing bundled file
        # degrades that dimension of scoring gracefully.
        gh_dict = None
        eqtl_dict = None
        if use_genehancer:
            try:
                gh_path = resources.genehancer_hg38
            except Exception:
                gh_path = None
            if gh_path is not None:
                logger.info("Loading GeneHancer connections from: %s", gh_path)
                gh_dict = load_genehancer(gh_path)
            for attr, loader, name in (
                ("eqtl_hg38", load_eqtl, "eQTL"),
            ):
                try:
                    _p = getattr(resources, attr, None)
                except Exception:
                    _p = None
                if _p is not None:
                    logger.info("Loading %s track from: %s", name, _p)
                    if name == "eQTL":
                        eqtl_dict = loader(_p)

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
            # Layer 2 annotation-representative: when the loader has
            # attached ``annot_pos`` / ``annot_snp`` columns, use
            # them for the actual gene lookup so LD-block members
            # sitting in richer regulatory context get to define the
            # annotation while the statistical lead's identity
            # remains in the CHR / POS / SNP / P columns.
            _ap = row.get("annot_pos", row["POS"])
            annotation = _annotate_variant(
                chrom=row["CHR"],
                pos=int(_ap) if _ap is not None and not pd.isna(_ap) else int(row["POS"]),
                genes_dict=genes_dict,
                window=window,
                genehancer=gh_dict,
                eqtl=eqtl_dict,
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

    if "SNP" in locus_table.columns and "CHR" in locus_table.columns:
        locus_table = locus_table.sort_values(
            ["CHR", "POS"], key=natsort.natsort_keygen()
        ).drop_duplicates(
            subset=["CHR", "SNP"], keep="first"
        )

    #if table_out is not None:
    #    outpath = table_out.replace(" ", "_").lower() + '.tsv'
    #    locus_table.to_csv(outpath, index=False, sep="\t", na_rep="None")
    #    logger.info("Locus summary written to: %s", outpath)
    
    ## Lead-SNP clumping window is *independent* of the annotation
    ## search window.  Default 250 kb matches PLINK's ``--clump-kb``
    ## convention and aligns with typical European-population LD
    ## extent (~100-200 kb) with a modest safety margin.  For
    ## African-ancestry cohorts LD is shorter (~30-50 kb) and users
    ## can pass a smaller value.  (``get_lead_snps`` in
    ## ``pycmplot.stats`` already clumps at the same window upstream
    ## of this call; this is a defensive second pass on the
    ## annotated table.)
    #
    _clumped = _clump_by_distance(locus_table, window_kb=clump_window_kb)
    #
    ## Collapse (CHR, SNP) duplicates that survive the distance
    ## clumper.  These arise when the same rsID appears at slightly
    ## different POS values across build-mixed tracks (e.g. rs123
    ## pre- vs post-liftover) — distance-based clumping keeps both
    ## because their POS separation exceeds the window, but they
    ## represent one physical variant.  Left un-deduped this
    ## clutters the cached hits overlay (``hits.<group>.tsv``) with
    ## phantom rsID duplicates that propagate into the circular
    ## plot (which consumes the overlay directly).  Keep the first
    ## occurrence per (CHR, SNP); the ``keep="first"`` here mirrors
    ## the greedy pick already done by the clumper.
    #
    if "SNP" in _clumped.columns and "CHR" in _clumped.columns:
        _clumped = _clumped.drop_duplicates(subset=["CHR", "SNP"], keep="first")

    # No gene-based row dropping here.  The old
    # ``drop_duplicates(["CHR", "nearest_gene"])`` removed rows from the
    # table itself, which also cost SNP-mode plots their independent
    # rsIDs.  Repeated gene labels are now collapsed only at plot time,
    # and only when the plot is labelled by gene (see
    # :func:`prepare_annotation_labels`).  The per-row display label and
    # the window it was assigned from are recorded for that step.
    if not _clumped.empty and (
        "nearest_gene" in _clumped.columns or "top_gene" in _clumped.columns
    ):
        _clumped = _clumped.copy()
        _clumped[GENE_LABEL_COL] = [
            _gene_label(r) or (str(r["SNP"]) if "SNP" in r else None)
            for _, r in _clumped.iterrows()
        ]
        _clumped[ANNOT_WINDOW_COL] = int(window_kb)
    if not _clumped.empty:
        # Locus half-width used for ``in_locus``; read back by the
        # plot-time highlight filter (``resolve_highlight_window_kb``).
        _clumped = _clumped.copy()
        _clumped[CLUMP_WINDOW_COL] = int(clump_window_kb)

    # Functional annotation of each lead (reporting only; does not
    # change which lead or gene was chosen).  Runs with the same switch
    # that enables the functional tracks for selection.  Annotates the
    # variant used for the gene lookup: ``annot_pos`` when the LD-block
    # tie-break chose a representative, else the lead itself.
    if use_genehancer and not _clumped.empty:
        logger.info("Adding GeneHancer / eQTL / ENCODE / CpG annotations to hits ...")
        _tracks = load_functional_tracks(resources)
        _recs = []
        for _, _row in _clumped.iterrows():
            _p = _row.get("annot_pos", _row["POS"])
            _p = int(_p) if _p is not None and not pd.isna(_p) else int(_row["POS"])
            _recs.append(annotate_functional(str(_row["CHR"]), _p, _tracks))
        _func = pd.DataFrame(_recs, index=_clumped.index, columns=list(FUNCTIONAL_COLUMNS))
        _clumped = pd.concat([_clumped.drop(columns=[c for c in FUNCTIONAL_COLUMNS if c in _clumped.columns]), _func], axis=1)

    if table_out is not None:
        outpath = table_out.replace(" ", "_").lower() + '.tsv'
        _clumped.to_csv(outpath, index=False, sep="\t", na_rep="None")
        logger.info("Locus summary written to: %s", outpath)
    
    #return locus_table
    return _clumped


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
                # One column for every row.  (Previously the column was
                # re-chosen per row inside a loop, so the LAST hit's
                # genic status decided the label column for all hits.)
                # ``gene_label`` already holds the per-row choice:
                # containing gene when genic, else ``top_gene``.
                for _c in (GENE_LABEL_COL, "top_gene", "nearest_gene"):
                    if _c in hits_table.columns:
                        label_clm = _c
                        break
                else:
                    logger.warning(
                        "No gene annotation columns in hits table; "
                        "falling back to 'SNP'."
                    )

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


CLUMP_WINDOW_COL = "clump_window_kb"
DEFAULT_CLUMP_WINDOW_KB = 250  # loader default (``clump_window_kb``)


def resolve_highlight_window_kb(highlight_window_kb=None, hits_table=None) -> int:
    """Locus half-width (kb) for the plot-time highlight filter.

    Order: explicit *highlight_window_kb*; else the ``clump_window_kb``
    the loader recorded in *hits_table*; else the loader default
    (250 kb).  This is the window the loader used to mark ``in_locus``.
    """
    if highlight_window_kb is not None:
        return int(highlight_window_kb)
    if (hits_table is not None and isinstance(hits_table, pd.DataFrame)
            and CLUMP_WINDOW_COL in hits_table.columns):
        _v = pd.to_numeric(hits_table[CLUMP_WINDOW_COL], errors="coerce").dropna()
        if not _v.empty:
            return int(_v.max())
    return DEFAULT_CLUMP_WINDOW_KB


def filter_in_locus_by_threshold(sig, threshold, window: int,
                                 chr_col: str = "CHR", pos_col: str = "POS"):
    """Plot-time, per-locus highlight filter for ONE track.

    *sig* is that track's highlighted (``in_locus``) variants.  A row is
    kept when it lies within *window* bp, on the same chromosome, of a
    row of *sig* that passes *threshold*: ``P_UNSIGNED >= threshold``
    for signed statistics (the loader adds ``P_UNSIGNED``), otherwise
    ``P <= threshold``.

    Each locus is judged by its own lead.  A locus is every in-locus
    variant within *window* (the loader's ``clump_window_kb``) of the
    track's lead, and the lead is its most significant variant, so a
    passing locus keeps all of its variants and a failing locus loses
    all of them.  If nothing passes, nothing is kept.  ``threshold=None``
    returns *sig* unchanged.
    """
    import numpy as np

    if sig is None or not isinstance(sig, pd.DataFrame) or sig.empty or threshold is None:
        return sig
    thr = float(threshold)
    if "P_UNSIGNED" in sig.columns:
        passing = pd.to_numeric(sig["P_UNSIGNED"], errors="coerce").to_numpy() >= thr
    else:
        passing = pd.to_numeric(sig["P"], errors="coerce").to_numpy() <= thr
    if not passing.any():
        return sig.iloc[0:0]

    chr_arr = sig[chr_col].astype(str).to_numpy()
    pos_arr = sig[pos_col].to_numpy(dtype=np.int64)
    keep = np.zeros(len(sig.index), dtype=bool)
    w = int(window)
    for c in np.unique(chr_arr[passing]):
        on_c = chr_arr == c
        anchors = np.sort(pos_arr[on_c & passing])
        p = pos_arr[on_c]
        # distance from each row to its nearest passing anchor
        i = np.searchsorted(anchors, p)
        left = np.abs(p - anchors[np.clip(i - 1, 0, len(anchors) - 1)])
        right = np.abs(anchors[np.clip(i, 0, len(anchors) - 1)] - p)
        keep[on_c] = np.minimum(left, right) <= w
    return sig[keep]


def filter_hits_by_signif(hits_table, signif_threshold, p_col: str = "P"):
    """Return the subset of *hits_table* meeting *signif_threshold*.

    Auto-detects signed statistics: if any value in ``hits_table[p_col]``
    is negative, filters by ``|value| >= threshold`` (both tails).
    Otherwise treats the column as a p-value and filters by
    ``value <= threshold``.

    ``signif_threshold=None`` (or an empty / column-less hits table)
    returns the input unchanged.
    """
    import numpy as np
    import pandas as pd

    if (
        hits_table is None
        or not isinstance(hits_table, pd.DataFrame)
        or hits_table.empty
        or signif_threshold is None
        or p_col not in hits_table.columns
    ):
        return hits_table
    vals = hits_table[p_col].to_numpy(dtype=float)
    thresh = float(signif_threshold)
    if np.any(vals < 0):
        mask = np.abs(vals) >= thresh
    else:
        mask = vals <= thresh
    return hits_table[mask].reset_index(drop=True)


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
