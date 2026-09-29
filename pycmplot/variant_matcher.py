"""
pycmplot.variant_matcher
=========================

Cross-source variant identity resolver.

GWAS summary-stats files and LD reference panels frequently name the
same physical variant with different identifier conventions:

* ``rs12345`` — dbSNP rsID
* ``1:636975:A:G`` — ``CHR:POS:REF:ALT`` (Ensembl / PLINK-modern)
* ``chr1:636975:A:G`` — same with ``chr`` prefix (UCSC)
* ``1:636975:G:A`` — same variant, alleles swapped
* ``1:636975`` — position-only (lossy)

``VariantMatcher`` normalises every registered identifier to a canonical
key so that lookups across heterogeneous sources always resolve to the
same variant.  It's the read-time harmonisation layer that closes the
"missing from LD reference" fallback in :func:`pycmplot.stats.clump`
whenever the mismatch is a naming issue rather than a biological one.

Canonical key
-------------

* Full form:  ``(CHR, POS, frozenset({REF, ALT}))`` — allele-orderless
  so ``"1:636975:A:G"`` and ``"1:636975:G:A"`` map to the same key.
* Position-only fallback: ``(CHR, POS, None)`` when alleles are unknown.

Typical use
-----------

Map the significant hits onto the LD graph's *existing* node names
(the graph is never relabelled):

    from pycmplot.variant_matcher import ld_id_map

    ld_names = ld_id_map(sig_df, graph.snp_names, graph.snp_index,
                         snp_col="SNP", chr_col="CHR", pos_col="POS",
                         ld_chrom=graph.snp_chrom, ld_pos=graph.snp_pos)
    # {sumstats_id: ld_node_name} -- look r² up with the mapped name

:func:`pycmplot.stats.clump` does this internally when called with
``harmonize_variants=True``.  The older graph-relabelling route
(``add_many`` + :meth:`VariantMatcher.rename_map` +
``LDGraph.remap_ids``) still works but touches every LD node and
mutates the graph in place.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from typing import Iterable, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# CHR:POS:REF:ALT and variants (chr prefix, upper/lower case) — the
# canonical parse for structured SNP ID strings.  Captures CHR, POS,
# REF, ALT groups.  Separators between the four fields can be any
# of ``:``, ``_``, ``/`` or ``|`` (and mixes thereof — e.g.
# ``chr1:63748574_G_A`` is a real bcftools-lifted format), so the
# same variant identifier survives round-trips through different
# genotyping / imputation pipelines.  ``-`` and ``.`` are
# deliberately excluded from the separator set because they double
# as deletion sentinels (``-``) and version dots in allele strings.
_SEP = r"[:_/|]"
_CPRA_RE = re.compile(
    rf"^(?:chr)?([\dXYMTxymt]+){_SEP}(\d+){_SEP}"
    rf"([ACGTNacgtn\-]+){_SEP}([ACGTNacgtn\-]+)$"
)
# CHR:POS — position-only fallback (lossy, no alleles).  Same
# separator flexibility.
_CP_RE = re.compile(rf"^(?:chr)?([\dXYMTxymt]+){_SEP}(\d+)$")

# CHR:POS:TYPE — GWAS-Catalog / meta-analysis style where the third
# field is a variant-class label instead of alleles.  Recognised
# tokens (case-insensitive): SNV, SNP, MNV, MNP, INDEL, DEL, INS,
# CNV, DELINS, SUB.  Treated as position-only for matching — the
# class label is not part of the canonical key so ``1:7535638:SNV``
# and ``1:7535638:SNP`` and ``1:7535638`` all resolve together.
_VAR_TYPE_TOKENS = r"SNV|SNP|MNV|MNP|INDEL|DEL|INS|CNV|DELINS|SUB"
_CPT_RE = re.compile(
    rf"^(?:chr)?([\dXYMTxymt]+){_SEP}(\d+){_SEP}(?:{_VAR_TYPE_TOKENS})$",
    re.IGNORECASE,
)


@dataclass
class VariantRec:
    """One canonical variant with every registered alias."""
    chrom: str
    pos: int
    ref: Optional[str]  = None      # None when only CHR:POS known
    alt: Optional[str]  = None
    aliases: set[str] = field(default_factory=set)

    @property
    def canonical_id(self) -> str:
        """Human-readable canonical string: CHR:POS:{REF/ALT} or CHR:POS."""
        if self.ref is not None and self.alt is not None:
            # Sort alleles so the canonical string is orderless.
            a, b = sorted((self.ref.upper(), self.alt.upper()))
            return f"{self.chrom}:{self.pos}:{a}:{b}"
        return f"{self.chrom}:{self.pos}"


class VariantMatcher:
    """Cross-source variant identity resolver — see module docstring."""

    def __init__(self) -> None:
        # canonical tuple -> VariantRec
        self._canon: dict[tuple, VariantRec] = {}
        # alias string -> canonical tuple
        self._alias_to_canon: dict[str, tuple] = {}
        # (chrom, pos) -> canonical tuple (first-hit wins), for
        # position-only queries against full-alleles records.
        self._pos_index: dict[tuple[str, int], tuple] = {}

    # ------------------------------------------------------------------
    # Building
    # ------------------------------------------------------------------

    # Chromosome aliases — collapses every synonym for a sex or
    # organelle chromosome to a single canonical name so
    # ``chrX == chr23 == X == 23`` all resolve to the same variant.
    # Follows the PLINK numeric convention: 23=X, 24=Y, 25=XY
    # (pseudo-autosomal), 26=MT.
    _CHR_ALIASES: dict[str, str] = {
        "23": "X",  "X": "X",  "CHRX": "X",  "CHR23": "X",
        "24": "Y",  "Y": "Y",  "CHRY": "Y",  "CHR24": "Y",
        "25": "XY", "XY": "XY", "CHRXY": "XY", "CHR25": "XY",
        "26": "M",  "M": "M",  "MT": "M",  "CHRM": "M", "CHRMT": "M", "CHR26": "M",
    }

    @classmethod
    def _normalise_chr(cls, chrom: str) -> str:
        s = str(chrom).strip()
        s_up = s.upper()
        # Consult the alias map first (covers "chr23" -> "X" without
        # stripping the prefix twice).
        if s_up in cls._CHR_ALIASES:
            return cls._CHR_ALIASES[s_up]
        if s_up.startswith("CHR"):
            s = s[3:]
            s_up = s.upper()
            if s_up in cls._CHR_ALIASES:
                return cls._CHR_ALIASES[s_up]
        return s_up

    @classmethod
    def _parse_identifier(cls, ident: str) -> Optional[dict]:
        """Try to parse a SNP-ID string into structured fields.

        Returns ``None`` for rsIDs and free-form identifiers that
        don't match the CHR:POS[:REF:ALT] template — those are still
        registrable via :meth:`add` with explicit ``chr``/``pos`` args.
        """
        if ident is None:
            return None
        s = str(ident).strip()
        if not s:
            return None
        m = _CPRA_RE.match(s)
        if m:
            return {
                "chrom": cls._normalise_chr(m.group(1)),
                "pos":   int(m.group(2)),
                "ref":   m.group(3).upper(),
                "alt":   m.group(4).upper(),
            }
        m = _CP_RE.match(s)
        if m:
            return {
                "chrom": cls._normalise_chr(m.group(1)),
                "pos":   int(m.group(2)),
                "ref":   None,
                "alt":   None,
            }
        # CHR:POS:TYPE — GWAS-Catalog style with a variant-class token
        # (SNV / SNP / MNV / MNP / INDEL / DEL / INS / CNV / DELINS /
        # SUB) in the third field.  Falls back to position-only
        # matching so a hit at ``1:7535638:SNV`` resolves to the
        # same variant as one at ``1:7535638:A:G`` after a subsequent
        # allele-carrying registration upgrades the record.
        m = _CPT_RE.match(s)
        if m:
            return {
                "chrom": cls._normalise_chr(m.group(1)),
                "pos":   int(m.group(2)),
                "ref":   None,
                "alt":   None,
            }
        return None

    @staticmethod
    def _canon_key(chrom: str, pos: int,
                   ref: Optional[str], alt: Optional[str]) -> tuple:
        """Canonical hashable key: alleles are sorted so ``A/G`` and
        ``G/A`` collide, and alleles are set to ``None`` when unknown
        so a position-only lookup can still find a full-form record.
        """
        if ref is not None and alt is not None:
            a, b = sorted((ref.upper(), alt.upper()))
            return (chrom, int(pos), a, b)
        return (chrom, int(pos), None, None)

    def add(self, ident: Optional[str] = None, *,
            chrom: Optional[str] = None, pos: Optional[int] = None,
            ref: Optional[str] = None, alt: Optional[str] = None,
            aliases: Iterable[str] = ()) -> Optional[str]:
        """Register a variant.

        Fields are inferred from ``ident`` when it matches
        ``CHR:POS[:REF:ALT]``; otherwise pass ``chrom``/``pos`` (+
        optional ``ref``/``alt``) explicitly.  ``ident`` and every
        ``aliases`` entry become searchable via :meth:`match`.

        Returns the canonical ID string, or ``None`` if not enough
        information was supplied to place the variant on the genome.
        """
        parsed = None
        if ident is not None:
            parsed = self._parse_identifier(ident)
        if parsed is not None:
            chrom = chrom or parsed["chrom"]
            pos   = pos   if pos   is not None else parsed["pos"]
            ref   = ref   or parsed["ref"]
            alt   = alt   or parsed["alt"]
        if chrom is None or pos is None:
            return None
        chrom_n = self._normalise_chr(chrom)
        # Prefer a full-alleles canonical key; if a position-only key
        # already exists, fold the position-only aliases into the
        # full-alleles record when we later learn REF/ALT.
        full_key = self._canon_key(chrom_n, pos, ref, alt)
        pos_key  = self._canon_key(chrom_n, pos, None, None)

        rec = self._canon.get(full_key)
        if rec is None:
            rec = self._canon.get(pos_key)
            if rec is not None and ref is not None and alt is not None:
                # Upgrade the pos-only record to a full-alleles one.
                self._canon.pop(pos_key)
                rec.ref, rec.alt = ref.upper(), alt.upper()
                self._canon[full_key] = rec
                for alias in rec.aliases:
                    self._alias_to_canon[alias] = full_key
            else:
                rec = VariantRec(chrom=chrom_n, pos=int(pos),
                                 ref=(ref.upper() if ref else None),
                                 alt=(alt.upper() if alt else None))
                self._canon[full_key if ref and alt else pos_key] = rec
        canonical = full_key if ref and alt else pos_key

        if ident is not None:
            s = str(ident)
            rec.aliases.add(s)
            self._alias_to_canon[s] = canonical
            # Also register the canonical-form string (alleles sorted)
            # and the "chr"-prefixed version so lookup handles both.
            canon_str = rec.canonical_id
            rec.aliases.add(canon_str)
            self._alias_to_canon[canon_str] = canonical
            if not s.startswith("chr"):
                self._alias_to_canon.setdefault("chr" + s, canonical)
        for a in aliases:
            a = str(a)
            rec.aliases.add(a)
            self._alias_to_canon[a] = canonical
        # Maintain a position-only index so ``match("chr:pos")`` can
        # find full-alleles records too.  First-hit wins per (chrom,
        # pos) — multi-allelic sites keep the first-added variant.
        self._pos_index.setdefault((chrom_n, int(pos)), canonical)
        return rec.canonical_id

    def add_many(self, idents: Iterable[str]) -> int:
        """Register a batch of ``CHR:POS[:REF:ALT]`` identifiers.

        Skips (silently) any that don't parse; returns the count of
        successful additions.
        """
        n = 0
        for i in idents:
            if self.add(i) is not None:
                n += 1
        return n

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame,
                       snp_col: str = "SNP",
                       chr_col: str = "CHR",
                       pos_col: str = "POS",
                       ref_col: Optional[str] = None,
                       alt_col: Optional[str] = None) -> "VariantMatcher":
        """Build a matcher from a summary-stats DataFrame.

        Uses the *snp_col* value as the primary alias and enriches
        with structured fields from *chr_col*/*pos_col*/*ref_col*/
        *alt_col* when present.
        """
        vm = cls()
        cols = df.columns
        has_ref = ref_col is not None and ref_col in cols
        has_alt = alt_col is not None and alt_col in cols
        for row in df.itertuples(index=False):
            row_d = row._asdict() if hasattr(row, "_asdict") else dict(zip(df.columns, row))
            ident = row_d.get(snp_col)
            vm.add(
                ident,
                chrom=row_d.get(chr_col),
                pos=row_d.get(pos_col),
                ref=row_d.get(ref_col) if has_ref else None,
                alt=row_d.get(alt_col) if has_alt else None,
            )
        return vm

    # ------------------------------------------------------------------
    # Querying
    # ------------------------------------------------------------------

    def match(self, ident: str) -> Optional[str]:
        """Return the canonical ID string for *ident*, or ``None``.

        Tries in order:
          1. direct alias lookup (fast — the string was registered)
          2. structural parse + full-alleles key
          3. structural parse + position-only key (allele-blind)
        """
        s = str(ident)
        canon = self._alias_to_canon.get(s)
        if canon is None:
            parsed = self._parse_identifier(s)
            if parsed is not None:
                full = self._canon_key(parsed["chrom"], parsed["pos"],
                                       parsed["ref"], parsed["alt"])
                if full in self._canon:
                    canon = full
                else:
                    pos_only = self._canon_key(parsed["chrom"],
                                               parsed["pos"], None, None)
                    if pos_only in self._canon:
                        canon = pos_only
                    else:
                        # Position-only query against a full-alleles
                        # record: consult the (chrom, pos) index.
                        canon = self._pos_index.get(
                            (parsed["chrom"], parsed["pos"])
                        )
        if canon is None:
            return None
        return self._canon[canon].canonical_id

    def rename_map(self, source: Iterable[str],
                   target: Iterable[str]) -> dict[str, str]:
        """Return ``{source_id: target_id}`` for every pair where both
        source and target IDs resolve to the same canonical variant.

        Useful for :meth:`pycmplot.ld.LDGraph.remap_ids`.  For clumping
        prefer :func:`ld_id_map`, which scales with the number of hits
        rather than the size of the LD reference and leaves the graph
        untouched.
        """
        target_canon: dict[str, str] = {}
        for t in target:
            c = self.match(t)
            if c is not None:
                target_canon.setdefault(c, str(t))
        out: dict[str, str] = {}
        for s in source:
            c = self.match(s)
            if c is None:
                continue
            t = target_canon.get(c)
            if t is not None and str(s) != t:
                out[str(s)] = t
        return out

    # ------------------------------------------------------------------
    # Introspection
    # ------------------------------------------------------------------

    def __len__(self) -> int:
        return len(self._canon)

    def __repr__(self) -> str:
        return (f"<VariantMatcher n_variants={len(self._canon):,} "
                f"n_aliases={len(self._alias_to_canon):,}>")

    def __contains__(self, ident: str) -> bool:
        return self.match(ident) is not None


# ----------------------------------------------------------------------
# Hits -> LD-graph name mapping (graph left untouched)
# ----------------------------------------------------------------------

# Leading ``CHR<sep>POS`` of any structured ID; used for a vectorised
# position pre-filter over the LD node names.
_CP_HEAD = rf"^(?:chr)?([\dXYMTxymt]+){_SEP}(\d+)"


def ld_id_map(sig: pd.DataFrame,
              ld_names: Iterable[str],
              ld_index: dict,
              snp_col: str = "SNP",
              chr_col: str = "CHR",
              pos_col: str = "POS",
              *,
              ld_chrom: Optional[np.ndarray] = None,
              ld_pos: Optional[np.ndarray] = None) -> dict[str, str]:
    """Map each significant hit's ID onto an LD-graph node name.

    Returns ``{sumstats_id: ld_node_name}``.  IDs already present in
    *ld_index* map to themselves; the rest are resolved through a
    :class:`VariantMatcher` built from the unresolved hits only.

    *ld_chrom* / *ld_pos* are the graph's stored position key
    (:attr:`pycmplot.ld.LDGraph.snp_chrom` / ``snp_pos``).  When given,
    candidate LD nodes are selected with a numpy ``isin`` on positions
    and no LD name is regex-parsed except the few candidates; this also
    lets rsID-named nodes match by position.  Without them, the LD
    names are scanned once with a vectorised regex (slower, and rsID
    nodes can only match by exact string).

    Matching is allele-strict on the LD side: an LD node whose name
    carries alleles matches a hit with the same (orderless) alleles or
    a hit with unknown alleles, never a hit with a different ALT at the
    same position.  LD nodes without alleles (position-only names,
    rsIDs) match by position; at a multi-allelic hit site the first
    registered hit wins.

    Hits that cannot be resolved are simply absent from the result.
    """
    ids = sig[snp_col].astype(str)
    out: dict[str, str] = {s: s for s in ids.unique() if s in ld_index}
    todo = sig.loc[~ids.isin(out.keys())]
    todo = todo.dropna(subset=[chr_col, pos_col])
    if todo.empty:
        return out

    vm = VariantMatcher.from_dataframe(todo, snp_col=snp_col,
                                       chr_col=chr_col, pos_col=pos_col)
    hit_chr = [VariantMatcher._normalise_chr(c) for c in todo[chr_col]]
    hit_pos = todo[pos_col].astype(np.int64).to_numpy()
    want = set(zip(hit_chr, hit_pos.tolist()))

    if not isinstance(ld_names, (np.ndarray, pd.Series, list, tuple)):
        ld_names = list(ld_names)
    names = np.asarray(ld_names).astype(str)

    if ld_chrom is not None and ld_pos is not None:
        ld_chrom = np.asarray(ld_chrom)
        ld_pos = np.asarray(ld_pos)
        cand_idx = np.flatnonzero(np.isin(ld_pos, np.unique(hit_pos)))
        cands = [(names[i], str(ld_chrom[i]), int(ld_pos[i])) for i in cand_idx
                 if (str(ld_chrom[i]), int(ld_pos[i])) in want]
    else:
        ext = pd.Series(names, dtype=str).str.extract(_CP_HEAD)
        chr_map = {u: VariantMatcher._normalise_chr(u)
                   for u in ext[0].dropna().unique()}
        keys = ext[0].map(chr_map) + ":" + ext[1].str.lstrip("0")
        want_s = {f"{c}:{p}" for c, p in want}
        cands = []
        for n in names[keys.isin(want_s).to_numpy()]:
            p = VariantMatcher._parse_identifier(n)
            cands.append((n, p["chrom"], p["pos"]))

    canon_to_ld: dict[tuple, str] = {}
    for n, chrom, pos in cands:
        p = VariantMatcher._parse_identifier(n)
        # Alleles only count when the name's own position agrees with
        # the stored key; otherwise treat the node as position-only.
        if p is not None and p["chrom"] == chrom and p["pos"] == pos:
            ref, alt = p["ref"], p["alt"]
        else:
            ref = alt = None
        k = VariantMatcher._canon_key(chrom, pos, ref, alt)
        if k not in vm._canon:
            k = VariantMatcher._canon_key(chrom, pos, None, None)
            if k not in vm._canon:
                if ref is not None:
                    continue            # different ALT -> different variant
                k = vm._pos_index.get((chrom, pos))
                if k is None:
                    continue
        canon_to_ld.setdefault(k, n)

    for s in todo[snp_col].astype(str):
        k = vm._alias_to_canon.get(s)
        if k is not None and k in canon_to_ld:
            out[s] = canon_to_ld[k]
    return out
