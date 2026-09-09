"""
pycmplot.constants
==================

Genome-level constants shared across pycmplot modules.

Contents
--------
hg38_chr_lengths : dict
    GRCh38 chromosome lengths in base-pairs for chromosomes 1–22, X, and Y.
    Used to compute Circos sector sizes when summary statistics do not
    cover the full chromosome.

BIOTYPE_WEIGHTS : dict
    Numeric priority weights for Ensembl gene biotypes used by
    :func:`~pycmplot.annotation._annotate_variant` when ranking
    candidate genes at intergenic loci.  Values are grouped into four
    tiers reflecting the current Ensembl biotype hierarchy (release
    116, June 2026):

    * **Protein-coding & immune-receptor genes** (weight 0.85 – 1.00) —
      ``protein_coding``, the eight ``IG_*_gene`` / ``TR_*_gene``
      classes, plus ``protein_coding_LoF`` and
      ``protein_coding_CDS_not_defined`` which are still
      protein-producing in some individuals/isoforms.
    * **Non-coding RNAs** (0.55 – 0.75) — small ncRNAs (miRNA,
      piRNA, siRNA, snRNA, snoRNA, scaRNA, tRNA, ribozyme, vaultRNA,
      rRNA, miscRNA, plus their ``Mt_`` mitochondrial counterparts)
      and long ncRNAs (lncRNA / lincRNA and every ``antisense`` /
      ``sense_*`` / ``retained_intron`` / ``macro_lncRNA`` subtype
      Ensembl classifies as a long-ncRNA sibling).
    * **Pseudogenes** (0.20 – 0.45) — ordered by evidence of
      expression: ``transcribed_*`` and ``translated_*`` subtypes
      score above plain processed / unprocessed / IG / TR pseudogenes.
    * **Speculative / provisional** (0.15 – 0.35) — ``TEC``,
      ``readthrough``, ``artifact``.

    ``antisense`` moved from 0.30 to 0.65 in 0.4.x — Ensembl
    classifies it as a long-ncRNA subtype, not a pseudogene-tier
    biotype.  Aliases for legacy spellings (``vault_RNA`` ↔
    ``vaultRNA``, ``misc_RNA`` ↔ ``miscRNA``,
    ``3prime_overlapping_ncRNA`` ↔ ``3_prime_overlapping_ncRNA``,
    ``lincRNA`` ↔ ``long_intergenic_ncRNA``) are present so a
    biotype coming in with either spelling resolves to the same
    weight.  Reference:
    https://jun2026.archive.ensembl.org/info/genome/genebuild/biotypes.html

CHROM_ORDER : list of str
    Standard chromosome ordering for autosomes 1–22 followed by X, Y, and
    MT.  Used for natural-sort validation and display ordering.

Notes
-----
``hg38_chr_lengths`` reflects the GRCh38 primary assembly (GCA_000001405).
Values may differ slightly from builds that include alternate contigs or
patches.
"""

# ---------------------------------------------------------------------------
# hg38 chromosome lengths (GRCh38)
# ---------------------------------------------------------------------------
hg38_chr_lengths: dict[str, int] = {
    "chr1":  249698942,
    "chr2":  242508799,
    "chr3":  198450956,
    "chr4":  190424264,
    "chr5":  181630948,
    "chr6":  170805979,
    "chr7":  159345973,
    "chr8":  145138636,
    "chr9":  138688728,
    "chr10": 133797422,
    "chr11": 135186938,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 108136338,
    "chr15": 102439437,
    "chr16":  92211104,
    "chr17":  83836422,
    "chr18":  80373285,
    "chr19":  58617616,
    "chr20":  64444167,
    "chr21":  46709983,
    "chr22":  51857516,
    "chrX":  156040895,
    "chrY":   57264655,
}

# ---------------------------------------------------------------------------
# Gene biotype weights used for nearest-gene prioritisation
# ---------------------------------------------------------------------------
BIOTYPE_WEIGHTS: dict[str, float] = {
    # ---- Protein-coding + immune-receptor genes ----------------------
    # ``gene`` catches unqualified entries in older references; treated
    # as protein_coding by default.  ``IG_*`` and ``TR_*`` genes
    # undergo somatic recombination in immune tissues but are
    # functionally protein-producing.  ``protein_coding_LoF`` /
    # ``protein_coding_CDS_not_defined`` are still coding for some
    # individuals or isoforms — down-weighted slightly relative to
    # canonical protein_coding but well above ncRNA/pseudogene tiers.
    "gene":                                   1.00,
    "protein_coding":                         1.00,
    "protein_coding_LoF":                     0.90,
    "protein_coding_CDS_not_defined":         0.85,
    "IG_C_gene":                              1.00,
    "IG_D_gene":                              1.00,
    "IG_J_gene":                              1.00,
    "IG_V_gene":                              1.00,
    "TR_C_gene":                              1.00,
    "TR_D_gene":                              1.00,
    "TR_J_gene":                              1.00,
    "TR_V_gene":                              1.00,
    # NMD transcripts still contain an ORF but are flagged for
    # targeted degradation — score below coding, above ncRNA.
    "nonsense_mediated_decay":                0.55,
    "non_stop_decay":                         0.55,

    # ---- Small non-coding RNAs (regulatory) --------------------------
    "miRNA":                                  0.75,
    "piRNA":                                  0.70,
    "siRNA":                                  0.70,
    "ribozyme":                               0.70,
    "snRNA":                                  0.65,
    "snoRNA":                                 0.65,
    "scaRNA":                                 0.65,
    "tRNA":                                   0.60,
    "Mt_tRNA":                                0.60,
    "vaultRNA":                               0.60,
    "vault_RNA":                              0.60,   # legacy spelling
    "rRNA":                                   0.55,
    "Mt_rRNA":                                0.55,
    "miscRNA":                                0.55,
    "misc_RNA":                               0.55,   # legacy spelling

    # ---- Long non-coding RNAs (lncRNA subfamily) --------------------
    # Ensembl 116 classifies antisense, sense_intronic, sense_overlapping,
    # 3'-overlapping ncRNA, macro_lncRNA and retained_intron as lncRNA
    # subtypes.  Weighted together with lincRNA (long intergenic ncRNA).
    "lncRNA":                                 0.70,
    "lincRNA":                                0.70,
    "long_intergenic_ncRNA":                  0.70,   # legacy spelling
    "ncRNA":                                  0.70,
    "antisense":                              0.65,   # 0.30 -> 0.65 (see docstring)
    "antisense_RNA":                          0.65,
    "3prime_overlapping_ncRNA":               0.65,
    "3_prime_overlapping_ncRNA":              0.65,   # legacy spelling
    "macro_lncRNA":                           0.65,
    "non_coding":                             0.65,
    "sense_intronic":                         0.60,
    "sense_overlapping":                      0.60,
    "retained_intron":                        0.55,
    # ``processed_transcript`` is Ensembl's parent category for
    # anything without an ORF that doesn't fit a more specific slot;
    # kept below the specific lncRNA/ncRNA subtypes.
    "processed_transcript":                   0.50,

    # ---- Pseudogenes -------------------------------------------------
    # Ordered by evidence of expression: transcribed / translated
    # subtypes score above plain processed / unprocessed.
    "transcribed_processed_pseudogene":       0.45,
    "transcribed_unitary_pseudogene":         0.40,
    "translated_processed_pseudogene":        0.40,
    "transcribed_unprocessed_pseudogene":     0.35,
    "translated_unprocessed_pseudogene":      0.35,
    "polymorphic_pseudogene":                 0.30,
    "processed_pseudogene":                   0.30,
    "unitary_pseudogene":                     0.25,
    "unprocessed_pseudogene":                 0.20,
    "pseudogene":                             0.20,
    "IG_pseudogene":                          0.20,
    "IG_C_pseudogene":                        0.20,
    "IG_J_pseudogene":                        0.20,
    "IG_V_pseudogene":                        0.20,
    "TR_pseudogene":                          0.20,
    "TR_J_pseudogene":                        0.20,
    "TR_V_pseudogene":                        0.20,

    # ---- Speculative / provisional biotypes -------------------------
    # These are low-confidence classifications; kept above pseudogene
    # (they might still be real genes) but below any confirmed ncRNA.
    "stop_codon_readthrough":                 0.55,
    "readthrough":                            0.35,
    "TEC":                                    0.30,
    "artifact":                               0.15,
}

# ---------------------------------------------------------------------------
# Standard chromosome order (autosomes + sex + MT)
# ---------------------------------------------------------------------------
CHROM_ORDER: list[str] = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
