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
    :func:`~pycmplot.annotation._annotate_and_prioritize_variant` when
    ranking candidate genes at intergenic loci.  Protein-coding genes and
    general ``'gene'`` entries receive the highest weight (1.00); various
    pseudogene classes receive the lowest (0.20–0.30).

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
    "gene":                                   1.00,
    "protein_coding":                         1.00,
    "miRNA":                                  0.75,
    "lncRNA":                                 0.70,
    "ncRNA":                                  0.70,
    "lincRNA":                                0.70,
    "ribozyme":                               0.70,
    "snRNA":                                  0.65,
    "snoRNA":                                 0.65,
    "scaRNA":                                 0.65,
    "vault_RNA":                              0.60,
    "antisense":                              0.30,
    "rRNA":                                   0.55,
    "processed_transcript":                   0.50,
    "transcribed_processed_pseudogene":       0.45,
    "transcribed_unitary_pseudogene":         0.40,
    "transcribed_unprocessed_pseudogene":     0.35,
    "processed_pseudogene":                   0.30,
    "pseudogene":                             0.20,
    "unprocessed_pseudogene":                 0.20,
}

# ---------------------------------------------------------------------------
# Standard chromosome order (autosomes + sex + MT)
# ---------------------------------------------------------------------------
CHROM_ORDER: list[str] = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
