Changelog
=========

All notable changes to **pycmplot** are documented here.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_
and this project adheres to `Semantic Versioning <https://semver.org/>`_.

---


0.4.3 - 2026-09-30
------------------------------------------------------------------------------

**Added**

- **``pycmplot.LDGraph`` — CSR-encoded sparse LD reference.**
  Optional companion to :func:`pycmplot.clump` for LD-based
  clumping.  Stores pairwise r² as three contiguous numpy arrays
  (``indptr`` int32, ``indices`` int32, ``data`` float32) plus a
  fixed-width Unicode SNP-name table — ~24 bytes per edge vs
  ~1 kB per dict-of-dict entry, so a 1 M SNP × 10 neighbour
  reference lives in ~240 MB rather than ~10 GB of Python heap.

  Key features:

    * ``LDGraph.from_plink_ld(path, r2_threshold=0.1)`` — build
      from a PLINK ``.ld`` file, filtering below-threshold edges
      at build time (so ``absence in graph == independent`` by
      construction; query paths never need a threshold comparison).
    * ``LDGraph.save(path)`` / ``LDGraph.load(path)`` — round-trip
      to a compressed ``.npz`` with a JSON metadata sidecar; build
      once, reload in milliseconds across runs.
    * ``LDGraph.r2(a, b)`` / ``.neighbors(snp, r2_min=None)`` /
      ``.any_linked(snp, others, r2_min=None)`` — O(log k) sorted-
      slice lookups.  ``any_linked`` is the fast inner loop for
      greedy clumping; ``neighbors`` is the primitive future
      regional-plot LD-tier colouring builds on.

  ``pycmplot.clump()`` now accepts an ``LDGraph`` instance
  directly (preferred for reuse), a ``.npz`` path (auto-loads),
  a PLINK ``.ld`` path (auto-builds), a dict-of-dict (legacy),
  or ``None`` (distance-only fallback).  Regressions:
  ``test_ld_graph_from_plink_ld_roundtrip`` and
  ``test_ld_clump_accepts_ldgraph``.

- **LDGraph on-disk footprint reductions.**  ``save()`` gains three
  compression optimisations stackable via the ``compact=True``
  convenience flag or individually via ``symmetric_only``,
  ``quantize_r2``, ``compression`` / ``compression_level``:

    * **Symmetric-only storage** — only canonical edges
      (src < dst) written; reverse edges reconstructed at load
      time so in-memory queries stay O(log k).  ~40% reduction
      of the ``indices`` + ``data`` arrays.
    * **uint16 r² quantisation** — 65,536 levels, precision
      ~1.5e-5 (well below LD estimation noise at ~0.01).  Halves
      the ``data`` array; auto-dequantised to float32 on load.
    * **zstd container (``.ldz``)** — optional custom binary
      format written via the ``zstandard`` package.  Default
      level 19 yields ~15% additional shrink over zlib but
      ~13× slower to write; users can dial down with
      ``compression_level=9`` for speed.  Zlib remains the
      default because at moderate zstd levels it matches or
      beats zstd on numeric LD arrays.

  Measured on chr22 1401-AFR reference (101,641 SNPs,
  786,274 edges):

  ============================================  =======  =========
  variant                                       size     save
  ============================================  =======  =========
  default (dense + float32 + zlib)              3108 KB  294 ms
  compact (symm + uint16 + zlib) — recommended  1702 KB  304 ms
  compact + zstd L19 — max compression          1464 KB  3.9 s
  ============================================  =======  =========

  Scaled to genome-wide ~25 M-edge references, ``compact=True``
  is expected to drop a ~350 MB .npz to ~130 MB with no query-
  time cost (in-memory representation stays dense, expanded on
  load).  Regression:
  ``test_ld_graph_compact_roundtrip_preserves_queries``.

- **LDGraph genome-build metadata + LD-clumping sanity check.**
  ``LDGraph`` now carries a ``build`` attribute (``hg18`` /
  ``hg19`` / ``hg38``, normalised from user-supplied aliases like
  ``GRCh37`` or ``b37``).  ``LDGraph.from_plink_ld(path,
  build='hg19')`` declares it at build time; the value is
  persisted in the ``.npz`` / ``.ldz`` and its JSON sidecar so
  pre-built graphs travel with their build.

  ``pycmplot.io.load()`` uses this to guard the LD-clumping step
  against silent build mismatches — a common failure mode when a
  hg19 LD reference meets hg38 sumstats (or GeneHancer auto-lifts
  hg19 sumstats to hg38 before clumping).  If the LD reference
  build and the track's effective post-liftover build disagree,
  the loader raises with a clear message rather than degrading to
  100% distance-only fallback::

      ValueError: LD clumping build mismatch for track 'HbF':
        LD reference build = hg19
        track effective build after liftover = ['hg38']
      Options: (a) rebuild the LD reference in the track's build,
      (b) provide a same-build LD reference, or (c) disable
      GeneHancer (--no_genehancer) so the track stays in native
      build.

  Falls back to a warning + skipped check when the LD reference
  is undeclared (older graphs with no build metadata; use
  ``ld_reference_build='hg19'`` on the loader as an override).

  New CLI flag ``-ldrb`` / ``--ld_reference_build`` declares the
  build when the user is pointing at a raw PLINK ``.ld`` file
  (rather than a pre-built LDGraph ``.npz`` / ``.ldz`` that
  already knows).  Regression:
  ``test_ld_graph_build_metadata_roundtrip``.

- **``pycmplot.clump()`` — LD-based greedy clumping with an
  external reference panel.**  Optional standalone alternative to
  the built-in distance-based
  :func:`~pycmplot.stats.get_lead_snps`.  Accepts a PLINK
  ``--r2`` (``.ld``) file path or a pre-parsed
  ``{snp_a: {snp_b: r²}}`` dict.  Greedy algorithm matches PLINK
  ``--clump``: within a ``kb``-distance window, a candidate is
  clumped away only if its r² with any already-accepted lead is
  ``>= r2``.  Variants absent from the reference are treated as
  r²=0 (independent), matching PLINK's default behaviour.

  Usage::

      import pycmplot
      leads = pycmplot.clump(
          df_of_significant_snps,
          ld_reference="ref.ld",   # PLINK .ld output
          r2=0.1,
          kb=500,
          p_threshold=5e-8,
      )

  The distance-only fallback (``ld_reference=None``) makes this
  drop-in-safe: it degrades cleanly to independent-if-not-linked
  semantics without failing.  Not yet wired into the plotting
  pipeline (users invoke it directly to build a lead-SNP table
  before calling :func:`pycmplot.io.load`); pipeline integration
  is a follow-up.  Regressions:
  ``test_ld_clump_greedy_semantics``,
  ``test_ld_clump_no_ld_reference_falls_back_to_distance``.

- **``pycmplot.clump()`` distance-based fallback for missing
  reference overlap.**  Previously, when a variant was absent
  from the LD reference, ``r2()`` returned 0.0 → the pair was
  silently treated as "independent".  On a low-overlap reference
  this over-selected independent leads without warning.

  New behaviour: when either the lead or candidate isn't in the
  LD reference, ``clump()`` falls back to distance-only inside the
  same ``kb`` window — biology-safe conservative default
  (unknown LD → assume possibly linked → clump).  Variants
  beyond the ``kb`` window from every lead stay independent
  regardless of reference status, so the behaviour matches
  :func:`pycmplot.stats.get_lead_snps` exactly when
  ``ld_reference=None``.

  The clumping run logs a summary at INFO level::

      clump: 42,318 pair checks, 3,241 (7.7%) fell back to
      distance-only because at least one variant was absent
      from the LD reference

  A WARNING fires when the fallback fraction exceeds 20%,
  pointing users at reference / summary-stats naming
  reconciliation as a likely cause.  Regression:
  ``test_ld_clump_missing_snp_falls_back_to_distance``.

- **``pycmplot.VariantMatcher`` — cross-source variant identity
  resolver.**  New class that reconciles SNP identifiers across
  formats:

    * ``rs12345`` — dbSNP rsID
    * ``1:636975:A:G`` — CHR:POS:REF:ALT (Ensembl / PLINK-modern)
    * ``chr1:636975:A:G`` — same with UCSC ``chr`` prefix
    * ``1:636975:G:A`` — allele-order-agnostic
    * ``1:636975`` — position-only (allele-blind fallback)
    * mixed-separator forms with ``:``, ``_``, ``/``, ``|``
      (e.g. ``chr1_63748574_A_G``, ``1:63748574_G_A``,
      ``1_63748574:A|G``) as seen in bcftools-lifted files,
      TOPMed-style dbSNP dumps, PGS-Catalog exports, and various
      genotyping pipelines.  Deliberately excludes ``-`` and ``.``
      from the separator set to avoid ambiguity with deletion
      sentinels and version dots in allele strings.
    * chromosome aliases across PLINK numeric / UCSC ``chr`` /
      Ensembl bare conventions:
      ``X ↔ 23 ↔ chrX ↔ chr23``, ``Y ↔ 24 ↔ chrY ↔ chr24``,
      ``M ↔ MT ↔ 26 ↔ chrM ↔ chrMT``, plus ``XY ↔ 25`` for the
      pseudo-autosomal region.  So ``chr23:5737583`` and
      ``chrX:5737583`` resolve to the same variant.
    * GWAS-Catalog / meta-analysis style ``CHR:POS:TYPE`` where
      TYPE is a variant-class token
      (``SNV`` / ``SNP`` / ``MNV`` / ``MNP`` / ``INDEL`` / ``DEL``
      / ``INS`` / ``CNV`` / ``DELINS`` / ``SUB``), case-
      insensitive.  Class tokens are treated as position-only
      aliases: ``1:7535638:SNV``, ``1:7535638:SNP``,
      ``1:7535638:INDEL`` all resolve to the same canonical
      variant as ``1:7535638:A:G`` (regardless of registration
      order, thanks to the position-only-record → full-alleles
      upgrade path).

  Canonical key is ``(CHR, POS, frozenset({REF, ALT}))`` when
  alleles are known, else ``(CHR, POS)``.  Builds a bidirectional
  index (alias → canonical, plus (chrom, pos) → canonical) so
  every registered alias — including deliberately-swapped alleles
  — resolves to the same variant record.

  Integration points:

    * ``pycmplot.LDGraph.remap_ids(rename)`` — in-place relabel of
      graph nodes using a matcher-produced rename table.
    * ``pycmplot.clump(..., variant_matcher=vm)`` — auto-relabels
      the LD graph before greedy clumping so cross-source naming
      differences don't force the missing-SNP fallback.
    * ``VariantMatcher.from_dataframe(df)`` — builds a matcher
      from a summary-stats DataFrame in one line.

  Regressions: ``test_variant_matcher_basic_forms``,
  ``test_variant_matcher_rename_map``,
  ``test_variant_matcher_clump_integration``.

- **CLI: LD-based clumping + variant harmonisation.**  Three new
  flags plumb LD-based clumping and cross-source ID matching all
  the way from the shell:

    * ``-ldr`` / ``--ld_reference PATH`` — path to a PLINK ``.ld``
      file or an ``LDGraph`` ``.npz`` / ``.ldz`` (auto-loaded).
      When set, ``pycmplot`` switches from distance-only clumping
      to :func:`pycmplot.stats.clump`.
    * ``-ldr2`` / ``--ld_r2 FLOAT`` — r² threshold (default 0.1,
      PLINK ``--clump-r2`` convention).
    * ``--harmonize_variants`` — build a ``VariantMatcher`` from
      the loaded sumstats and relabel LD-graph SNP IDs to match.
      Turn on when your sumstats use rsIDs and the LD reference
      uses ``CHR:POS:REF:ALT`` (or vice versa).

- **UCSC functional annotation tracks + LD-block-aware lead
  selection.**  Five new bundled tracks (hg38) drive the priority
  score and the tiebreaker used when perfect-LD variants share the
  same P-value:

    * ``gtex_eqtl_caviar.tsv.gz`` (4.4 MB) — GTEx CAVIAR
      fine-mapped eQTL variant→gene→CPP.  Adds per-gene
      ``eqtl_bonus = max(CPP) × 3`` to the priority score; a
      CPP-1.0 fine-mapped eQTL contributes +3 to that gene,
      dominating all other components.  The strongest possible
      signal that a variant regulates a specific target.
    * ``ucsc_ccre.tsv.gz`` (8.9 MB) — ENCODE candidate cis-
      Regulatory Elements (pELS / dELS / PLS / CTCF-bound).
    * ``ucsc_dnase.tsv.gz`` (5.6 MB, score ≥ 250 filter) — DNase
      hypersensitivity clusters (open chromatin).
    * ``ucsc_tfbs.tsv.gz`` (33 MB, score ≥ 500 filter) — TF
      binding sites (with TF name).
    * ``ucsc_cpg.tsv.gz`` (259 KB) — CpG islands.

  Total bundled: ~52 MB.  Preparation script
  ``scripts/prep_functional_tracks.py`` converts UCSC raw
  downloads to bundled derivatives with per-track score filtering
  so the ``absence == no evidence`` semantic is clean.

  New public helper :func:`pycmplot.annotation.position_informativeness_score`
  aggregates evidence across all tracks into a single per-position
  score used exclusively for tiebreaking within LD blocks.  cCRE,
  CpG, DNase, and TFBS are used only for this tiebreaker; they
  don't directly contribute to the per-gene priority score
  (position-only tracks don't tell us *which* gene a SNP regulates).

- **GeneHancer enhancer-gene interactions in the priority score.**
  New ``genehancer_bonus`` term in
  :func:`~pycmplot.annotation._annotate_variant`, computed from the
  bundled ``genehancer.tsv.gz`` (218,117 elements, 818,358 gene
  connections; hg38).  For each SNP, GH element(s) that *contain*
  the SNP position (``start <= pos <= end``) contribute
  ``score / 5`` (uncapped, no distance decay) to every connected
  gene of that element.  A GH connection score of 20 → +4.0 bonus
  (rare, dominates all other components); score 10 → +2.0
  (comparable to ``genic_bonus``); score 1 → +0.2 (nudge only).
  Multiple containing elements linking the same gene collapse to
  the ``max`` bonus.  Strict overlap: SNPs outside every GH element
  receive zero bonus and the algorithm falls back to
  geometry+strand+biotype.

  Priority formula becomes::

      priority = biotype_weight
                 · (distance_score + 2·genic + 1·upstream
                    + 2·promoter + genehancer_bonus)

  Hits table gains two columns: ``genehancer_bonus`` (the bonus
  actually applied to the winning gene, 0.0 if none) and
  ``gh_score`` (the raw GH connection score, ``None`` if there was
  no GH support).  Turn off with ``-no_gh`` / ``--no_genehancer``
  to fall back to pure geometry+strand+biotype scoring.

  Solves the AC110792.1 / PDGFRA case at
  chr4:54,768,663 (hg19-lifted to hg38): the region is a validated
  GeneHancer-predicted PDGFRA enhancer, and GH-based scoring
  correctly overrides the tiny orphan Ensembl PC in favour of
  PDGFRA when the SNP falls inside a GH element.  Regressions:
  ``test_genehancer_bonus_boosts_gh_connected_gene`` (verifies
  bonus math + PDGFRA lift) and
  ``test_genehancer_no_op_when_position_uncovered`` (identity
  fallback for GH-uncovered positions).

  Custom GeneHancer paths via ``PYCMPLOT_GENEHANCER_HG38``
  environment variable or :class:`~pycmplot.resources.ResourceConfig`
  ``genehancer_hg38`` attribute.  Preparation script
  ``scripts/prep_genehancer.py`` converts the source semicolon-CSV
  into the bundled TSV format.

  **hg19 auto-liftover under GH:** since GH is hg38-only, a
  pure-hg19 loader group used to skip the liftover step entirely
  (neither ``_needs_lift_file`` nor ``_needs_lift_group`` fired
  without mixed builds).  A third liftover trigger,
  ``_needs_lift_gh``, now fires when ``use_genehancer=True`` and
  the track's declared build is hg18/hg19, so GH coordinates always
  match the annotated variant positions.  Turn off with
  ``--no_genehancer`` if you want to keep native-build hg19
  coordinates.

- **Two-layer annotation architecture (2026-09-15).**  Clumping
  and functional prioritisation are now cleanly decoupled:

    * **Layer 1 — Independent lead identification** (statistical):
      :func:`~pycmplot.stats.get_lead_snps` /
      :func:`~pycmplot.stats.get_highlight_snps` /
      :func:`pycmplot.stats.clump` operate on P and distance / LD
      alone.  Which variant is called a "lead" depends only on
      statistics and reference-panel LD, never on the annotation
      panel.  Reproducible across annotation-source changes.
    * **Layer 2 — Annotation representative selection**
      (functional): after Layer 1, ``pycmplot.io.load`` scans each
      lead's LD-block members (variants tied at the same P within
      ``clump_window_kb`` on the same chromosome) for the position
      with the highest
      :func:`~pycmplot.annotation.position_informativeness_score`.
      That position drives the ``_annotate_variant`` call; the
      statistical lead's rsID / POS / P remain the record's
      identity.  New ``annot_pos`` / ``annot_snp`` columns record
      the representative when it differs from the lead.

  ``get_lead_snps`` retains its optional ``tiebreak_score_col``
  parameter for direct-API callers who want mixed statistical +
  functional clumping in one pass, but the pipeline no longer
  uses it — the two layers are separate.

  Motivating example: chr4 hg19 54,768,663 (rs111878819) — one
  of three tied-significant variants (rs111341649 at 54,774,323
  and rs6833117 at 54,775,640).  Pre-fix greedy pick chose the
  leftmost (rs111878819) as lead → geometric fallback to
  ``LNX1-GSX2``.  With decoupled layers, rs111878819 remains the
  statistical lead (deterministic POS-sort), but ``annot_pos``
  is swapped to rs111341649 (inside a PDGFRA-linked GH element
  at 53,907,897-53,909,741 hg38) → annotation returns
  ``LNX1-PDGFRA``.  The hits table now reports both:
  ``SNP=rs111878819`` (statistical lead) and
  ``annot_snp=rs111341649`` (annotation representative).
  Regression: ``test_ld_block_tiebreak_prefers_functional_variant``.

- **Distance-based clumping window is now a first-class knob
  (``-cw`` / ``--clump_window_kb``, default 250 kb).**  Previously
  ``get_hits_summary_table`` piped its ``window_kb`` argument into
  both the gene-candidate search AND
  :func:`~pycmplot.annotation._clump_by_distance` (independent-lead
  selector), so bumping the annotation window silently dropped
  nearby leads.  Clumping is now decoupled from annotation search
  and exposed as its own parameter with a **new default of 250 kb**,
  matching PLINK's ``--clump-kb`` convention and aligning with
  typical European-population LD extent (~100-200 kb) plus a
  modest safety margin.  Old default was 500 kb.

  The primary lead-selection clumping in
  :func:`pycmplot.stats.get_lead_snps` (and
  :func:`~pycmplot.stats.get_highlight_snps`) also drops from
  500 kb to 250 kb so that the primary pass and the defensive
  second pass agree.  African-ancestry or admixed cohorts where
  LD is shorter (~30-50 kb) can go tighter (``-cw 100``);
  conservative locus-definition workflows can widen
  (``-cw 500``).  Cache invalidates automatically when
  ``clump_window_kb`` changes.  Regression:
  ``test_clumping_window_decoupled_from_annotation_window``.

- **Configurable annotation search window (``-aw`` /
  ``--annotation_window_kb``).**  The pipeline was hard-coded to a
  2 Mb window in :func:`~pycmplot.io.load`, which diverged from the
  500 kb value stated in ``gene_selection_algorithm.md``.  Search
  radius is now a top-level knob, defaults to **500** kb (matching
  the spec), and cascades through :func:`~pycmplot.annotation.get_hits_summary_table`
  and the hits-overlay cache key so the auto-hits table regenerates
  cleanly when the window changes.  Wider windows can be requested
  per-run (e.g. ``-aw 2000``) when the local neighbourhood is
  gene-poor and TAD-scale PC candidates are wanted.

  Motivating case: ``chr6:115443568`` (hg38) at the pipeline's old
  2 Mb default resolved to ``HS3ST5-NT5DC1`` (an intergenic_PC
  join across 1.1 Mb + 0.66 Mb PC flankers).  With the 500 kb
  spec-aligned default it resolves to ``FRK`` (the closest PC at
  487 kb on one side).

- **Custom highlight-legend size and padding.**  Two new
  user-facing knobs on both plotters:

    * ``highlight_legend_size`` (``-hll_size`` /
      ``--highlight_legend_size``) — overrides the legend font
      size for both entries and title.  Defaults preserve
      the pre-existing look: ``annotation_size`` on the linear
      plot, ``track_label_size`` on the circular plot.
    * ``highlight_legend_pad`` (``-hll_pad`` /
      ``--highlight_legend_pad``) — padding between the legend
      frame and the axes edge, passed to matplotlib as
      ``borderaxespad``.  Default: matplotlib's built-in 0.5.
      Increase to push the legend further from the plot;
      negative values allowed for slight overlap.

  Both threaded through :func:`~pycmplot.plotting.linear.plot_linearm`,
  :func:`~pycmplot.plotting.linear.plot_linear`, and
  :func:`~pycmplot.plotting.circular.circular`.

- **Liftover-unmapped report: one file for all tracks.**  Variants that
  cannot be placed in hg38 were previously dropped silently.  The loader
  now writes them to ``<output prefix>.liftover_unmapped.tsv`` (next to
  the hits table; ``pycmplot.liftover_unmapped.tsv`` when no prefix is
  set).  Columns: ``LABEL`` (track) first, then the original-build
  coordinates and fields, plus ``UNMAPPED_REASON``, ``LIFTED_CHR``,
  ``LIFTED_POS`` and ``SIGNIFICANT``.  Reasons are
  ``no_chain_mapping``, ``maps_to_other_chromosome`` (includes alt
  contigs such as ``chr8_KI270821v1_alt``) and
  ``beyond_hg38_chromosome_length``.  Rows are grouped by track in the
  order the tracks were given.  Tracks served from the cache do not
  re-run liftover, so their rows are carried over from the existing
  file; if that file has been moved or deleted a warning names the
  affected tracks.  Each track that is lifted over in the current run
  logs its count per reason at INFO; cached tracks do not repeat it.

- **Warning when significant variants are lost in liftover.**  If any
  unmapped variant passes the loader's significance threshold
  (``P <= threshold``, or ``|P| >= threshold`` for signed statistics),
  a WARNING lists up to 10 of them and states that they are missing
  from the plot and the hits table.  The advice depends on why
  liftover ran: when it was triggered only because GeneHancer needs
  hg38, it suggests re-running with ``--no_genehancer`` /
  ``use_genehancer=False``; for mixed builds, it suggests plotting the
  builds as separate plots in their native coordinates.  The summary
  is stored in the track's cache entry, so this warning is repeated on
  cached runs.

- **``liftover_position(..., return_unmapped=True)``.**  Returns
  ``(clean_df, unmapped_df)`` with the columns described above.  The
  default (``False``) still returns a single DataFrame.

- **``gene_label`` and ``annot_window_kb`` columns in the hits table.**
  ``gene_label`` holds each hit's gene label (the containing gene,
  ``nearest_gene``, for genic hits; ``top_gene`` otherwise; the rsID
  when neither exists).  ``annot_window_kb`` records the annotation
  search window.  Both are used by the plot-time gene-label collapse
  described under **Changed**.

- **``pycmplot.annotation.prepare_annotation_labels``** (and
  ``labels_to_draw``).  Resolves the label column for a plot and marks
  which hits get a label (``label_shown``).  Used by both
  :func:`~pycmplot.plotting.linear.linear` and
  :func:`~pycmplot.plotting.circular.circular`; the input table is not
  modified.

- **Plot-time, per-locus highlight filter** (``-pht`` /
  ``--plot_highlight_thresh``; ``highlight_thresh`` on
  :func:`~pycmplot.plotting.linear.linear` and
  :func:`~pycmplot.plotting.circular.circular`).  Applied to each track
  separately: a highlighted variant keeps its highlight only if it lies
  within the locus window of a highlighted variant in the same track
  that passes the cutoff (``P <= threshold``, or ``|P| >= threshold``
  for signed statistics).  The window is the loader's
  ``clump_window_kb`` (``-cw``), the same window that defines the
  highlighted loci, so each locus is judged by its own lead: passing
  loci keep their whole highlight, failing loci are drawn in the
  background colour, and if nothing passes nothing is highlighted.  A
  locus merged into another track's lead in the hits table is still
  judged by its own track's lead.  The loader records
  ``clump_window_kb`` in the hits table and the CLI passes ``-cw``
  explicitly; Python users can override it with
  ``highlight_window_kb``.  The loaded data and the gene labels are
  unaffected (use ``-psig`` / ``--plot_signif_threshold`` for labels),
  so a run can load broadly and choose what to highlight at plot time.
  New helpers :func:`pycmplot.annotation.filter_in_locus_by_threshold`
  and :func:`pycmplot.annotation.resolve_highlight_window_kb`.  (0.4.1
  removed an unused ``highlight_thresh`` parameter from
  :func:`~pycmplot.plotting.circular.circular`; the name now returns
  with this meaning.)

- **Functional-annotation columns in the hits table and overlay.**
  With GeneHancer enabled (the default), each lead's hits-table and
  ``hits.<group>.tsv`` row now reports the functional evidence at the
  variant used for its gene assignment (``annot_pos`` when the
  LD-block tie-break chose a representative, else the lead):
  ``gh_id`` / ``gh_feature`` (GeneHancer element), ``eqtl`` (GTEx
  CAVIAR fine-mapped eQTLs as ``GENE_TISSUE``, highest CPP first),
  ``tfbs`` (ENCODE TF clusters as ``TF_SCORE``), ``ccre_id`` / ``ccre``
  / ``ccre_class`` (ENCODE cCRE accession, signature description and
  class), ``dnase_score`` / ``dnase_sources`` (ENCODE DNase cluster)
  and ``cpg_island``.  Multi-valued fields are comma-separated;
  ``None`` means no overlap.  The columns are reporting only; the same
  tracks also drive the LD-block tie-break (see **Fixed**).  New helpers
  :func:`pycmplot.annotation.load_functional_tracks` and
  :func:`pycmplot.annotation.annotate_functional`.

- **Bundled functional tracks rebuilt from the raw downloads**
  with the new ``scripts/prep_functional_tracks.py`` (sources in
  ``annotation/``).  Every column used for scoring is identical to the
  previous files (verified row by row for all six tracks); the files
  gain ``GH_ID`` (GeneHancer), ``CCRE_ID`` / ``CCRE_CLASS`` /
  ``DESCRIPTION`` (cCRE), ``SOURCES`` (DNase) and ``NAME`` /
  ``CPG_COUNT`` (CpG).  ``gtex_eqtl_caviar.tsv.gz`` now holds one row
  per variant-gene-tissue (adding ``SNP`` and ``TISSUE``, CPP to 6
  decimals) instead of one per variant-gene;
  :func:`~pycmplot.annotation.load_eqtl` derives the scoring value
  (maximum CPP across tissues, rounded to 4 decimals), which
  reproduces the previous file exactly, and still accepts a per-gene
  file without ``TISSUE`` supplied through ``PYCMPLOT_EQTL_HG38``.
  Bundled data grows by about 15 MB.  The hits-overlay cache fingerprint now also covers the
  GeneHancer, eQTL and ENCODE / CpG track files, so swapping any of
  them regenerates cached hits tables (previously only the gene-info
  files were covered).

**Changed**

- **Switched ``pyliftover`` dependency to ``liftover`` for performance.**
  Prior versions used ``pyliftover`` to convert hg18/hg19 coordinates to 
  hg38. ``liftover`` appears to offer better performance with its c++ 
  implementation, and it is more actively maintained. Notably, both tools 
  produce identical results, with a ``liftover`` substantially faster for 
  large files.

  Verified on 200,000 random positions per chain file (hg19→hg38 and
  hg18→hg38): identical output for every position.  Chain files load
  16–35× faster and single lookups run 2–3× faster.  ``pyliftover``
  is no longer imported anywhere and was removed from the
  ``pyproject.toml`` and ``setup.cfg`` dependencies; the README and the
  :mod:`pycmplot.liftover` docstrings now refer to
  ``liftover.ChainFile``.  :func:`~pycmplot.liftover.liftover_position`
  also resolves each chain file once per call instead of once per row,
  making it about 3× faster overall (9.6 s → 3.4 s on 1 M hg19 rows).

- **Loader runtime: heavy resources hoisted out of the per-track loop.**
  Prior to this release, every non-cached track that requested LD-based
  clumping re-loaded the LD reference graph from disk inside its own
  iteration of :func:`pycmplot.load`.  For a 31 M-pair reference panel
  (~30 s parse), an N-track run paid that fixed cost N times.  The same
  N-multiplication applied to GeneHancer + eQTL bundle loads inside the
  Layer-2 annotation-representative swap when ``use_genehancer=True``.

  The loader now installs two lazy-load closures above the per-track
  loop.  The first non-cached track that needs each resource pays the
  parse; every subsequent iteration reuses the cached object.  A
  fully-cached run (all tracks hit the on-disk cache and ``continue``)
  never triggers a load at all — the closures stay dormant.  For a
  five-track run against a full 1000G reference the wall-clock
  improvement is roughly an order of magnitude on the LD-clumping
  branch.

  Semantics are otherwise unchanged: ``in_locus`` is still computed
  per-track, and the per-track effective-
  build sanity check (which can differ between tracks after per-track
  liftover) still fires inside the loop against the graph's declared
  ``build`` metadata.  The build check itself no longer costs a full
  graph parse — it consults the metadata carried on the already-loaded
  object.

  One semantic difference: lead extraction now runs **once, after the
  per-track loop**, on the significant variants pooled from every
  track (Layer 1 clumping, then the Layer-2 annotation-representative
  swap).  A locus significant in several tracks therefore yields one
  lead, the most significant across tracks.  The plotters do not use
  the track of origin when drawing hits, so plots are unaffected apart
  from the fixes listed under **Fixed** below.

- **``pycmplot.stats.get_signif_snps`` (public API).**  New helper that
  extracts significant variants from a summary-statistics DataFrame
  without applying clumping or annotation-representative selection.
  Preserves the calling convention of :func:`~pycmplot.stats.clump` and
  :func:`~pycmplot.stats.get_lead_snps` (``logp`` / ``score_col`` /
  ``ascending`` auto-detection for signed statistics) and returns a
  minimal DataFrame with only pipeline-relevant columns (``CHR``,
  ``POS``, ``SNP``, ``P``, ``logP``, ``OLD_POS``, ``OLD_BUILD``,
  ``BUILD``, ``LABEL``).  Exported from the top-level ``pycmplot``
  namespace.  Regression:
  ``benchmark/tests/test_priority_algorithm.py::test_ld_reference_loaded_once_across_tracks``.

- **``pycmplot.stats.clump`` accepts a live ``LDGraph`` object.**  The
  ``ld_reference`` parameter's type hint now includes ``LDGraph``
  alongside ``str`` / ``dict`` / ``None``.  Runtime dispatch already
  handled this case; the annotation update makes it discoverable and
  is what the hoisted loader uses to hand the pre-loaded graph to each
  per-track ``clump`` call without re-parsing from disk.

- **Biologically-meaningful gene-assignment algorithm.**  The
  priority score in :func:`pycmplot.annotation._annotate_variant`
  was reworked to encode the design in ``gene_selection_algorithm.md``:

  * ``distance_score = exp(-|d| / 100 kb)`` — nearby genes dominate
    more sharply than the previous ``1/log10(d+10)`` decay.
  * Explicit strand-aware ``upstream_bonus`` when the SNP sits 5'
    of the gene body (TSS side), equal to ``distance_score`` (at most
    +1.0; see *Upstream bonus decays with distance* below).
    Upstream = ``pos < start`` for ``+`` strand, ``pos > end`` for
    ``-`` strand; genic positions are neither upstream nor downstream.
  * ``genic_bonus`` (+2.0) for containing gene body.
  * ``promoter_bonus`` (+2.0) when ``|pos - TSS| <= 2 kb``, strand-
    aware (TSS = START for ``+``, END for ``-``).  Fires regardless
    of whether the SNP is inside the gene body.
  * The whole component sum is **multiplied by ``biotype_weight``**,
    which is the load-bearing choice: a distant protein-coding gene
    now beats a nearby pseudogene the SNP happens to sit inside,
    and an antisense/divergent transcript containing the SNP loses
    to a promoter-proximal PC gene 500 bp away.  ``top_gene`` is
    always the highest-scoring candidate — no early-out on
    "containing gene" or on the intergenic "LEFT-RIGHT" pair.

  New fields on every hits-table row: ``snp_position``
  (``"genic"|"upstream"|"downstream"``), ``upstream_of_gene``,
  ``promoter_proximal``, ``dist_to_tss``, ``upstream_bonus``,
  ``genic_bonus``, ``strand`` (of the winning gene).  Regression:
  ``benchmark/tests/test_priority_algorithm.py``.

- **PC-preference selection in ``_annotate_variant``.**  ``top_gene``
  is no longer a pure argmax over the priority score.  A cascade
  runs against the candidate pool:

    1. No candidates in window -> ``top_gene = None``.
    2. No PC-like candidate in window -> highest-priority candidate
       overall (may be a lncRNA or pseudogene; nothing better
       exists in the region).
    3. A PC-like candidate contains the SNP -> highest-priority
       PC-like *genic* candidate wins.
    4. Otherwise -> take the highest-priority PC-like candidate on
       each positional side of the SNP.  If both sides have one,
       report as joined ``"LEFT-RIGHT"`` label with
       ``biotype = "intergenic_PC"``; if only one side has a PC,
       report that single PC.

  ``PC_LIKE_BIOTYPES`` (module constant in ``pycmplot.annotation``)
  defines the high-confidence coding set: ``protein_coding``,
  ``protein_coding_LoF``, ``protein_coding_CDS_not_defined``,
  ``nonsense_mediated_decay``, and the IG/TR gene biotypes.  Rules
  out lncRNAs (including divergent transcripts like ``*-DT``),
  pseudogenes, antisense, sense_intronic/overlapping, and misc RNAs
  as ``top_gene`` when a PC candidate exists.  Motivating case:
  ``chr6:135820013`` (hg38), where the SNP sits inside ``AHI1-DT``
  but is 31 kb upstream of ``PDE7B``; the cascade now correctly
  returns ``PDE7B`` as ``top_gene``.  Regression
  ``benchmark/tests/test_priority_algorithm.py::test_ahi1dt_pde7b_regression_hg38``.

- **``lncRNA`` / ``lincRNA`` / ``ncRNA`` biotype weights lowered
  from 0.70 to 0.55.**  Complements the PC-preference cascade: when
  no PC candidate exists in the window and the highest-priority
  non-PC candidate wins, the lower weight keeps lncRNAs from
  dominating over closer pseudogenes / antisense in edge cases.
  ``3prime_overlapping_ncRNA`` and ``macro_lncRNA`` dropped from
  0.65 to 0.60; anti-sense biotypes kept at 0.65 because they
  carry more information about the sense PC gene than generic
  lncRNAs do (regulatory antisense RNAs such as PDE7B-AS1).

- **Ensembl 116 GRCh38 + Ensembl 87 GRCh37 geneinfo aliases in
  ``BIOTYPE_WEIGHTS``.**  Added ``3prime_overlapping_ncrna`` (GRCh37
  lowercase spelling, matches ``3prime_overlapping_ncRNA`` at 0.65)
  and ``rRNA_pseudogene`` (GRCh38.116, matches pseudogene subtypes
  at 0.20) so the full vocabulary of both lightweight reference
  files gets a defined weight.

- **CLI ``--ld_reference`` now requires a pre-built LDGraph.**
  Raw PLINK ``.ld`` paths are rejected at :func:`pycmplot.io.load`
  with instructions to convert to an ``.npz`` / ``.ldz`` first
  via ``LDGraph.from_plink_ld(path, build='hg19').save('ref.npz')``.
  Rationale: the pipeline's build sanity check needs a genome
  build to compare against, and a raw ``.ld`` file has no place
  to declare one.  Making the CLI insist on pre-built graphs
  forces build declaration exactly once — at graph-build time —
  and the metadata travels with the file thereafter.  Removes
  the ``ld_reference_build`` parameter on
  :func:`pycmplot.io.load` and the corresponding
  ``-ldrb`` / ``--ld_reference_build`` CLI flag (introduced
  earlier in this same release).  ``pycmplot.stats.clump()``
  Python API still accepts raw ``.ld`` paths for scripted use;
  only the loader-driven pipeline is stricter.

- **Per-track cache key now covers LD clumping + clump window.**
  On a cache hit the cached DataFrame + leads bypass the entire
  clumping path — including the LD lookup — so flipping
  ``--ld_reference``, ``--ld_r2``, ``--harmonize_variants``, or
  ``--clump_window_kb`` between runs previously served silently
  stale distance-clumped leads.  These four params are now
  included in the per-track cache key so any change invalidates
  the affected entries automatically.  LD params only contribute
  to the key when ``ld_reference`` is set, so users on the
  distance-only path keep their historical cache behaviour
  unchanged.  Regression:
  ``test_cache_key_includes_ld_and_clump_params``.

- **Upstream bonus decays with distance.**  During 0.4.3 development
  ``upstream_bonus`` was a flat +1.0 for any gene whose 5' end faced
  the SNP anywhere in the 500 kb window, so a distant gene could
  outrank an adjacent one.  Example (hg38): a SNP at chr2:60,450,000,
  520 bp outside BCL11A (− strand, SNP at its 3' end), was assigned
  PAPOLG, 306 kb away (score 1.047 vs 0.995).  The bonus is now
  ``distance_score`` (``exp(-d / 100 kb)``: 1.0 adjacent, 0.37 at
  100 kb, 0.05 at 300 kb), and the locus is assigned BCL11A.  Compared
  with the flat bonus, on 20,000 random hg38 positions ``top_gene``
  changed at 14.6% of those with a gene in the window; 2,693 of the
  2,701 changes moved to a closer gene, and the median distance of
  the chosen gene(s) fell from 333 kb to 99 kb.  The remaining 8 are gene deserts with no
  protein-coding gene in range.  HBB, HBG2, HBS1L, HBS1L-MYB,
  NPRL3-HBA2 and intragenic BCL11A assignments are unchanged.  A
  protein-coding gene 20 kb upstream still beats a pseudogene
  containing the SNP (1.64 vs 0.6).  ``upstream_of_gene`` is now set
  from ``snp_position`` rather than from a non-zero bonus.  **Hits
  tables regenerated with this release may label some loci
  differently from earlier runs.**

- **Repeated gene labels are collapsed at plot time, not removed from
  the hits table.**  A hits-table step added during 0.4.3 development
  that dropped rows sharing ``(CHR, nearest_gene)`` is removed.  It compared ``nearest_gene``
  although intergenic hits are labelled with ``top_gene``, treated all
  hits with no gene on a chromosome as duplicates, merged independent
  signals at any distance, kept the leftmost row rather than the most
  significant one, and removed rows that SNP-labelled plots needed.
  The hits table now keeps every independent lead.  When a plot is
  labelled by gene, hits sharing a gene label within the annotation
  window (default 500 kb) on the same chromosome get one label, on the
  most significant hit; the others keep their highlighted points,
  guide lines and per-locus colours.  Hits with no gene label are
  never merged.  SNP-labelled (and other non-gene) plots label every
  lead.

- **Hits-overlay cache key is versioned** (``HITS_LOGIC_VERSION = 4``
  in :mod:`pycmplot.cache`).  Cached ``hits.<group>.tsv`` overlays
  built by earlier hits-table logic regenerate once; rows added by
  hand (``source="user"``) are kept.

**Fixed**

- **(CHR, SNP) duplicate collapse across annotation, cached
  overlay, and circular plot.**  Same rsID at slightly different
  POS values across build-mixed tracks (e.g. ``rs123`` in hg19
  vs hg38, or pre- vs post-liftover) was surviving
  :func:`~pycmplot.annotation._clump_by_distance` because POS
  separation exceeded the window.  The duplicates then
  propagated to the cached ``hits.<group>.tsv`` overlay, the
  circular plot (which consumes the overlay directly), and the
  linear plot's SNP-labelled annotation track.  Fixed at three
  layers:

    1. :func:`~pycmplot.annotation.get_hits_summary_table` now
       drops ``(CHR, SNP)`` duplicates at write time so fresh
       overlays are always clean.
    2. :func:`pycmplot.cache.read_hits_overlay` dedups on read
       so users with pre-fix overlays don't need to invalidate
       their cache.
    3. :func:`~pycmplot.plotting.linear.plot_linearm` picks the
       dedup key based on ``label_col``: ``(CHR, SNP)`` when
       annotating by SNP; ``(CHR, POS, label)`` when annotating
       by gene / custom column (so distinct variants at the
       same locus still keep separate labels).

  Regressions:
  ``test_hits_overlay_dedup_by_snp_at_source``,
  ``test_hits_overlay_dedup_on_read_defends_pre_fix_overlays``,
  ``test_annotation_dedup_by_snp_when_snp_labeled``.

- **Gene annotation missing for all but the last track.**  After the
  loader runtime hoist above, the hits table was built from a stale
  per-track variable holding only the last freshly loaded track's
  significant variants, and the pooled, clumped leads were discarded.
  Loci from every other track got no label; the hits table was empty
  (no labels at all) whenever the last track had no significant
  variants; and a run where every track came from the cache failed
  with ``KeyError: 'CHR'``.  The hits table is now built from the
  pooled leads of all tracks, and tracks served from the cache feed
  the same pooled clumping as freshly loaded ones.  Affected
  unreleased 0.4.3 code only (0.4.2 collected leads per track).

- **Layer-2 annotation-representative swap scanned the wrong track.**
  Tied-P neighbours for each lead were looked up in the last freshly
  loaded track's data instead of the lead's own track.  Each lead is
  now matched against its own track.

- **No significant variants in any track raised ``KeyError: 'CHR'``.**
  The loader now returns an empty hits table.

- **``annotate="GENE"`` chose one label column for all hits from the
  last hit.**  :func:`~pycmplot.annotation.get_annotation_column`
  re-chose the column for every row in a loop, so the last hit's
  genic status decided the column for every hit: if it was genic, all
  intergenic hits were labelled with ``nearest_gene`` instead of their
  ``LEFT-RIGHT`` ``top_gene`` label, and vice versa.  Each hit now gets
  its own label via ``gene_label``; cached overlays that predate that
  column get it filled in at plot time.  Present in released 0.4.2.

- **Variants lifted onto a different chromosome kept their old
  ``CHR``.**  The liftover helpers took the mapped position but
  ignored the mapped chromosome, so these variants were plotted at a
  coordinate from another chromosome.  They are now treated as
  unmapped (``maps_to_other_chromosome``) and reported as above:
  about 0.05% of mapped positions in a random sample.  Present in
  released 0.4.2.

- **LD build-mismatch error named the wrong track.**  The message
  referenced a leftover loop variable; it now refers to the loaded
  tracks as a group.

- **``ax`` and ``chrom_label_size`` now work on the linear
  plotters.**  Both have been accepted and documented since 0.4.1 but
  had no effect: a loop variable overwrote ``ax``, and the chromosome
  labels were always drawn at matplotlib's default size.  Passing
  ``ax=`` to :func:`~pycmplot.plotting.linear.linear` /
  :func:`~pycmplot.plotting.linear.plot_linearm` now draws the plot
  (annotation panel plus one row per track) inside that axes' slot,
  so it can share a figure with other panels; the slot's placeholder
  axes is removed, the rest of the figure is left alone, and no file
  is saved (the caller saves the figure), matching
  :func:`~pycmplot.plotting.circular.circular`.  Axes created with
  ``fig.add_axes`` work as well as gridspec slots.
  ``chrom_label_size`` now sets the chromosome label size.  **Its
  default of 6 makes linear chromosome labels smaller than before
  (matplotlib's default, 10)**; pass ``chrom_label_size=10`` /
  ``-cl_size 10`` for the old look.  The CLI's ``-cl_size`` /
  ``--chrom_label_size`` and ``-tl_size`` / ``--track_label_size``
  now apply to linear plots as well (previously listed as
  circular-only and not passed to the linear plotter).

- **Stale highlight docstrings.**  :func:`~pycmplot.plotting.linear.linear`
  documented ``highlight_thresh`` with a default of ``5e-8`` (it is
  ``None``) and did not describe what it does, and
  :func:`~pycmplot.plotting.circular.circular` did not document it.
  Both now describe the per-locus filter above, and
  :func:`~pycmplot.plotting.linear.plot_linearm` documents its
  ``highlight_thresh`` / ``highlight_window`` parameters.

- **``pycmplot.qq_single`` was missing, and ``from pycmplot import *``
  failed.**  ``qq_single`` (one QQ plot drawn onto an ``Axes`` you
  supply, for custom and multi-panel figures) was listed in
  ``pycmplot.__all__`` but never imported, so ``pycmplot.qq_single``
  raised ``AttributeError`` (although ``plot_qq_single``'s deprecation
  warning points users to it) and any ``from pycmplot import *``
  failed.  It is now exported; the deprecated ``plot_qq_single`` is
  listed once, among the deprecated aliases.
  :mod:`pycmplot.plotting` now also exports the current names
  (``linear``, ``circular``, ``qq_single``, ``qq_combined``,
  ``qq_separate``, ``qq_overlay``) alongside the deprecated ``plot_*``
  aliases, and the API reference describes ``qq_single``.  Present in
  released 0.4.2.

- **Circular plot crashed when no chromosome's data started near its
  first base.**  The track labels in the spacer sector were placed at
  ``(end - start) / 6``, a width rather than a position, so pycirclize
  raised ``ValueError: ... is invalid range of 'Spacer1' sector``
  whenever that sector did not start at 0 (e.g. regional or heavily
  trimmed summary statistics).  Labels are now placed one sixth into
  the sector; plots whose spacer starts at 0, as with genome-wide
  data, are unchanged.  Present in released 0.4.2.

- **cCRE, CpG, DNase and TFBS tracks now used in the LD-block
  tie-break.**  The UCSC functional-annotation entry above describes
  them breaking ties between variants with identical P-values when
  choosing each lead's annotation representative (``annot_pos``), but
  they were never loaded: only GeneHancer and GTEx eQTL scored.  All six
  tracks now contribute, with the weights in
  :func:`~pycmplot.annotation.position_informativeness_score`: eQTL
  3.0 and GeneHancer 2.0 when they link to a candidate gene of the
  locus (a protein-coding gene within the annotation window of the
  lead, ``-aw``, default 500 kb; new helper
  :func:`~pycmplot.annotation.candidate_gene_pool`), otherwise 1.5 and
  1.0; cCRE 1.5, CpG 1.0, DNase 1.0, TFBS 0.5.  Overlaps use each
  file's coordinate convention (BED for the UCSC tracks).  Only tied variants are affected; the statistical lead
  (``SNP`` / ``POS``) never changes.  Each track file is now read once
  per run and shared between the tie-break, gene scoring and the
  hits-table columns (a GeneHancer-enabled load of one track went from
  about 16 s to 11 s in testing).  Affected unreleased 0.4.3 code only.

---


0.4.2 - 2026-09-12
------------------------------------------------------------------------------

**Added**

- **Signed-statistic support in the loader (iHS / XP-EHH / Fay & Wu's
  H / Tajima's D, etc.).**  When ``logp=False`` and the score column
  carries negative values, the loader now:

  * Adds a companion ``P_UNSIGNED = |P|`` column to each track's
    DataFrame.
  * Routes lead-SNP extraction and highlight-window selection through
    ``P_UNSIGNED`` (via new ``score_col`` / ``ascending`` parameters
    on :func:`~pycmplot.stats.get_lead_snps` and
    :func:`~pycmplot.stats.get_highlight_snps`), so both tails of the
    distribution contribute leads.  The original signed column stays
    in ``P`` and drives the y-axis unchanged.
  * Requires an explicit ``signif_threshold`` — the p-value fallback
    ``max(0.05/N, 5e-8)`` is meaningless on ``|value|``.  Loader raises
    a clear :class:`ValueError` pointing at the missing threshold.
  * Records a mirrored ``genome_neg = -signif_threshold`` (and
    ``suggestive_neg`` where applicable) in the per-track
    ``signif_lines`` dict.  Linear and circular plotters draw both
    ``+threshold`` and ``-threshold`` dashed reference lines when
    those keys are present, so both selection tails are annotated.
  * Clamps each reference-line y-value to the observed data range
    (``min(threshold, max(P))`` on the positive tail,
    ``max(-threshold, min(P))`` on the negative tail) so a
    threshold beyond the data extremes still draws at the plot edge
    instead of floating off-screen.

  Before this change, ``signif_threshold=4`` on signed data selected
  ``iHS <= 4`` (kept all negatives, missed positive-selection hits) and
  drew a single reference line — the semantic was inherited from the
  p-value path.  Regression: ``benchmark/tests/test_signed_stats.py``.

- **Plot-time significance filter on the hits table.**  New
  ``signif_threshold`` parameter on
  :func:`~pycmplot.plotting.linear.plot_linear` and
  :func:`~pycmplot.plotting.circular.circular`, and matching CLI flag
  ``-psig`` / ``--plot_signif_threshold``.  Loci whose lead SNP fails
  this cutoff are dropped from gene-label annotations without changing
  the loaded data, highlighted points, or reference lines.  Enables a
  "load broadly, annotate strictly" workflow: run the loader with a
  permissive ``--signif_threshold`` / ``--highlight_thresh`` to keep a
  rich hits table (useful for hand-editing the ``hits.<group>.tsv``
  overlay), then tighten annotation stringency at plot time.  Backed
  by a new :func:`pycmplot.annotation.filter_hits_by_signif` helper
  that auto-detects signed vs unsigned by inspecting the hits table's
  ``P`` column.

- **``highlight_legend`` on/off toggle on both plotters.**  A
  boolean parameter (default ``True``) that suppresses the
  "Highlighted Categories" legend entirely when set to ``False``.
  Designed for the multi-panel figure workflow where the same
  categories apply to every panel — rendering the legend on the
  first panel only lets it read as a shared legend for the whole
  figure without every panel drawing its own copy.  Exposed on
  :func:`~pycmplot.plotting.linear.plot_linear`,
  :func:`~pycmplot.plotting.linear.plot_linearm`,
  :func:`~pycmplot.plotting.circular.plot_circular`, and on the CLI
  as ``-no_hll / --no_highlight_legend`` (store_true).

- **``highlight_legend_loc`` parameter on both plotters.**  The
  "Highlighted Categories" legend can now be positioned anywhere on
  the plot without editing pycmplot's source.  Accepts any
  matplotlib ``loc`` string (``"upper center"`` — the new default —
  ``"upper right"``, ``"lower center"``, ``"best"``, etc.) or a
  2-tuple ``(x, y)`` for :func:`~matplotlib.axes.Axes.legend`
  ``bbox_to_anchor`` placement (useful for anchoring outside the
  polar frame, e.g. ``(0.5, -0.05)``).  Exposed on
  :func:`~pycmplot.plotting.linear.plot_linear`,
  :func:`~pycmplot.plotting.linear.plot_linearm`,
  :func:`~pycmplot.plotting.circular.plot_circular`, and on the CLI
  as ``-hll_loc / --highlight_legend_loc``.  The default was
  formerly ``"upper right"`` (linear) and lower-left ``bbox_to_anchor``
  (circular); both are now ``"upper center"`` for consistency and to
  avoid overlapping the annotation panel on the right side.  Shared
  translation helper:
  :func:`~pycmplot.annotation.resolve_highlight_legend_placement`.

- **Per-file ``col:<name>`` (or bare ``col`` / ``C``) token in
  ``--build`` / ``build_list=``.**
  Individual entries can now be either a literal build name
  (``hg18``/``hg19``/``hg38`` and their aliases — the original
  behaviour) *or* a column reference of the form ``col:<colname>``
  (equivalently ``C:<colname>``), meaning "for this file, source
  per-row builds from the named column".  This is the escape hatch
  for the awkward case where most files declare a single literal
  build but one file happens to carry a per-row build column with
  a non-standard name that pycmplot's auto-detector can't find.
  Example: ``--build hg19,hg38,col:my_build`` treats file 1 as
  entirely hg19, file 2 as entirely hg38, and reads file 3's builds
  row-by-row from its ``my_build`` column.  Column lookup is
  case-insensitive; a missing column errors upfront naming the
  offending file and requested column.  Combines transparently with
  the group-level mixed-build detection (see *Fixed (liftover)*
  below), so any hg18/hg19 rows surfaced via the referenced column
  are still lifted to hg38.  For the common case where the user
  knows a file has *some* build column but doesn't want to look up
  its exact name, a bare marker — ``col``, ``column``, ``C``,
  ``col:``, ``C:`` — also works: pycmplot then falls back to the
  standard ``BUILD`` / ``Genome`` / ``Genome_Build`` /
  ``Genome-build`` candidate list, and only errors out (with a
  friendly hint) if none of them appears in that file's header.

- **Batch editing across all three overlay backends.**  For overlays
  with dozens to hundreds of loci, single-cell editing was tedious.
  Three new mechanisms, all using the same pandas
  :meth:`~pandas.DataFrame.query` grammar so users learn one syntax:

  * **TUI**: ``/`` opens a filter modal (grid re-renders to matching
    rows); ``Space`` toggles row selection; ``Ctrl+A`` selects every
    filtered row; ``Ctrl+E`` opens a batch-edit modal with a
    column-toggle (``highlight_color`` / ``category``), value input,
    and live colour swatch — one value applied to every selected
    row.  Selection is keyed on the *original* dataframe index so
    it survives filtering and refresh; ``Esc`` clears both filter
    and selection in one keystroke.  When no rows are explicitly
    selected, batch-edit falls back to "the current filtered view"
    — so ``/ … Ctrl+E`` is a two-step batch flow.

  * **New** ``pycmplot hits`` **CLI**.  Scriptable batch operations
    over the same overlay files, no GUI required::

        pycmplot hits list  --cache_dir X [--where EXPR] [--columns CHR,POS] [--tsv]
        pycmplot hits set   --cache_dir X --where EXPR [--color VAL] [--category VAL] [--dry-run]

    ``--dry-run`` prints the N rows that would change and exits
    without writing.  Invalid colours are rejected pre-write; a bad
    ``--where`` prints the offending expression cleanly rather than
    a pandas traceback.  The CLI needs no extras — ideal for
    headless CI / Makefile steps that reproducibly apply a colouring
    policy.

- **Terminal-UI editor for the hits overlay TSV.**  A second
  backend for the ``pycmplot edit`` subcommand, selected via the
  ``--tui`` flag, renders a Textual-based spreadsheet directly in
  the terminal — no browser or port forwarding needed.  Ideal for
  SSH-only cluster sessions where the browser backend isn't
  usable.  Feature parity with the Streamlit backend: cell edits
  via an F2/Enter modal with a live truecolor swatch preview,
  category autocomplete-by-hint, add/delete user rows, save via
  Ctrl+S (writes through
  :func:`~pycmplot.cache.write_hits_overlay`), and a Ctrl+P
  preview that renders a matplotlib PNG to
  ``<cache_dir>/preview.png``.  Install with
  ``pip install "pycmplot[editor-tui]"`` — the two backends are
  independent extras, so users can install just the one they need
  (or both).

- **Browser-based editor for the hits overlay TSV.**  A new
  ``pycmplot edit --cache_dir <path>`` subcommand launches a local
  Streamlit app that presents the ``hits.<group>.tsv`` as a typed
  data-editor grid — ``source`` as an ``auto``/``user`` dropdown,
  ``category`` as an autocomplete-with-add selectbox, and
  ``highlight_color`` as a text field backed by a live colour-swatch
  preview strip (valid colours render at their true colour, ``auto``
  shows as a dashed grey chip, invalid values render red).  Save
  delegates to :func:`~pycmplot.cache.write_hits_overlay` so the
  atomic-write and inheritance-across-regeneration guarantees hold
  identically to the text-editor workflow; a "Preview plot" button
  re-renders a linear Manhattan against the cached tracks and the
  in-memory overlay so users see changes reflected before saving.
  Streamlit is an **optional extra** — install with
  ``pip install "pycmplot[editor]"`` — so headless / CI pipelines
  don't pay the install cost.  The ``edit`` subcommand is dispatched
  via a preflight in :func:`~pycmplot._core.main` (rather than an
  argparse subparser) so existing ``pycmplot --sum_stats …``
  invocations remain unchanged.

**BREAKING (public API) — soft-deprecation, not removal**

- **Seven core public-API functions renamed to short, memorable
  forms.**  Before pycmplot leaves the pre-release window the old
  descriptive-but-verbose names have been retired in favour of
  concise verbs that read naturally with a module prefix
  (``import pycmplot as pcm``):

  ================================================================  ============================
  Old name (deprecated, still importable)                           New short name
  ================================================================  ============================
  ``prep_pycmplot_input_info``                                       ``prep``
  ``get_sumstats_and_merged_sector_list``                            ``load``
  ``plot_linear``                                                    ``linear``
  ``plot_circular``                                                  ``circular``
  ``plot_qq_combined``                                               ``qq_combined``
  ``plot_qq_overlay``                                                ``qq_overlay``
  ``plot_qq_separate``                                               ``qq_separate``
  ================================================================  ============================

  Each old name remains importable for **one release cycle** and
  its first call per Python process emits a single
  ``DeprecationWarning`` naming the new short form.  A shared shim
  in :mod:`pycmplot._deprecation` preserves ``__doc__``,
  ``__wrapped__`` and ``__signature__`` so IDE hints and Sphinx
  autofunction resolve to the real function.  Every internal
  caller (``_core``, ``editor``, ``editor_tui``, ``hits_cli``,
  ``benchmark/bench_python.py``) has been updated to the new
  names in the same commit; only user scripts pinning the old
  names will see the warning.

  New idiomatic pattern in the tutorial / README::

      import pycmplot as pcm
      fi     = pcm.prep(sum_stats=files, labels=labels)
      bundle = pcm.load(sum_stats=files, labels=labels, file_info=fi,
                        logp=True, trim_pval=0.01, highlight=True)
      pcm.linear(sumstats_loaded=bundle["dfs"],
                 hits_table=bundle["annot"],
                 signif_lines=bundle["lines"], output_dir="./out")
      pcm.circular(sumstats_loaded=bundle["dfs"],
                   sector_sizes=bundle["sectors"],
                   hits_table=bundle["annot"])

  Internal engines (``plot_linearm``, ``plot_circosm``) keep their
  ``m`` suffix — they signal "give me an axes and I'll render" and
  aren't the recommended entry point.

**Changed (annotation)**

- **``BIOTYPE_WEIGHTS`` reweighted against Ensembl 116 hierarchy.**
  The prioritisation table in :data:`pycmplot.constants.BIOTYPE_WEIGHTS`
  has been aligned with the current Ensembl biotype classification
  (`Ensembl release 116, June 2026
  <https://jun2026.archive.ensembl.org/info/genome/genebuild/biotypes.html>`_).
  Two visible changes and a large silent one:

  * **``antisense`` raised from 0.30 → 0.65.**  Ensembl 116 classifies
    ``antisense`` as an lncRNA subtype (transcripts that overlap the
    genomic span of a protein-coding locus on the opposite strand),
    not a pseudogene-tier biotype.  The old weight ordered it below
    ``processed_pseudogene`` (0.30), which was inconsistent with the
    sibling lncRNA subtypes (``lncRNA`` / ``lincRNA`` at 0.70).
  * **33 new biotype entries added.**  The pre-0.4.x table had 20
    entries; the new one has 64.  Previously-missing categories that
    used to silently default to weight 0 now score correctly:

    - **Immune-receptor genes** — ``IG_C_gene`` / ``IG_D_gene`` /
      ``IG_J_gene`` / ``IG_V_gene`` / ``TR_C_gene`` / ``TR_D_gene`` /
      ``TR_J_gene`` / ``TR_V_gene`` at 1.00 (protein-coding
      equivalents).
    - ``protein_coding_LoF`` (0.90) and
      ``protein_coding_CDS_not_defined`` (0.85) — still coding for
      some individuals or isoforms.
    - ``nonsense_mediated_decay`` / ``non_stop_decay`` at 0.55.
    - Small ncRNAs previously missing: ``piRNA`` (0.70), ``siRNA``
      (0.70), ``tRNA`` (0.60), ``Mt_tRNA`` (0.60), ``Mt_rRNA``
      (0.55), ``miscRNA`` (0.55).
    - lncRNA subtypes: ``macro_lncRNA``, ``non_coding``,
      ``3prime_overlapping_ncRNA``, ``sense_intronic``,
      ``sense_overlapping``, ``retained_intron``.
    - Pseudogene subtypes: ``polymorphic_pseudogene``,
      ``translated_processed_pseudogene``,
      ``translated_unprocessed_pseudogene``, ``unitary_pseudogene``,
      and the eight ``IG_*_pseudogene`` / ``TR_*_pseudogene``
      subtypes.
    - Speculative / provisional: ``TEC`` (0.30), ``readthrough``
      (0.35), ``stop_codon_readthrough`` (0.55), ``artifact`` (0.15).

  * **Legacy-spelling aliases** so ``vault_RNA`` / ``vaultRNA``,
    ``misc_RNA`` / ``miscRNA``, ``long_intergenic_ncRNA`` /
    ``lincRNA``, and ``3prime_overlapping_ncRNA`` /
    ``3_prime_overlapping_ncRNA`` all resolve to the same weight —
    no more silent zero-scores when the reference file uses a
    non-canonical spelling.

  Only the ``priority_score`` component of the intergenic tier is
  affected; positional flanker selection and the ``nearest_gene``
  columns are unchanged (they don't use the weight table).  The
  annotation-flowchart, schematic, and priority-score panel scripts
  reflect the new table.

**Changed (annotation, internal)**

- **Two annotation passes merged into a single window walk.**
  Before this release, ``_annotate_variant`` (positional nearest-gene
  detection) and ``_annotate_and_prioritize_variant`` (biotype-weighted
  priority scoring) walked every candidate gene in the ±500 kb window
  independently, then had their output dicts merged inside
  ``get_hits_summary_table``.  Since both computed the same distances,
  the same genic checks, and the same gene_density, the redundancy was
  wasted work — and any drift between the two passes would silently
  produce inconsistent columns in the same row.  Both passes are now
  consolidated inside a single ``_annotate_variant`` that emits the
  union of the previous two dicts in one traversal, roughly halving the
  annotation-step cost and giving one source of truth for the shared
  semantics.  ``_build_genes_dict`` now includes ``BIOTYPE`` in the
  interval tuples (5-tuple ``(start, end, strand, gene, biotype)``) so
  the priority scoring never needs to re-touch the source DataFrame.
  Column names + values on the hits summary table are unchanged; every
  pre-merge test scenario still passes when routed through the merged
  function.

**Fixed**

- **``signif_threshold`` now propagates into the drawn reference
  line.**  The loader was re-initialising ``resolved_signif_line``
  from ``max(0.05/N, 5e-8)`` on every iteration regardless of the
  supplied ``signif_threshold``.  For unsigned data the two values
  coincidentally agreed; for signed data
  ``signif_threshold=4`` was set correctly for lead-picking but the
  reference line silently rendered at ``5e-8`` (a p-value scale) instead
  of at ``4``.  ``resolved_signif_line`` now initialises from
  ``signif_threshold`` and is only overridden by an explicit
  ``signif_line=<float>``.  Bonus: on unsigned data, a hand-picked
  ``signif_threshold`` (e.g. ``1e-6``) that used to disagree with the
  drawn line now matches by default; pass ``signif_line=<value>``
  to keep them different.

- **``prep()`` column resolution: additions to ``pvl_candidates``
  are now honored.**  The per-file inner loop was rebuilding
  ``pvl_cands`` from a hard-coded literal list, silently discarding
  any additions made to ``pvl_candidates`` at function scope (e.g.
  ``IHS``, ``RSB``, ``LOGP``).  The rebuild is removed; the outer
  list is now the single source of truth.

- **``prep()`` user-supplied column hints now take priority.**
  Two-pass resolution: (1) case-insensitive match against the file
  header for ``chrom`` / ``pos`` / ``snp`` / ``pcol`` if the user
  supplied one; (2) leftmost header column matching any built-in
  candidate.  Fixes the bug where a header like
  ``SNP CHR POSITION IHS LOGPVALUE P BH_adj_P Bonf`` resolved the
  p-value column to ``IHS`` (leftmost header match against the
  candidate *set*) even when the user explicitly passed ``pcol="P"``.
  A hint that names a column not present in the header now errors
  out clearly rather than silently falling back to a different
  column.

**Fixed (safety)**

- **``--clear_cache`` no longer removes the parent ``--cache_dir``
  itself.**  A serious footgun: earlier releases did
  ``shutil.rmtree(cache_dir)``, which meant that pointing
  ``--cache_dir`` at a directory that also held analysis notes,
  plot outputs, or any other unrelated files would silently blow
  everything away.  The clear routine now targets only the three
  artefacts pycmplot writes — ``metadata.json``, ``tracks/``, and
  ``annotations/`` — and inside those subdirectories only deletes
  files whose names match the pycmplot filename conventions
  (parquet + leads + pvals sidecars under ``tracks/``, and
  ``hits.<group>.tsv`` / ``.meta.json`` under ``annotations/``).
  Foreign files are kept and reported at INFO level.  The parent
  ``--cache_dir`` is never removed.  Before deleting
  ``metadata.json`` the routine also verifies it looks like a
  pycmplot manifest (has ``cache_version`` or ``tracks`` keys) —
  a stray ``metadata.json`` from another tool is left in place.

**Fixed (liftover)**

- **Cross-file build detection.**  Before this release, the liftover
  trigger only inspected each file's own ``BUILD`` column and fired
  liftover when that single file mixed hg19 and hg38 rows.  A loader
  group like ``build_list=['hg19', 'hg38']`` — where each file was
  internally single-build but the group as a whole was mixed —
  silently left every track in its native coordinate system, so
  cross-track highlight lines drew at different genomic positions
  (bug reported 2026-09-05).  The loader now runs a group-level
  build scan **before** the per-track loop: it collects each track's
  declared build from ``build_list=`` (or ``--build``) and, when
  absent, from a light-touch ``pd.read_csv(nrows=2000)`` peek at the
  file's ``BUILD`` column.  When more than one distinct build is
  detected across the group, every hg18/hg19 track is lifted over to
  hg38 so that shared coordinates and highlight lines line up.  The
  logger reports which builds were detected and why liftover fired
  (``same-file mixed builds`` vs.
  ``cross-file mixed builds in the loader group``).  Single-build
  groups (``[hg19, hg19]``, ``[hg38]``, etc.) are unchanged.

**Fixed (annotation)**

- **``nearest_upstream_gene`` / ``nearest_downstream_gene`` now use
  positional semantics.**  The previous release defined "upstream"
  strand-aware, meaning "the variant is 5' of the gene's TSS
  (respecting strand)".  Each candidate gene got assigned to exactly
  one of the two buckets based on its strand, so a very close gene
  sitting to the left of the variant but whose 3' end faced the
  variant landed in the ``downstream`` slot, and the reported
  ``nearest_upstream_gene`` was often not the positionally nearest
  gene at all.  From this release, "upstream" means **lower genomic
  coordinate** (left of the variant) and "downstream" means **higher
  coordinate** (right of the variant), matching intuition and the
  convention used by most external tools.  The strand-aware
  definition is retained only for :data:`promoter_upstream_flag`
  where it's genuinely biological.

- **Intergenic** ``top_gene`` **now actually flanks the variant.**
  Previously the intergenic branch of the prioritiser did
  ``candidates.head(2)`` on a priority-sorted list, which could
  return two genes both sitting on the same side of the variant
  (two adjacent coding genes to the left would beat a farther coding
  gene on the right).  Now one gene is picked from each positional
  side by base-pair distance to the gene body; the joined
  ``LEFT-RIGHT`` label always references genes that actually
  flank the variant.  When only one side has a gene in the search
  window, a single symbol is returned rather than a joined pair.

- **New fields** ``nearest_gene`` **and** ``nearest_gene_distance``
  in the hits table — the actual closest gene by base-pair distance
  regardless of side.  This is what users typically mean by
  "nearest".  ``get_annotation_column("GENE")`` now prefers
  ``nearest_gene`` over ``nearest_upstream_gene`` for genic labels
  (with a fallback to ``nearest_upstream_gene`` when reading legacy
  cached hits tables that predate the new column).

**Docs**

- Tutorial note under :ref:`cli-tut-linear` explaining signed-stat
  behavior for iHS / XP-EHH and how ``signif_threshold`` /
  ``highlight_thresh`` are applied on ``|value|``.

---


0.4.1 - 2026-08-23
------------------------------------------------------------------------------

**Added**

- **``track_label_size`` and ``chrom_label_size`` on the linear
  plotters.**  :func:`~pycmplot.plotting.linear.plot_linear` and
  :func:`~pycmplot.plotting.linear.plot_linearm` gain both parameters
  (default ``6``).  ``track_label_size`` sets the track-label font
  size, which was previously fixed at ``10``, so track labels are
  smaller by default.  ``chrom_label_size`` is accepted but not yet
  applied.

- **``ax`` parameter on the linear plotters** for multi-panel
  figures.  Accepted and documented, but not yet applied: the linear
  plotter still creates its own figure.

**Changed**

- **Linear plot guide lines** (``highlight_line=True``) are now drawn
  once, at each position in the hits table, above the data points
  (``zorder=4``).  Previously a line was drawn at every highlighted
  variant, separately for each track, beneath the data (``zorder=0``).

- **Circular highlight legend moved to the upper right**
  (``loc="upper right"``, ``bbox_to_anchor=(1.05, 1.1)``); previously
  upper left (``bbox_to_anchor=(-0.05, 0.0)``).

**Removed**

- **Unused ``highlight_thresh`` parameter** of
  :func:`~pycmplot.plotting.circular.circular`.  It was accepted but
  had no effect; the CLI no longer passes it.

**Fixed**

- **Hits overlay writes missing values as ``NA``** instead of empty
  fields, so hand-edited ``hits.<group>.tsv`` files round-trip
  cleanly.


---


0.4.0 - 2026-08-13
------------------------------------------------------------------------------

**Changed**

- **QQ plot legend no longer includes a "95% CI" entry.**  The CI
  band is visually self-evident (it's the shaded region hugging the
  diagonal), and the extra legend entry crowded the top corner
  alongside the track label and the ``λ`` annotation.  Legend now
  shows only the track name(s); band still renders identically.
  Callers who prefer the old behaviour can add an
  overlay/annotation of their own.

- **Pvals sidecar switched from uncompressed ``.npy`` to compressed
  ``.npz`` with sorted ``float32`` payload.**  The QQ plotter sorts
  the array internally, so pre-sorting is semantically transparent;
  ``float32`` preserves ``median(-log10 P)`` to ~9 decimal digits
  (well beyond any GWAS precision concern).  Result: **~5× smaller
  on-disk footprint** — 8 MB → 1.7 MB on a 1 M-variant sidecar,
  scaling proportionally to 80 MB → ~17 MB at 10 M variants.
  Parquet was evaluated for consistency with the other cache files
  but performed *worse* (117–71% of raw) due to per-column overhead
  on random-ish float payloads; ``np.savez_compressed`` reaches
  ~20–30% of raw.  Read path is backwards-compatible: legacy
  ``.npy`` and short-lived ``.parquet`` sidecars are still readable
  transparently.

- **Hits overlay filename is now group-scoped.**  The file layout
  changes from a single ``<cache_dir>/annotations/hits.tsv`` to
  ``<cache_dir>/annotations/hits.<group_key>.tsv``, where
  ``group_key`` is a short SHA-256 of the sorted list of per-track
  cache_keys in the loader call.  Because each cache_key already
  encodes the raw file SHA-256 plus every Stage-1 parameter, each
  distinct ``(files, parameters)`` combination gets its own hits
  overlay.  This mirrors the per-track cache invalidation model and
  matches the reproducibility mental model: different parameters →
  different analysis → different cached artefacts.  Multi-panel
  workflows are safe (two panels with different sumstats never
  clobber each other's overlays), and warm re-runs with the same
  ``(files, parameters)`` reuse and update the same file — so
  user-authored rows persist within one parameter setting.  Changing
  a Stage-1 parameter (e.g. ``highlight_thresh`` or ``trim_pval``)
  spawns a fresh overlay under a new ``group_key``; the earlier
  overlay lingers on disk unchanged (users who want to carry a
  hand-edit forward can copy rows manually).

**BREAKING (Python API)**

- **Default for** ``compute_pvals`` **flipped from** ``True`` **to**
  ``False`` on :func:`~pycmplot.io.get_sumstats_and_merged_sector_list`.
  The bundled p-value array is expensive to materialise (~80 MB per
  track at 10 M variants) and is only needed for QQ plots, so opting
  in matches what the CLI already does (``-qq/--qq_plot`` sets it
  automatically).

  **Migration** — Python-API scripts that call the loader with
  defaults and then feed ``bundle['pvals']`` into a QQ plotter must
  now pass ``compute_pvals=True`` explicitly:

  .. code-block:: python

     # Before 0.4.0 — bundle['pvals'] populated implicitly
     bundle = get_sumstats_and_merged_sector_list(..., logp=True)
     plot_qq_combined(bundle['pvals'], output_path='qq.png')

     # 0.4.0 onward — opt in
     bundle = get_sumstats_and_merged_sector_list(
         ..., logp=True, compute_pvals=True,
     )
     plot_qq_combined(bundle['pvals'], output_path='qq.png')

  The four public QQ plotters (:func:`~pycmplot.plotting.qq.plot_qq_single`,
  :func:`~pycmplot.plotting.qq.plot_qq_combined`,
  :func:`~pycmplot.plotting.qq.plot_qq_separate`,
  :func:`~pycmplot.plotting.qq.plot_qq_overlay`) now raise
  ``ValueError`` with a directive message when handed ``None`` p-values,
  so callers that missed this migration get a clear signpost rather
  than a mystery downstream ``TypeError``.  Manhattan / circular
  workflows are unaffected.

**Added**

- **Per-locus categories driving a custom highlight legend.**  The
  auto-generated hits table now carries a ``category`` column,
  populated with the sentinel ``"significant"`` at cache-write time.
  Users hand-edit rows to give each locus a label
  (e.g. ``"novel"``, ``"replicated"``, ``"MHC"``); **both the linear
  and circular plotters** then render a "Highlighted Categories"
  legend, with one entry per unique category in first-appearance
  order — so re-ordering rows in the TSV re-orders the legend.  The
  linear plotter anchors the legend in the top-right of the topmost
  data axes; the circular plotter anchors it below the polar sectors
  (which have no natural interior real-estate).  When everything is
  left at the defaults (all rows ``category=significant`` and
  ``highlight_color=auto``), **no legend is added** — the pre-feature
  layout is preserved.  As with the ``highlight_color`` column,
  user-set categories are inherited across cache regenerations by
  ``(CHR, POS)`` lookup (no need to also change ``source='auto'`` to
  ``user``).  Both plotters delegate the entry-building logic to the
  shared :func:`~pycmplot.annotation.build_highlight_legend_entries`
  helper (deduplicates on category; warns when the same category
  appears with multiple colours; tolerates legacy caches without the
  ``category`` / ``highlight_color`` columns).  Companion helper:
  :func:`~pycmplot.annotation.resolve_highlight_categories`.

- **Per-locus highlight colors via the hits overlay.**  The
  auto-generated hits table now carries a ``highlight_color`` column,
  populated with the sentinel ``"auto"`` at cache-write time.  Users
  can hand-edit any row's value to any matplotlib-parseable color
  (name, hex ``#rrggbb``, RGB tuple) to give that locus its own
  highlight colour; ``"auto"`` / blank / NaN / invalid entries fall
  back to the plot-time ``highlight_color`` argument.  Applied by
  both :func:`~pycmplot.plotting.linear.plot_linear` and
  :func:`~pycmplot.plotting.circular.plot_circular` — each
  ``in_locus`` variant is looked up against its nearest lead in the
  overlay (within 500 kb) and coloured accordingly.  Invalid colours
  emit a ``logger.warning`` naming the offending value so typos in
  the TSV are easy to fix.

  User-set colours **survive regeneration** even when they were made
  to auto rows (without changing ``source`` to ``user``): the writer
  looks up each new auto row's ``(CHR, POS)`` in the previous overlay
  and inherits any non-``"auto"`` colour it finds.  Users don't need
  to know about the ``source`` column just to recolour a locus.

  A pre-existing bug in the overlay reader — pandas'
  ``comment='#'`` truncating hex-colour values mid-cell — was fixed
  as part of this work.  Header comments are now stripped only when
  they appear at the top of the file, so hex codes are read
  faithfully.

- **User-editable hits overlay** (``<cache_dir>/annotations/hits.tsv``).
  When ``--cache`` is enabled, the auto-generated hits table is
  persisted as a plain TSV with a ``source`` column
  (``"auto" | "user"``).  Regeneration (triggered when the leads or
  the annotation resource files change) only replaces rows tagged
  ``"auto"``; rows tagged ``"user"`` survive every invalidation.
  Users can therefore edit the TSV directly to add custom loci
  (e.g. an *MHC region* label; a meta-analysis lead absent from any
  input file), correct auto-picked gene names, or suppress
  false-positive leads — all without touching Python.  A header
  comment in the TSV documents the schema.  Skipping the annotation
  pass on a cache hit saves an additional 1–5 s per run on top of the
  per-track cache.

- **Per-track Stage-1 cache and resume** (``pycmplot.cache``, new
  ``--cache`` / ``--cache_dir`` / ``--no_resume`` / ``--clear_cache``
  CLI flags; ``cache=``, ``cache_dir=``, ``resume=`` kwargs on
  :func:`~pycmplot.io.get_sumstats_and_merged_sector_list`).  When
  ``cache=True``, each input file's post-load / post-liftover /
  post-thinning DataFrame is written to
  ``<cache_dir>/tracks/<label>.<short_key>.parquet`` alongside its
  lead-SNP table and (optionally) its raw p-value array.  A JSON
  metadata file records a per-track ``cache_key`` computed from the
  raw file's SHA-256, the pycmplot version, and every Stage-1
  parameter that affects the cached data (``trim_pval``,
  ``auto_thin*``, ``highlight*``, ``signif_threshold``, ``build``,
  ``logp``).  Subsequent runs skip Stage 1 entirely for tracks whose
  cache_key matches — a 1.9×–10× speedup that grows with input size.
  A parameter change or file rewrite invalidates the affected entries
  silently, so users don't need to remember to clear the cache; the
  ``--clear_cache`` flag is provided for explicit resets.  Multi-panel
  callers that reuse a label (e.g. ``"Hb"`` in two different panels
  pointing at different sumstats) get two coexisting cache entries,
  disambiguated by the ``<short_key>`` suffix in the filename.  The
  "resume" semantics (as described in ``to-do.md``) fall out for
  free: every completed track is atomically committed to
  ``metadata.json`` before the next iteration begins, so a run that
  crashes on track N leaves tracks 1…N−1 cached and the next
  ``pycmplot --cache`` invocation continues from track N.

- **``--version`` / ``-V`` flag** on the CLI — prints the installed
  pycmplot version and exits.

- **Directive-error guards in the QQ plotters.**  When any of
  :func:`~pycmplot.plotting.qq.plot_qq_single`,
  :func:`~pycmplot.plotting.qq.plot_qq_combined`,
  :func:`~pycmplot.plotting.qq.plot_qq_separate`, or
  :func:`~pycmplot.plotting.qq.plot_qq_overlay` receives ``None``
  p-values (either a bare ``None`` or a ``pval_dict`` containing
  ``None`` values), a clear ``ValueError`` is raised naming the
  offending track and pointing the user at the ``compute_pvals=True``
  loader kwarg.  Complements the compute_pvals default flip above
  so migration hiccups produce a directive message rather than a
  downstream ``TypeError``.

- **Multi-Panel Canvas Support (`plot_circular`)**:
  Added an optional `ax` parameter to `plot_circular()`, enabling users to render
  circular Manhattan plots onto existing Matplotlib polar axes (`projection='polar'`).
  Allows embedding `pyCirclize` figures into complex, multi-panel layouts using
  `matplotlib.figure.SubFigure`, `GridSpec`, or standard subplots.
  Bypasses internal auto-saving (`fig.savefig()`) when `ax` is provided, delegating
  layout control and rendering pipeline management to the user.

**Fixed**

- QQ plotters
  (:func:`~pycmplot.plotting.qq.plot_qq_single` +
  :func:`~pycmplot.plotting.qq.plot_qq_overlay` and the wrappers that
  forward to them) no longer raise
  ``TypeError: float() argument must be a real number, not a 'list'``
  when a scalar keyword argument (``fontsize``, ``point_size``) arrives
  as a single-element list, 0-d ndarray, or 1-element ``pandas.Series``.
  A new :func:`~pycmplot.plotting.qq._to_scalar_float` helper coerces
  these to plain floats at the top of each entry point; multi-element
  inputs still raise, but now with a directive message naming the
  offending kwarg (rather than surfacing from deep inside
  ``matplotlib``'s numeric parsing at an unhelpful call site).
  :func:`~pycmplot.plotting.qq._compute_lambda` was hardened at the
  same time to always return a Python ``float`` (never a NumPy scalar
  or 0-d array), so ``f"λ = {lam:.4f}"`` formatting is safe regardless
  of upstream pvals shape.

- Guarded ``df["BUILD"].unique()`` in the loader with an explicit
  ``if "BUILD" in df.columns`` check.  The previous version dereferenced
  ``df["BUILD"]`` before the guard and raised ``KeyError: 'BUILD'`` when
  the input file had no BUILD column.

- ``_atomic_write`` no longer produces ``FileNotFoundError`` when the
  writer auto-appends its own extension (notably ``numpy.save``, which
  silently appends ``.npy``).  The temp filename is now
  ``<stem>.__tmp__<ext>`` (e.g. ``Hb.pvals.__tmp__.npy``) so any
  auto-append lands on the exact path the subsequent ``os.replace``
  expects.

- ``compute_pvals`` no longer participates in the per-track cache key,
  so toggling QQ requirements between runs doesn't invalidate the
  main cache.  When a subsequent run wants pvals but no sidecar
  exists from earlier calls, Stage 1 re-runs *just for that track* to
  populate the sidecar; every run after that hits the cache with
  pvals included.  Fixes the awkward interaction between the
  ``compute_pvals=True`` Python-API default and the ``--cache`` flag.

- **Multi-panel cache safety.**  Two loader calls sharing the same
  ``cache_dir`` and reusing a track label (e.g. ``"Hb"`` in panel A
  and ``"Hb"`` in panel B, pointing at different sumstats) no longer
  clobber each other's cached parquet / metadata.  Cache metadata is
  now keyed on the full 64-char cache_key rather than the label alone,
  and per-track filenames include an 8-char slice of the key
  (``Hb.a1b2c3d4.parquet``) so multiple entries for the same label
  coexist safely.  Cache-layout version bumped to ``2``; older
  ``.pycmplot_cache/`` directories are transparently wiped on first
  access.  Documentation added for ``cache``, ``cache_dir``,
  ``resume``, ``compute_pvals``, ``auto_thin*`` in the loader
  docstring.

- Track iteration in the loader now follows the user-supplied
  ``labels`` list (order-preserving) rather than the non-deterministic
  set intersection ``sumstats.keys() & file_info.keys()``.  Without
  this, ``signif_threshold``'s auto-computation from the
  first-processed track's SNP count fed the second track's cache_key
  differently on cold vs warm runs, breaking the cache for every
  track after the first.  The cache-HIT path also now mirrors the
  cold-path side-effect on ``signif_threshold`` so downstream tracks
  see identical state.

**Documentation**

- Added NumPy-style parameter blocks for ``compute_pvals``,
  ``auto_thin``, ``auto_thin_threshold``, ``auto_thin_max_below``,
  ``cache``, ``cache_dir``, and ``resume`` on
  :func:`~pycmplot.io.get_sumstats_and_merged_sector_list`.  The
  ``resume`` entry is honestly documented as reserved-for-future — the
  actual resume semantics are inherent in the atomic per-track commit
  and do not require the flag.

- Rewrote the :mod:`pycmplot.io` module docstring to cover the
  density-aware sub-sampling algorithm, the per-track cache layout
  (with pointers to :mod:`pycmplot.cache`), the multi-panel
  content-vs-label keying, and a public-function summary listing
  ``auto_thin_for_manhattan`` and ``get_output_paths`` alongside the
  two headliners.

- New :mod:`pycmplot.cache` module docstring covering the on-disk
  layout (``metadata.json``, ``tracks/``, ``annotations/hits.tsv``),
  the ``cache_key`` construction, atomic-write semantics, and the
  user-editable hits overlay design.



----


0.3.1 - 2026-07-31
------------------------------------------------------------------------------

**Fixed**

- **LiftOver hg18/hg19 to hg38**

When hg38 was absent, hg18-only and hg19-only summary statistics were converted to
hg38. This behaviour has now been changed such that only hg18-only files are converted 
while hg19-only files remain in hg19 coordinates and the bundled ENSEMBL GFF3 file in 
GRCh37 is used for annotation.

- **Significance thresholds**

The genome-wide significant and highlighting thresholds defaults were set to ``5e-08``.
Suggestive line was only included when significance line was enabled. Significant and 
highlighting thresholds have now been explicitly set to be caulculated from number of 
markers using ``0.05 / number of markers``. A marker/snp count dictionary is now generated 
upfront before auto-thinning and pvalue trimming are applied to use in this calculation.
This effectively enables sumstat-specific thresholds and lines.
Suggestive line inclusion now depends only on whether ``-sug/--suggest_threshold`` is 
specified.

----


0.3.0 - 2026-06-01
------------------------------------------------------------------------------

**Fixed**

- Bug fix:
  - linear plotting ``t_heights`` local variable access failure.


**Changed**

- Suggestive line color from ``lightblue`` to ``navy`` in circular plotting
- Significance line color from ``red`` to ``orangered`` in linear plotting to match
  circular plotting
- Suggestive line color from ``blue`` to ``navy`` in linear plotting  to match 
  circular plotting


**Added**

- ``ylabel``: optional ylabel text in circular plotting to match linear plotting
 

----


0.2.8 - 2026-05-30
------------------------------------------------------------------------------

**Added**

- **Dual annotation renderer architecture**

Two complementary annotation functions now handle sparse and dense
annotation scenarios independently:

- ``_draw_annotation_arrows`` — sparse annotation renderer with tiered
  label placement, chromosome-boundary spreading, cumulative-distance
  stacking, and straight ``arc3`` arrows (curvature fixed at zero for
  visual clarity in low-density contexts).

- ``_draw_annotation_arrows_multirail`` — dense annotation renderer
  implementing a three-step layout pipeline (see below) with curved
  ``arc`` arrows and adaptive ``ylim``.

- **Three-step dense annotation layout pipeline** (``_draw_annotation_arrows_multirail``)

1. *Relaxation pass* — bidirectional ``min_sep`` enforcement starting
   from ``x_signal`` positions.  Labels in dense regions drift further
   from their signals than labels in sparse regions, producing a
   natural density signal with no explicit cluster detection.
2. *Drift-based rail assignment* — each label's relaxation drift is
   binned into a rail index using
   ``rail_stride = rail_width / max_rails``.  Denser regions
   automatically receive higher rail indices proportionally across the
   full rail range.  No per-rail queue processing or ``max_drift``
   threshold is required.
3. *linspace rank-reassignment* — labels are sorted by ``x_signal``
   and assigned evenly-spaced ``x_text`` slots via
   ``np.linspace(rail_start, rail_end, n)``.  This guarantees
   ``x_text`` rank equals ``x_signal`` rank (no arrow crossings by
   construction) and full rail coverage regardless of ``rail_frac`` or
   signal distribution.

- **Auto char_width from axes geometry**

For vertical text (``rotation=90``), the horizontal label footprint is
one character wide regardless of string length.  ``char_width`` is now
derived from the axes pixel extent and figure DPI at draw time::

    px_per_bp  = ax_bbox.width / (xmax - xmin)
    char_width = 0.6 * fsize * (fig.dpi / 72.0) / px_per_bp

The ``char_width_factor`` parameter has been removed from
``_draw_annotation_arrows_multirail``; ``char_width`` is computed
automatically and scales correctly with figure size, DPI, and font
size.

- **Proportional space budgeting and rail_frac awareness**

Rail width is derived from ``rail_frac`` as
``rail_width = genome_width * rail_frac``, centred on the genome
midpoint.  ``rail_stride`` and slot spacing scale proportionally with
``rail_frac``, ensuring even label distribution at any rail fraction
without choking at rail boundaries.

- **Layout table** (``pd.DataFrame``)

Placement, relaxation, and rendering are now cleanly separated via a
layout table with columns ``label``, ``x_signal``, ``x_text``,
``rail_id``.  ``rail_id`` is written during placement and not read
again until the rendering pass, enforcing strict separation of layout
and rendering concerns.

- **Chromosome-boundary detection** (``_draw_annotation_arrows``)

For each adjacent chromosome pair, the inter-chromosome gap is
computed.  If the gap is narrower than ``spread_width``, both boundary
annotations receive an ``x_bound`` value encoding direction and
magnitude, used downstream to push boundary labels apart before
general spreading.

- **Cumulative x-position porting from tracks**

Annotation cumulative x positions are now ported directly from track
DataFrames via a three-column merge on ``(chr_col, pos_col, LABEL)``
rather than being recomputed independently, guaranteeing exact
consistency between annotation and track coordinates.

- **track_heights sanity check and y-label positioning**

``track_heights`` is validated against the expected count
(``n_tracks + 1`` when annotating, ``n_tracks`` otherwise) with
explicit ``ValueError`` and ``TypeError`` messages.  The y-label
position (``-log10(P)``) is computed from actual height ratios
accounting for top-to-bottom track orientation::

    y_lab_pos = data_total / (2 * total_height)


**Changed**

- ``_draw_annotation_arrows``: ``max_rad`` parameter removed; curvature
  is intentionally fixed at zero (straight arrows) for sparse
  annotation contexts.  Dense annotation curvature is handled
  exclusively by ``_draw_annotation_arrows_multirail``.

- Annotation deduplication now occurs at the top of both renderers
  via ``drop_duplicates(subset=[chr_col, "x", label_col])`` to prevent
  replicated arrows when ``annot_df`` is a merged multi-track table.

- Chromosome order in boundary detection now uses ``natsorted`` instead
  of ``set`` to guarantee correct genomic ordering.

- ``x_bound`` is only set when the inter-chromosome gap is
  ``<= spread_width`` (previously unconditional), preventing spurious
  boundary constraints between well-separated chromosomes.


**Fixed**


- Arrow crossings eliminated unconditionally by the linspace
  rank-reassignment step: ``x_text`` rank is guaranteed equal to
  ``x_signal`` rank for all labels across all rails.

- Annotation spill past genome right boundary resolved: ``rail_end``
  acts as a hard clamp during relaxation; labels cannot exceed it
  regardless of local density.

- Higher-rail priority inversion fixed: the drift-based rail assignment
  correctly places the densest labels (largest drift) on higher rails,
  not the labels nearest the rail boundary.

- ``x_texts`` sort-order mismatch resolved: cumulative-scaled positions
  are now mapped back to original signal order via ``np.argsort``
  before use, preventing label-to-wrong-position assignment.

- ``char_width`` underestimation fixed: replacing the hardcoded
  ``8e6`` fallback with axes-geometry derivation corrects the ~2×
  underestimate that caused stacking to never fire for typical figure
  sizes at ``fontsize=6``.

- ``natsorted`` applied to chromosome order throughout to prevent
  incorrect pairing of chromosomes (e.g. chr3 with chr17) caused by
  ``set`` iteration order.


0.2.7 - 2026-04-27
------------------

**Added**

- **Default-on density-aware auto-thinning** for Manhattan / circular
  rendering, inspired by ``gwaslab`` and applied on top of (i.e. in
  addition to) the existing ``--trim_pval``.  A new helper
  :func:`~pycmplot.io.auto_thin_for_manhattan` keeps **every** variant
  whose "interestingness" signal is at or above ``--auto_thin_threshold``
  and uniformly sub-samples the dense bulk to at most
  ``--auto_thin_max_below`` rows per track (default ``200 000``).  Lead
  SNPs are still extracted from the *full* unthinned data, so peak
  annotations are unaffected.

  Two modes, switched by ``--logp``:

  * **P-value mode** (``--logp`` set, the GWAS default).  Signal is
    ``-log10(P)``; ``--auto_thin_threshold`` is in ``-log10(P)`` units
    (default ``2.0`` => ``P <= 0.01``).  Every suggestive /
    genome-wide-significant variant survives untouched.
  * **Raw-statistic mode** (``--logp`` off).  The ``P`` column is
    interpreted as a raw test statistic and the signal becomes
    ``|value|``, so the same machinery works for selection scans like
    **iHS, XP-EHH, F_ST, Fay & Wu's H, Tajima's D**, etc.  The default
    threshold of ``2.0`` works for the standardised \|iHS\| / \|XP-EHH\|
    scans; override (e.g. ``--auto_thin_threshold 0.05``) for F_ST.

  Negative extremes are preserved as well as positive ones, so for
  signed statistics (iHS, XP-EHH) both tails of the distribution
  survive intact.

  New CLI flags:

  ============================== ================================================
  Flag                           Description
  ============================== ================================================
  ``--no_auto_thin``             Disable auto-thinning entirely.
  ``--auto_thin_threshold``      ``-log10(P)`` floor above which every variant
                                 is kept (default 2.0).
  ``--auto_thin_max_below``      Cap on background variants per track
                                 (default 200 000).
  ``--no_qq_thin``               Counterpart for QQ log-uniform thinning,
                                 which is now ON by default.
  ============================== ================================================

  Combined with the rendering and data-prep optimisations from earlier
  in this release, this brings pycmplot's untrimmed timings to:

  +-------+-------------------+--------------+----------------+
  | Size  | manhattan (s)     | qq (s)       | circular (s)   |
  +=======+===================+==============+================+
  | 500K  | 4.4 (was 32.6)    | 4.1 (19.0)   | 18.5 (119)     |
  +-------+-------------------+--------------+----------------+
  | 1M    | 5.1 (was 63.7)    | 4.9 (37)     | 19.6 (235)     |
  +-------+-------------------+--------------+----------------+
  | 2M    | 6.6 (was 127)     | 6.4 (75)     | 21.3 (469)     |
  +-------+-------------------+--------------+----------------+
  | 5M    | 12.7 (was 317)    | 11.7 (191)   | 28.7 (1169)    |
  +-------+-------------------+--------------+----------------+

  i.e. circular plotting at 5 M variants is now **41x faster** than the
  pre-0.2.7 untrimmed path, and projects to ~38 s at 10 M variants
  (down from ~38 minutes — and faster than CMplot's circular path).

**Performance**

- Linear Manhattan rendering switched from ``ax.scatter`` (one ``PathCollection``
  carrying a path-per-point with per-point ``should_simplify`` checks) to
  one ``ax.plot(..., marker='.', linestyle='none')`` per chromosome
  (a single ``Line2D`` whose marker-draw loop is dramatically cheaper).
  Visually identical rasterised output; on a 1 M-variant single-track plot
  this alone shrinks ``plot_linearm`` from ~6 s to ~0.5 s.

- QQ plots (``plot_qq_single`` and ``plot_qq_combined``) make the same
  scatter → plot switch for the observed points.

- Chromosome-name normalisation in
  :func:`~pycmplot.io.get_sumstats_and_merged_sector_list` is now applied
  to the **categories** of the CHR ``Categorical`` (≤25 distinct values)
  rather than to the underlying N-row code array.  The result is stored
  as a ``Categorical`` ordered by ``CHROM_ORDER`` so downstream code can
  derive ``chr_idx`` from ``cat.codes`` directly.
- Linear-plot ``_prep`` recognises the canonical Categorical CHR column
  produced by the loader and skips the redundant ``str.replace +
  str.upper + replace`` pass that was running on every plot call.
- Optional CSV reader switched to ``engine='pyarrow'`` with safe fallback
  to the default C engine when pyarrow is unavailable.
- New ``compute_pvals`` parameter on
  :func:`~pycmplot.io.get_sumstats_and_merged_sector_list` (default
  ``True``); ``_core.py`` now sets it to ``False`` when no QQ plot is
  requested, skipping an ~80 MB-at-10 M-variants p-value-array copy that
  was unused on Manhattan- or circular-only runs.

Combined effect (measured, single-track untrimmed, fresh subprocess):

==========  ===========  ==========  ========
plot_type   500K before  500K after  speed-up
==========  ===========  ==========  ========
manhattan   32.6 s       4.6 s       7.1x
qq          19.0 s       6.7 s       2.8x
circular    119.0 s      39.9 s      3.0x
==========  ===========  ==========  ========

==========  ==========  ==========  ========
plot_type   1M before   1M after    speed-up
==========  ==========  ==========  ========
manhattan   63.7 s      6.0 s       10.6x
qq          37.1 s      10.2 s      3.6x
circular    235.3 s     73.3 s      3.2x
==========  ==========  ==========  ========

**Fixed**

- ``POS`` is now stored as plain ``int64`` after a ``to_numeric +
  dropna`` pass, rather than the nullable ``Int64`` that leaked ``pd.NA``
  into reductions like ``groupby(...).max()`` and caused
  ``TypeError: boolean value of NA is ambiguous`` further down the
  pipeline.
- ``plot_linearm``'s ``df.groupby(CHR)[POS].max()`` now passes
  ``observed=True`` so categorical chromosomes with no rows in a
  particular track produce no entry (``s.get(c, 0)`` handles the missing
  case), avoiding the ``NA``-propagation crash described above.
- Stripped 5 288 stray ``NUL`` bytes that had been appended to the end
  of ``pycmplot/plotting/linear.py`` (filesystem-level corruption from a
  partial overwrite — the file imported only after the trailing zeros
  were removed).

----

0.2.5 - 2026-04-20
------------------

**Fixed**

- Chromosome-22 positions falling outside hg38 chr22 limits after liftover
  no longer crash circular plotting.  The liftover post-filter now guards
  against unknown chromosomes with an informative warning.
- ``prep_pycmplot_input_info`` now resolves and stores column mappings
  **per file** rather than collapsing everything onto the last file.  This
  fixes incorrect column resolution when the input summary statistics files
  use different header names.
- ``io.get_file_header`` now correctly honours the ``delim`` argument when
  reading the header line.
- ``stats.get_highlight_snps`` now forwards ``logp`` through to
  ``get_lead_snps`` instead of hard-coding it to ``False`` — highlighting
  works correctly when plotting on the −log₁₀(p) axis.
- ``_core.py`` annotation resolution now uses the value of ``--annotate``
  (not the column name) when checking whether the requested annotation
  column exists in the hits table, and falls back to ``SNP`` with a warning
  when it does not.
- Chromosome-length sort (``--sort_track chrom_len``) now actually sorts by
  the number of chromosomes (most chromosomes first) rather than by track
  label.
- ``resources.ResourceConfig.require`` now imports ``as_file`` from
  ``importlib.resources`` so the bundled-resource fallback no longer raises
  ``NameError``. The fallback also now verifies that the resolved file
  actually exists before returning, rather than silently returning a
  phantom path.
- ``prep_pycmplot_input_info`` no longer emits a spurious "no build column
  detected" warning when the input files contain a ``BUILD`` column. The
  check previously inspected the length of the top-level info list, which
  only distinguishes the ``--build`` path from the no-build path; the fix
  also checks whether a build column was appended to ``old_cols``.
- Linear Manhattan plot: per-track labels and the shared
  ``-log₁₀(p-value)`` y-axis label no longer overlap in the left margin.
  Track labels are now rendered as a right-aligned sub-title above each
  axes (``ax.set_title(..., loc='right')``), which keeps them out of the
  data region entirely — so labels remain legible for dense null tracks,
  iHS/F_ST/XP-EHH panels, or any other plot where data can reach the
  upper-right corner.  The figure also reserves an explicit left strip
  for the shared y-label via ``fig.subplots_adjust`` instead of relying
  on ``tight_layout`` (which was incompatible with the shared-x gridspec
  and silently emitted a matplotlib warning).
- Linear Manhattan plot: the ``df = df[df[p_col] >= 0]`` sanity filter is
  now only applied when plotting ``-log₁₀(p)``. For non-p-value
  statistics (iHS, XP-EHH, Fay & Wu's H) negative values are legitimate
  and are preserved.  The filter was also previously applied *after*
  ``color_cycle`` was constructed, which caused a latent
  ``ValueError: 'c' argument has N elements, which is inconsistent with
  'x' and 'y'`` whenever the filter actually dropped rows.

- Annotation in circular plotting when **GENE** selected but **SNP** 
  annotated.


**Added**

- ``--ylabel`` / ``-yl`` flag (and ``ylabel=`` kwarg on
  :func:`~pycmplot.plotting.linear.plot_linear` and
  :func:`~pycmplot.plotting.linear.plot_linearm`) for overriding the
  shared y-axis label on linear Manhattan plots.  Intended for non-p-value
  statistics, e.g. ``--ylabel 'iHS'`` or ``--ylabel 'F_ST'``.
- All QQ-plotting functions (:func:`~pycmplot.plotting.qq.plot_qq_single`,
  :func:`~pycmplot.plotting.qq.plot_qq_combined`,
  :func:`~pycmplot.plotting.qq.plot_qq_separate`,
  :func:`~pycmplot.plotting.qq.plot_qq_overlay`) are now re-exported at the
  top level (``from pycmplot import plot_qq_combined``) and through the
  :mod:`pycmplot.plotting` subpackage.
- **hg18 → hg38 liftover.** ``BUILD`` column values of ``hg18`` (or
  ``--build hg18``) now trigger direct hg18 → hg38 coordinate conversion
  via a bundled UCSC chain file
  (``pycmplot/data/hg18ToHg38.over.chain.gz``). A new
  :func:`~pycmplot.liftover.liftover_hg18_to_hg38` helper and
  ``ResourceConfig.chain_hg18_hg38`` attribute (overridable via
  ``PYCMPLOT_CHAIN_HG18_HG38``) are exposed alongside the existing
  hg19 → hg38 path. Together these cover virtually all publicly available
  GWAS summary statistics.
- ``python -m pycmplot`` entry point (via a new ``__main__.py``).
- New Jupyter notebook demonstrating QQ-plotting workflows.
- All module-, class-, and function-level docstrings now use the bare
  ``"""..."""`` form so that Sphinx autodoc / numpydoc and
  :func:`help` render them correctly.
- Information about the sumstats printed to screen now includes number 
  of variants pre and post trimming, memory usage, and progress bar.


**Changed**

- Enhanced memory efficiency by changing **CHR** and **BUILD** columns 
  dtypes from ``str`` to ``category`` in ``io.py``

- Licence changed to MIT Licence.

----

0.2.2 - 2026-04-18
------------------

**Added**

QQ plots (:mod:`pycmplot.plotting.qq`):

- :func:`~pycmplot.plotting.qq.plot_qq_single` — single QQ panel on a
  provided axes, with 95% CI band, null diagonal, optional genome-wide
  line, and λ annotation.
- :func:`~pycmplot.plotting.qq.plot_qq_combined` — all sumstats as
  per-panel grid with configurable column count.
- :func:`~pycmplot.plotting.qq.plot_qq_separate` — one file per sumstat.
- :func:`~pycmplot.plotting.qq.plot_qq_overlay` — all sumstats on one
  shared axes, with λ in legend entries.
- :func:`~pycmplot.plotting.qq.thin_pvals` — log-uniform p-value thinning
  helper that preserves tail density while sparsifying the bulk, with no
  hard threshold seam.

CLI flags for QQ plotting:

=====================================  =======================================
Flag                                   Description
=====================================  =======================================
``-qq`` / ``--qq_plot``                Generate QQ plot(s) alongside the Manhattan plot.
``-qq_sep`` / ``--qq_separate``        Save one file per sumstat instead of a combined figure.
``-qq_ov`` / ``--qq_overlay``          Overlay all sumstats on a single QQ axes.
``-qq_cols`` / ``--qq_ncols``          Number of columns in the combined grid (default 3).
``-qq_max_pts`` / ``--qq_max_points``  Maximum points per track after thinning (default 50 000).
``-qq_thin`` / ``--qq_thin``           Enable log-uniform p-value thinning (off by default).
``-thin_below`` / ``--thin_below``     P-value floor below which all points are kept (default 0.01).
=====================================  =======================================

**Performance**

- Log-uniform thinning reduces a 10 M-SNP dataset to ≤ 50 000 plotted
  points with no perceptible visual difference.
- Scatter points are rasterised inside PDF/SVG output
  (``rasterized=True``), reducing file sizes from hundreds of MB to a
  few MB for large datasets.

**Fixed**

- ``_qq_arrays``: removed an erroneous reverse on the ``observed`` array
  that paired the largest expected quantile with the smallest observed
  p-value, breaking the diagonal.
- ``thin_pvals``: replaced the two-region split that could produce a
  zero bulk budget (silently dropping the diagonal below
  −log₁₀(p) = 2) with seamless log-uniform thinning.
- ``_plot_circularm``: increased padding between the first and last
  tracks to improve visibility of track labels and y-axis ticks.
- ``--build_column`` detection no longer fails when the flag is omitted.

----

0.2.1 — 2026-04-16
------------------

**Added**

- ``--build`` option for supplying per-file genome builds when the
  summary statistics files do not carry a ``BUILD`` column.
- ``--build`` and ``--build_column`` are both optional; plotting
  proceeds without genome-build information when neither is supplied.

**Changed**

- Expanded ``--annotate`` choices from ``snp``/``gene`` to any column in
  the hits table (and any column in a user-supplied annotation table in
  the Python API).

**Caveat**

- When multiple summary statistics files use different coordinate
  systems and ``--annotate`` is set, annotation defaults to hg38
  coordinates, which may mis-annotate hg19 variants.  Supplying correct
  builds avoids this.

----

0.1.9 — 2026-04-14
------------------

**Fixed**

- Column name auto-detection now covers both lower- and upper-case
  variants of every built-in candidate.
- ``build`` parameter of
  :func:`~pycmplot.io.prep_pycmplot_input_info` is now consistent with
  the CLI equivalent (required instead of optional).

----

0.1.8 — 2026-04-14
------------------

**Added**

- ``--highlight_color`` and ``--highlight_line_color`` options.
- Short form for ``--colors``.
- Long forms for ``-r_min``, ``-r_max``, ``-t_space``, ``-pad``.

**Fixed**

- ``from __future__ import annotations`` import bug.
- Short form for ``--highlight_line``.

----

0.1.0 — 2026-04-18
------------------

Initial release.

**Added**

Package structure:

- Installable Python package with ``src/`` layout, ``pyproject.toml``,
  ``setup.cfg``, and a ``setup.py`` compatibility shim.
- Console script ``pycmplot`` (also runnable as ``python -m pycmplot``).

Modules:

- :mod:`pycmplot.constants` — hg38 chromosome lengths, biotype priority
  weights, standard chromosome order.
- :mod:`pycmplot.resources` — :class:`~pycmplot.resources.ResourceConfig`
  dataclass for reference-file paths, configurable via environment
  variables (``PYCMPLOT_CHAIN_HG19_HG38``, ``PYCMPLOT_GENEINFO_HG38``,
  ``PYCMPLOT_GENEINFO_HG19``).
- :mod:`pycmplot.liftover` — lazy-initialised hg19 → hg38 coordinate
  conversion.
- :mod:`pycmplot.stats` — :func:`~pycmplot.stats.get_lead_snps` and
  :func:`~pycmplot.stats.get_highlight_snps`.
- :mod:`pycmplot.io` — summary statistics loader with auto-detection of
  delimiters and column names.
- :mod:`pycmplot.annotation` — strand-aware nearest-gene annotation
  with biotype-weighted prioritisation, promoter flagging, and
  :func:`~pycmplot.annotation.get_hits_summary_table`.
- :mod:`pycmplot.plotting.linear` — multi-track stacked linear
  Manhattan plotter.
- :mod:`pycmplot.plotting.circular` — multi-track Circos-style circular
  Manhattan plotter.
- :mod:`pycmplot.cli` — ``argparse`` CLI.
- :mod:`pycmplot._core` — the ``main()`` orchestration function.

**Fixed** (relative to the original monolithic script):

- Module-level ``LiftOver(hardcoded_path)`` call replaced by a lazy
  singleton; ``import pycmplot`` no longer raises ``FileNotFoundError``.
- Hardcoded ``/vast/awonkam1/...`` resourc
