Changelog
=========

All notable changes to **pycmplot** are documented here.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_
and this project adheres to `Semantic Versioning <https://semver.org/>`_.

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

**Docs**

- Tutorial note under :ref:`cli-tut-linear` explaining signed-stat
  behavior for iHS / XP-EHH and how ``signif_threshold`` /
  ``highlight_thresh`` are applied on ``|value|``.


0.4.1 - 2026-08-23
------------------------------------------------------------------------------

**Added**

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

**Added**

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
