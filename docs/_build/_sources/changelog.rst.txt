Changelog
=========

All notable changes to **pycmplot** are documented here.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_
and this project adheres to `Semantic Versioning <https://semver.org/>`_.

----

0.2.7 — 2026-04-27
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

0.2.5 — 2026-04-20
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

0.2.2 — 2026-04-18
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