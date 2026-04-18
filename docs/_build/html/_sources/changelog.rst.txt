.. include:: ./changelog.rst

Changelog
=========

All notable changes to **pycmplot** are documented here.
Format follows `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_.
Versioning follows `Semantic Versioning <https://semver.org/>`_.

----

0.2.2 — 2026-04-18
-------------------

Added


**QQ plots** (``pycmplot.plotting.qq``)

``plot_qq_single()``
    Draws a single QQ plot onto a provided Matplotlib ``Axes``.
    Includes a 95% confidence interval band (derived from the beta
    distribution), diagonal null line, optional genome-wide significance
    line, and genomic inflation factor λ annotation.

``plot_qq_combined()``
    Plots all sumstats as individual panels arranged in a configurable
    column grid in a single figure.

``plot_qq_separate()``
    Saves one QQ plot file per sumstat, named
    ``{stem}_{label}.{format}``.

``plot_qq_overlay()``
    Plots all sumstats on a single shared axes, each coloured by label
    with λ values embedded in the legend entries.

``thin_pvals()``
    Log-uniform p-value thinning helper for fast QQ plotting.
    Selects up to ``max_points`` evenly-spaced positions along the
    −log₁₀(p) axis so the null bulk is sparse and the significant tail
    is dense, with no hard threshold boundary and no visible seam.
    Lambda is always computed on the full unfiltered array before
    thinning.

**CLI flags**

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Flag
     - Description
   * - ``-qq`` / ``--qq_plot``
     - Generate QQ plot(s) alongside the Manhattan plot.
   * - ``-qq_sep`` / ``--qq_separate``
     - Save one file per sumstat instead of a combined figure.
   * - ``-qq_ov`` / ``--qq_overlay``
     - Overlay all sumstats on a single QQ axes.
   * - ``-qq_cols`` / ``--qq_ncols``
     - Number of columns in the combined grid (default: 3).
   * - ``-qq_max_pts`` / ``--qq_max_points``
     - Maximum points per track after thinning (default: 50 000).
   * - ``-qq_no_thin`` / ``--qq_no_thin``
     - Disable thinning and plot all points (slow for large datasets).

**Performance**

- **Log-uniform p-value thinning** reduces a 10 M-SNP dataset to
  ≤ 50 000 plotted points with no perceptible visual difference.
  Thinning is on by default for all QQ functions and can be disabled
  with ``thin=False`` (API) or ``--qq_no_thin`` (CLI).
- **Rasterised scatter** (``rasterized=True``) renders the point cloud
  as a bitmap inside PDF/SVG output, reducing file sizes from hundreds
  of MB to a few MB for large datasets.

Fixed

``_qq_arrays`` — reversed observed array
    The ``observed`` array was reversed with ``[::-1]``, pairing the
    largest expected quantile with the smallest observed value and
    breaking the diagonal entirely.  Removed the reversal so rank *i*
    correctly maps to the *i*-th smallest p-value.

``thin_pvals`` — zero bulk budget
    When the tail region alone exceeded ``max_points``, the bulk budget
    became zero and all null-region points were silently dropped, leaving
    the diagonal invisible below −log₁₀(p) = 2.  Replaced the two-region
    split with seamless log-uniform thinning, eliminating both the
    zero-budget bug and the visible seam at the threshold boundary.

  ``_plot_circularm`` — track label space
    Increased space between first and last tracks to improve visibility of 
    track labels, and y-axis and its label

----


Version 0.2.1 (2026-04-16)
==========================

Added


- Added option to supply genome builds for summary stats files ``--build``` if 'BUILD' column is not in the files.
      This is a revertion to earlier versions of the plotting script. 
      Also made ``--build`` and ``--build_column`` optional allowing plotting to still proceed without genome build information.
      However, caution must be taken when multiple summary stats files are provided with different coordinate systems.
      For example, if ``--annotate`` is set, hits table generation will default to `hg38` coordinate, potentially leading to 
      in accurate annotations for variants in different coordinates.


Changed


- Changed ``--annotate`` choices:
      Expnaded choices from `snp` and `gene` to include other columns in hits table.
      Also allowed for other columns in user supplied annotation table (available in python API only).


Version 0.1.9 (2026-04-14)
==========================

Fixed


- Fixed column name auto-detection:
   - Expanded candidates list by adding lower and upper case versions for existing condidates.

  Fixed ``build`` option for ``prep_pycmplot_input_info`` function.
   - Updated it from optional to required parameter to be consistent with command line version.


Version 0.1.8 (2026-04-14)
==========================

Added


- Added options to specify color of highlighted positions ``highlight_color``
  and line running through highlighted positions ``highlight_line_color``.

  Added command-line short forms for ``--colors``.

  Added command-line long forms for ``-r_min``, ``-r_max``, ``-t_space``, and
  ``-pad`` (`#1 <https://github.com/esohkevin/pycmplot/issues/1>`_)


Fixed


- Fixed bug with __future__ import.

  Fixed command-line short form for ``--highlight_line``. (`#2 <https://github.com/esohkevin/pycmplot/issues/2>`_)


----


0.1.0 — 2026-04-18
-------------------

Initial release.

Added


**Package structure**

- Installable Python package with ``src/`` layout, ``pyproject.toml``,
  ``setup.cfg``, and a ``setup.py`` compatibility shim for older
  setuptools.
- Console script entry point: ``pycmplot`` (also runnable as
  ``python -m pycmplot``).

**Modules**

``constants``
    hg38 chromosome lengths, biotype priority weights, and standard
    chromosome order.

``resources``
    ``ResourceConfig`` dataclass for configuring external reference file
    paths (liftover chain, gene info files) via constructor arguments or
    environment variables (``PYCMPLOT_CHAIN_HG19_HG38``,
    ``PYCMPLOT_GENEINFO_HG38``, ``PYCMPLOT_GENEINFO_HG19``,
    ``PYCMPLOT_FEATURESINFO``).

``liftover``
    Lazy-initialised hg19 → hg38 coordinate conversion.  The
    ``LiftOver`` object is created on first use rather than at import
    time, preventing ``FileNotFoundError`` on systems where the chain
    file is not configured.

``stats``
    ``get_lead_snps()`` (greedy distance clumping) and
    ``get_highlight_snps()`` (locus window flagging).

``io``
    Summary statistics loading, auto-detection of file delimiters, and
    ``get_sumstats_and_merged_sector_list()`` for Circos sector-size
    computation.

``annotation``
    Strand-aware nearest-gene annotation with biotype-weighted
    prioritisation, promoter flagging, and ``get_hits_summary_table()``.

``plotting.linear``
    ``multi_track_linear_manhattan()``: multi-track linear Manhattan
    plot with per-track significance lines, locus highlighting, and
    gene/SNP annotation arrows.

``plotting.circular``
    ``plot_circosm()``: per-chromosome Circos track plotter.
    ``compute_track_radii_dict()``: dynamic track radius calculator.

``cli``
    Full ``argparse`` CLI with required and optional argument groups.

``_core``
    ``main()`` orchestration function wiring CLI → data loading →
    plotting.

Fixed


Relative to the original monolithic script:

- Module-level ``LiftOver(hardcoded_path)`` call replaced with a lazy
  singleton so importing the package no longer raises
  ``FileNotFoundError``.
- Three hardcoded ``/vast/awonkam1/...`` resource paths replaced with
  ``ResourceConfig``.
- ``highlight`` was used as a free variable inside
  ``get_sumstats_and_merged_sector_list`` — promoted to an explicit
  parameter.
- ``geneinfo`` was used as an implicit closure variable inside
  ``build_locus_summary`` — now passed explicitly.
- ``suggest_line`` was used before assignment in certain code paths —
  initialisation reordered.



