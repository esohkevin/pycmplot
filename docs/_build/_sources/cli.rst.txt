.. _cli:

Command-Line Interface
======================

pycmplot exposes a full command-line interface that mirrors the Python API.
After installation, the ``pycmplot`` command is available in your PATH.

.. code-block:: bash

   pycmplot -h    # print full help

Input / output options
----------------------

.. list-table::
   :widths: 30 15 55
   :header-rows: 1

   * - Flag
     - Default
     - Description
   * - ``-s``, ``--sum_stats``
     - *required*
     - Comma-separated list of summary statistics files (gzip supported).
   * - ``-l``, ``--labels``
     - *required*
     - Comma-separated track labels, one per file in ``--sum_stats``.
   * - ``-od``, ``--output_dir``
     - ``.``
     - Output directory. Created if it does not exist.
   * - ``-of``, ``--output_format``
     - ``png``
     - Output image format: ``png``, ``pdf``, ``svg``, or ``jpg``.
   * - ``--dpi``
     - ``300``
     - Resolution in dots per inch (raster formats only).

Column auto-detection options
------------------------------

pycmplot infers column names automatically. Use these flags only if your
column names fall outside the recognised defaults.

.. list-table::
   :widths: 30 15 55
   :header-rows: 1

   * - Flag
     - Default
     - Description
   * - ``-chr``, ``--chrom_column``
     - *auto*
     - Name of the chromosome column.
   * - ``-pos``, ``--pos_column``
     - *auto*
     - Name of the base-pair position column.
   * - ``-snp``, ``--snp_column``
     - *auto*
     - Name of the variant/marker ID column.
   * - ``-p``, ``--pval_column``
     - *auto*
     - Name of the p-value (or test statistic) column.
   * - ``-d``, ``--delim``
     - *auto*
     - Field delimiter (``tab``, ``space``, ``comma``, ``colon``,
       ``semi-colon``). Auto-detected when omitted.
   * - ``-bc``, ``--build_column``
     - *auto*
     - Name of an in-file genome build column (values: ``hg18``, ``hg19``,
       ``hg38``).
   * - ``-b``, ``--build``
     - *none*
     - Comma-separated list of per-file genome builds in the same order as
       ``--sum_stats``, used when the files have no build column
       (e.g. ``hg19,hg38,hg18``). ``hg18`` and ``hg19`` coordinates are
       lifted to ``hg38`` automatically before plotting.

.. _cli_liftover:

Plotting behaviour options
--------------------------

.. list-table::
   :widths: 30 15 55
   :header-rows: 1

   * - Flag
     - Default
     - Description
   * - ``-m``, ``--mode``
     - ``lm``
     - Plot mode: ``lm`` (linear) or ``cm`` (circular).
   * - ``--logp``
     - off
     - Plot –log₁₀(p) instead of raw p-values.
   * - ``-sig``, ``--signif_threshold``
     - off
     - Draw a genome-wide significance threshold line (auto 0.05/N).
   * - ``-sigl``, ``--signif_line``
     - ``5e-8``
     - Explicit significance threshold value; overrides ``--signif_threshold``.
   * - ``-sug``, ``--suggest_threshold``
     - off
     - Draw a suggestive significance threshold line.
   * - ``-hl``, ``--highlight``
     - off
     - Highlight significant loci with a distinct colour.
   * - ``-hc``, ``--highlight_color``
     - ``brown``
     - Colour for highlighted variants.
   * - ``-hll``, ``--highlight_line``
     - off
     - Draw a vertical dashed line through each highlighted locus.
   * - ``-hlc``, ``--highlight_line_color``
     - ``grey``
     - Colour of the highlight line.
   * - ``-ht``, ``--highlight_thresh``
     - ``5e-8``
     - P-value threshold for highlighting.
   * - ``-a``, ``--annotate``
     - ``SNP``
     - Annotate significant hits with ``SNP`` IDs, nearest ``GENE`` names,
       or any column available in the hits table (e.g. ``top_gene``,
       ``nearest_upstream_gene``).
   * - ``-tp``, ``--trim_pval``
     - off
     - Exclude variants with p-value above this threshold (e.g. ``0.01``)
       before plotting. **Strongly recommended for large files.**
   * - ``-st``, ``--sort_track``
     - input order
     - Sort tracks by ``label`` (alphabetical) or ``chrom_len``
       (chromosome length, longest first).
   * - ``--plot_title``
     - *(none)*
     - Title string printed above the plot.

QQ plot options
---------------

.. list-table::
   :widths: 30 15 55
   :header-rows: 1

   * - Flag
     - Default
     - Description
   * - ``-qq``, ``--qq_plot``
     - off
     - Generate QQ plot(s) alongside the Manhattan plot.
   * - ``-qq_sep``, ``--qq_separate``
     - off
     - Save one QQ-plot file per summary statistics file.
   * - ``-qq_ov``, ``--qq_overlay``
     - off
     - Overlay all sumstats on a single shared QQ axes, with λ in the legend.
   * - ``-qq_cols``, ``--qq_ncols``
     - ``3``
     - Number of columns in the combined QQ-plot grid.
   * - ``-qq_thin``, ``--qq_thin``
     - off
     - Enable log-uniform p-value thinning for fast QQ plotting.
   * - ``-thin_below``, ``--thin_below``
     - ``0.01``
     - P-value threshold below which all points are kept (points above are
       downsampled when thinning is enabled).
   * - ``-qq_max_pts``, ``--qq_max_points``
     - ``50000``
     - Maximum points plotted per QQ track after thinning.

Example commands
----------------

Single-trait linear Manhattan plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pycmplot \
     --sum_stats HbF.tsv.gz \
     --labels HbF \
     --logp \
     --signif_line \
     --highlight \
     --annotate GENE \
     --output_dir ./results \
     --output_format png \
     --dpi 300

Three-trait stacked linear Manhattan plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pycmplot \
     --sum_stats HbF.tsv.gz,MCV.txt.gz,MCH.tsv.gz \
     --labels HbF,MCV,MCH \
     --logp \
     --signif_line \
     --highlight \
     --annotate GENE \
     --trim_pval 0.01 \
     --output_dir ./results

Circular Manhattan plot
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pycmplot \
     --sum_stats HbF.tsv.gz,MCV.tsv.gz \
     --labels HbF,MCV \
     --mode cm \
     --trim_pval 0.01 \
     --logp \
     --signif_threshold \
     --plot_title "RBC Traits" \
     --output_dir ./results

Manhattan plot with companion QQ plots
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pycmplot \
     --sum_stats HbF.tsv.gz,MCV.tsv.gz,MCH.tsv.gz \
     --labels HbF,MCV,MCH \
     --logp \
     --signif_line \
     --annotate GENE \
     --qq_plot \
     --qq_thin \
     --qq_ncols 3 \
     --output_dir ./results

Supplying per-file genome builds
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When summary statistics files do not carry a ``BUILD`` column, supply the
builds in the same order as ``--sum_stats``:

.. code-block:: bash

   pycmplot \
     --sum_stats hg19_study.tsv.gz,hg38_study.tsv.gz \
     --labels Study_A,Study_B \
     --build hg19,hg38 \
     --logp \
     --annotate GENE \
     --output_dir ./results
