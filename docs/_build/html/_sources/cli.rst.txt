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
   * - ``-b``, ``--build_column``
     - *auto*
     - Name of the genome build column (values: ``hg19``, ``hg38``).

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
   * - ``-a``, ``--annotate``
     - ``SNP``
     - Annotate significant hits with ``SNP`` IDs or nearest ``GENE`` names.
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
   * - ``-qq``, ``--qq_plot``
     - off
     - Generate a QQ-plot alongside the Manhattan plot *(coming soon)*.

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
