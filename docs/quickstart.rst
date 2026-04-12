.. _quickstart:

Quickstart
==========

This page walks through the most common workflows in a few lines each.
For a detailed walkthrough with real data, see :ref:`python_api_notebook`.

Input file format
-----------------

pycmplot accepts any whitespace- or comma-delimited summary statistics file,
including gzip-compressed (``.gz``) files. Required columns (auto-detected from
common names) are:

- Chromosome (e.g. ``CHR``, ``CHROM``, ``#CHROM``)
- Base-pair position (e.g. ``BP``, ``POS``, ``pos``)
- Variant identifier (e.g. ``SNP``, ``RSID``, ``MarkerName``)
- P-value or test statistic (e.g. ``P``, ``pvalue``, ``Wald_P``)

Optionally, a genome-build column (``hg19`` / ``hg38``) enables automatic
liftover (see :ref:`cli_liftover`).

.. tip::
   For large summary statistics files, always pass ``--trim_pval 0.01`` to
   discard variants with p > 0.01 before plotting. This can reduce memory usage
   by an order of magnitude.

Linear Manhattan plot (single trait)
-------------------------------------

**Command line**

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

**Python API**

.. code-block:: python

   from pycmplot import plot_linear

   plot_linear(
       sum_stats=["HbF.tsv.gz"],
       labels=["HbF"],
       logp=True,
       signif_line=True,
       highlight=True,
       annotate="GENE",
       output_dir="./results",
       output_format="png",
       dpi=300,
   )

Multi-track linear Manhattan plot
----------------------------------

Compare signals across three RBC traits in a single stacked figure:

.. code-block:: python

   from pycmplot import plot_linear

   plot_linear(
       sum_stats=["HbF.tsv.gz", "MCV.txt.gz", "MCH.tsv.gz"],
       labels=["HbF", "MCV", "MCH"],
       logp=True,
       signif_line=True,
       highlight=True,
       annotate="GENE",
       output_dir="./results",
   )

Circular (Circos-style) Manhattan plot
----------------------------------------

.. code-block:: python

   from pycmplot import plot_circular

   plot_circular(
       sum_stats=["HbF.tsv.gz", "MCV.tsv.gz"],
       labels=["HbF", "MCV"],
       trim_pval=0.01,
       logp=True,
       signif_threshold=True,
       plot_title="RBC Traits",
       output_dir="./results",
   )

Mixed genome builds with liftover
-----------------------------------

If your summary statistics were generated on different reference panels (hg19
and hg38), add a build column to each file, concatenate them, and pass the
column name:

.. code-block:: bash

   pycmplot \
     --sum_stats combined.tsv.gz \
     --labels MyTrait \
     --build_column BUILD \
     --logp \
     --output_dir ./results

pycmplot will liftover all hg19 positions to hg38 before plotting.

Next steps
----------

- Full CLI reference: :ref:`cli`
- Complete Python API: :ref:`api`
- Interactive notebook tutorial: :ref:`python_api_notebook`
