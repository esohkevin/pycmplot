.. _tutorial:

Tutorial
========

This tutorial is an end-to-end walkthrough of every user-facing feature
in **pycmplot**, working from the same synthetic dataset throughout so
that every snippet is directly runnable.  It's longer and more granular
than :ref:`quickstart` — read the quickstart first if you just want to
see the shape of the API in a page.

The tutorial is organised as a progression:

#. :ref:`tut-setup` — a synthetic sumstats file you can paste-and-run.
#. :ref:`tut-load` — the three-step load / prep / plot mental model.
#. :ref:`tut-linear` — single- and multi-track linear Manhattan plots.
#. :ref:`tut-circular` — the Circos-style circular layout.
#. :ref:`tut-highlight` — thresholds, colours, and highlight controls.
#. :ref:`tut-annotate` — automatic nearest-gene labels and the hits table.
#. :ref:`tut-qq` — QQ plots and the ``compute_pvals`` opt-in.
#. :ref:`tut-cache` — the per-track cache and warm-resume behaviour.
#. :ref:`tut-overlay` — the user-editable hits-overlay TSV.
#. :ref:`tut-colors` — per-locus highlight colours.
#. :ref:`tut-categories` — per-locus categories and the custom legend.
#. :ref:`tut-multipanel` — multiple sumstats groups on one canvas.
#. :ref:`tut-liftover` — mixed genome builds.
#. :ref:`tut-cli` — CLI equivalents for every step above.
#. :ref:`tut-faq` — troubleshooting and common gotchas.


.. _tut-setup:

Setup: a synthetic dataset
--------------------------

Every snippet below assumes the following synthetic sumstats file on
disk.  Save it once and reuse throughout the tutorial:

.. code-block:: python

   import numpy as np
   import pandas as pd

   rng = np.random.default_rng(0)
   n = 200_000
   df = pd.DataFrame({
       "CHR": rng.choice([str(i) for i in range(1, 23)], size=n),
       "BP":  rng.integers(1, 200_000_000, size=n),
       "SNP": [f"rs{i}" for i in range(n)],
       "P":   rng.uniform(1e-12, 1.0, size=n),
   })
   # Plant ~20 genome-wide-significant hits
   hits_idx = rng.choice(n, 20, replace=False)
   df.loc[hits_idx, "P"] = 10 ** -rng.uniform(9, 15, 20)
   df.to_csv("hb.tsv", sep="\t", index=False)

   # A second trait, for multi-track examples
   df2 = df.sample(frac=1.0, random_state=1).reset_index(drop=True)
   df2["P"] = rng.uniform(1e-12, 1.0, size=n)
   df2.loc[rng.choice(n, 15, replace=False), "P"] = \
       10 ** -rng.uniform(9, 15, 15)
   df2.to_csv("mcv.tsv", sep="\t", index=False)


.. _tut-load:

The three-step mental model
---------------------------

pycmplot separates loading from rendering so that expensive per-file
work (I/O, liftover, lead extraction, hits-table construction) happens
exactly once even when you produce multiple plot types (linear +
circular + QQ) from the same data.  Every workflow in this tutorial
uses the same three steps:

#. :func:`~pycmplot.io.prep_pycmplot_input_info` — resolve column
   names and delimiters for each input file.
#. :func:`~pycmplot.io.get_sumstats_and_merged_sector_list` — load,
   trim, lift over (if needed), extract leads, build the hits table.
   Returns a ``bundle`` dict.
#. A plotter (:func:`~pycmplot.plotting.linear.plot_linear`,
   :func:`~pycmplot.plotting.circular.plot_circular`, or one of the
   QQ plotters) — consumes fields from ``bundle``.

Here is the canonical shape:

.. code-block:: python

   from pycmplot import (
       prep_pycmplot_input_info,
       get_sumstats_and_merged_sector_list,
   )

   files, labels = ["hb.tsv"], ["Hb"]
   file_info = prep_pycmplot_input_info(sum_stats=files, labels=labels)
   bundle = get_sumstats_and_merged_sector_list(
       sum_stats=files, labels=labels, file_info=file_info,
       logp=True, trim_pval=0.01, signif_threshold=5e-8,
   )
   # bundle has: 'dfs', 'sectors', 'annot', 'lines', 'pvals'
   print(list(bundle))


.. _tut-linear:

Linear Manhattan plots
----------------------

Single-track
~~~~~~~~~~~~

.. code-block:: python

   from pycmplot import plot_linear

   plot_linear(
       sumstats_loaded=bundle["dfs"],
       signif_lines=bundle["lines"],
       hits_table=bundle["annot"],
       logp=True,
       colors=["steelblue", "silver"],
       plot_title="Hb",
       output_dir="./out", output_format="png", dpi=300,
   )

Multi-track
~~~~~~~~~~~

Passing more than one sumstats file stacks the tracks vertically, one
axes per file, sharing the chromosomal x-axis:

.. code-block:: python

   files, labels = ["hb.tsv", "mcv.tsv"], ["Hb", "MCV"]
   file_info = prep_pycmplot_input_info(sum_stats=files, labels=labels)
   bundle = get_sumstats_and_merged_sector_list(
       sum_stats=files, labels=labels, file_info=file_info,
       logp=True, trim_pval=0.01, signif_threshold=5e-8,
   )
   plot_linear(
       sumstats_loaded=bundle["dfs"],
       signif_lines=bundle["lines"],
       hits_table=bundle["annot"],
       logp=True,
       colors=["steelblue", "silver"],
       output_dir="./out",
   )

Colour scheme is applied per chromosome (alternating), not per track,
matching classic Manhattan convention.


.. _tut-circular:

Circular (Circos) plots
-----------------------

Same ``bundle``, different plotter — plus one extra input,
``sector_sizes``, which controls the sector layout:

.. code-block:: python

   from pycmplot import plot_circular

   plot_circular(
       sumstats_loaded=bundle["dfs"],
       sector_sizes=bundle["sectors"],
       signif_lines=bundle["lines"],
       hits_table=bundle["annot"],
       logp=True,
       colors=["steelblue", "silver"],
       plot_title="RBC Traits",
       output_dir="./out",
   )

If you plan to render both linear and circular from the same data,
call the loader **once** and pass the same ``bundle`` to both — the
loader is the expensive step.


.. _tut-highlight:

Highlighting
------------

Set ``highlight=True`` at both load and plot time to colour signals
above a threshold:

.. code-block:: python

   bundle = get_sumstats_and_merged_sector_list(
       sum_stats=files, labels=labels, file_info=file_info,
       logp=True, trim_pval=0.01,
       highlight=True,             # extract hits during load
       highlight_thresh=5e-8,      # p-value cutoff
       signif_threshold=5e-8,
   )
   plot_linear(
       sumstats_loaded=bundle["dfs"],
       hits_table=bundle["annot"],
       signif_lines=bundle["lines"],
       highlight=True,             # render coloured
       highlight_color="brown",    # default colour
       logp=True, output_dir="./out",
   )

The default ``highlight_color`` applies to every hit unless a specific
locus has been given a custom colour in the hits overlay
(see :ref:`tut-colors`).


.. _tut-annotate:

Automatic gene annotation
-------------------------

Add ``annotate="GENE"`` to label each lead SNP with its nearest gene:

.. code-block:: python

   plot_linear(
       sumstats_loaded=bundle["dfs"],
       hits_table=bundle["annot"],
       signif_lines=bundle["lines"],
       highlight=True, logp=True,
       annotate="GENE", label_col="top_gene",
       output_dir="./out",
   )

The hits table produced by the loader (``bundle["annot"]``) is a
tidy summary you can inspect or export directly:

.. code-block:: python

   print(bundle["annot"].head())
   # CHR   POS       SNP    P              LABEL   top_gene   source
   # 1     14822344  rs123  3.4e-11        Hb      FTO        auto
   # ...

Columns you'll see include ``CHR``, ``POS``, ``SNP``, ``P``, ``LABEL``
(track name), ``top_gene``, plus the overlay columns
``source``, ``highlight_color``, and ``category`` (introduced by the
overlay system — see :ref:`tut-overlay`).


.. _tut-qq:

QQ plots and ``compute_pvals``
------------------------------

.. important::

   As of pycmplot 0.4.0, the loader **does not** materialise the full
   sorted p-value array by default.  Pass ``compute_pvals=True`` if
   you plan to make a QQ plot:

.. code-block:: python

   from pycmplot import plot_qq_combined, plot_qq_overlay, plot_qq_separate

   bundle = get_sumstats_and_merged_sector_list(
       sum_stats=files, labels=labels, file_info=file_info,
       logp=True, trim_pval=0.01,
       compute_pvals=True,          # <-- opt in
   )

   plot_qq_combined(
       pval_dict=bundle["pvals"],
       thin=True, max_points=50_000,
       ncols=2, title="RBC",
       output_path="./out/rbc_qq", fig_format="png",
   )

   plot_qq_overlay(
       pval_dict=bundle["pvals"],
       thin=True, max_points=50_000,
       title="RBC", output_path="./out/rbc_qq_overlay",
   )

   plot_qq_separate(
       pval_dict=bundle["pvals"], base_name="RBC",
       thin=True, max_points=50_000,
       output_path="./out/rbc_qq",
   )

Every QQ plotter draws the 95% CI band around the diagonal and
annotates each track with its genomic inflation factor λ.  The CI
band is not in the legend — it's visually self-evident, and removing
it from the legend keeps the top corner uncluttered.


.. _tut-cache:

Caching and warm resume
-----------------------

Loading is the expensive step (I/O + trim + liftover + lead
extraction).  Turn on caching and every re-run of the same
``(files, parameters)`` combination completes in milliseconds:

.. code-block:: python

   bundle = get_sumstats_and_merged_sector_list(
       sum_stats=files, labels=labels, file_info=file_info,
       logp=True, trim_pval=0.01, highlight=True,
       cache=True,                # enable
       cache_dir="./.pycmplot",   # where to store artefacts
       resume=True,               # reuse existing cache entries
   )

**How keying works.** Each track gets a cache entry keyed on
``SHA-256(raw_file_sha256 + cache_version + Stage-1 params)``.  If any
of those change — you re-genotype and re-run, or you bump
``highlight_thresh`` — the affected tracks silently regenerate and
the rest are reused.  There are no stale-cache accidents.

**Cache layout** (rooted at ``cache_dir``):

.. code-block:: text

   .pycmplot/
   ├── tracks/
   │   ├── Hb.<key>.parquet         # loaded rows
   │   ├── Hb.<key>.leads.parquet   # per-track hits
   │   ├── Hb.<key>.pvals.npz       # sorted float32 pvals (only if compute_pvals=True)
   │   └── Hb.<key>.meta.json
   └── annotations/
       └── hits.<group_key>.tsv     # merged hits overlay (see next section)

**Disabling resume.** ``resume=False`` forces regeneration but still
writes the fresh entries to disk (useful for benchmarking or after
you've bumped ``CACHE_VERSION`` upstream).

**Clearing the cache.** Delete ``cache_dir`` on the filesystem, or
use ``pycmplot --clear_cache --cache_dir ./.pycmplot`` from the CLI.


.. _tut-overlay:

The hits overlay TSV
--------------------

When caching is enabled, the loader writes the hits table to a
group-scoped TSV that **you're expected to hand-edit**.  The path is:

.. code-block:: text

   <cache_dir>/annotations/hits.<group_key>.tsv

where ``group_key`` is a short SHA-256 of the sorted list of per-track
cache keys — so distinct ``(files, parameters)`` combinations get
distinct overlay files and never clobber each other, even on the same
canvas.

Every row starts with ``source='auto'`` (populated from the sumstats)
and three overlay columns you can edit:

.. code-block:: text

   # pycmplot hits overlay
   # Edit any row's highlight_color / category to customise the plot.
   # source=user rows are added by you; source=auto rows come from the loader.
   CHR  POS       SNP    P        LABEL  top_gene  source  highlight_color  category
   1    14822344  rs123  3.4e-11  Hb     FTO       auto    auto             significant
   2    91344001  rs456  2.1e-10  Hb     MYADM     auto    auto             significant
   ...

You can:

- **Edit ``highlight_color``** to give any locus its own colour
  (see :ref:`tut-colors`).
- **Edit ``category``** to give any locus a legend label
  (see :ref:`tut-categories`).
- **Add rows with ``source=user``** to force annotation of loci that
  didn't make the automatic cutoff.

**Persistence across cache regenerations.**  When Stage-1 parameters
change and auto rows are re-derived, your custom values are inherited
by ``(CHR, POS)`` lookup — you don't have to also change
``source='auto'`` to ``user`` to keep your edits.


Example batch edit ``hits.<group_key>.tsv`` files using ``awk`` in commandline:

.. code-block:: bash

   cachedir=/your/cachedir
   
   for i in ${cachdir}/annotations/hits.*.tsv; do 
      awk '
         OFS="\t" 
         {
            if($1 ~ /^#/) {print $0} 
            else{
               if($1 ~ /^source/) {print $0} 
               else{
                  if($5 >= 5e-08) {$26="orange"; $27="marginally significant (P < 1e-07)"} 
                  else {$27="significant (P < 5e-08)"} {print $0}
               }
            }
         }' ${i} > ${i}.bak

      mv ${i}.bak ${i}
   done

- This highlights all signals with ``P < 5e-08`` with the defaul ``brown`` color and all signals
not reaching the genome-wide significance threshold but have P < 1e-07 with orange.
- It updates the categories for both to ``significant (P < 5e-08)`` and ``marginally significant (P < 1e-07)``
respectively. This would be used to add a custom legend to the Manhattan plots, making it self explanatory.

.. _tut-colors:

Per-locus highlight colours
---------------------------

Set the ``highlight_color`` column on any row to a matplotlib-parseable
colour (name, ``#rrggbb`` hex, or an ``(r, g, b)`` tuple):

.. code-block:: text

   CHR  POS       ...  source  highlight_color  category
   1    14822344  ...  auto    red              significant
   2    91344001  ...  auto    #00cc44          significant
   3    50000000  ...  auto    auto             significant  # falls back to default

The sentinel ``auto`` (or blank / NaN / invalid) falls back to the
plot-time ``highlight_color`` argument.  Invalid colours emit a
warning naming the offending value so typos are easy to fix.

Matching is done by nearest-lead within 500 kb on the same
chromosome, so re-running with a different thinning or trim never
loses your colour choices.


.. _tut-categories:

Per-locus categories and the custom legend
------------------------------------------

The ``category`` column lets you group loci in the legend:

.. code-block:: text

   CHR  POS       ...  highlight_color  category
   1    14822344  ...  red              novel
   2    91344001  ...  #00cc44          replicated
   3    50000000  ...  auto             significant

**Both the linear and circular plotters** render a "Highlighted
Categories" legend with one entry per unique category in
first-appearance order — reorder your rows in the TSV to reorder the
legend.  When everything is left at defaults (all rows
``category=significant`` and ``highlight_color=auto``), **no legend
is added**, preserving the pre-feature layout.

Same-category-different-colour is deduplicated: the first colour
wins and a warning names the loser so you can fix the conflict.

Missing columns (e.g. loading a legacy cache) are tolerated — the
helper simply returns no legend entries in that case.


.. _tut-multipanel:

Multi-panel canvas
------------------

To place multiple *groups* of sumstats on the same figure — each
group is its own stacked linear (or circular) sub-plot — pass an
explicit matplotlib ``Axes`` (or ``SubFigure``) via ``ax=``:

.. code-block:: python

   import matplotlib.pyplot as plt

   fig = plt.figure(figsize=(14, 8), constrained_layout=True)
   sub_top, sub_bot = fig.subfigures(2, 1)

   for sub, group_files, group_labels in [
       (sub_top, ["hb.tsv"],  ["Hb"]),
       (sub_bot, ["mcv.tsv"], ["MCV"]),
   ]:
       fi = prep_pycmplot_input_info(
           sum_stats=group_files, labels=group_labels,
       )
       b = get_sumstats_and_merged_sector_list(
           sum_stats=group_files, labels=group_labels, file_info=fi,
           logp=True, highlight=True, cache=True, cache_dir="./.pycmplot",
       )
       plot_linear(
           sumstats_loaded=b["dfs"],
           hits_table=b["annot"], signif_lines=b["lines"],
           logp=True, highlight=True,
           ax=sub,  # <-- render into this subfigure
       )

   fig.savefig("./out/two_panel.png", dpi=300)

Every cache file, hits overlay, per-locus colour, and category legend
is *group-scoped*, so two panels with different sumstats never
clobber each other.


.. _tut-liftover:

Mixed genome builds
-------------------

If your files were generated on different reference panels, pycmplot
can lift over hg18 and hg19 coordinates to hg38 before plotting.
Supply the builds either through a ``BUILD`` column in the files or
by passing ``build_list=`` (Python) / ``--build`` (CLI):

.. code-block:: python

   bundle = get_sumstats_and_merged_sector_list(
       sum_stats=["study_hg18.tsv", "study_hg19.tsv", "study_hg38.tsv"],
       labels=["A", "B", "C"],
       build_list=["hg18", "hg19", "hg38"],
       logp=True, trim_pval=0.01,
       file_info=prep_pycmplot_input_info(
           sum_stats=["study_hg18.tsv", "study_hg19.tsv", "study_hg38.tsv"],
           labels=["A", "B", "C"],
       ),
   )

Liftover is cached alongside the loaded rows, so subsequent runs
(with ``cache=True``) skip it entirely.


.. _tut-cli:

CLI equivalents
---------------

Everything above has a CLI equivalent.  A representative one-liner:

.. code-block:: bash

   pycmplot \
     --sum_stats hb.tsv,mcv.tsv \
     --labels Hb,MCV \
     --logp --signif_line --highlight \
     --annotate GENE \
     --trim_pval 0.01 \
     --cache --cache_dir ./.pycmplot \
     --output_dir ./out

Cache flags (all mirror the Python API):

- ``--cache`` — enable caching + hits overlay.
- ``--cache_dir PATH`` — where to store artefacts (default ``./.pycmplot``).
- ``--no_resume`` — regenerate everything but still write fresh entries.
- ``--clear_cache`` — delete the cache tree and exit.
- ``-qq / --qq_plot`` — implies ``compute_pvals=True`` under the hood.
- ``-V / --version`` — print version and exit.

See :ref:`cli` for the full reference.


.. _tut-faq:

Troubleshooting
---------------

**QQ plot fails with** ``bundle['pvals']`` **being** ``None``
   You forgot to pass ``compute_pvals=True`` to the loader (see
   :ref:`tut-qq`).  The default flipped from ``True`` to ``False``
   in 0.4.0 to avoid the ~80 MB-per-track memory cost when you're
   not making QQ plots.

**Cache says HIT for track A but MISS for track B on a warm re-run**
   Almost always means Stage-1 parameters differ between the two
   tracks (e.g. ``signif_threshold`` was auto-computed for one but
   explicitly set for the other).  Explicitly pass every Stage-1
   parameter you care about to lock in reproducibility.

**Hex colours in my TSV are being truncated / dropped**
   Don't use ``pd.read_csv(comment='#')`` to inspect the overlay —
   it treats ``#`` mid-cell as a comment marker and eats hex-color
   values.  The pycmplot reader strips only leading ``#`` header
   lines; if you inspect the file yourself, strip leading ``#`` lines
   manually and then read as a plain TSV.

**Highlight legend didn't appear**
   :func:`~pycmplot.annotation.build_highlight_legend_entries` returns
   an empty list when everything is at defaults
   (``category=significant`` for all rows AND ``highlight_color=auto``
   for all rows) — the plotter then skips the legend to preserve the
   pre-feature layout.  Edit at least one row's ``category`` or
   ``highlight_color`` to trigger it.

**"Category X appears with multiple colours" warning**
   You set two loci to the same ``category`` but gave them different
   ``highlight_color`` values.  The first colour wins in the legend;
   the warning names the value that lost so you can fix the TSV.

**Multi-panel run overwrote another panel's hits overlay**
   This shouldn't happen from 0.4.0 onward — each ``(files,
   parameters)`` group gets its own ``hits.<group_key>.tsv``.  If
   you see it, check that you're really passing distinct
   ``sum_stats`` lists to each panel's loader call.

**Hits overlay TSV got corrupted (crash mid-save, truncated line, etc.)**
   The loader logs ``"Hits overlay unreadable (…); ignoring."`` and
   regenerates the hits table from the cached per-track leads
   (``<label>.<key>.leads.parquet``) — no raw sumstats reload, no
   per-track cache invalidation.  The plot renders fine.  **User
   edits in the corrupted file are not recovered**, though: any
   ``source=user`` rows you added, ``highlight_color`` overrides,
   and ``category`` labels need to be re-applied in the freshly
   written TSV.  If you make heavy manual edits, keep a copy of the
   overlay TSV under version control alongside your analysis
   scripts.


Next steps
----------

- :ref:`cli` — full CLI reference.
- :ref:`api` — complete Python API reference.
- :ref:`python_api_notebook` — an executable Jupyter walkthrough.
- :doc:`changelog` — what changed in each release.
