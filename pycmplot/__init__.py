"""
pycmplot
========
Multi-track circular and linear Manhattan plot generation for GWAS summary statistics.

Quickstart
----------
Command-line::

    pycmplot -s file1.gz,file2.gz -l HbF,MCV --logp --mode lm

Python API (0.4.x — short names)::

    import pycmplot as pcm

    fi     = pcm.prep(sum_stats=files, labels=labels)
    bundle = pcm.load(sum_stats=files, labels=labels, file_info=fi,
                      logp=True, trim_pval=0.01, highlight=True)
    pcm.linear(sumstats_loaded=bundle["dfs"], hits_table=bundle["annot"],
               signif_lines=bundle["lines"], output_dir="./out")
    pcm.circular(sumstats_loaded=bundle["dfs"],
                 sector_sizes=bundle["sectors"], hits_table=bundle["annot"])

Backwards-compat (old names still work with a soft ``DeprecationWarning``)::

    from pycmplot import (
        prep_pycmplot_input_info,          # -> prep
        get_sumstats_and_merged_sector_list,  # -> load
        plot_linear, plot_circular,         # -> linear, circular
        plot_qq_combined, plot_qq_overlay, plot_qq_separate,  # -> qq_*
    )

Public surface
--------------
"""

# ---------------------------------------------------------------------------
# New short names (0.4.x)
# ---------------------------------------------------------------------------
from pycmplot.io import prep, load
from pycmplot.plotting.linear import linear
from pycmplot.plotting.circular import circular, compute_track_radii_dict
from pycmplot.plotting.qq import (
    qq_combined, qq_overlay, qq_separate, plot_qq_single,
)
from pycmplot.stats import get_lead_snps, get_highlight_snps
from pycmplot.annotation import get_hits_summary_table
from pycmplot.constants import hg38_chr_lengths, BIOTYPE_WEIGHTS
from pycmplot.resources import ResourceConfig

# ---------------------------------------------------------------------------
# Backwards-compat aliases (deprecated in 0.4.x — will be removed in a
# future release).  Each of these names is the same object exposed by
# its home module's deprecation wrapper, so importing either the module
# path or the top-level attribute both emit the warning uniformly.
# ---------------------------------------------------------------------------
from pycmplot.io import (
    prep_pycmplot_input_info,
    get_sumstats_and_merged_sector_list,
)
from pycmplot.plotting.linear import plot_linear
from pycmplot.plotting.circular import plot_circular
from pycmplot.plotting.qq import (
    plot_qq_combined, plot_qq_overlay, plot_qq_separate, plot_qq_single,
)

__all__ = [
    # New short names
    "prep",
    "load",
    "linear",
    "circular",
    "qq_single",
    "qq_combined",
    "qq_overlay",
    "qq_separate",
    "plot_qq_single",
    "compute_track_radii_dict",
    "get_lead_snps",
    "get_highlight_snps",
    "get_hits_summary_table",
    "hg38_chr_lengths",
    "BIOTYPE_WEIGHTS",
    "ResourceConfig",
    # Deprecated aliases (kept for one release cycle)
    "prep_pycmplot_input_info",
    "get_sumstats_and_merged_sector_list",
    "plot_linear",
    "plot_circular",
    "plot_qq_single",
    "plot_qq_combined",
    "plot_qq_overlay",
    "plot_qq_separate",
]

__version__ = "0.4.2"
