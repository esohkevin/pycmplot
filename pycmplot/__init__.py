"""
pycmplot
========
Multi-track circular and linear Manhattan plot generation for GWAS summary statistics.

Quickstart
----------
Command-line::

    pycmplot -s file1.gz,file2.gz -l HbF,MCV --logp --mode lm

Python API::

    from pycmplot.io import prep_pycmplot_input_info, get_sumstats_and_merged_sector_list
    from pycmplot.plotting import plot_linear, plot_circular, plot_qq_single, plot_qq_separate, plot_qq_overlay, plot_qq_combined
    from pycmplot.stats import get_lead_snps
    from pycmplot.annotation import get_hits_summary_table

Public surface
--------------
"""

from pycmplot.plotting.linear import plot_linear
from pycmplot.plotting.circular import plot_circular, compute_track_radii_dict
from pycmplot.plotting.qq import plot_qq_single, plot_qq_separate, plot_qq_overlay, plot_qq_combined
from pycmplot.stats import get_lead_snps, get_highlight_snps
from pycmplot.io import prep_pycmplot_input_info, get_sumstats_and_merged_sector_list
from pycmplot.annotation import get_hits_summary_table
from pycmplot.constants import hg38_chr_lengths, BIOTYPE_WEIGHTS
from pycmplot.resources import ResourceConfig

__all__ = [
    "plot_linear",
    "plot_circular",
    "plot_qq_single",
    "plot_qq_separate",
    "plot_qq_overlay",
    "plot_qq_combined",
    "compute_track_radii_dict",
    "get_lead_snps",
    "get_highlight_snps",
    "prep_pycmplot_input_info",
    "get_sumstats_and_merged_sector_list",
    "get_hits_summary_table",
    "hg38_chr_lengths",
    "BIOTYPE_WEIGHTS",
    "ResourceConfig",
]

__version__ = "0.2.8"
