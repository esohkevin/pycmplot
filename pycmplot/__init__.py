"""
pycmplot
========
Multi-track circular and linear Manhattan plot generation for GWAS summary statistics.

Quickstart
----------
Command-line::

    pycmplot -s file1.gz,file2.gz -l HbF,MCV --logp --mode lm

Python API::

    from pycmplot.plotting import multi_track_linear_manhattan, plot_circosm
    from pycmplot.stats import get_lead_snps
    from pycmplot.annotation import get_hits_summary_table

Public surface
--------------
"""

from pycmplot.plotting.linear import multi_track_linear_manhattan
from pycmplot.plotting.circular import plot_circosm, compute_track_radii_dict
from pycmplot.stats import get_lead_snps, get_highlight_snps
from pycmplot.io import get_sumstats_and_merged_sector_list
from pycmplot.annotation import get_hits_summary_table
from pycmplot.constants import hg38_chr_lengths, BIOTYPE_WEIGHTS
from pycmplot.resources import ResourceConfig

__all__ = [
    "multi_track_linear_manhattan",
    "plot_circosm",
    "compute_track_radii_dict",
    "get_lead_snps",
    "get_highlight_snps",
    "get_sumstats_and_merged_sector_list",
    "get_hits_summary_table",
    "hg38_chr_lengths",
    "BIOTYPE_WEIGHTS",
    "ResourceConfig",
]

__version__ = "0.1.0"
