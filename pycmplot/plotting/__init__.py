"""
pycmplot.plotting
=================

Plotting sub-package for pycmplot.  Exposes the three public entry points
used by the CLI and by downstream scripts:

* :func:`plot_linear` — multi-track stacked linear Manhattan plot.
* :func:`plot_circular` — Circos-style circular Manhattan plot.
* :func:`plot_qq_combined`, :func:`plot_qq_separate`,
  :func:`plot_qq_overlay`, :func:`plot_qq_single` — QQ plotting helpers.
"""

from pycmplot.plotting.linear import plot_linear
from pycmplot.plotting.circular import plot_circular, compute_track_radii_dict
from pycmplot.plotting.qq import (
    plot_qq_single,
    plot_qq_combined,
    plot_qq_separate,
    plot_qq_overlay,
)

__all__ = [
    "plot_linear",
    "plot_circular",
    "compute_track_radii_dict",
    "plot_qq_single",
    "plot_qq_combined",
    "plot_qq_separate",
    "plot_qq_overlay",
]
