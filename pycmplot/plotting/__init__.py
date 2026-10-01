"""
pycmplot.plotting
=================

Plotting sub-package for pycmplot.  Public entry points:

* :func:`linear` — multi-track stacked linear Manhattan plot.
* :func:`circular` — Circos-style circular Manhattan plot.
* :func:`qq_single` — one QQ plot drawn onto an Axes you supply (for
  custom figures and multi-panel layouts).
* :func:`qq_combined`, :func:`qq_separate`, :func:`qq_overlay` — QQ
  plots for several traits: a grid, one file per trait, or one shared
  panel.

The ``plot_*`` names (:func:`plot_linear`, :func:`plot_circular`,
:func:`plot_qq_single`, :func:`plot_qq_combined`,
:func:`plot_qq_separate`, :func:`plot_qq_overlay`) are deprecated
aliases kept for backwards compatibility; they emit a
``DeprecationWarning`` when called.
"""

from pycmplot.plotting.linear import linear, plot_linear
from pycmplot.plotting.circular import (
    circular, plot_circular, compute_track_radii_dict,
)
from pycmplot.plotting.qq import (
    qq_single,
    qq_combined,
    qq_separate,
    qq_overlay,
    plot_qq_single,
    plot_qq_combined,
    plot_qq_separate,
    plot_qq_overlay,
)

__all__ = [
    "linear",
    "circular",
    "compute_track_radii_dict",
    "qq_single",
    "qq_combined",
    "qq_separate",
    "qq_overlay",
    # Deprecated aliases
    "plot_linear",
    "plot_circular",
    "plot_qq_single",
    "plot_qq_combined",
    "plot_qq_separate",
    "plot_qq_overlay",
]
