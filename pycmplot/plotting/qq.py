from __future__ import annotations

"""
pycmplot.plotting.qq
====================
QQ (quantile-quantile) plots for GWAS p-values.

Speed notes
-----------
GWAS summary statistics often contain millions of SNPs.  Most of those points
lie near the null diagonal and are visually redundant.  Two optimisations are
applied by default:

1. **P-value thinning** (``thin_below`` / ``max_points``):
   All points above a -log10(p) tail threshold are kept in full; the bulk
   of null-like points below that threshold are randomly downsampled to at
   most ``max_points`` total.  Lambda (λ) is always computed on the *full*
   unfiltered array before thinning, so the statistic is never affected.

2. **Rasterised scatter** (``rasterized=True``):
   The scatter layer is rendered as a bitmap inside vector formats (PDF/SVG),
   dramatically reducing file size and save time for large point clouds.

Public functions
----------------
thin_pvals          Downsample null-like p-values for fast plotting.
qq_single      Draw one QQ plot onto a given Axes.
plot_qq_combined    All QQ plots in a single figure (grid layout).
plot_qq_separate    One output file per sumstat.
plot_qq_overlay     All sumstats overlaid on one axes, coloured by label.
"""

import logging
import math
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from scipy.stats import beta as beta_dist
from pycmplot.io import get_output_paths

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

_COMPUTE_PVALS_HINT = (
    "This usually means the loader did not materialise per-track p-value "
    "arrays.  As of pycmplot 0.4.0 the default of "
    "``compute_pvals`` on ``get_sumstats_and_merged_sector_list`` is "
    "``False``; Python-API callers that plan to draw a QQ figure must "
    "pass ``compute_pvals=True`` explicitly.  On the CLI, ``--qq`` "
    "sets it automatically."
)


def _reject_missing_pvals(pvals, label: Optional[str] = None) -> None:
    """Raise a clear ``ValueError`` when a QQ input is missing.

    Used at the entry point of every public ``plot_qq_*`` function so a
    caller who forgot ``compute_pvals=True`` on the loader (or who
    supplied a ``pval_dict`` with ``None`` values) gets a directive
    error message instead of a downstream ``TypeError`` when the array
    is later dereferenced.
    """
    if pvals is None:
        who = f" for track {label!r}" if label else ""
        raise ValueError(
            f"QQ plot requested but p-values are None{who}.  "
            + _COMPUTE_PVALS_HINT
        )


def _to_scalar_float(value, name: str = "value") -> float:
    """Coerce a value to a Python ``float`` with a clear error path.

    Accepts:
      * plain ``int`` / ``float``
      * NumPy 0-d arrays and 1-element 1-d arrays / lists / tuples /
        pandas Series (unwrapped to their single element)

    Raises :class:`TypeError` naming *value* and its offending type
    when it can't reasonably be coerced -- much more actionable than
    matplotlib's downstream
    ``"float() argument must be a real number, not a 'list'"``.
    """
    if isinstance(value, (int, float)):
        return float(value)
    if hasattr(value, "item"):
        try:
            return float(value.item())
        except (ValueError, TypeError):
            pass
    if hasattr(value, "__len__"):
        try:
            if len(value) == 1:
                return float(next(iter(value)))
        except (ValueError, TypeError):
            pass
    raise TypeError(
        f"{name!r} must be a scalar number; received {type(value).__name__} "
        f"with value {value!r}.  If this came from a QQ plotter, check that "
        f"you passed a scalar for {name} (not a list / ndarray / Series)."
    )


def _validate_pval_dict(pval_dict) -> None:
    """Validate a ``pval_dict`` input to a multi-track QQ plotter."""
    if pval_dict is None:
        raise ValueError(
            "QQ plot requested but pval_dict is None.  "
            + _COMPUTE_PVALS_HINT
        )
    if not len(pval_dict):
        raise ValueError("pval_dict is empty.")
    missing = [k for k, v in pval_dict.items() if v is None]
    if missing:
        raise ValueError(
            f"QQ plot requested but p-values are None for tracks: "
            f"{missing!r}.  " + _COMPUTE_PVALS_HINT
        )


# ---------------------------------------------------------------------------
# Thinning helper
# ---------------------------------------------------------------------------

def thin_pvals(
    pvals: np.ndarray,
    tail_threshold: float = 0.01,
    max_points: int = 50_000,
    seed: int = 42,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Downsample p-values for faster QQ plotting with no visible breaks.
 
    Rather than splitting into tail / bulk regions with different sampling
    strategies (which produces a visible seam at the threshold), this function
    uses a single **log-uniform** thinning pass over all p-values:
 
    1. Sort p-values ascending and convert to −log₁₀ scale.
    2. Pick ``max_points`` evenly-spaced indices along the −log₁₀ axis.
       Because −log₁₀ compresses large p-values and expands small ones, this
       automatically gives dense coverage in the interesting tail and sparse
       coverage in the null bulk — with no hard boundary.
 
    Parameters
    ----------
    pvals:
        Full array of raw p-values.
    tail_threshold:
        Kept for API compatibility; no longer used as a hard split point.
        All points above −log₁₀(tail_threshold) are always represented because
        the log-uniform spacing naturally keeps them.
    max_points:
        Maximum number of points to return (default 50 000).
    seed:
        Unused (kept for API compatibility — log-uniform selection is
        deterministic).
 
    Returns
    -------
    (kept_pvals, kept_ranks, n_full)
        *kept_pvals*  — thinned p-values in ascending order.
        *kept_ranks*  — 1-based ranks in the full sorted array.
        *n_full*      — total SNP count before thinning (for expected quantiles).
 
    Notes
    -----
    Lambda (λ) must be computed on the full *pvals* array **before** calling
    this function — thinning changes the empirical distribution.
    """
    pvals = np.asarray(pvals, dtype=float)
    pvals = pvals[np.isfinite(pvals) & (pvals > 0) & (pvals <= 1)]
    n_full = len(pvals)
 
    if n_full <= max_points:
        # Nothing to thin
        sort_idx = np.argsort(pvals)
        return pvals[sort_idx], np.arange(1, n_full + 1), n_full
 
    # Sort ascending; full_ranks[i] = i+1
    pvals_sorted = np.sort(pvals)
    full_ranks   = np.arange(1, n_full + 1)
 
    # Work in −log10 space so spacing is proportional to visual separation
    logp = -np.log10(pvals_sorted)          # ascending p → descending logp
    logp_min = logp[0]                      # smallest logp (bulk end)
    logp_max = logp[-1]                     # largest logp (tail end)
 
    # Evenly-spaced target positions along the logp axis
    targets = np.linspace(logp_min, logp_max, max_points)
 
    # For each target, pick the closest actual point (searchsorted on
    # the reversed array since logp is descending)
    logp_desc = logp[::-1]                  # descending for searchsorted
    idx_desc  = np.searchsorted(logp_desc, targets, side="left")
    idx_desc  = np.clip(idx_desc, 0, n_full - 1)
 
    # Convert back to ascending-p indices and deduplicate
    idx_asc = (n_full - 1 - idx_desc)
    idx_asc = np.unique(idx_asc)            # sorted, no duplicates
 
    kept_pvals = pvals_sorted[idx_asc]
    kept_ranks = full_ranks[idx_asc]
 
    n_kept = len(kept_pvals)
    logger.debug(
        "QQ thinning: %d → %d points (%.1f%% retained)",
        n_full, n_kept, 100 * n_kept / n_full,
    )
 
    return kept_pvals, kept_ranks, n_full


# ---------------------------------------------------------------------------
# Core array builder
# ---------------------------------------------------------------------------

def _qq_arrays(
    pvals: np.ndarray,
    ranks: Optional[np.ndarray] = None,
    n_full: Optional[int] = None,
    ci: float = 0.95,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (expected, observed, ci_lower, ci_upper) in −log10 scale.

    Parameters
    ----------
    pvals:
        Sorted (ascending) p-values to plot — may be a thinned subset.
    ranks:
        1-based ranks of *pvals* in the full distribution.  If ``None``,
        assumes *pvals* is the complete set and ranks are 1..n.
    n_full:
        Total number of SNPs in the full (pre-thinning) dataset.  Used to
        compute correct expected quantiles.  Defaults to ``len(pvals)``.
    ci:
        Confidence interval level.
    """
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)

    if n_full is None:
        n_full = n
    if ranks is None:
        ranks = np.arange(1, n + 1)

    # Expected −log10(p): rank i → expected p = i/(n_full+1)
    expected = -np.log10(ranks / (n_full + 1))

    # Observed −log10(p): rank i paired with the i-th smallest p-value
    observed = -np.log10(pvals)

    # CI from the beta distribution (uses original ranks in full dataset)
    alpha = 1.0 - ci
    ci_lo = -np.log10(beta_dist.ppf(1 - alpha / 2, ranks, n_full - ranks + 1))
    ci_hi = -np.log10(beta_dist.ppf(    alpha / 2, ranks, n_full - ranks + 1))

    # Sort by expected ascending for clean polygon fill
    order = np.argsort(expected)
    return expected[order], observed[order], ci_lo[order], ci_hi[order]


# ---------------------------------------------------------------------------
# Lambda
# ---------------------------------------------------------------------------

def _compute_lambda(pvals: np.ndarray) -> float:
    """Genomic inflation factor λ = median(χ²_obs) / median(χ²_expected).

    Always returns a plain Python ``float`` (never a NumPy scalar or
    0-d array) so downstream f-string formatting like ``f"{lam:.4f}"``
    can't accidentally raise the opaque ``float()``-argument TypeError
    when a caller passed an unusual pvals shape.
    """
    from scipy.stats import chi2
    pvals = np.asarray(pvals, dtype=float).ravel()
    pvals = pvals[np.isfinite(pvals) & (pvals > 0) & (pvals <= 1)]
    if len(pvals) == 0:
        return float("nan")
    obs_median_chi2 = chi2.ppf(1 - float(np.median(pvals)), df=1)
    expected_median_chi2 = chi2.ppf(0.5, df=1)   # ≈ 0.4549
    return round(float(obs_median_chi2 / expected_median_chi2), 4)


# ---------------------------------------------------------------------------
# Get Legend Best Location
# ---------------------------------------------------------------------------

def get_legend_loc_string(legend, ax):
    from matplotlib.legend import Legend
    import random
    """
    Returns the string name of the location selected by loc='best'.
    """
    # 1. Force a layout update to compute the final geometry
    fig = ax.get_figure()
    fig.canvas.draw()
    
    # 2. Get bounding boxes in pixels
    leg_bbox = legend.get_window_extent()
    ax_bbox = ax.bbox
    
    # 3. Calculate normalized center position of the legend inside the axes (0 to 1)
    leg_center_x = (leg_bbox.x0 + leg_bbox.x1) / 2
    leg_center_y = (leg_bbox.y0 + leg_bbox.y1) / 2
    
    norm_x = (leg_center_x - ax_bbox.x0) / ax_bbox.width
    norm_y = (leg_center_y - ax_bbox.y0) / ax_bbox.height
    
    # 4. Map the normalized positions to text positions
    if norm_y > 0.66:
        y_str = "upper"
    elif norm_y < 0.33:
        y_str = "lower"
    else:
        y_str = "center"
        
    if norm_x > 0.66:
        x_str = "right"
    elif norm_x < 0.33:
        x_str = "left"
    else:
        x_str = "center"

    # Standardize output for center alignments
    if y_str == "center" and x_str == "center":
        best_loc = "center"
    elif x_str == "center":
        best_loc = f"{y_str} center"
    elif y_str == "center":
        best_loc = x_str # matches 'right' or 'left'
    else:
        best_loc = f"{y_str} {x_str}"

    # 5. Get a list of all legend locations and exclude it and other undesirable locations
    # best_loc, 'center', 'upper center', 'upper left'
    # Genomic control inflation factor (lambda) already occupies 'upper left'
    remove_locs = [best_loc, 'center', 'upper center', 'upper left']
    leg_locs = list(Legend.codes.keys())
    filtered_locs = list(filter(lambda x: x not in remove_locs, leg_locs))
    leg_loc = random.choice(filtered_locs)

    return leg_loc

# ---------------------------------------------------------------------------
# Get Legend Best Location
# ---------------------------------------------------------------------------

def set_legend_loc_string():
    from matplotlib.legend import Legend
    import random
    """
    Returns the string name of the location selected by loc='best'.
    """
    # Get a list of all legend locations and exclude it and other undesirable locations
    # 'center', 'upper center', 'upper left'
    # Genomic control inflation factor (lambda) already occupies 'upper left'
    remove_locs = ['center', 'upper center', 'upper left']
    leg_locs = list(Legend.codes.keys())
    filtered_locs = list(filter(lambda x: x not in remove_locs, leg_locs))
    leg_loc = random.choice(filtered_locs)

    return leg_loc

# ---------------------------------------------------------------------------
# Single-axis QQ plot
# ---------------------------------------------------------------------------

def qq_single(
    pvals: np.ndarray | pd.Series,
    ax: plt.Axes,
    label: Optional[str] = None,
    color: str = "steelblue",
    point_size: float = 8,
    ci: float = 0.95,
    ci_alpha: float = 0.15,
    signif_threshold: Optional[float] = 5e-8,
    show_lambda: bool = True,
    title: Optional[str] = None,
    # --- speed options ---
    thin: bool = True,
    thin_below: float = 0.01,
    max_points: int = 50_000,
    fontsize: float = 8,
    rasterized: bool = True,
) -> plt.Axes:
    """Draw a single QQ plot onto *ax*.

    Parameters
    ----------
    pvals:
        Array or Series of raw p-values (not −log10).
    ax:
        Matplotlib Axes to draw on.
    label:
        Legend label for the scatter points.
    color:
        Colour for points and CI fill.
    point_size:
        Scatter point size.
    ci:
        Confidence interval level (default 0.95).
    ci_alpha:
        Transparency of the CI band.
    signif_threshold:
        If given, draw a horizontal dashed line at −log10(threshold).
    show_lambda:
        Annotate the plot with the genomic inflation factor λ.
    title:
        Axes title.
    thin:
        Enable p-value thinning for speed (default ``True``).
    thin_below:
        P-value threshold below which all points are always kept.
        Points above this threshold are downsampled.
    max_points:
        Maximum number of points to plot after thinning (default 50 000).
    rasterized:
        Render the scatter as a bitmap inside vector output formats —
        greatly reduces PDF/SVG file size (default ``True``).

    Returns
    -------
    plt.Axes
    """

    # Guard against the common mistake of passing a numpy array of Axes
    # (e.g. from plt.subplots(1, 2)) instead of a single Axes object.
    if not hasattr(ax, "fill_between"):
        raise TypeError(
            "'ax' must be a single Matplotlib Axes object, but received "
            f"{type(ax).__name__}.\n"
            "If you created the figure with plt.subplots(nrows, ncols), "
            "index the returned array, e.g.:\n"
            "  fig, axes = plt.subplots(1, 2)\n"
            "  qq_single(pvals, ax=axes[0])"
        )

    # Directive error message when pvals is None — usually a symptom of
    # the loader default ``compute_pvals=False`` (see _reject_missing_pvals).
    _reject_missing_pvals(pvals, label=label)

    # Coerce numeric kwargs to plain floats once, at the top.  Every
    # ax.<something>(..., fontsize=fontsize, point_size=...) call below
    # is now safe against callers who accidentally hand in a
    # single-element list / 0-d ndarray / pandas Series — matplotlib
    # otherwise raises the opaque
    # "TypeError: float() argument must be a real number, not a 'list'"
    # from deep inside its own numeric parsing.
    fontsize = _to_scalar_float(fontsize, name="fontsize")
    point_size = _to_scalar_float(point_size, name="point_size")

    pvals_full = np.asarray(pvals, dtype=float)
    pvals_full = pvals_full[np.isfinite(pvals_full) & (pvals_full > 0) & (pvals_full <= 1)]

    # Lambda always on the full array
    lam = _compute_lambda(pvals_full)

    if thin and len(pvals_full) > max_points:
        plot_pvals, plot_ranks, n_full = thin_pvals(
            pvals_full, tail_threshold=thin_below, max_points=max_points
        )
    else:
        plot_pvals = np.sort(pvals_full)
        plot_ranks = np.arange(1, len(plot_pvals) + 1)
        n_full = len(plot_pvals)

    expected, observed, ci_lo, ci_hi = _qq_arrays(
        plot_pvals, ranks=plot_ranks, n_full=n_full, ci=ci
    )

    # CI band -- deliberately unlabelled to keep the legend uncluttered.
    # The shaded region hugging the diagonal is self-evidently a
    # confidence band; adding a "95% CI" legend entry only crowds the
    # top corner alongside the track label and the ``lambda``
    # annotation.  The CI level is captured in the docstring / caller;
    # if downstream users want a legend entry back, they can supply
    # ``label=`` on their own overlay.
    ax.fill_between(
        expected, ci_lo, ci_hi,
        color=color, alpha=ci_alpha, linewidth=0,
    )

    # Diagonal null line
    max_val = max(expected.max(), observed.max()) * 1.05
    ax.plot([0, max_val], [0, max_val], color="grey", linewidth=0.8,
            linestyle="--", zorder=1)

    # Observed points — use ``ax.plot`` with marker-only style rather than
    # ``ax.scatter``.  At 1 M+ points this swaps a giant ``PathCollection``
    # (one path + ``should_simplify`` check per point) for a single
    # ``Line2D`` whose marker draw loop is far cheaper, while producing
    # visually identical rasterised output.  ``point_size`` is in scatter
    # units (points²); convert to plot's ``markersize`` (points) via
    # ``sqrt`` so existing user-supplied values keep their meaning.
    _ms = float(np.sqrt(point_size)) if point_size else 2.0
    ax.plot(
        expected, observed,
        marker="o", linestyle="none",
        color=color, alpha=0.85,
        markersize=_ms, markeredgecolor="none",
        label=label, zorder=2,
        rasterized=rasterized,
    )

    '''"""
    # Significance line
    if signif_threshold is not None:
        sig_logp = -np.log10(signif_threshold)
        ax.axhline(sig_logp, color="red", linewidth=0.7, linestyle="--",
                   label=f"p={signif_threshold:.0e}")
    """'''

    # Lambda annotation.  Defensive scalar coercion for both ``lam`` and
    # ``fontsize`` -- when a caller (or a wrapper that forwarded
    # kwargs) accidentally hands in a list / 1-d ndarray, the f-string
    # ``.4f`` spec on ``lam`` and matplotlib's numeric parsing on
    # ``fontsize`` both raise the opaque
    # ``TypeError: float() argument must be a real number, not a 'list'``
    # at this call site.  Coerce here so the error surfaces earlier with
    # the offending argument named, and so the common accident
    # (accidental list-wrapping) just Works.
    _lam = _to_scalar_float(lam, name="lam")
    _fs = _to_scalar_float(fontsize, name="fontsize")
    if show_lambda and not math.isnan(_lam):
        ax.text(
            0.05, 0.95,
            f"λ = {_lam:.4f}",
            transform=ax.transAxes,
            va="top", ha="left",
            fontsize=_fs, fontstyle="italic",
            color="black",
        )

    ax.set_xlabel("Expected −log₁₀(p)", fontsize=10)
    ax.set_ylabel("Observed −log₁₀(p)", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if title:
        ax.set_title(title, fontsize=10, pad=6)
    else:
        if label:
            ax.set_title(label, fontsize=10, pad=6)
            #leg_loc = set_legend_loc_string()
            #ax.legend(fontsize=fontsize, frameon=False, loc=leg_loc)
            #ax.legend(fontsize=fontsize, frameon=False, loc="best")

    ax.set_xlim(0, max(expected)+1)
    ax.set_ylim(0, max(observed)+2)

    return ax


# ---------------------------------------------------------------------------
# Combined multi-panel figure
# ---------------------------------------------------------------------------

def qq_combined(
    pval_dict: dict[str, np.ndarray | pd.Series],
    colors: Optional[list[str]] = None,
    point_size: float = 8,
    ci: float = 0.95,
    signif_threshold: Optional[float] = 5e-8,
    show_lambda: bool = True,
    ncols: int = 3,
    figsize: Optional[tuple] = None,
    dpi: int = 300,
    title: Optional[str] = None,
    output_path: Optional[str] = None,
    fig_format: str = "png",
    thin: bool = True,
    thin_below: float = 0.01,
    max_points: int = 50_000,
    fontsize: float = 8,
    rasterized: bool = True,
) -> tuple[plt.Figure, list[plt.Axes]]:
    """Plot all QQ plots in a single figure arranged in a grid.

    Parameters
    ----------
    pval_dict:
        Ordered dict of ``{label: p_value_array}``.
    colors:
        List of colours, one per track.  Cycles if fewer than tracks.
    ncols:
        Number of columns in the subplot grid (default 3).
    figsize:
        Figure size.  Auto-calculated from *ncols* and number of tracks
        if ``None``.
    output_path:
        If given, save the figure here.
    thin, thin_below, max_points, rasterized:
        See :func:`qq_single`.

    Returns
    -------
    (fig, axes)
    """
    _validate_pval_dict(pval_dict)
    n = len(pval_dict)

    nrows = math.ceil(n / ncols)

    cmap = plt.get_cmap("tab10")
    colors = [mcolors.to_hex(cmap(i % 10)) for i in range(n)]
    #if colors is None:
    #    cmap = plt.get_cmap("tab10")
    #    colors = [mcolors.to_hex(cmap(i % 10)) for i in range(n)]
    #elif len(colors) < n:
    #    colors = [colors[i % len(colors)] for i in range(n)]

    if figsize is None:
        figsize = (ncols * 4.5, nrows * 4.5)

    fig, axes_grid = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False, layout="constrained")
    axes_flat = axes_grid.flatten()

    for idx, (label, pvals) in enumerate(pval_dict.items()):
        qq_single(
            pvals=pvals,
            ax=axes_flat[idx],
            label=label,
            color=colors[idx],
            point_size=point_size,
            ci=ci,
            signif_threshold=signif_threshold,
            show_lambda=show_lambda,
            title=label,
            thin=thin,
            fontsize=fontsize,
            thin_below=thin_below,
            max_points=max_points,
            rasterized=rasterized,
        )

    for ax in axes_flat[n:]:
        ax.set_visible(False)

    if title:
        # Compute the left and right bounds of the active subplot grid
        bboxes = [ax.get_position() for ax in axes_flat[:n] if ax.get_visible()]
        left = min(b.x0 for b in bboxes)
        right = max(b.x1 for b in bboxes)
        center_x = (left + right) / 2.0

        fig.suptitle(title, fontsize=13, x=center_x, y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    #if title:
    #    fig.suptitle(title, fontsize=13)
    #plt.tight_layout()

    if output_path:
        fmt = fig_format or Path(output_path).suffix.lstrip(".") or "png"
        fig.savefig(f"{output_path}.{fmt}", format=fmt, dpi=dpi, bbox_inches="tight")
        logger.info("Saved combined QQ plot: %s", f"{output_path}.{fmt}")

    return fig, list(axes_flat[:n])


# ---------------------------------------------------------------------------
# Separate figures — one file per sumstat
# ---------------------------------------------------------------------------

def qq_separate(
    pval_dict: dict[str, np.ndarray | pd.Series],
    base_name: str = None,
    output_path: str = ".",
    colors: Optional[list[str]] = None,
    point_size: float = 8,
    ci: float = 0.95,
    signif_threshold: Optional[float] = 5e-8,
    show_lambda: bool = True,
    figsize: tuple = (5, 5),
    dpi: int = 300,
    fig_format: str = "png",
    thin: bool = True,
    thin_below: float = 0.01,
    max_points: int = 50_000,
    fontsize: float = 8,
    rasterized: bool = True,
) -> list[str]:
    """Save one QQ plot per sumstat as individual files.

    Parameters
    ----------
    pval_dict:
        Ordered dict of ``{label: p_value_array}``.
    output_dir:
        Directory to save files in.
    file_stem:
        Prefix for output filenames.
    colors:
        List of colours, one per track.
    thin, thin_below, max_points, rasterized:
        See :func:`qq_single`.

    Returns
    -------
    List of output file paths.
    """

    #labels = pval_dict.keys()

    ## plot name
    #(
    #    plt_name, 
    #    table_out,
    #    plt_base,
    #) = get_output_paths(
    #    labels = labels,
    #    mode='qq',
    #    output_dir=output_path, 
    #    plot_title=base_name, 
    #    output_format=fig_format
    #)

    _validate_pval_dict(pval_dict)
    n = len(pval_dict)

    cmap = plt.get_cmap("tab10")
    colors = [mcolors.to_hex(cmap(i % 10)) for i in range(n)]
    #if colors is None:
    #    cmap = plt.get_cmap("tab10")
    #    colors = [mcolors.to_hex(cmap(i % 10)) for i in range(n)]
    #elif len(colors) < n:
    #    colors = [colors[i % len(colors)] for i in range(n)]

    saved: list[str] = []

    for idx, (label, pvals) in enumerate(pval_dict.items()):
        fig, ax = plt.subplots(figsize=figsize)

        qq_single(
            pvals=pvals,
            ax=ax,
            label=label,
            color=colors[idx],
            point_size=point_size,
            ci=ci,
            signif_threshold=signif_threshold,
            show_lambda=show_lambda,
            title=label,
            thin=thin,
            fontsize=fontsize,
            thin_below=thin_below,
            max_points=max_points,
            rasterized=rasterized,
        )

        plt.tight_layout()

        safe_label = label.replace(" ", "_").replace("/", "-")
        out_path = f"{output_path}_{safe_label}.{fig_format}"
        fig.savefig(out_path, format=fig_format, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        logger.info("Saved QQ plot: %s", out_path)
        saved.append(out_path)

    return saved


# ---------------------------------------------------------------------------
# Overlay — all sumstats on one axes
# ---------------------------------------------------------------------------

def qq_overlay(
    pval_dict: dict[str, np.ndarray | pd.Series],
    colors: Optional[list[str]] = None,
    point_size: float = 8,
    ci: float = 0.95,
    ci_alpha: float = 0.10,
    signif_threshold: Optional[float] = 5e-8,
    show_lambda: bool = True,
    figsize: tuple = (6, 6),
    dpi: int = 300,
    title: Optional[str] = None,
    output_path: Optional[str] = None,
    fig_format: str = "png",
    thin: bool = True,
    thin_below: float = 0.01,
    max_points: int = 50_000,
    fontsize: float = 8,
    rasterized: bool = True,
) -> tuple[plt.Figure, plt.Axes]:
    """Plot all sumstats on a single QQ axes, each coloured differently.

    Lambda (λ) values appear in the legend label for each sumstat.

    Parameters
    ----------
    pval_dict:
        Ordered dict of ``{label: p_value_array}``.
    colors:
        List of colours, one per sumstat.  Defaults to ``tab10`` palette.
    ci_alpha:
        Transparency of CI bands (default 0.10 — lower than single-panel
        default to keep overlapping bands readable).
    show_lambda:
        Append λ to each legend entry.
    thin, thin_below, max_points, rasterized:
        See :func:`qq_single`.

    Returns
    -------
    (fig, ax)
    """

    #labels = pval_dict.keys()

    ## plot name
    #(
    #    plt_name, 
    #    table_out,
    #    plt_base,
    #) = get_output_paths(
    #    labels = labels,
    #    mode='qq',
    #    output_dir=output_path,
    #    plot_title=title, 
    #    output_format=fig_format
    #)

    _validate_pval_dict(pval_dict)
    n = len(pval_dict)

    # Coerce numeric kwargs once (see qq_single for rationale).
    point_size = _to_scalar_float(point_size, name="point_size")

    cmap = plt.get_cmap("tab10")
    colors = [mcolors.to_hex(cmap(i % 10)) for i in range(n)]
    #if colors is None:
    #    cmap = plt.get_cmap("tab10")
    #    colors = [mcolors.to_hex(cmap(i % 10)) for i in range(n)]
    #elif len(colors) < n:
    #    colors = [colors[i % len(colors)] for i in range(n)]

    fig, ax = plt.subplots(figsize=figsize)
    global_max = 0.0
    max_expected_global = 0.0
    max_observed_global = 0.0

    for idx, (label, pvals) in enumerate(pval_dict.items()):
        pvals_full = np.asarray(pvals, dtype=float)
        pvals_full = pvals_full[np.isfinite(pvals_full) & (pvals_full > 0) & (pvals_full <= 1)]

        # Lambda on full array before any thinning
        lam = _compute_lambda(pvals_full)

        if thin and len(pvals_full) > max_points:
            plot_pvals, plot_ranks, n_full = thin_pvals(
                pvals_full, tail_threshold=thin_below, max_points=max_points
            )
        else:
            plot_pvals = np.sort(pvals_full)
            plot_ranks = np.arange(1, len(plot_pvals) + 1)
            n_full = len(plot_pvals)

        expected, observed, ci_lo, ci_hi = _qq_arrays(
            plot_pvals, ranks=plot_ranks, n_full=n_full, ci=ci
        )

        color = colors[idx]
        _lam_i = _to_scalar_float(lam, name="lam")
        legend_label = f"{label}  (λ={_lam_i:.4f})" if show_lambda else label

        ax.fill_between(
            expected, ci_lo, ci_hi,
            color=color, alpha=ci_alpha, linewidth=0,
        )
        # Use ax.plot (Line2D) instead of ax.scatter (PathCollection) for
        # the same draw-time reasons as in qq_single — much faster at
        # large N with identical rasterised output.
        _ms = float(np.sqrt(point_size)) if point_size else 2.0
        ax.plot(
            expected, observed,
            marker="o", linestyle="none",
            color=color, alpha=0.85,
            markersize=_ms, markeredgecolor="none",
            label=legend_label, zorder=2 + idx,
            rasterized=rasterized,
        )

        global_max = max(global_max, expected.max(), observed.max())
        max_expected_global = max(max_expected_global, float(expected.max()))
        max_observed_global = max(max_observed_global, float(observed.max()))

    ax.plot(
        [0, global_max * 1.05], [0, global_max * 1.05],
        color="grey", linewidth=0.8, linestyle="--", zorder=1,
    )

    '''"""
    if signif_threshold is not None:
        ax.axhline(
            -np.log10(signif_threshold),
            color="red", linewidth=0.7, linestyle="--",
            label=f"p = {signif_threshold:.0e}",
        )
    """'''

    ax.set_xlabel("Expected −log₁₀(p)", fontsize=11)
    ax.set_ylabel("Observed −log₁₀(p)", fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.legend(
        fontsize=8, 
        frameon=False, 
        framealpha=0.7, 
        edgecolor="lightgrey", 
        loc="upper left",
    )

    if title:
        ax.set_title(title, fontsize=11, pad=8)

    #plt.xlim(0, max(expected)+2)
    #plt.ylim(0, max(observed)+1)
    ax.set_xlim(0, max_expected_global + 1)
    ax.set_ylim(0, max_observed_global + 2)

    plt.tight_layout()

    if output_path:
        fmt = fig_format or Path(output_path).suffix.lstrip(".") or "png"
        fig.savefig(f"{output_path}.{fmt}", format=fmt, dpi=dpi, bbox_inches="tight")
        logger.info("Saved overlay QQ plot: %s", f"{output_path}.{fmt}")

    return fig, ax


# ---------------------------------------------------------------------------
# Backwards-compatible aliases (deprecated in 0.4.1)
# ---------------------------------------------------------------------------
# The three QQ entry points shed the ``plot_`` prefix so the API reads as
# ``pcm.qq_combined(...)``, ``pcm.qq_overlay(...)``, ``pcm.qq_separate(...)``.
# Old names remain importable for one release cycle with a soft warning.
from pycmplot._deprecation import _deprecated_alias as _da
plot_qq_single = _da(
    qq_single, old_name="plot_qq_single", new_name="qq_single",
)
plot_qq_combined = _da(
    qq_combined, old_name="plot_qq_combined", new_name="qq_combined",
)
plot_qq_overlay = _da(
    qq_overlay, old_name="plot_qq_overlay", new_name="qq_overlay",
)
plot_qq_separate = _da(
    qq_separate, old_name="plot_qq_separate", new_name="qq_separate",
)