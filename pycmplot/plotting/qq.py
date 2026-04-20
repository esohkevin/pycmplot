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
plot_qq_single      Draw one QQ plot onto a given Axes.
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
    """Genomic inflation factor λ = median(χ²_obs) / median(χ²_expected)."""
    from scipy.stats import chi2
    pvals = pvals[np.isfinite(pvals) & (pvals > 0) & (pvals <= 1)]
    if len(pvals) == 0:
        return float("nan")
    obs_median_chi2 = chi2.ppf(1 - np.median(pvals), df=1)
    expected_median_chi2 = chi2.ppf(0.5, df=1)   # ≈ 0.4549
    return round(float(obs_median_chi2 / expected_median_chi2), 4)


# ---------------------------------------------------------------------------
# Single-axis QQ plot
# ---------------------------------------------------------------------------

def plot_qq_single(
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
    thin: bool = False,
    thin_below: float = 0.01,
    max_points: int = 50_000,
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
            "  plot_qq_single(pvals, ax=axes[0])"
        )
 

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

    # CI band
    ax.fill_between(
        expected, ci_lo, ci_hi,
        color=color, alpha=ci_alpha, linewidth=0,
        label=f"{int(ci * 100)}% CI",
    )

    # Diagonal null line
    max_val = max(expected.max(), observed.max()) * 1.05
    ax.plot([0, max_val], [0, max_val], color="grey", linewidth=0.8,
            linestyle="--", zorder=1)

    # Observed points
    ax.scatter(
        expected, observed,
        s=point_size, color=color, alpha=0.85,
        label=label, zorder=2, edgecolors="none",
        rasterized=rasterized,
    )

    '''"""
    # Significance line
    if signif_threshold is not None:
        sig_logp = -np.log10(signif_threshold)
        ax.axhline(sig_logp, color="red", linewidth=0.7, linestyle="--",
                   label=f"p={signif_threshold:.0e}")
    """'''

    # Lambda annotation
    if show_lambda and not math.isnan(lam):
        ax.text(
            0.05, 0.95,
            f"λ = {lam:.4f}",
            transform=ax.transAxes,
            va="top", ha="left",
            fontsize=9, fontstyle="italic",
            color="black",
        )

    ax.set_xlabel("Expected −log₁₀(p)", fontsize=10)
    ax.set_ylabel("Observed −log₁₀(p)", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if title:
        ax.set_title(title, fontsize=10, pad=6)
    if label:
        ax.legend(fontsize=8, frameon=False, loc="lower right")

    return ax


# ---------------------------------------------------------------------------
# Combined multi-panel figure
# ---------------------------------------------------------------------------

def plot_qq_combined(
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
    thin: bool = False,
    thin_below: float = 0.01,
    max_points: int = 50_000,
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
        See :func:`plot_qq_single`.

    Returns
    -------
    (fig, axes)
    """
    n = len(pval_dict)
    if n == 0:
        raise ValueError("pval_dict is empty.")

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

    fig, axes_grid = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    axes_flat = axes_grid.flatten()

    for idx, (label, pvals) in enumerate(pval_dict.items()):
        plot_qq_single(
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
            thin_below=thin_below,
            max_points=max_points,
            rasterized=rasterized,
        )

    for ax in axes_flat[n:]:
        ax.set_visible(False)

    if title:
        fig.suptitle(title, fontsize=13, y=1.01)

    plt.tight_layout()

    if output_path:
        fmt = fig_format or Path(output_path).suffix.lstrip(".") or "png"
        fig.savefig(f"{output_path}.{fmt}", format=fmt, dpi=dpi, bbox_inches="tight")
        logger.info("Saved combined QQ plot: %s", f"{output_path}.{fmt}")

    return fig, list(axes_flat[:n])


# ---------------------------------------------------------------------------
# Separate figures — one file per sumstat
# ---------------------------------------------------------------------------

def plot_qq_separate(
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
    thin: bool = False,
    thin_below: float = 0.01,
    max_points: int = 50_000,
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
        See :func:`plot_qq_single`.

    Returns
    -------
    List of output file paths.
    """

    labels = pval_dict.keys()

    # plot name
    (
        plt_name, 
        table_out,
        plt_base,
    ) = get_output_paths(
        labels = labels,
        mode='qq',
        output_dir=output_path, 
        plot_title=base_name, 
        output_format=fig_format
    )

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

        plot_qq_single(
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
            thin_below=thin_below,
            max_points=max_points,
            rasterized=rasterized,
        )

        plt.tight_layout()

        safe_label = label.replace(" ", "_").replace("/", "-")
        out_path = f"{plt_base}_{safe_label}.{fig_format}"
        fig.savefig(out_path, format=fig_format, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        logger.info("Saved QQ plot: %s", out_path)
        saved.append(out_path)

    return saved


# ---------------------------------------------------------------------------
# Overlay — all sumstats on one axes
# ---------------------------------------------------------------------------

def plot_qq_overlay(
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
    thin: bool = False,
    thin_below: float = 0.01,
    max_points: int = 50_000,
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
        See :func:`plot_qq_single`.

    Returns
    -------
    (fig, ax)
    """

    labels = pval_dict.keys()

    # plot name
    (
        plt_name, 
        table_out,
        plt_base,
    ) = get_output_paths(
        labels = labels,
        mode='qq',
        output_dir=output_path,
        plot_title=title, 
        output_format=fig_format
    )

    n = len(pval_dict)
    if n == 0:
        raise ValueError("pval_dict is empty.")


    cmap = plt.get_cmap("tab10")
    colors = [mcolors.to_hex(cmap(i % 10)) for i in range(n)]
    #if colors is None:
    #    cmap = plt.get_cmap("tab10")
    #    colors = [mcolors.to_hex(cmap(i % 10)) for i in range(n)]
    #elif len(colors) < n:
    #    colors = [colors[i % len(colors)] for i in range(n)]

    fig, ax = plt.subplots(figsize=figsize)
    global_max = 0.0

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
        legend_label = f"{label}  (λ={lam:.4f})" if show_lambda else label

        ax.fill_between(
            expected, ci_lo, ci_hi,
            color=color, alpha=ci_alpha, linewidth=0,
        )
        ax.scatter(
            expected, observed,
            s=point_size, color=color, alpha=0.85,
            label=legend_label, zorder=2 + idx, edgecolors="none",
            rasterized=rasterized,
        )

        global_max = max(global_max, expected.max(), observed.max())

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
        fontsize=8, frameon=True, framealpha=0.7,
        edgecolor="lightgrey", loc="lower right",
    )

    if title:
        ax.set_title(title, fontsize=11, pad=8)

    plt.tight_layout()

    if output_path:
        fmt = fig_format or Path(output_path).suffix.lstrip(".") or "png"
        fig.savefig(f"{plt_base}.{fmt}", format=fmt, dpi=dpi, bbox_inches="tight")
        logger.info("Saved overlay QQ plot: %s", f"{plt_base}.{fmt}")

    return fig, ax
