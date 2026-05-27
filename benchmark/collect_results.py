#!/usr/bin/env python3
"""
collect_results.py
Aggregates benchmark CSVs and produces a publication-ready summary.

Usage:
    python collect_results.py --resultsdir results/ --out summary.csv
    python collect_results.py --resultsdir results/ --out summary.csv --plot
"""

import argparse
import os
import sys
import warnings

import pandas as pd
import numpy as np


SIZE_ORDER = ["500K", "1M", "2M", "5M", "10M"]
TOOL_ORDER = ["pycmplot", "gwaslab", "qmplot", "CMplot", "qqman"]

# Approximate variant counts for each size label (matches generate_sumstats.py)
SIZE_N = {"500K": 500_000, "1M": 1_000_000, "2M": 2_000_000,
          "5M": 5_000_000, "10M": 10_000_000}


# ---------------------------------------------------------------------------
# Plot-series styling
# ---------------------------------------------------------------------------
# Plotted series (one entry per tool key found in the CSVs).  pycmplot
# variants share a common colour and are distinguished by linestyle so that
# the trim-related speedups read as an ablation on top of pycmplot, rather
# than as separate competing tools.  Order here determines legend order.
SERIES_STYLE = {
    # --- pycmplot family ---------------------------------------------------
    "pycmplot":
        dict(label="pycmplot",                              color="#534AB7", ls="-",  marker="o"),
    "pycmplot_trimmed_p0.01":
        dict(label="pycmplot (--trim_p 0.01)",              color="#534AB7", ls="--", marker="o"),
    "pycmplot_trimmed_p0.001":
        dict(label="pycmplot (--trim_p 0.001)",             color="#534AB7", ls=":",  marker="o"),
    # multitrack variants share pycmplot's colour but use their own keys so
    # they only show on multitrack_* plot_types (where no untrimmed key
    # exists; the multitrack column in bench_python.csv uses the
    # "pycmplot_multitrack" tool label).
    "pycmplot_multitrack":
        dict(label="pycmplot (multitrack)",                 color="#534AB7", ls="-",  marker="o"),
    "pycmplot_multitrack_trimmed_p0.01":
        dict(label="pycmplot multitrack (--trim_p 0.01)",   color="#534AB7", ls="--", marker="o"),
    # --- Other Python tools ------------------------------------------------
    "gwaslab":
        dict(label="gwaslab",                               color="#1D9E75", ls="-",  marker="s"),
    "qmplot":
        dict(label="qmplot",                                color="#378ADD", ls="-",  marker="^"),
    # --- R tools -----------------------------------------------------------
    "CMplot":
        dict(label="CMplot",                                color="#D85A30", ls="-",  marker="D"),
    "CMplot_multitrack":
        dict(label="CMplot (multitrack)",                   color="#D85A30", ls="-",  marker="D"),
    "qqman":
        dict(label="qqman",                                 color="#888780", ls="-",  marker="v"),
}


# Friendly display names for plot titles
PLOT_TYPE_TITLE = {
    "manhattan":            "Manhattan plot",
    "circular":             "Circular Manhattan plot",
    "qq":                   "QQ plot",
    "multitrack_manhattan": "Multi-track Manhattan plot",
    "multitrack_circular":  "Multi-track circular Manhattan plot",
}


def load_results(resultsdir: str) -> pd.DataFrame:
    frames = []
    for fname in (
        "bench_python.csv", "bench_r.csv",
        #"pycmplot_trimmed_bench_python.csv",
        #"pycmplot_multitrack_multi_bench_python.csv",
        #"pycmplot_multitrack_multi_trimmed_bench_python.csv",
        #"gwaslab_bench_python.csv",
        #"qmplot_bench_python.csv",
        #"CMplot_bench_r.csv",
        #"CMplot_multitrack_bench_r.csv",
        #"qqman_bench_r.csv",
    ):
        fpath = os.path.join(resultsdir, fname)
        if os.path.exists(fpath):
            df = pd.read_csv(fpath)
            frames.append(df)
        else:
            print(f"[warn] {fpath} not found — skipping")

    if not frames:
        sys.exit("No result files found. Run benchmarks first.")

    df = pd.concat(frames, ignore_index=True)

    # Drop failed rows
    n_before = len(df)
    df = df[df["wall_time_s"].astype(str) != "ERROR"].copy()
    n_dropped = n_before - len(df)
    if n_dropped:
        print(f"[warn] Dropped {n_dropped} failed rows")

    df["wall_time_s"] = pd.to_numeric(df["wall_time_s"])
    df["peak_mem_mb"] = pd.to_numeric(df["peak_mem_mb"])
    df["out_file_kb"] = pd.to_numeric(df["out_file_kb"])
    df["n_variants"]  = pd.to_numeric(df["n_variants"])

    # Normalise multitrack labelling.  Some tools (e.g. ``CMplot_multitrack``)
    # tag their rows with ``plot_type='manhattan'`` / ``'circular'`` while
    # pycmplot's multitrack rows use ``'multitrack_manhattan'`` /
    # ``'multitrack_circular'``.  Retag the former on load so cross-tool
    # multitrack comparisons land on the multitrack figures rather than
    # contaminating the single-track ones.
    mt_mask = (
        df["tool"].astype(str).str.contains("_multitrack", na=False)
        & ~df["plot_type"].astype(str).str.startswith("multitrack_")
    )
    if mt_mask.any():
        df.loc[mt_mask, "plot_type"] = (
            "multitrack_" + df.loc[mt_mask, "plot_type"].astype(str)
        )

    # Drop any exact duplicate (tool, plot_type, size_label, replicate) rows
    # that could sneak in when the same benchmark was logged to several CSVs
    # (e.g. a per-tool CSV *and* the aggregated ``bench_python.csv``).
    # Without this, ``summarize`` keeps both as separate replicates and
    # ``pivot_table`` later fails with "Index contains duplicate entries".
    before = len(df)
    df = df.drop_duplicates(
        subset=["tool", "plot_type", "size_label", "replicate"],
        keep="last",
    )
    dropped = before - len(df)
    if dropped:
        print(f"[info] Dropped {dropped} duplicate replicate row(s)")

    return df


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute mean ± SD per (tool, plot_type, size_label).
    """
    grp = df.groupby(["tool", "plot_type", "size_label"])

    summary = grp.agg(
        n_replicates      = ("replicate",    "count"),
        wall_time_mean    = ("wall_time_s",  "mean"),
        wall_time_sd      = ("wall_time_s",  "std"),
        wall_time_min     = ("wall_time_s",  "min"),
        wall_time_max     = ("wall_time_s",  "max"),
        peak_mem_mean_mb  = ("peak_mem_mb",  "mean"),
        peak_mem_sd_mb    = ("peak_mem_mb",  "std"),
        out_file_kb_mean  = ("out_file_kb",  "mean"),
    ).reset_index()

    # Categorical ordering for clean display.  Canonical tools sort first;
    # trimmed / multitrack variants come after but are *kept* (otherwise the
    # restricted-category cast would silently turn them into NaN, which
    # later breaks pivot_table with "Index contains duplicate entries").
    all_tools = list(TOOL_ORDER) + [
        t for t in summary["tool"].unique() if t not in TOOL_ORDER
    ]
    summary["tool"]       = pd.Categorical(summary["tool"],       categories=all_tools, ordered=True)
    summary["size_label"] = pd.Categorical(summary["size_label"], categories=SIZE_ORDER, ordered=True)
    summary = summary.sort_values(["plot_type", "tool", "size_label"])

    # Formatted string for paper table: "mean ± sd"
    summary["time_str"] = summary.apply(
        lambda r: f"{r.wall_time_mean:.1f} ± {r.wall_time_sd:.1f}", axis=1
    )
    summary["mem_str"] = summary.apply(
        lambda r: f"{r.peak_mem_mean_mb:.0f} ± {r.peak_mem_sd_mb:.0f}", axis=1
    )

    return summary


def pivot_table(summary: pd.DataFrame, plot_type: str, metric: str = "time_str") -> pd.DataFrame:
    """
    Produce a tools × sizes pivot table for a given plot type and metric.
    Suitable for direct inclusion in a manuscript table.
    """
    sub = summary[summary["plot_type"] == plot_type].copy()
    pivot = sub.pivot(index="tool", columns="size_label", values=metric)
    pivot = pivot.reindex(index=[t for t in TOOL_ORDER if t in pivot.index],
                          columns=[s for s in SIZE_ORDER if s in pivot.columns])
    return pivot


def speedup_table(summary: pd.DataFrame, baseline: str = "CMplot",
                  plot_type: str = "manhattan") -> pd.DataFrame:
    """
    Compute speedup of each tool relative to a baseline (default CMplot).
    Useful for the paper's key claim of Python speed advantage.
    """
    sub = summary[(summary["plot_type"] == plot_type)].copy()
    base = sub[sub["tool"] == baseline][["size_label", "wall_time_mean"]].rename(
        columns={"wall_time_mean": "base_time"})
    merged = sub.merge(base, on="size_label")
    merged["speedup_vs_CMplot"] = merged["base_time"] / merged["wall_time_mean"]
    merged["speedup_vs_CMplot"] = merged["speedup_vs_CMplot"].round(1)

    pivot = merged.pivot(index="tool", columns="size_label", values="speedup_vs_CMplot")
    pivot = pivot.reindex(index=[t for t in TOOL_ORDER if t in pivot.index],
                          columns=[s for s in SIZE_ORDER if s in pivot.columns])
    return pivot


def _plot_metric(
    summary: pd.DataFrame,
    outdir: str,
    *,
    mean_col: str,
    sd_col: str,
    ylabel: str,
    title_prefix: str,
    file_suffix: str,
):
    """Internal: render one scaling plot per plot_type for a chosen metric.

    Parameters
    ----------
    summary
        Aggregated benchmark summary as produced by :func:`summarize`.
    outdir
        Output directory; created if missing.
    mean_col, sd_col
        Column names in *summary* for the metric's mean and standard
        deviation (e.g. ``"wall_time_mean"`` / ``"wall_time_sd"``).
    ylabel
        Axis label for the metric.
    title_prefix
        Prefix for the per-figure title (e.g. ``"Runtime scaling"`` or
        ``"Peak memory scaling"``).  Combined with the plot-type
        display name from :data:`PLOT_TYPE_TITLE`.
    file_suffix
        Suffix appended to each output PDF name (empty string for the
        canonical wall-time plots; ``"_memory"`` for memory plots).
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[warn] matplotlib not available — skipping plots")
        return

    os.makedirs(outdir, exist_ok=True)
    fallback_style = dict(label=None, color="#888780", ls="-", marker="o")

    for plot_type in summary["plot_type"].unique():
        sub = summary[summary["plot_type"] == plot_type].copy()
        sub["n_approx"] = sub["size_label"].map(SIZE_N)
        sub = sub.dropna(subset=["n_approx", mean_col])
        if sub.empty:
            continue

        # Iterate in SERIES_STYLE order for stable legend ordering, then
        # append any tools not yet styled (so new entries are not silently
        # dropped).
        tools_in_data = list(sub["tool"].unique())
        ordered_tools = [t for t in SERIES_STYLE if t in tools_in_data] + \
                        [t for t in tools_in_data if t not in SERIES_STYLE]
        if not ordered_tools:
            continue

        fig, ax = plt.subplots(figsize=(7, 4.4))

        for tool in ordered_tools:
            tdf = sub[sub["tool"] == tool].sort_values("n_approx")
            if tdf.empty or tdf[mean_col].isna().all():
                continue
            style = SERIES_STYLE.get(tool, dict(fallback_style, label=tool))
            ax.plot(
                tdf["n_approx"], tdf[mean_col],
                marker=style["marker"],
                linestyle=style["ls"],
                color=style["color"],
                label=style["label"] or tool,
                linewidth=1.8,
                markersize=5,
            )
            ax.fill_between(
                tdf["n_approx"],
                tdf[mean_col] - tdf[sd_col].fillna(0),
                tdf[mean_col] + tdf[sd_col].fillna(0),
                alpha=0.12, color=style["color"], linewidth=0,
            )

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Number of variants", fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        title = PLOT_TYPE_TITLE.get(plot_type, plot_type)
        ax.set_title(f"{title_prefix} — {title}", fontsize=11)
        ax.legend(frameon=False, fontsize=8, loc="upper left")
        ax.grid(True, which="both", linestyle="--", linewidth=0.4, alpha=0.5)
        fig.tight_layout()

        out_path = os.path.join(outdir, f"benchmark_{plot_type}{file_suffix}.pdf")
        fig.savefig(out_path, dpi=200)
        plt.close(fig)
        print(f"  Saved: {out_path}")


def plot_results(summary: pd.DataFrame, outdir: str):
    """
    Generate scaling plots (one per plot_type) for both wall-clock time
    and peak memory.

    Wall-time plots are written to ``benchmark_<plot_type>.pdf`` and the
    matching memory plots to ``benchmark_<plot_type>_memory.pdf``.  Both
    share :data:`SERIES_STYLE`, so pycmplot's trimmed variants stay
    visually grouped with canonical pycmplot in both views.
    """
    _plot_metric(
        summary, outdir,
        mean_col="wall_time_mean",
        sd_col="wall_time_sd",
        ylabel="Wall-clock time (s)",
        title_prefix="Runtime scaling",
        file_suffix="",
    )
    _plot_metric(
        summary, outdir,
        mean_col="peak_mem_mean_mb",
        sd_col="peak_mem_sd_mb",
        ylabel="Peak memory (MB)",
        title_prefix="Peak memory scaling",
        file_suffix="_memory",
    )


def main():
    parser = argparse.ArgumentParser(description="Aggregate benchmark results")
    parser.add_argument("--resultsdir", default="results")
    parser.add_argument("--out",        default="results/summary.csv")
    parser.add_argument("--plot",       action="store_true",
                        help="Generate scaling plots (requires matplotlib)")
    parser.add_argument("--plotdir",    default="results/plots")
    args = parser.parse_args()

    print("Loading results...")
    df = load_results(args.resultsdir)
    print(f"  {len(df)} replicate rows across {df['tool'].nunique()} tools")

    summary = summarize(df)
    summary.to_csv(args.out, index=False)
    print(f"\nFull summary saved to: {args.out}")

    # Print paper-ready tables
    for plot_type in sorted(summary["plot_type"].unique()):
        print(f"\n=== {plot_type.upper()} — Wall-clock time (s, mean ± SD) ===")
        tbl = pivot_table(summary, plot_type=plot_type, metric="time_str")
        print(tbl.to_string())

        print(f"\n=== {plot_type.upper()} — Peak memory (MB, mean ± SD) ===")
        tbl2 = pivot_table(summary, plot_type=plot_type, metric="mem_str")
        print(tbl2.to_string())

    print(f"\n=== SPEEDUP vs CMplot (manhattan) ===")
    if "CMplot" in summary["tool"].values:
        sp = speedup_table(summary, baseline="CMplot", plot_type="manhattan")
        print(sp.to_string())
    else:
        print("  CMplot results not yet available")

    if args.plot:
        print(f"\nGenerating scaling plots -> {args.plotdir}/")
        plot_results(summary, args.plotdir)


if __name__ == "__main__":
    main()
