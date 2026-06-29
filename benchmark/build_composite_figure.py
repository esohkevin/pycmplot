#!/usr/bin/env python3
"""Build the composite scaling figure used in the manuscript.

Output: ``benchmark_composite_scaling.png`` (and an accompanying PDF), written
to the directory passed via ``--outdir`` (default: same directory as this
script).

Layout (2 rows x 3 columns):

  Top row    -- single-track wall-time scaling at 500K-10M variants:
    (A) Manhattan
    (B) Circular Manhattan
    (C) QQ

  Bottom row -- multi-track wall-time scaling at 1M and 2M variants
                (multi-track benchmarks were only collected at 1M / 2M):
    (D) Multi-track Manhattan          (pycmplot vs CMplot)
    (E) Multi-track circular Manhattan (pycmplot vs CMplot)
    (F) Shared legend

The figure reuses :data:`collect_results.SERIES_STYLE` so styling stays
consistent with the per-plot-type PDFs produced by
``collect_results.py --plot``.

Usage
-----
::

    cd benchmark
    python collect_results.py --resultsdir results --out results/summary.csv
    python build_composite_figure.py                      # writes alongside this script
    python build_composite_figure.py --outdir ../figures  # custom output dir
    python build_composite_figure.py --summary path/to/summary.csv
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# Resolve ``collect_results`` relative to this script so the figure builder
# works regardless of the current working directory.
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from collect_results import SERIES_STYLE, SIZE_N  # noqa: E402


PANELS = [
    # (row, col, letter, plot_type, title)
    (0, 0, "A", "manhattan",            "Manhattan plot"),
    (0, 1, "B", "circular",             "Circular Manhattan plot"),
    (0, 2, "C", "qq",                   "QQ plot"),
    (1, 0, "D", "multitrack_manhattan", "Multi-track Manhattan plot"),
    (1, 1, "E", "multitrack_circular",  "Multi-track circular Manhattan plot"),
]


def build(summary_path: Path, outdir: Path, basename: str = "benchmark_composite_scaling") -> Path:
    """Render the composite figure and return the path to the PNG output."""
    summary = pd.read_csv(summary_path)

    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5))
    plt.subplots_adjust(
        top=0.94, bottom=0.07, left=0.06, right=0.98, wspace=0.30, hspace=0.45,
    )

    all_handles: list = []
    all_labels: list = []
    seen: set = set()

    for r, c, letter, pt, title in PANELS:
        ax = axes[r, c]
        sub = summary[summary["plot_type"] == pt].copy()
        sub["n_approx"] = sub["size_label"].map(SIZE_N)
        sub = sub.dropna(subset=["n_approx"])

        tools_in_data = list(sub["tool"].unique())
        ordered_tools = (
            [t for t in SERIES_STYLE if t in tools_in_data]
            + [t for t in tools_in_data if t not in SERIES_STYLE]
        )

        for tool in ordered_tools:
            tdf = sub[sub["tool"] == tool].sort_values("n_approx")
            if tdf.empty:
                continue
            s = SERIES_STYLE.get(
                tool, dict(label=tool, color="#888780", ls="-", marker="o"),
            )
            line, = ax.plot(
                tdf["n_approx"], tdf["wall_time_mean"],
                marker=s["marker"], linestyle=s["ls"], color=s["color"],
                label=s["label"] or tool,
                linewidth=1.6, markersize=4.5,
            )
            ax.fill_between(
                tdf["n_approx"],
                tdf["wall_time_mean"] - tdf["wall_time_sd"].fillna(0),
                tdf["wall_time_mean"] + tdf["wall_time_sd"].fillna(0),
                alpha=0.12, color=s["color"], linewidth=0,
            )
            label = s["label"] or tool
            if label not in seen:
                all_handles.append(line)
                all_labels.append(label)
                seen.add(label)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Number of variants", fontsize=10)
        if c == 0:
            ax.set_ylabel("Wall-clock time (s)", fontsize=10)
        ax.set_title(f"({letter}) {title}", fontsize=11, loc="left")
        ax.grid(True, which="both", linestyle="--", linewidth=0.3, alpha=0.5)
        ax.tick_params(labelsize=9)

    # Bottom-right cell (1, 2) holds the shared legend so each axes stays
    # uncluttered.
    legend_ax = axes[1, 2]
    legend_ax.axis("off")
    legend_ax.legend(
        all_handles, all_labels,
        loc="center", ncol=1, fontsize=10, frameon=False,
        title="Tool", title_fontsize=11,
    )

    outdir.mkdir(parents=True, exist_ok=True)
    png_path = outdir / f"{basename}.png"
    pdf_path = outdir / f"{basename}.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote {png_path}")
    print(f"Wrote {pdf_path}")
    return png_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build the composite scaling figure for the manuscript.",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=None,
        help=(
            "Path to summary.csv produced by collect_results.py.  "
            "Default: <script_dir>/summary.csv, then <script_dir>/results/summary.csv."
        ),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=SCRIPT_DIR,
        help="Output directory for the figure files (default: this script's directory).",
    )
    parser.add_argument(
        "--basename",
        type=str,
        default="benchmark_composite_scaling",
        help="Output file basename (without extension).  Default: benchmark_composite_scaling.",
    )
    args = parser.parse_args()

    if args.summary is None:
        candidates = [SCRIPT_DIR / "summary.csv", SCRIPT_DIR / "results" / "summary.csv"]
        for c in candidates:
            if c.exists():
                args.summary = c
                break
        else:
            sys.exit(
                "Could not locate summary.csv.  Pass --summary or run "
                "`python collect_results.py --resultsdir results --out results/summary.csv` first."
            )

    if not args.summary.exists():
        sys.exit(f"summary file not found: {args.summary}")

    build(args.summary, args.outdir, basename=args.basename)


if __name__ == "__main__":
    main()
