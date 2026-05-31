#!/usr/bin/env python3
"""
bench_python.py
Benchmarks Python GWAS visualization tools: pycmplot, gwaslab, qmplot.

Usage:
    python bench_python.py --tool pycmplot --input data/sumstats_1M.tsv \
        --size 1M --replicates 5 --outdir results/

Writes one CSV row per replicate to results/<tool>_trimmed_bench_python.csv.
"""

import argparse
import csv
import gc
import os
import sys
import time
import tracemalloc
from pathlib import Path

RESULT_COLS = [
    "tool", "plot_type", "size_label", "n_variants",
    "replicate", "wall_time_s", "peak_mem_mb", "out_file_kb"
]


def _record(writer, row: dict):
    writer.writerow({k: row.get(k, "") for k in RESULT_COLS})


# ---------------------------------------------------------------------------
# Individual tool wrappers
# Each wrapper must:
#   1. Load data from disk (include I/O in timing)
#   2. Produce a PNG to out_path (or nearby with an appended suffix)
#   3. Return nothing
# ---------------------------------------------------------------------------

def run_pycmplot(input_path: str, out_path: str, plot_type: str = "manhattan", trim: float = False):
    """
    pycmplot benchmark.

    pycmplot appends an informative suffix to whatever stem is passed as
    plot_title, so the file written to disk will differ from out_path.
    benchmark_one() resolves the actual output via a directory snapshot.
    """
    import pycmplot

    sumstat = [input_path]
    label = ["Trait"]

    out_path_base = out_path.replace(".png","")
    out_parts = out_path_base.rsplit('/', 1)
    out_dir   = out_parts[0]
    out_file_parts = ["pycmplot"] + out_parts[1].split("_")[2:6]
    out_file = "_".join(out_file_parts)

    sumstats_info_dict = pycmplot.prep_pycmplot_input_info(
        sum_stats=sumstat,
        labels=label,
        chrom="CHR",
        pos="BP",
        pcol="P",
        snp="SNP",
    )

    if trim:
        pycmplot_dict = pycmplot.get_sumstats_and_merged_sector_list(
            sum_stats=sumstat,
            labels=label,
            file_info=sumstats_info_dict,
            logp=True,
            trim_pval=trim,
        )
    else:
        pycmplot_dict = pycmplot.get_sumstats_and_merged_sector_list(
            sum_stats=sumstat,
            labels=label,
            file_info=sumstats_info_dict,
            logp=True,
            #trim_pval=0.001,
        )        

    if plot_type == "manhattan":
        pycmplot.plot_linear(
            sumstats_loaded=pycmplot_dict["dfs"],
            plot_title=out_file,
            output_dir=out_dir,
            output_format='png',
            dpi=150,
            logp=True,
        )
    elif plot_type == "circular":
        pycmplot.plot_circular(
            sumstats_loaded=pycmplot_dict["dfs"],
            sector_sizes=pycmplot_dict["sectors"],
            logp=True,
            plot_title=out_file,
            output_dir=out_dir,
            dpi=150,
            output_format='png',
            pad=1,
            r_min=50,
            r_max=100,
        )
    elif plot_type == "qq":
        if trim:
            pycmplot.plot_qq_overlay(
                pval_dict=pycmplot_dict["pvals"],
                thin=True,
                thin_below=0.001,
                title=out_file,
                output_path=out_dir,
            )
        else:
            pycmplot.plot_qq_overlay(
                pval_dict=pycmplot_dict["pvals"],
                #thin=True,
                #thin_below=0.001,
                title=out_file,
                output_path=out_dir,
            )



def run_gwaslab(input_path: str, out_path: str, plot_type: str = "manhattan", trim: float = False):
    """
    gwaslab benchmark.
    https://cloufield.github.io/gwaslab/
    """
    import pandas as pd
    import gwaslab as gl
    import matplotlib
    matplotlib.use("Agg")

    df = pd.read_csv(input_path, sep="\t")

    mysumstats = gl.Sumstats(
        df,
        snpid="SNP",
        chrom="CHR",
        pos="BP",
        p="P",
        beta="BETA",
        se="SE",
        ea="A1",
        nea="A2",
    )

    if plot_type == "manhattan":
        mysumstats.plot_mqq(
            mode="m",
            save=out_path,
            save_kwargs={"dpi": 150},
            verbose=False,
        )
    elif plot_type == "qq":
        mysumstats.plot_mqq(
            mode="qq",
            save=out_path,
            save_kwargs={"dpi": 150},
            verbose=False,
        )


def run_qmplot(input_path: str, out_path: str, plot_type: str = "manhattan", trim: float = False):
    """
    qmplot benchmark.
    https://github.com/ShujiaHuang/qmplot
    """
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from qmplot import manhattanplot, qqplot

    df = pd.read_csv(input_path, sep="\t")

    fig, ax = plt.subplots(figsize=(12, 4))

    if plot_type in ("manhattan", "circular"):
        manhattanplot(
            data=df,
            chrom="CHR",
            pos="BP",
            pv="P",
            snp="SNP",
            ax=ax,
        )
    elif plot_type == "qq":
        qqplot(
            data=df["P"],
            ax=ax,
        )

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)



def run_pycmplot_multitrack(manifest_path: str, out_path: str, plot_type: str = "multitrack_manhattan", trim: float = False):
    """
    pycmplot multi-track benchmark.

    Reads trait TSV paths from a manifest file (one path per line) and passes
    them as a list to pycmplot — no merging required on the user's part.
    This is the key input-format advantage over CMplot, which requires a
    column-wise merge performed manually before plotting (see bench_r.R).

    Timing therefore covers: manifest parse + N file loads + plot + save.
    """
    import pycmplot

    with open(manifest_path) as fh:
        sumstat = [line.strip() for line in fh if line.strip()]
    n_traits = len(sumstat)
    labels   = [f"Trait{k+1}" for k in range(n_traits)]

    out_path = out_path.replace(".png","")
    out_parts = out_path.rsplit('/', 1)
    out_dir   = out_parts[0]
    out_file_parts = ["pycmplot"] + out_parts[1].split("_")[2:6]
    out_file = "_".join(out_file_parts)
    

    sumstats_info_dict = pycmplot.prep_pycmplot_input_info(
        sum_stats=sumstat,
        labels=labels,
        chrom="CHR",
        pos="BP",
        pcol="P",
        snp="SNP"
    )

    pycmplot_dict = pycmplot.get_sumstats_and_merged_sector_list(
        sum_stats=sumstat,
        labels=labels,
        file_info=sumstats_info_dict,
        logp=True,
        trim_pval=0.001,
    )

    if plot_type == "multitrack_manhattan":
        pycmplot.plot_linear(
            sumstats_loaded=pycmplot_dict["dfs"],
            plot_title=out_file,
            output_dir=out_dir,
            output_format='png',
            dpi=150,
            logp=True,
            signif_lines=pycmplot_dict['lines'],
        )
    elif plot_type == "multitrack_circular":
        pycmplot.plot_circular(
            sumstats_loaded=pycmplot_dict["dfs"],
            sector_sizes=pycmplot_dict["sectors"],
            logp=True,
            signif_lines=pycmplot_dict['lines'],
            plot_title=out_file,
            output_dir=out_dir,
            dpi=150,
            output_format='png',
            pad=1,
            r_min=40,
            r_max=100,
        )

# ---------------------------------------------------------------------------
# Timing + memory harness
# ---------------------------------------------------------------------------

TOOL_RUNNERS = {
    "pycmplot":            run_pycmplot,
    "pycmplot_multitrack": run_pycmplot_multitrack,
    "gwaslab":            run_gwaslab,
    "qmplot":             run_qmplot,
}


def benchmark_one(tool: str, input_path: str, out_path: str, plot_type: str, trim: float):
    """
    Run one timed, memory-tracked benchmark call.

    pycmplot appends an informative suffix to the stem passed as plot_title,
    so the file written to disk differs from out_path. We snapshot the output
    directory before and after the call to find whatever was actually written
    and return that path so the caller can measure its size correctly.

    Tools that honour out_path exactly (gwaslab, qmplot) are unaffected —
    the pre/post snapshot simply confirms out_path appeared and returns it.

    Returns
    -------
    tuple[float, float, str]
        (wall_time_seconds, peak_memory_mb, actual_output_path)
    """
    runner = TOOL_RUNNERS[tool]
    out_dir = os.path.dirname(os.path.abspath(out_path))

    # Snapshot PNGs present before the call
    before = set(
        os.path.join(out_dir, f)
        for f in os.listdir(out_dir)
        if f.lower().endswith(".png")
    )

    # Force a full GC cycle before each run so prior allocations don't inflate
    gc.collect()

    tracemalloc.start()
    t0 = time.perf_counter()

    runner(input_path, out_path, plot_type, trim)

    t1 = time.perf_counter()
    _, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    wall_time   = t1 - t0
    peak_mem_mb = peak_bytes / 1024 / 1024

    # Resolve the actual output file written
    after     = set(
        os.path.join(out_dir, f)
        for f in os.listdir(out_dir)
        if f.lower().endswith(".png")
    )
    new_files = after - before

    if os.path.exists(out_path):
        # Tool honoured out_path exactly (gwaslab, qmplot, pycmplot qq mode)
        actual_out = out_path
    elif len(new_files) == 1:
        # Tool wrote a single new file with its own naming (pycmplot manhattan/circular)
        actual_out = next(iter(new_files))
    elif len(new_files) > 1:
        # Multiple new files — prefer the one whose name starts with out_path's stem
        stem = os.path.splitext(os.path.basename(out_path))[0]
        matches = [f for f in new_files if os.path.basename(f).startswith(stem)]
        actual_out = matches[0] if matches else next(iter(new_files))
    else:
        # Nothing written — return out_path so caller's exists() check gives False
        actual_out = out_path

    return wall_time, peak_mem_mb, actual_out


def main():
    parser = argparse.ArgumentParser(description="Benchmark Python GWAS visualization tools")
    parser.add_argument("--tool",        required=True, choices=list(TOOL_RUNNERS.keys()))
    parser.add_argument("--input",       default=None,  help="Path to sumstats TSV (single-trait runs)")
    parser.add_argument("--manifest",    default=None,
                        help="Path to manifest file listing trait TSV paths, one per line "
                             "(multi-track runs — used with pycmplot_multitrack)")
    parser.add_argument("--size",        required=True, help="Dataset size label (e.g. 1M)")
    parser.add_argument("--plot-type",   default="manhattan",
                        choices=["manhattan", "circular", "qq",
                                 "multitrack_manhattan", "multitrack_circular"],
                        help="Plot type to benchmark")
    parser.add_argument("--replicates",  type=int, default=5)
    parser.add_argument("--trim", nargs='?', type=float, default=False, const=0.1, required=False)
    parser.add_argument("--outdir",      default="results", help="Directory for CSV results")
    parser.add_argument("--figdir",      default="figures", help="Directory for generated figures")
    args = parser.parse_args()

    # Validate: multi-track plot types require --manifest; single-trait require --input
    is_multitrack = args.plot_type.startswith("multitrack")
    if is_multitrack and not args.manifest:
        parser.error("--manifest is required for multitrack plot types")
    if not is_multitrack and not args.input:
        parser.error("--input is required for single-trait plot types")

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.figdir, exist_ok=True)

    if args.trim:
        bench_res = args.tool + f"_bench_python_trim{args.trim}.csv"
    else:
        bench_res = args.tool + "_bench_python.csv"

    #bench_res = args.tool + "_bench_python.csv"
    #bench_res = args.tool + "_trimmed_p0.001_bench_python.csv"
    #bench_res = args.tool + "_multi_trimmed_bench_python.csv"
    csv_path = os.path.join(args.outdir, bench_res)
    write_header = not os.path.exists(csv_path)

    # Resolve input path and variant count
    # Multi-track: use manifest path as the "input" passed to the runner;
    # report variant count from the first trait file in the manifest.
    if is_multitrack:
        input_for_bench = args.manifest
        with open(args.manifest) as mf:
            first_trait = next(l.strip() for l in mf if l.strip())
        n_variants = sum(1 for _ in open(first_trait)) - 1
        with open(args.manifest) as mf:
            n_traits = sum(1 for l in mf if l.strip())
    else:
        input_for_bench = args.input
        n_variants = sum(1 for _ in open(args.input)) - 1
        n_traits = 1

    print(f"\n[bench] tool={args.tool}  size={args.size}  "
          f"n={n_variants:,}  traits={n_traits}  plot={args.plot_type}  reps={args.replicates}")

    with open(csv_path, "a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=RESULT_COLS)
        if write_header:
            writer.writeheader()

        for rep in range(1, args.replicates + 1):
            if args.trim:
                out_fig = os.path.join(
                    args.figdir,
                    f"{args.tool}_{args.plot_type}_{args.size}_trim{args.trim}_rep{rep}.png"
                )
            else:                
                out_fig = os.path.join(
                    args.figdir,
                    f"{args.tool}_{args.plot_type}_{args.size}_rep{rep}.png"
                )

            try:
                wall, mem, actual_fig = benchmark_one(
                    args.tool, input_for_bench, out_fig, args.plot_type, args.trim
                )
                out_kb = os.path.getsize(actual_fig) / 1024 if os.path.exists(actual_fig) else 0

                row = dict(
                    tool=args.tool,
                    plot_type=args.plot_type,
                    size_label=args.size,
                    n_variants=n_variants,
                    replicate=rep,
                    wall_time_s=round(wall, 3),
                    peak_mem_mb=round(mem, 2),
                    out_file_kb=round(out_kb, 1),
                )
                writer.writerow(row)
                fh.flush()  # write incrementally in case of OOM on large runs
                print(f"  rep {rep}/{args.replicates}  "
                      f"time={wall:.2f}s  mem={mem:.0f}MB  fig={out_kb:.0f}KB")

            except Exception as e:
                print(f"  rep {rep} FAILED: {e}", file=sys.stderr)
                writer.writerow(dict(
                    tool=args.tool, plot_type=args.plot_type,
                    size_label=args.size, n_variants=n_variants,
                    replicate=rep, wall_time_s="ERROR",
                    peak_mem_mb="ERROR", out_file_kb="ERROR"
                ))
                fh.flush()


if __name__ == "__main__":
    main()

# #!/usr/bin/env python3
# """
# bench_python.py
# Benchmarks Python GWAS visualization tools: pycmplot, gwaslab, qmplot.
# 
# Usage:
#     python bench_python.py --tool pycmplot --input data/sumstats_1M.tsv \
#         --size 1M --replicates 5 --outdir results/
# 
# Writes one CSV row per replicate to results/<tool>_trimmed_bench_python.csv.
# """
# 
# import argparse
# import csv
# import gc
# import os
# import sys
# import time
# import tracemalloc
# from pathlib import Path
# 
# RESULT_COLS = [
#     "tool", "plot_type", "size_label", "n_variants",
#     "replicate", "wall_time_s", "peak_mem_mb", "out_file_kb"
# ]
# 
# 
# def _record(writer, row: dict):
#     writer.writerow({k: row.get(k, "") for k in RESULT_COLS})
# 
# 
# # ---------------------------------------------------------------------------
# # Individual tool wrappers
# # Each wrapper must:
# #   1. Load data from disk (include I/O in timing)
# #   2. Produce a PNG to out_path (or nearby with an appended suffix)
# #   3. Return nothing
# # ---------------------------------------------------------------------------
# 
# def run_pycmplot(input_path: str, out_path: str, plot_type: str = "manhattan"):
#     """
#     pycmplot benchmark.
# 
#     pycmplot appends an informative suffix to whatever stem is passed as
#     plot_title, so the file written to disk will differ from out_path.
#     benchmark_one() resolves the actual output via a directory snapshot.
#     """
#     import pycmplot
# 
#     sumstat = [input_path]
#     label = ["Trait"]
# 
#     out_parts = out_path.rsplit('/', 1)
#     out_dir  = out_parts[0]
#     out_file = out_parts[1]
# 
#     sumstats_info_dict = pycmplot.prep_pycmplot_input_info(
#         sum_stats=sumstat,
#         labels=label,
#         chrom="CHR",
#         pos="BP",
#         pcol="P",
#         snp="SNP"
#     )
# 
#     pycmplot_dict = pycmplot.get_sumstats_and_merged_sector_list(
#         sum_stats=sumstat,
#         labels=label,
#         file_info=sumstats_info_dict,
#         logp=True,
#         #table_out="pycmplot_bench", 
#         #trim_pval=0.01, 
#         #signif_threshold=1e-7        
#     )
# 
#     if plot_type == "manhattan":
#         pycmplot.plot_linear(
#             sumstats_loaded=pycmplot_dict["dfs"],
#             plot_title=out_file,
#             output_dir=out_dir,
#             output_format='png',
#             dpi=150,
#             logp=True,
#             #hits_table=hits_table, 
#             #chr_spacing=9e6,
#             #point_size=5, 
#             #linear_track_spacing=0.01,
#             #signif_lines=signif_lines,
#             #colors=['steelblue','orange'], 
#             #figsize=(10, 4), 
#             #trim_pval=0.01,            
#         )
#     elif plot_type == "circular":
#         pycmplot.plot_circular(
#             sumstats_loaded=pycmplot_dict["dfs"],
#             sector_sizes=pycmplot_dict["sectors"],
#             logp=True,
#             plot_title=out_file,
#             output_dir=out_dir,
#             dpi=150,
#             output_format='png',
#             pad=1,
#             r_min=50,
#             r_max=100,
#             #hits_table=hits_table,
#             #annotate="GENE",
#             #plot_title_size=9,
#             #label_col='top_gene',
#             #annotation_size=8,
#             #highlight=True,
#             #highlight_thresh=1e-7,
#             #highlight_line=True,
#             #signif_lines=signif_lines,
#             #colors=['steelblue','lightblue'],
#             #figsize=(10, 4),            
#         )
#     elif plot_type == "qq":
#         pycmplot.plot_qq_overlay(
#             pval_dict=pycmplot_dict["pvals"],
#             thin=True,
#             output_path=out_dir,
#             #thin_below=0.01
#         )
# 
# 
# def run_gwaslab(input_path: str, out_path: str, plot_type: str = "manhattan"):
#     """
#     gwaslab benchmark.
#     https://cloufield.github.io/gwaslab/
#     """
#     import pandas as pd
#     import gwaslab as gl
#     import matplotlib
#     matplotlib.use("Agg")
# 
#     df = pd.read_csv(input_path, sep="\t")
# 
#     mysumstats = gl.Sumstats(
#         df,
#         snpid="SNP",
#         chrom="CHR",
#         pos="BP",
#         p="P",
#         beta="BETA",
#         se="SE",
#         ea="A1",
#         nea="A2",
#     )
# 
#     if plot_type == "manhattan":
#         mysumstats.plot_mqq(
#             mode="m",
#             save=out_path,
#             save_kwargs={"dpi": 150},
#             verbose=False,
#         )
#     elif plot_type == "qq":
#         mysumstats.plot_mqq(
#             mode="qq",
#             save=out_path,
#             save_kwargs={"dpi": 150},
#             verbose=False,
#         )
# 
# 
# def run_qmplot(input_path: str, out_path: str, plot_type: str = "manhattan"):
#     """
#     qmplot benchmark.
#     https://github.com/ShujiaHuang/qmplot
#     """
#     import pandas as pd
#     import matplotlib
#     matplotlib.use("Agg")
#     import matplotlib.pyplot as plt
#     from qmplot import manhattanplot, qqplot
# 
#     df = pd.read_csv(input_path, sep="\t")
# 
#     fig, ax = plt.subplots(figsize=(12, 4))
# 
#     if plot_type in ("manhattan", "circular"):
#         manhattanplot(
#             data=df,
#             chrom="CHR",
#             pos="BP",
#             pv="P",
#             snp="SNP",
#             ax=ax,
#         )
#     elif plot_type == "qq":
#         qqplot(
#             data=df["P"],
#             ax=ax,
#         )
# 
#     fig.savefig(out_path, dpi=150, bbox_inches="tight")
#     plt.close(fig)
# 
# 
# # ---------------------------------------------------------------------------
# # Timing + memory harness
# # ---------------------------------------------------------------------------
# 
# TOOL_RUNNERS = {
#     "pycmplot": run_pycmplot,
#     #"gwaslab":  run_gwaslab,
#     #"qmplot":   run_qmplot,
# }
# 
# 
# def benchmark_one(tool: str, input_path: str, out_path: str, plot_type: str):
#     """
#     Run one timed, memory-tracked benchmark call.
# 
#     pycmplot appends an informative suffix to the stem passed as plot_title,
#     so the file written to disk differs from out_path. We snapshot the output
#     directory before and after the call to find whatever was actually written
#     and return that path so the caller can measure its size correctly.
# 
#     Tools that honour out_path exactly (gwaslab, qmplot) are unaffected —
#     the pre/post snapshot simply confirms out_path appeared and returns it.
# 
#     Returns
#     -------
#     tuple[float, float, str]
#         (wall_time_seconds, peak_memory_mb, actual_output_path)
#     """
#     runner = TOOL_RUNNERS[tool]
#     out_dir = os.path.dirname(os.path.abspath(out_path))
# 
#     # Snapshot PNGs present before the call
#     before = set(
#         os.path.join(out_dir, f)
#         for f in os.listdir(out_dir)
#         if f.lower().endswith(".png")
#     )
# 
#     # Force a full GC cycle before each run so prior allocations don't inflate
#     gc.collect()
# 
#     tracemalloc.start()
#     t0 = time.perf_counter()
# 
#     runner(input_path, out_path, plot_type)
# 
#     t1 = time.perf_counter()
#     _, peak_bytes = tracemalloc.get_traced_memory()
#     tracemalloc.stop()
# 
#     wall_time   = t1 - t0
#     peak_mem_mb = peak_bytes / 1024 / 1024
# 
#     # Resolve the actual output file written
#     after     = set(
#         os.path.join(out_dir, f)
#         for f in os.listdir(out_dir)
#         if f.lower().endswith(".png")
#     )
#     new_files = after - before
# 
#     if os.path.exists(out_path):
#         # Tool honoured out_path exactly (gwaslab, qmplot, pycmplot qq mode)
#         actual_out = out_path
#     elif len(new_files) == 1:
#         # Tool wrote a single new file with its own naming (pycmplot manhattan/circular)
#         actual_out = next(iter(new_files))
#     elif len(new_files) > 1:
#         # Multiple new files — prefer the one whose name starts with out_path's stem
#         stem = os.path.splitext(os.path.basename(out_path))[0]
#         matches = [f for f in new_files if os.path.basename(f).startswith(stem)]
#         actual_out = matches[0] if matches else next(iter(new_files))
#     else:
#         # Nothing written — return out_path so caller's exists() check gives False
#         actual_out = out_path
# 
#     return wall_time, peak_mem_mb, actual_out
# 
# 
# def main():
#     parser = argparse.ArgumentParser(description="Benchmark Python GWAS visualization tools")
#     parser.add_argument("--tool",        required=True, choices=list(TOOL_RUNNERS.keys()))
#     parser.add_argument("--input",       required=True, help="Path to sumstats TSV")
#     parser.add_argument("--size",        required=True, help="Dataset size label (e.g. 1M)")
#     parser.add_argument("--plot-type",   default="manhattan",
#                         choices=["manhattan", "circular", "qq"],
#                         help="Plot type to benchmark")
#     parser.add_argument("--replicates",  type=int, default=5)
#     parser.add_argument("--outdir",      default="results", help="Directory for CSV results")
#     parser.add_argument("--figdir",      default="figures", help="Directory for generated figures")
#     args = parser.parse_args()
# 
#     os.makedirs(args.outdir, exist_ok=True)
#     os.makedirs(args.figdir, exist_ok=True)
# 
#     #bench_res = args.tool + "_trimmed_bench_python.csv"
#     bench_res = args.tool + "_bench_python.csv"
#     csv_path = os.path.join(args.outdir, bench_res)
#     write_header = not os.path.exists(csv_path)
# 
#     # Count variants in input
#     import pandas as pd
#     n_variants = sum(1 for _ in open(args.input)) - 1  # subtract header
# 
#     print(f"\n[bench] tool={args.tool}  size={args.size}  "
#           f"n={n_variants:,}  plot={args.plot_type}  reps={args.replicates}")
# 
#     with open(csv_path, "a", newline="") as fh:
#         writer = csv.DictWriter(fh, fieldnames=RESULT_COLS)
#         if write_header:
#             writer.writeheader()
# 
#         for rep in range(1, args.replicates + 1):
#             out_fig = os.path.join(
#                 args.figdir,
#                 f"{args.tool}_{args.plot_type}_{args.size}_rep{rep}.png"
#             )
# 
#             try:
#                 wall, mem, actual_fig = benchmark_one(
#                     args.tool, args.input, out_fig, args.plot_type
#                 )
#                 out_kb = os.path.getsize(actual_fig) / 1024 if os.path.exists(actual_fig) else 0
# 
#                 row = dict(
#                     tool=args.tool,
#                     plot_type=args.plot_type,
#                     size_label=args.size,
#                     n_variants=n_variants,
#                     replicate=rep,
#                     wall_time_s=round(wall, 3),
#                     peak_mem_mb=round(mem, 2),
#                     out_file_kb=round(out_kb, 1),
#                 )
#                 writer.writerow(row)
#                 fh.flush()  # write incrementally in case of OOM on large runs
#                 print(f"  rep {rep}/{args.replicates}  "
#                       f"time={wall:.2f}s  mem={mem:.0f}MB  fig={out_kb:.0f}KB")
# 
#             except Exception as e:
#                 print(f"  rep {rep} FAILED: {e}", file=sys.stderr)
#                 writer.writerow(dict(
#                     tool=args.tool, plot_type=args.plot_type,
#                     size_label=args.size, n_variants=n_variants,
#                     replicate=rep, wall_time_s="ERROR",
#                     peak_mem_mb="ERROR", out_file_kb="ERROR"
#                 ))
#                 fh.flush()
# 
# 
# if __name__ == "__main__":
#     main()
