#!/usr/bin/env python3
"""
generate_sumstats.py
Generates synthetic GWAS summary statistics for benchmarking.

Usage:
    python generate_sumstats.py --n 1000000 --out data/sumstats_1M.tsv
    python generate_sumstats.py --all --outdir data/
    python generate_sumstats.py --n 100_000_000 --out data/sumstats_100M.tsv
    python generate_sumstats.py --n 1_000_000 --build hg19 --out data/sumstats_1M_hg19.tsv

Notes on scale
--------------
The 50M and 100M paths write **one chromosome at a time**: build a per-chrom
DataFrame, append it to the output TSV, drop it.  Peak memory is therefore
bounded by chr1's share (~10 % of total) rather than the full-genome frame,
which lets a 64 GB node handle 100M variants without OOMing.  This matters
for the extended benchmark:

    500K  →  ~40 MB TSV
    100M  →  ~8-9 GB TSV, ~3-4 GB peak Python memory

Runtime for --all now takes ~5-10 min per size at 100M, so we skip files that
already exist unless --force is passed.
"""

import numpy as np
import pandas as pd
import argparse
import os

# hg38 (GRCh38) chromosome sizes in bp (chr1–22).  These are the reference
# coordinates the sumstats live in unless --build hg19 is specified.
HG38_CHROM_SIZES = {
    1: 248956422, 2: 242193529, 3: 198295559, 4: 190214555,
    5: 181538259, 6: 170805979, 7: 159345973, 8: 145138636,
    9: 138394717, 10: 133797422, 11: 135086622, 12: 133275309,
    13: 114364328, 14: 107043718, 15: 101991189, 16: 90338345,
    17: 83257441,  18: 80373285,  19: 58617616,  20: 64444167,
    21: 46709983,  22: 50818468
}

# hg19 (GRCh37) chromosome sizes in bp (chr1–22).  Used when --build hg19
# is set so the generated positions are actually plausible hg19 coordinates
# (roughly the same lengths, but distinct enough that a liftover pass will
# produce genuinely different values).
HG19_CHROM_SIZES = {
    1: 249250621, 2: 243199373, 3: 198022430, 4: 191154276,
    5: 180915260, 6: 171115067, 7: 159138663, 8: 146364022,
    9: 141213431, 10: 135534747, 11: 135006516, 12: 133851895,
    13: 115169878, 14: 107349540, 15: 102531392, 16: 90354753,
    17: 81195210,  18: 78077248,  19: 59128983,  20: 63025520,
    21: 48129895,  22: 51304566
}

# Backwards compatible alias — earlier revisions of this script exposed
# ``CHROM_SIZES`` at module scope, and downstream code (bench_python.py,
# various analysis notebooks) imports it directly.
CHROM_SIZES = HG38_CHROM_SIZES

DATASET_SIZES = {
    "500K":  500_000,
    "1M":  1_000_000,
    "2M":  2_000_000,
    "5M":  5_000_000,
    "10M": 10_000_000,
    # -- Added for the 100M-variant scaling demo --------------------
    "25M":  25_000_000,
    "50M":  50_000_000,
    #"100M": 100_000_000,
}


def _chrom_sizes_for_build(build: str) -> dict[int, int]:
    b = (build or "hg38").lower()
    if b in ("hg19", "grch37"):
        return HG19_CHROM_SIZES
    if b in ("hg38", "grch38"):
        return HG38_CHROM_SIZES
    raise ValueError(f"Unsupported build {build!r}; use 'hg19' or 'hg38'.")


def _distribute_variants(n_variants: int,
                         chrom_sizes: dict[int, int]) -> dict[int, int]:
    """Split n_variants across chr1..22 proportional to length.

    Guarantees at least 1 variant per chromosome and that the counts
    sum to exactly ``n_variants`` (the tail chromosome absorbs any
    rounding remainder).
    """
    total_bp = sum(chrom_sizes.values())
    chroms = list(chrom_sizes.keys())
    counts: dict[int, int] = {}
    remaining = n_variants
    for chrom in chroms[:-1]:
        count = int(n_variants * chrom_sizes[chrom] / total_bp)
        counts[chrom] = max(count, 1)
        remaining -= counts[chrom]
    counts[chroms[-1]] = max(remaining, 1)
    return counts


def _chrom_frame(chrom: int, count: int, snp_offset: int,
                 max_bp: int, rng: np.random.Generator,
                 build_label: str | None) -> pd.DataFrame:
    """Materialise the rows for a single chromosome.

    Kept as a function so both the whole-genome path (small sizes,
    for API compatibility) and the streaming path (large sizes) share
    exactly the same row shape.
    """
    positions = np.sort(rng.integers(10_000, max_bp - 10_000, size=count))
    df = pd.DataFrame({
        "CHR":  chrom,
        "SNP":  [f"rs{snp_offset + i:09d}" for i in range(count)],
        "BP":   positions,
        "A1":   rng.choice(["A", "C", "G", "T"], size=count),
        "A2":   rng.choice(["A", "C", "G", "T"], size=count),
        "BETA": rng.normal(0, 0.05, size=count),
        "SE":   np.abs(rng.normal(0.02, 0.005, size=count)),
        "P":    rng.uniform(0, 1, size=count),
    })
    if build_label:
        df["BUILD"] = build_label
    return df


def generate_sumstats(n_variants: int, n_signals: int = 30, seed: int = 42,
                      build: str | None = None) -> pd.DataFrame:
    """Generate synthetic GWAS summary statistics (single-shot, in-memory).

    Parameters
    ----------
    n_variants : int
        Total number of variants.
    n_signals : int
        Number of simulated association signals (p < 5e-8).
    seed : int
        Random seed for reproducibility.
    build : {"hg19", "hg38", None}
        When set, an additional ``BUILD`` column is emitted and the
        chromosome-size table is picked accordingly.  ``None``
        preserves the pre-liftover behaviour of the script (no
        BUILD column, hg38 coordinates).

    Returns
    -------
    pd.DataFrame
        DataFrame with columns CHR, SNP, BP, A1, A2, BETA, SE, P
        (and BUILD when ``build`` is passed).

    Notes
    -----
    This path holds every row in memory.  For 50M+ variants,
    :func:`stream_sumstats_to_tsv` is used by ``main()`` instead;
    this function is still exposed for programmatic use and small
    sizes (kept unchanged for API compatibility with callers that
    imported it before the 100M extension).
    """
    rng = np.random.default_rng(seed)
    chrom_sizes = _chrom_sizes_for_build(build) if build else HG38_CHROM_SIZES
    counts = _distribute_variants(n_variants, chrom_sizes)

    rows = []
    snp_offset = 0
    for chrom, count in counts.items():
        rows.append(_chrom_frame(
            chrom, count, snp_offset,
            max_bp=chrom_sizes[chrom],
            rng=rng, build_label=build,
        ))
        snp_offset += count

    df = pd.concat(rows, ignore_index=True)

    # Inject association signals
    signal_indices = rng.choice(len(df), size=n_signals, replace=False)
    for idx in signal_indices:
        df.loc[idx, "P"] = 10 ** rng.uniform(-50, -8)

    # Ensure no p=0 or p>1
    df["P"] = df["P"].clip(1e-300, 1.0)

    return df


def stream_sumstats_to_tsv(out_path: str, n_variants: int,
                           n_signals: int = 30, seed: int = 42,
                           build: str | None = None,
                           chunk_bytes_hint: int = 1 << 30) -> None:
    """Write a synthetic sumstats TSV chromosome-by-chromosome.

    This is the streaming counterpart to :func:`generate_sumstats`,
    used automatically by ``main()`` for the 25M / 50M / 100M sizes
    where the full-genome frame would exhaust node memory.  Peak
    Python memory is bounded by the largest single-chromosome frame
    (chr1 ≈ 10 % of total) rather than the whole file.

    Signals are injected on the fly: we sample ``n_signals`` global
    row indices at the start, then per chromosome check which of them
    fall inside the current row range and rewrite those P values
    before the chunk is flushed.
    """
    rng = np.random.default_rng(seed)
    chrom_sizes = _chrom_sizes_for_build(build) if build else HG38_CHROM_SIZES
    counts = _distribute_variants(n_variants, chrom_sizes)
    total = sum(counts.values())

    # Pre-sample signal row indices in the *global* row space so their
    # distribution matches the in-memory ``generate_sumstats`` path.
    signal_global_idx = np.sort(
        rng.choice(total, size=n_signals, replace=False)
    )
    signal_p = 10 ** rng.uniform(-50, -8, size=n_signals)

    snp_offset = 0
    header_written = False
    with open(out_path, "w") as fh:
        for chrom, count in counts.items():
            df = _chrom_frame(
                chrom, count, snp_offset,
                max_bp=chrom_sizes[chrom],
                rng=rng, build_label=build,
            )
            # Apply signals whose global index lies in this chunk.
            lo, hi = snp_offset, snp_offset + count
            in_chunk = (signal_global_idx >= lo) & (signal_global_idx < hi)
            if in_chunk.any():
                local_positions = signal_global_idx[in_chunk] - lo
                df.iloc[local_positions,
                        df.columns.get_loc("P")] = signal_p[in_chunk]
            df["P"] = df["P"].clip(1e-300, 1.0)

            df.to_csv(fh, sep="\t", index=False, header=not header_written)
            header_written = True
            snp_offset += count
            # Free frame between chromosomes so peak stays flat.
            del df


def main():
    parser = argparse.ArgumentParser(description="Generate synthetic GWAS sumstats for benchmarking")
    parser.add_argument("--n", type=int, help="Number of variants")
    parser.add_argument("--out", type=str, help="Output file path")
    parser.add_argument("--all", action="store_true", help="Generate all benchmark sizes")
    parser.add_argument("--outdir", type=str, default="data", help="Output directory (used with --all)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument(
        "--build", choices=["hg19", "hg38"], default=None,
        help=(
            "Emit a BUILD column and use build-specific chromosome "
            "sizes.  Used by the multi-trait liftover benchmark to "
            "generate genuinely mixed-coordinate inputs."
        ),
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Overwrite existing files (default: skip if present).",
    )
    parser.add_argument(
        "--stream-threshold", type=int, default=20_000_000,
        help=(
            "Row count above which the streaming (per-chromosome) "
            "writer is used instead of the in-memory generator.  "
            "Default 20M, which keeps peak memory under ~2 GB for "
            "any single-chromosome chunk."
        ),
    )
    args = parser.parse_args()

    def _emit_one(n_variants: int, out_path: str, build: str | None) -> None:
        if os.path.exists(out_path) and not args.force:
            print(f"[skip] {out_path} already exists (use --force to overwrite)")
            return
        os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
        stream = n_variants > args.stream_threshold
        tag = "streaming" if stream else "in-memory"
        print(f"Generating {n_variants:,} variants "
              f"(build={build or 'hg38'}, mode={tag}) -> {out_path} ...")
        if stream:
            stream_sumstats_to_tsv(out_path, n_variants,
                                   seed=args.seed, build=build)
        else:
            df = generate_sumstats(n_variants, seed=args.seed, build=build)
            df.to_csv(out_path, sep="\t", index=False)
            del df
        size_mb = os.path.getsize(out_path) / 1e6
        print(f"  -> {out_path}  ({size_mb:.1f} MB)")

    if args.all:
        os.makedirs(args.outdir, exist_ok=True)
        for label, n in DATASET_SIZES.items():
            suffix = f"_{args.build}" if args.build else ""
            out_path = os.path.join(args.outdir,
                                    f"sumstats_{label}{suffix}.tsv")
            _emit_one(n, out_path, args.build)
    elif args.n and args.out:
        _emit_one(args.n, args.out, args.build)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
