#!/usr/bin/env python3
"""
generate_sumstats.py
Generates synthetic GWAS summary statistics for benchmarking.

Usage:
    python generate_sumstats.py --n 1000000 --out data/sumstats_1M.tsv
    python generate_sumstats.py --all --outdir data/
"""

import numpy as np
import pandas as pd
import argparse
import os

# Approximate hg38 chromosome sizes in bp (chr1–22)
CHROM_SIZES = {
    1: 248956422, 2: 242193529, 3: 198295559, 4: 190214555,
    5: 181538259, 6: 170805979, 7: 159345973, 8: 145138636,
    9: 138394717, 10: 133797422, 11: 135086622, 12: 133275309,
    13: 114364328, 14: 107043718, 15: 101991189, 16: 90338345,
    17: 83257441,  18: 80373285,  19: 58617616,  20: 64444167,
    21: 46709983,  22: 50818468
}

DATASET_SIZES = {
    "500K":  500_000,
    "1M":  1_000_000,
    "2M":  2_000_000,
    "5M":  5_000_000,
    "10M": 10_000_000,
}


def generate_sumstats(n_variants: int, n_signals: int = 30, seed: int = 42) -> pd.DataFrame:
    """
    Generate synthetic GWAS summary statistics.

    Parameters
    ----------
    n_variants : int
        Total number of variants.
    n_signals : int
        Number of simulated association signals (p < 5e-8).
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns CHR, SNP, BP, P, BETA, SE, A1, A2.
    """
    rng = np.random.default_rng(seed)
    total_bp = sum(CHROM_SIZES.values())
    chroms = list(CHROM_SIZES.keys())

    # Distribute variants proportionally by chromosome length
    chrom_counts = {}
    remaining = n_variants
    for chrom in chroms[:-1]:
        count = int(n_variants * CHROM_SIZES[chrom] / total_bp)
        chrom_counts[chrom] = max(count, 1)
        remaining -= chrom_counts[chrom]
    chrom_counts[22] = max(remaining, 1)

    rows = []
    snp_offset = 0
    for chrom in chroms:
        count = chrom_counts[chrom]
        positions = np.sort(
            rng.integers(10_000, CHROM_SIZES[chrom] - 10_000, size=count)
        )
        pvals = rng.uniform(0, 1, size=count)
        betas = rng.normal(0, 0.05, size=count)
        ses = np.abs(rng.normal(0.02, 0.005, size=count))

        rows.append(pd.DataFrame({
            "CHR":  chrom,
            "SNP":  [f"rs{snp_offset + i:09d}" for i in range(count)],
            "BP":   positions,
            "A1":   rng.choice(["A", "C", "G", "T"], size=count),
            "A2":   rng.choice(["A", "C", "G", "T"], size=count),
            "BETA": betas,
            "SE":   ses,
            "P":    pvals,
        }))
        snp_offset += count

    df = pd.concat(rows, ignore_index=True)

    # Inject association signals
    signal_indices = rng.choice(len(df), size=n_signals, replace=False)
    for idx in signal_indices:
        # Cluster nearby variants around each signal for realism
        df.loc[idx, "P"] = 10 ** rng.uniform(-50, -8)

    # Ensure no p=0 or p>1
    df["P"] = df["P"].clip(1e-300, 1.0)

    return df


def main():
    parser = argparse.ArgumentParser(description="Generate synthetic GWAS sumstats for benchmarking")
    parser.add_argument("--n", type=int, help="Number of variants")
    parser.add_argument("--out", type=str, help="Output file path")
    parser.add_argument("--all", action="store_true", help="Generate all benchmark sizes")
    parser.add_argument("--outdir", type=str, default="data", help="Output directory (used with --all)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    if args.all:
        os.makedirs(args.outdir, exist_ok=True)
        for label, n in DATASET_SIZES.items():
            out_path = os.path.join(args.outdir, f"sumstats_{label}.tsv")
            if os.path.exists(out_path):
                print(f"[skip] {out_path} already exists")
                continue
            print(f"Generating {label} ({n:,} variants)...")
            df = generate_sumstats(n, seed=args.seed)
            df.to_csv(out_path, sep="\t", index=False)
            print(f"  -> {out_path}  ({os.path.getsize(out_path) / 1e6:.1f} MB)")
    elif args.n and args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        df = generate_sumstats(args.n, seed=args.seed)
        df.to_csv(args.out, sep="\t", index=False)
        print(f"Generated {args.n:,} variants -> {args.out}")
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
