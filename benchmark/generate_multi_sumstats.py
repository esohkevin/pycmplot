#!/usr/bin/env python3
"""
generate_multi_sumstats.py
Generates multi-trait GWAS summary statistics for benchmarking multi-track plots.

All traits share the same variant scaffold (CHR, BP, SNP, A1, A2) — mimicking
a real study where multiple traits are measured in the same cohort. Each trait
gets independent p-values and effect sizes drawn from different seeds.

Outputs per (size, n_traits):
  - N separate TSV files  (pycmplot format): sumstats_<SIZE>_trait<K>.tsv
  - 1 manifest file       (lists TSV paths): sumstats_<SIZE>_<N>traits.manifest
  The column-wise CMplot merge is done inside bench_r.R as part of its timed
  section, since that merge is required workflow for CMplot users.

Usage:
  # Generate 3-trait files at 1M and 2M variants
  python generate_multi_sumstats.py --sizes 1M 2M --n-traits 3 --outdir data/

  # Generate 2-, 3-, and 5-trait files at 1M
  python generate_multi_sumstats.py --sizes 1M --n-traits 2 3 5 --outdir data/

  # Generate all benchmark sizes, 3 traits
  python generate_multi_sumstats.py --all --n-traits 3 --outdir data/
"""

import argparse
import os

import numpy as np
import pandas as pd

# Approximate hg38 chromosome sizes in bp (chr1–22) — same as generate_sumstats.py
CHROM_SIZES = {
    1: 248956422, 2: 242193529, 3: 198295559, 4: 190214555,
    5: 181538259, 6: 170805979, 7: 159345973, 8: 145138636,
    9: 138394717, 10: 133797422, 11: 135086622, 12: 133275309,
    13: 114364328, 14: 107043718, 15: 101991189, 16: 90338345,
    17: 83257441,  18: 80373285,  19: 58617616,  20: 64444167,
    21: 46709983,  22: 50818468
}

DATASET_SIZES = {
    "500K": 500_000,
    "1M":   1_000_000,
    "2M":   2_000_000,
    "5M":   5_000_000,
}


def _build_variant_scaffold(n_variants: int, seed: int = 42) -> pd.DataFrame:
    """
    Build the shared CHR / BP / SNP / A1 / A2 scaffold for all traits.
    Identical to the scaffold in generate_sumstats.py so single- and
    multi-trait files are directly comparable in benchmarks.
    """
    rng = np.random.default_rng(seed)
    total_bp = sum(CHROM_SIZES.values())
    chroms = list(CHROM_SIZES.keys())

    chrom_counts: dict[int, int] = {}
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
        rows.append(pd.DataFrame({
            "CHR": chrom,
            "SNP": [f"rs{snp_offset + i:09d}" for i in range(count)],
            "BP":  positions,
            "A1":  rng.choice(["A", "C", "G", "T"], size=count),
            "A2":  rng.choice(["A", "C", "G", "T"], size=count),
        }))
        snp_offset += count

    return pd.concat(rows, ignore_index=True)


def _add_trait_columns(
    scaffold: pd.DataFrame,
    n_signals: int = 25,
    seed: int = 1,
) -> pd.DataFrame:
    """
    Add BETA, SE, P columns to a scaffold copy using an independent seed.
    Each trait gets its own signals at randomly chosen loci.
    """
    rng = np.random.default_rng(seed)
    n = len(scaffold)
    df = scaffold.copy()
    df["BETA"] = rng.normal(0, 0.05, size=n)
    df["SE"]   = np.abs(rng.normal(0.02, 0.005, size=n))
    df["P"]    = rng.uniform(0, 1, size=n).clip(1e-300, 1.0)

    # Inject association signals
    signal_idx = rng.choice(n, size=n_signals, replace=False)
    for idx in signal_idx:
        df.loc[idx, "P"] = 10 ** rng.uniform(-50, -8)

    return df


def generate_multi_sumstats(
    n_variants: int,
    n_traits: int,
    size_label: str,
    outdir: str,
    scaffold_seed: int = 42,
    trait_seed_offset: int = 100,
    n_signals: int = 25,
    force: bool = False,
):
    """
    Generate n_traits sumstats TSV files sharing the same variant scaffold.

    Parameters
    ----------
    n_variants : int
    n_traits : int
    size_label : str
        Label used in filenames, e.g. "1M".
    outdir : str
        Output directory.
    scaffold_seed : int
        Seed for variant position / allele generation (shared across traits).
    trait_seed_offset : int
        Trait k gets seed = trait_seed_offset + k, ensuring independence.
    n_signals : int
        Simulated association signals per trait.
    force : bool
        Overwrite existing files if True.

    Returns
    -------
    list[str]
        Paths to the generated trait TSV files (in trait order).
    """
    os.makedirs(outdir, exist_ok=True)

    # Check if all files already exist
    trait_paths = [
        os.path.join(outdir, f"sumstats_{size_label}_trait{k+1}.tsv")
        for k in range(n_traits)
    ]
    manifest_path = os.path.join(outdir, f"sumstats_{size_label}_{n_traits}traits.manifest")

    if not force and all(os.path.exists(p) for p in trait_paths):
        print(f"[skip] All {n_traits} trait files for {size_label} already exist")
        return trait_paths

    print(f"Building variant scaffold: {n_variants:,} variants (seed={scaffold_seed})...")
    scaffold = _build_variant_scaffold(n_variants, seed=scaffold_seed)

    for k in range(n_traits):
        trait_seed = trait_seed_offset + k
        label = f"Trait{k+1}"
        out_path = trait_paths[k]

        if not force and os.path.exists(out_path):
            print(f"  [skip] {os.path.basename(out_path)} already exists")
            continue

        print(f"  Generating {label} (seed={trait_seed})...")
        df = _add_trait_columns(scaffold, n_signals=n_signals, seed=trait_seed)
        df.to_csv(out_path, sep="\t", index=False)
        size_mb = os.path.getsize(out_path) / 1e6
        print(f"    -> {os.path.basename(out_path)}  ({size_mb:.1f} MB)")

    # Write manifest — one absolute path per line
    with open(manifest_path, "w") as fh:
        for p in trait_paths:
            fh.write(os.path.abspath(p) + "\n")
    print(f"  Manifest: {manifest_path}")

    return trait_paths


def main():
    parser = argparse.ArgumentParser(
        description="Generate multi-trait GWAS sumstats for benchmarking"
    )
    parser.add_argument(
        "--sizes", nargs="+", default=["1M"],
        choices=list(DATASET_SIZES.keys()),
        help="Dataset sizes to generate (default: 1M)"
    )
    parser.add_argument(
        "--n-traits", nargs="+", type=int, default=[3],
        help="Number of traits per dataset (default: 3). Multiple values allowed."
    )
    parser.add_argument(
        "--all", action="store_true",
        help="Generate all benchmark sizes (500K, 1M, 2M, 5M)"
    )
    parser.add_argument("--outdir", default="data")
    parser.add_argument(
        "--force", action="store_true",
        help="Overwrite existing files"
    )
    args = parser.parse_args()

    sizes = list(DATASET_SIZES.keys()) if args.all else args.sizes

    for size_label in sizes:
        n_variants = DATASET_SIZES[size_label]
        for n_traits in args.n_traits:
            print(f"\n=== {size_label} | {n_traits} traits ({n_variants:,} variants each) ===")
            generate_multi_sumstats(
                n_variants=n_variants,
                n_traits=n_traits,
                size_label=size_label,
                outdir=args.outdir,
                force=args.force,
            )

    print("\nDone.")
    print("Pass the .manifest file to bench scripts via --manifest <path>")


if __name__ == "__main__":
    main()
