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

Mixed-build (liftover benchmark)
--------------------------------
Pass ``--liftover`` (or equivalently ``--builds hg19,hg19,hg38,hg38``) at the
1M size to produce a 4-trait set where two traits carry hg19 coordinates and
two carry hg38 coordinates.  Each trait's TSV gets its own single-build
scaffold *and* a ``BUILD`` column so downstream tools can dispatch the
liftover.  A companion ``sumstats_<SIZE>_<N>traits.builds.txt`` file lists
the per-trait builds in the same order as the manifest — this is what
``bench_python.py`` reads to pass ``build_list=`` to pycmplot's loader,
so the liftover happens inside the timed section.

Usage:
  # Standard multi-trait (all hg38, single scaffold)
  python generate_multi_sumstats.py --sizes 1M 2M --n-traits 3 --outdir data/

  # Liftover demo: 4 traits at 1M, two hg19 + two hg38
  python generate_multi_sumstats.py --sizes 1M --liftover --outdir data/

  # Scaling demo: multi-trait 50M / 100M
  python generate_multi_sumstats.py --sizes 50M 100M --n-traits 3 --outdir data/

  # All benchmark sizes, 3 traits
  python generate_multi_sumstats.py --all --n-traits 3 --outdir data/
"""

import argparse
import os
from typing import Iterable, List, Optional

import numpy as np
import pandas as pd

# We deliberately re-use the size + chrom-size tables from generate_sumstats.py
# so the two scripts can never drift.  If that import fails (e.g. the two
# scripts have been vendored to different paths), fall back to the local
# definitions below.
try:  # pragma: no cover — import-order defence
    from generate_sumstats import (  # type: ignore[import-not-found]
        HG19_CHROM_SIZES, HG38_CHROM_SIZES,
        DATASET_SIZES, stream_sumstats_to_tsv,
        _chrom_sizes_for_build, _distribute_variants, _chrom_frame,
    )
except Exception:  # pragma: no cover
    import sys as _sys
    _sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from generate_sumstats import (  # type: ignore[import-not-found]
        HG19_CHROM_SIZES, HG38_CHROM_SIZES,
        DATASET_SIZES, stream_sumstats_to_tsv,
        _chrom_sizes_for_build, _distribute_variants, _chrom_frame,
    )

# Streaming threshold for the per-trait writer.  Above this row count we
# emit per-chromosome and drop each chunk between writes; below it we
# build the whole trait DataFrame in memory (which is faster for small
# sizes).  20M keeps peak memory ~2 GB per trait for the largest chunk.
STREAM_THRESHOLD = 20_000_000


def _build_variant_scaffold(n_variants: int, seed: int = 42,
                            build: str = "hg38") -> pd.DataFrame:
    """Build the shared CHR / BP / SNP / A1 / A2 scaffold for all traits.

    Identical to the scaffold in generate_sumstats.py so single- and
    multi-trait files are directly comparable in benchmarks.
    """
    rng = np.random.default_rng(seed)
    chrom_sizes = _chrom_sizes_for_build(build)
    counts = _distribute_variants(n_variants, chrom_sizes)

    rows = []
    snp_offset = 0
    for chrom, count in counts.items():
        # Position + alleles only; per-trait BETA/SE/P are added later.
        df = _chrom_frame(
            chrom, count, snp_offset,
            max_bp=chrom_sizes[chrom],
            rng=rng, build_label=None,
        )
        rows.append(df[["CHR", "SNP", "BP", "A1", "A2"]])
        snp_offset += count

    return pd.concat(rows, ignore_index=True)


def _add_trait_columns(
    scaffold: pd.DataFrame,
    n_signals: int = 25,
    seed: int = 1,
) -> pd.DataFrame:
    """Add BETA, SE, P columns to a scaffold copy using an independent seed.

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


def _write_trait_stream(out_path: str, n_variants: int, seed: int,
                        build: Optional[str], n_signals: int) -> None:
    """Write one trait TSV without ever materialising the full frame.

    Delegates to :func:`generate_sumstats.stream_sumstats_to_tsv` so the
    large-scale path shares its per-chromosome logic (and the same signal
    injection semantics) with the single-trait generator.
    """
    stream_sumstats_to_tsv(
        out_path=out_path,
        n_variants=n_variants,
        n_signals=n_signals,
        seed=seed,
        build=build,
    )


def generate_multi_sumstats(
    n_variants: int,
    n_traits: int,
    size_label: str,
    outdir: str,
    scaffold_seed: int = 42,
    trait_seed_offset: int = 100,
    n_signals: int = 25,
    force: bool = False,
    builds: Optional[List[str]] = None,
):
    """Generate n_traits sumstats TSV files sharing the same variant scaffold.

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
    builds : list[str] or None
        Per-trait build labels (``"hg19"`` / ``"hg38"``).  ``None`` means
        "all hg38, single shared scaffold, no BUILD column" — the original
        behaviour.  When set, each trait's scaffold is regenerated under
        its declared build (so hg19 traits get real hg19 coordinates) and
        the trait TSVs carry a ``BUILD`` column.  A companion
        ``sumstats_<SIZE>_<N>traits.builds.txt`` file is written alongside
        the manifest, listing one build per line in trait order.

    Returns
    -------
    list[str]
        Paths to the generated trait TSV files (in trait order).
    """
    os.makedirs(outdir, exist_ok=True)
    if builds is not None:
        if len(builds) != n_traits:
            raise ValueError(
                f"--builds gave {len(builds)} entries but --n-traits is "
                f"{n_traits}; they must match one-to-one."
            )
        # Filename suffix tags the mixed-build set so it doesn't collide
        # with the plain hg38 output at the same size.
        suffix = "_mixedbuild"
    else:
        suffix = ""

    trait_paths = [
        os.path.join(outdir,
                     f"sumstats_{size_label}{suffix}_trait{k+1}.tsv")
        for k in range(n_traits)
    ]
    manifest_path = os.path.join(
        outdir, f"sumstats_{size_label}{suffix}_{n_traits}traits.manifest",
    )
    builds_path = os.path.join(
        outdir, f"sumstats_{size_label}{suffix}_{n_traits}traits.builds.txt",
    )

    if not force and all(os.path.exists(p) for p in trait_paths):
        print(f"[skip] All {n_traits} trait files for {size_label}"
              f"{suffix} already exist")
        # Re-emit manifest + builds anyway so downstream scripts don't
        # depend on the caller having them from a previous run.
        _write_manifest_and_builds(manifest_path, builds_path,
                                   trait_paths, builds)
        return trait_paths

    # ---- Streaming (large-N) path ---------------------------------
    if n_variants > STREAM_THRESHOLD:
        print(f"Streaming mode: peak memory bounded by chr1 per trait "
              f"(~{n_variants // 22 * 40 / 1e9:.1f} GB per chunk).")
        for k in range(n_traits):
            trait_seed = trait_seed_offset + k
            trait_build = builds[k] if builds else None
            label = f"Trait{k+1}"
            out_path = trait_paths[k]
            if not force and os.path.exists(out_path):
                print(f"  [skip] {os.path.basename(out_path)} already exists")
                continue
            print(f"  Streaming {label} "
                  f"(seed={trait_seed}, build={trait_build or 'hg38'})...")
            _write_trait_stream(
                out_path=out_path, n_variants=n_variants,
                seed=trait_seed, build=trait_build,
                n_signals=n_signals,
            )
            size_mb = os.path.getsize(out_path) / 1e6
            print(f"    -> {os.path.basename(out_path)}  ({size_mb:.1f} MB)")
        _write_manifest_and_builds(manifest_path, builds_path,
                                   trait_paths, builds)
        return trait_paths

    # ---- In-memory (small-N) path --------------------------------
    if builds is None:
        # Shared scaffold, no BUILD column — original behaviour.
        print(f"Building variant scaffold: "
              f"{n_variants:,} variants (seed={scaffold_seed})...")
        scaffold = _build_variant_scaffold(n_variants, seed=scaffold_seed)

        for k in range(n_traits):
            trait_seed = trait_seed_offset + k
            label = f"Trait{k+1}"
            out_path = trait_paths[k]
            if not force and os.path.exists(out_path):
                print(f"  [skip] {os.path.basename(out_path)} already exists")
                continue
            print(f"  Generating {label} (seed={trait_seed})...")
            df = _add_trait_columns(scaffold, n_signals=n_signals,
                                    seed=trait_seed)
            df.to_csv(out_path, sep="\t", index=False)
            size_mb = os.path.getsize(out_path) / 1e6
            print(f"    -> {os.path.basename(out_path)}  ({size_mb:.1f} MB)")
    else:
        # Mixed-build path.  Because each build has its own chromosome
        # sizes we don't attempt to share a single scaffold — each
        # trait gets a freshly-built scaffold under its own build, with
        # a per-build scaffold seed so hg19 traits share positions and
        # hg38 traits share positions (a realistic pattern when two
        # cohorts were genotyped on the same panel).
        by_build_scaffold: dict[str, pd.DataFrame] = {}
        for k in range(n_traits):
            trait_build = builds[k]
            if trait_build not in by_build_scaffold:
                print(f"Building {trait_build} scaffold: {n_variants:,} variants ...")
                by_build_scaffold[trait_build] = _build_variant_scaffold(
                    n_variants,
                    seed=scaffold_seed + hash(trait_build) % 10_000,
                    build=trait_build,
                )
            scaffold = by_build_scaffold[trait_build]
            trait_seed = trait_seed_offset + k
            label = f"Trait{k+1}"
            out_path = trait_paths[k]
            if not force and os.path.exists(out_path):
                print(f"  [skip] {os.path.basename(out_path)} already exists")
                continue
            print(f"  Generating {label} "
                  f"(seed={trait_seed}, build={trait_build})...")
            df = _add_trait_columns(scaffold, n_signals=n_signals,
                                    seed=trait_seed)
            df["BUILD"] = trait_build
            df.to_csv(out_path, sep="\t", index=False)
            size_mb = os.path.getsize(out_path) / 1e6
            print(f"    -> {os.path.basename(out_path)}  ({size_mb:.1f} MB)")

    _write_manifest_and_builds(manifest_path, builds_path,
                               trait_paths, builds)
    return trait_paths


def _write_manifest_and_builds(manifest_path: str, builds_path: str,
                               trait_paths: Iterable[str],
                               builds: Optional[List[str]]) -> None:
    """Write the manifest (always) and the .builds.txt sidecar (mixed only).

    Split out so both the streaming and in-memory paths share the same
    logic — the two file formats are tiny but downstream scripts depend
    on them, so a single writer is safer than duplicating a couple of
    ``open()`` blocks in each path.
    """
    with open(manifest_path, "w") as fh:
        for p in trait_paths:
            fh.write(os.path.abspath(p) + "\n")
    print(f"  Manifest: {manifest_path}")

    if builds is not None:
        with open(builds_path, "w") as fh:
            for b in builds:
                fh.write(f"{b}\n")
        print(f"  Builds:   {builds_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate multi-trait GWAS sumstats for benchmarking",
        formatter_class=argparse.RawDescriptionHelpFormatter,
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
        help="Generate all benchmark sizes (500K, 1M, 2M, 5M, 10M, 25M, 50M, 100M)"
    )
    parser.add_argument(
        "--builds", type=str, default=None,
        help=(
            "Comma-separated per-trait builds (hg19|hg38) for the "
            "liftover benchmark.  E.g. --builds hg19,hg19,hg38,hg38.  "
            "Length must match --n-traits.  Emits a BUILD column and a "
            "companion .builds.txt sidecar."
        ),
    )
    parser.add_argument(
        "--liftover", action="store_true",
        help=(
            "Convenience flag for the manuscript's mixed-build "
            "benchmark: shorthand for `--sizes 1M --n-traits 4 "
            "--builds hg19,hg19,hg38,hg38`."
        ),
    )
    parser.add_argument("--outdir", default="data")
    parser.add_argument(
        "--force", action="store_true",
        help="Overwrite existing files"
    )
    args = parser.parse_args()

    if args.liftover:
        # Override any user-supplied conflicting flags with a warning.
        if args.sizes != ["1M"]:
            print(f"[warn] --liftover overrides --sizes; using ['1M'].")
        if args.n_traits != [3]:
            print(f"[warn] --liftover overrides --n-traits; using [4].")
        if args.builds:
            print(f"[warn] --liftover overrides --builds; using "
                  f"hg19,hg19,hg38,hg38.")
        args.sizes = ["1M"]
        args.n_traits = [4]
        args.builds = "hg19,hg19,hg38,hg38"

    builds: Optional[List[str]] = None
    if args.builds:
        builds = [b.strip() for b in args.builds.split(",") if b.strip()]
        bad = [b for b in builds if b.lower() not in ("hg19", "hg38")]
        if bad:
            parser.error(f"--builds contains unsupported entries: {bad}")

    sizes = list(DATASET_SIZES.keys()) if args.all else args.sizes

    for size_label in sizes:
        n_variants = DATASET_SIZES[size_label]
        for n_traits in args.n_traits:
            if builds is not None and len(builds) != n_traits:
                parser.error(
                    f"--builds has {len(builds)} entries but "
                    f"--n-traits is {n_traits}; they must match."
                )
            hdr = f"{size_label} | {n_traits} traits ({n_variants:,} variants each)"
            if builds:
                hdr += f" | builds={','.join(builds)}"
            print(f"\n=== {hdr} ===")
            generate_multi_sumstats(
                n_variants=n_variants,
                n_traits=n_traits,
                size_label=size_label,
                outdir=args.outdir,
                force=args.force,
                builds=builds,
            )

    print("\nDone.")
    print("Pass the .manifest file to bench scripts via --manifest <path>")
    if builds:
        print("For mixed-build runs, the sibling .builds.txt is auto-consumed "
              "by bench_python.py for the liftover benchmark.")


if __name__ == "__main__":
    main()
