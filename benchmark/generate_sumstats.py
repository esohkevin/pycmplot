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
from typing import Optional

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


# ---------------------------------------------------------------------------
# Target spike-in loci — well-known GWAS peaks in both coordinate systems.
# ---------------------------------------------------------------------------
# The synthetic sumstats should look biologically plausible, so instead of
# planting random genome-wide-significant hits everywhere we plant *peaks*
# at these six loci — one lead SNP per locus at the exact target position
# with p = min_p_value, plus ~20 supporting SNPs distributed within
# ±WINDOW_HALF_SIZE around the lead.  Supporting p-values decay from the
# lead as ``p(d) = min_p · 10^(|d|/tau)`` with ``tau = 100 kb``, so at
# 100 kb from the lead the p-value is 10× weaker, at 200 kb it's 100×
# weaker.  This produces the tall lead + diffuse cloud shape typical of
# a real GWAS Manhattan peak.
#
# The hg19 and hg38 lists are pre-paired — each row has the same locus 
# identity across builds so liftover round-trips can be validated by 
# comparing (build, position) tuples.
HG19_TARGET_SPIKES: list[tuple[str, int, float]] = [
    ("1",  172_058_886, 7e-12),   # rs186897198
    ("2",  36_733_328,  3e-08),   # rs2030645 - hg19-specific
    ("7",  18_786_817,  7e-23),   # rs727851
    ("10", 31_127_166,  7e-08),   # rs12413361
    ("11", 12_879_123,  5e-09),   # rs546512774
    ("11", 2_802_090,   9e-08),   # rs234886
    ("15", 22_791_431,  1e-08),   # rs6606792 - hg19-specific
    ("17", 65_854_602,  2e-12),   # rs55931203
]

HG38_TARGET_SPIKES: list[tuple[str, int, float]] = [
    ("1",  172_089_746, 7e-12),   # rs186897198
    ("7",  18_747_194,  7e-23),   # rs727851
    ("10", 30_838_237,  7e-08),   # rs12413361
    ("11", 12_857_576,  5e-09),   # rs546512774
    ("11", 2_780_860,   9e-08),   # rs234886
    ("12", 122_933_684, 7e-09),   # rs73230017 - hg38-specific
    ("16", 69_181_056,  1e-08),   # rs12444184 - hg38-specific
    ("17", 67_858_486,  2e-12),   # rs55931203
]

# Half-width of the peak window in bp — variants *outside* every target
# window are guaranteed to have p ≥ 5e-8 (no accidental noise looks
# significant).  Variants *inside* a window may carry the peak-shape
# p-values described above.
WINDOW_HALF_SIZE: int = 250_000

# Decay constant for the exponential peak shape (in bp).  ``p(d) = min_p
# · 10^(|d|/PEAK_TAU_BP)`` — at 100 kb the p-value is 10× weaker than the
# lead, at 200 kb it's 100× weaker.
PEAK_TAU_BP: int = 100_000

# Number of supporting SNPs to plant around each lead.  Real GWAS peaks
# typically show 10–50 sub-threshold supporting SNPs from LD tagging;
# 20 gives a visually convincing peak without inflating file size.
N_SUPPORT_PER_PEAK: int = 20


def _targets_for_build(build: Optional[str], override: Optional[str] = None
                       ) -> list[tuple[str, int, float]]:
    """Return the spike list matching *build* (or the *override* argument).

    * ``override='hg19' / 'hg38'`` picks that build regardless of *build*.
    * ``override='off'`` returns an empty list (falls back to legacy
      random-signal behaviour).
    * ``override='auto'`` (or ``None``) uses *build*, defaulting to hg38
      when no build is specified.
    """
    src = (override or "auto").lower()
    if src == "off":
        return []
    if src == "hg19":
        return HG19_TARGET_SPIKES
    if src == "hg38":
        return HG38_TARGET_SPIKES
    # auto
    if build is None:
        return HG38_TARGET_SPIKES
    b = str(build).strip().lower()
    if b in ("hg19", "grch37", "b37", "19"):
        return HG19_TARGET_SPIKES
    return HG38_TARGET_SPIKES


def _clip_random_noise(df, target_spikes, rng):
    """Push every non-target-region variant to p >= 5e-8.

    Guarantees the only genome-wide-significant hits are the injected
    spikes.  Rows inside a ±WINDOW_HALF_SIZE window around any target
    on the matching chromosome are left alone (so peak-shape supporting
    SNPs keep their sub-threshold p-values).
    """
    if not target_spikes:
        return df
    near = pd.Series(False, index=df.index)
    for chrom, pos, _ in target_spikes:
        m = ((df["CHR"].astype(str) == str(chrom))
             & (df["BP"].sub(pos).abs() <= WINDOW_HALF_SIZE))
        near = near | m
    accidental = (df["P"] < 5e-8) & (~near)
    n_bad = int(accidental.sum())
    if n_bad:
        df.loc[accidental, "P"] = rng.uniform(5e-8, 1.0, size=n_bad)
    return df


def _inject_target_spikes(df, target_spikes, rng, build_label=None,
                          n_support: int = N_SUPPORT_PER_PEAK):
    """Add a lead SNP + ``n_support`` supporting SNPs at each target locus.

    Peak shape:
      * Lead SNP at exactly ``(chr, pos)`` with ``p = min_p_value``.
      * Supporting SNPs at Gaussian-sampled offsets around the lead
        (σ = WINDOW_HALF_SIZE / 3, clipped to ±WINDOW_HALF_SIZE) with
        ``p(d) = min_p · 10^(|d|/PEAK_TAU_BP)`` × a small multiplicative
        jitter in [0.5, 2.0].  Result: dense significant cloud right
        at the lead, sub-threshold tail out to the window edge.
    """
    if not target_spikes:
        return df
    new_rows = []
    for chrom, pos, min_p in target_spikes:
        # Lead SNP at exact position.
        new_rows.append(_make_spike_row(
            chrom, pos, min_p, kind="lead",
            build_label=build_label,
            snp_id=f"rs_lead_{chrom}_{pos}",
        ))
        # Supporting SNPs — sample offsets, then repeatedly top up any
        # that fell outside the ±WINDOW_HALF_SIZE window.
        sigma = WINDOW_HALF_SIZE / 3.0
        _offsets = rng.normal(0.0, sigma, size=n_support * 2)
        _offsets = _offsets[np.abs(_offsets) <= WINDOW_HALF_SIZE][:n_support]
        while len(_offsets) < n_support:
            _more = rng.normal(0.0, sigma, size=n_support)
            _more = _more[np.abs(_more) <= WINDOW_HALF_SIZE]
            _offsets = np.concatenate([_offsets, _more])
        _offsets = _offsets[:n_support]

        _support_p = min_p * np.power(10.0, np.abs(_offsets) / PEAK_TAU_BP)
        _support_p = _support_p * rng.uniform(0.5, 2.0, size=n_support)
        _support_p = np.clip(_support_p, 1e-300, 1.0)

        for i, (off, p) in enumerate(zip(_offsets, _support_p)):
            _bp = int(pos + off)
            if _bp < 1:
                continue
            new_rows.append(_make_spike_row(
                chrom, _bp, float(p), kind="support",
                build_label=build_label,
                snp_id=f"rs_supp_{chrom}_{pos}_{i}",
            ))

    if not new_rows:
        return df
    spike_df = pd.DataFrame(new_rows)
    # Match the target frame's column set (and fill any missing).
    for col in df.columns:
        if col not in spike_df.columns:
            spike_df[col] = df[col].iloc[0] if len(df) else None
    return pd.concat([df, spike_df[df.columns]], ignore_index=True)


def _make_spike_row(chrom, pos, p, kind, build_label, snp_id):
    """Construct a single spike-in row dict matching the emitted schema."""
    row = {
        "CHR":  str(chrom),
        "SNP":  snp_id,
        "BP":   int(pos),
        "A1":   "A", "A2": "G",
        "BETA": 0.15 if kind == "lead" else 0.05,
        "SE":   0.02,
        "P":    float(p),
    }
    if build_label:
        row["BUILD"] = build_label
    return row

DATASET_SIZES = {
    "500K":  500_000,
    "1M":  1_000_000,
    "2M":  2_000_000,
    "5M":  5_000_000,
    "10M": 10_000_000,
    # -- Added for the 100M-variant scaling demo --------------------
    "25M":  25_000_000,
    "50M":  50_000_000,
    "100M": 100_000_000,
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
                      build: str | None = None,
                      targets: str = "auto") -> pd.DataFrame:
    """Generate synthetic GWAS summary statistics (single-shot, in-memory).

    Parameters
    ----------
    n_variants : int
        Total number of variants.
    n_signals : int
        **Legacy parameter, ignored when ``targets != 'off'``.**  When
        ``targets='off'`` the pre-0.4.x random-signal behaviour is used
        (``n_signals`` variants scattered at random with p ∈ [1e-50,
        5e-8]).  Kept in the signature so callers pinning the old
        keyword still import.
    seed : int
        Random seed for reproducibility.
    build : {"hg19", "hg38", None}
        When set, an additional ``BUILD`` column is emitted and the
        chromosome-size table is picked accordingly.  ``None``
        preserves the pre-liftover behaviour of the script (no
        BUILD column, hg38 coordinates).
    targets : {"auto", "hg19", "hg38", "off"}, optional
        Whether to plant peak-shape spikes at the well-known target
        loci in :data:`HG19_TARGET_SPIKES` / :data:`HG38_TARGET_SPIKES`.

        * ``"auto"`` (default) — pick the list matching *build*,
          defaulting to hg38 when no build is set.
        * ``"hg19"`` / ``"hg38"`` — force a specific list.
        * ``"off"`` — no target spikes; fall back to the legacy
          random-signal behaviour driven by *n_signals*.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns CHR, SNP, BP, A1, A2, BETA, SE, P
        (and BUILD when ``build`` is passed).

    Notes
    -----
    In the target-spike mode every ``P < 5e-8`` row lies inside a
    ±:data:`WINDOW_HALF_SIZE` window of one of the target loci, and
    each locus carries one lead SNP at the exact target position with
    ``p = min_p_value``.
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

    _spikes = _targets_for_build(build, override=targets)
    if _spikes:
        # Target-spike mode: guarantee no accidental noise looks
        # significant, then plant peak-shape signals at each target.
        df = _clip_random_noise(df, _spikes, rng)
        df = _inject_target_spikes(df, _spikes, rng=rng, build_label=build)
    else:
        # Legacy random-signal mode (``targets='off'``).
        signal_indices = rng.choice(len(df), size=n_signals, replace=False)
        for idx in signal_indices:
            df.loc[idx, "P"] = 10 ** rng.uniform(-50, -8)

    # Ensure no p=0 or p>1
    df["P"] = df["P"].clip(1e-300, 1.0)
    return df


def stream_sumstats_to_tsv(out_path: str, n_variants: int,
                           n_signals: int = 30, seed: int = 42,
                           build: str | None = None,
                           targets: str = "auto",
                           chunk_bytes_hint: int = 1 << 30) -> None:
    """Write a synthetic sumstats TSV chromosome-by-chromosome.

    Streaming counterpart to :func:`generate_sumstats`, used
    automatically by ``main()`` for the 25M / 50M / 100M sizes where
    the full-genome frame would exhaust node memory.  Peak Python
    memory is bounded by the largest single-chromosome frame
    (chr1 ≈ 10 % of total) rather than the whole file.

    Target-spike injection (``targets != 'off'``) is applied
    per-chromosome: each chunk's frame carries any target loci whose
    chromosome matches, and the peak-shape supporting SNPs are
    generated inside that same chunk so they fit alongside the random
    background for that chromosome.  The ``targets='off'`` path
    preserves the pre-target legacy behaviour (random ``n_signals``
    hits sampled from the global index space).
    """
    rng = np.random.default_rng(seed)
    chrom_sizes = _chrom_sizes_for_build(build) if build else HG38_CHROM_SIZES
    counts = _distribute_variants(n_variants, chrom_sizes)
    total = sum(counts.values())

    _spikes = _targets_for_build(build, override=targets)

    # Only the legacy ``targets='off'`` path uses the global random-
    # signal sampling; target-spike mode ignores ``n_signals`` and
    # plants a fixed peak per matching target chromosome.
    if not _spikes:
        signal_global_idx = np.sort(
            rng.choice(total, size=n_signals, replace=False)
        )
        signal_p = 10 ** rng.uniform(-50, -8, size=n_signals)
    else:
        signal_global_idx = np.array([], dtype=int)
        signal_p = np.array([], dtype=float)

    snp_offset = 0
    header_written = False
    with open(out_path, "w") as fh:
        for chrom, count in counts.items():
            df = _chrom_frame(
                chrom, count, snp_offset,
                max_bp=chrom_sizes[chrom],
                rng=rng, build_label=build,
            )
            if _spikes:
                # Target-spike mode: clip accidental noise on this
                # chromosome, then plant the chunk's peaks.
                _chr_targets = [t for t in _spikes if str(t[0]) == str(chrom)]
                df = _clip_random_noise(df, _chr_targets, rng)
                if _chr_targets:
                    df = _inject_target_spikes(
                        df, _chr_targets, rng=rng, build_label=build,
                    )
            else:
                # Legacy random-signal path.
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
    parser.add_argument(
        "--targets", choices=["auto", "hg19", "hg38", "off"], default="auto",
        help=(
            "Which target-locus list to spike into the output.  "
            "'auto' (default) matches --build (falling back to hg38 "
            "when no build is set); 'hg19'/'hg38' force a specific "
            "list; 'off' falls back to the legacy random-signal "
            "generator that scatters n_signals hits across the "
            "genome with no controlled positions.  See "
            ":data:`HG19_TARGET_SPIKES` / :data:`HG38_TARGET_SPIKES` "
            "for the coordinate lists (10 body-height associated gwas-catalog " 
            "loci: rs186897198, rs727851, rs2030645, rs6606792, rs73230017, rs12444184"
            "rs12413361, rs546512774, rs234886, rs55931203)."
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
                                   seed=args.seed, build=build,
                                   targets=args.targets)
        else:
            df = generate_sumstats(n_variants, seed=args.seed, build=build,
                                   targets=args.targets)
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
