"""Batch CLI for editing the pycmplot hits overlay TSV.

This is the scriptable companion to :mod:`pycmplot.editor` (browser
GUI) and :mod:`pycmplot.editor_tui` (terminal UI).  It's built for
Makefiles, pipelines, and reproducible analysis notebooks — a step
like "colour every genome-wide-significant novel hit red" belongs in
a script, not a GUI.

Two subcommands::

    pycmplot hits list                        \\
        --cache_dir ./.pycmplot_cache         \\
        [--group KEY]                         \\
        [--where 'P < 5e-8']                  \\
        [--columns CHR,POS,SNP,category]

    pycmplot hits set                         \\
        --cache_dir ./.pycmplot_cache         \\
        [--group KEY]                         \\
        --where 'P < 5e-8'                    \\
        [--color red]                         \\
        [--category "genome-wide"]            \\
        [--dry-run]

The ``--where`` expression uses pandas ``DataFrame.query()`` grammar,
so it accepts any of the usual operators — ``<``, ``<=``, ``==``,
``!=``, ``in``, ``and``/``or``, ``str.contains(...)`` on string
columns, backtick-quoted column names for names that contain
special characters, and so on.  This is the same grammar the TUI's
inline filter uses, so users learn one syntax.

Writes go through :func:`pycmplot.cache.write_hits_overlay`, so the
atomic-write and inheritance-across-regeneration guarantees hold
identically to the GUI backends and to hand-editing the TSV.

Dispatch happens in :mod:`pycmplot._core` via the same preflight
that routes ``pycmplot edit`` — a leading positional ``hits`` is
peeled off before argparse sees it, so existing
``pycmplot --sum_stats …`` invocations are unaffected.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Public entry point (called from pycmplot._core.main preflight)
# ---------------------------------------------------------------------------

def dispatch(argv: list[str]) -> None:
    """Parse ``argv`` (the post-``hits`` slice of ``sys.argv``) and run.

    ``argv[0]`` is expected to be the subcommand name (``list`` or
    ``set``); anything else is a friendly error.  We keep the parser
    hand-rolled rather than using argparse subparsers because the
    outer preflight already peels off the top-level positional, and
    threading two argparse layers through here gains nothing.
    """
    if not argv or argv[0] in ("-h", "--help"):
        _print_top_help()
        return
    sub = argv[0]
    if sub == "list":
        _cmd_list(argv[1:])
    elif sub == "set":
        _cmd_set(argv[1:])
    else:
        raise SystemExit(
            f"Unknown `pycmplot hits` subcommand: {sub!r}. "
            "Expected one of: list, set."
        )


def _print_top_help() -> None:
    print(
        "Usage: pycmplot hits <subcommand> [options]\n"
        "\n"
        "Subcommands:\n"
        "  list   Print rows from the hits overlay (with optional --where filter)\n"
        "  set    Apply --color / --category to every row matching --where\n"
        "\n"
        "Both subcommands accept:\n"
        "  --cache_dir PATH  (default: .pycmplot_cache)\n"
        "  --group KEY       (auto-detected when only one exists)\n"
        "\n"
        "Run `pycmplot hits list -h` or `pycmplot hits set -h` for details."
    )


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _resolve(cache_dir: str, group: Optional[str]) -> tuple[str, str]:
    """Return ``(absolute_cache_dir, group_key)`` or raise SystemExit."""
    if not Path(cache_dir).exists():
        raise SystemExit(f"cache_dir does not exist: {cache_dir}")
    # Reuse the browser backend's group-picker so semantics match.
    from pycmplot.editor import _resolve_group
    return str(Path(cache_dir).resolve()), _resolve_group(cache_dir, group)


def _load_overlay(cache_dir: str, group_key: str):
    """Read the overlay + sidecar; SystemExit if the TSV is unreadable."""
    from pycmplot.cache import read_hits_overlay
    df, meta = read_hits_overlay(cache_dir, group_key)
    if df is None:
        raise SystemExit(
            f"Could not read hits.{group_key}.tsv. "
            "Re-run pycmplot with --cache to regenerate."
        )
    return df, meta


def _apply_where(df, where: Optional[str]):
    """Apply a df.query() expression; return the filtered view.

    ``df.query()`` errors are turned into a crisp ``SystemExit`` so
    the user sees the offending column name / operator rather than a
    pandas traceback.
    """
    if not where:
        return df
    try:
        return df.query(where, engine="python")
    except Exception as exc:
        raise SystemExit(
            f"--where expression is not a valid pandas query: {where!r}\n"
            f"  ({exc.__class__.__name__}: {exc})\n"
            "See https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html"
        )


# ---------------------------------------------------------------------------
# `pycmplot hits list`
# ---------------------------------------------------------------------------

def _cmd_list(argv: list[str]) -> None:
    p = argparse.ArgumentParser(
        prog="pycmplot hits list",
        description=(
            "Print rows from the hits overlay.  Supports the same "
            "--where expression syntax as `pycmplot hits set` (pandas "
            "DataFrame.query grammar)."
        ),
    )
    p.add_argument("--cache_dir", default=".pycmplot_cache")
    p.add_argument("--group", default=None)
    p.add_argument("--where", default=None,
                   help="Pandas query expression, e.g. 'P < 5e-8'.")
    p.add_argument(
        "--columns", default=None,
        help=(
            "Comma-separated subset of columns to print.  Default: "
            "every column.  Useful for `awk`-friendly output."
        ),
    )
    p.add_argument(
        "--tsv", action="store_true",
        help=(
            "Emit the result as TSV (headerless-friendly for piping "
            "into `awk` or `cut`).  Default is a padded human-readable "
            "table."
        ),
    )
    args = p.parse_args(argv)

    cache_dir, group_key = _resolve(args.cache_dir, args.group)
    df, _ = _load_overlay(cache_dir, group_key)
    filtered = _apply_where(df, args.where)

    if args.columns:
        wanted = [c.strip() for c in args.columns.split(",") if c.strip()]
        missing = [c for c in wanted if c not in filtered.columns]
        if missing:
            raise SystemExit(
                f"--columns includes unknown name(s): {missing}. "
                f"Available: {list(filtered.columns)}"
            )
        filtered = filtered[wanted]

    if args.tsv:
        # to_csv with sep='\t' keeps everything script-parseable.
        sys.stdout.write(filtered.to_csv(sep="\t", index=False))
    else:
        # Rely on pandas' default renderer — it handles wide/narrow
        # sensibly and truncates rather than wrapping.
        import pandas as pd
        with pd.option_context("display.max_rows", None,
                               "display.max_columns", None,
                               "display.width", 200):
            print(filtered.to_string(index=False))

    # Summary line to stderr so it doesn't interfere with TSV piping.
    print(f"({len(filtered)} row(s) of {len(df)})", file=sys.stderr)


# ---------------------------------------------------------------------------
# `pycmplot hits set`
# ---------------------------------------------------------------------------

def _cmd_set(argv: list[str]) -> None:
    p = argparse.ArgumentParser(
        prog="pycmplot hits set",
        description=(
            "Apply --color and/or --category to every row matching "
            "--where.  Writes atomically through write_hits_overlay so "
            "the next loader call cache-HITs and inherits the edits."
        ),
        epilog=(
            "Example:\n"
            "  pycmplot hits set --cache_dir ./.pycmplot_cache \\\n"
            "      --where 'P < 5e-8 and category == \"significant\"' \\\n"
            "      --color '#00cc44' --category 'genome-wide'"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--cache_dir", default=".pycmplot_cache")
    p.add_argument("--group", default=None)
    p.add_argument(
        "--where", required=True,
        help=(
            "Pandas query expression selecting the rows to mutate.  "
            "Use `--where 'True'` to match every row."
        ),
    )
    p.add_argument(
        "--color", default=None, metavar="VALUE",
        help=(
            "Value to write into the highlight_color column.  Any "
            "matplotlib-parseable value (name, `#rrggbb`, RGB tuple, "
            "or the sentinel `auto`).  Invalid values are rejected."
        ),
    )
    p.add_argument(
        "--category", default=None, metavar="VALUE",
        help="Value to write into the category column.",
    )
    p.add_argument(
        "--dry-run", action="store_true",
        help=(
            "Print the rows that would change but don't write.  "
            "Recommended before running a broad --where."
        ),
    )
    args = p.parse_args(argv)

    if args.color is None and args.category is None:
        raise SystemExit(
            "Nothing to do: pass at least one of --color / --category."
        )
    if args.color is not None:
        # Fail fast on a typo — same validation logic the plotter uses.
        from matplotlib.colors import is_color_like
        v = args.color.strip()
        if v.lower() not in {"auto", ""} and not is_color_like(v):
            raise SystemExit(
                f"--color {args.color!r} is not a matplotlib-parseable "
                "colour and not the sentinel `auto`.  Refusing to write."
            )

    cache_dir, group_key = _resolve(args.cache_dir, args.group)
    df, meta = _load_overlay(cache_dir, group_key)
    if not meta or "auto_key" not in meta:
        raise SystemExit(
            "Overlay sidecar is missing `auto_key`; refusing to save "
            "(would break the next cache-HIT check).  Re-run pycmplot "
            "with --cache to regenerate the sidecar."
        )

    filtered = _apply_where(df, args.where)
    if filtered.empty:
        print(f"--where matched 0 row(s); nothing to change.")
        return

    # Preview so the user can see what changed, in both dry-run and
    # real-write flows.  Only prints the columns we're actually
    # touching + a few identifiers, to keep the output scannable.
    _preview_cols = [
        c for c in ("CHR", "POS", "SNP", "P", "LABEL",
                    "highlight_color", "category")
        if c in filtered.columns
    ]
    print(f"Matching {len(filtered)} row(s) of {len(df)}:")
    import pandas as pd
    with pd.option_context("display.max_rows", 40,
                           "display.max_columns", None,
                           "display.width", 200):
        print(filtered[_preview_cols].to_string(index=False))

    if args.dry_run:
        _msg = []
        if args.color is not None:
            _msg.append(f"highlight_color = {args.color!r}")
        if args.category is not None:
            _msg.append(f"category = {args.category!r}")
        print(f"\n[dry-run] would set {' and '.join(_msg)} on those rows.")
        return

    # Mutate in place on the full frame (not the filtered view) so
    # index alignment is unambiguous.
    idx = filtered.index
    if args.color is not None:
        df.loc[idx, "highlight_color"] = args.color
    if args.category is not None:
        df.loc[idx, "category"] = args.category

    # Save via write_hits_overlay — auto rows first, user rows below,
    # atomic replace, sidecar meta preserved.
    from pycmplot.cache import AUTO_TAG, SOURCE_COL, write_hits_overlay
    auto_df = df[df[SOURCE_COL].astype(str) == AUTO_TAG].copy()
    written = write_hits_overlay(
        cache_dir, auto_df,
        auto_key=meta["auto_key"],
        group_key=group_key,
        preserve_user_from=df,
    )
    _changed = []
    if args.color is not None:
        _changed.append(f"highlight_color = {args.color!r}")
    if args.category is not None:
        _changed.append(f"category = {args.category!r}")
    print(
        f"\nSet {' and '.join(_changed)} on {len(filtered)} row(s).\n"
        f"Wrote: {written}"
    )
