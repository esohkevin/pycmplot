"""Streamlit-based GUI for editing the pycmplot hits overlay.

The overlay TSV that :func:`~pycmplot.io.get_sumstats_and_merged_sector_list`
writes to ``<cache_dir>/annotations/hits.<group>.tsv`` is designed to
be hand-edited (see :ref:`tut-overlay`).  This module wraps the same
file in a lightweight browser UI so users can pick highlight colours
visually, autocomplete categories, and add/remove ``source=user`` rows
without touching a text editor.

The UI is delivered as an **optional extra** — Streamlit is not a
runtime dependency of pycmplot.  Install it with::

    pip install "pycmplot[editor]"

then launch::

    pycmplot edit --cache_dir ./.pycmplot_cache

The subcommand dispatch lives in :mod:`pycmplot._core`; this module
provides two entry points:

* :func:`launch` — CLI-side helper that boots a Streamlit server on
  this file.  Called from ``pycmplot edit``.
* :func:`_run_app` — the Streamlit script itself, invoked by
  Streamlit's runtime when the server starts.  Not intended for
  direct import.

Design notes
------------
* Editing is done via :func:`streamlit.data_editor` with per-column
  ``ColumnConfig`` (color as text with a live swatch, category as a
  selectbox with type-to-add, source as an ``auto``/``user`` dropdown).
* Save goes through :func:`pycmplot.cache.write_hits_overlay`, so the
  same atomic-write + inheritance-across-regeneration guarantees the
  Python API relies on hold here too.  No parallel write path.
* Live plot preview re-invokes :func:`pycmplot.plotting.linear.plot_linear`
  (or :func:`~pycmplot.plotting.circular.plot_circular`) on the
  in-memory edited overlay and renders the figure via
  :func:`streamlit.pyplot`.  This means the plotters must be able to
  consume the edited overlay directly — which they already do, so
  there's no plotter-side change.
* Multi-group caches (multiple ``hits.<group>.tsv`` files under one
  ``cache_dir``) get a group-picker sidebar; single-group caches skip
  it entirely.
"""

from __future__ import annotations

import argparse
import glob
import logging
import os
import re
import sys
from pathlib import Path
from typing import List, Optional, Tuple

logger = logging.getLogger(__name__)

# Environment variables used to pass state from ``launch()`` (the CLI
# entry) into ``_run_app()`` (the Streamlit-side entry).  Streamlit
# spawns the script in a subprocess without argv, so a small
# environment handshake is the cleanest way to pass configuration.
_ENV_CACHE_DIR = "PYCMPLOT_EDIT_CACHE_DIR"
_ENV_GROUP = "PYCMPLOT_EDIT_GROUP"


# ---------------------------------------------------------------------------
# Group discovery helpers (usable without streamlit installed)
# ---------------------------------------------------------------------------

def _list_groups(cache_dir: str | os.PathLike) -> List[str]:
    """Return the group_keys for every hits overlay under *cache_dir*.

    Filenames look like ``hits.<8-hex-chars>.tsv``.  Anything not
    matching that pattern is ignored so stray files in the annotations
    folder don't corrupt the picker.
    """
    ann_dir = Path(cache_dir) / "annotations"
    if not ann_dir.exists():
        return []
    pat = re.compile(r"^hits\.([0-9a-f]{4,64})\.tsv$")
    groups: List[str] = []
    for p in sorted(ann_dir.glob("hits.*.tsv")):
        m = pat.match(p.name)
        if m:
            groups.append(m.group(1))
    return groups


def _resolve_group(cache_dir: str, explicit_group: Optional[str]) -> str:
    """Pick a group_key, either from the CLI arg or the sole discovered one.

    Raises ``SystemExit`` with a user-facing message when nothing can
    be resolved — the editor is useless without a target file, and
    Streamlit swallows plain exceptions inside ``_run_app``.
    """
    if explicit_group:
        return explicit_group
    groups = _list_groups(cache_dir)
    if not groups:
        raise SystemExit(
            f"No hits overlays found under {cache_dir}/annotations. "
            "Run pycmplot with --cache first to generate one."
        )
    if len(groups) == 1:
        return groups[0]
    # Multiple groups: don't pick blindly; ask the user.
    listing = "\n  ".join(groups)
    raise SystemExit(
        f"Multiple hits overlays found under {cache_dir}/annotations.\n"
        f"Re-run with --group <key>. Available:\n  {listing}"
    )


# ---------------------------------------------------------------------------
# CLI-side entry — boots streamlit on THIS file
# ---------------------------------------------------------------------------

def launch(cache_dir: str,
           group: Optional[str] = None,
           host: str = "localhost",
           port: int = 8501) -> None:
    """Boot a Streamlit server rendering this module's editor.

    Called from ``pycmplot edit ...``.  Fails fast with an actionable
    message when Streamlit isn't installed (users see the extras hint
    rather than a raw ``ModuleNotFoundError``).
    """
    # Sanity-check the cache_dir + group before boot so users get a
    # crisp error at the terminal rather than a Streamlit "Something
    # went wrong" toast.
    if not Path(cache_dir).exists():
        raise SystemExit(f"cache_dir does not exist: {cache_dir}")
    resolved_group = _resolve_group(cache_dir, group)

    # Environment handshake — see the _ENV_* constants above.
    os.environ[_ENV_CACHE_DIR] = str(Path(cache_dir).resolve())
    os.environ[_ENV_GROUP] = resolved_group

    # IMPORTANT: set Streamlit config env vars BEFORE importing
    # ``streamlit.web.bootstrap``.  Streamlit's config system reads
    # its env vars at first-import; setting them afterwards is a
    # silent no-op on server.port / server.address (the reason we
    # were previously seeing every editor session bind to 8501
    # regardless of --port).
    os.environ["STREAMLIT_SERVER_PORT"] = str(port)
    os.environ["STREAMLIT_SERVER_ADDRESS"] = host
    os.environ["STREAMLIT_SERVER_HEADLESS"] = "true"
    os.environ["STREAMLIT_BROWSER_GATHER_USAGE_STATS"] = "false"

    try:
        # ``streamlit.web.bootstrap`` is the internal API Streamlit
        # itself uses when you run ``streamlit run script.py``.  We
        # call it directly so ``pycmplot edit`` is a single command
        # rather than a shell recipe wrapping ``streamlit run``.
        import streamlit.web.bootstrap as _boot  # type: ignore
    except ImportError as exc:
        raise SystemExit(
            "The 'edit' subcommand requires Streamlit. Install it with:\n"
            "    pip install \"pycmplot[editor]\"\n"
            f"(underlying error: {exc})"
        )

    # ``flag_options`` is a belt-and-braces fallback; the env-vars set
    # above are the primary channel (see comment).
    flag_options = {
        "server.address": host,
        "server.port": port,
        "server.headless": True,
        "browser.gatherUsageStats": False,
    }
    script_path = str(Path(__file__).resolve())

    print(
        f"pycmplot editor: launching Streamlit at http://{host}:{port}\n"
        f"  cache_dir : {os.environ[_ENV_CACHE_DIR]}\n"
        f"  group     : {resolved_group}\n"
        "  Ctrl-C to stop."
    )
    # ``bootstrap.run`` sets up config *watchers* but doesn't apply
    # ``flag_options`` immediately — the ``Server`` it constructs
    # otherwise reads server.port / server.address from the already-
    # loaded defaults, silently ignoring our overrides.  Applying
    # ``load_config_options`` explicitly first fixes this.
    try:
        _boot.load_config_options(flag_options)
    except Exception:  # pragma: no cover — very old Streamlit
        pass

    # Signature has shifted slightly across Streamlit versions; the
    # positional args are stable but newer builds accept a fifth
    # ``args`` list.  Try the modern signature first, fall back to the
    # older one.
    try:
        _boot.run(script_path, is_hello=False, args=[], flag_options=flag_options)
    except TypeError:  # pragma: no cover — older Streamlit
        _boot.run(script_path, "", [], flag_options)


# ---------------------------------------------------------------------------
# Streamlit-side entry — the actual UI
# ---------------------------------------------------------------------------

def _run_app() -> None:
    """The Streamlit script.  Invoked by ``streamlit run`` on this file."""
    # Local imports keep pycmplot importable when streamlit is absent.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    import streamlit as st

    from pycmplot.annotation import (
        CATEGORY_COL, CATEGORY_DEFAULT,
        HIGHLIGHT_COLOR_COL, HIGHLIGHT_COLOR_AUTO,
    )
    from pycmplot.cache import (
        AUTO_TAG, SOURCE_COL, USER_TAG,
        read_hits_overlay, write_hits_overlay,
    )

    st.set_page_config(page_title="pycmplot hits editor", layout="wide")

    cache_dir = os.environ.get(_ENV_CACHE_DIR)
    group_key = os.environ.get(_ENV_GROUP)
    if not cache_dir or not group_key:
        st.error(
            "Editor was launched without the required environment. "
            "Use `pycmplot edit --cache_dir <path>` rather than "
            "`streamlit run editor.py`."
        )
        return

    st.title("pycmplot hits overlay editor")
    st.caption(
        f"Editing `hits.{group_key}.tsv` in `{cache_dir}/annotations`. "
        "Changes are held in memory until you press **Save**."
    )

    df, meta = read_hits_overlay(cache_dir, group_key)
    if df is None:
        st.error(
            f"Could not read hits.{group_key}.tsv. If it's corrupted, "
            "run pycmplot with --cache once to regenerate it from the "
            "cached leads (user edits will need to be re-applied)."
        )
        return

    # Existing categories, used as autocomplete options for the
    # SelectboxColumn.  Include the sentinel default so users can
    # revert an edited row cleanly.
    known_cats = sorted({str(v) for v in df.get(CATEGORY_COL, []).astype(str)
                         if v and str(v).lower() not in {"nan", "none"}})
    if CATEGORY_DEFAULT not in known_cats:
        known_cats = [CATEGORY_DEFAULT, *known_cats]

    # Column configs mirror the schema documented in the tutorial.
    # Anything not listed here is rendered with sensible defaults.
    col_cfg = {
        SOURCE_COL: st.column_config.SelectboxColumn(
            "source", options=[AUTO_TAG, USER_TAG],
            help="`auto` rows come from the loader; `user` rows are yours to add.",
            required=True,
        ),
        HIGHLIGHT_COLOR_COL: st.column_config.TextColumn(
            "highlight_color",
            help=(
                "Matplotlib colour name, `#rrggbb`, or `auto` to fall "
                "back to the plotter's default. Invalid values fall "
                "back too (with a warning)."
            ),
            default=HIGHLIGHT_COLOR_AUTO,
        ),
        CATEGORY_COL: st.column_config.SelectboxColumn(
            "category", options=known_cats,
            help=(
                "Legend label. Rows sharing a category share a legend "
                "entry. Type a new value to add it."
            ),
            required=True,
            default=CATEGORY_DEFAULT,
        ),
    }

    st.subheader("Rows")
    edited = st.data_editor(
        df,
        column_config=col_cfg,
        num_rows="dynamic",   # allow adding user rows inline
        width="stretch",
        key="hits_editor",
    )

    # ---- Colour-swatch preview (below the editor, since ``data_editor``
    # doesn't support inline cell rendering as of Streamlit 1.32).
    st.subheader("Colour preview")
    _render_swatches(st, edited)

    # ---- Actions ---------------------------------------------------------
    st.subheader("Actions")
    col_save, col_reload, col_preview = st.columns(3)

    if col_save.button("Save", type="primary", width="stretch"):
        _do_save(st, cache_dir, group_key, edited, df, meta)

    if col_reload.button("Discard & reload", width="stretch"):
        # Rerun re-reads from disk on the next tick; clearing the
        # data_editor's session key drops the in-memory buffer.
        st.session_state.pop("hits_editor", None)
        st.rerun()

    if col_preview.button("Preview plot", width="stretch"):
        _do_preview(st, cache_dir, group_key, edited)

    with st.expander("Overlay metadata"):
        st.json(meta or {"note": "no sidecar metadata found"})


# ---------------------------------------------------------------------------
# Small render/action helpers, split out to keep _run_app readable
# ---------------------------------------------------------------------------

def _render_swatches(st, df) -> None:
    """Render one small chip per row so users can spot colour typos.

    Uses a plain HTML table (via ``st.markdown(unsafe_allow_html=True)``)
    to keep the swatch grid dense.  ``matplotlib.colors.is_color_like``
    tells us whether a value would round-trip through the plotter.
    """
    from matplotlib.colors import is_color_like, to_hex
    from pycmplot.annotation import (
        CATEGORY_COL, HIGHLIGHT_COLOR_AUTO, HIGHLIGHT_COLOR_COL,
    )

    if df is None or df.empty:
        st.caption("Nothing to preview.")
        return
    if HIGHLIGHT_COLOR_COL not in df.columns:
        st.caption("No `highlight_color` column in this overlay.")
        return

    parts = ["<div style='display:flex;flex-wrap:wrap;gap:6px;'>"]
    for _, row in df.iterrows():
        raw = str(row.get(HIGHLIGHT_COLOR_COL, HIGHLIGHT_COLOR_AUTO)).strip()
        label = str(row.get(CATEGORY_COL, "")) or "—"
        # Sentinel/blank/invalid → grey with a dashed border, so users
        # can see at a glance which rows will fall back to the default.
        if not raw or raw.lower() in {"nan", "none", HIGHLIGHT_COLOR_AUTO}:
            css = "background:#eee;border:1px dashed #999;color:#666"
            tip = "auto → falls back to plotter default"
        elif not is_color_like(raw):
            css = "background:#fdd;border:1px solid #c33;color:#900"
            tip = f"invalid: {raw}"
        else:
            css = f"background:{to_hex(raw)};border:1px solid #333;color:#fff"
            tip = raw
        parts.append(
            f"<span title='{tip}' "
            f"style='display:inline-block;padding:4px 8px;border-radius:4px;"
            f"font-size:12px;{css};'>{label}</span>"
        )
    parts.append("</div>")
    st.markdown("".join(parts), unsafe_allow_html=True)


def _do_save(st, cache_dir, group_key, edited, original, meta) -> None:
    """Delegate to write_hits_overlay so we inherit its guarantees."""
    from pycmplot.cache import write_hits_overlay
    try:
        # Split auto vs user so write_hits_overlay's auto-block +
        # preserve-user-below layout is preserved.
        from pycmplot.cache import SOURCE_COL, AUTO_TAG
        auto_df = edited[edited[SOURCE_COL].astype(str) == AUTO_TAG].copy()
        # ``auto_key`` from the previous sidecar is the right one to
        # carry forward — the *auto payload identity* hasn't changed
        # (only user-facing columns have), so the next loader call
        # will still cache-HIT and inherit these edits.
        auto_key = (meta or {}).get("auto_key")
        if not auto_key:
            st.error(
                "No `auto_key` in the sidecar metadata; refusing to "
                "save (would break the next cache-HIT check). Re-run "
                "pycmplot with --cache to regenerate the sidecar."
            )
            return
        written = write_hits_overlay(
            cache_dir, auto_df,
            auto_key=auto_key,
            group_key=group_key,
            preserve_user_from=edited,
        )
        st.success(f"Saved: {written}")
    except Exception as exc:
        st.error(f"Save failed: {exc}")


def _do_preview(st, cache_dir, group_key, edited) -> None:
    """Render a linear Manhattan preview using the in-memory overlay.

    Preview is intentionally best-effort: it reuses whatever per-track
    cache entries already exist under *cache_dir* and skips raw
    sumstats reloading.  If the cache is incomplete for any reason,
    the preview simply refuses rather than triggering an expensive
    full re-load — the user can re-run their normal pycmplot command
    to regenerate first.
    """
    import matplotlib.pyplot as plt
    from pycmplot.cache import TrackCache
    from pycmplot import __version__ as _pcm_version

    tc = TrackCache(cache_dir, _pcm_version)
    # ``_metadata['tracks']`` maps cache_key -> entry.  We walk it to
    # collect one (label, cache_key) pair per label.  When a label
    # appears under more than one cache_key (rare — happens if the
    # user has cached the same label with different Stage-1 params in
    # the same cache_dir), we pick the most recently written entry so
    # the preview reflects the current parameter set as closely as
    # possible.
    tracks = getattr(tc, "_metadata", {}).get("tracks", {}) or {}
    by_label: dict = {}
    for cache_key, entry in tracks.items():
        lbl = entry.get("label")
        if not lbl:
            continue
        prev = by_label.get(lbl)
        if prev is None or entry.get("written", 0) > prev[1].get("written", 0):
            by_label[lbl] = (cache_key, entry)
    if not by_label:
        st.warning(
            "Preview needs at least one cached track under this "
            "cache_dir. Run pycmplot with --cache once first."
        )
        return

    # ``plot_linear`` expects sumstats_loaded to be a dict of
    # ``label -> [df, n_chroms]`` (list, indexed positionally in the
    # plotter).  A previous version of this preview passed raw
    # DataFrames and hit ``KeyError(0)``.
    dfs = {}
    for lbl, (ckey, _entry) in by_label.items():
        cached = tc.get(lbl, ckey)
        if cached is not None and getattr(cached, "df", None) is not None:
            n_chroms = int(getattr(cached, "extras", {}).get(
                "n_chroms",
                cached.df["CHR"].nunique() if "CHR" in cached.df.columns else 22,
            ))
            dfs[lbl] = [cached.df, n_chroms]
    if not dfs:
        st.warning(
            "Cached track metadata exists but the parquet files are "
            "missing or unreadable. Re-run pycmplot with --cache."
        )
        return

    fig, _ = plt.subplots(figsize=(12, 3 * max(1, len(dfs))))
    from pycmplot.plotting.linear import linear
    linear(
        sumstats_loaded=dfs, logp=True, highlight=True,
        hits_table=edited, annotate=False, ax=fig,
    )
    st.pyplot(fig, clear_figure=True)


# ---------------------------------------------------------------------------
# Dispatch: when Streamlit runs this file as a script, call _run_app().
# When it's imported as a normal module (tests, `pycmplot edit` launcher),
# do nothing.
# ---------------------------------------------------------------------------

def _is_streamlit_runtime() -> bool:
    """Detect whether this process is being executed by Streamlit.

    Streamlit sets ``STREAMLIT_SERVER_PORT`` and exposes an internal
    runtime module; either signal alone is enough to be confident.
    """
    if os.environ.get("STREAMLIT_SERVER_PORT"):
        return True
    try:
        import streamlit.runtime as _rt  # type: ignore
        return _rt.exists()
    except Exception:
        return False


if _is_streamlit_runtime():
    _run_app()
