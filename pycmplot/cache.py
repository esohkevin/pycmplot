"""
pycmplot.cache
==============

Per-track cache for the Stage-1 output of
:func:`~pycmplot.io.get_sumstats_and_merged_sector_list`.

Motivation
----------
Loading a multi-million-variant summary statistics file involves several
expensive steps: raw CSV read, dtype coercion, chromosome-name
normalisation, hg18/hg19 -> hg38 liftover, density-aware auto-thinning,
lead-SNP extraction, per-track statistics.  Together these dominate the
runtime of a multi-track plot call.  Stage-2 (matplotlib rendering) is
comparatively cheap.

When Stage 2 fails (SSH drop, notebook restart, kernel crash), the user
currently has to re-run Stage 1 for every input file from scratch.  This
module provides an on-disk, per-track cache so subsequent runs read
completed tracks straight from parquet and re-execute Stage 1 only for
tracks that are still missing or whose parameters have changed.

Layout
------
::

    <cache_dir>/
        metadata.json                          # cache_key + status + params per track
        tracks/
            <label>.parquet                    # main sumstats DataFrame
            <label>.leads.parquet              # lead SNP DataFrame (optional)
            <label>.pvals.npy                  # full raw p-values (optional)

``metadata.json`` is written atomically (temp file + rename) after every
per-track write, so a crash mid-run never leaves the cache in a
half-updated state.

Cache validation
----------------
Every cache entry records a ``cache_key`` computed from:

* SHA-256 of the raw input file (so any real content change invalidates)
* pycmplot version
* All Stage-1 parameters that would change the cached DataFrame
  (``trim_pval``, ``auto_thin``, ``auto_thin_threshold``,
  ``auto_thin_max_below``, ``logp``, ``highlight``, ``highlight_thresh``,
  ``build``, ...)

A cache hit requires an exact key match.  A mismatch triggers silent
regeneration (see :func:`compute_cache_key`).

Public API
----------
The module is used internally by
:func:`~pycmplot.io.get_sumstats_and_merged_sector_list`; direct use is
supported for advanced pipelines.

>>> from pycmplot.cache import TrackCache, compute_cache_key
>>> cache = TrackCache("/tmp/mycache")
>>> key = compute_cache_key(raw_path="/data/track1.tsv.gz",
...                         version="0.3.2",
...                         trim_pval=0.01, auto_thin=True)
>>> hit = cache.get("track1", key)         # None on miss / stale
>>> if hit is None:
...     df, leads, pvals, extras = ...     # run Stage 1
...     cache.put("track1", key, df, leads=leads, pvals=pvals, extras=extras)
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import shutil
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Hashing
# ---------------------------------------------------------------------------

_HASH_CHUNK = 1 << 20  # 1 MiB


def sha256_file(path: str | os.PathLike) -> str:
    """Return the hex SHA-256 digest of a file's byte contents.

    Reads the file in 1 MiB chunks so memory usage is bounded regardless
    of file size.  Works transparently on gzip-compressed files (the
    digest is over the compressed bytes on disk, which is what we want:
    the on-disk fingerprint changes whenever the file changes).
    """
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(_HASH_CHUNK)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def compute_cache_key(**params: Any) -> str:
    """Deterministic short cache key from Stage-1 parameters.

    Parameters
    ----------
    **params
        Arbitrary keyword arguments.  All values must be JSON-serialisable
        (``str`` fallback is applied for anything else).  Order is
        irrelevant: keys are sorted before hashing.

    Returns
    -------
    str
        Hex SHA-256 digest of the sorted, JSON-serialised parameter dict.
        Callers should include everything that affects the cached
        DataFrame content (raw file hash, version, thinning, trimming,
        liftover build, highlight thresholds, ...).  Anything left out
        risks a false cache hit after a parameter change.
    """
    payload = json.dumps(params, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


# ---------------------------------------------------------------------------
# TrackCache
# ---------------------------------------------------------------------------

@dataclass
class CachedTrack:
    """Handle to a cache-hit result returned by :meth:`TrackCache.get`."""

    df: pd.DataFrame
    leads: Optional[pd.DataFrame]
    pvals: Optional[np.ndarray]
    extras: dict[str, Any]  # arbitrary scalars: n_chroms, snp_count, ...


class TrackCache:
    """On-disk per-track cache for pycmplot Stage-1 output.

    Multi-panel safety
    ------------------
    Users routinely call the loader multiple times in one process (e.g.
    to render several panels of a multi-panel canvas), and different
    calls can legitimately share a track label — say ``"Hb"`` from
    Cameroon GWAS in panel 1 and ``"Hb"`` from Nigeria GWAS in panel 2.
    Metadata therefore keys on the full ``cache_key`` (not the label
    alone), and per-track filenames include the first 8 chars of the
    cache_key so two entries for the same label but different content
    coexist as separate files (``Hb.a1b2c3d4.parquet`` and
    ``Hb.e5f6g7h8.parquet``).  The label survives inside the metadata
    entry so cache dumps stay browsable.

    Not thread-safe: assume one process at a time.  A crashed process
    that leaves partial entries behind is handled cleanly on the next
    run (partial entries are treated as misses).
    """

    METADATA_FILENAME = "metadata.json"
    TRACKS_SUBDIR = "tracks"
    # Bump when the on-disk layout changes in a backwards-incompatible
    # way.  Older caches are transparently wiped on first access.
    CACHE_VERSION = 3
    KEY_SHORT_LEN = 8

    def __init__(self, cache_dir: str | os.PathLike, version: str) -> None:
        self.cache_dir = Path(cache_dir)
        self.tracks_dir = self.cache_dir / self.TRACKS_SUBDIR
        self.metadata_path = self.cache_dir / self.METADATA_FILENAME
        self.version = version
        self._metadata = self._load_metadata()

    # ---- public API -----------------------------------------------------

    def get(self, label: str, cache_key: str) -> Optional[CachedTrack]:
        """Return a :class:`CachedTrack` for (label, cache_key) if valid.

        Returns ``None`` when:
        * no metadata entry exists for *cache_key*,
        * the entry's ``status`` is not ``"complete"``,
        * the recorded ``label`` disagrees with the queried label
          (extremely unlikely — same key across different labels would
          require a SHA-256 collision), or
        * the on-disk parquet file has gone missing.

        A hit reads all sidecar files (leads, pvals) from disk and
        returns them alongside the main DataFrame.
        """
        entry = self._metadata.get("tracks", {}).get(cache_key)
        if entry is None:
            return None
        if entry.get("status") != "complete":
            return None
        if entry.get("label") not in (None, label):
            # Guard against accidental cross-label reuse.
            return None

        main_path = self.tracks_dir / entry.get("cache_file", f"{label}.parquet")
        if not main_path.exists():
            return None

        try:
            df = pd.read_parquet(main_path)
        except Exception as exc:  # pragma: no cover - corrupt parquet
            logger.warning("Cache read failed for %r (%s); treating as miss.",
                           label, exc)
            return None

        leads = None
        leads_name = entry.get("leads_file")
        if leads_name:
            leads_path = self.tracks_dir / leads_name
            if leads_path.exists():
                try:
                    leads = pd.read_parquet(leads_path)
                except Exception as exc:  # pragma: no cover
                    logger.warning("Cache leads read failed for %r (%s).",
                                   label, exc)

        pvals = None
        pvals_name = entry.get("pvals_file")
        if pvals_name:
            pvals_path = self.tracks_dir / pvals_name
            if pvals_path.exists():
                try:
                    # New format (0.4.0+): np.savez_compressed with a
                    # sorted float32 payload under key "P" (see put()
                    # for the compression rationale).  Older caches
                    # used raw ``.npy`` (uncompressed float64) or a
                    # short-lived parquet layout; both are read
                    # transparently.
                    name = str(pvals_path)
                    if name.endswith(".npz"):
                        with np.load(pvals_path) as npz:
                            pvals = np.asarray(npz["P"], dtype=float)
                    elif name.endswith(".parquet"):
                        pvals = pd.read_parquet(pvals_path)["P"].to_numpy()
                    else:  # .npy fallback
                        pvals = np.load(pvals_path)
                except Exception as exc:  # pragma: no cover
                    logger.warning("Cache pvals read failed for %r (%s).",
                                   label, exc)

        extras = dict(entry.get("extras", {}))
        logger.info("Cache HIT  %s  (%d rows)", label, len(df))
        return CachedTrack(df=df, leads=leads, pvals=pvals, extras=extras)

    def put(
        self,
        label: str,
        cache_key: str,
        df: pd.DataFrame,
        *,
        leads: Optional[pd.DataFrame] = None,
        pvals: Optional[np.ndarray] = None,
        extras: Optional[dict[str, Any]] = None,
        raw_source: Optional[str] = None,
    ) -> None:
        """Persist a completed Stage-1 result for the (label, cache_key)."""
        self.tracks_dir.mkdir(parents=True, exist_ok=True)
        short = cache_key[: self.KEY_SHORT_LEN]
        safe_label = _safe_filename_component(label)

        main_name = f"{safe_label}.{short}.parquet"
        _atomic_write(
            self.tracks_dir / main_name,
            lambda p: df.to_parquet(p),
        )

        leads_name: Optional[str] = None
        if leads is not None and not leads.empty:
            leads_name = f"{safe_label}.{short}.leads.parquet"
            _atomic_write(
                self.tracks_dir / leads_name,
                lambda p: leads.to_parquet(p),
            )

        pvals_name: Optional[str] = None
        if pvals is not None:
            # Compressed .npz with sorted float32 payload:
            #
            # * float32 preserves median(-log10 P) to ~9 decimal digits
            #   (well past any GWAS precision concern) and halves the
            #   raw byte footprint;
            # * pre-sorting improves DEFLATE compression by ~2x on
            #   heavy-tailed p-value distributions with no semantic
            #   loss (the QQ plotter sorts internally anyway);
            # * np.savez_compressed reaches ~29% of raw float64 size
            #   on typical GWAS scans -- much better than parquet's
            #   per-column overhead achieves on this data shape.
            pvals_name = f"{safe_label}.{short}.pvals.npz"
            pvals_sorted_f32 = np.sort(
                np.asarray(pvals, dtype=np.float32)
            )
            _atomic_write(
                self.tracks_dir / pvals_name,
                lambda p: np.savez_compressed(p, P=pvals_sorted_f32),
            )

        entry: dict[str, Any] = {
            "status": "complete",
            "label": label,
            "cache_key": cache_key,
            "cache_file": main_name,
            "n_rows": int(len(df)),
            "written": time.time(),
            "extras": dict(extras or {}),
        }
        if leads_name:
            entry["leads_file"] = leads_name
        if pvals_name:
            entry["pvals_file"] = pvals_name
        if raw_source:
            entry["raw_source"] = str(raw_source)

        self._metadata.setdefault("tracks", {})[cache_key] = entry
        self._save_metadata()

    def mark_failed(self, label: str, cache_key: str, reason: str) -> None:
        """Record a Stage-1 failure so a later resume can retry cleanly."""
        self._metadata.setdefault("tracks", {})[cache_key] = {
            "status": "failed",
            "label": label,
            "cache_key": cache_key,
            "reason": str(reason)[:500],
            "written": time.time(),
        }
        self._save_metadata()

    def status(self, cache_key: str) -> Optional[str]:
        """Return the current status for *cache_key* (or ``None`` if unknown)."""
        entry = self._metadata.get("tracks", {}).get(cache_key)
        return entry.get("status") if entry else None

    def entries_for_label(self, label: str) -> list[dict[str, Any]]:
        """Return all metadata entries whose recorded label matches.

        Useful for inspection: a multi-panel plot with two different
        ``"Hb"`` sumstats will have two entries here, each with its own
        ``cache_key``.
        """
        return [
            entry
            for entry in self._metadata.get("tracks", {}).values()
            if entry.get("label") == label
        ]

    # Filename patterns we consider "pycmplot-owned" under the two
    # subdirectories.  Anything not matching one of these is left
    # alone by :meth:`clear` so users who point ``--cache_dir`` at a
    # directory that also holds unrelated files can't lose data.
    _TRACKS_FILE_PATTERNS = (
        re.compile(r"^.+\.[0-9a-f]{4,64}\.parquet$"),                # main frame
        re.compile(r"^.+\.[0-9a-f]{4,64}\.leads\.parquet$"),         # leads
        re.compile(r"^.+\.[0-9a-f]{4,64}\.pvals\.(npz|parquet|npy)$"),  # pvals
        # Historical/legacy layouts we still recognise as ours.
        re.compile(r"^.+\.parquet$"),
        re.compile(r"^.+\.pvals\.(npz|parquet|npy)$"),
    )
    _ANNOTATIONS_FILE_PATTERNS = (
        re.compile(r"^hits\.[0-9a-f]{4,64}\.tsv$"),
        re.compile(r"^hits\.[0-9a-f]{4,64}\.meta\.json$"),
        # Legacy single-group layout, before 0.4.x
        re.compile(r"^hits\.tsv$"),
        re.compile(r"^hits\.meta\.json$"),
    )

    def clear(self) -> None:
        """Remove pycmplot cache artefacts and reset the in-memory metadata.

        **Safety model.**  This method never removes ``cache_dir``
        itself.  It targets only the three artefacts pycmplot writes:

        * ``metadata.json`` — the top-level manifest.
        * ``tracks/`` — per-track parquet + leads + pvals sidecars.
        * ``annotations/`` — group-scoped ``hits.<group>.tsv`` overlays
          and their ``.meta.json`` sidecars.

        Before deleting ``metadata.json`` we verify it *looks like a
        pycmplot metadata file* (has ``cache_version`` + ``tracks``
        keys, or is empty/JSON-invalid — both cases are safely ours to
        drop).  Inside each subdirectory we only delete files that
        match the pycmplot filename conventions (see
        :attr:`_TRACKS_FILE_PATTERNS` and :attr:`_ANNOTATIONS_FILE_PATTERNS`);
        anything foreign is left untouched and logged at INFO level.
        Empty subdirectories are then removed; non-empty ones (because
        the user dropped their own files inside) stay behind so
        nothing is silently lost.

        Motivation: an earlier release did
        ``shutil.rmtree(cache_dir)``, which would blow away
        *everything* under whatever the user pointed ``--cache_dir`` at
        — a serious footgun when the same directory happened to hold
        analysis notes, plot outputs, or other unrelated files.
        """
        cache_dir = self.cache_dir
        if not cache_dir.exists():
            self._metadata = self._empty_metadata()
            return

        # ---- 1. metadata.json ------------------------------------------
        if self.metadata_path.exists():
            _looks_ours = False
            try:
                _blob = json.loads(self.metadata_path.read_text() or "{}")
                if isinstance(_blob, dict) and (
                    "cache_version" in _blob or "tracks" in _blob
                    or not _blob  # empty dict was written by a very old build
                ):
                    _looks_ours = True
            except json.JSONDecodeError:
                # Corrupt / non-JSON — safe to treat as ours; the file
                # is inside our own subdir already.
                _looks_ours = True
            except Exception as exc:  # pragma: no cover — I/O race
                logger.warning(
                    "Could not inspect %s (%s); leaving it in place.",
                    self.metadata_path, exc,
                )
            if _looks_ours:
                try:
                    self.metadata_path.unlink()
                except Exception as exc:  # pragma: no cover
                    logger.warning(
                        "Failed to remove %s: %s",
                        self.metadata_path, exc,
                    )
            else:
                logger.info(
                    "Leaving %s in place — file does not look like a "
                    "pycmplot metadata manifest (missing 'cache_version' "
                    "and 'tracks' keys).", self.metadata_path,
                )

        # ---- 2. tracks/ + annotations/ ---------------------------------
        self._prune_subdir(
            self.tracks_dir,
            self._TRACKS_FILE_PATTERNS,
            label="tracks",
        )
        # Import lazily to avoid a top-of-file circular dep on the
        # module-level ``ANNOTATIONS_SUBDIR`` constant (defined later
        # in this file).
        annotations_dir = cache_dir / ANNOTATIONS_SUBDIR
        self._prune_subdir(
            annotations_dir,
            self._ANNOTATIONS_FILE_PATTERNS,
            label="annotations",
        )

        # ---- 3. Reset in-memory state ----------------------------------
        self._metadata = self._empty_metadata()

    @staticmethod
    def _prune_subdir(subdir: Path,
                      patterns: tuple,
                      label: str) -> None:
        """Delete files matching *patterns* under *subdir*, then rmdir if empty.

        Foreign files (anything not matching one of the pycmplot
        filename patterns) are left in place and the user is told at
        INFO level.  A non-empty subdir at the end is kept — we never
        force-remove a directory the user might have added their own
        files to.
        """
        if not subdir.exists() or not subdir.is_dir():
            return
        _kept: list[str] = []
        for entry in subdir.iterdir():
            if not entry.is_file():
                # Nested directory the user put here — untouched.
                _kept.append(entry.name)
                continue
            if any(pat.match(entry.name) for pat in patterns):
                try:
                    entry.unlink()
                except Exception as exc:  # pragma: no cover
                    logger.warning("Failed to remove %s: %s", entry, exc)
                    _kept.append(entry.name)
            else:
                _kept.append(entry.name)
        if _kept:
            logger.info(
                "Kept %d non-pycmplot file(s) under %s/: %s",
                len(_kept), label,
                ", ".join(sorted(_kept)[:6])
                + (" …" if len(_kept) > 6 else ""),
            )
            return
        # Directory is now empty — drop it.
        try:
            subdir.rmdir()
        except Exception as exc:  # pragma: no cover
            logger.warning("Failed to rmdir %s: %s", subdir, exc)

    # ---- metadata helpers ----------------------------------------------

    def _empty_metadata(self) -> dict[str, Any]:
        return {
            "cache_version": self.CACHE_VERSION,
            "pycmplot_version": self.version,
            "created": time.time(),
            "tracks": {},
        }

    def _load_metadata(self) -> dict[str, Any]:
        if not self.metadata_path.exists():
            return self._empty_metadata()
        try:
            data = json.loads(self.metadata_path.read_text())
        except Exception as exc:
            logger.warning(
                "Cache metadata unreadable (%s); starting fresh.", exc,
            )
            return self._empty_metadata()

        # Bail out on any cache-layout version mismatch.  Older caches
        # used label-keyed metadata (cache_version 1 / unset) which is
        # incompatible with the current cache_key-keyed layout — trying
        # to read them silently would return misleading misses forever.
        seen_version = data.get("cache_version")
        if seen_version != self.CACHE_VERSION:
            logger.info(
                "Cache layout changed (found version %r, expected %r); "
                "resetting %s.",
                seen_version, self.CACHE_VERSION, self.cache_dir,
            )
            # Wipe on-disk data too so stale parquet files don't linger.
            if self.tracks_dir.exists():
                try:
                    shutil.rmtree(self.tracks_dir)
                except OSError as exc:
                    logger.warning("Cache reset: could not remove %s (%s).",
                                   self.tracks_dir, exc)
            return self._empty_metadata()

        # Preserve unrelated top-level keys but ensure structure exists.
        data.setdefault("tracks", {})
        data.setdefault("pycmplot_version", self.version)
        return data

    def _save_metadata(self) -> None:
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._metadata["cache_version"] = self.CACHE_VERSION
        self._metadata["pycmplot_version"] = self.version
        self._metadata["updated"] = time.time()
        payload = json.dumps(self._metadata, indent=2, sort_keys=True)
        _atomic_write(self.metadata_path, lambda p: p.write_text(payload))


def _safe_filename_component(name: str) -> str:
    """Sanitise a label into a filename-safe component.

    Track labels can contain arbitrary characters (users have been seen
    to pass e.g. ``"HbF (Kenya)"``), but not every filesystem tolerates
    ``/``, ``\\``, ``:`` etc.  We strip those to underscores and cap
    the result to a reasonable length; uniqueness is guaranteed
    elsewhere by the cache-key suffix.
    """
    if not name:
        return "unnamed"
    safe = re.sub(r"[^A-Za-z0-9._+-]+", "_", name)
    return safe[:64] or "unnamed"


# ---------------------------------------------------------------------------
# Hits table cache with user-editable overlay
# ---------------------------------------------------------------------------
#
# The hits table is a small (usually <500-row) DataFrame produced by
# ``pycmplot.annotation.get_hits_summary_table``.  Regenerating it costs a
# few seconds (Ensembl GFF3 load + per-lead nearest-gene lookup), which
# is a noticeable fraction of the warm-cache runtime.  Beyond the speed
# win, exposing the cached table as a **user-editable TSV** turns the
# annotation layer into a hand-tunable overlay:
#
#   * Users can hand-add rows (custom loci, curated names, meta-analysis
#     hits absent from the input files) by editing the TSV directly.
#   * Users can hand-edit auto-generated rows (override auto-picked gene
#     name, correct biotype, suppress a false-positive lead).
#   * User-authored rows survive any invalidation because they carry
#     ``source="user"`` and the regeneration pass only touches rows
#     tagged ``source="auto"``.
#
# The overlay lives at ``<cache_dir>/annotations/hits.tsv`` and starts
# with a header comment describing the schema.

ANNOTATIONS_SUBDIR = "annotations"
HITS_FILENAME_STEM = "hits"       # actual filenames get a .<group>.tsv suffix
SOURCE_COL = "source"
AUTO_TAG = "auto"
USER_TAG = "user"
HITS_GROUP_SHORT_LEN = 8

HITS_HEADER_COMMENT = (
    "# pycmplot hits overlay -- edit rows with source='user' to add or\n"
    "# override annotations; auto-generated rows (source='auto') are\n"
    "# regenerated whenever the underlying leads / signif_threshold / GFF3\n"
    "# resource change.  Blank-line-separated blocks are preserved as-is.\n"
    "# Required columns depend on the plot mode: at minimum CHR, POS, and\n"
    "# whichever column is used for label text (SNP by default; --annotate\n"
    "# GENE for gene-name annotation).\n"
    "#\n"
    "# Per-locus highlight color: the 'highlight_color' column takes any\n"
    "# matplotlib-recognised color (name, hex, rgb tuple).  The sentinel\n"
    "# 'auto' -- or a blank / NaN -- falls back to the plot-time\n"
    "# --highlight_color argument.  Invalid values are rejected with a\n"
    "# warning and treated as 'auto'.  Edit specific rows to give\n"
    "# individual loci their own color (e.g. grey out an MHC row, red for\n"
    "# a novel hit).\n"
    "#\n"
    "# Per-locus category: the 'category' column gives each locus a\n"
    "# label that drives an optional custom legend (e.g. 'novel',\n"
    "# 'replicated', 'MHC').  When every row is left at the default\n"
    "# 'significant' AND every highlight_color is 'auto', no legend is\n"
    "# added.  As soon as any row is edited (colour and/or category),\n"
    "# the plotter renders a legend with one entry per unique category\n"
    "# in first-appearance order -- reorder rows in this file to\n"
    "# reorder the legend.\n"
    "#\n"
    "# This file is scoped to ONE loader call's group of tracks.  A\n"
    "# different loader call with a different set of tracks writes to a\n"
    "# separate hits.<group_key>.tsv in this directory, so multi-panel\n"
    "# workflows don't clobber each other's annotations.\n"
)


def compute_hits_group_key(track_cache_keys: list[str]) -> str:
    """Deterministic short key naming a group of tracks in one loader call.

    The key is derived from the sorted set of **per-track cache keys**
    supplied to the loader.  Because each cache_key already covers the
    raw file SHA-256, the pycmplot version, and every Stage-1 parameter
    that affects the cached data (``trim_pval``, ``auto_thin*``,
    ``highlight*``, ``signif_threshold``, ``build``, ``logp``), the
    group key naturally scopes the hits overlay to a specific
    ``(files, parameters)`` combination.  This mirrors the per-track
    cache invalidation model:

    * Two loader calls that share ``cache_dir`` but supply different
      raw files (multi-panel pattern) get different group keys, so
      their hits overlays land in separate files and don't clobber
      each other.
    * Two loader calls with the same raw files but different Stage-1
      parameters (e.g. ``highlight_thresh=5e-8`` vs ``1e-6``) also get
      different group keys, so each parameter setting has its own
      overlay.  This matches the reproducibility mental model that
      "different parameters → different analysis → different cached
      artefacts".  User-authored rows are preserved *within* one
      ``(files, parameters)`` group across warm re-runs; orphan
      overlays from earlier parameter settings are left on disk (they
      are not consulted, but a user who wants to carry a hand-edit
      forward can copy rows manually).
    """
    if not track_cache_keys:
        return "default"
    payload = "|".join(sorted(str(k) for k in track_cache_keys)).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()[:HITS_GROUP_SHORT_LEN]


def _hits_paths(cache_dir: str | os.PathLike,
                group_key: str) -> tuple[Path, Path]:
    base = Path(cache_dir) / ANNOTATIONS_SUBDIR
    fname = f"{HITS_FILENAME_STEM}.{group_key}"
    return base / f"{fname}.tsv", base / f"{fname}.meta.json"


def hits_auto_key(leads_df: pd.DataFrame, resources_signature: str,
                  window_kb: int = 500) -> str:
    """Cache key for the auto-generated hits table.

    Includes a fingerprint of the leads DataFrame (rows the annotation
    pass would consume) plus the GFF3 resource fingerprint and the
    window size, so any real content or resource change invalidates.
    """
    if leads_df is None or leads_df.empty:
        leads_sig = "empty"
    else:
        # Hash of a stable CSV projection — fast at O(n_leads), which is
        # tiny compared to the annotation pass itself.
        try:
            core = leads_df[["CHR", "POS", "P"]].sort_values(["CHR", "POS"])
        except KeyError:
            core = leads_df.sort_index()
        blob = core.to_csv(index=False).encode("utf-8")
        leads_sig = hashlib.sha256(blob).hexdigest()
    return compute_cache_key(
        leads_sig=leads_sig,
        resources=resources_signature,
        window_kb=int(window_kb),
    )


def read_hits_overlay(cache_dir: str | os.PathLike,
                      group_key: str
                      ) -> tuple[Optional[pd.DataFrame], Optional[dict]]:
    """Read the hits overlay TSV and its sidecar metadata, if present.

    Parameters
    ----------
    cache_dir : path-like
        Cache root directory.
    group_key : str
        Short group identifier from :func:`compute_hits_group_key`.
        Different sets of tracks land in different files so multi-panel
        callers don't clobber each other's overlays.

    Returns
    -------
    (df, meta)
        Either or both may be ``None`` when the files don't exist.
        Lines beginning with ``#`` are treated as comments (compatible
        with ``pd.read_csv(comment='#')``).
    """
    hits_path, meta_path = _hits_paths(cache_dir, group_key)
    df = None
    if hits_path.exists():
        try:
            # Skip only *header* lines beginning with ``#`` — passing
            # ``comment='#'`` to pd.read_csv would strip mid-field ``#``
            # characters too, which breaks legitimate hex-color values
            # like ``#00cc44`` in the ``highlight_color`` column.
            import io as _io
            with open(hits_path, "r", encoding="utf-8") as fh:
                lines = []
                header_done = False
                for line in fh:
                    if not header_done and line.lstrip().startswith("#"):
                        continue
                    header_done = True
                    lines.append(line)
            df = pd.read_csv(_io.StringIO("".join(lines)), sep="\t")
            if SOURCE_COL not in df.columns:
                df[SOURCE_COL] = AUTO_TAG
        except Exception as exc:
            logger.warning("Hits overlay unreadable (%s); ignoring.", exc)
            df = None
    meta: Optional[dict] = None
    if meta_path.exists():
        try:
            meta = json.loads(meta_path.read_text())
        except Exception as exc:
            logger.warning("Hits overlay metadata unreadable (%s).", exc)
            meta = None
    return df, meta


def write_hits_overlay(
    cache_dir: str | os.PathLike,
    auto_df: pd.DataFrame,
    *,
    auto_key: str,
    group_key: str,
    preserve_user_from: Optional[pd.DataFrame] = None,
) -> Path:
    """Write the merged hits overlay to disk atomically.

    The output has ``source`` as the first column so users can spot
    auto- vs user-rows at a glance.  If *preserve_user_from* is
    supplied, any rows tagged ``source="user"`` there are appended
    below the auto rows (source order preserved), so user edits carry
    over across regenerations.

    ``group_key`` is included in the filename so multi-panel workflows
    that share a ``cache_dir`` but supply different track sets get
    separate overlays — see :func:`compute_hits_group_key`.
    """
    hits_path, meta_path = _hits_paths(cache_dir, group_key)
    hits_path.parent.mkdir(parents=True, exist_ok=True)

    auto = auto_df.copy() if auto_df is not None else pd.DataFrame()
    if not auto.empty:
        auto = auto.copy()
        auto[SOURCE_COL] = AUTO_TAG
        # Ensure the highlight_color + category columns exist with
        # sentinel defaults so users can hand-edit specific rows.  Both
        # columns are inherited from the previous overlay by
        # ``(CHR, POS)`` lookup below, so a user edit to an auto row
        # survives regeneration without them having to also change
        # ``source`` to 'user'.
        from pycmplot.annotation import (
            HIGHLIGHT_COLOR_COL, HIGHLIGHT_COLOR_AUTO,
            CATEGORY_COL, CATEGORY_DEFAULT,
        )
        if HIGHLIGHT_COLOR_COL not in auto.columns:
            auto[HIGHLIGHT_COLOR_COL] = HIGHLIGHT_COLOR_AUTO
        if CATEGORY_COL not in auto.columns:
            auto[CATEGORY_COL] = CATEGORY_DEFAULT

        # Inherit user-set colours AND categories from the previous
        # overlay for matching (CHR, POS) rows.  This lets an analyst
        # hand-edit either column on an auto row (without also
        # changing ``source='auto'`` to 'user') and have that edit
        # survive the next regeneration.  Only values that differ
        # from the respective sentinels are inherited, so a user
        # who wanted to reset an edit can just put 'auto' /
        # 'significant' back.
        if (preserve_user_from is not None
                and not preserve_user_from.empty
                and {"CHR", "POS"}.issubset(preserve_user_from.columns)
                and {"CHR", "POS"}.issubset(auto.columns)):
            # Normalise dtypes before merge -- the cached auto DF
            # stores CHR as a categorical, while the TSV-parsed
            # overlay reads it back as int64/object depending on
            # content; pandas refuses to merge across those.
            def _norm_keys(df: "pd.DataFrame") -> "pd.DataFrame":
                out = df.copy()
                out["CHR"] = out["CHR"].astype(str)
                out["POS"] = out["POS"].astype("int64")
                return out

            auto_key_df = _norm_keys(auto[["CHR", "POS"]])

            for col, sentinels in [
                (HIGHLIGHT_COLOR_COL,
                 {HIGHLIGHT_COLOR_AUTO, "", "nan", "none", "na"}),
                (CATEGORY_COL,
                 {CATEGORY_DEFAULT, "", "nan", "none", "na"}),
            ]:
                if col not in preserve_user_from.columns:
                    continue
                prev = _norm_keys(preserve_user_from[["CHR", "POS", col]])
                prev[col] = prev[col].astype(str)
                # Filter out sentinel / blank rows so we don't overwrite
                # a fresh non-default with a stale sentinel.
                prev = prev[
                    ~prev[col].str.strip().str.lower().isin(sentinels)
                ]
                if prev.empty:
                    continue
                # Deduplicate on (CHR, POS) — keep the last user edit.
                prev = prev.drop_duplicates(
                    subset=["CHR", "POS"], keep="last",
                )
                inherited = auto_key_df.merge(
                    prev, on=["CHR", "POS"], how="left",
                )
                mask = inherited[col].notna()
                auto.loc[mask.values, col] = \
                    inherited.loc[mask, col].values

    user_rows: Optional[pd.DataFrame] = None
    if preserve_user_from is not None and not preserve_user_from.empty:
        if SOURCE_COL in preserve_user_from.columns:
            user_rows = preserve_user_from[
                preserve_user_from[SOURCE_COL].astype(str) == USER_TAG
            ].copy()
        if user_rows is not None and user_rows.empty:
            user_rows = None

    if user_rows is not None:
        # Align columns: use union so hand-added rows can carry extra
        # columns (e.g. a "notes" column the user adds).
        all_cols = list(auto.columns) + [
            c for c in user_rows.columns if c not in auto.columns
        ]
        auto = auto.reindex(columns=all_cols)
        user_rows = user_rows.reindex(columns=all_cols)
        merged = pd.concat([auto, user_rows], ignore_index=True)
    else:
        merged = auto

    # Move ``source`` to the front for readability.
    if SOURCE_COL in merged.columns:
        cols = [SOURCE_COL] + [c for c in merged.columns if c != SOURCE_COL]
        merged = merged[cols]

    def _write(p: Path) -> None:
        with open(p, "w", encoding="utf-8") as fh:
            fh.write(HITS_HEADER_COMMENT)
            merged.to_csv(fh, sep="\t", index=False, na_rep='NA')

    _atomic_write(hits_path, _write)

    meta = {
        "auto_key": auto_key,
        "group_key": group_key,
        "written": time.time(),
        "n_auto": int((merged[SOURCE_COL] == AUTO_TAG).sum())
                  if SOURCE_COL in merged.columns else int(len(merged)),
        "n_user": int((merged[SOURCE_COL] == USER_TAG).sum())
                  if SOURCE_COL in merged.columns else 0,
    }
    _atomic_write(meta_path,
                  lambda p: p.write_text(json.dumps(meta, indent=2, sort_keys=True)))
    return hits_path


def resources_fingerprint(resources) -> str:
    """Short fingerprint of the annotation resource files.

    Combines the file paths *and* their sizes so a swap of the bundled
    GFF3 for a fresh release invalidates the cache.  Full SHA-256 of the
    ~40 MB gene-info file would be more principled but is disproportionate
    given the annotation pass itself is only a few seconds; size + path
    is a good practical compromise.
    """
    from pycmplot.resources import ResourceConfig, default_resources
    r = resources or default_resources
    parts: list[str] = []
    for attr in ("geneinfo_hg38", "geneinfo_hg19"):
        try:
            path = getattr(r, attr, None)
            if path is None:
                parts.append(f"{attr}:none")
                continue
            p = Path(path)
            if p.exists():
                parts.append(f"{attr}:{p.name}:{p.stat().st_size}")
            else:
                parts.append(f"{attr}:{p.name}:missing")
        except Exception:
            parts.append(f"{attr}:err")
    return "|".join(parts)


# ---------------------------------------------------------------------------
# Atomic write helper
# ---------------------------------------------------------------------------

def _atomic_write(target: Path, writer) -> None:
    """Write to *target* atomically via a temp file + rename.

    *writer* is a callable ``writer(path)`` that persists the payload to
    the supplied temporary path.  The temp file is renamed over *target*
    on success (rename is atomic on POSIX and Windows for same-directory
    moves), and cleaned up on failure so we never leave garbage next to
    the real cache file.

    Temp-file naming
    ----------------
    The temp path is constructed to preserve the original extension so
    library writers that auto-append an extension (notably
    :func:`numpy.save`, which appends ``.npy`` when the path doesn't
    already end in ``.npy``) still write to the exact path we expect.
    Concretely, ``foo.pvals.npy`` becomes ``foo.pvals.__tmp__.npy`` —
    not ``foo.pvals.npy.tmp`` (which used to yield
    ``foo.pvals.npy.tmp.npy`` on disk and a ``FileNotFoundError`` on the
    subsequent :func:`os.replace`).
    """
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    # Insert ``.__tmp__`` just before the final extension so any
    # auto-append behaviour in the writer sees the same extension the
    # target has.
    if target.suffix:
        tmp = target.with_name(target.stem + ".__tmp__" + target.suffix)
    else:
        tmp = target.with_name(target.name + ".__tmp__")
    try:
        if isinstance(writer, (str, bytes, bytearray)):
            tmp.write_text(writer if isinstance(writer, str)
                           else writer.decode("utf-8"))
        else:
            writer(tmp)
        os.replace(tmp, target)
    finally:
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass
