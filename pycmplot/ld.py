"""
pycmplot.ld
============

CSR-encoded sparse graph representation of pairwise LD (r²) from an
external reference panel (e.g. a PLINK ``--r2`` ``.ld`` file).

Design rationale
----------------

The naive ``{snp_a: {snp_b: r²}}`` dict-of-dict is convenient but
carries ~1 kB of Python-object overhead per entry.  At the scale of
1 M SNPs × 10 neighbours on average that's ~10 GB of heap for what
should be a ~240 MB numpy structure.  Compressed Sparse Row (CSR)
storage collapses this by:

* filtering weak edges (r² < ``threshold``) at build time — the
  saved graph then satisfies ``absence in graph == independent``
  by construction, so query paths need no threshold comparison;
* storing per-SNP neighbour lists as contiguous ``np.int32``
  slices of a single ``indices`` array — cache-friendly scans,
  no hashing on lookup;
* keeping the SNP name table separate from the graph so per-edge
  memory is exactly 8 bytes (int32 neighbour + float32 r²).

Reference-panel LD graphs built once at pipeline setup can be
serialised (``.npz``) and reloaded across runs in milliseconds.

The primary consumer is :func:`pycmplot.stats.clump`, which uses
:meth:`LDGraph.any_linked` for its greedy-clumping inner loop, but
the graph also supports :meth:`LDGraph.neighbors` for downstream
regional-plot-style LD-tier colouring.
"""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class LDGraph:
    """CSR-encoded sparse LD reference.

    Attributes
    ----------
    snp_names : numpy.ndarray of str
        Node ``i`` in the graph corresponds to ``snp_names[i]``.
    snp_index : dict[str, int]
        Inverse lookup, ``snp_name -> node_index``.
    snp_chrom, snp_pos : numpy.ndarray, shape (N,)
        Per-node genomic position key used by
        :func:`pycmplot.variant_matcher.ld_id_map` to harmonise IDs.
        ``snp_chrom`` holds normalised chromosome names (``"1"`` ..
        ``"22"``, ``"X"``, ``"Y"``, ``"XY"``, ``"M"``; ``""`` when
        unknown) and ``snp_pos`` int32 base-pair positions (``-1`` when
        unknown).  Filled from the ``CHR``/``BP`` columns at build time
        when available, otherwise parsed once from the SNP names; saved
        with the graph so the work is never repeated.  Graphs saved
        before this field existed compute it lazily on first access.
    indptr : numpy.ndarray of int32, shape (N+1,)
        Standard CSR row-pointer array.  Neighbours of node ``i``
        live in ``indices[indptr[i]:indptr[i+1]]``.
    indices : numpy.ndarray of int32, shape (M,)
        Neighbour node indices, sorted ascending within each node's
        slice so ``searchsorted`` can find pairs in O(log k).
    data : numpy.ndarray of float32, shape (M,)
        r² value on each edge, aligned with ``indices``.
    threshold : float
        The r² threshold applied at build time.  Edges with
        ``r² < threshold`` are absent from the graph.
    """

    def __init__(self, snp_names: np.ndarray, indptr: np.ndarray,
                 indices: np.ndarray, data: np.ndarray,
                 threshold: float, has_r2: bool = True,
                 build: Optional[str] = None,
                 snp_chrom: Optional[np.ndarray] = None,
                 snp_pos: Optional[np.ndarray] = None):
        self.snp_names = np.asarray(snp_names)
        self.indptr = np.asarray(indptr, dtype=np.int32)
        self.indices = np.asarray(indices, dtype=np.int32)
        self.data = np.asarray(data, dtype=np.float32)
        self.threshold = float(threshold)
        # False when the graph was built from a reduced-format .ld
        # (no R2 column) — every ``data`` value is a sentinel 1.0.
        # Downstream regional-plot code should refuse to run when
        # this is False.
        self.has_r2 = bool(has_r2)
        # Genome build the reference panel was computed on.  Stored
        # on the graph so ``pycmplot.io.load`` can sanity-check that
        # the sumstats / GeneHancer coordinate system agrees with
        # the LD reference before clumping.  ``None`` = unknown
        # (older graphs or user chose not to declare).
        self.build: Optional[str] = self._normalise_build(build)
        self.snp_index: dict[str, int] = {
            str(name): int(i) for i, name in enumerate(self.snp_names)
        }
        # Position key (see class docstring).  Both arrays or neither;
        # a length mismatch is a corrupt input, not a lazy-fill case.
        self._snp_chrom: Optional[np.ndarray] = None
        self._snp_pos: Optional[np.ndarray] = None
        if snp_chrom is not None and snp_pos is not None:
            _c = np.asarray(snp_chrom).astype("U")
            _p = np.asarray(snp_pos, dtype=np.int32)
            if _c.shape[0] != self.snp_names.shape[0] or _p.shape[0] != self.snp_names.shape[0]:
                raise ValueError(
                    f"snp_chrom/snp_pos length ({_c.shape[0]}/{_p.shape[0]}) "
                    f"does not match snp_names ({self.snp_names.shape[0]})."
                )
            self._snp_chrom, self._snp_pos = _c, _p

    # ------------------------------------------------------------------
    # Position key (built once, persisted by save())
    # ------------------------------------------------------------------

    @staticmethod
    def _positions_from_names(names: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Vectorised ``CHR<sep>POS`` parse of SNP names.

        Names that do not start with a structured position (rsIDs,
        free-form labels) get ``("", -1)``.
        """
        from pycmplot.variant_matcher import VariantMatcher, _CP_HEAD
        n = int(np.asarray(names).shape[0])
        if n == 0:
            return np.array([], dtype="U1"), np.array([], dtype=np.int32)
        ext = pd.Series(np.asarray(names), dtype=str).str.extract(_CP_HEAD)
        chr_map = {u: VariantMatcher._normalise_chr(u)
                   for u in ext[0].dropna().unique()}
        chrom = ext[0].map(chr_map).fillna("").to_numpy(dtype="U")
        pos = pd.to_numeric(ext[1], errors="coerce").fillna(-1)
        return chrom, pos.to_numpy(dtype=np.int32)

    def _ensure_positions(self) -> None:
        if self._snp_chrom is None or self._snp_pos is None:
            logger.info(
                "LDGraph: no stored position key; parsing %s SNP names "
                "once (re-save the graph to persist it).",
                f"{self.n_snps:,}",
            )
            self._snp_chrom, self._snp_pos = self._positions_from_names(self.snp_names)

    @property
    def snp_chrom(self) -> np.ndarray:
        self._ensure_positions()
        return self._snp_chrom

    @property
    def snp_pos(self) -> np.ndarray:
        self._ensure_positions()
        return self._snp_pos

    @staticmethod
    def _normalise_build(build: Optional[str]) -> Optional[str]:
        if build is None:
            return None
        s = str(build).strip().lower()
        _MAP = {
            "hg38": "hg38", "grch38": "hg38", "38": "hg38",
            "hg19": "hg19", "grch37": "hg19", "b37": "hg19", "19": "hg19",
            "hg18": "hg18", "ncbi36": "hg18", "18": "hg18",
        }
        if s not in _MAP:
            raise ValueError(
                f"Unsupported LD reference build {build!r}; "
                f"supported: hg18, hg19, hg38."
            )
        return _MAP[s]

    # ------------------------------------------------------------------
    # Introspection helpers
    # ------------------------------------------------------------------

    @property
    def n_snps(self) -> int:
        return int(self.snp_names.shape[0])

    @property
    def n_edges(self) -> int:
        """Number of directed edges (each pair is counted twice)."""
        return int(self.indices.shape[0])

    def __len__(self) -> int:
        return self.n_snps

    def __contains__(self, snp: str) -> bool:
        return str(snp) in self.snp_index

    def __repr__(self) -> str:
        return (f"<LDGraph n_snps={self.n_snps:,} "
                f"n_edges={self.n_edges:,} "
                f"threshold={self.threshold}>")

    # ------------------------------------------------------------------
    # Query API
    # ------------------------------------------------------------------

    def r2(self, snp_a: str, snp_b: str) -> float:
        """Return r² between *snp_a* and *snp_b*, or 0.0 if the pair
        is not in the graph (i.e. below the build-time threshold or
        the SNP is missing from the reference panel).
        """
        ia = self.snp_index.get(str(snp_a))
        if ia is None:
            return 0.0
        ib = self.snp_index.get(str(snp_b))
        if ib is None:
            return 0.0
        s, e = int(self.indptr[ia]), int(self.indptr[ia + 1])
        if s == e:
            return 0.0
        row_indices = self.indices[s:e]
        j = int(np.searchsorted(row_indices, ib))
        if j < row_indices.shape[0] and int(row_indices[j]) == ib:
            return float(self.data[s + j])
        return 0.0

    def neighbors(self, snp: str,
                  r2_min: Optional[float] = None) -> list[tuple[str, float]]:
        """Return ``[(snp, r²)]`` for every LD partner of *snp*.

        Passing *r2_min* filters the return set to r² >= r2_min
        (must be >= the graph's build-time threshold).  Neighbours
        are returned in ascending node-index order (deterministic
        but not r²-sorted).
        """
        ia = self.snp_index.get(str(snp))
        if ia is None:
            return []
        s, e = int(self.indptr[ia]), int(self.indptr[ia + 1])
        if s == e:
            return []
        names = self.snp_names[self.indices[s:e]]
        r2s = self.data[s:e]
        if r2_min is not None:
            mask = r2s >= float(r2_min)
            names = names[mask]
            r2s = r2s[mask]
        return [(str(n), float(r)) for n, r in zip(names, r2s)]

    def any_linked(self, snp: str, others: Iterable[str],
                   r2_min: Optional[float] = None) -> bool:
        """Fast batch test: is *snp* linked to ANY member of *others*
        at r² >= *r2_min* (defaults to the graph's threshold)?

        Used by :func:`pycmplot.stats.clump` in its inner loop.
        """
        ia = self.snp_index.get(str(snp))
        if ia is None:
            return False
        s, e = int(self.indptr[ia]), int(self.indptr[ia + 1])
        if s == e:
            return False
        row_indices = self.indices[s:e]
        row_data = self.data[s:e]
        thresh = self.threshold if r2_min is None else float(r2_min)
        for other in others:
            ib = self.snp_index.get(str(other))
            if ib is None:
                continue
            j = int(np.searchsorted(row_indices, ib))
            if (j < row_indices.shape[0]
                    and int(row_indices[j]) == ib
                    and float(row_data[j]) >= thresh):
                return True
        return False

    # ------------------------------------------------------------------
    # Cross-source relabelling
    # ------------------------------------------------------------------

    def remap_ids(self, rename: dict) -> None:
        """Relabel graph nodes in-place using ``rename[old] -> new``.

        Nodes not in *rename* keep their existing name.  The graph's
        edge topology and weights are untouched — only the SNP-name
        table and the ``snp_index`` reverse map are rebuilt.  The
        position key (``snp_chrom``/``snp_pos``) is kept as is, since
        node order and physical positions do not change.  Used
        primarily to harmonise the graph's SNP naming with an
        external DataFrame's naming via
        :class:`pycmplot.variant_matcher.VariantMatcher`.
        """
        new_names = np.array(
            [str(rename.get(str(n), str(n))) for n in self.snp_names],
            dtype="U",
        )
        self._ensure_positions()  # pin positions to the pre-rename names
        self.snp_names = new_names
        self.snp_index = {str(n): int(i) for i, n in enumerate(new_names)}

    # ------------------------------------------------------------------
    # Serialisation
    # ------------------------------------------------------------------

    def save(self, path: str | Path, *,
             compact: bool = False,
             symmetric_only: Optional[bool] = None,
             quantize_r2: Optional[bool] = None,
             compression: Optional[str] = None,
             compression_level: Optional[int] = None) -> None:
        """Save as a ``.npz`` (with metadata sidecar).

        Layouts:

        * **Default (dense + float32 + zlib)** — every pair stored
          in both directions, r² as float32, zlib-compressed.  Portable
          and self-describing via the JSON sidecar.
        * **Compact** (``compact=True``) — enables symmetric-only
          storage (canonical A→B where idx(A)<idx(B); reverse edges
          recomputed at load time) plus uint16 quantisation of r²
          (65,536 levels, precision ~1.5e-5).  Compression stays at
          zlib by default because on typical numeric LD arrays zlib
          matches or beats zstd at moderate levels.

        Compression:

        * ``compression='zlib'`` (default) — numpy's own ``.npz``.
        * ``compression='zstd'`` (opt-in) — writes a ``.ldz``
          container using the ``zstandard`` library.  Falls back
          to zlib with a warning if zstandard isn't installed.
          Only worth using at ``compression_level >= 15``; at
          default level ~10 it typically loses to zlib on LD data.
        * ``compression_level`` — passed through to zstd only
          (ignored for zlib).  ``19`` gives ~15% smaller than zlib
          but is ~13× slower to write; useful for build-once,
          load-many workflows.

        Individual flags override ``compact``: pass
        ``symmetric_only=False`` to keep dense storage etc.
        """
        # Resolve individual flags from the compact convenience.
        if symmetric_only is None:
            symmetric_only = compact
        if quantize_r2 is None:
            quantize_r2 = compact
        if compression is None:
            compression = "zlib"  # zlib wins for numeric LD arrays
                                   # at moderate zstd levels; users
                                   # opt into zstd explicitly for
                                   # max compression at high levels.
        if compression_level is None:
            compression_level = 19  # only used when compression='zstd'

        path = Path(path)

        # ---- 1. Symmetric-only reduction --------------------------------
        # Keep only canonical edges (src < dst).  On load we expand
        # back to dense so all in-memory queries stay O(log k).
        if symmetric_only and self.indices.size > 0:
            _keep_mask = np.empty(self.indices.size, dtype=bool)
            _n = self.snp_names.shape[0]
            # For each row src, the neighbours are indices[indptr[src]:indptr[src+1]].
            # An edge src->dst is canonical iff src < dst.
            for src in range(_n):
                s, e = int(self.indptr[src]), int(self.indptr[src + 1])
                _keep_mask[s:e] = self.indices[s:e] > src
            _indices_c = self.indices[_keep_mask]
            _data_c = self.data[_keep_mask]
            # Rebuild indptr from the filtered edges.
            _sources = np.repeat(np.arange(_n, dtype=np.int32),
                                 np.diff(self.indptr))
            _sources_c = _sources[_keep_mask]
            _counts_c = np.bincount(_sources_c, minlength=_n)
            _indptr_c = np.concatenate([[0], np.cumsum(_counts_c)]).astype(np.int32)
        else:
            _indices_c = self.indices
            _data_c = self.data
            _indptr_c = self.indptr

        # ---- 2. r² quantisation to uint16 -------------------------------
        # 0.0 -> 0, 1.0 -> 65535; reconstruct as float32 / 65535 on load.
        # Precision ~1.5e-5 is trivial vs. the underlying LD estimation
        # noise (bootstrap CIs of r² are typically ±0.01 at N=1000).
        if quantize_r2:
            _data_saved = np.clip(np.round(_data_c * 65535.0), 0, 65535).astype(np.uint16)
        else:
            _data_saved = _data_c.astype(np.float32)

        arrays = {
            "indptr": _indptr_c,
            "indices": _indices_c,
            "data": _data_saved,
            "snp_names": self.snp_names,
            "threshold": np.float32(self.threshold),
            "has_r2": np.bool_(self.has_r2),
            "symmetric_only": np.bool_(bool(symmetric_only)),
            "quantized_r2": np.bool_(bool(quantize_r2)),
            # Build stored as a short unicode string so .npz can round-
            # trip with allow_pickle=False.  Empty string = unset.
            "build": np.array(self.build or "", dtype="U8"),
            # Position key: computed once (at build or here) and stored
            # so harmonised clumping never re-parses SNP names.
            "snp_chrom": self.snp_chrom,
            "snp_pos": self.snp_pos,
        }

        # ---- 3. Container / compression --------------------------------
        used_compression = compression
        if compression == "zstd":
            try:
                import zstandard  # noqa: F401  (checked for availability)
                _write_ldz(path, arrays, level=compression_level)
            except ImportError:
                logger.warning(
                    "zstandard not installed; falling back to zlib "
                    "(.npz).  ``pip install zstandard`` for the "
                    "``.ldz`` option."
                )
                used_compression = "zlib"
                np.savez_compressed(path, **arrays)
        else:
            np.savez_compressed(path, **arrays)

        meta = {
            "version": 3,          # 3: adds snp_chrom / snp_pos
            "n_snps": self.n_snps,
            "n_edges": self.n_edges,
            "threshold": self.threshold,
            "has_r2": self.has_r2,
            "symmetric_only": bool(symmetric_only),
            "quantized_r2": bool(quantize_r2),
            "compression": used_compression,
            "build": self.build,
        }
        path.with_suffix(path.suffix + ".json").write_text(
            json.dumps(meta, indent=2)
        )
        logger.info(
            "Saved LDGraph -> %s (%s SNPs, %s edges [%s in-file], "
            "threshold %.2f; symmetric_only=%s quantized_r2=%s "
            "compression=%s, %.1f KB)",
            path, f"{self.n_snps:,}", f"{self.n_edges:,}",
            f"{_indices_c.size:,}", self.threshold,
            symmetric_only, quantize_r2, used_compression,
            path.stat().st_size / 1024,
        )

    @classmethod
    def load(cls, path: str | Path) -> "LDGraph":
        """Load a graph previously written by :meth:`save`.

        Auto-detects the container format:

        * ``.npz`` — numpy's own zipped format (default).
        * ``.ldz`` — custom zstd-compressed binary written by
          ``save(compact=True)`` when ``zstandard`` is available.

        Also auto-detects and reverses the on-disk optimisations
        (symmetric-only expansion, uint16 → float32 r² dequantisation)
        so downstream query methods always see the dense in-memory
        representation.
        """
        path = Path(path)

        # ---- Container ------------------------------------------------
        if path.suffix.lower() == ".ldz":
            arrays = _read_ldz(path)
        else:
            with np.load(path, allow_pickle=False) as z:
                arrays = {k: z[k] for k in z.files}

        # ---- Metadata (with backward-compat defaults) -----------------
        # Scalars are unwrapped with reshape(-1)[0]: ``.ldz`` files
        # written before the shape fix in ``_write_ldz`` store them as
        # 1-element arrays, which NumPy >= 2 refuses to pass to float().
        def _scalar(a):
            return np.asarray(a).reshape(-1)[0]
        threshold = float(_scalar(arrays["threshold"]))
        has_r2 = bool(_scalar(arrays["has_r2"])) if "has_r2" in arrays else True
        symmetric_only = (bool(_scalar(arrays["symmetric_only"]))
                          if "symmetric_only" in arrays else False)
        quantized_r2 = (bool(_scalar(arrays["quantized_r2"]))
                        if "quantized_r2" in arrays else False)
        # ``build`` comes back as either a 0-d unicode array or a
        # 1-elem unicode array depending on the container; unwrap
        # both cleanly.
        if "build" in arrays:
            _ba = arrays["build"]
            if hasattr(_ba, "item"):
                _build = str(_ba.item()) if _ba.shape == () else str(_ba[0])
            else:
                _build = str(_ba)
        else:
            _build = ""
        build = _build.strip() or None

        snp_names = arrays["snp_names"]
        indptr = arrays["indptr"]
        indices = arrays["indices"]
        data = arrays["data"]

        # ---- 1. Dequantise r² ----------------------------------------
        if quantized_r2:
            data = data.astype(np.float32) / 65535.0

        # ---- 2. Expand symmetric-only back to dense -------------------
        # Iterate each canonical edge src->dst (src < dst) and emit
        # both directions so the in-memory graph is symmetric and
        # query methods stay O(log k).
        if symmetric_only and indices.size > 0:
            _n = snp_names.shape[0]
            _sources = np.repeat(np.arange(_n, dtype=np.int32),
                                 np.diff(indptr))
            src_full = np.concatenate([_sources, indices])
            dst_full = np.concatenate([indices, _sources])
            val_full = np.concatenate([data, data])
            _order = np.lexsort((dst_full, src_full))
            src_full = src_full[_order]
            dst_full = dst_full[_order]
            val_full = val_full[_order]
            _counts = np.bincount(src_full, minlength=_n)
            indptr = np.concatenate([[0], np.cumsum(_counts)]).astype(np.int32)
            indices = dst_full
            data = val_full

        # Position key: present in v3+ files; older files leave it
        # unset and the graph fills it lazily on first access.
        return cls(
            snp_names=snp_names,
            indptr=indptr,
            indices=indices,
            data=data,
            threshold=threshold,
            has_r2=has_r2,
            build=build,
            snp_chrom=arrays.get("snp_chrom"),
            snp_pos=arrays.get("snp_pos"),
        )

    # ------------------------------------------------------------------
    # Builder
    # ------------------------------------------------------------------

    @classmethod
    def from_plink_ld(cls, path: str | Path,
                      r2_threshold: float = 0.1,
                      build: Optional[str] = None) -> "LDGraph":
        """Build a graph from a PLINK ``--r2`` (``.ld``) file.

        Two source layouts are supported:

        * **Standard** — the seven-column PLINK output
          (``CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2``).  Edges with
          ``R2 < r2_threshold`` are dropped and the exact r² is
          stored in ``data``.
        * **Reduced** — a trimmed presence-only layout with
          ``SNP_A`` + ``SNP_B`` (and optionally ``CHR_A`` / ``BP_A``
          for readability) but no ``R2`` column.  Every row is
          treated as an edge at r² == 1.0 sentinel — meaning
          "we know r² ≥ the PLINK-side filter that generated this
          file".  The graph's ``has_r2`` flag is set ``False`` so
          downstream code can gate r²-tier features (regional
          plots) that need real numbers.  Roughly 2× more compact
          than the standard format at the source.

        Pairs are always stored symmetrically (A→B and B→A) so
        ``any_linked`` / ``neighbors`` work regardless of insertion
        order.

        The per-node position key (``snp_chrom``/``snp_pos``) is taken
        from ``CHR_A``/``BP_A`` and ``CHR_B``/``BP_B`` when present, so
        rsID-named panels still carry positions.  Nodes without a
        column value fall back to parsing their SNP name.
        """
        path = Path(path)
        logger.info("Building LDGraph from PLINK .ld: %s "
                    "(r² threshold %.2f)", path, r2_threshold)
        df = pd.read_csv(path, sep=r"\s+", engine="python")
        # SNP_A + SNP_B are always required.  R2 is optional
        # (reduced-format files omit it — the row's presence is the
        # LD statement, honouring the PLINK-side r² filter that
        # produced the file).
        for col in ("SNP_A", "SNP_B"):
            if col not in df.columns:
                raise ValueError(
                    f"PLINK .ld file at {path} missing required "
                    f"column {col!r}.  Header seen: {list(df.columns)}"
                )
        has_r2 = "R2" in df.columns
        if has_r2:
            # Filter early — everything below the graph threshold is dropped.
            df = df[df["R2"] >= float(r2_threshold)].copy()
        else:
            logger.info(
                "  reduced format detected (no R2 column); "
                "every row treated as an edge at r² == 1.0 sentinel"
            )
        if df.empty:
            logger.warning(
                "LDGraph.from_plink_ld: no pairs survived r² >= %s; "
                "returning empty graph.", r2_threshold,
            )
            return cls(
                snp_names=np.array([], dtype="U1"),
                indptr=np.array([0], dtype=np.int32),
                indices=np.array([], dtype=np.int32),
                data=np.array([], dtype=np.float32),
                threshold=r2_threshold,
                has_r2=has_r2,
                build=build,
                snp_chrom=np.array([], dtype="U1"),
                snp_pos=np.array([], dtype=np.int32),
            )

        # Assign integer indices to every SNP seen.
        snp_a_arr = df["SNP_A"].astype(str).to_numpy()
        snp_b_arr = df["SNP_B"].astype(str).to_numpy()
        if has_r2:
            r2_arr = df["R2"].astype(np.float32).to_numpy()
        else:
            # Reduced-format sentinel: every edge equals 1.0.  The
            # graph's ``has_r2`` flag records that these are
            # sentinels, not real measurements.
            r2_arr = np.ones(len(df.index), dtype=np.float32)

        unique = np.unique(np.concatenate([snp_a_arr, snp_b_arr]))
        # np.unique on string input returns an object-dtype array;
        # cast to fixed-width Unicode so the .npz save path can use
        # ``allow_pickle=False`` on load (safer + faster).
        unique = unique.astype("U")
        snp_index = {str(name): i for i, name in enumerate(unique)}

        # Emit each pair as two directed edges.
        n_pairs = len(df.index)
        n_edges = 2 * n_pairs
        src = np.empty(n_edges, dtype=np.int32)
        dst = np.empty(n_edges, dtype=np.int32)
        val = np.empty(n_edges, dtype=np.float32)
        idx_a = np.fromiter((snp_index[s] for s in snp_a_arr),
                            dtype=np.int32, count=n_pairs)
        idx_b = np.fromiter((snp_index[s] for s in snp_b_arr),
                            dtype=np.int32, count=n_pairs)
        src[:n_pairs] = idx_a; dst[:n_pairs] = idx_b; val[:n_pairs] = r2_arr
        src[n_pairs:] = idx_b; dst[n_pairs:] = idx_a; val[n_pairs:] = r2_arr

        # Sort primarily by src, secondarily by dst so ``indices``
        # is sorted within each row for ``searchsorted`` lookups.
        order = np.lexsort((dst, src))
        src = src[order]; dst = dst[order]; val = val[order]

        # Build CSR indptr via bincount cumulative sum.
        counts = np.bincount(src, minlength=len(unique))
        indptr = np.concatenate([[0], np.cumsum(counts)]).astype(np.int32)

        snp_chrom, snp_pos = cls._positions_from_names(unique)
        _parts = []
        for side in ("A", "B"):
            if {f"SNP_{side}", f"CHR_{side}", f"BP_{side}"} <= set(df.columns):
                _parts.append(pd.DataFrame({
                    "SNP": df[f"SNP_{side}"].astype(str).to_numpy(),
                    "CHR": df[f"CHR_{side}"].to_numpy(),
                    "BP": pd.to_numeric(df[f"BP_{side}"], errors="coerce").to_numpy(),
                }))
        if _parts:
            from pycmplot.variant_matcher import VariantMatcher
            _cols = (pd.concat(_parts, ignore_index=True)
                       .dropna(subset=["CHR", "BP"])
                       .drop_duplicates("SNP"))
            _node = np.fromiter((snp_index[x] for x in _cols["SNP"]),
                                dtype=np.int64, count=len(_cols.index))
            _cmap = {u: VariantMatcher._normalise_chr(u)
                     for u in pd.unique(_cols["CHR"].astype(str))}
            _cchr = _cols["CHR"].astype(str).map(_cmap).to_numpy(dtype="U")
            if _cchr.dtype.itemsize > snp_chrom.dtype.itemsize:
                snp_chrom = snp_chrom.astype(_cchr.dtype)
            snp_chrom[_node] = _cchr
            snp_pos[_node] = _cols["BP"].to_numpy(dtype=np.int64).astype(np.int32)

        graph = cls(
            snp_names=unique,
            indptr=indptr,
            indices=dst,
            data=val,
            threshold=r2_threshold,
            has_r2=has_r2,
            build=build,
            snp_chrom=snp_chrom,
            snp_pos=snp_pos,
        )
        logger.info(
            "Built %s SNPs, %s directed edges (%s pairs kept from "
            "%s input rows)",
            f"{graph.n_snps:,}", f"{graph.n_edges:,}",
            f"{n_pairs:,}", f"{len(df.index) + (r2_arr < r2_threshold).sum():,}",
        )
        return graph


# ---------------------------------------------------------------------------
# .ldz zstd container (used when save(compact=True) can import zstandard)
# ---------------------------------------------------------------------------
#
# On-disk layout::
#
#     magic       : 4 bytes b"LDZ1"
#     header_len  : uint32 (bytes)
#     header      : JSON blob { "arrays": [ { "name", "dtype", "shape",
#                                             "n_bytes" }, ... ] }
#     arrays      : concatenated zstd-compressed frames, in header order
#
# One-shot compression per array (not chunked) — the whole graph fits
# comfortably in RAM at load time so streaming buys nothing.
# ---------------------------------------------------------------------------

_LDZ_MAGIC = b"LDZ1"


def _write_ldz(path: "Path", arrays: dict, level: int = 19) -> None:
    import struct
    import zstandard  # required at call time by save(compression='zstd')
    # Default level 19 — for numeric LD arrays, zstd only beats zlib
    # at high levels.  Users on time-sensitive workflows can lower
    # via ``save(compression_level=9)`` etc.
    cctx = zstandard.ZstdCompressor(level=int(level))

    header_arrays: list[dict] = []
    compressed_frames: list[bytes] = []
    for name, arr in arrays.items():
        # Record the true shape: ascontiguousarray promotes 0-d scalars
        # to shape (1,), which used to be written into the header.
        _arr = np.asarray(arr)
        raw = np.ascontiguousarray(_arr).tobytes()
        blob = cctx.compress(raw)
        header_arrays.append({
            "name": name,
            "dtype": _arr.dtype.str,
            "shape": list(_arr.shape),
            "n_bytes_compressed": len(blob),
        })
        compressed_frames.append(blob)

    header = json.dumps({"arrays": header_arrays}).encode("utf-8")
    with open(path, "wb") as fh:
        fh.write(_LDZ_MAGIC)
        fh.write(struct.pack("<I", len(header)))
        fh.write(header)
        for blob in compressed_frames:
            fh.write(blob)


def _read_ldz(path: "Path") -> dict:
    import struct
    import zstandard
    dctx = zstandard.ZstdDecompressor()

    with open(path, "rb") as fh:
        magic = fh.read(4)
        if magic != _LDZ_MAGIC:
            raise ValueError(
                f"{path} is not an .ldz file (magic {magic!r} != "
                f"{_LDZ_MAGIC!r})"
            )
        (header_len,) = struct.unpack("<I", fh.read(4))
        header = json.loads(fh.read(header_len).decode("utf-8"))
        arrays: dict = {}
        for spec in header["arrays"]:
            blob = fh.read(spec["n_bytes_compressed"])
            raw = dctx.decompress(blob)
            arr = np.frombuffer(raw, dtype=np.dtype(spec["dtype"]))
            if spec["shape"]:
                arr = arr.reshape(spec["shape"])
            else:
                # 0-d scalar: reshape to () and take .item() as needed.
                arr = arr.reshape(())
            arrays[spec["name"]] = arr
    return arrays
