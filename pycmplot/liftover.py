"""
pycmplot.liftover
=================

Genome coordinate liftover utilities (hg18 → hg38 and hg19 → hg38).

The :class:`liftover.ChainFile` objects are initialised **lazily** — they
are created on first use and cached in a module-level dictionary, so
importing this module never triggers a file-not-found error even if the
chain files have not been configured yet.

Supported conversions
---------------------
pycmplot harmonises input coordinates to GRCh38. Two source assemblies are
supported:

* ``hg19`` / GRCh37 → GRCh38 (default, bundled chain file)
* ``hg18`` / NCBI36 → GRCh38 (bundled chain file; used when input rows
  carry a ``hg18`` build label)

Resource configuration
----------------------
Chain file paths are resolved through
:class:`~pycmplot.resources.ResourceConfig`.  By default, bundled chain
files are used (``pycmplot/data/hg19ToHg38.over.chain.gz`` and
``pycmplot/data/hg18ToHg38.over.chain.gz``).  They can be overridden by
setting the environment variables:

.. code-block:: bash

    export PYCMPLOT_CHAIN_HG19_HG38=/path/to/hg19ToHg38.over.chain.gz
    export PYCMPLOT_CHAIN_HG18_HG38=/path/to/hg18ToHg38.over.chain.gz
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd

from pycmplot.resources import ResourceConfig, default_resources
from pycmplot.constants import hg38_chr_lengths

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Lazy singleton — one LiftOver object per chain file path
# ---------------------------------------------------------------------------
_lo_cache: dict[str, object] = {}


def _get_liftover(chain_path: str):
    """Return a cached :class:`liftover.ChainFile` for *chain_path*.

    Loads the chain file on first call and stores the resulting
    :class:`liftover.ChainFile` instance in a module-level dict.  Subsequent
    calls with the same *chain_path* return the cached object without re-reading
    the file.

    Parameters
    ----------
    chain_path : str
        Absolute path to a UCSC-format ``.over.chain`` (or ``.over.chain.gz``)
        file.

    Returns
    -------
    liftover.ChainFile
        A ready-to-use liftover object for the specified chain file.
    """

    if chain_path not in _lo_cache:
        from liftover import ChainFile  # deferred import

        logger.info("Loading LiftOver chain file: %s", chain_path)
        _lo_cache[chain_path] = ChainFile(chain_path)
    return _lo_cache[chain_path]


# Reasons recorded in the ``UNMAPPED_REASON`` column of the unmapped report.
REASON_NO_MAPPING = "no_chain_mapping"
REASON_OTHER_CHROM = "maps_to_other_chromosome"
REASON_BEYOND_LENGTH = "beyond_hg38_chromosome_length"


def _convert_detail(lo, chrom, pos):
    """Lift one position; return ``(new_pos, reason, lifted_chrom, lifted_pos)``.

    ``new_pos`` is ``None`` when the variant cannot be placed on the same
    chromosome in the target build, in which case ``reason`` says why and
    ``lifted_chrom`` / ``lifted_pos`` carry the off-chromosome hit (if any)
    so it can be reported.

    ``liftover.ChainFile.convert_coordinate`` returns a list of
    ``(chrom, pos, strand)`` tuples (empty when unmapped).  Hits that
    land on a *different* chromosome are treated as unmapped: callers
    only receive a position, so keeping one would silently pair the
    source ``CHR`` with a coordinate from another chromosome.
    """
    results = lo.convert_coordinate(f"chr{chrom}", pos)
    if not results:
        return None, REASON_NO_MAPPING, None, None
    new_chrom, new_pos, _strand = results[0]
    new_chrom = str(new_chrom).removeprefix("chr")
    if new_chrom != str(chrom).removeprefix("chr"):
        return None, REASON_OTHER_CHROM, new_chrom, int(new_pos)
    return new_pos, None, None, None


def _convert(lo, chrom, pos) -> Optional[int]:
    """Lift one position with *lo*; ``None`` if unmapped or off-chromosome."""
    return _convert_detail(lo, chrom, pos)[0]


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------

def liftover_hg19_to_hg38(
    chrom: str,
    pos: int,
    resources: Optional[ResourceConfig] = None,
) -> Optional[int]:
    """Convert a single hg19 position to its hg38 equivalent.

    Uses a lazily loaded and cached :class:`liftover.ChainFile` object backed
    by the chain file specified in *resources*.  When multiple hg38 mappings
    exist for a given position, the one with the highest chain score is returned.

    Parameters
    ----------
    chrom : str
        Chromosome name **without** the ``'chr'`` prefix (e.g. ``'1'``,
        ``'X'``).  The prefix is added internally before querying liftover.
    pos : int
        0-based hg19 position, as expected by :class:`liftover.ChainFile`.
    resources : ResourceConfig | Target Build Version, optional
        :class:`~pycmplot.resources.ResourceConfig` instance.  Falls back to
        :data:`~pycmplot.resources.default_resources` when ``None``.

    Returns
    -------
    int or None
        Corresponding 0-based hg38 position, or ``None`` if the position
        could not be mapped (unmapped region, chromosome gap, or deleted
        sequence).

    Notes
    -----
    liftover (like pyliftover) uses **0-based** coordinates (BED convention).  GWAS summary
    statistics files typically use **1-based** coordinates (VCF/Ensembl
    convention).  The caller (:func:`liftover_position`) is responsible for any
    coordinate-system adjustment.

    See Also
    --------
    liftover_position :
        Applies :func:`liftover_hg19_to_hg38` row-wise to a full DataFrame.

    Examples
    --------
    >>> from pycmplot.liftover import liftover_hg19_to_hg38
    >>> new_pos = liftover_hg19_to_hg38("11", 5246695)
    >>> new_pos
    5225465
    """

    if resources is None:
        resources = default_resources

    chain_path = resources.require("chain_hg19_hg38")
    lo = _get_liftover(chain_path)

    return _convert(lo, chrom, pos)


def liftover_hg18_to_hg38(
    chrom: str,
    pos: int,
    resources: Optional[ResourceConfig] = None,
) -> Optional[int]:
    """Convert a single hg18 (NCBI36) position to its hg38 equivalent.

    Uses a lazily loaded and cached :class:`liftover.ChainFile` object
    backed by the hg18→hg38 chain file specified in *resources*.  When
    multiple hg38 mappings exist for a given position, the one with the
    highest chain score is returned.

    Parameters
    ----------
    chrom : str
        Chromosome name **without** the ``'chr'`` prefix (e.g. ``'1'``,
        ``'X'``).  The prefix is added internally before querying
        liftover.
    pos : int
        0-based hg18 position, as expected by :class:`liftover.ChainFile`.
    resources : ResourceConfig, optional
        :class:`~pycmplot.resources.ResourceConfig` instance.  Falls back
        to :data:`~pycmplot.resources.default_resources` when ``None``.

    Returns
    -------
    int or None
        Corresponding 0-based hg38 position, or ``None`` if the position
        could not be mapped (unmapped region, chromosome gap, or deleted
        sequence).

    See Also
    --------
    liftover_hg19_to_hg38 :
        Equivalent helper for hg19 coordinates.
    liftover_position :
        Applies the appropriate per-row dispatcher to a full DataFrame.
    """

    if resources is None:
        resources = default_resources

    chain_path = resources.require("chain_hg18_hg38")
    lo = _get_liftover(chain_path)

    return _convert(lo, chrom, pos)


def liftover_position(
    df: pd.DataFrame,
    hg38_chr_limits: dict = None,
    resources: Optional[ResourceConfig] = None,
    return_unmapped: bool = False,
):
    """Liftover all hg18/hg19 rows in *df* to hg38 coordinates.

    Iterates over every row in *df* and dispatches to
    :func:`liftover_hg19_to_hg38` for rows whose ``BUILD`` column equals
    ``'hg19'`` or to :func:`liftover_hg18_to_hg38` for rows whose ``BUILD``
    column equals ``'hg18'``.  Rows with any other build value are passed
    through unchanged.  Rows that cannot be placed in hg38 are dropped:
    no chain mapping, a mapping onto a different chromosome, or a lifted
    position beyond the hg38 chromosome length.  Pass
    ``return_unmapped=True`` to get those rows back as a second table.

    Two provenance columns are added to the returned DataFrame so that the
    original coordinates remain accessible:

    * ``OLD_POS`` — the pre-liftover base-pair position.
    * ``OLD_BUILD`` — the original build value (``'hg19'``).

    After processing, the ``BUILD`` column is updated to ``'hg38'`` for all
    rows.

    Parameters
    ----------
    df : pandas.DataFrame
        Summary statistics DataFrame with canonical columns ``CHR``, ``POS``,
        and ``BUILD``.  The ``POS`` column is coerced to ``int`` before
        processing.
    resources : ResourceConfig, optional
        :class:`~pycmplot.resources.ResourceConfig` instance supplying the
        chain file path.  Falls back to
        :data:`~pycmplot.resources.default_resources` when ``None``.
    return_unmapped : bool, optional
        If ``True``, return ``(clean_df, unmapped_df)``.  ``unmapped_df``
        holds the dropped rows in their original coordinates (``POS`` /
        ``BUILD`` untouched) plus ``UNMAPPED_REASON`` (one of
        ``'no_chain_mapping'``, ``'maps_to_other_chromosome'``,
        ``'beyond_hg38_chromosome_length'``), and ``LIFTED_CHR`` /
        ``LIFTED_POS`` for the hit the chain file did return, if any.
        Default ``False`` returns only ``clean_df`` (backward compatible).

    Returns
    -------
    pandas.DataFrame or tuple of pandas.DataFrame
        A copy of *df* with:

        * ``POS`` replaced by hg38 coordinates for all hg19 rows.
        * ``BUILD`` set to ``'hg38'`` for all rows.
        * ``OLD_POS`` and ``OLD_BUILD`` columns added.
        * Rows with unmappable positions (new ``POS == 0``) removed.

    See Also
    --------
    liftover_hg19_to_hg38 :
        Single-position conversion function called internally.

    Examples
    --------
    >>> from pycmplot.liftover import liftover_position
    >>> df_hg38 = liftover_position(df)
    >>> df_hg38["BUILD"].unique()
    array(['hg38'], dtype=object)
    >>> "OLD_POS" in df_hg38.columns
    True
    """


    if resources is None:
        resources = default_resources

    if hg38_chr_limits is None:
        hg38_chr_limits = {k.replace("chr",""): v for k, v in hg38_chr_lengths.items()}
        

    df = df.copy()
    df["POS"] = df["POS"].astype(int)

    # Resolve each chain once, not once per row.
    _builds = set(df["BUILD"].unique())
    lifters = {}
    if "hg19" in _builds:
        lifters["hg19"] = _get_liftover(resources.require("chain_hg19_hg38"))
    if "hg18" in _builds:
        lifters["hg18"] = _get_liftover(resources.require("chain_hg18_hg38"))

    n = len(df.index)
    new_positions: list[Optional[int]] = [None] * n
    reasons: list[Optional[str]] = [None] * n
    lifted_chrom: list[Optional[str]] = [None] * n
    lifted_pos: list[Optional[int]] = [None] * n
    for i, (chrom, pos, build) in enumerate(zip(df["CHR"], df["POS"], df["BUILD"])):
        lo = lifters.get(build)
        if lo is None:
            new_positions[i] = pos
            continue
        new_positions[i], reasons[i], lifted_chrom[i], lifted_pos[i] = (
            _convert_detail(lo, chrom, pos)
        )

    reasons_s = pd.Series(reasons, index=df.index, dtype=object)
    new_pos_s = pd.Series(new_positions, index=df.index, dtype="float64")

    # Range check against hg38 chromosome lengths.
    chr_str = df["CHR"].astype(str)
    limits = chr_str.map(hg38_chr_limits)
    for chrom in chr_str[limits.isna()].unique():
        logger.warning(
            "Chromosome %r not in hg38 chromosome-length table; "
            "keeping all variants without range check.", chrom,
        )
    beyond = reasons_s.isna() & limits.notna() & (new_pos_s > limits)
    reasons_s[beyond] = REASON_BEYOND_LENGTH
    lifted_pos_s = pd.Series(lifted_pos, index=df.index, dtype="Int64")
    lifted_pos_s[beyond] = new_pos_s[beyond].astype("Int64")
    lifted_chrom_s = pd.Series(lifted_chrom, index=df.index, dtype=object)
    lifted_chrom_s[beyond] = chr_str[beyond]
    # A lifted position of 0 was historically treated as unmapped.
    zero = reasons_s.isna() & (new_pos_s == 0)
    reasons_s[zero] = REASON_NO_MAPPING

    dropped = reasons_s.notna()

    unmapped_df = df[dropped].copy()
    unmapped_df["UNMAPPED_REASON"] = reasons_s[dropped]
    unmapped_df["LIFTED_CHR"] = lifted_chrom_s[dropped]
    unmapped_df["LIFTED_POS"] = lifted_pos_s[dropped]

    df["OLD_POS"] = df["POS"]
    df["OLD_BUILD"] = df["BUILD"]
    df["BUILD"] = "hg38"
    df["POS"] = new_pos_s.fillna(0).astype(int)
    kept = df[~dropped]

    # Keep the historical row layout: grouped by chromosome in
    # first-appearance order, fresh RangeIndex.
    if kept.empty:
        clean_df = kept.reset_index(drop=True)
    else:
        clean_df = pd.concat(
            [kept[kept["CHR"] == c] for c in kept["CHR"].unique()],
            axis=0, ignore_index=True,
        )

    if dropped.any():
        counts = unmapped_df["UNMAPPED_REASON"].value_counts()
        logger.info(
            "Liftover: %s of %s variants could not be placed in hg38 (%s).",
            int(dropped.sum()), n,
            ", ".join(f"{k}: {v}" for k, v in counts.items()),
        )

    if return_unmapped:
        return clean_df, unmapped_df.reset_index(drop=True)
    return clean_df
