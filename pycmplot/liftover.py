"""
pycmplot.liftover
=================

Genome coordinate liftover utilities (hg18 → hg38 and hg19 → hg38).

The :class:`pyliftover.LiftOver` objects are initialised **lazily** — they
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
    """Return a cached :class:`~pyliftover.LiftOver` for *chain_path*.

    Loads the chain file on first call and stores the resulting
    :class:`~pyliftover.LiftOver` instance in a module-level dict.  Subsequent
    calls with the same *chain_path* return the cached object without re-reading
    the file.

    Parameters
    ----------
    chain_path : str
        Absolute path to a UCSC-format ``.over.chain`` (or ``.over.chain.gz``)
        file.

    Returns
    -------
    pyliftover.LiftOver
        A ready-to-use liftover object for the specified chain file.
    """

    if chain_path not in _lo_cache:
        from pyliftover import LiftOver  # deferred import

        logger.info("Loading LiftOver chain file: %s", chain_path)
        _lo_cache[chain_path] = LiftOver(chain_path)
    return _lo_cache[chain_path]


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------

def liftover_hg19_to_hg38(
    chrom: str,
    pos: int,
    resources: Optional[ResourceConfig] = None,
) -> Optional[int]:
    """Convert a single hg19 position to its hg38 equivalent.

    Uses a lazily loaded and cached :class:`~pyliftover.LiftOver` object backed
    by the chain file specified in *resources*.  When multiple hg38 mappings
    exist for a given position, the one with the highest chain score is returned.

    Parameters
    ----------
    chrom : str
        Chromosome name **without** the ``'chr'`` prefix (e.g. ``'1'``,
        ``'X'``).  The prefix is added internally before querying pyliftover.
    pos : int
        0-based hg19 position, as expected by :class:`pyliftover.LiftOver`.
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
    pyliftover uses **0-based** coordinates (BED convention).  GWAS summary
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

    results = lo.convert_coordinate(f"chr{chrom}", pos)
    if not results:
        return None
    # pyliftover returns sorted by chain score; take the best hit
    _new_chrom, new_pos, _strand, _score = results[0]
    return new_pos


def liftover_hg18_to_hg38(
    chrom: str,
    pos: int,
    resources: Optional[ResourceConfig] = None,
) -> Optional[int]:
    """Convert a single hg18 (NCBI36) position to its hg38 equivalent.

    Uses a lazily loaded and cached :class:`~pyliftover.LiftOver` object
    backed by the hg18→hg38 chain file specified in *resources*.  When
    multiple hg38 mappings exist for a given position, the one with the
    highest chain score is returned.

    Parameters
    ----------
    chrom : str
        Chromosome name **without** the ``'chr'`` prefix (e.g. ``'1'``,
        ``'X'``).  The prefix is added internally before querying
        pyliftover.
    pos : int
        0-based hg18 position, as expected by :class:`pyliftover.LiftOver`.
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

    results = lo.convert_coordinate(f"chr{chrom}", pos)
    if not results:
        return None
    _new_chrom, new_pos, _strand, _score = results[0]
    return new_pos


def liftover_position(
    df: pd.DataFrame,
    hg38_chr_limits: dict = None,
    resources: Optional[ResourceConfig] = None,
) -> pd.DataFrame:
    """Liftover all hg18/hg19 rows in *df* to hg38 coordinates.

    Iterates over every row in *df* and dispatches to
    :func:`liftover_hg19_to_hg38` for rows whose ``BUILD`` column equals
    ``'hg19'`` or to :func:`liftover_hg18_to_hg38` for rows whose ``BUILD``
    column equals ``'hg18'``.  Rows with any other build value are passed
    through unchanged.  Rows for which liftover returns ``None`` or ``0``
    (unmappable positions) are silently dropped.

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

    Returns
    -------
    pandas.DataFrame
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

    new_positions: list[Optional[int]] = []
    for chrom, pos, build in zip(df["CHR"], df["POS"], df["BUILD"]):
        if build == "hg19":
            new_positions.append(liftover_hg19_to_hg38(chrom, pos, resources))
        elif build == "hg18":
            new_positions.append(liftover_hg18_to_hg38(chrom, pos, resources))
        else:
            new_positions.append(pos)

    df["OLD_POS"] = df["POS"]
    df["OLD_BUILD"] = df["BUILD"]
    df["BUILD"] = "hg38"
    df["POS"] = new_positions
    df["POS"] = df["POS"].fillna(0).astype(int)

    clean_frames: list[pd.DataFrame] = []
    for chrom in df["CHR"].unique():
        chr_df = df[df["CHR"] == chrom]
        chr_limit = hg38_chr_limits.get(str(chrom))
        if chr_limit is not None:
            chr_df = chr_df[chr_df["POS"] <= chr_limit]
        else:
            logger.warning(
                "Chromosome %r not in hg38 chromosome-length table; "
                "keeping all variants without range check.", chrom,
            )
        clean_frames.append(chr_df)

    if not clean_frames:
        return df.iloc[0:0]

    clean_df = pd.concat(clean_frames, axis=0, ignore_index=True)
    return clean_df[clean_df["POS"] != 0]
