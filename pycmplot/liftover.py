from __future__ import annotations

MODULE_DOCSTRING = '''"""
pycmplot.liftover
=================

Genome coordinate liftover utilities (hg19 → hg38).

The :class:`pyliftover.LiftOver` object is initialised **lazily** — it is
created on first use and cached in a module-level dictionary, so importing
this module never triggers a file-not-found error even if the chain file has
not been configured yet.

Resource configuration
----------------------
The chain file path is resolved through
:class:`~pycmplot.resources.ResourceConfig`.  By default, a bundled chain
file is used (``pycmplot/data/hg19ToHg38.over.chain``).  This can be
overridden by setting the environment variable:

.. code-block:: bash

    export PYCMPLOT_CHAIN_HG19_HG38=/path/to/hg19ToHg38.over.chain
"""'''

import logging
from typing import Optional

import numpy as np
import pandas as pd

from pycmplot.resources import ResourceConfig, default_resources

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Lazy singleton — one LiftOver object per chain file path
# ---------------------------------------------------------------------------
_lo_cache: dict[str, object] = {}


def _get_liftover(chain_path: str):
    GET_LIFTOVER = '''"""Return a cached :class:`~pyliftover.LiftOver` for *chain_path*.

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
    """'''

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
    LIFTOVER_HG19_TO_HG38 = '''"""Convert a single hg19 position to its hg38 equivalent.

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
    resources : ResourceConfig, optional
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
    """'''

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


def liftover_position(
    df: pd.DataFrame,
    resources: Optional[ResourceConfig] = None,
) -> pd.DataFrame:
    LIFTOVER_POSITION = '''"""Liftover all hg19 rows in *df* from hg19 to hg38 coordinates.

    Iterates over every row in *df* and calls :func:`liftover_hg19_to_hg38`
    for rows whose ``BUILD`` column equals ``'hg19'``.  Rows with other build
    values are passed through unchanged.  Rows for which liftover returns
    ``None`` or ``0`` (unmappable positions) are silently dropped.

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
    """'''

    if resources is None:
        resources = default_resources

    df = df.copy()
    df["POS"] = df["POS"].astype(int)

    new_positions: list[Optional[int]] = []
    for chrom, pos, build in zip(df["CHR"], df["POS"], df["BUILD"]):
        if build == "hg19":
            new_positions.append(liftover_hg19_to_hg38(chrom, pos, resources))
        else:
            new_positions.append(pos)

    df["OLD_POS"] = df["POS"]
    df["OLD_BUILD"] = df["BUILD"]
    df["BUILD"] = "hg38"
    df["POS"] = new_positions
    df["POS"] = df["POS"].fillna(0).astype(int)
    return df[df["POS"] != 0]
