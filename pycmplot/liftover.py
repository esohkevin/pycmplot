"""
pycmplot.liftover
=================
Genome coordinate liftover utilities (hg19 → hg38).

The :class:`pyliftover.LiftOver` object is initialised **lazily** — it is
created only when ``liftover_position`` is first called, so importing this
module never raises a :class:`FileNotFoundError` even if the chain file has
not been configured yet.
"""

from __future__ import annotations

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
    """Return a cached :class:`~pyliftover.LiftOver` for *chain_path*."""
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
    """Convert a single hg19 coordinate to hg38.

    Parameters
    ----------
    chrom:
        Chromosome name **without** the ``chr`` prefix (e.g. ``"1"``, ``"X"``).
    pos:
        0-based position (as expected by pyliftover).
    resources:
        :class:`~pycmplot.resources.ResourceConfig` instance.  Falls back to
        the module-level :data:`~pycmplot.resources.default_resources`.

    Returns
    -------
    int or None
        New hg38 position, or ``None`` if liftover failed for that coordinate.
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


def liftover_position(
    df: pd.DataFrame,
    resources: Optional[ResourceConfig] = None,
) -> pd.DataFrame:
    """Liftover all hg19 rows in *df* to hg38, in place.

    Expects columns ``CHR``, ``POS``, and ``BUILD``.  Rows whose ``BUILD``
    is ``'hg19'`` are lifted; others are left unchanged.  Rows that fail
    liftover (new position == 0 or ``None``) are dropped.

    Returns the modified DataFrame with two additional columns:
    ``OLD_POS`` and ``OLD_BUILD``.
    """
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
