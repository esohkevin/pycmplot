"""Deprecation shim used by the 0.4.x public-API renames.

The pre-0.4.x names — ``prep_pycmplot_input_info``,
``get_sumstats_and_merged_sector_list``, ``plot_linear``,
``plot_circular``, and the three ``plot_qq_*`` variants — remain
importable for one release cycle so existing user scripts don't
break, but each call emits a single ``DeprecationWarning`` naming
the new short name and immediately delegates to it.

Usage
-----
::

    def load(sum_stats, labels, *, file_info=None, ...):
        ...  # the real function

    get_sumstats_and_merged_sector_list = _deprecated_alias(
        load,
        old_name="get_sumstats_and_merged_sector_list",
        new_name="load",
    )

The wrapper preserves ``__doc__``, ``__wrapped__`` (so Sphinx
autofunction picks up the real docstring), and ``__signature__``
(so IDEs still show the right parameters).  A module-level set of
already-warned names ensures the warning fires only once per
Python process, even when the deprecated name is looped over.
"""

from __future__ import annotations

import functools
import warnings
from typing import Callable, TypeVar

_WARNED: set[str] = set()

F = TypeVar("F", bound=Callable)


def _deprecated_alias(target: F, *, old_name: str, new_name: str) -> F:
    """Return a wrapper that warns once and forwards to *target*.

    Parameters
    ----------
    target : callable
        The renamed function.  Called with the wrapper's args/kwargs
        untouched.
    old_name : str
        The name being deprecated (used in the warning text).
    new_name : str
        The recommended replacement (used in the warning text).
    """

    @functools.wraps(target)
    def _wrapper(*args, **kwargs):
        if old_name not in _WARNED:
            _WARNED.add(old_name)
            warnings.warn(
                f"pycmplot: {old_name!r} is deprecated and will be removed "
                f"in a future release; use {new_name!r} instead. "
                "(This warning is emitted only once per Python process.)",
                DeprecationWarning,
                stacklevel=2,
            )
        return target(*args, **kwargs)

    # Cosmetic overrides so introspection points at the real function.
    _wrapper.__name__ = old_name
    _wrapper.__qualname__ = old_name
    _wrapper.__wrapped__ = target  # type: ignore[attr-defined]
    return _wrapper  # type: ignore[return-value]
