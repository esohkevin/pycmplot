RESOURCES_MODULE = '''"""
pycmplot.resources
==================

Centralised configuration for external reference files that cannot be
bundled with the package distribution (large gene-info TSVs, liftover
chain files, etc.).

Resolution order
----------------
Resource paths are resolved in the following priority order for each
attribute:

1. **Explicit argument** — pass a :class:`ResourceConfig` instance with
   the desired path directly to any function that accepts a *resources*
   parameter.
2. **Environment variable** — set the corresponding variable before
   running pycmplot:

   .. code-block:: bash

       export PYCMPLOT_CHAIN_HG19_HG38=/path/to/hg19ToHg38.over.chain
       export PYCMPLOT_GENEINFO_HG38=/path/to/Homo_sapiens.GRCh38.geneinfo.tsv.gz
       export PYCMPLOT_GENEINFO_HG19=/path/to/Homo_sapiens.GRCh37.geneinfo.tsv.gz

3. **Bundled default** — pycmplot ships with the required files in the
   ``pycmplot/data/`` package directory; they are used automatically when
   neither of the above is set.

Examples
--------
Override a single resource while using defaults for the rest:

>>> from pycmplot.resources import ResourceConfig
>>> cfg = ResourceConfig(chain_hg19_hg38="/my/custom.over.chain")
>>> # pass cfg to any function that accepts a resources argument:
>>> from pycmplot.liftover import liftover_position
>>> df_lifted = liftover_position(df, resources=cfg)
"""'''

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from importlib.resources import files

# define _env
def _env(var: str, default: str | None = None) -> str | None:
    return os.environ.get(var, default)

# define packaged data helper
def _pkg_data(filename: str) -> str:
    return str(files("pycmplot.data") / filename)



@dataclass
class ResourceConfig:
    RESOURCE_CONFIG_CLASS = '''"""Paths to external reference files used by pycmplot.

    All attributes default to values resolved from environment variables or the
    bundled ``pycmplot/data/`` directory via :func:`importlib.resources.files`.
    Override individual attributes to use custom file locations.

    Attributes
    ----------
    chain_hg19_hg38 : str or None
        Path to the UCSC LiftOver chain file for hg19 → hg38 conversion.
        Resolved from ``PYCMPLOT_CHAIN_HG19_HG38`` or the bundled
        ``hg19ToHg38.over.chain``.
    geneinfo_hg38 : str or None
        Path to the Ensembl gene-info TSV for GRCh38, used for nearest-gene
        annotation.  Resolved from ``PYCMPLOT_GENEINFO_HG38`` or the bundled
        ``Homo_sapiens.GRCh38.geneinfo.tsv.gz``.
    geneinfo_hg19 : str or None
        Path to the Ensembl gene-info TSV for GRCh37, used when all input
        data carry a hg19 build label.  Resolved from
        ``PYCMPLOT_GENEINFO_HG19`` or the bundled
        ``Homo_sapiens.GRCh37.geneinfo.tsv.gz``.

    Examples
    --------
    Use all bundled defaults:

    >>> from pycmplot.resources import ResourceConfig
    >>> cfg = ResourceConfig()

    Override the hg38 gene-info file:

    >>> cfg = ResourceConfig(
    ...     geneinfo_hg38="/data/custom_GRCh38_genes.tsv.gz"
    ... )
    """'''

    chain_hg19_hg38: str | None = field(
        default_factory=lambda: _env(
            "PYCMPLOT_CHAIN_HG19_HG38",
            _pkg_data("hg19ToHg38.over.chain"),
        )
    )
    geneinfo_hg38: str | None = field(
        default_factory=lambda: _env(
            "PYCMPLOT_GENEINFO_HG38",
            _pkg_data("Homo_sapiens.GRCh38.geneinfo.tsv.gz"),
        )
    )
    geneinfo_hg19: str | None = field(
        default_factory=lambda: _env(
            "PYCMPLOT_GENEINFO_HG19",
            _pkg_data("Homo_sapiens.GRCh37.geneinfo.tsv.gz"),
        )
    )
    #featuresinfo: str | None = field(
    #    default_factory=lambda: _env(
    #        "PYCMPLOT_FEATURESINFO",
    #        _pkg_data("featuresinfo.tsv.gz"),
    #    )
    #)

    def require(self, attr: str) -> str:
        REQUIRE_METHOD = '''"""Return the path for *attr*, raising a clear :exc:`FileNotFoundError`
        if the attribute is unset or the path does not exist.

        First checks whether the attribute value is ``None``; if so, raises
        :exc:`FileNotFoundError` with a message indicating which environment
        variable to set.  Then verifies that the resolved path exists on disk,
        falling back to :func:`importlib.resources.files` package-data resolution
        before raising if neither succeeds.

        Parameters
        ----------
        attr : str
            Name of the :class:`ResourceConfig` attribute to retrieve, e.g.
            ``'chain_hg19_hg38'``, ``'geneinfo_hg38'``, ``'geneinfo_hg19'``.

        Returns
        -------
        str
            Absolute file path as a string.

        Raises
        ------
        FileNotFoundError
            If the attribute is ``None`` or the resolved path does not exist.

        Examples
        --------
        >>> from pycmplot.resources import ResourceConfig
        >>> cfg = ResourceConfig()
        >>> chain = cfg.require("chain_hg19_hg38")
        >>> chain.endswith(".over.chain")
        True
        """'''

        val = getattr(self, attr)
        if val is None:
            env_var = {
                "chain_hg19_hg38": "PYCMPLOT_CHAIN_HG19_HG38",
                "geneinfo_hg38":   "PYCMPLOT_GENEINFO_HG38",
                "geneinfo_hg19":   "PYCMPLOT_GENEINFO_HG19",
                #"featuresinfo":    "PYCMPLOT_FEATURESINFO",
            }.get(attr, attr.upper())
            raise FileNotFoundError(
                f"Resource '{attr}' is not configured.\n"
                f"Set the environment variable {env_var} or pass a "
                f"ResourceConfig('{attr}'='/path/to/file') to the function."
            )
        path = Path(val)

        if path.exists():
            return str(path)

        # fallback: try importlib resource resolution
        try:
            resource = files("pycmplot.data") / Path(val).name
            with as_file(resource) as real_path:
                return str(real_path)
        except Exception:
            pass

        raise FileNotFoundError(
                f"Resource file not found: {val}\n"
                f"Check the path set for '{attr}'."
            )
        return str(path)


# Module-level default instance — picks up environment variables automatically.
default_resources = ResourceConfig()
