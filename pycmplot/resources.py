"""
pycmplot.resources
==================
Centralised configuration for external resource files that cannot be bundled
with the package (large reference files, chain files, etc.).

Users can supply paths in three ways, in order of priority:

1. Pass a :class:`ResourceConfig` instance directly to functions that need it.
2. Set environment variables before running:

   .. code-block:: bash

       export PYCMPLOT_CHAIN_HG19_HG38=/path/to/hg19ToHg38.over.chain
       export PYCMPLOT_GENEINFO_HG38=/path/to/Homo_sapiens.GRCh38.geneinfo.tsv.gz
       export PYCMPLOT_GENEINFO_HG19=/path/to/Homo_sapiens.GRCh37.geneinfo.tsv.gz
       export PYCMPLOT_FEATURESINFO=/path/to/Homo_sapiens.GRCh38.features.tsv.gz

3. Edit the defaults in this module for a site-wide installation.
"""

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
    """Paths to external reference files used by pycmplot.

    Attributes
    ----------
    chain_hg19_hg38 :
        LiftOver chain file for hg19 → hg38 conversion.
    geneinfo_hg38 :
        Tab-delimited gene info file for GRCh38 (used for nearest-gene annotation).
    geneinfo_hg19 :
        Tab-delimited gene info file for GRCh37 (fallback when data is hg19).
    featuresinfo :
        Extended features info file (all biotypes) for GRCh38.
    """

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
        """Return the path for *attr*, raising a clear error if it is unset."""
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
