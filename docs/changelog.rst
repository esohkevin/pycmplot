.. include:: ../CHANGELOG.rst

Changelog
=========

All notable changes to pycmplot are documented here.
This project follows `Semantic Versioning <https://semver.org/>`_.


Version 0.1.0 (2026)
---------------------

Initial public release.

**Added**

- Multi-track stacked linear Manhattan plots (``--mode lm``).
- Circos-style circular Manhattan plots (``--mode cm``).
- Automatic column detection for chromosome, position, SNP ID, p-value,
  and genome-build columns.
- hg19 → hg38 coordinate liftover via ``pyliftover``.
- Nearest-gene annotation and structured hits summary table output.
- Cluster-aware label spreading and intelligent arrow-angle calculation
  for dense genomic regions.
- ``--trim_pval`` option for memory-efficient processing of large files.
- ``--sort_track`` option to order tracks by label or chromosome length.
- Full Python API mirroring the command-line interface.
- Interactive Jupyter notebook tutorial.
