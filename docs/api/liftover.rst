.. _api_liftover:

pycmplot.liftover
=================

Lazy hg18 → hg38 and hg19 → hg38 coordinate conversion powered by
`pyliftover <https://github.com/konstantint/pyliftover>`_.
Conversion is triggered only when a genome-build column is detected (or
explicitly named via ``--build_column``), and only for rows annotated as
``hg18`` or ``hg19``. All other rows are passed through unchanged.

.. currentmodule:: pycmplot.liftover

.. automodule:: pycmplot.liftover
   :members:
   :undoc-members:
   :show-inheritance:
