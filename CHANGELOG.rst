Version 0.2.1 (2026-04-16)
==========================

Added
-----

- Added option to supply genome builds for summary stats files ``--build``` if 'BUILD' column is not in the files.
      This is a revertion to earlier versions of the plotting script. 
      Also made ``--build`` and ``--build_column`` optional allowing plotting to still proceed without genome build information.
      However, caution must be taken when multiple summary stats files are provided with different coordinate systems.
      For example, if ``--annotate`` is set, hits table generation will default to `hg38` coordinate, potentially leading to 
      in accurate annotations for variants in different coordinates.


Changed
-------

- Changed ``--annotate`` choices:
      Expnaded choices from `snp` and `gene` to include other columns in hits table.
      Also allowed for other columns in user supplied annotation table (available in python API only).


Version 0.1.9 (2026-04-14)
==========================

Fixed
-----

- Fixed column name auto-detection:
   - Expanded candidates list by adding lower and upper case versions for existing condidates.

  Fixed ``build`` option for ``prep_pycmplot_input_info`` function.
   - Updated it from optional to required parameter to be consistent with command line version.


Version 0.1.8 (2026-04-14)
==========================

Added
-----

- Added options to specify color of highlighted positions ``highlight_color``
  and line running through highlighted positions ``highlight_line_color``.

  Added command-line short forms for ``--colors``.

  Added command-line long forms for ``-r_min``, ``-r_max``, ``-t_space``, and
  ``-pad`` (`#1 <https://github.com/esohkevin/pycmplot/issues/1>`_)


Fixed
-----

- Fixed bug with __future__ import.

  Fixed command-line short form for ``--highlight_line``. (`#2 <https://github.com/esohkevin/pycmplot/issues/2>`_)
