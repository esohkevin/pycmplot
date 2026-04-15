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
