.. _contributing:

Contributing
============

Bug reports, feature requests, and pull requests are welcome on the
`GitHub repository <https://github.com/esohkevin/pycmplot>`_.

Reporting issues
----------------

Please open a `GitHub issue <https://github.com/esohkevin/pycmplot/issues>`_
and include:

- Your pycmplot version (``pip show pycmplot``).
- Python version and OS.
- A minimal reproducible example, including the command or Python code
  used and the full error traceback.
- A small example input file if the issue is data-related.

Development setup
-----------------

.. code-block:: bash

   git clone https://github.com/esohkevin/pycmplot.git
   cd pycmplot
   python -m venv .venv
   source .venv/bin/activate
   pip install -e ".[dev]"

Building the documentation locally
------------------------------------

.. code-block:: bash

   cd docs
   pip install -r requirements.txt
   make html
   # open _build/html/index.html in your browser

Code style
----------

pycmplot follows `PEP 8 <https://peps.python.org/pep-0008/>`_. Please run
``flake8`` and ``black`` before submitting a pull request:

.. code-block:: bash

   black pycmplot/
   flake8 pycmplot/

Docstrings
----------

All public functions should use
`NumPy-style docstrings <https://numpydoc.readthedocs.io/en/latest/format.html>`_.
See the existing modules for examples.
