Installation
------------

The latest release of the ``chemiscope`` Python package can be installed from
the Python package index, using

.. code-block::

    pip install chemiscope

This installs the core package, which includes utilities for creating a Jupyter widget
and  generating JSON input files. For advanced features, such as automatic dataset
exploration with featurization, install the optional dependencies:

.. code-block::

    bash pip install chemiscope[explore]

To install a development version, first ensure you have Node.js and npm installed (for
compiling the TypeScript library). Then, clone the GitHub repository and install
locally:

.. code-block::

    git clone https://github.com/lab-cosmo/chemiscope.git
    cd chemiscope
    pip install .
