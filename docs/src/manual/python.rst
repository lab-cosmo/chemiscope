.. _python-module:

Python module
=============

The `chemiscope` package provides utilities to prepare JSON input files, 
and to visualize them in a :ref:`jupyter environment <jupyter>`.
The package also provides a small set of utility functions to convert `ase`_
structure data, and a command-line command to convert structure files to
a JSON input. 


Installation
------------

The latest release of the `chemiscope` Python package can be installed from 
the Python package index, using 

.. code-block::

    pip install chemiscope
    
To install a development version, one should instead clone the
`github repository <https://github.com/lab-cosmo/chemiscope>`_,
make sure that all the dependencies (including those needed to compile 
the typescript library) are present and then run 

.. code-block::

    pip install .
    
in the main folder. 



``chemiscope`` functions reference
----------------------------------

.. autofunction:: chemiscope.write_input

.. autofunction:: chemiscope.create_input

.. autofunction:: chemiscope.quick_settings

.. autofunction:: chemiscope.extract_properties

.. autofunction:: chemiscope.composition_properties

.. autofunction:: chemiscope.all_atomic_environments

.. autofunction:: chemiscope.librascal_atomic_environments

.. autofunction:: chemiscope.ellipsoid_from_tensor

.. autofunction:: chemiscope.arrow_from_vector

.. autofunction:: chemiscope.ase_vectors_to_arrows

.. autofunction:: chemiscope.ase_tensors_to_ellipsoids

.. _chemiscope-input-cli:

``chemiscope-input`` command line interface
-------------------------------------------

.. sphinx_argparse_cli::
   :module: chemiscope.main
   :func: _chemiscope_input_parser
   :prog: chemiscope-input
   :title:
   :group_title_prefix:
