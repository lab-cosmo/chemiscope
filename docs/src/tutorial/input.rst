.. _input-format:

Creating chemiscope input files
===============================

When using the default chemiscope interface, all the structures and properties
in a dataset are loaded from a single JSON file. This sections describe how to
generate such JSON file, either using a pre-existing python script that does
most of the work for you, or by writing the JSON file directly. Since the
resulting JSON file can be quite large and thus harder to share with
collaborators, the default chemiscope interface also allows to load JSON files
compressed with gzip.

tl;dr if you would like to generate a simple chemiscope for your dataset, we
have a `Google Colab notebook <https://colab.research.google.com/drive/1NU0gjtaHcB5Oc3NbFZiQYtctY2190hDu>`_
that can help!

Tools able to create chemiscope input
-------------------------------------

``chemiscope`` Python module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The easiest way to create a JSON input file is to use the ``chemiscope`` Python
module. Install the package with ``pip install chemiscope``, and use
:py:func:`chemiscope.write_input` or :py:func:`chemiscope.create_input` in your
own script to generate the JSON file.

If all the properties you want to include into chemiscope are already stored in
a file `ase`_ can read, the ``chemiscope`` python package also install a
`chemiscope-input <chemiscope-input-cli_>`_ command line script.

Note that chemiscope does not compute structural representations or
dimensionality reduction, and you need to do this yourself or use another
package such as ASAP.

``ASAP``
^^^^^^^^

The `ASAP`_ structural analysis package is another tool that can directly
generate an output in chemiscope format.

``chemiscope`` functions reference
----------------------------------

.. autofunction:: chemiscope.write_input

.. autofunction:: chemiscope.create_input

.. autofunction:: chemiscope.extract_properties

.. autofunction:: chemiscope.extract_shapes

.. autofunction:: chemiscope.composition_properties

.. autofunction:: chemiscope.all_atomic_environments

.. autofunction:: chemiscope.librascal_atomic_environments

.. _ase: https://wiki.fysik.dtu.dk/ase/index.html
.. _ASAP: https://github.com/BingqingCheng/ASAP


.. _chemiscope-input-cli:

``chemiscope-input`` command line interface
-------------------------------------------

.. sphinx_argparse_cli::
   :module: chemiscope.main
   :func: _chemiscope_input_parser
   :prog: chemiscope-input
   :title:
   :group_title_prefix:
