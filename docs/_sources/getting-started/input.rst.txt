.. _input:

Creating inputs with Python
===========================

Chemiscope loads datasets from a single JSON file containing structures, properties, and
optional metadata. For large files, use gzip compression (e.g., ``output.json.gz``) to
ease sharing and loading.

This section explains how to create chemiscope JSON input files using Python. For
details on the structure and format of the JSON input, see the :ref:`JSON file format
<json-format>` page.

The most common way to create and save a chemiscope dataset from Python is
:py:func:`chemiscope.write_input`, which creates the dictionary with chemiscope data
(structures, properties, metadata, etc) and writes the JSON file. Internally, this
function builds the Chemiscope JSON with :py:func:`chemiscope.create_input` that can be
used to create the input file and to manipulate the data structure before saving.

For quick interactive visualization, :py:func:`chemiscope.show` displays a dataset
directly in a Jupyter notebook without creating any files. If you already have a saved
chemiscope JSON file, :py:func:`chemiscope.show_input` loads and displays it in the same
way. It is also possible to display the saved Chemiscope files in
https://chemiscope.org/.

See the :ref:`Python module documentation <python-module>` for complete API details.
