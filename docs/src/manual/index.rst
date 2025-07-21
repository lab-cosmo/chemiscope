.. _user-manual:

User manual
===========

This manual presents an overview of how to use the Chemiscope viewer in its standard
implementations, including the web tool and the Jupyter widget. It does not discuss in
detail how to build a low-dimensional representation of a chemical dataset - Chemiscope
is primarily a visualization tool, but it also offers functionality to for dataset
exploration based on automatically computed structural representations.

This section starts introducing the concept of structural and physical properties,
before describing how to use the different panels in the standard visualization. It
continues by presenting how you can generate a Chemiscope input file to load on
https://chemiscope.org, as well as within a standalone HTML viewer that works offline.
See the :ref:`Python module documentation <python-module>` for how to interact with
Chemiscope in a script, or to explore a dataset directly inside a Jupyter notebook.

.. _chemiscope: https://chemiscope.org

.. toctree::
    :maxdepth: 2

    properties
    panels
    input
    sharing
