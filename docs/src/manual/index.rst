.. _user-manual:

User manual
===========

This manual presents an overview of how to interact with the chemiscope viewer
in its "standard" implementation - including the web tool and the jupyter widget. 
It does not discuss in detail how to build a low-dimensional representation of 
a chemical dataset - chemiscope is solely a viewer, and there are many tools available
to perform this kind of analyses. 

This section starts introducing the concept of structural and physical properties, 
before describing how to use the different panels in the standard visualization. 
It continues by presenting how you can generate a chemiscope input file to load on
https://chemiscope.org, as well as within a standalone HTML viewer which does
not require internet connectivity. 
See the :ref:`Python module documentation <python-module>`  
for how to interact with chemiscope in a script, or to explore a dataset directly 
inside a jupyter notebook.

.. _chemiscope: https://chemiscope.org

.. toctree::
    :maxdepth: 2

    properties
    panels
    input
    sharing
