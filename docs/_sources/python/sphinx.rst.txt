.. _sphinx:

Chemiscope in `sphinx`
======================

The chemiscope python module also provides a `sphinx` extension that can process
a custom RST directive into an interacting viewer, that can be embedded in the 
generated HTML documentation.

Configuration
^^^^^^^^^^^^^

Assuming a project is already configured to use `sphinx` to build its documentation,
in order to be able to use the chemiscope directive one has to first
add the `chemiscope.sphinx` extension to the `docs/conf.py` file:

.. code-block:: python

   extensions = ["chemiscope.sphinx"]

Chemiscope directive
^^^^^^^^^^^^^^^^^^^^

In order to include a chemiscope viewer into a RST documentation page one can use a 
``chemiscope`` directive. When compiled with ``sphinx-build``, the directive will
be translated into an HTML embedding of and interactive widget, that loads a JSON 
file and displays it. The directive requires a mandatory argument, that indicates 
the path (relative to the RST file) where the chemiscope JSON is to be found, 
and an optional ``mode`` parameter that specifies the type of visualization to be used
(``default``, ``structure``, or ``map`` - ``default`` being unimaginatively the default).

For instance, to show a viewer combining structure and property panels you can 
simply use the directive
 
.. code-block:: rst

    .. chemiscope:: ../datasets/showcase.json.gz

Once compiled, this will show as this widget
 
.. chemiscope:: ../datasets/showcase.json.gz
    
The ``structure`` mode will show the structures only

.. chemiscope:: ../datasets/showcase.json.gz 
    :mode: structure
 
and the ``map`` mode only the property map

.. chemiscope:: ../datasets/showcase.json.gz
    :mode: map 
    
Viewing a documentation locally
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Because of the need to dynamically load the chemiscope module and the 
desired file, it is not possible to load the widget by viewing the
the raw HTML pages. Instead, one needs to run a local HTTP server. 
For instance, using the port 8765,

.. code-block:: bash

    cd docs/build/html
    python3 -m http.server 8765

You can then view the documentation in a browser by loading the URL
``localhost:8765``. 

