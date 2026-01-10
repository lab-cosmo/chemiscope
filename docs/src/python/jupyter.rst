.. _jupyter:

Jupyter notebooks
=================

Chemiscope can be used as a widget in Jupyter notebooks, that should work in
both Jupyter classic and JupyterLab. The widget can be created in ``default``
mode (showing both a structure and a map panel), or used to display only
structures or only properties.

Once created, it is possible to interact with the widget using a traitlet
interface, modeled after `Jupyter widgets <http://ipywidgets.readthedocs.io>`_.

Creating a chemiscope widget
----------------------------

.. autofunction:: chemiscope.show
.. autofunction:: chemiscope.show_input

Widget properties and methods
-----------------------------

The widget object returned by :func:`chemiscope.show` provides traitlets to
interact with the visualization state programmatically. For more information
on using traitlets and widget events (like ``observe``), please refer to 
the `ipywidgets documentation
<https://ipywidgets.readthedocs.io/en/stable/examples/Widget%20Events.html>`_.

.. py:attribute:: settings

    A dictionary containing the current visualization settings. 
    This traitlet is synchronized between Python and the JavaScript 
    frontend. You can update it to change visualization options
    (e.g. map ranges, structure representation).

    By changing the ``pinned`` setting, you can control the number
    of viewers shown in the grid, as well as which structures are
    displayed in each viewer.

    .. code-block:: python

        # Enable space filling representation for the active viewer
        widget.settings = {"structure": [{"spaceFilling": True}]}

        # Set up 4 viewers with the specified structures
        widget.settings = {"pinned": [0, 1, 12, 5]

.. py:attribute:: selected_ids

    A dictionary describing the currently selected structure or 
    environment. It contains a ``structure`` index and optionally
    an ``atom`` index to indicate the active environment. 

    .. code-block:: python

        # Select the structure with index 3 and the environment 8
        widget.selected_ids = {"structure": 3, "atom": 8}

.. py:attribute:: active_viewer

    An integer indicating the index of the currently active viewer
    in the grid (0-based).


Saving the widget state as a standalone dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The "save" method allows exporting the current state of the widget 
as a standalone JSON viewer that can be opened in the web app 
(or loaded with :func:`chemiscope.show_input`).

.. automethod:: chemiscope.jupyter.ChemiscopeWidgetBase.save


Exporting images
~~~~~~~~~~~~~~~~

The chemiscope widget provides methods to capture snapshots of the map and 
structure panels programmatically. These methods are asynchronous, as they 
require communication with the browser's JavaScript engine to render and 
capture the data.

There are two types of methods: ``save_*`` methods, which write the image data 
directly to a file given the path, and ``get_*`` methods, which return the raw 
image bytes (PNG formatted).

.. note::
    These methods return ``asyncio.Future`` objects. Depending on the environment,
    environments, you may be able to ``await`` directly in a cell to wait for the 
    result. If ``await`` hangs in your environment, you can use 
    catch the ``get_*`` return value, and then in a separate cell get `.result()` 
    to access the actual return value.

Basic usage example:

.. code-block:: python

    # Cell 1: Create and display the widget
    import chemiscope
    widget = chemiscope.show(structures, properties)
    widget

.. code-block:: python

    # Cell 2: Save a snapshot of the current map to a file
    widget.save_map_image("current_map.png")

.. code-block:: python

    # Cell 3: Capture structure image data as raw PNG bytes
    img_future = widget.get_structure_image()
    
.. code-block:: python

    # Display image
    from IPython.display import Image
    img_data = img_future.result()
    display(Image(img_data))


Capturing sequences
~~~~~~~~~~~~~~~~~~~

You can also capture a sequence of structure snapshots (e.g., for a trajectory 
animation). The ``indices`` parameter can be a list of structure indices or 
a list of dictionaries specifying both structure and atom indices. 
An optional ``settings`` parameter allows applying specific visualization 
settings to each frame.

.. code-block:: python

    indices = [0, 10, 20]
    paths = ["frame_0.png", "frame_10.png", "frame_20.png"]
    
    # Capture and save a sequence of frames
    widget.save_structure_sequence(indices, paths)

.. automethod:: chemiscope.jupyter.ChemiscopeWidgetBase.save_map_image
.. automethod:: chemiscope.jupyter.ChemiscopeWidgetBase.save_structure_image
.. automethod:: chemiscope.jupyter.ChemiscopeWidgetBase.get_map_image
.. automethod:: chemiscope.jupyter.ChemiscopeWidgetBase.get_structure_image
.. automethod:: chemiscope.jupyter.ChemiscopeWidgetBase.save_structure_sequence
.. automethod:: chemiscope.jupyter.ChemiscopeWidgetBase.get_structure_sequence

Dataset exploration
-------------------

.. autofunction:: chemiscope.explore
.. autofunction:: chemiscope.get_featurizer
.. autofunction:: chemiscope.metatomic_featurizer
