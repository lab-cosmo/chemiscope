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
    (e.g. map ranges, structure representation, camera location). See the
    :ref:`json-settings` section for a list of available settings.

    By changing the ``pinned`` setting, you can control the number
    of viewers shown in the grid, as well as which structures are
    displayed in each viewer.

    .. code-block:: python

        # Enable space filling representation for the active viewer
        widget.settings = {"structure": [{"spaceFilling": True}]}

        # Set up 4 viewers with the specified structures
        widget.settings = {"pinned": [0, 1, 12, 5]}

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

.. py:attribute:: warning_timeout

    An integer controlling how warning messages are displayed (in ms). Set to a positive
    value for auto dismissing warnings after the specified time, ``0`` for warnings that
    must be manually closed or ``-1`` to disable warnings entirely.

    .. code-block:: python

        # Auto dismiss warnings after 5 seconds
        widget.warning_timeout = 5000


Saving the dataset and settings as a standalone file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The "save" method allows exporting the current state of the widget 
as a standalone JSON file that can be opened in the web app 
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
    you may be able to ``await`` directly in a cell to wait for the 
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

.. automethod:: chemiscope.jupyter.ChemiscopeWidgetBase.save_map_image
.. automethod:: chemiscope.jupyter.ChemiscopeWidgetBase.save_structure_image
.. automethod:: chemiscope.jupyter.ChemiscopeWidgetBase.get_map_image
.. automethod:: chemiscope.jupyter.ChemiscopeWidgetBase.get_structure_image

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

.. automethod:: chemiscope.jupyter.ChemiscopeWidgetBase.save_structure_sequence
.. automethod:: chemiscope.jupyter.ChemiscopeWidgetBase.get_structure_sequence

Headless widget
---------------

If you want to use chemiscope to generate screenshots or inspect datasets programmatically without running a Jupyter notebook, you can use the headless widget. This widget uses a headless browser (Playwright) to render the chemiscope interface.

.. autofunction:: chemiscope.headless

The headless widget provides the same ``settings`` and ``selected_ids`` traitlets as the main widget, and exposes the same methods for saving screenshots and sequences.

.. code-block:: python

    import chemiscope
    from chemiscope import headless

    # Load your data
    structures = ...
    properties = ...

    # Create a headless widget instance
    # This automatically downloads and configures the required browser
    widget = headless(structures=structures, properties=properties, mode="structure")

    # Modify settings programmatically
    widget.settings = {'structure': [{'spaceFilling': True}]}

    # Save a snapshot of the active structure
    widget.save_structure_image("snapshot.png")

    # Save a sequence of images
    indices = [0, 1, 2]
    paths = ["frame_0.png", "frame_1.png", "frame_2.png"]
    widget.save_structure_sequence(indices, paths)

    # Save the dataset to a JSON file
    widget.save("dataset.json")

    # Clean up resources
    widget.close()

.. note::
    The headless widget requires the `playwright` python package to be installed. You can install it with ``pip install playwright`` and then run ``playwright install`` to download the browser binaries.

Dataset exploration
-------------------

.. autofunction:: chemiscope.explore
.. autofunction:: chemiscope.get_featurizer
.. autofunction:: chemiscope.metatomic_featurizer
