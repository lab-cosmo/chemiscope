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

The widget object returned by :func:`chemiscope.show` and 
:func:`chemiscope.show_input` provides traitlets 
to interact with the visualization state programmatically.

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


The "save" method allows exporting the current state of the widget 
as a standalone JSON viewer that can be opened in the web app 
(or loaded with :func:`chemiscope.show_input`).

.. automethod:: chemiscope.jupyter.ChemiscopeWidgetBase.save

Dataset exploration
-------------------

.. autofunction:: chemiscope.explore
.. autofunction:: chemiscope.get_featurizer
.. autofunction:: chemiscope.metatomic_featurizer
