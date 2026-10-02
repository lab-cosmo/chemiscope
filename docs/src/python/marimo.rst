.. _marimo:

Chemiscope in ``marimo``
========================

The chemiscope Python module can embed the notebook widget in a
`marimo <https://marimo.io>`_ notebook. Marimo does not display a raw anywidget,
so use :py:func:`chemiscope.marimo.viewer` instead of :py:func:`chemiscope.show`.


Installation
^^^^^^^^^^^^

.. code-block:: bash

    pip install chemiscope[marimo]


Basic usage
^^^^^^^^^^^

:py:func:`chemiscope.marimo.viewer` takes the same arguments as
:py:func:`chemiscope.show`. The last expression of the cell has to be the viewer
(or a layout that contains it):

.. code-block:: python

    import ase.io
    import chemiscope

    structures = ase.io.read("structures.xyz", ":")
    viewer = chemiscope.marimo.viewer(
        structures,
        properties=chemiscope.extract_properties(structures),
        settings=chemiscope.quick_settings(x="energy", y="volume"),
    )
    viewer

Run the notebook with:

.. code-block:: bash

    marimo edit app.py

An example that loads the showcase dataset lives at
``python/marimo/example/app.py``.


Reading the selection
^^^^^^^^^^^^^^^^^^^^^

The object ``viewer`` returns is a marimo UI element. Cells that read
``viewer.selected_ids``, ``viewer.settings``, or ``viewer.active_viewer`` re-run
when that value changes. ``viewer.save`` and the other widget methods are
available on the same object.

.. code-block:: python

    viewer.selected_ids

Assign a complete dictionary to ``viewer.settings``. Changing a nested key in
place does not update the widget.

.. autofunction:: chemiscope.marimo.viewer
.. autofunction:: chemiscope.marimo.viewer_input
