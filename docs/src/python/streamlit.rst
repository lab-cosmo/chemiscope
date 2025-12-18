.. _streamlit:

Chemiscope in ``streamlit``
===========================

The chemiscope python module provides a ``streamlit`` component that allows to
embed interactive chemiscope viewers directly into Streamlit applications.


Installation
^^^^^^^^^^^^

To use the chemiscope Streamlit component, install chemiscope with the streamlit
extra:

.. code-block:: bash

    pip install chemiscope[streamlit]


Basic Usage
^^^^^^^^^^^

The component is available through the :py:func:`chemiscope.streamlit.viewer`
function. Here's a minimal example:

.. code-block:: python

    import ase.io
    import chemiscope
    import streamlit as st

    # Read structures
    structures = ase.io.read("structures.xyz", ":")

    # Create or load your chemiscope dataset
    dataset = chemiscope.create_input(structures)

    # Display the viewer
    chemiscope.streamlit.viewer(dataset, mode="structure")

Run the application:

.. code-block:: bash

    streamlit run your_app.py

Viewer parameters
+++++++++++++++++

The :py:func:`~chemiscope.streamlit.viewer` function accepts several optional parameters
to customize the viewer appearance and behavior:

.. code-block:: python

    chemiscope.streamlit.viewer(
        dataset,                      # Chemiscope dataset
        settings=settings,            # Custom settings dict
        mode="default",               # Visualization mode: "default", "map", or "structure"
        no_info_panel=False,          # Hide the info panel
        key="my_viewer",              # Unique widget key
        width="stretch",              # Width: "stretch" or pixels
        height=600,                   # Height in pixels
        selected_index=sel_id,        # Initial selected structure index
        on_select=callback_function   # Callback for selection changes
        on_settings_change=callback_function   # Callback for settings changes
    )

Visualization modes
+++++++++++++++++++

The ``mode`` parameter controls which panels are displayed:

``default`` shows both structure and property panels, ``map`` shows only the
property map, ``structure`` shows only structure panel.

.. code-block:: python

    chemiscope.streamlit.viewer(dataset, mode=mode)

Interactive structure selection
+++++++++++++++++++++++++++++++

You can respond to structure selection changes using the ``on_select`` callback:

.. code-block:: python

    import streamlit as st
    import chemiscope

    def handle_selection(selection_id):
        print(f"Selected structure: {selection_id}")
        # You can set state variables to update automatically the text
        st.session_state.selected_id = selection_id

    if "selected_id" not in st.session_state:
        st.session_state["selected_id"] = None

    chemiscope.streamlit.viewer(
        dataset,
        key="viewer",
        on_select=handle_selection
    )

    st.text(f"Selected structure ID: {st.session_state.selected_id}")

The `on_select` callback allows to modify other parts of the app whenever the
structure is changed in the viewer.  By setting the ``selected_index`` to a
Streamlit session state variable, you can also programmatically change the
selected structure from other parts of the app.

.. code-block:: python

    import streamlit as st
    import chemiscope

    st.slider(
            "Select structure by index",
            min_value=0,
            max_value=len(structures) - 1,
            key="selected_id",
        )

    def handle_selection(selection_id):
        print(f"Selected structure: {selection_id}")
        # You can set state variables to update automatically the text
        st.session_state.selected_id = selection_id

    if "selected_id" not in st.session_state:
        st.session_state["selected_id"] = None

    chemiscope.streamlit.viewer(
        dataset,
        selected_index=st.session_state.selected_id,
        on_select=handle_selection
    )

    st.text(f"Selected structure ID: {st.session_state.selected_id}")

Given that the chemiscope viewer has its own mechanism to change the selected
structure, it is essential to set up a callback that updates the same state
variable that controls the selected index through the ``selected_index`` option.
If you want to use the chemiscope component as a "passive" structure viewer, you
can use ``mode="structure"`` and ``no_info_panel=True``.


Custom settings
+++++++++++++++

You can customize the viewer appearance and behavior using the ``settings``
parameter, that can contain any parameter that can be used for the chemiscope
options, as obtained for instance using :py:func:`chemiscope.quick_settings`:

.. code-block:: python

   settings = chemiscope.quick_settings(
       x="property_1",
       y="property_2",
       color="energy",
       structure_settings={
           "bonds": True,
           "unitCell": True,
       }
   )

   chemiscope.streamlit.viewer(
       dataset,
       settings=settings,
       mode="default"
   )

Settings can also be manipulated in a bi-directional manner using the
``on_settings_change`` to define a callback. You can see the app example for an
example of how this can be used to manipulate the visualization and map style
from the app. Simple example:

.. code-block:: python

    def handle_settings_change(settings):
        print("Visualization settings changed:", settings)
        # store settings in session state, e.g. to update other components
        st.session_state.viewer_settings = settings

    chemiscope.streamlit.viewer(
        dataset,
        settings=settings,
        on_settings_change=handle_settings_change
    )

API reference
^^^^^^^^^^^^^

.. autofunction:: chemiscope.streamlit.viewer
