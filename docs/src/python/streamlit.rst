.. _streamlit:

Chemiscope in ``streamlit``
===========================

The chemiscope python module provides a ``streamlit`` component that allows to embed
interactive chemiscope viewers directly into Streamlit applications.


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
    frames = ase.io.read("structures.xyz", ":")

    # Create or load your chemiscope dataset
    dataset = chemiscope.create_input(frames)

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
       mode="default",               # Visualization mode
       key="my_viewer",              # Unique widget key
       width="stretch",              # Width: "stretch" or pixels
       height=600,                   # Height in pixels
       selected_index=0,             # Initially selected structure
       on_select=callback_function   # Selection callback
   )

Visualization modes
+++++++++++++++++++

The ``mode`` parameter controls which panels are displayed:

``default`` - shows both structure and property panels, ``map`` shows only the property
map, ``structure`` shows only structure panel.

.. code-block:: python

   chemiscope.streamlit.viewer(dataset, mode=mode)

Interactive structure selection
+++++++++++++++++++++++++++++++

You can respond to structure selection changes using the ``on_select`` callback:

.. code-block:: python

   import streamlit as st
   import chemiscope

   def handle_selection():
        index = st.state.get("viewer")
        st.write(f"Selected structure: {index}")

   chemiscope.streamlit.viewer(
       dataset,
       key="viewer",
       on_select=handle_selection
   )

The `on_select` callback allows to modify other parts of the app whenever the 
structure is changed in the viewer.  By setting the ``selected_index`` to a 
Streamlit session state variable, you can also programmatically change the selected
structure from other parts of the app.

.. code-block:: python

   import streamlit as st
   import chemiscope

    st.slider(
            "Select structure by index",
            min_value=0,
            max_value=len(frames) - 1,
            key="selected_index",
        )

    chemiscope.streamlit.viewer(
        dataset,
        key="viewer",
        selected_index=st.session_state.selected_index,
    )


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


API reference
^^^^^^^^^^^^^

.. autofunction:: chemiscope.streamlit.viewer
