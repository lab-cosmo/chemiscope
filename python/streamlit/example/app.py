import streamlit as st
import ase.io
import chemiscope
from io import StringIO

from chemiscope.streamlit import viewer  # or chemiscope.streamlit_viewer if you re-exported

def _mark_slider_change():
    st.session_state["selection_source"] = "slider"

st.set_page_config(page_title="Chemiscope + Streamlit demo", layout="wide")
st.title("Chemiscope inside Streamlit (plain TypeScript component)")

uploaded = st.file_uploader("Upload an extended XYZ file", type=["xyz"])

if uploaded is not None:
    # read the file content as text
    text = uploaded.getvalue().decode("utf-8")
    file_obj = StringIO(text)

    # read from the uploaded file-like object
    frames = ase.io.read(file_obj, ":", format="extxyz")

    # very simple example: extract all properties chemiscope knows about
    properties = chemiscope.extract_properties(frames)

    # default meta/settings
    dataset = chemiscope.create_input(
        frames=frames,
        properties=properties,
        meta={"name": "Streamlit demo"},
    )

    # interactive selection with a slider, and two-way communication with the viewer
    n_structures = len(frames)

    # Central state for selection
    if "selected_structure" not in st.session_state:
        st.session_state.selected_structure = 0

    # Track where the last change came from: "slider" or "viewer"
    if "selection_source" not in st.session_state:
        st.session_state.selection_source = "viewer"

    pending_key = "pending_selected_structure"
    if pending_key in st.session_state:
        new_idx = st.session_state[pending_key]
        # clamp + validate
        if isinstance(new_idx, int) and 0 <= new_idx < n_structures:
            st.session_state.selected_structure = new_idx
        del st.session_state[pending_key]

    # UI controls
    col1, col2 = st.columns([3, 1])

    with col1:
        st.slider(
            "Structure index",
            min_value=0,
            max_value=n_structures - 1,
            key="selected_structure",
            on_change=_mark_slider_change,
        )

    with col2:
        st.markdown(
            f"**Selected structure:** {st.session_state.selected_structure}"
        )        

    selection_from_viewer = viewer(dataset, selected_index=st.session_state.selected_structure, height=700, key="chemiscope")

    if selection_from_viewer is not None:
        try:
            new_idx = int(selection_from_viewer)
        except (TypeError, ValueError):
            new_idx = st.session_state.selected_structure

        if ( st.session_state.selection_source != "slider" and
            0 <= new_idx < n_structures and 
            new_idx != st.session_state.selected_structure):
            st.session_state[pending_key] = new_idx            
            st.rerun()

        st.session_state.selection_source = "viewer"


else:
    st.info("Upload an .xyz file to see the Chemiscope viewer.")

