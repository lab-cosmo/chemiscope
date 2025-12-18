from io import StringIO

import ase.io
import streamlit as st

import chemiscope


def display_selected_structure(structures):
    """Displays raw info about the currently selected structure"""
    selected_index = st.session_state.get("selected_index")
    st.subheader("Selected Structure")

    if selected_index is None:
        st.text("No structure selected")
        return

    st.text(f"Index: {selected_index}")
    st.code(str(structures[selected_index]))


def build_settings(
    palette: str, opacity: int, size: int, show_bonds: bool, space_filling: bool
):
    """
    Constructs the configuration dict required by chemiscope. The structure of this
    dictionary determines how the visualizer looks (colors, sizes) based on the current
    'mode'. See :py:func:`chemiscope.create_input` for details
    """
    settings = {}
    s_state = st.session_state.settings

    if s_state["mode"] in ["default", "map"]:
        settings["map"] = {
            "color": {"palette": palette, "opacity": opacity},
            "size": {"factor": size},
        }

    if s_state["mode"] in ["default", "structure"]:
        settings["structure"] = [{"bonds": show_bonds, "spaceFilling": space_filling}]

    return settings


def create_sidebar_widgets(uploaded: bool):
    """Generates the sidebar setting and returns the current settings"""
    st.header("Viewer Settings on Load")

    mode = st.radio(
        "Viewer Mode", ["default", "structure", "map"], key="mode", disabled=uploaded
    )

    st.subheader("Display Settings")
    height = st.slider("Height", 100, 1200, step=50, key="height")

    st.subheader("Map Settings")

    # disable map controls if we are in structure-only mode
    is_structure_only = mode == "structure"
    opacity = st.slider(
        "Opacity", 0, 100, step=10, key="opacity", disabled=is_structure_only
    )
    size = st.slider(
        "Point size", 1, 100, step=10, key="size", disabled=is_structure_only
    )

    palette_options = ["viridis", "magma", "plasma", "inferno", "cividis"]
    palette = st.selectbox(
        "Palette", palette_options, key="palette", disabled=is_structure_only
    )

    st.subheader("Structure Settings")

    # disable structure controls if we are in map-only mode
    is_map_only = mode == "map"
    show_bonds = st.checkbox("Show Bonds", key="show_bonds", disabled=is_map_only)
    space_filling = st.checkbox(
        "Space Filling", key="space_filling", disabled=is_map_only
    )

    viewer_settings = build_settings(palette, opacity, size, show_bonds, space_filling)
    return viewer_settings, height, mode


@st.cache_data(show_spinner="Loading and processing structures...")
def process_uploaded_file(file_bytes: bytes, file_name: str):
    """
    Parses the uploaded file bytes into ASE atoms and generates the chemiscope input.
    It uses @st.cache_data so we don't re-parse the potentially large XYZ file
    every time a slider is moved.
    """
    text = file_bytes.decode("utf-8")
    structures = ase.io.read(StringIO(text), ":", format="extxyz")
    properties = chemiscope.extract_properties(structures)

    dataset = chemiscope.create_input(
        structures=structures, properties=properties, meta={"name": file_name}
    )

    return dataset, structures
