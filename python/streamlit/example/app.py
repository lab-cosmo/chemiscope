import json
from io import StringIO
from typing import Any, Dict, Tuple

import ase.io
import streamlit as st

import chemiscope
from chemiscope.streamlit import viewer


DEFAULT_SETTINGS: Dict[str, Any] = {
    "structure_mode": True,
    "map_mode": True,
    "height": 500,
    "opacity": 100,
    "size": 30,
    "palette": "viridis",
    "show_bonds": True,
    "space_filling": False,
}

if "settings" not in st.session_state:
    st.session_state.settings = DEFAULT_SETTINGS.copy()


def build_settings(
    palette: str,
    opacity: int,
    factor: int,
    show_bonds: bool,
    space_filling: bool,
) -> Dict[str, Any]:
    structure_settings_list = [
        {
            "bonds": show_bonds,
            "spaceFilling": space_filling,
        }
    ]
    map_settings_dict = {
        "color": {"palette": palette, "opacity": opacity},
        "size": {"factor": factor},
    }

    final_settings = {}
    if st.session_state.settings["map_mode"]:
        final_settings["map"] = map_settings_dict
    if st.session_state.settings["structure_mode"]:
        final_settings["structure"] = structure_settings_list

    return final_settings


@st.cache_data(show_spinner="Loading and processing structures...")
def process_uploaded_file(file_bytes: bytes, file_name: str):
    text = file_bytes.decode("utf-8")
    file_obj = StringIO(text)

    frames = ase.io.read(file_obj, ":", format="extxyz")
    properties = chemiscope.extract_properties(frames)

    dataset = chemiscope.create_input(
        frames=frames,
        properties=properties,
        meta={"name": file_name},
    )

    return dataset, frames


def on_structure_select(selected_id: int):
    if selected_id is not None:
        st.session_state.selected_index = selected_id


def on_settings_change(new_settings: Dict[str, Any]):
    if new_settings:
        st.session_state.settings.update(new_settings)


def display_selected_structure():
    selected_index = st.session_state.get("selected_index")
    st.subheader("Selected Structure")
    if selected_index is not None:
        frames = st.session_state.get("uploaded_frames")
        st.text(f"Index: {selected_index}")
        st.code(str(frames[selected_index]))
    else:
        st.text("No structure selected")


def create_sidebar_widgets(uploaded: bool) -> Tuple[Dict[str, Any], str, int, str]:
    st.header("Viewer Settings on Load")
    settings_locked = uploaded
    s = st.session_state.settings

    st.subheader("Viewer Modes")
    structure_mode = st.checkbox(
        "Structure", value=s["structure_mode"], disabled=settings_locked
    )
    map_mode = st.checkbox("Map", value=s["map_mode"], disabled=settings_locked)

    modes = []
    if structure_mode:
        modes.append("structure")
    if map_mode:
        modes.append("map")
    if not modes:
        modes = ["default"]

    is_structure_only_mode = modes == ["structure"]
    is_map_only_mode = modes == ["map"]

    st.subheader("Display Settings")
    height = st.slider(
        "Height (px)", 100, 1200, s["height"], 50, disabled=settings_locked
    )

    st.subheader("Structure Settings")
    show_bonds = st.checkbox(
        "Show Bonds", value=s["show_bonds"], disabled=is_map_only_mode
    )
    space_filling = st.checkbox(
        "Space Filling", value=s["space_filling"], disabled=is_map_only_mode
    )

    st.subheader("Map Settings")
    opacity = st.slider(
        "Opacity", 0, 100, s["opacity"], 10, disabled=is_structure_only_mode
    )
    size = st.slider(
        "Point size", 1, 100, s["size"], 10, disabled=is_structure_only_mode
    )
    palette = st.selectbox(
        "Palette",
        ["viridis", "magma", "plasma", "inferno", "cividis"],
        index=["viridis", "magma", "plasma", "inferno", "cividis"].index(s["palette"]),
        disabled=is_structure_only_mode,
    )

    viewer_settings = build_settings(palette, opacity, size, show_bonds, space_filling)
    mode_display = modes[0] if len(modes) == 1 else "default"

    return viewer_settings, height, mode_display


st.set_page_config(page_title="Chemiscope + Streamlit", layout="wide")
st.title("Chemiscope inside Streamlit")

uploaded = st.file_uploader(
    "Upload extended XYZ file", type=["xyz", "extxyz"], key="static_uploader_key"
)

with st.sidebar:
    settings, height, mode_display = create_sidebar_widgets(uploaded is not None)

if uploaded:
    uploaded_bytes = uploaded.getvalue()
    file_name = uploaded.name

    dataset, frames = process_uploaded_file(uploaded_bytes, file_name)
    st.session_state["uploaded_frames"] = frames

    if "selected_index" not in st.session_state:
        st.session_state.selected_index = 0

    st.slider(
        "Select structure by index",
        min_value=0,
        max_value=len(frames) - 1,
        key="selected_index",
    )

    viewer(
        dataset,
        height=height,
        mode=mode_display,
        settings=settings,
        key="chemiscope_viewer",
        selected_index=st.session_state.selected_index,
        on_select=on_structure_select,
        on_settings_change=on_settings_change,
    )

    st.markdown("---")
    display_selected_structure()

    st.subheader("Export")
    download_dataset = dataset.copy()
    download_dataset["settings"] = settings
    json_data = json.dumps(download_dataset, indent=2)
    st.download_button(
        label="Download JSON dataset",
        data=json_data,
        file_name=f"{file_name.split('.')[0]}.chemiscope.json",
        mime="application/json",
        use_container_width=True,
    )
