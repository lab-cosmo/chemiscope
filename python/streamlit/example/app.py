import json
from io import StringIO
from typing import Any, Dict

import ase.io
import streamlit as st

import chemiscope
from chemiscope.streamlit import viewer


DEFAULT_SETTINGS: Dict[str, Any] = {
    "mode": "default",
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
    palette: str, opacity: int, size: int, show_bonds: bool, space_filling: bool
):
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


@st.cache_data(show_spinner="Loading and processing structures...")
def process_uploaded_file(file_bytes: bytes, file_name: str):
    text = file_bytes.decode("utf-8")
    frames = ase.io.read(StringIO(text), ":", format="extxyz")
    properties = chemiscope.extract_properties(frames)

    dataset = chemiscope.create_input(
        frames=frames, properties=properties, meta={"name": file_name}
    )

    return dataset, frames


def on_structure_select(selected_id):
    if selected_id is not None:
        st.session_state.selected_index = selected_id


def on_settings_change(new_settings):
    if not new_settings:
        return

    s_state = st.session_state.settings

    if "map" in new_settings:
        color = new_settings["map"]["color"]
        s_state.update({"palette": color["palette"], "opacity": color["opacity"]})
        s_state["size"] = new_settings["map"]["size"]["factor"]

    if "structure" in new_settings:
        struct = new_settings["structure"][0]
        s_state.update(
            {"show_bonds": struct["bonds"], "space_filling": struct["spaceFilling"]}
        )


def display_selected_structure():
    selected_index = st.session_state.get("selected_index")
    st.subheader("Selected Structure")

    if selected_index is None:
        st.text("No structure selected")
        return

    frames = st.session_state.get("uploaded_frames")
    st.text(f"Index: {selected_index}")
    st.code(str(frames[selected_index]))


def create_sidebar_widgets(uploaded: bool):
    st.header("Viewer Settings on Load")
    s = st.session_state.settings

    mode = st.radio("Viewer Mode", ["default", "structure", "map"], disabled=uploaded)

    st.subheader("Display Settings")
    height = st.slider("Height", 100, 1200, s["height"], 50)

    st.subheader("Structure Settings")
    is_map_only = mode == "map"
    show_bonds = st.checkbox("Show Bonds", value=s["show_bonds"], disabled=is_map_only)
    space_filling = st.checkbox(
        "Space Filling", value=s["space_filling"], disabled=is_map_only
    )

    st.subheader("Map Settings")
    is_structure_only = mode == "structure"
    opacity = st.slider("Opacity", 0, 100, s["opacity"], 10, disabled=is_structure_only)
    size = st.slider("Point size", 1, 100, s["size"], 10, disabled=is_structure_only)
    palette_options = ["viridis", "magma", "plasma", "inferno", "cividis"]
    palette = st.selectbox(
        "Palette",
        palette_options,
        index=palette_options.index(s["palette"]),
        disabled=is_structure_only,
    )

    viewer_settings = build_settings(palette, opacity, size, show_bonds, space_filling)

    return viewer_settings, height, mode


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

    st.slider("Select structure by index", 0, len(frames) - 1, key="selected_index")

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
    st.download_button(
        "Download JSON dataset",
        data=json.dumps(download_dataset, indent=2),
        file_name=f"{file_name.split('.')[0]}.chemiscope.json",
        mime="application/json",
        use_container_width=True,
    )
