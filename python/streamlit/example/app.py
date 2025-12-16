"""
Chemiscope Streamlit example app

This app demonstrates an integration of the chemiscope viewer inside Streamlit-built web
application. It provides a complete workflow for uploading and visualizing structures
from XYZ files, and demonstrates a setup for bidirectional synchronization between the
Python state (Streamlit widgets) and JavaScript state (chemiscope viewer).

Run the application with: `streamlit run app.py`.
"""

import json
from typing import Any, Dict

import streamlit as st
from ui import (
    create_sidebar_widgets,
    display_selected_structure,
    process_uploaded_file,
)

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


# Map specific widget keys to the main defaults. These keys must match the `key`
# arguments in slider and checkpoint components in ui.py
WIDGET_DEFAULTS = {
    "mode": DEFAULT_SETTINGS["mode"],
    "height": DEFAULT_SETTINGS["height"],
    "opacity": DEFAULT_SETTINGS["opacity"],
    "size": DEFAULT_SETTINGS["size"],
    "palette": DEFAULT_SETTINGS["palette"],
    "show_bonds": DEFAULT_SETTINGS["show_bonds"],
    "space_filling": DEFAULT_SETTINGS["space_filling"],
}


def on_structure_select(selected_id):
    """Callback when a user clicks a point on the map or selects a structure"""
    if selected_id is not None:
        st.session_state.selected_index = selected_id


def on_settings_change(new_settings):
    """
    Callback when settings are changed inside the chemiscope widget. We capture these
    changes to update the streamlit sidebar sliders to match.
    """
    if not new_settings:
        return

    # update internal settings store
    st.session_state.visualizer_settings = new_settings

    # sync streamlit state with chemiscope map settings
    if "map" in new_settings:
        color = new_settings["map"]["color"]
        st.session_state["palette"] = color["palette"]
        st.session_state["opacity"] = int(color["opacity"])
        st.session_state["size"] = int(new_settings["map"]["size"]["factor"])

    # sync streamlit state with chemiscope structure settings
    if "structure" in new_settings and isinstance(new_settings["structure"], dict):
        struct = new_settings["structure"]
        st.session_state["show_bonds"] = bool(struct["bonds"])
        st.session_state["space_filling"] = bool(struct["spaceFilling"])


# Initialize session state variables if they don't exist. Without this, widgets would
# reset to default every time the user interacts with the app (triggering a re-run)
for k, v in WIDGET_DEFAULTS.items():
    if k not in st.session_state:
        st.session_state[k] = v

if "settings" not in st.session_state:
    st.session_state.settings = DEFAULT_SETTINGS.copy()


st.set_page_config(page_title="Chemiscope + Streamlit", layout="wide")
st.title("Chemiscope inside Streamlit")

uploaded = st.file_uploader(
    "Upload XYZ file", type=["xyz", "extxyz"], key="static_uploader_key"
)

with st.sidebar:
    settings, height, mode_display = create_sidebar_widgets(uploaded is not None)

if uploaded:
    uploaded_bytes = uploaded.getvalue()
    file_name = uploaded.name
    dataset, frames = process_uploaded_file(uploaded_bytes, file_name)

    if "selected_index" not in st.session_state:
        st.session_state.selected_index = 0

    st.slider("Select structure by index", 0, len(frames) - 1, key="selected_index")

    viewer(
        dataset,
        height=height,
        mode=mode_display,
        settings=settings,
        key="chemiscope_viewer",
        selected_index=st.session_state.selected_index,  # bind python selection
        on_select=on_structure_select,  # sync js selection -> python
        on_settings_change=on_settings_change,  # sync js settings -> python
    )

    st.markdown("---")
    display_selected_structure(frames)

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
