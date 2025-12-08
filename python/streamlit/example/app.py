import json
from io import StringIO

import ase.io
import streamlit as st

import chemiscope
from chemiscope.streamlit import viewer


def build_settings(palette="inferno", opacity=100, factor=30, join_points=False):
    return {
        "map": {
            "color": {"palette": palette, "opacity": opacity},
            "size": {"factor": factor},
            "joinPoints": join_points,
        }
    }


st.set_page_config(page_title="Chemiscope + Streamlit demo", layout="wide")
st.title("Chemiscope inside Streamlit")

with st.sidebar:
    st.header("Visualization settings")
    st.subheader("Viewer Modes")

    structure_mode = st.checkbox("Structure", value=True)
    map_mode = st.checkbox("Map", value=True)

    available_modes = []
    if structure_mode:
        available_modes.append("structure")
    if map_mode:
        available_modes.append("map")
    if not available_modes:
        st.warning("At least one mode must be selected")
        available_modes = ["default"]

    height = st.slider(
        "Custom height (px)", min_value=100, max_value=2000, value=700, step=50
    )
    width_option = st.selectbox("Width mode", ["stretch", "custom"], index=0)
    if width_option == "custom":
        width_value = st.slider(
            "Custom width (px)", min_value=200, max_value=2000, value=1200, step=50
        )
        width = width_value
    else:
        width = width_option

    opacity = st.slider("Opacity", min_value=0, max_value=100, value=100, step=10)
    size = st.slider("Points size", min_value=1, max_value=100, value=30, step=10)
    palette = st.selectbox(
        "Palette",
        ["viridis", "magma", "plasma"],
    )
    join_points = st.checkbox("Link points", value=False)

    settings = build_settings(palette, opacity, size, join_points)

uploaded = st.file_uploader("Upload an extended XYZ file", type=["xyz"])

if uploaded is not None:
    with st.spinner("Loading and processing structures..."):
        try:
            text = uploaded.getvalue().decode("utf-8")
            file_obj = StringIO(text)

            progress_text = "Reading structures from file..."
            progress_bar = st.progress(0, text=progress_text)

            frames = ase.io.read(file_obj, ":", format="extxyz")

            if len(frames) == 0:
                st.error("No structures found in the uploaded file")
                st.stop()

            progress_bar.progress(33, text="Extracting properties...")

            properties = chemiscope.extract_properties(frames)
            progress_bar.progress(66, text="Creating Chemiscope dataset...")

            dataset = chemiscope.create_input(
                frames=frames,
                properties=properties,
                meta={"name": uploaded.name},
            )

            progress_bar.progress(100, text="Ready!")

        except Exception as e:
            st.error(f"Error processing file: {str(e)}")
            st.stop()

    mode_display = "default" if len(available_modes) > 1 else available_modes[0]

    viewer(
        dataset,
        height=height,
        width=width,
        mode=mode_display,
        settings=settings,
        key="chemiscope",
    )

    st.divider()
    st.subheader("Export")

    json_data = json.dumps(dataset, indent=2)
    st.download_button(
        label="Download JSON dataset",
        data=json_data,
        file_name=f"{uploaded.name.split('.')[0]}.chemiscope.json",
        mime="application/json",
        use_container_width=True,
    )
