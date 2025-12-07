import json
from io import StringIO

import ase.io
import chemiscope
import streamlit as st
from chemiscope.streamlit import viewer

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
        "Viewer height", min_value=100, max_value=2000, value=700, step=100
    )

    st.divider()

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

    _component = viewer(
        dataset,
        height=height,  # TODO
        mode=mode_display,
        #  settings=settings,  TODO
        key="chemiscope",
    )

    st.divider()
    st.subheader("Export")

    col1, col2 = st.columns(2)
    with col1:
        if st.button("Refresh Viewer", use_container_width=True):
            st.rerun()

    with col2:
        json_data = json.dumps(dataset, indent=2)
        st.download_button(
            label="Download JSON dataset",
            data=json_data,
            file_name=f"{uploaded.name.split('.')[0]}.chemiscope.json",
            mime="application/json",
            use_container_width=True,
        )
