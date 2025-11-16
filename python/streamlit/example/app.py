import streamlit as st
import ase.io
import chemiscope
from io import StringIO

from chemiscope.streamlit import viewer  # or chemiscope.streamlit_viewer if you re-exported

st.set_page_config(page_title="Chemiscope + Streamlit demo", layout="wide")
st.title("Chemiscope inside Streamlit (plain TypeScript component)")

uploaded = st.file_uploader("Upload an extended XYZ file", type=["xyz"])

if uploaded is not None:
    text = uploaded.getvalue().decode("utf-8")  # or the encoding your files use
    file_obj = StringIO(text)

    # read from the uploaded file-like object
    frames = ase.io.read(file_obj, ":", format="extxyz")

    # very simple example: extract all properties chemiscope knows about
    properties = chemiscope.extract_properties(frames)

    # default meta/settings; you can customize as usual
    dataset = chemiscope.create_input(
        frames=frames,
        properties=properties,
        meta={"name": "Streamlit demo"},
    )

    n_structures = len(frames)  # or dataset["structures"]["size"]

    # Slider controlling which structure is shown
    idx = st.slider("Structure index", min_value=0, max_value=n_structures - 1, value=0)

    viewer(dataset, selected_index=idx, height=700, key="chemiscope")

else:
    st.info("Upload an .xyz file to see the Chemiscope viewer.")

