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
    st.write("DEBUG: frames loaded =", len(frames))

    # very simple example: extract all properties chemiscope knows about
    properties = chemiscope.extract_properties(frames)

    st.write("DEBUG: properties")

    # default meta/settings; you can customize as usual
    dataset = chemiscope.create_input(
        frames=frames,
        properties=properties,
        meta={"name": "Streamlit demo"},
    )

    st.write("DEBUG: created dataset")
    viewer(dataset, height=700)
else:
    st.info("Upload an .xyz file to see the Chemiscope viewer.")

