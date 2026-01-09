"""
Chemiscope Explore Streamlit example app

This app demonstrates an integration of chemiscope.explore inside a Streamlit-built web
application. It provides a workflow for uploading atomic structures and automatically
generating low-dimensional visualizations using featurization and dimensionality
reduction.

Run the application with: `streamlit run explore_app.py`.

This requires additional dependencies. Install them with:
    pip install chemiscope[explore]
"""

import json
from io import StringIO

import streamlit as st

import ase.io

import chemiscope
from chemiscope.streamlit import viewer
from chemiscope.structures import extract_properties


FEATURIZER_NAME = "pet-mad-1.0"


def on_structure_select(selected_id):
    if selected_id is not None:
        st.session_state.selected_index = selected_id


@st.cache_data(show_spinner=False)
def load_structures(file_bytes: bytes, file_name: str):
    text = file_bytes.decode("utf-8")
    structures = ase.io.read(StringIO(text), ":", format="extxyz")
    return structures


@st.cache_resource(show_spinner=False)
def load_featurizer(featurizer_name: str, batch_size: int = 1):
    if featurizer_name == "pet-mad-1.0":
        try:
            import torch
            from pet_mad.explore import PETMADFeaturizer
        except ImportError as e:
            raise ImportError(
                "Required package not found. Please install the dependencies with "
                "`pip install chemiscope[explore]`."
            ) from e

        device = "cuda" if torch.cuda.is_available() else "cpu"
        return PETMADFeaturizer(version="1.0.0", device=device, batch_size=batch_size)
    else:
        raise ValueError(f"Unknown featurizer: {featurizer_name}")


@st.cache_data(show_spinner=False)
def create_explore_dataset(
    _structures,
    featurizer_name: str,
    file_name: str,
    batch_size: int = 1,
):
    properties = extract_properties(_structures, environments=None)

    featurizer = load_featurizer(featurizer_name, batch_size)
    features = featurizer(_structures, None)
    properties["features"] = features

    color_property = None
    for prop_name in properties:
        if prop_name != "features":
            color_property = prop_name
            break
    if color_property is None:
        color_property = "features[3]"

    settings = chemiscope.quick_settings(
        x="features[1]", y="features[2]", map_color=color_property
    )

    dataset = chemiscope.create_input(
        structures=_structures,
        properties=properties,
        settings=settings,
        metadata={"name": file_name},
    )

    return dataset


st.set_page_config(page_title="Chemiscope Explore", layout="wide")
st.title("🔬 Chemiscope Explore")

st.markdown(
    "Upload atomic structures to automatically generate a low-dimensional "
    "visualization using the **PET-MAD** featurizer."
)

with st.sidebar:
    st.header("Settings")

    batch_size = st.slider(
        "Batch size",
        min_value=1,
        max_value=128,
        value=32,
        key="batch_size",
    )

    height = st.slider("Viewer height", 400, 1000, value=600, step=50, key="height")

uploaded = st.file_uploader(
    "Upload XYZ file", type=["xyz", "extxyz"], key="explore_uploader"
)

if not uploaded:
    if "processed_file" in st.session_state:
        del st.session_state.processed_file
    if "computing" in st.session_state:
        del st.session_state.computing
    st.stop()

uploaded_bytes = uploaded.getvalue()
file_name = uploaded.name

cache_key = (file_name, batch_size)
is_new_file = st.session_state.get("processed_file") != cache_key

if is_new_file:
    st.session_state.computing = True

    with st.status("Processing...", expanded=True) as status:
        st.write("Loading structures...")
        structures = load_structures(uploaded_bytes, file_name)

        st.write("Loading featurizer model...")
        _ = load_featurizer(FEATURIZER_NAME, batch_size)

        st.write(f"Computing features for {len(structures)} structures...")
        dataset = create_explore_dataset(
            structures,
            FEATURIZER_NAME,
            file_name,
            batch_size=batch_size,
        )
        status.update(label="Done!", state="complete", expanded=False)

    if not uploaded:
        del st.session_state.computing
        st.rerun()

    st.session_state.processed_file = cache_key
    del st.session_state.computing
    st.rerun()
else:
    structures = load_structures(uploaded_bytes, file_name)
    dataset = create_explore_dataset(
        structures,
        FEATURIZER_NAME,
        file_name,
        batch_size=batch_size,
    )

st.caption(f"Loaded {len(structures)} structures from {file_name}")

try:
    if "selected_index" not in st.session_state:
        st.session_state.selected_index = 0

    viewer(
        dataset,
        height=height,
        mode="default",
        settings=dataset.get("settings"),
        key="chemiscope_explore_viewer",
        selected_index=st.session_state.selected_index,
        on_select=on_structure_select,
    )

    with st.sidebar:
        st.markdown("---")
        st.subheader("Export")
        st.download_button(
            "Download dataset",
            data=json.dumps(dataset, indent=2),
            file_name=f"{file_name.rsplit('.', 1)[0]}.chemiscope.json",
            mime="application/json",
            use_container_width=True,
        )

except ImportError as e:
    st.error(
        f"Missing dependencies: {e}\n\n"
        "Please install with: `pip install chemiscope[explore]`"
    )
except Exception as e:
    st.error(f"Error: {e}")
