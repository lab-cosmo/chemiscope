from __future__ import annotations

import os
from typing import Any, Mapping

import streamlit as st  # debug
import streamlit.components.v1 as components

# flip this to False when using the dev server (npm start)
_RELEASE = True


def _get_build_path() -> str:
    here = os.path.dirname(__file__)
    return os.path.join(here, "stcomponent")


if _RELEASE:
    _component_func = components.declare_component(
        "chemiscope_viewer",
        path=_get_build_path(),
    )
else:
    _component_func = components.declare_component(
        "chemiscope_viewer",
        url="http://localhost:3001",
    )


def viewer(
    dataset: Mapping[str, Any],
    *,
    key: str | None = None,
    height: int = 600,
    selected_index: int | None = None,
) -> Any:
    """
    Render a Chemiscope viewer inside a Streamlit app.

    Parameters
    ----------
    dataset
        Chemiscope dataset (same dict used in write_input/create_input).
    key
        Optional Streamlit widget key.
    height
        Height in pixels for the viewer.
    """

    return _component_func(dataset=dataset, height=height, selected_index=selected_index,
 key=key)

