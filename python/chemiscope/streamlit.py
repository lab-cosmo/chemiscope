from __future__ import annotations

import os
from typing import Any, Mapping


# flip this to False when using the dev server (npm start)
_RELEASE = True

# Lazily declared component function. We avoid importing streamlit at
# module import time so the package can be imported in environments where
# streamlit is not installed (e.g. building docs).
_component_func = None


def _get_build_path() -> str:
    here = os.path.dirname(__file__)
    # prefer a local development build next to the component source
    src_build = os.path.normpath(os.path.join(here, "..", "streamlit", "build"))
    if os.path.exists(src_build):
        return src_build
    # fallback to the packaged stcomponent folder
    return os.path.join(here, "stcomponent")


def _ensure_component():
    """Ensure `_component_func` is declared by importing Streamlit and
    calling `declare_component`. Raises a helpful ImportError if Streamlit
    is missing.
    """
    global _component_func
    if _component_func is not None:
        return

    try:
        import streamlit.components.v1 as components
    except Exception as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "Streamlit is required to use chemiscope.streamlit.viewer. "
            "Install it with: pip install 'chemiscope[streamlit]'"
        ) from exc

    if _RELEASE:
        _component_func = components.declare_component(
            "chemiscope_viewer", path=_get_build_path()
        )
    else:
        _component_func = components.declare_component(
            "chemiscope_viewer", url="http://localhost:3001"
        )


def viewer(
    dataset: Mapping[str, Any],
    *,
    key: str | None = None,
    width: str | int = "stretch",
    height: int | None = None,
    selected_index: int | None = None,
    mode: str = "default",
    settings: dict | None = None,
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
    selected_index
        Optional initial structure index to select.
    mode
        Optional visualization mode: "default", "structure", or "map".
    settings
        Optional settings dictionary (from quick_settings or custom dict) to configure
        the viewer appearance and behavior.
    """

    _ensure_component()
    return _component_func(
        key=key,
        width=width,
        height=height,
        dataset=dataset,
        selected_index=selected_index,
        mode=mode,
        settings=settings,
    )
