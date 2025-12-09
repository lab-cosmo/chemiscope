from __future__ import annotations

import os
from typing import Any, Callable, Mapping, Optional


_component_func = None


def _get_build_path() -> str:
    here = os.path.dirname(__file__)
    src_build = os.path.normpath(os.path.join(here, "..", "build"))

    if os.path.exists(src_build):
        return src_build
    return os.path.join(here, "stcomponent")


def viewer(
    dataset: Mapping[str, Any],
    *,
    width: str | int = "stretch",
    height: int = 550,
    mode: str = "default",
    key: Optional[str] = None,
    selected_index: Optional[int] = None,
    settings: Optional[dict] = None,
    on_select: Optional[Callable[..., Any]] = None,
) -> Any:
    """
    Render a Chemiscope viewer inside a Streamlit app.

    Parameters
    ----------
    dataset
        Chemiscope dataset (same dict used in write_input/create_input).
    width
        Width in pixels or "stretch" for full width.
    height
        Height in pixels for the viewer (or None to auto-size).
    mode
        Visualization mode: "default", "structure", or "map".
    key
        Optional Streamlit widget key.
    selected_index
        Optional initial structure index to select.
    settings
        Optional settings dictionary to configure the viewer appearance and behavior.
    on_select
        Optional callback function called when the structure selection changes.
        Receives the new selected index or None if unselected.
    """
    global _component_func

    if _component_func is None:
        try:
            import streamlit.components.v1 as components
        except ImportError as exc:
            raise ImportError(
                "Streamlit is required to use chemiscope.streamlit.viewer. "
                "Install it with: pip install 'chemiscope[streamlit]'"
            ) from exc

        _component_func = components.declare_component(
            "chemiscope_viewer", path=_get_build_path()
        )

    return _component_func(
        key=key,
        width=width,
        height=height,
        dataset=dataset,
        selected_index=selected_index,
        mode=mode,
        settings=settings,
        on_change=on_select,
    )
