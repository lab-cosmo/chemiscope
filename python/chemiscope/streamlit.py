from __future__ import annotations

import os
import uuid
from typing import Any, Callable, Mapping, Optional


_component_func = None


def viewer(
    dataset: Mapping[str, Any],
    *,
    width: str | int = "stretch",
    height: int = 550,
    mode: str = "default",
    key: Optional[str] = None,
    selected_index: Optional[int] = None,
    settings: Optional[dict] = None,
    on_select: Optional[Callable[[Optional[int]], Any]] = None,
    on_settings_change: Optional[Callable[[dict], Any]] = None,
) -> Any:
    """
    Render a Chemiscope viewer inside a Streamlit app.

    Parameters
    ----------
    dataset: str
        Chemiscope dataset (same dict used in write_input/create_input).
    width: str | int, default is "stretch"
        Width in pixels or "stretch" for full width.
    height: int, default is 500
        Height in pixels for the viewer (or None to auto-size).
    mode: str
        Visualization mode: "default", "structure", or "map".
    key: is
        Optional Streamlit widget key.
    selected_index: int
        Optional initial structure index to select.
    settings: dict
        Optional settings dictionary to configure the viewer appearance and behavior.
    on_select: callable
        Optional callback function called when the structure selection changes.
        Receives the new selected index or None if unselected.
    on_settings_change: callable
        Called when settings change
    """
    global _component_func

    try:
        import streamlit as st
        import streamlit.components.v1 as components
    except ImportError:
        raise ImportError(
            "Streamlit is required to use chemiscope.streamlit.viewer. "
            "Install it with: pip install 'chemiscope[streamlit]'"
        )

    if _component_func is None:
        build_path = os.path.join(os.path.dirname(__file__), "stcomponent")
        _component_func = components.declare_component(
            "chemiscope_viewer", path=build_path
        )

    if key is None:
        key = f"chemiscope_viewer_{uuid.uuid4().hex[:8]}"

    state_key = f"{key}_state"
    if state_key not in st.session_state:
        st.session_state[state_key] = {
            "selected_index": selected_index,
            "settings": settings or {},
            "last_update": None,
        }

    def on_change():
        current_state = st.session_state[state_key]
        component_value = st.session_state[key]

        if component_value is None:
            return

        new_selection = component_value.get("selected_id")
        new_settings = component_value.get("settings", {})

        selection_changed = new_selection != current_state["selected_index"]
        settings_changed = new_settings != current_state["settings"]

        if selection_changed:
            current_state["selected_index"] = new_selection
            if on_select:
                on_select(new_selection)

        if settings_changed:
            current_state["settings"] = new_settings
            if on_settings_change:
                on_settings_change(new_settings)

        if selection_changed or settings_changed:
            current_state["last_update"] = "component"

    current_state = st.session_state[state_key]

    external_selection_changed = (
        selected_index is not None
        and selected_index != current_state["selected_index"]
        and current_state["last_update"] != "streamlit"
    )

    external_settings_changed = (
        settings is not None
        and settings != current_state["settings"]
        and current_state["last_update"] != "streamlit"
    )

    if external_selection_changed:
        current_state["selected_index"] = selected_index
        current_state["last_update"] = "streamlit"

    if external_settings_changed:
        current_state["settings"] = settings
        current_state["last_update"] = "streamlit"

    if current_state["last_update"] == "streamlit":
        current_state["last_update"] = None

    return _component_func(
        key=key,
        width=width,
        height=height,
        dataset=dataset,
        selected_index=current_state["selected_index"],
        mode=mode,
        settings=current_state["settings"],
        default=None,
        on_change=on_change,
    )
