from typing import Any, Optional


def viewer(
    structures=None,
    *,
    properties=None,
    metadata=None,
    environments=None,
    shapes=None,
    settings=None,
    parameters=None,
    mode: str = "default",
    warning_timeout: int = 10000,
    cache_structures: bool = True,
    # for backward compatibility with chemiscope.show
    frames=None,
    meta=None,
) -> Any:
    """
    Show a dataset inside a `marimo <https://marimo.io>`_ notebook.

    Arguments match :py:func:`chemiscope.show`. The returned object is a marimo UI
    element wrapping the chemiscope anywidget, so marimo re-runs cells that read
    ``selected_ids``, ``settings``, or ``active_viewer`` when the user changes the
    view. Methods such as ``save`` are proxied onto the underlying widget.

    .. code-block:: python

        import chemiscope

        viewer = chemiscope.marimo.viewer(
            structures,
            properties=properties,
            settings=chemiscope.quick_settings(x="energy", y="volume"),
        )
        viewer

    In another cell, ``viewer.selected_ids`` is the current selection.

    :param structures: list of atomic structures, as in :py:func:`chemiscope.show`
    :param dict properties: dictionary of properties
    :param dict metadata: optional metadata of the dataset
    :param list environments: optional list of ``(structure id, atom id, cutoff)``
    :param shapes: optional dictionary of shapes
    :param settings: optional dictionary of visualization settings
    :param dict parameters: optional dictionary of parameters for multidimensional
        properties
    :param str mode: ``"default"``, ``"structure"``, or ``"map"``
    :param int warning_timeout: timeout (in ms) for warning messages
    :param bool cache_structures: cache structure data on the Python side
    """
    try:
        import marimo as mo
    except ImportError as exc:
        raise ImportError(
            "marimo is required to use chemiscope.marimo.viewer. "
            "Install it with: pip install 'chemiscope[marimo]'"
        ) from exc

    from .widget import show

    widget = show(
        structures,
        properties=properties,
        metadata=metadata,
        environments=environments,
        shapes=shapes,
        settings=settings,
        parameters=parameters,
        mode=mode,
        warning_timeout=warning_timeout,
        cache_structures=cache_structures,
        frames=frames,
        meta=meta,
    )
    return mo.ui.anywidget(widget)


def viewer_input(
    path,
    *,
    settings: Optional[dict] = None,
    mode: str = "default",
    warning_timeout: int = 10000,
    cache_structures: bool = True,
) -> Any:
    """
    Load a chemiscope JSON file and show it inside a marimo notebook.

    Arguments match :py:func:`chemiscope.show_input`.

    :param path: path to a ``.json`` or ``.json.gz`` dataset, or a file-like object
    :param dict settings: optional settings merged over those stored in the file
    :param str mode: ``"default"``, ``"structure"``, or ``"map"``
    :param int warning_timeout: timeout (in ms) for warning messages
    :param bool cache_structures: cache structure data on the Python side
    """
    try:
        import marimo as mo
    except ImportError as exc:
        raise ImportError(
            "marimo is required to use chemiscope.marimo.viewer_input. "
            "Install it with: pip install 'chemiscope[marimo]'"
        ) from exc

    from .widget import show_input

    widget = show_input(
        path,
        settings=settings,
        mode=mode,
        warning_timeout=warning_timeout,
        cache_structures=cache_structures,
    )
    return mo.ui.anywidget(widget)
