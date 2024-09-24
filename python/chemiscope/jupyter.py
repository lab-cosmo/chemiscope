# -*- coding: utf-8 -*-
import gzip
import json
import warnings

import ipywidgets
from traitlets import Bool, Dict, Unicode

from .input import create_input
from .version import __version__

# this needs to match the version/name defined in
# python/jupyter/src/labextension.ts
PACKAGE_NAME = "chemiscope"
PACKAGE_VERSION = __version__


class ChemiscopeWidgetBase(ipywidgets.DOMWidget, ipywidgets.ValueWidget):
    _view_module = Unicode(PACKAGE_NAME).tag(sync=True)
    _view_module_version = Unicode(PACKAGE_VERSION).tag(sync=True)
    value = Unicode().tag(sync=True)
    has_metadata = Bool().tag(sync=True)

    # synchronized settings from the JS side. You can assign to this field to
    # change what's being displayed by chemiscope, but you need to assign a full
    # dictionary (`widget.settings["map"]["x"]["property"] = "foo"` will not
    # work, but `widget.settings = updated_settings` will).
    settings = Dict().tag(sync=True)
    # switch to disable automatic update of settings
    _settings_sync = Bool().tag(sync=True)

    def __init__(self, data, has_metadata):
        super().__init__()
        self.value = json.dumps(data)
        self.has_metadata = has_metadata
        if "settings" in data:
            self.settings = data["settings"]
        self._settings_sync = True

    def save(self, path):
        """
        Save the dataset displayed by this widget as JSON to the given ``path``.
        If ``path`` ends with ``.gz``, the file is written as gzip compressed
        JSON string.

        :param str path: where to save the dataset.
        """
        if path.endswith(".gz"):
            file = gzip.open(path, "w", 9)
        else:
            file = open(path, "wb")

        # update the settings in the data to the latest value
        data = json.loads(self.value)
        data["settings"] = self.settings

        file.write(json.dumps(data).encode("utf8"))
        file.close()

    def __repr__(self, max_length=64):
        # string representation of the chemiscope widget, outputs that are too large
        class_name = self.__class__.__name__
        string_repr = ""

        for key, value in json.loads(
            self.value
        ).items():  # loops over meta, structures, settings, ...
            value_repr = repr(value)
            truncated_repr = (
                (value_repr[:max_length] + "...}")
                if len(value_repr) > max_length
                else value_repr
            )
            string_repr += f"{key}={truncated_repr}, "

        return f"<{class_name}({string_repr})>"

    def _repr_html_(self):
        if _is_running_in_notebook():
            return self.__repr__()
        else:
            return ""


@ipywidgets.register
class ChemiscopeWidget(ChemiscopeWidgetBase):
    _view_name = Unicode("ChemiscopeView").tag(sync=True)

    def __init__(self, data, has_metadata):
        super().__init__(data, has_metadata)


@ipywidgets.register
class StructureWidget(ChemiscopeWidgetBase):
    _view_name = Unicode("StructureView").tag(sync=True)

    def __init__(self, data, has_metadata):
        super().__init__(data, has_metadata)


@ipywidgets.register
class MapWidget(ChemiscopeWidgetBase):
    _view_name = Unicode("MapView").tag(sync=True)

    def __init__(self, data, has_metadata):
        super().__init__(data, has_metadata)


def show_input(path, mode="default"):
    """
    Loads and shows the chemiscope widget in ``path``.
    If ``path`` ends with ``.gz``, the file is loaded as a gzip compressed
    JSON string.

    :param str path: load the chemiscope widget from path.

    :param str mode: widget mode, either ``default``, ``structure`` or ``map``.

    .. code-block:: python

        import chemiscope

        widget = chemiscope.show_input("dataset.json")

    ..
    """

    if mode not in ["default", "structure", "map"]:
        raise ValueError(
            f"invalid mode '{mode}' in chemiscope.load, expected one of "
            "'default', 'structure' or 'map'"
        )

    if mode == "default":
        widget_class = ChemiscopeWidget

    elif mode == "structure":
        widget_class = StructureWidget

    elif mode == "map":
        widget_class = MapWidget

    if path.endswith(".gz"):
        with gzip.open(path, "rt") as f:
            dict_input = json.load(f)
    elif path.endswith(".json"):
        with open(path, "rt") as f:
            dict_input = json.load(f)
    else:
        raise ValueError(
            "invalid file format in chemiscope.load,\
                         expected .json or .json.gz"
        )

    try:
        meta = dict_input["meta"]

        if meta == {"name": " "}:
            has_metadata = False
        else:
            has_metadata = True

    except KeyError:
        raise ValueError("missing metadata in chemiscope.load")

    return widget_class(dict_input, has_metadata=has_metadata)


def show(
    frames=None,
    properties=None,
    meta=None,
    environments=None,
    shapes=None,
    settings=None,
    mode="default",
):
    """
    Show the dataset defined by the given ``frames`` and ``properties``
    (optionally ``meta``, ``environments`` and ``shapes`` as well) using an embedded
    chemiscope visualizer inside a Jupyter notebook. These parameters have the same
    meaning as in the :py:func:`chemiscope.create_input` function.

    The ``mode`` keyword also allows overriding the default two-panels
    visualization to show only a structure panel (``mode = "structure"``) or the
    map panel (``mode = "map"``). These modes also make it possible to view a
    dataset for which properties (or frames) are not available.

    When inside a jupyter notebook, the returned object will create a new
    chemiscope visualizer displaying the dataset. The object exposes a
    `settings` traitlet, that allows to modify the visualization options
    (possibly even linking the parameters to another widget). Printing the value
    of the `settings` property is also a good way to see a full list of the
    available options.

    The returned object also have a ``save`` function that can be used to save
    the dataset to a ``.json`` or ``.json.gz`` file to load it in the main website
    later. The visualization options will be those used in the active widget, so
    this is also a good way to tweak the appearence of the visualization before
    saving it.

    .. code-block:: python

        import chemiscope
        from sklearn.decomposition import PCA
        import ase.io

        pca = PCA(n_components=3)

        frames = ase.io.read(...)
        properties = {
            "PCA": pca.fit_transform(some_data),
        }

        widget = chemiscope.show(frames, properties)
        # display the dataset in a chemiscope visualizer inside the notebook
        widget
        # ...

        # NB: due to how traitlet work, you should always set the value of
        # the `settings` property. Only the properties that are explicitly
        # indicated will be modified.
        widget.settings = {"map": {"palette": "seismic"}}
        widget.settings["map"]["palette"] = "viridis"  # << does nothing!

        # Save the file for later use
        widget.save("dataset.json")

    .. _ase.Atoms: https://wiki.fysik.dtu.dk/ase/ase/atoms.html
    """
    if not (_is_running_in_notebook() or _is_running_in_sphinx_gallery()):
        warnings.warn(
            "chemiscope.show only works in a jupyter notebook or a sphinx build",
            stacklevel=2,
        )

    has_metadata = meta is not None
    if not has_metadata:
        meta = {"name": " "}

    if mode == "default":
        widget_class = ChemiscopeWidget

    elif mode == "structure":
        if properties is None:
            properties = {}
        if "index" not in properties:
            # also adds an index property to have something to show in the info panel
            properties["index"] = {
                "target": "structure",
                "values": list(range(len(frames))),
            }

        widget_class = StructureWidget

    elif mode == "map":
        if properties is None:
            properties = {}

        widget_class = MapWidget

    else:
        raise ValueError(
            f"invalid mode '{mode}' in chemiscope.show, expected one of "
            "'default', 'structure' or 'map'"
        )

    dict_input = create_input(
        frames=frames,
        properties=properties,
        meta=meta,
        environments=environments,
        shapes=shapes,
        settings=settings,
    )

    if mode != "structure":
        # if there is a map, we need two properties, otherwise there will be no map
        # and error is only visible in the console
        if len(dict_input["properties"]) < 2:
            raise ValueError("Need at least two properties to visualize a map widget")

    return widget_class(dict_input, has_metadata=has_metadata)


def _is_running_in_notebook():
    """
    Check whether the python interpreter is running for a jupyter notebook or
    not. Taken from https://stackoverflow.com/a/39662359/4692076.
    """

    # apparently get_ipython is lost when this gets called from a callback of
    # an ipython widget. See https://github.com/jupyter/jupyter/issues/299
    try:
        from IPython import get_ipython
    except ImportError:
        return False

    try:
        shell = get_ipython().__class__.__name__
        # ZMQInteractiveShell is the standard Jupyter Kernel
        # Interpreter is used by pyiodide
        if shell in ["ZMQInteractiveShell", "Interpreter"]:
            return True
        elif shell == "TerminalInteractiveShell":
            return False
        else:
            return False
    except NameError:
        return False


def _is_running_in_sphinx_gallery():
    """
    Returns true if a file is being executed by sphinx-gallery.
    """

    # This is very hacky: we are relying on the fact that `sphinx_gallery.gen_gallery`
    # should only be imported when actually generating the gallery, and is not imported
    # otherwise.
    try:
        import sphinx_gallery

        sphinx_build = hasattr(sphinx_gallery, "gen_gallery")
    except ImportError:
        sphinx_build = False

    return sphinx_build
