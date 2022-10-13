# -*- coding: utf-8 -*-
import gzip
import json

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
        self.settings = {}
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


def show(
    frames=None,
    properties=None,
    meta=None,
    environments=None,
    settings=None,
    mode="default",
):
    """
    Show the dataset defined by the given ``frames`` and ``properties``
    (optionally ``meta`` and ``environments`` as well) using a embedded chemiscope
    visualizer inside a Jupyter notebook. These parameters have the same meaning
    as in the :py:func:`chemiscope.create_input` function.

    The ``mode`` keyword also allows overriding the default two-panels
    visualization to show only a structure panel (``mode = "structure"``) or the
    map panel (``mode = "map"``). These modes also make it possible to view a
    dataset for which properties (or frames) are not available.

    When inside a jupyter notebook, the returned object will create a new
    chemiscope visualizer displaying the dataset. The returned object also have
    a ``save`` function that can be used to save the dataset to a ``.json`` or
    ``.json.gz`` file to load it in the main website later.

    .. code-block:: python

        import chemiscope
        from sklearn.decomposition import PCA
        import ase.io

        pca = PCA(n_components = 3)

        frames = ase.io.read(...)
        properties = {
            "PCA": pca.fit_transform(some_data)
        }

        widget = chemiscope.show(frames, properties)
        # display the dataset in a chemiscope visualizer inside the notebook
        widget
        # ...


        # Save the file for later use
        widget.save("dataset.json")

    .. _ase.Atoms: https://wiki.fysik.dtu.dk/ase/ase/atoms.html
    """
    if not _is_running_in_notebook():
        raise Exception("chemiscope.show only works inside a jupyter notebook")

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
        if shell == "ZMQInteractiveShell":
            return True
        elif shell == "TerminalInteractiveShell":
            return False
        else:
            return False
    except NameError:
        return False
