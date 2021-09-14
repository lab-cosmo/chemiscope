# -*- coding: utf-8 -*-
import json
import gzip
import ipywidgets
from traitlets import Unicode

from .input import create_input


@ipywidgets.register
class ChemiscopeWidget(ipywidgets.DOMWidget, ipywidgets.ValueWidget):
    _view_name = Unicode("ChemiscopeView").tag(sync=True)
    _view_module = Unicode("chemiscope-widget").tag(sync=True)

    data = Unicode().tag(sync=True)

    def __init__(self, data):
        super().__init__()
        self.data = json.dumps(data)

    def save(self, path):
        """
        Save the dataset displayed by this :py:class:`ChemiscopeWidget` as JSON
        to the given ``path``. If ``path`` ends with ``.gz``, the file is
        written as gzip compressed JSON string.

        :param str path: where to save the dataset.
        """
        if path.endswith(".gz"):
            file = gzip.open(path, "w", 9)
        else:
            file = open(path, "wb")

        file.write(self.data.encode("utf8"))
        file.close()


@ipywidgets.register
class StructureWidget(ipywidgets.DOMWidget, ipywidgets.ValueWidget):
    _view_name = Unicode("StructureView").tag(sync=True)
    _view_module = Unicode("chemiscope-widget").tag(sync=True)

    data = Unicode().tag(sync=True)

    def __init__(self, data):
        super().__init__()
        self.data = json.dumps(data)


def show(frames, properties=None, meta=None, cutoff=None, mode="default"):
    """
    Show the dataset defined by the given ``frames`` and ``properties``
    (optionally ``meta`` and ``cutoff`` as well) using a embedded chemiscope
    visualizer inside a Jupyter notebook. These parameters have the same meaning
    as in the :py:func:`chemiscope.create_input` function.
    The ``mode`` keyword also allows overriding the default two-panels visualization
    to show only a structure panel or the map panel. These modes also make it possible
    to view a dataset for which properties (or frames) are not available.

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

    if meta is None:
        meta = {"name": " "}

    if mode == "default":
        dict_input = create_input(
            frames=frames, properties=properties, meta=meta, cutoff=cutoff
        )
        return ChemiscopeWidget(dict_input)
    elif mode == "structure":
        if properties is None:
            properties = {}
        if "index" not in properties:
            # also adds an index property to have something to show in the info panel
            properties["index"] = list(range(len(frames)))

        dict_input = create_input(
            frames=frames, properties=properties, meta=meta, cutoff=cutoff
        )
        return StructureWidget(dict_input)
    else:
        raise ValueError(
            f"invalid mode '{mode}' in chemiscope.show, expected one of 'default' or 'structure'"
        )


def _is_running_in_notebook():
    """
    Check whether the python interpreter is running for a jupyter notebook or
    not. Taken from https://stackoverflow.com/a/39662359/4692076.
    """
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
