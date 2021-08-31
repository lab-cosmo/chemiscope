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


def show(frames, properties, meta={"name": " "}, cutoff=None):
    """
    Show the dataset defined by the given ``frames`` and ``properties`` using a
    embedded chemiscope visualizer inside a Jupyter notebook.

    :param list frames: list of atomic structures. For now, only `ase.Atoms`_
                        objects are supported
    :param dict properties: optional dictionary of additional properties, see below
    :param dict meta: optional metadata of the dataset, see below
    :param float cutoff: optional. If present, will be used to generate
                         atom-centered environments

    :returns: a :py:class:`ChemiscopeWidget` that will display itself
    """
    if not _is_running_in_notebook():
        raise Exception("chemiscope.show only works inside a jupyter notebook")

    dict_input = create_input(
        frames=frames, properties=properties, meta=meta, cutoff=cutoff
    )
    return ChemiscopeWidget(dict_input)


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
