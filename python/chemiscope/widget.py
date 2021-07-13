# -*- coding: utf-8 -*-
import json

from ipywidgets import DOMWidget, ValueWidget, register
from traitlets import Unicode

from .input import create_input


module_name = "chemiscope-widget"
# TODO: figure out where the module version is defined (where is it used)
module_version = "^0.1.0"


@register
class ChemiscopeWidget(DOMWidget, ValueWidget):
    _view_name = Unicode("ChemiscopeView").tag(sync=True)
    _view_module = Unicode(module_name).tag(sync=True)
    _view_module_version = Unicode(module_version).tag(sync=True)

    data = Unicode().tag(sync=True)

    def __init__(self, data):
        super().__init__()
        self.data = json.dumps(data)


def display(frames, properties, meta={"name": " "}, cutoff=None):
    """
    Displays an instance of the chemiscope visualizer using the frames and
    properties given as parameters. It uses the `create_input` function.

    :param list frames: list of atomic structures. For now, only `ase.Atoms`_
                        objects are supported
    :param dict properties: optional dictionary of additional properties, see below
    :param dict meta: optional metadata of the dataset, see below
    :param float cutoff: optional. If present, will be used to generate
                         atom-centered environments
    """
    dict_input = create_input(
        frames=frames, properties=properties, meta=meta, cutoff=cutoff
    )
    return ChemiscopeWidget(dict_input)
