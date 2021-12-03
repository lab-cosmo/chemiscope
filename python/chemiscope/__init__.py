# -*- coding: utf-8 -*-
from .input import create_input, write_input  # noqa

try:
    # only import the chemiscope.show function if we have ipywidgets installed.
    import ipywidgets  # noqa
    from .jupyter import show  # noqa
except ImportError:
    pass

__version__ = "0.3.4"
