# -*- coding: utf-8 -*-
from .input import create_input, write_input  # noqa
from .structures import all_atomic_environments, librascal_atomic_environments  # noqa

try:
    # only import the chemiscope.show function if we have ipywidgets installed.
    import ipywidgets  # noqa

    from .jupyter import show  # noqa
except ImportError:
    pass

__version__ = "0.4.1"
