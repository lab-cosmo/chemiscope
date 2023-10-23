# -*- coding: utf-8 -*-
from .input import create_input, write_input  # noqa
from .structures import (  # noqa
    all_atomic_environments,
    arrow_from_vector,
    ase_tensors_to_ellipsoids,
    ase_vectors_to_arrows,
    center_shape,
    composition_properties,
    ellipsoid_from_tensor,
    extract_properties,
    librascal_atomic_environments,
)
from .version import __version__  # noqa

try:
    # only import the chemiscope.show function if we have ipywidgets installed.
    import ipywidgets  # noqa

    from .jupyter import show  # noqa
except ImportError:
    pass
