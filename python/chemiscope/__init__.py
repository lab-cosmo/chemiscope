# -*- coding: utf-8 -*-
from .input import create_input, quick_settings, write_input  # noqa: F401
from .structures import (  # noqa: F401
    all_atomic_environments,
    arrow_from_vector,
    ase_merge_pi_frames,
    ase_tensors_to_ellipsoids,
    ase_vectors_to_arrows,
    center_shape,
    composition_properties,
    ellipsoid_from_tensor,
    extract_properties,
    librascal_atomic_environments,
)
from .version import __version__  # noqa: F401

try:
    # only import the chemiscope.show function if we have ipywidgets installed.
    import ipywidgets  # noqa

    from .jupyter import show, show_input  # noqa
except ImportError:
    pass
