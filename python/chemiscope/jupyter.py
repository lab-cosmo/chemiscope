# -*- coding: utf-8 -*-
"""Deprecated alias for :mod:`chemiscope.widget`.

The widget module was renamed from ``jupyter`` to ``widget``. Import from
``chemiscope.widget`` (or use the top-level ``chemiscope.show``) instead.
"""

import warnings

from .widget import (  # noqa: F401
    ChemiscopeWidget,
    ChemiscopeWidgetBase,
    MapWidget,
    StructureWidget,
    show,
    show_input,
)


warnings.warn(
    "chemiscope.jupyter is deprecated and will be removed in a future release; "
    "use chemiscope.widget instead",
    DeprecationWarning,
    stacklevel=2,
)
