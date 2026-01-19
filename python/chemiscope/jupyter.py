# -*- coding: utf-8 -*-
import asyncio
import base64
import gzip
import itertools
import json
import warnings
from pathlib import Path

import ipywidgets
from traitlets import Bool, Dict, Int, Unicode

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
    # timeout for warning messages
    warning_timeout = Int(10000).tag(sync=True)
    # currently selected environment/structure/atom
    selected_ids = Dict().tag(sync=True)
    # index of the currently active viewer in the grid
    active_viewer = Int(0).tag(sync=True)

    def __init__(
        self, data, has_metadata, warning_timeout=10000, cache_structures=True
    ):
        super().__init__()
        self.has_metadata = has_metadata
        if "settings" in data:
            self.settings = data["settings"]

        # timeout for warning messages (ms). 0 to make persistent, -1 to disable
        self.warning_timeout = warning_timeout

        self._request_counter = itertools.count()
        self._pending_requests = {}

        # hold structures on the python side to save js memory
        self._structures_cache = {}
        if cache_structures:
            data["structures"] = self._cache_structures(data["structures"])

        # pass data to javascript through a traitlet
        self.value = json.dumps(data)

        # listen for custom messages from the JS side
        self.on_msg(self._handle_js_msg)

    def _cache_structures(self, structures):
        """Cache full structures and return structure references"""
        if not structures or "data" in structures[0]:
            return structures

        structure_refs = []
        for i, structure in enumerate(structures):
            if "data" in structure:
                raise ValueError("Got a mix of explicit and external structures")

            cache_key = f"structure-{i}"
            self._structures_cache[cache_key] = structure
            structure_refs.append({"size": structure["size"], "data": cache_key})

        return structure_refs

    def get_map_image(self):
        """
        Request a snapshot of the map. Returns a Future that resolves to the
        image data (PNG formatted).

        This method is asynchronous. In a Jupyter notebook, you should ``await``
        it to get the data.

        .. code-block:: python

            # Get raw image data
            data = await widget.get_map_image()

            # Display it
            from IPython.display import Image

            display(Image(data))

        :return: A Future that resolves to the image data as bytes.
        """

        if self._view_name == "StructureView":
            raise RuntimeError(
                "Cannot retrieve map image: this widget is a structure-only viewer."
            )

        return self._request_screenshot("map")

    def get_structure_image(self):
        """
        Request a snapshot of the active structure viewer.
        Returns a Future that resolves to the image data (PNG formatted).

        This method is asynchronous. In a Jupyter notebook,
        you should ``await`` it to get the data.

        .. code-block:: python

            data = await widget.get_structure_image()

        :return: A Future that resolves to the image data as bytes.
        """

        if self._view_name == "MapView":
            raise RuntimeError(
                "Cannot retrieve structure image: this widget is a map-only viewer."
            )

        return self._request_screenshot("structure")

    def save_map_image(self, path):
        """
        Save a snapshot of the map to a file.

        This method starts a background task to save the
        image. You can ``await`` it if you need to ensure the
        file is written before proceeding.

        .. code-block:: python

            widget.save_map_image("map.png")

        :param str path: Path where the image will be saved.
        :return: A Future that resolves when the file is written.
        """

        if self._view_name == "StructureView":
            raise RuntimeError(
                "Cannot retrieve map image: this widget is a structure-only viewer."
            )

        return asyncio.ensure_future(self._save_image_to_file(path, "map"))

    def save_structure_image(self, path):
        """
        Save a snapshot of the active structure viewer to a file.

        This method starts a background task to save the image.
        You can ``await`` it if you need to ensure the file is
        written before proceeding.

        .. code-block:: python

            widget.save_structure_image("structure.png")

        :param str path: Path where the image will be saved.
        :return: A Future that resolves when the file is written.
        """

        if self._view_name == "MapView":
            raise RuntimeError(
                "Cannot save structure image: this widget is a map-only viewer."
            )

        return asyncio.ensure_future(self._save_image_to_file(path, "structure"))

    def get_structure_sequence(self, indices, settings=None):
        """
        Request a sequence of structure snapshots. Returns a Future that
        resolves to a list of image data (PNG formatted).

        This allows rendering multiple frames efficiently without blocking
        the UI. The sequence is processed in the browser.

        The ``indices`` list can contain integers (structure index) or
        dictionaries specifying ``structure`` and ``atom`` indices
        (for environments).

        The ``settings`` list, if provided, must have the same length
        as ``indices``. Each element is a dictionary of structure
        settings (e.g. ``{"spaceFilling": True}``) to apply for
        that specific frame.

        .. code-block:: python

            indices = [0, 1, 2]
            settings = [{"spaceFilling": True}, {}, {"spaceFilling": False}]
            data_list = await widget.get_structure_sequence(indices, settings)

        :param list indices: List of indices (int or dict) to render.
        :param list settings: Optional list of settings dicts to apply
            for each frame.
        :return: A Future that resolves to a list of image data bytes.
        """

        if self._view_name == "MapView":
            raise RuntimeError(
                "Cannot save structure image: this widget is a map-only viewer."
            )

        return asyncio.ensure_future(
            self._process_structure_sequence(indices, settings)
        )

    def save_structure_sequence(self, indices, paths, settings=None):
        """
        Save a sequence of structure snapshots to files.

        This method acts as a wrapper around :meth:`get_structure_sequence`
        and writes the results to files.

        .. code-block:: python

            indices = [0, 10, 20]
            paths = ["frame_0.png", "frame_10.png", "frame_20.png"]
            await widget.save_structure_sequence(indices, paths)

        :param list indices: List of indices (int or dict) to render.
        :param list paths: List of file paths where images will be saved.
            Must match length of ``indices``.
        :param list settings: Optional list of settings dicts to apply
            to each frame.
        :return: A Future that resolves when all files are written.
        """
        if len(indices) != len(paths):
            raise ValueError("indices and paths must have the same length")
        if settings is not None and len(indices) != len(settings):
            raise ValueError("indices and settings must have the same length")
        if self._view_name == "MapView":
            raise RuntimeError(
                "Cannot save structure image: this widget is a map-only viewer."
            )

        async def _save_impl():
            data_list = await self.get_structure_sequence(indices, settings)
            for path, data in zip(paths, data_list, strict=True):
                with open(path, "wb") as f:
                    f.write(data)

        return asyncio.ensure_future(_save_impl())

    def _request_screenshot(self, target):
        request_id = next(self._request_counter)

        loop = asyncio.get_running_loop()
        future = loop.create_future()
        self._pending_requests[request_id] = future

        self.send({"type": "save-image", "target": target, "requestId": request_id})
        return future

    async def _save_image_to_file(self, path, target):
        data = await self._request_screenshot(target)
        with open(path, "wb") as f:
            f.write(data)

    async def _process_structure_sequence(self, indices, settings=None):
        request_id = next(self._request_counter)

        loop = asyncio.get_running_loop()
        future = loop.create_future()

        self._pending_requests[request_id] = {
            "future": future,
            "results": [None] * len(indices),
            "type": "sequence",
        }

        msg = {
            "type": "save-structure-sequence",
            "indices": indices,
            "requestId": request_id,
        }
        if settings is not None:
            msg["settings"] = settings

        self.send(msg)

        return await future

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
        if len(self._structures_cache) > 0:
            # there are cached structures, add them to the dump
            data["structures"] = list(self._structures_cache.values())

        file.write(json.dumps(data).encode("utf8"))
        file.close()

    def _handle_js_msg(self, _, content, buffers):
        """Handle custom messages sent from the JS widget.
        The pattern is that a message is sent from another
        function, and the response is collected here, handling
        a successful response (`msg_type="base-mst-result"`)
        or and error (`msg_type="base-mst-result"`). See below
        for a brief explanation of the different types of messages.

        Currently implemented messages:

        * "load-structure": Initiated on the JS side, used to
            retrieve an atomic geometry from disk or from a
            python-side cache.

            {"type": "load-structure",
             "requestId": int,
             "data": "filename-or-id"}

        * "save-image": Initiated on the Python side, used to get
            a snapshot of the map or of the active structure

            {"type": "save-image",
             "requestId": int,
             "target": "map" | "structure"}

        * "save-structure-sequence": Initiated on the Python side, used to
            generate snapshots for several atomic configurations.
            Should provide indices (in the same format as for `selected_ids`)
            and optionally structure settings as a list of dicts.

            {"type": "save-structure-sequence",
             "requestId": int,
             "indices": [{"structure": id1[, "atom": id2]}, ...],
             ["settings": [settings_1_dict, ...]]
             }

        Expected messages:

        - {"type": "load-structure",
           "requestId": int,
           "filename": "relative/path.json[.gz]"}
        """
        msg_type = content.get("type", None)
        request_id = content.get("requestId")

        if msg_type == "save-image-result":
            pending = self._pending_requests.pop(request_id, None)
            if pending and isinstance(pending, asyncio.Future):
                data = content.get("data")
                try:
                    _header, encoded = data.split(",", 1)
                    data = base64.b64decode(encoded)
                    pending.set_result(data)
                except Exception as e:
                    pending.set_exception(e)

        elif msg_type == "save-image-error":
            pending = self._pending_requests.pop(request_id, None)
            if pending and isinstance(pending, asyncio.Future):
                pending.set_exception(RuntimeError(content.get("error")))

        elif msg_type == "save-structure-sequence-result":
            pending = self._pending_requests.get(request_id)
            if pending and isinstance(pending, dict):
                step = content.get("step")
                data = content.get("data")
                results = pending["results"]

                try:
                    _header, encoded = data.split(",", 1)
                    decoded = base64.b64decode(encoded)
                    results[step] = decoded
                except Exception as e:
                    print(f"Error decoding frame {step}: {e}")

        elif msg_type == "save-structure-sequence-done":
            pending = self._pending_requests.pop(request_id, None)
            if pending and isinstance(pending, dict):
                pending["future"].set_result(pending["results"])

        elif msg_type == "save-structure-sequence-error":
            pending = self._pending_requests.pop(request_id, None)
            if pending and isinstance(pending, dict):
                pending["future"].set_exception(RuntimeError(content.get("error")))

        if msg_type == "load-structure":
            # Reads a structure file and sends it back to the JS side
            data = content.get("data")

            if data is None:
                self.send(
                    {
                        "type": "load-structure-error",
                        "requestId": request_id,
                        "error": "missing filename in load-structure message",
                    }
                )
                return
            if data in self._structures_cache:
                # structure is cached
                self.send(
                    {
                        "type": "load-structure-result",
                        "requestId": request_id,
                        "structure": self._structures_cache[data],
                        "data": data,  # also sends back data to cross-check
                    }
                )
            else:
                # structure is an external file
                try:
                    # read structure from file (including relative path)
                    path = Path(data)

                    # Allow gzip-compressed structure files as well
                    if str(path).endswith(".gz"):
                        with gzip.open(path, "rt") as f:
                            structure = json.load(f)
                    else:
                        with path.open("rt") as f:
                            structure = json.load(f)
                    # `structure` must be a plain JSON-serialisable object that matches
                    # chemiscope's `Structure` interface on the JS side.
                    self.send(
                        {
                            "type": "load-structure-result",
                            "requestId": request_id,
                            "structure": structure,
                            "data": data,  # also sends back data to cross-check
                        }
                    )
                except Exception as exc:
                    self.send(
                        {
                            "type": "load-structure-error",
                            "requestId": request_id,
                            "error": f"failed to load structure from '{data}': {exc}",
                        }
                    )

    def __repr__(self, max_length=64):
        # string representation of the chemiscope widget, outputs that are too large
        class_name = self.__class__.__name__
        string_repr = ""

        # loops over metadata, structures, settings, ...
        for key, value in json.loads(self.value).items():
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

    def __init__(self, data, has_metadata, warning_timeout, cache_structures):
        super().__init__(data, has_metadata, warning_timeout, cache_structures)


@ipywidgets.register
class StructureWidget(ChemiscopeWidgetBase):
    _view_name = Unicode("StructureView").tag(sync=True)

    def __init__(self, data, has_metadata, warning_timeout, cache_structures):
        super().__init__(data, has_metadata, warning_timeout, cache_structures)


@ipywidgets.register
class MapWidget(ChemiscopeWidgetBase):
    _view_name = Unicode("MapView").tag(sync=True)

    def __init__(self, data, has_metadata, warning_timeout, cache_structures):
        super().__init__(data, has_metadata, warning_timeout, cache_structures)


def show_input(
    path, *, settings=None, mode="default", warning_timeout=10000, cache_structures=True
):
    """
    Loads and shows the chemiscope input in ``path``.

    If ``path`` ends with ``.gz``, the file is loaded as a gzip compressed JSON string.
    If ``path`` is a file-like object, it is read as JSON input.

    :param str | Path | file-like path: load the chemiscope input from this path or
        file-like object

    :param dict settings: override the default settings in the input
    :param str mode: widget mode, either ``default``, ``structure`` or ``map``.
    :param float warning_timeout: timeout (in ms) for warnings. Set to a negative value
        to disable warnings, and to zero to make them persistent.
    :param bool cache_structures: whether to cache structure data on the Python side
        to reduce the JScript memory footprint

    .. code-block:: python

        import chemiscope

        widget = chemiscope.show_input("dataset.json")

        # or

        with open("dataset.json", "r") as f:
            widget = chemiscope.show_input(f)

    ..
    """

    if mode not in ["default", "structure", "map"]:
        raise ValueError(
            f"invalid mode '{mode}' in chemiscope.show_input, expected one of "
            "'default', 'structure' or 'map'"
        )

    if mode == "default":
        widget_class = ChemiscopeWidget

    elif mode == "structure":
        widget_class = StructureWidget

    elif mode == "map":
        widget_class = MapWidget

    if isinstance(path, Path):
        path = str(path)

    if isinstance(path, str):
        if path.endswith(".gz"):
            with gzip.open(path, "rt") as f:
                dict_input = json.load(f)
        elif path.endswith(".json"):
            with open(path, "rt") as f:
                dict_input = json.load(f)
        else:
            raise ValueError(
                "invalid file format in chemiscope.show_input, expected .json or "
                ".json.gz"
            )
    else:
        dict_input = json.load(path)

    try:
        metadata = dict_input["meta"]

        if metadata == {"name": " "}:
            has_metadata = False
        else:
            has_metadata = True

    except KeyError:
        raise ValueError("missing metadata in chemiscope.show_input")

    if settings is not None:
        if not isinstance(settings, dict):
            raise ValueError(
                f"expected 'settings' to be a dict, got {type(settings)} instead"
            )
        if "settings" not in dict_input:
            dict_input["settings"] = {}
        dict_input["settings"].update(settings)

    return widget_class(
        dict_input,
        has_metadata=has_metadata,
        warning_timeout=warning_timeout,
        cache_structures=cache_structures,
    )


def show(
    structures=None,
    *,
    properties=None,
    metadata=None,
    environments=None,
    shapes=None,
    settings=None,
    mode="default",
    warning_timeout=10000,
    cache_structures=True,
    # for backward compatibility
    frames=None,
    meta=None,
):
    """
    Show the dataset defined by the given ``structures`` and ``properties`` (optionally
    ``metadata``, ``environments`` and ``shapes`` as well) using an embedded chemiscope
    visualizer inside a Jupyter notebook. These parameters have the same meaning as in
    the :py:func:`chemiscope.create_input` function.

    The ``mode`` keyword also allows overriding the default two-panels visualization to
    show only a structure panel (``mode = "structure"``) or the map panel (``mode =
    "map"``). These modes also make it possible to view a dataset for which properties
    (or structures) are not available. The widget displays warning messages, that
    disappear after the specified ``warning_timeout`` (in ms). Set to a negative value
    to disable warnings, and to zero to make them persistent. ``cache_structures`` is a
    flag determining whether to cache structure data on the Python side to reduce the
    JScript memory footprint.

    When inside a jupyter notebook, the returned object will create a new chemiscope
    visualizer displaying the dataset. The object exposes a ``settings`` traitlet, that
    allows to modify the visualization options (possibly even linking the parameters to
    another widget). Printing the value of the ``settings`` property is also a good way
    to see a full list of the available options.

    The returned object also have a ``save`` function that can be used to save the
    dataset to a ``.json`` or ``.json.gz`` file to load it in the main website later.
    The visualization options will be those used in the active widget, so this is also a
    good way to tweak the appearance of the visualization before saving it.

    .. code-block:: python

        import chemiscope
        from sklearn.decomposition import PCA
        import ase.io

        pca = PCA(n_components=3)

        structures = ase.io.read(...)
        properties = {
            "PCA": pca.fit_transform(some_data),
        }

        widget = chemiscope.show(structures, properties)
        # display the dataset in a chemiscope visualizer inside the notebook
        widget
        # ...

        # NB: due to how traitlet work, you should always set the value of
        # the `settings` property. Only the properties that are explicitly
        # indicated will be modified.
        widget.settings = {"map": {"symbol": "tag"}}
        widget.settings["map"]["symbol"] = "tag"  # << does nothing!

        # Save the file for later use
        widget.save("dataset.json")

    .. _ase.Atoms: https://wiki.fysik.dtu.dk/ase/ase/atoms.html
    """
    if frames is not None:
        warnings.warn(
            "`frames` argument is deprecated, use `structures` instead",
            stacklevel=2,
        )
        if structures is not None:
            raise ValueError("cannot use both `structures` and `frames` arguments")

        structures = frames

    if meta is not None:
        warnings.warn(
            "`meta` argument is deprecated, use `metadata` instead",
            stacklevel=2,
        )
        if metadata is not None:
            raise ValueError("cannot use both `metadata` and `meta` arguments")

        metadata = meta

    if not (_is_running_in_notebook() or _is_running_in_sphinx_gallery()):
        warnings.warn(
            "chemiscope.show only works in a jupyter notebook or a sphinx build",
            stacklevel=2,
        )

    has_metadata = metadata is not None
    if not has_metadata:
        metadata = {"name": " "}

    if mode == "default":
        widget_class = ChemiscopeWidget

    elif mode == "structure":
        if properties is None:
            properties = {}

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
        structures=structures,
        properties=properties,
        metadata=metadata,
        environments=environments,
        shapes=shapes,
        settings=settings,
    )

    if mode == "structure" and "index" not in dict_input["properties"]:
        dict_input["properties"]["index"] = {
            "target": "structure",
            "values": list(range(len(dict_input["structures"]))),
        }

    if mode != "structure":
        # if there is a map, we need two properties, otherwise there will be no map
        # and error is only visible in the console
        if len(dict_input["properties"]) < 2:
            raise ValueError("Need at least two properties to visualize a map widget")

    return widget_class(
        dict_input,
        has_metadata=has_metadata,
        warning_timeout=warning_timeout,
        cache_structures=cache_structures,
    )


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
