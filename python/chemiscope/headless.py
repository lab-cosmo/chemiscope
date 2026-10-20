import asyncio
import atexit
import base64
import gzip
import http.server
import json
import math
import os
import secrets
import threading
import warnings
from pathlib import Path

from traitlets import Dict, HasTraits, observe

from .input import _read_input, create_input


try:
    from playwright.async_api import async_playwright

    PLAYWRIGHT_AVAILABLE = True
except ImportError:
    PLAYWRIGHT_AVAILABLE = False

# global PlaywrightServer instance for reuse across multiple
# ChemiscopeHeadless instances
_SERVER = None


class PlaywrightServer:
    """
    A dedicated thread running an asyncio loop to manage Playwright.
    This avoids conflicts with existing loops in the main thread.
    """

    def __init__(self):
        if not PLAYWRIGHT_AVAILABLE:
            raise ImportError(
                "chemiscope.headless requires the 'playwright' package. "
                "Please install it with `pip install playwright` and "
                "run `playwright install`."
            )

        self._running = True
        self._loop = asyncio.new_event_loop()
        self._thread = threading.Thread(target=self._run_loop, daemon=True)
        self._thread.start()

        self._playwright = self.run_sync(async_playwright().start())

        # enable software rendering to support headless environments without GPU
        args = ["--enable-unsafe-swiftshader"]
        self._browser = self.run_sync(self._playwright.chromium.launch(args=args))

    def _run_loop(self):
        asyncio.set_event_loop(self._loop)
        self._loop.run_forever()

    def run_sync(self, coro, timeout=10):
        """Run a coroutine in the server loop and wait for the result."""
        if not self._running:
            raise RuntimeError("PlaywrightServer has been stopped")
        future = asyncio.run_coroutine_threadsafe(coro, self._loop)
        # timeout prevents hanging if the browser process died
        return future.result(timeout=timeout)

    def stop(self):
        """Clean up resources and stop the loop."""
        if not self._running:
            return
        self._running = False
        # browser/playwright may already be dead; warn but don't raise
        for coro in [self._browser.close(), self._playwright.stop()]:
            try:
                future = asyncio.run_coroutine_threadsafe(coro, self._loop)
                future.result(timeout=5)
            except Exception as e:
                warnings.warn(
                    f"error during PlaywrightServer cleanup: {e}",
                    stacklevel=2,
                )
        self._loop.call_soon_threadsafe(self._loop.stop)
        self._thread.join(timeout=5)

    def new_page(self, **kwargs):
        """Create a new page in the browser."""
        return self.run_sync(self._browser.new_page(**kwargs))


def _get_server():
    global _SERVER
    if _SERVER is None:
        _SERVER = PlaywrightServer()
        # clean up server on exit
        atexit.register(_close_server)
    return _SERVER


def _close_server():
    """Shut down the server and discard it so a fresh one is created next time."""
    global _SERVER
    if _SERVER is not None:
        _SERVER.stop()
        _SERVER = None


class ChemiscopeHeadless(HasTraits):
    """
    A headless version of the Chemiscope widget, running in a managed browser
    instance. This allows for programmatic generation of screenshots without
    a Jupyter environment.
    """

    # synchronized settings
    settings = Dict()
    # currently selected environment/structure/atom
    selected_ids = Dict()

    def __init__(self, data, mode="default", width=None, height=None, **kwargs):
        super().__init__(**kwargs)

        # Handle dimensions with 4:3 ratio default
        if width is None and height is None:
            width = 800
            height = 600
        elif width is None:
            width = int(height * 4 / 3)
        elif height is None:
            height = int(width * 3 / 4)

        self._width = width
        self._height = height

        # Get a shared server, or create one. If the browser has crashed
        # (common in long-running Jupyter sessions), retry with a fresh server.
        self._server = _get_server()
        try:
            self._page = self._server.new_page(device_scale_factor=1)
        except Exception as e:
            warnings.warn(
                f"browser connection lost ({e}), restarting server",
                stacklevel=2,
            )
            _close_server()
            self._server = _get_server()
            self._page = self._server.new_page(device_scale_factor=1)
        self._run_sync(
            self._page.set_viewport_size({"width": self._width, "height": self._height})
        )

        if "settings" in data:
            self.settings = data["settings"]

        # Load chemiscope library
        self._load_chemiscope_library()

        # Initialize visualizer
        self._initialize_visualizer(data, mode)
        self._mode = mode
        self._data = data

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
        data = self._data.copy()
        data["settings"] = self.settings

        file.write(json.dumps(data).encode("utf8"))
        file.close()

    def _run_sync(self, coro, timeout=10):
        return self._server.run_sync(coro, timeout=timeout)

    def _load_chemiscope_library(self):
        # look for the library in the install path
        # (we reuse the `chemiscope.min.js` stored for chemiscope.sphinx)
        package_static = (
            Path(os.path.dirname(__file__)) / "sphinx" / "static" / "chemiscope.min.js"
        )

        content = None
        if package_static.exists():
            with open(package_static, "r") as f:
                content = f.read()
        else:
            raise ImportError(
                "Cannot locate `chemiscope.min.js`. Please check that the chemiscope"
                " package is installed and accessible."
            )

        # minimal chemiscope page with divs for the various components
        html = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <style>
        body {{ margin: 0; padding: 0; overflow: hidden; }}
        #chemiscope-structure, #chemiscope-map, #chemiscope-meta, #chemiscope-info {{
            width: {self._width}px; height: {self._height}px;
        }}
    </style>
</head>
<body>
    <div id="chemiscope-meta"></div>
    <div id="chemiscope-map"></div>
    <div id="chemiscope-structure"></div>
    <div id="chemiscope-info"></div>
</body>
</html>"""
        self._run_sync(self._page.set_content(html))

        if content:
            self._run_sync(self._page.add_script_tag(content=content))

    def _initialize_visualizer(self, data, mode):
        # determine visualizer class based on mode
        if mode == "default":
            vis_class = "DefaultVisualizer"
            config_keys = ["meta", "map", "structure", "info"]
        elif mode == "structure":
            vis_class = "StructureVisualizer"
            config_keys = ["meta", "structure", "info"]
        elif mode == "map":
            vis_class = "MapVisualizer"
            config_keys = ["meta", "map", "info"]
        else:
            raise ValueError(f"Invalid mode: {mode}")

        config = {key: f"chemiscope-{key}" for key in config_keys}

        # Keep the structures in Python and send them to the browser when they
        # are displayed, instead of sending the whole dataset upfront, which
        # would fail for large datasets
        self._structures = data["structures"]
        dataset = dict(data)
        dataset["structures"] = [
            {"size": structure["size"], "data": i}
            for i, structure in enumerate(self._structures)
        ]
        self._run_sync(
            self._page.expose_function("__chsp_load_structure", self._load_structure)
        )

        # Send the rest of the dataset as JSON, which is much faster to parse in
        # the browser than the equivalent JavaScript code
        dataset_json, has_non_finite = _to_json(dataset)
        server, url = _serve(dataset_json.encode("utf8"))

        script = f"""
            async () => {{
                const response = await fetch({json.dumps(url)});
                // non-finite values are sent as strings, see `_to_json`
                const nonFinite = new Map([
                    ['__chemiscope_nan__', NaN],
                    ['__chemiscope_inf__', Infinity],
                    ['__chemiscope_-inf__', -Infinity],
                ]);
                const reviver = (_, value) =>
                    nonFinite.has(value) ? nonFinite.get(value) : value;
                const dataset = JSON.parse(
                    await response.text(),
                    {json.dumps(has_non_finite)} ? reviver : undefined
                );

                const config = {json.dumps(config)};
                if (config.structure !== undefined) {{
                    config.loadStructure = (_, structure) =>
                        window.__chsp_load_structure(structure.data);
                }}
                window.visualizer = await Chemiscope.{vis_class}.load(config, dataset);
                return 'loaded';
            }}
        """

        try:
            # use a large timeout as it might take some time to load the dataset
            self._run_sync(self._page.evaluate(script), timeout=300)
        finally:
            server.shutdown()
            server.server_close()

        # ensure initial settings are applied
        if self.settings:
            self._apply_settings(self.settings)

    def _load_structure(self, index):
        structure = self._structures[index]
        if "data" in structure:
            # structure stored in an external file, see `write_external_structures`
            return _read_input(structure["data"], "chemiscope.headless")
        return structure

    # connect traitlets
    @observe("settings")
    def _on_settings_change(self, change):
        self._apply_settings(change.new)

    def _apply_settings(self, settings):
        script = """
            (settings) => {
                if (window.visualizer) {
                    window.visualizer.applySettings(settings);
                }
            }
        """
        self._run_sync(self._page.evaluate(script, settings))

    @observe("selected_ids")
    def _on_selected_ids_change(self, change):
        selected = change.new
        if not selected:
            return

        script = """
            ({structure, atom}) => {
                if (!window.visualizer) return;
                
                const indexer = window.visualizer.indexer;
                const target = window.visualizer.saveSettings().target;
                
                let indexes;
                if (atom !== undefined) {
                    indexes = indexer.fromStructureAtom(target, structure, atom);
                } else {
                    indexes = indexer.fromStructure(structure, target);
                }
                
                if (indexes) {
                    window.visualizer.select(indexes);
                }
            }
        """
        self._run_sync(self._page.evaluate(script, selected))

    def get_structure_image(self):
        """
        Request a snapshot of the active structure viewer.
        Returns the image data (PNG formatted) as bytes.
        """
        if self._mode == "map":
            raise RuntimeError(
                "Cannot retrieve structure image: this widget is a map-only viewer."
            )

        data_url = self._run_sync(
            self._page.evaluate(
                """
            async () => {
                if (!window.visualizer || !window.visualizer.structure) {
                    throw new Error("No structure visualizer available");
                }
                // Wait a frame to ensure rendering
                await new Promise(r => requestAnimationFrame(r));
                return window.visualizer.structure.exportActivePNG();
            }
        """
            )
        )
        return self._decode_base64(data_url)

    def get_map_image(self):
        """
        Request a snapshot of the map.
        Returns the image data (PNG formatted) as bytes.
        """
        if self._mode == "structure":
            raise RuntimeError(
                "Cannot retrieve map image: this widget is a structure-only viewer."
            )

        data_url = self._run_sync(
            self._page.evaluate(
                """
            async () => {
                if (!window.visualizer || !window.visualizer.map) {
                    throw new Error("No map visualizer available");
                }
                return await window.visualizer.map.exportPNG();
            }
        """
            )
        )
        return self._decode_base64(data_url)

    def save_structure_image(self, path):
        """
        Save a snapshot of the active structure viewer to a file.
        """

        data = self.get_structure_image()
        with open(path, "wb") as f:
            f.write(data)

    def save_map_image(self, path):
        """
        Save a snapshot of the map to a file.
        """
        data = self.get_map_image()
        with open(path, "wb") as f:
            f.write(data)

    def get_structure_sequence(self, indices, settings=None):
        """
        Request a sequence of structure snapshots.
        Returns a list of image data (PNG formatted) as bytes.

        The ``indices`` list can contain integers (structure index) or
        dictionaries specifying ``structure`` and ``atom`` indices
        (for environments).

        The ``settings`` list, if provided, must have the same length
        as ``indices``. Each element is a dictionary of structure
        settings (e.g. ``{"spaceFilling": True}``) to apply for
        that specific frame.
        """
        if settings is not None and len(indices) != len(settings):
            raise ValueError("indices and settings must have the same length")

        if self._mode == "map":
            raise RuntimeError(
                "Cannot retrieve structure image: this widget is a map-only viewer."
            )

        images = []
        for i, index in enumerate(indices):
            current_settings = settings[i] if settings else None

            # Prepare arguments
            args = {"index": index}
            if current_settings:
                args["settings"] = current_settings

            data_url = self._run_sync(
                self._page.evaluate(
                    """
                async ({index, settings}) => {
                    const vis = window.visualizer;
                    const indexer = vis.indexer;

                    if (settings) {
                        vis.structure.applySettings([settings]);
                    }

                    const target = vis.saveSettings().target;
                    let indexes;

                    // Handle index being int or object
                    let structIdx = typeof index === 'number' ? index : index.structure;
                    let atomIdx = typeof index === 'number' ? undefined : index.atom;

                    if (atomIdx !== undefined) {
                        indexes = indexer.fromStructureAtom(target, structIdx, atomIdx);
                    } else {
                        indexes = indexer.fromStructure(structIdx, target);
                    }

                    await vis.structure.show(indexes);
                    await new Promise(r => requestAnimationFrame(r));

                    return vis.structure.exportActivePNG();
                }
            """,
                    args,
                )
            )

            images.append(self._decode_base64(data_url))

        return images

    def save_structure_sequence(self, indices, paths, settings=None):
        """
        Save a sequence of structure snapshots to files.
        """

        if len(indices) != len(paths):
            raise ValueError("indices and paths must have the same length")

        images = self.get_structure_sequence(indices, settings)
        for path, data in zip(paths, images, strict=True):
            with open(path, "wb") as f:
                f.write(data)

    def _decode_base64(self, data_url):
        if "," in data_url:
            base64_data = data_url.split(",")[1]
        else:
            base64_data = data_url

        return base64.b64decode(base64_data)

    def close(self):
        """Close the browser page."""
        if hasattr(self, "_page") and self._page is not None:
            # skip if the server is already shut down (avoids hangs in __del__)
            if hasattr(self, "_server") and self._server._running:
                try:
                    self._run_sync(self._page.close())
                except Exception:
                    pass
            self._page = None

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass


def headless(
    structures=None,
    *,
    properties=None,
    metadata=None,
    environments=None,
    shapes=None,
    settings=None,
    parameters=None,
    mode="default",
    width=None,
    height=None,
    **kwargs,
):
    """
    Create a headless Chemiscope instance.
    Parameters match those of `chemiscope.show`.

    :param int width: width of the widget container in pixels. If not specified,
        it is calculated from the height (defaulting to 800 if both are missing)
        to keep a 4:3 aspect ratio. Note that the output image resolution might
        differ from this value (e.g. due to high-DPI scaling or internal padding).
    :param int height: height of the widget container in pixels. If not specified,
        it is calculated from the width (defaulting to 600 if both are missing)
        to keep a 4:3 aspect ratio.
    """
    data = create_input(
        structures=structures,
        properties=properties,
        metadata=metadata,
        environments=environments,
        shapes=shapes,
        settings=settings,
        parameters=parameters,
    )

    return ChemiscopeHeadless(data, mode=mode, width=width, height=height, **kwargs)


def headless_input(
    path, *, settings=None, mode="default", width=None, height=None, **kwargs
):
    """
    Create a headless Chemiscope instance from the chemiscope input in ``path``,
    which can be a ``.json`` or ``.json.gz`` file, or a file-like object.

    :param dict settings: override the default settings in the input, as in
        :py:func:`chemiscope.show_input`

    The other parameters match those of :py:func:`chemiscope.headless`.
    """
    data = _read_input(path, "chemiscope.headless_input")

    if settings is not None:
        if not isinstance(settings, dict):
            raise ValueError(
                f"expected 'settings' to be a dict, got {type(settings)} instead"
            )
        data.setdefault("settings", {}).update(settings)

    return ChemiscopeHeadless(data, mode=mode, width=width, height=height, **kwargs)


def _serve(body):
    """
    Serve ``body`` from a local HTTP server, and return the server and the URL
    of the data. This is much faster than sending large data to the browser
    through Playwright, e.g. with ``page.route``.
    """
    path = "/" + secrets.token_urlsafe(16)

    class Handler(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            if self.path != path:
                self.send_error(404)
                return
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Access-Control-Allow-Origin", "*")
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, *args):
            pass

    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server, f"http://127.0.0.1:{server.server_address[1]}{path}"


def _to_json(data):
    """
    Convert ``data`` to JSON. NaN and infinite values are not valid JSON, and are
    replaced by strings. The second returned value is ``True`` if this happened.
    """
    try:
        return json.dumps(data, allow_nan=False), False
    except ValueError:
        return json.dumps(_replace_non_finite(data), allow_nan=False), True


def _replace_non_finite(value):
    if isinstance(value, float) and not math.isfinite(value):
        # "nan", "inf" or "-inf"
        return f"__chemiscope_{value}__"
    elif isinstance(value, dict):
        return {key: _replace_non_finite(v) for key, v in value.items()}
    elif isinstance(value, (list, tuple)):
        return [_replace_non_finite(v) for v in value]
    return value
