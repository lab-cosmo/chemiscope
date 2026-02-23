import asyncio
import atexit
import base64
import gzip
import json
import os
import threading
from pathlib import Path

from traitlets import Dict, HasTraits, observe

from .input import create_input


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

    def run_sync(self, coro):
        """Run a coroutine in the server loop and wait for the result."""
        future = asyncio.run_coroutine_threadsafe(coro, self._loop)
        return future.result()

    def stop(self):
        """Clean up resources and stop the loop."""
        self.run_sync(self._browser.close())
        self.run_sync(self._playwright.stop())
        self._loop.call_soon_threadsafe(self._loop.stop)
        self._thread.join()

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
    global _SERVER
    if _SERVER:
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

        # Gets a shared server instance (or creates one if it doesn't exist)
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

    def _run_sync(self, coro):
        return self._server.run_sync(coro)

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

        script = f"""
            async (dataset) => {{
                const config = {json.dumps(config)};
                window.visualizer = await Chemiscope.{vis_class}.load(config, dataset);
                return 'loaded';
            }}
        """

        self._run_sync(self._page.evaluate(script, data))
        # ensure initial settings are applied
        if self.settings:
            self._apply_settings(self.settings)

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
            self._run_sync(self._page.close())
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
    )

    return ChemiscopeHeadless(data, mode=mode, width=width, height=height, **kwargs)
