# -*- coding: utf-8 -*-
import atexit
import base64
import gzip
import json
import os
from pathlib import Path

from traitlets import Dict, HasTraits, observe

from .input import create_input


try:
    from playwright.sync_api import sync_playwright

    PLAYWRIGHT_AVAILABLE = True
except ImportError:
    PLAYWRIGHT_AVAILABLE = False

# global browser instance for reuse across multiple
# ChemiscopeHeadless instances
_PLAYWRIGHT = None
_BROWSER = None


def _get_browser():
    global _PLAYWRIGHT, _BROWSER
    if _PLAYWRIGHT is None:
        if not PLAYWRIGHT_AVAILABLE:
            raise ImportError(
                "chemiscope.headless requires the 'playwright' package. "
                "Please install it with `pip install playwright` and "
                "run `playwright install`."
            )
        _PLAYWRIGHT = sync_playwright().start()

        args = []
        if _should_disable_sandbox():
            args.append("--no-sandbox")
            args.append("--disable-setuid-sandbox")

        print("Launching headless browser...")
        _BROWSER = _PLAYWRIGHT.chromium.launch(args=args)
        print("Browser launched successfully.")
        # clean up browser on exit
        atexit.register(_close_browser)
    return _BROWSER


def _should_disable_sandbox():
    """
    Determine if the browser sandbox should be disabled.
    This is necessary in CI environments, containers, and when running as root.
    """
    # CI environments
    if os.environ.get("CI"):
        return True

    # Running as root (Chrome refuses to run as root with sandbox)
    if hasattr(os, "getuid") and os.getuid() == 0:
        return True

    # Running in a container (Docker/Podman)
    # Check for .dockerenv
    if os.path.exists("/.dockerenv"):
        return True

    # Check cgroups for docker/containerd strings
    try:
        with open("/proc/self/cgroup", "r") as f:
            if any("docker" in line or "containerd" in line for line in f):
                return True
    except (OSError, IOError):
        pass

    return False


def _close_browser():
    global _PLAYWRIGHT, _BROWSER
    if _BROWSER:
        _BROWSER.close()
        _BROWSER = None
    if _PLAYWRIGHT:
        _PLAYWRIGHT.stop()
        _PLAYWRIGHT = None


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

    def __init__(self, data, mode="default", **kwargs):
        super().__init__(**kwargs)

        # Gets a shared browser instance (or creates one if it doesn't exist)
        self._browser = _get_browser()
        self._page = self._browser.new_page()

        if "settings" in data:
            self.settings = data["settings"]

        # Load chemiscope library
        self._load_chemiscope_library()

        # Initialize visualizer
        self._initialize_visualizer(data, mode)
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
                "Cannot locate `chemiscope.min.js`. Please chech that the chemiscope"
                " package is installed and accessible."
            )

        # minimal chemiscope page with divs for the various components
        html = """
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <style>
        #chemiscope-structure, #chemiscope-map, #chemiscope-meta, #chemiscope-info {
            width: 800px; height: 600px;
        }
    </style>
</head>
<body>
    <div id="chemiscope-meta"></div>
    <div id="chemiscope-map"></div>
    <div id="chemiscope-structure"></div>
    <div id="chemiscope-info"></div>
</body>
</html>"""
        self._page.set_content(html)

        if content:
            self._page.add_script_tag(content=content)

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

        self._page.evaluate(script, data)
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
        self._page.evaluate(script, settings)

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
        self._page.evaluate(script, selected)

    def get_structure_image(self):
        """
        Request a snapshot of the active structure viewer.
        Returns the image data (PNG formatted) as bytes.
        """
        data_url = self._page.evaluate("""
            async () => {
                if (!window.visualizer || !window.visualizer.structure) {
                    throw new Error("No structure visualizer available");
                }
                // Wait a frame to ensure rendering
                await new Promise(r => requestAnimationFrame(r));
                return window.visualizer.structure.exportActivePNG();
            }
        """)
        return self._decode_base64(data_url)

    def get_map_image(self):
        """
        Request a snapshot of the map.
        Returns the image data (PNG formatted) as bytes.
        """
        data_url = self._page.evaluate("""
            async () => {
                if (!window.visualizer || !window.visualizer.map) {
                    throw new Error("No map visualizer available");
                }
                return await window.visualizer.map.exportPNG();
            }
        """)
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

        images = []
        for i, index in enumerate(indices):
            current_settings = settings[i] if settings else None

            # Prepare arguments
            args = {"index": index}
            if current_settings:
                args["settings"] = current_settings

            data_url = self._page.evaluate(
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
        if hasattr(self, "_page"):
            self._page.close()

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
    **kwargs,
):
    """
    Create a headless Chemiscope instance.
    Parameters match those of `chemiscope.show`.
    """
    data = create_input(
        structures=structures,
        properties=properties,
        metadata=metadata,
        environments=environments,
        shapes=shapes,
        settings=settings,
    )

    return ChemiscopeHeadless(data, mode=mode, **kwargs)
