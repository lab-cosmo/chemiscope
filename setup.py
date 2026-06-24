#!/usr/bin/env python3
import glob
import os
import shutil
import subprocess
import sys

from setuptools import setup
from setuptools.command.bdist_egg import bdist_egg


class bdist_egg_disabled(bdist_egg):
    """Disabled version of bdist_egg

    Prevents setup.py install performing setuptools' default easy_install,
    which it should never ever do.
    """

    def run(self):
        sys.exit(
            "Aborting implicit building of eggs. "
            + "Use `pip install .` or `python setup.py bdist_wheel && pip "
            + "uninstall chemiscope -y && pip install dist/chemiscope-*.whl` "
            + "to install from source."
        )


NPM_BUILD_INPUT = [
    # source files
    *glob.glob("src/**/*", recursive=True),
    *glob.glob("python/widget/src/**/*", recursive=True),
    # dependencies
    "package.json",
    "package-lock.json",
    # build configuration
    "tsconfig.json",
    "webpack.config.ts",
    # streamlit component
    *glob.glob("python/streamlit/src/**/*", recursive=True),
]

NPM_BUILD_OUTPUT = [
    # anywidget ES module for the notebook widget
    "python/chemiscope/static/chemiscope-widget.mjs",
    # sphinx extension
    "python/chemiscope/sphinx/static/chemiscope.min.js",
    # streamlit component
    "python/chemiscope/streamlit/index.html",
    "python/chemiscope/streamlit/chemiscope-streamlit.min.js",
]


def needs_npm_build():
    """
    Determine if we need to re-build the JS code with npm by checking the modification
    time for inputs & outputs of this build.
    """
    last_output_modification_time = -1
    for file in NPM_BUILD_OUTPUT:
        if os.path.exists(file):
            last_output_modification_time = max(
                os.path.getmtime(file), last_output_modification_time
            )
        else:
            # if the file does not exist, we need to build
            return True

    last_input_modification_time = -1
    for file in NPM_BUILD_INPUT:
        if os.path.exists(file):
            last_input_modification_time = max(
                os.path.getmtime(file), last_input_modification_time
            )
        # if the file does not exist, we might be building a sdist which already
        # contains the pre-built extensions

    return last_output_modification_time < last_input_modification_time


def run_npm_build():
    """Build the JavaScript code required by the notebook widget using npm"""
    root = os.path.dirname(os.path.realpath(__file__))

    # This file does not exist when building a wheel (for pip installation)
    # from a sdist. This is fine since the sdist already contains the required
    # js files.
    if os.path.exists(os.path.join(root, "package.json")):
        subprocess.run("npm ci", check=True, shell=True)

        # build the main code, used by the sphinx extension
        subprocess.run("npm run build", check=True, shell=True)
        shutil.copyfile(
            src="dist/chemiscope.min.js",
            dst="python/chemiscope/sphinx/static/chemiscope.min.js",
        )

        # build the anywidget ES module for the notebook widget
        subprocess.run("npm run build:anywidget", check=True, shell=True)

        subprocess.run("npm run build:streamlit", check=True, shell=True)

        streamlit_output_dir = os.path.join("python", "chemiscope", "streamlit")
        os.makedirs(streamlit_output_dir, exist_ok=True)
        shutil.copyfile(
            os.path.join("python", "streamlit", "build", "chemiscope-streamlit.min.js"),
            os.path.join(streamlit_output_dir, "chemiscope-streamlit.min.js"),
        )

        shutil.copyfile(
            os.path.join("python", "streamlit", "src", "index.html"),
            os.path.join(streamlit_output_dir, "index.html"),
        )


if __name__ == "__main__":
    if needs_npm_build():
        print("==== rebuilding the JS code with npm, this might take a while ...")
        run_npm_build()

    # The anywidget ES module is loaded directly by anywidget through the `_esm`
    # trait, so it only needs to be installed as package data (no notebook
    # frontend-extension registration is required).
    setup(
        cmdclass={
            "bdist_egg": bdist_egg if "bdist_egg" in sys.argv else bdist_egg_disabled,
        },
        package_data={
            "chemiscope": [
                "static/*",
                "sphinx/static/*",
                "streamlit/**/*",
            ]
        },
    )
