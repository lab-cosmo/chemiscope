#!/usr/bin/env python3
import glob
import json
import os
import subprocess
import sys
import warnings

from setuptools import setup

# customize dist directory if the user did not already request one, and if we
# are running setup.py directly (i.e. not from pip)
if not sys.argv[0].endswith("setup.py"):
    change_dist_dir = False

change_dist_dir = True
for arg in sys.argv:
    if "--dist-dir" in arg:
        change_dist_dir = False


def run_npm_build():
    """Build the JavaScript code required by the jupyter widget using npm"""
    root = os.path.dirname(os.path.realpath(__file__))

    # This file does not exists when building a wheel (for pip installation)
    # from a sdist. This is fine since the sdist already contains the required
    # js files.
    if os.path.exists(os.path.join(root, "package.json")):
        subprocess.run("npm ci", check=True, shell=True)
        # build the labextension first, since it triggers a rebuild of the main
        # package code
        try:
            import jupyterlab  # noqa

            subprocess.run("npm run build:labextension", check=True, shell=True)
            with open(
                os.path.join(root, "python", "jupyter", "labextension", "install.json"),
                "w",
            ) as fd:
                install_json = {
                    "__comment": "metadata for chemiscope jupyterlab extension",
                    "packageManager": "python",
                    "packageName": "chemiscope",
                    "uninstallInstructions": "Use your Python package manager "
                    + "(pip, conda, etc.) to uninstall chemiscope",
                }
                json.dump(install_json, fd)
        except ImportError:
            # skip building jupyterlab if it is not installed on the developer
            # machine
            warnings.warn(
                "skipping the lab extension build since jupyterlab is not installed"
            )

        subprocess.run("npm run build:nbextension", check=True, shell=True)


if __name__ == "__main__":
    # we need to run npm build outside of the call to setup to be able to get
    # the full list of files to install for the labextension. Building the
    # labextension creates unpredictable file names, and `data_file` wants to
    # explicitly list all files one by one.

    # TODO: this currently runs more than once because pip runs this function
    # multiple times. Can we cache the outcome?
    run_npm_build()

    setup(
        data_files=[
            # this is what `jupyter nbextension install --sys-prefix` does
            (
                "share/jupyter/nbextensions/chemiscope",
                [
                    "python/jupyter/extension.js",
                    "python/jupyter/nbextension/chemiscope.min.js",
                ],
            ),
            # this is what `jupyter nbextension enable --sys-prefix` does
            (
                "etc/jupyter/nbconfig/notebook.d",
                [
                    "python/jupyter/chemiscope.json",
                ],
            ),
            # install the labextension
            (
                "share/jupyter/labextensions/chemiscope",
                list(
                    filter(
                        lambda f: os.path.isfile(f),
                        glob.glob("python/jupyter/labextension/*"),
                    )
                ),
            ),
            (
                "share/jupyter/labextensions/chemiscope/static",
                glob.glob("python/jupyter/labextension/static/*"),
            ),
        ],
    )
