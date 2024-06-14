#!/usr/bin/env python3
import glob
import json
import os
import shutil
import subprocess
import sys
import warnings

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
    *glob.glob("python/jupyter/src/**/*", recursive=True),
    # dependencies
    "package.json",
    "package-lock.json",
    # build configuration
    "tsconfig.json",
    "webpack.config.ts",
]

NPM_BUILD_OUTPUT = [
    "python/jupyter/nbextension/chemiscope.min.js",
    "python/jupyter/labextension/package.json",
    "python/chemiscope/sphinx/static/chemiscope.min.js",
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
    """Build the JavaScript code required by the jupyter widget using npm"""
    root = os.path.dirname(os.path.realpath(__file__))

    # This file does not exist when building a wheel (for pip installation)
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
                "skipping the lab extension build since jupyterlab is not installed",
                stacklevel=1,
            )

            # still build the main code for the sphinx extension
            subprocess.run("npm run build", check=True, shell=True)

        shutil.copyfile(
            src="dist/chemiscope.min.js",
            dst="python/chemiscope/sphinx/static/chemiscope.min.js",
        )

        subprocess.run("npm run build:nbextension", check=True, shell=True)


if __name__ == "__main__":
    # we need to run npm build outside of the call to setup to be able to get
    # the full list of files to install for the labextension. Building the
    # labextension creates unpredictable file names, and `data_file` wants to
    # explicitly list all files one by one.

    if needs_npm_build():
        print("==== rebuilding the JS code with npm, this might take a while ...")
        run_npm_build()

    setup(
        cmdclass={
            "bdist_egg": bdist_egg if "bdist_egg" in sys.argv else bdist_egg_disabled,
        },
        package_data={"chemiscope": ["chemiscope/sphinx/static/*"]},
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
