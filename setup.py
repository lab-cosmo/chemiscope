#!/usr/bin/env python3
import os
import subprocess
import sys
from distutils.command import install, sdist

from setuptools import setup
from wheel import bdist_wheel

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
        subprocess.run("npm run build:nbextension", check=True, shell=True)


class sdist_with_npm(sdist.sdist):
    def run(self):
        if change_dist_dir:
            self.dist_dir = os.path.join("dist", "python")
        run_npm_build()
        super().run()


class bdist_wheel_with_npm(bdist_wheel.bdist_wheel):
    def run(self):
        run_npm_build()
        if change_dist_dir:
            self.dist_dir = os.path.join("dist", "python")
        super().run()


class install_with_npm(install.install):
    def run(self):
        run_npm_build()
        super().run()


if __name__ == "__main__":
    setup(
        cmdclass={
            # we need to override multiple classes here to cover all possible
            # use cases. `sdist` is used to generate the source upload to PyPI,
            # `bdist_wheel` is used both to generate the wheel from the sdist
            # (pip install from PyPI) and the wheel from the full repository
            # (pip install with local path); and `install` adds support for
            # `python setup.py install`. Notablely missing here is `python
            # setup.py develop`/`pip install -e <local/path>`.
            "sdist": sdist_with_npm,
            "bdist_wheel": bdist_wheel_with_npm,
            "install": install_with_npm,
        },
        data_files=[
            # this is what `jupyter nbextension install --sys-prefix` does
            (
                "share/jupyter/nbextensions/chemiscope-widget",
                [
                    "python/nbextension/extension.js",
                    "python/nbextension/build/chemiscope-widget.min.js",
                ],
            ),
            # this is what `jupyter nbextension enable --sys-prefix` does
            (
                "etc/jupyter/nbconfig/notebook.d",
                [
                    "python/nbextension/chemiscope-widget.json",
                ],
            ),
        ],
    )
