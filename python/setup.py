#!/usr/bin/env python3
from setuptools import setup
from distutils.command.build import build
from distutils import log
import json
import os
import subprocess


class CustomBuild(build):
    def run(self):
        log.info("running npm --prefix ../ run build")
        print("path", os.path.abspath(os.getcwd()))
        print("files")
        files = [f for f in os.listdir("../") if os.path.isfile(f)]
        for f in files:
            print(f)
        try:
            os.system("npm --prefix ../ run build")
        except subprocess.CalledProcessError as e:
            # os.system("npm run build")
            os.system("npm --prefix ../../../ run build")
        build.run(self)


if __name__ == "__main__":
    # get the version directly from the package.json file
    package_json = os.path.join(os.path.dirname(__file__), "package.json")
    version = json.load(open(package_json))["version"]
    setup(
        version=version,
        cmdclass={"build": CustomBuild},
        data_files=[
            # like `jupyter nbextension install --sys-prefix`
            # TODO: fix the above, does not seem to work and has to be done manually
            # cf Chiheb's and Michele's computer
            (
                "share/jupyter/nbextensions/chemiscope-widget",
                [
                    "chemiscope/nbextension/extension.js",
                    "chemiscope/nbextension/chemiscope-widget.min.js",
                ],
            ),
            # like `jupyter nbextension enable --sys-prefix`
            # TODO: fix the above, does not seem to work and has to be done manually
            # cf Chiheb's and Michele's computer
            (
                "etc/jupyter/nbconfig/notebook.d",
                ["chemiscope/nbextension/chemiscope-widget.json"],
            ),
        ],
    )
