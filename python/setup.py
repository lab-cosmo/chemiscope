#!/usr/bin/env python3
from setuptools import setup
import json
import os

if __name__ == "__main__":
    # get the version directly from the package.json file
    package_json = os.path.join(os.path.dirname(__file__), "package.json")
    version = json.load(open(package_json))["version"]
    setup(version=version)
