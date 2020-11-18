# -*- coding: utf-8 -*-
from .input import write_input


# Get the version from setuptools metadata (which took it from package.json)
# cf https://stackoverflow.com/a/17638236/4692076
from pkg_resources import get_distribution, DistributionNotFound
import os

try:
    dist = get_distribution("chemiscope")
    # Normalize case for Windows systems
    dist_loc = os.path.normcase(dist.location)
    here = os.path.normcase(__file__)
    if not here.startswith(os.path.join(dist_loc, "chemiscope")):
        # not installed, but there is another version that *is*
        raise DistributionNotFound
except DistributionNotFound:
    __version__ = "dev"
else:
    __version__ = dist.version
