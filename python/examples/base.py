"""
Simple chemiscope input
=======================

This example demonstrates the basic usage of the 
chemiscope package, reading structures and properties
from an ASE package and preparing a chemiscope file
to visualize them. Use `chemiscope.show` in a 
Jupyter notebook for interactive visualization

"""

import numpy as np
import ase.io
import chemiscope

# load structures
frames = ase.io.read("data/showcase.xyz", ":")

# write the chemiscope input file. use `chemiscope.show` to view directly in a Jupyter environment
chemiscope.write_input(
    path="showcase.json.gz",
    frames=frames,
    # quickly extract properties from the ASE frames
    properties=chemiscope.extract_properties(frames, only=["dipole_ccsd", "ccsd_pol"]),
    meta=dict(name="Dipole and polarizability"),
    settings={
        "map": {
            "x": {"property": "ccsd_pol[1]"},
            "y": {"property": "ccsd_pol[2]"},
            "color": {"property": "dipole_ccsd[1]"},
        }
    },
)
