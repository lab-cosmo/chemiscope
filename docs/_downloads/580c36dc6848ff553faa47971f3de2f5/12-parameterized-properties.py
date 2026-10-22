"""
Parameterized properties
========================

This example demonstrates the use of parameterized (multidimensional) properties
in chemiscope. These are displayed as curves in the info panel.

We compute the radial distribution function (RDF) for each structure and display
it as a parameterized property alongside the structure viewer.
"""

# %%
#

import ase.io
import numpy as np
from ase.geometry.rdf import get_rdf

import chemiscope


# %%
#
# Load structures

structures = ase.io.read("data/explore_c-gap-20u.xyz", ":")

# %%
#
# Compute the RDF for each structure using ASE

rmax = 3.5
nbins = 50
rdf_values = []
for structure in structures:
    rdf = get_rdf(structure, rmax=rmax, nbins=nbins, no_dists=True)
    rdf_values.append(rdf)

# The r grid corresponding to the RDF bins
dr = rmax / nbins
r_grid = np.linspace(dr / 2, rmax - dr / 2, nbins)

# %%
#
# Create a chemiscope widget with the RDF as a parameterized property

chemiscope.show(
    structures=structures,
    properties={
        "RDF": {
            "target": "structure",
            "values": np.array(rdf_values),
            "parameters": ["r"],
            "description": "Radial distribution function g(r)",
        },
        "config_type": [s.info["config_type"] for s in structures],
    },
    parameters={
        "r": {
            "values": r_grid,
            "units": "Angstrom",
        },
    },
    mode="structure",
    metadata=dict(
        name="RDF of selected C structures",
        description="Radial distribution functions computed for carbon structures "
        "from the C-GAP-20U dataset, displayed as parameterized properties.",
    ),
    settings=chemiscope.quick_settings(periodic=True),
)
