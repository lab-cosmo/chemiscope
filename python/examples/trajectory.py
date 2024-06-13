"""
Trajectory plotting
===================

This example demonstrates the visualization of trajectory data. It also includes
visualization of atomic forces.

The same parameters can be used with `chemiscope.show` to visualize an interactive
widget in a Jupyter notebook.
"""

# %%
#

import ase.io
import numpy as np

import chemiscope

# %%
#
# Load structures and properties
frames = ase.io.read("data/trajectory.xyz", ":")

properties = {
    # concise definition of a property, with just an array and the type
    # inferred by the size
    "index": np.arange(len(frames)),
    # an example of the verbose definition
    "energy": {
        "target": "structure",
        "values": [frame.info["energy"] for frame in frames],
        "units": "eV",
        "description": "potential energy, computed with DFTB+",
    },
}


# %%
#
# Uses a chemiscope widget to visualize structures and properties.
# Note how to display the forces as vectorial properties

cs = chemiscope.show(
    # dataset metadata can also be included, to provide a self-contained description
    # of the data, authors and references
    meta={
        "name": "Allyl alcohol MD trajectory.",
        "description": (
            "This dataset contains data from a DFTB+ trajectory of allyl alcohol."
        ),
        "authors": ["The chemiscope developers"],
        "references": [
            (
                "G. Fraux, R. Cersonsky, and M. Ceriotti, "
                '"Chemiscope: interactive structure-property explorer for materials '
                'and molecules," JOSS 5(51), 2117 (2020).'
            )
        ],
    },
    frames=frames,
    properties=properties,
    # visualize forces as vectors
    shapes={
        "forces": chemiscope.ase_vectors_to_arrows(
            frames, "forces", scale=1, radius=0.15
        )
    },
    settings={  # these are reasonable settings for trajectory visualization
        "structure": [
            {
                "keepOrientation": True,
                "playbackDelay": 100,
                "shape": "forces",  # visualize force vectors
            }
        ]
    },
)
