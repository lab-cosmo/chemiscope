"""
Extended shape visualization in chemiscope
==========================================

This example demonstrates how to visualize structure and atomic 
properties in the structure panel, using different types of 
visualizations (ellipsoids for tensors, arrows for vectors).
The example also shows how to define custom shapes. 

Note that the same parameters can be used with `chemiscope.show`
to visualize an interactive widget in a Jupyter notebook.
"""

import ase.io as aseio
import chemiscope
import numpy as np

# loads a dataset of structures
frames = aseio.read("data/alpha-mu.xyz", ":")

# converts the arrays from the format they are stored in to an array
# format that can be processed by the ASE utilities
for a in frames:
    a.arrays["alpha"] = np.array(
        [
            [axx, ayy, azz, axy, axz, ayz]
            for (axx, ayy, azz, axy, axz, ayz) in zip(
                a.arrays["axx"],
                a.arrays["ayy"],
                a.arrays["azz"],
                a.arrays["axy"],
                a.arrays["axz"],
                a.arrays["ayz"],
            )
        ]
    )

# here we define shapes that will later be used to create the input;
# input generation can also be achieved as a single call, but in practice
# it is wise to define separate entities for better readability

# cubes with smooth shading, centered on atoms. these are created as
# "custom" shapes and then are just scaled to atom-dependent sizes
atom_sizes = {"O": 0.4, "N": 0.5, "C": 0.45, "H": 0.2}
smooth_cubes = dict(
    kind="custom",
    parameters={
        "global": {
            "vertices": [  # this list defines the vertices
                [-1, -1, -1],
                [1, -1, -1],
                [-1, 1, -1],
                [1, 1, -1],
                [-1, -1, 1],
                [1, -1, 1],
                [-1, 1, 1],
                [1, 1, 1],
            ],
            "simplices": [  # and these are the indices of the vertices that form the triangular mesh
                [0, 2, 1],
                [1, 2, 3],
                [4, 6, 5],
                [5, 6, 7],
                [0, 1, 4],
                [0, 4, 2],
                [1, 3, 5],
                [2, 6, 3],
                [1, 5, 4],
                [2, 4, 6],
                [3, 7, 5],
                [3, 6, 7],
            ],
        },
        # the cube is used at each atomic position, only difference being the scaling
        "atom": [
            {"scale": atom_sizes[label]}
            for label in (list(frames[0].symbols) + list(frames[1].symbols))
        ],
    },
)

# structure-based shapes. also demonstrates how to achieve sharp-edge shading.
# requires defining multiple times the same vertices
sharp_cubes = dict(
    kind="custom",
    parameters={
        "global": {
            "vertices": [  # in order to get "sharp" edges, you need to define separate vertices for each facet
                [0, 0, 0],
                [0, 1, 0],
                [1, 1, 0],
                [1, 0, 0],
                [0, 0, 0],
                [0, 0, 1],
                [0, 1, 1],
                [0, 1, 0],
                [0, 1, 0],
                [0, 1, 1],
                [1, 1, 1],
                [1, 1, 0],
                [1, 1, 0],
                [1, 1, 1],
                [1, 0, 1],
                [1, 0, 0],
                [1, 0, 0],
                [1, 0, 1],
                [0, 0, 1],
                [0, 0, 0],
                [0, 0, 1],
                [1, 0, 1],
                [1, 1, 1],
                [0, 1, 1],
            ],
            "simplices": [  # simplices defining the mesh - two triangles per facet
                [0, 1, 2],
                [2, 3, 0],
                [4, 5, 6],
                [6, 7, 4],
                [8, 9, 10],
                [10, 11, 8],
                [12, 13, 14],
                [14, 15, 12],
                [16, 17, 18],
                [18, 19, 16],
                [20, 21, 22],
                [22, 23, 20],
            ],
        },
        "structure": [  # structure positioning is relative to the origin of the axes
            {"position": [0, 4, 0], "color": 0xFF0000},
            {
                "position": [3, 2, 1],
                "color": 0x00FF00,
                "orientation": [0.2, 0.4, 0.1, 1],
            },
        ],
    },
)
# dipole moments visualized as arrows. this is just to demonstrate manual insertion,
# see below to extract directly from the ASE info
dipoles_manual = (
    dict(
        kind="arrow",
        parameters={
            "global": {
                "base_radius": 0.2,
                "head_radius": 0.3,
                "head_length": 0.5,
                "color": 0xFF00B0,
            },
            "structure": [
                {"vector": frames[0].info["dipole_ccsd"].tolist()},
                {"vector": frames[0].info["dipole_ccsd"].tolist()},
            ],
        },
    ),
)

dipoles_auto = chemiscope.extract_vectors_from_ase(frames, "dipole_ccsd", scale=0.5)
# one can always update the defaults created by these automatic functions
dipoles_auto["parameters"]["global"] = {
    "base_radius": 0.2,
    "head_radius": 0.3,
    "head_length": 0.5,
    "color": 0xFF00B0,
}

chemiscope.write_input(
    "alphamu-shapes.json.gz",
    frames=frames,
    properties=chemiscope.extract_properties(frames, only=["alpha"]),
    shapes={
        # cubes with smooth shading, centered on atoms
        "smooth_cubes": smooth_cubes,
        # demonstrates showing a "global" shape for each structure
        "cube": sharp_cubes,
        # (molecular) electric dipole
        "dipole": dipoles_auto,
        # atomic decomposition of the polarizability as ellipsoids. use utility to extract from the ASE frames
        "alpha": chemiscope.extract_tensors_from_ase(
            frames, "alpha", force_positive=True, scale=0.2
        ),
    },
    settings={  # the write_input function also allows defining the default visualization settings
        "map": {
            "x": {"property": "alpha[1]"},
            "y": {"property": "alpha[2]"},
            "z": {"property": "alpha[3]"},
            "palette": "seismic",
            "color": {"property": ""},
        },
        "structure": [
            {
                "spaceFilling": False,
                "atomLabels": False,
                "atoms": False,
                "shape": "alpha,dipole",  # multiple shapes can be visualized at the same time!
                "axes": "off",
                "keepOrientation": False,
                "playbackDelay": 700,
                "environments": {
                    "activated": True,
                    "bgColor": "CPK",
                    "bgStyle": "licorice",
                    "center": False,
                    "cutoff": 0.5,
                },
            }
        ],
    },
    environments=chemiscope.all_atomic_environments(frames),
)
