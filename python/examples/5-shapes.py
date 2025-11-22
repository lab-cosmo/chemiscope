"""
Shape visualization
===================

This example demonstrates how to visualize structure and atomic
properties in the structure panel, using different types of
predefined shapes (ellipsoids for tensors, arrows for vectors).
The example also shows how to define custom shapes.

Note that the same parameters can be used with :py:func:`chemiscope.show` to visualize
an interactive widget in a Jupyter notebook.
"""

# %%
#

import ase.io
import numpy as np
from scipy.spatial.transform import Rotation

import chemiscope


# %%
#
# Loads a dataset of structures

frames = ase.io.read("data/alpha-mu.xyz", ":")

quaternions = []
# converts the arrays from the format they are stored in to an array
# format that can be processed by the ASE utilities
for a in frames:
    a.positions += a.cell.diagonal() * 0.5
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
                strict=True,
            )
        ]
    )

    # interatomic separations (used to orient "stuff" later on)
    dists = a.get_all_distances(mic=True)
    np.fill_diagonal(dists, np.max(dists))
    for i in range(len(a)):
        nneigh = dists[i].argmin()
        vec = a.get_distance(i, nneigh, vector=True, mic=True)
        vec /= np.linalg.norm(vec)
        quat = Rotation.align_vectors([np.array([0, 0, 1])], [vec])[0].as_quat()
        quaternions.append(quat)

# %%
#
# Shapes generation
# -----------------
#
# Here we define shapes that will later be used to create the input.
# Input generation can also be achieved as a single call, but in practice
# it is wise to define separate entities for better readability.

# %%
#
# Cubes with smooth shading, centered on atoms. These are created as
# "custom" shapes and then are just scaled to atom-dependent sizes:

atom_sizes = {"O": 0.4, "N": 0.5, "C": 0.45, "H": 0.2}
smooth_cubes = dict(
    kind="custom",
    parameters={
        "global": {
            # this list defines the vertices
            "vertices": [
                [-1, -1, -1],
                [1, -1, -1],
                [-1, 1, -1],
                [1, 1, -1],
                [-1, -1, 1],
                [1, -1, 1],
                [-1, 1, 1],
                [1, 1, 1],
            ],
            # and these are the indices of the vertices that form the triangular mesh
            "simplices": [
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

# %%
#
# Structure-based shapes. Also demonstrates how to achieve sharp-edge shading.
# It requires defining multiple times the same vertices:

sharp_cubes = dict(
    kind="custom",
    parameters={
        "global": {
            # in order to get "sharp" edges, you need to define separate vertices for
            # each facet
            "vertices": [
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
            # simplices defining the mesh - two triangles per facet
            "simplices": [
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
        # structure positioning is relative to the origin of the axes
        "structure": [
            {"position": [12, 14, 12], "color": 0xFF0000},
            {
                "position": [15, 14, 12],
                "color": 0x00FF00,
                "orientation": [0.5, -0.5, 0, 1 / np.sqrt(2)],
            },
        ],
    },
)

# %%
#
# Load vertices and simplices from an external file
# and use it to draw a very irreverent molecule
irreverent_dict = np.load("data/irreverence.npz")
irreverent_shape = dict(
    kind="custom",
    parameters={
        "global": {
            "vertices": irreverent_dict["vertices"].tolist(),
            "simplices": irreverent_dict["simplices"].tolist(),
            "scale": 0.02,
        },
        "atom": [
            {"orientation": quaternions[i].tolist()}
            for i in range(len(frames[0]) + len(frames[1]))
        ],
    },
)

# %%
#
# Dipole moments visualized as arrows. This is just to demonstrate manual insertion,
# see below to extract directly from the ASE info
dipoles_manual = (
    dict(
        kind="arrow",
        parameters={
            "global": {
                "baseRadius": 0.2,
                "headRadius": 0.3,
                "headLength": 0.5,
                "color": 0xFF00B0,
            },
            "structure": [
                {
                    "position": [12, 12, 12],
                    "vector": frames[0].info["dipole_ccsd"].tolist(),
                },
                {
                    "position": [12, 12, 12],
                    "vector": frames[1].info["dipole_ccsd"].tolist(),
                },
            ],
        },
    ),
)


dipoles_auto = chemiscope.ase_vectors_to_arrows(frames, "dipole_ccsd", scale=0.5)
# one can always update the defaults created by these automatic functions
dipoles_auto["parameters"]["global"] = {
    "baseRadius": 0.2,
    "headRadius": 0.3,
    "headLength": 0.5,
    "color": 0xFF00B0,
}
for d in dipoles_auto["parameters"][
    "structure"
]:  # center the dipoles close to the molecule
    d["position"] = [11, 11, 11]

# %%
#
# Create a visualization and save it as a file that can be viewed at
# `<chemiscope.org>`_:

chemiscope.write_input(
    "shapes-example.json.gz",
    frames=frames,
    properties=chemiscope.extract_properties(frames, only=["alpha"]),
    shapes={
        # cubes with smooth shading, centered on atoms
        "smooth_cubes": smooth_cubes,
        # demonstrates showing a "global" shape for each structure
        "cube": sharp_cubes,
        # (molecular) electric dipole
        "dipole": dipoles_auto,
        # atomic decomposition of the polarizability as ellipsoids. use utility to
        # extract from the ASE frames
        "alpha": chemiscope.ase_tensors_to_ellipsoids(
            frames, "alpha", force_positive=True, scale=0.2
        ),
        # shapes with a bit of flair
        "irreverence": irreverent_shape,
    },
    # the write_input function also allows defining the default visualization settings
    settings={
        "map": {
            "x": {"property": "alpha[1]"},
            "y": {"property": "alpha[2]"},
            "z": {"property": "alpha[3]"},
            "color": {"property": "", "palette": "seismic"},
        },
        "structure": [
            {
                "spaceFilling": False,
                "atomLabels": False,
                "atoms": False,
                # multiple shapes can be visualized at the same time!
                "shape": "alpha,dipole",
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


# %%
#
# The file can also be viewed in a notebook. Use `chemiscope.show` above to bypass the
# creation of a JSON file and directly create a viewer.

chemiscope.show_input("shapes-example.json.gz", mode="structure")
