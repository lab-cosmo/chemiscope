"""
Atom property coloring in chemiscope
====================================

This example demonstrates how to color atoms based on scalar
properties.

Note that the same parameters can be used with `chemiscope.show`
to visualize an interactive widget in a Jupyter notebook.
"""

import ase.io
import chemiscope
import numpy as np

# loads a dataset of structures
frames = ase.io.read("data/alpha-mu.xyz", ":")

# converts the arrays from the format they are stored in to an array
# format that can be processed by the ASE utilities
for a in frames:
    a.arrays["polarizability"] = np.array(
        [
            (axx + ayy + azz) / 3
            for (axx, ayy, azz) in zip(
                a.arrays["axx"], a.arrays["ayy"], a.arrays["azz"]
            )
        ]
    )

    # one possible measure of anisotropy...
    a.arrays["alpha_eigenvalues"] = np.array(
        [
            np.linalg.eigvalsh([[axx, axy, axz], [axy, ayy, ayz], [axz, ayz, azz]])
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

    a.arrays["anisotropy"] = (
        a.arrays["alpha_eigenvalues"][:, 2] - a.arrays["alpha_eigenvalues"][:, 0]
    )

chemiscope.write_input(
    "colors-example.json.gz",
    frames=frames,
    properties=chemiscope.extract_properties(
        frames, only=["polarizability", "anisotropy", "alpha_eigenvalues"]
    ),
    settings={  # the write_input function also allows defining the default visualization settings
        "map": {
            "x": {"property": "alpha_eigenvalues[1]"},
            "y": {"property": "alpha_eigenvalues[2]"},
            "z": {"property": "alpha_eigenvalues[3]"},
            "palette": "seismic",
            "color": {"property": "anisotropy"},
        },
        "structure": [
            {
                "spaceFilling": False,
                "atomLabels": False,
                "atoms": True,
                "axes": "off",
                "keepOrientation": False,
                "playbackDelay": 700,
                "environments": {
                    "activated": True,
                    "bgColor": "CPK",
                    "bgStyle": "licorice",
                    "center": False,
                    "cutoff": 4,
                },
                "color": {"property": "anisotropy", "min": 1, "max": 15},
            }
        ],
    },
    environments=chemiscope.all_atomic_environments(frames),
)
