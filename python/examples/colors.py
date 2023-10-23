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
# format, and computes scalar quantities to display as atom coloring
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


# now we just write the chemiscope datafile
chemiscope.write_input(
    "colors-example.json.gz",
    frames=frames,
    # properties can be extracted from the ASE.Atoms frames
    properties=chemiscope.extract_properties(
        frames, only=["polarizability", "anisotropy", "alpha_eigenvalues"]
    ),
    # the write_input function also allows defining the default visualization settings
    settings={  
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
                # these are the color visualization options
                "color": {"property": "anisotropy", "min": 1, "max": 15},
            }
        ],
    },
    # these are atomic properties - in order to view them in the map panel, we must 
    # indicate that all atoms are environment centers
    environments=chemiscope.all_atomic_environments(frames),
)
