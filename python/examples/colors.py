"""
Atom property coloring
======================

This example demonstrates how to color atoms based on scalar properties.

Note that the same parameters can be used with `chemiscope.show` to visualize an
interactive widget in a Jupyter notebook.
"""

import ase.io
import numpy as np

import chemiscope

# loads a dataset of structures
frames = ase.io.read("data/alpha-mu.xyz", ":")

# compute some scalar quantities to display as atom coloring
polarizability = []
alpha_eigenvalues = []
anisotropy = []

for frame in frames:
    # center in the box
    frame.positions += frame.cell.diagonal()*0.5
    for axx, ayy, azz, axy, axz, ayz in zip(
        frame.arrays["axx"],
        frame.arrays["ayy"],
        frame.arrays["azz"],
        frame.arrays["axy"],
        frame.arrays["axz"],
        frame.arrays["ayz"],
    ):
        polarizability.append((axx + ayy + azz) / 3)

        # one possible measure of anisotropy...
        eigenvalues = np.linalg.eigvalsh(
            [[axx, axy, axz], [axy, ayy, ayz], [axz, ayz, azz]]
        )
        alpha_eigenvalues.append(eigenvalues)

        anisotropy.append(eigenvalues[2] - eigenvalues[0])


# now we just write the chemiscope input file
chemiscope.write_input(
    "colors-example.json.gz",
    frames=frames,
    # properties can be extracted from the ASE.Atoms frames
    properties={
        "polarizability": np.vstack(polarizability),
        "anisotropy": np.vstack(anisotropy),
        "alpha_eigenvalues": np.vstack(alpha_eigenvalues),
    },
    # the write_input function also allows defining the default visualization settings
    settings={
        "map": {
            "x": {"property": "alpha_eigenvalues[1]"},
            "y": {"property": "alpha_eigenvalues[2]"},
            "z": {"property": "alpha_eigenvalues[3]"},
            "color": {"property": "anisotropy"},
        },
        "structure": [
            {
                "color": {"property": "anisotropy"},
            }
        ],
    },
    # the properties we want to visualise are atomic properties - in order to view them
    # in the map panel, we must indicate that all atoms are environment centers
    environments=chemiscope.all_atomic_environments(frames),
)
