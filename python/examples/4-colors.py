"""
Atom property coloring
======================

This example demonstrates how to color atoms based on scalar properties.

Note that the same parameters can be used with :py:func:`chemiscope.show` to visualize
an interactive widget in a Jupyter notebook.
"""

# %%
#

import ase.io
import numpy as np

import chemiscope


# %%
#
# Load a dataset of structures containing polarizability and dipole data

frames = ase.io.read("data/alpha-mu.xyz", ":")

# %%
#
# Compute some scalar quantities to display as atom coloring

polarizability = []
alpha_eigenvalues = []
anisotropy = []

for frame in frames:
    # center molecule in the box
    frame.positions += frame.cell.diagonal() * 0.5
    for axx, ayy, azz, axy, axz, ayz in zip(
        frame.arrays["axx"],
        frame.arrays["ayy"],
        frame.arrays["azz"],
        frame.arrays["axy"],
        frame.arrays["axz"],
        frame.arrays["ayz"],
        strict=True,
    ):
        polarizability.append((axx + ayy + azz) / 3)

        # one possible measure of anisotropy...
        eigenvalues = np.linalg.eigvalsh(
            [[axx, axy, axz], [axy, ayy, ayz], [axz, ayz, azz]]
        )
        alpha_eigenvalues.append(eigenvalues)

        anisotropy.append(eigenvalues[2] - eigenvalues[0])

# %%
#
# Create a visualization and save it as a file that can be viewed at
# `<chemiscope.org>`_:

chemiscope.write_input(
    "colors-example.json.gz",
    frames=frames,
    # properties can also be extracted from the ASE.Atoms frames
    properties={
        "polarizability": np.vstack(polarizability),
        "anisotropy": np.vstack(anisotropy),
        "alpha_eigenvalues": np.vstack(alpha_eigenvalues),
    },
    # it is also possible to define the default visualization settings, e.g. map axes,
    # color property and palette, optionally
    settings={
        "map": {
            "x": {"property": "alpha_eigenvalues[1]"},
            "y": {"property": "alpha_eigenvalues[2]"},
            "z": {"property": "alpha_eigenvalues[3]"},
            "color": {"property": "anisotropy", "palette": "inferno"},
        },
        "structure": [
            {
                "color": {"property": "anisotropy", "palette": "bwr"},
            }
        ],
    },
    # the properties we want to visualise are atomic properties - in order to view them
    # in map panel we must indicate the list of environments (all atoms in this case)
    environments=chemiscope.all_atomic_environments(frames),
)


# %%
#
# The file can also be viewed in a notebook. Use `chemiscope.show` above to bypass the
# creation of a JSON file and directly create a viewer.

chemiscope.show_input("colors-example.json.gz")
