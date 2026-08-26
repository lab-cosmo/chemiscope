"""
Atom property coloring and labels
=================================

This example demonstrates how to color atoms based on scalar properties, and how to
visualize property values as labels for the atoms.

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

structures = ase.io.read("data/alpha-mu.xyz", ":")

# %%
#
# Compute some scalar quantities to display as atom coloring

polarizability = []
alpha_eigenvalues = []
anisotropy = []

for structure in structures:
    # center molecule in the box
    structure.positions += structure.cell.diagonal() * 0.5
    for axx, ayy, azz, axy, axz, ayz in zip(
        structure.arrays["axx"],
        structure.arrays["ayy"],
        structure.arrays["azz"],
        structure.arrays["axy"],
        structure.arrays["axz"],
        structure.arrays["ayz"],
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
    structures=structures,
    # properties can also be extracted from the ASE.Atoms structures
    properties={
        "polarizability": np.vstack(polarizability),
        "anisotropy": np.vstack(anisotropy),
        "alpha_eigenvalues": np.vstack(alpha_eigenvalues),
    },
    # it is also possible to define the default visualization settings, e.g. map axes,
    # color property and palette, and to indicate that we want to show atom labels
    # with the anisotropy value
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
                "atomLabels": True,
                "labelsProperty": "anisotropy",
            }
        ],
    },
    # the properties we want to visualise are atomic properties - in order to view them
    # in map panel we must indicate the list of environments (all atoms in this case)
    environments=chemiscope.all_atomic_environments(structures),
)


# %%
#
# The file can also be viewed in a notebook. Use `chemiscope.show` above to bypass the
# creation of a JSON file and directly create a viewer.

chemiscope.show_input("colors-example.json.gz")
