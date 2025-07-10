"""
Structure-property maps
=======================

This example demonstrates the visualization of structures (or environments) using
data-driven descriptors of their geometry to cluster together similar motifs. Here, the
geometric descriptors have been computed by PCA starting from SOAP representations but
are provided as text files to avoid external dependencies for the example.

The same parameters demonstrated on this example can be used with ``chemiscope.show`` to
visualize an interactive widget in a Jupyter notebook.
"""

# %%
#

import ase.io
import numpy as np

import chemiscope


# %%
# Load structures:

frames = ase.io.read("data/trajectory.xyz", ":")

# %%
#
# Load the SOAP-PCA descriptors. To obtain them yourself, check the :ref:`examples
# <explore-advanced-example>` of feature calculation and dimentionality reduction.
# Chemiscope also provides ``chemiscope.explore`` function to obtain low-dimentional
# representation based on PETMAD features.

pca_atom = np.loadtxt("data/trajectory-pca_atom.dat")
pca_structure = np.loadtxt("data/trajectory-pca_structure.dat")

# %%
#
# When both environments and structure properties are present, a toggle allows you to
# switch between both modes.
#
# .. note::
#
#     if there are properties stored in the ASE frames, you can extract them with
#     chemiscope.extract_properties(frames)

properties = {
    # concise definition of a property, with just an array and the type
    # inferred by the size
    "structure PCA": pca_structure,
    "atom PCA": pca_atom,
    # an example of the verbose definition
    "energy": {
        "target": "structure",
        "values": [frame.info["dftb_energy_eV"] for frame in frames],
        "units": "eV",
        "description": "potential energy, computed with DFTB+",
    },
}

# %%
#
# Environment descriptors have only been computed for C and O atoms.
environments = []
cutoff = 4.0
for frame_i, frame in enumerate(frames):
    for atom_i, atom in enumerate(frame.numbers):
        if atom == 6 or atom == 8:
            environments.append((frame_i, atom_i, cutoff))


# %%
#
# Create a visualization and save it as a file that can be viewed at `chemiscope.org
# <https://chemiscope.org>`_:

chemiscope.write_input(
    "trajectory-pca.json.gz",
    # dataset metadata can also be included to provide a self-contained description
    # of the data, authors, and references
    meta={
        "name": "Allyl alcohol PCA map",
        "description": (
            "This dataset contains a PCA map of the C and O environments "
            "from a few frames out of a MD simulation of allyl alcohol, C3H5OH."
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
    environments=environments,
    settings={  # these are reasonable settings for trajectory visualization
        "structure": [{"keepOrientation": True, "playbackDelay": 100}]
    },
)

# %%
#
# The file can also be viewed in a notebook. Use `chemiscope.show` above to bypass the
# creation of a JSON file and directly create a viewer.

chemiscope.show_input("trajectory-pca.json.gz")
