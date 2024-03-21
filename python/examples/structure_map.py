"""
Structure-property maps in chemiscope
=====================================

This example demonstrates the visualization of structures 
(or environments) using data-driven descriptors of their 
geometry, to cluster together similar motifs. 
Here the geometric descriptors have been computed by 
PCA starting from SOAP representations, but are provided
as text files to avoid external dependencies for the 
example. 

The same parameters can be used with `chemiscope.show`
to visualize an interactive widget in a Jupyter notebook.
"""

import numpy as np
import ase.io
import chemiscope

# load structures
frames = ase.io.read("data/trajectory.xyz", ":")

# load the SOAP-PCA descriptors. chemiscope does not provide
# analysis routines, but you can look up for instance
# scikit-matter as a package to do dimensionality reduction
# analyses.
pca_atom = np.loadtxt("data/trajectory-pca_atom.dat")
pca_struc = np.loadtxt("data/trajectory-pca_structure.dat")

# when both environments and structure property are present
# only environment properties are shown. still they can be stored,
# and future versions of chemiscope may allow switching between
# the two modes.
# NB: if there are properties stored in the ASE frames, you can extract
#     them with chemiscope.extract_properties(frames)
properties = {
    # concise definition of a property, with just an array and the type
    # inferred by the size
    "structure PCA": pca_struc,
    "atom PCA": pca_atom,
    # an example of the verbose definition
    "energy": {
        "target": "structure",
        "values": [frame.info["energy"] for frame in frames],
        "units": "eV",
        "description": "potential energy, computed with DFTB+",
    },
}

# environment descriptors have only been computed for C and O atoms.
# we use a mask and then a utility function to generate the proper
# list of environments
for frame in frames:
    frame_mask = np.zeros(len(frame))
    frame_mask[np.where((frame.numbers == 6) | (frame.numbers == 8))[0]] = 1
    frame.arrays["center_atoms_mask"] = frame_mask

environments = chemiscope.librascal_atomic_environments(frames, cutoff=4.0)

# write the chemiscope input file.
chemiscope.write_input(
    path="trajectory-pca.json.gz",
    # dataset metadata can also be included, to provide a self-contained description
    # of the data, authors and references
    meta={
        "name": "Allyl alcohol PCA map",
        "description": "This dataset contains a PCA map of the C and O environments from a few frames out of a MD simulation of allyl alcohol, C3H5OH.",
        "authors": ["The chemiscope developers"],
        "references": [
            'G. Fraux, R. Cersonsky, and M. Ceriotti, "Chemiscope: interactive structure-property explorer for materials and molecules," JOSS 5(51), 2117 (2020).'
        ],
    },
    frames=frames,
    properties=properties,
    environments=environments,
    settings={  # these are reasonable settings for trajectory visualization
        "structure": [{"keepOrientation": True, "playbackDelay": 100}]
    },
)
