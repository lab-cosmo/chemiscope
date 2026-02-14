"""
Visualize biomolecules with MDAnalysis or chemfiles
===================================================

This example shows how to visualize biomolecules in chemiscope with `MDAnalysis
<https://www.mdanalysis.org/>`_ or `chemfiles <https://chemfiles.org/>`_, and
how  to leverage the selection capabilities of these libraries to show only a
subset of the atoms.

Biomolecules often contain a large number of atoms, which makes the classical
ball-and-stick representation of molecules hard to be understand.
The "cartoon" representation focuses on the main structural elements, such as
backbone atoms, to highlight the secondary structure, and is often more readable.
"""

import urllib.request

import chemfiles
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.analysis.rms import RMSD
from MDAnalysis.tests.datafiles import GRO_MEMPROT, XTC_MEMPROT

import chemiscope


# %%
#
# Retrieving the PDB file from RCSB Protein Data Bank
# +++++++++++++++++++++++++++++++++++++++++++++++++++
#
# The `RCSB Protein Data Bank <https://www.rcsb.org/>`_ (RCSB PDB) is a database of
# crystal structures of proteins, nucleic acids and small molecules. To start with, we
# will retrieve a structure from the PDB database. Here we choose "10MH", a complex
# consisting of a protein, a nucleic acid, small molecules, and crystallographic water.

pdb_id = "10MH"
urllib.request.urlretrieve(
    f"https://files.rcsb.org/view/{pdb_id}.pdb", f"./{pdb_id}.pdb"
)

# %%
#
# Using MDAnalysis
# ++++++++++++++++
#
# We use `MDAnalysis <https://www.mdanalysis.org/>`_ to read the PDB file, interpreting
# also the metadata that describes the structure of the protein.

universe = mda.Universe(f"./{pdb_id}.pdb")

# The `chemiscope` takes a `MDAnalysis.AtomGroup
# <https://userguide.mdanalysis.org/stable/atomgroup.html>`_ as input. You can toggle
# the cartoon representation in the hamburger menu in the top-right corner of the
# widget. When the cartoon representation is off, the representation will automatically
# fall back to the ball-and-stick representation.

ag = universe.atoms
chemiscope.show(
    structures=ag,
    mode="structure",
    settings=chemiscope.quick_settings(structure_settings={"cartoon": True}),
)

# %%
#
# Selecting atoms of interest
# ---------------------------
#
# The crystallographic water in the structure is not of interest, so we can use the
# `select_atoms()` method to only show the complex for a cleaner visualization.

sol = universe.select_atoms("not water")
chemiscope.show(
    structures=sol,
    mode="structure",
    settings=chemiscope.quick_settings(structure_settings={"cartoon": True}),
)

# %%
#
# Using chemfiles
# +++++++++++++++
#
# Alternatively, we can use `chemfiles <https://chemfiles.org/>`_ to read the PDB file.

trajectory = chemfiles.Trajectory(f"./{pdb_id}.pdb")
frame = trajectory.read()

# Just like with MDAnalysis, we can pass the chemfiles frame directly to chemiscope.
chemiscope.show(
    structures=frame,
    mode="structure",
    settings=chemiscope.quick_settings(structure_settings={"cartoon": True}),
)

# %%
#
# Selecting atoms of interest
# ---------------------------
#
# We can also use chemfiles selections to filter atoms. The selection returns
# the indices of the atoms matching the query. We then remove the atoms we
# don't want from the frame.

selection = chemfiles.Selection("resname HOH")
to_remove = selection.evaluate(frame)

for i in sorted(to_remove, reverse=True):
    frame.remove(i)

chemiscope.show(
    structures=frame,
    mode="structure",
    settings=chemiscope.quick_settings(structure_settings={"cartoon": True}),
)


# %%
#
# Exploring the sampled conformational space
# ++++++++++++++++++++++++++++++++++++++++++
#
# We can use the map mode to explore the conformational space sampled by the MD
# simulation easily. Here, we use a protein-lipid system taken from the
# `MDAnalysisTests <https://userguide.mdanalysis.org/stable/testing.html>`_ as an
# example.

complx = mda.Universe(GRO_MEMPROT, XTC_MEMPROT)

# We describe the conformational space by two features: the z-axis distance between the
# geometric centers of protein and lipid, and the root mean square deviation (RMSD) of
# the atomic positions of protein with respect to its initial conformation.

# Distance calculation
distances = []
for _ in complx.trajectory:
    lipid_center = complx.select_atoms("resname POP*").center_of_geometry()
    protein_center = complx.select_atoms("protein").center_of_geometry()
    distances.append((protein_center - lipid_center)[2])
distances = np.abs(distances)

# RMSD calculation
ref = mda.Universe(GRO_MEMPROT)
R = RMSD(complx, ref, select="backbone")
R.run()
rmsd = R.results.rmsd.T[2]

# We can then use the map mode to visualize the sampled conformational space.

# Given that trajectories can be very large, we load the structures on disk to
# reduce the memory usage of the viewer.
# Note: write_external_structures works with both MDAnalysis AtomGroups and chemfiles structures.
external_structures = chemiscope.write_external_structures(
    complx.atoms, prefix="protein-rmsd"
)

chemiscope.show(
    structures=external_structures,
    metadata={
        "name": "Protein-Lipid Complex (MDAnalysis)",
        "description": (
            "Conformational space of a protein-lipid complex featurized "
            "by the protein-lipid z-axis distance and the protein RMSD"
        ),
    },
    properties={
        "Protein-Lipid Distance": {
            "target": "structure",
            "values": distances,
            "units": "Å",
            "description": (
                "Z-axis distance between the geometric centers of protein and lipid"
            ),
        },
        "Protein Backbone RMSD": {
            "target": "structure",
            "values": rmsd,
            "units": "Å",
            "description": (
                "RMSD of the atomic positions of protein with respect to "
                "its initial conformation"
            ),
        },
    },
)

# %%
#
# We can achieve the same visualization using chemfiles to load the trajectory.
# We will use the properties calculated above for consistency.

# We use the paths from MDAnalysisTests for this example
trajectory = chemfiles.Trajectory(XTC_MEMPROT)
trajectory.set_topology(GRO_MEMPROT)

# chemiscope.write_external_structures accepts a chemfiles Trajectory directly
external_structures_chemfiles = chemiscope.write_external_structures(
    trajectory, prefix="protein-rmsd-chemfiles"
)

chemiscope.show(
    structures=external_structures_chemfiles,
    metadata={
        "name": "Protein-Lipid Complex (chemfiles)",
        "description": (
            "Conformational space of a protein-lipid complex loaded with chemfiles"
        ),
    },
    properties={
        "Protein-Lipid Distance": {
            "target": "structure",
            "values": distances,
            "units": "Å",
            "description": (
                "Z-axis distance between the geometric centers of protein and lipid"
            ),
        },
        "Protein Backbone RMSD": {
            "target": "structure",
            "values": rmsd,
            "units": "Å",
            "description": (
                "RMSD of the atomic positions of protein with respect to "
                "its initial conformation"
            ),
        },
    },
)
