"""
Visualize biomolecules with MDAnalysis
======================================

This example shows how to visualize biomolecules in chemiscope with `MDAnalysis
<https://www.mdanalysis.org/>`_, and how to leverage the `select_atoms()
<https://userguide.mdanalysis.org/stable/selections.html>`_ method to only show the
atoms of interest.

Biomolecules are often consisting of a large number of atoms, which makes the classical
ball-and-stick repersentation of molecules hard to be understood. Cartoon representation
focuses on the main structural elements, such as backbone atoms, to highlight the
secondary structure, and is often more readable.
"""

import urllib.request

import MDAnalysis as mda

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
# Reading the PDB file and visualizing it in Chemiscope
# +++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# We can use `MDAnalysis <https://www.mdanalysis.org/>`_ to read the PDB file.

universe = mda.Universe(f"./{pdb_id}.pdb")

# %%
#
# The `chemiscope` takes a `MDAnalysis.AtomGroup
# <https://userguide.mdanalysis.org/stable/atomgroup.html>`_ as input. You can toggle
# the cartoon representation in the hamburger menu in the top-right corner of the
# widget. When the cartoon representation is off, the representation will automatically
# fall back to the ball-and-stick.

ag = universe.atoms
chemiscope.show(
    frames=ag,
    mode="structure",
)

# %%
#
# Selecting atoms of interest
# +++++++++++++++++++++++++++
#
# The crystallographic water in the structure is not of interest, so we can use the
# `select_atoms()` method to only show the complex for a cleaner visualization.

sol = universe.select_atoms("not water")
chemiscope.show(
    frames=sol,
    mode="structure",
)
