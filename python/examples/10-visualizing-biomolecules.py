"""
Visualize biomolecules with MDAnalysis, chemfiles or OpenMM
===========================================================

This example shows how to visualize biomolecules in chemiscope with `MDAnalysis
<https://www.mdanalysis.org/>`_, `chemfiles <https://chemfiles.org/>`_ or `OpenMM
<https://openmm.org/>`_, and how to leverage the selection capabilities of these
libraries to show only a subset of the atoms.

Biomolecules often contain a large number of atoms, which makes the classical
ball-and-stick representation of molecules hard to be understand.
The "cartoon" representation focuses on the main structural elements, such as
backbone atoms, to highlight the secondary structure, and is often more readable.
"""

import urllib.request

import chemfiles
import MDAnalysis as mda
import numpy as np
import openmm.app
import openmm.unit
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

# Just like with MDAnalysis, we can pass the chemfiles frame directly to chemiscope.
chemiscope.show(
    structures=trajectory,
    mode="structure",
    settings=chemiscope.quick_settings(structure_settings={"cartoon": True}),
)

# %%
#
# Selecting atoms of interest
# ---------------------------
#
# We can also use `chemfiles selections <chfl-select>`_ to filter atoms. The selection
# returns the indices of the atoms matching the query. We then remove the atoms we
# don't want from the frame.
#
# .. _chfl-select: https://chemfiles.org/chemfiles/latest/selections.html

# takes the first frame from the trajectory (in case there are many...)
frame = trajectory.read()
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
# Using OpenMM
# ++++++++++++
#
# `OpenMM <https://openmm.org/>`_ keeps a structure as a topology plus a set of
# positions, and chemiscope accepts anything carrying those two attributes:
# ``openmm.app.PDBFile``, ``openmm.app.PDBxFile``, ``openmm.app.Modeller``, and
# whatever else is built on top of them.

pdb = openmm.app.PDBFile(f"./{pdb_id}.pdb")

chemiscope.show(
    structures=pdb,
    mode="structure",
    settings=chemiscope.quick_settings(structure_settings={"cartoon": True}),
)

# %%
#
# Selecting atoms of interest
# ---------------------------
#
# OpenMM edits a structure through ``openmm.app.Modeller``, so dropping the
# crystallographic water means deleting those residues. The modeller itself is a
# valid chemiscope structure, so there is nothing to convert afterwards.

modeller = openmm.app.Modeller(pdb.topology, pdb.positions)
modeller.delete(
    [residue for residue in modeller.topology.residues() if residue.name == "HOH"]
)

chemiscope.show(
    structures=modeller,
    mode="structure",
    settings=chemiscope.quick_settings(structure_settings={"cartoon": True}),
)

# %%
#
# Showing a binding site as whole residues
# ----------------------------------------
#
# 10MH is a methyltransferase caught with its cofactor, S-adenosyl-homocysteine
# (residue ``SAH``). To look at the pocket around it, we put an environment on every
# cofactor atom and let chemiscope draw the surroundings within a cutoff.
#
# A spherical cutoff slices through whatever it meets, which leaves half-drawn side
# chains at the edge of the pocket. The ``byResidue`` environment setting keeps every
# residue that has at least one atom inside the cutoff whole, so the pocket is made of
# complete amino acids. The switch is in the hamburger menu, under "environments", and
# only appears for structures that carry residue information.

topology = modeller.topology
cofactor = [atom.index for atom in topology.atoms() if atom.residue.name == "SAH"]

# one environment per cofactor atom, each with a 5 A cutoff
environments = [(0, atom_index, 5.0) for atom_index in cofactor]

# environments need a property for each of them: here, the distance of each cofactor
# atom to the centre of the cofactor
positions = np.asarray(
    modeller.positions.value_in_unit(openmm.unit.angstrom), dtype=float
)
cofactor_positions = positions[cofactor]
distances = np.linalg.norm(cofactor_positions - cofactor_positions.mean(axis=0), axis=1)

chemiscope.show(
    structures=modeller,
    mode="structure",
    properties={
        "distance to cofactor centre": {
            "target": "atom",
            "values": distances,
            "units": "A",
        }
    },
    environments=environments,
    settings=chemiscope.quick_settings(
        structure_settings={
            "cartoon": True,
            "environments": {"cutoff": 5.0, "byResidue": True, "bgStyle": "cartoon"},
        }
    ),
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
# Note: write_external_structures works with both MDAnalysis AtomGroups
# and chemfiles structures.
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


# %%
#
# Nucleic acids and glycans
# +++++++++++++++++++++++++
#
# The cartoon representation is not limited to protein secondary structure. The DNA
# bound to the methyltransferase above is drawn as a backbone ribbon with the
# nucleobases as filled rings, coloured by base (A red, C yellow, G green, T blue,
# U cyan). Nucleotide names are recognised in the PDB (``DA``, ``A``), Amber
# (``DA5``, ``RA3``) and CHARMM (``ADE``) conventions.
#
# Showing the DNA on its own makes the rings easier to see than in the full complex.

dna = universe.select_atoms("nucleic")
chemiscope.show(
    structures=dna,
    mode="structure",
    settings=chemiscope.quick_settings(structure_settings={"cartoon": True}),
)

# %%
#
# Carbohydrates are drawn with the `Symbol Nomenclature for Glycans
# <https://www.ncbi.nlm.nih.gov/glycans/snfg.html>`_ instead: each monosaccharide
# becomes its SNFG symbol (blue cube for GlcNAc, green sphere for mannose, yellow
# sphere for galactose, red cone for fucose, ...) and the symbols are linked by
# sticks following the glycosidic bonds. This is on by default, under ``snfg``.
#
# 1IGY is an intact antibody, which carries the N-glycan of the Fc region, so it
# shows both representations at once: protein cartoon plus the glycan symbols.

pdb_id = "1IGY"
urllib.request.urlretrieve(
    f"https://files.rcsb.org/view/{pdb_id}.pdb", f"./{pdb_id}.pdb"
)

antibody = mda.Universe(f"./{pdb_id}.pdb")
chemiscope.show(
    structures=antibody.select_atoms("not water"),
    mode="structure",
    settings=chemiscope.quick_settings(
        structure_settings={"cartoon": True, "snfg": True}
    ),
)

# %%
#
# The glycan alone, with the protein removed, is a small tree of four to eight
# symbols; this is the view used to compare glycoforms.

glycan = antibody.select_atoms("resname NAG NDG BMA MAN GAL FUC")
chemiscope.show(
    structures=glycan,
    mode="structure",
    settings=chemiscope.quick_settings(structure_settings={"snfg": True}),
)
