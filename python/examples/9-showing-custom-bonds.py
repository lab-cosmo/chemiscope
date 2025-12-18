"""
Using stk structures and showing custom bonds
=============================================

This example demonstrates interfacing `stk <https://stk.readthedocs.io/en/stable/>`_
with chemiscope and the use of shapes to show custom bonding in a molecule.

`stk <https://stk.readthedocs.io/en/stable/>`_ comes with topology/bonding
information by default (using the cheminformatic software rdkit).

We use `stko <https://stko-docs.readthedocs.io/en/latest/>`_ to calculate
some rudimentary properties of `stk` molecules. `stko` can be installed with
``pip install stko``.

"""

import itertools as it
import tempfile

import ase.io
import stk
import stko
from rdkit.Chem import AllChem as rdkit

import chemiscope


# %%
#
# Interfacing chemiscope with stk molecules
# ++++++++++++++++++++++++++++++++++++++++++
#
# Generate a list of stk BuildingBlocks (representation of a molecule) with
# properties. We start by constructing a cage and host-guest complex with stk.

cage = stk.ConstructedMolecule(
    topology_graph=stk.cage.FourPlusSix(
        building_blocks=(
            stk.BuildingBlock(
                smiles="NCCN",
                functional_groups=[stk.PrimaryAminoFactory()],
            ),
            stk.BuildingBlock(
                smiles="O=CC(C=O)C=O",
                functional_groups=[stk.AldehydeFactory()],
            ),
        ),
        optimizer=stk.MCHammer(),
    ),
)
host_guest = stk.ConstructedMolecule(
    topology_graph=stk.host_guest.Complex(
        host=stk.BuildingBlock.init_from_molecule(cage),
        guests=stk.host_guest.Guest(
            building_block=stk.BuildingBlock("[Br][Br]"),
        ),
    ),
)

# %%
#
# Including using stk to interface with rdkit molecules.

rdkitmol = rdkit.MolFromSmiles("Cc1ccccc1")
rdkitmol = rdkit.AddHs(rdkitmol)
rdkit.Kekulize(rdkitmol)
params = rdkit.ETKDGv3()
params.randomSeed = 0xF00D
rdkit.EmbedMolecule(rdkitmol, params)

# %%
#
# We can put this into a list of stk.Molecule objects.


structures = [
    # A mostly optimised cage molecule.
    cage,
    # A host guest molecule.
    host_guest,
    # From rdkit.
    stk.BuildingBlock.init_from_rdkit_mol(rdkitmol),
    # A metal-organic cage.
    stk.ConstructedMolecule(
        stk.cage.M2L4Lantern(
            building_blocks=(
                stk.BuildingBlock(
                    smiles="[Pd+2]",
                    functional_groups=(
                        stk.SingleAtom(stk.Pd(0, charge=2)) for i in range(4)
                    ),
                    position_matrix=[[0.0, 0.0, 0.0]],
                ),
                stk.BuildingBlock(
                    smiles=("C1=NC=CC(C2=CC=CC(C3=CC=NC=C3)=C2)=C1"),
                    functional_groups=[
                        stk.SmartsFunctionalGroupFactory(
                            smarts="[#6]~[#7X2]~[#6]",
                            bonders=(1,),
                            deleters=(),
                        ),
                    ],
                ),
            ),
            # Ensure that bonds between the
            # GenericFunctionalGroups of the ligand and the
            # SingleAtom functional groups of the metal are
            # dative.
            reaction_factory=stk.DativeReactionFactory(
                stk.GenericReactionFactory(
                    bond_orders={
                        frozenset(
                            {
                                stk.GenericFunctionalGroup,
                                stk.SingleAtom,
                            }
                        ): 9,
                    },
                ),
            ),
        ),
    ),
]


# %%
#
# Here we use stko (https://stko-docs.readthedocs.io/en/latest/) to compute
# properties of each molecule in the list and format it into the dictionary
# required by chemiscope.

energy = stko.UFFEnergy(ignore_inter_interactions=False)
shape_calc = stko.ShapeCalculator()
properties = {
    "uffenergy": [energy.get_energy(molecule) for molecule in structures],
    "aspheriticty": [
        shape_calc.get_results(molecule).get_asphericity() for molecule in structures
    ],
}

# %%
#
# A chemiscope widget showing the result with standard bonding derived from the
# bonding topology in stk.

chemiscope.show(
    structures=structures,
    properties=properties,
    settings=chemiscope.quick_settings(
        x="aspheriticty",
        y="uffenergy",
        color="",
        structure_settings={
            "atoms": True,
            "bonds": True,
            "spaceFilling": False,
        },
    ),
)

# %%
#
# Writing to a json.gz file, again without added bonding.

chemiscope.write_input(
    path="noshape_example.json.gz",
    structures=structures,
    properties=properties,
    meta=dict(name="Standard stk bonding."),
    settings=chemiscope.quick_settings(
        x="aspheriticty",
        y="uffenergy",
        color="",
        structure_settings={
            "atoms": True,
            "bonds": True,
            "spaceFilling": False,
        },
    ),
)

# %%
#
# For comparison, we show the same molecule imported with ASE versus stk, where
# chemiscope will automatically try to detect the bonds based on the geometry
# (because ASE does not contain this information), which will not be right due
# to the far-from equilibrium nature of the structure.

with tempfile.NamedTemporaryFile(suffix=".xyz") as tmpfile:
    structures[0].write(tmpfile.name)
    chemiscope.show(
        structures=[ase.io.read(tmpfile.name)],
        properties={i: [properties[i][0]] for i in properties},
        settings=chemiscope.quick_settings(
            x="aspheriticty",
            y="uffenergy",
            color="",
            structure_settings={
                "atoms": True,
                "bonds": True,
                "spaceFilling": False,
            },
        ),
    )

# %%
#
# Showing custom bonds beyond the existing topology
# ++++++++++++++++++++++++++++++++++++++++++++++++++
#
# With the custom bond features, we can now overlay connections of interest on
# the structure and existing topology by using, the
# `convert_stk_bonds_to_shapes` function.
# For example, here, we show the cage topology graph of metal-organic cages
# as defined by the metal atoms only. We start by building a series of
# Pd_nL_n2 metal-organic cages, with n=2,3,4,6 and the same organic ligand,
# using stk.


structures = [
    # A metal-organic cage.
    stk.ConstructedMolecule(
        func(
            building_blocks=(
                stk.BuildingBlock(
                    smiles="[Pd+2]",
                    functional_groups=(
                        stk.SingleAtom(stk.Pd(0, charge=2)) for i in range(4)
                    ),
                    position_matrix=[[0.0, 0.0, 0.0]],
                ),
                stk.BuildingBlock(
                    smiles=("C1=NC=CC(C2=CC=CC(C3=CC=NC=C3)=C2)=C1"),
                    functional_groups=[
                        stk.SmartsFunctionalGroupFactory(
                            smarts="[#6]~[#7X2]~[#6]",
                            bonders=(1,),
                            deleters=(),
                        ),
                    ],
                ),
            ),
            # Ensure that bonds between the
            # GenericFunctionalGroups of the ligand and the
            # SingleAtom functional groups of the metal are
            # dative.
            reaction_factory=stk.DativeReactionFactory(
                stk.GenericReactionFactory(
                    bond_orders={
                        frozenset(
                            {
                                stk.GenericFunctionalGroup,
                                stk.SingleAtom,
                            }
                        ): 9,
                    },
                ),
            ),
        ),
    )
    for func in (
        stk.cage.M2L4Lantern,
        stk.cage.M3L6,
        stk.cage.M4L8,
        stk.cage.M6L12Cube,
    )
]


# %%
#
# Again, we write some properties.

energy = stko.UFFEnergy(ignore_inter_interactions=False)
shape_calc = stko.ShapeCalculator()
properties = {
    "uffenergy": [energy.get_energy(molecule) for molecule in structures],
    "aspheriticty": [
        shape_calc.get_results(molecule).get_asphericity() for molecule in structures
    ],
}


# %%
#
# Now, we use some stk features to extract the metal atoms, and create "fake"
# bonds between them.

metal_atoms = [
    # Atoms with atomic number 46 (Pd).
    [i.get_id() for i in molecule.get_atoms() if i.get_atomic_number() == 46]
    for molecule in structures
]

fake_bonds = [
    # Combinations of those atoms, this does not filter for nearest neighbors.
    [(a1id, a2id) for a1id, a2id in it.combinations(metal_atoms[i], 2)]
    for i, molecule in enumerate(structures)
]

# %%
#
# With these fake bonds, we can create a new list of stk structures only having
# the fake bonds, to extract the shape dictionary defined by only the fake
# bonds. Note that these fake bonds could be arbitrarily defined, as long as
# they map to the atoms in the original structure.

structures_with_pd_pd_bonds = [
    stk.BuildingBlock.init(
        atoms=struct.get_atoms(),
        # Only including fake bonds.
        bonds=tuple(
            stk.Bond(
                atom1=next(struct.get_atoms(a1id)),
                atom2=next(struct.get_atoms(a2id)),
                order=1,
            )
            for a1id, a2id in fake_bonds[i]
        ),
        position_matrix=struct.get_position_matrix(),
    )
    for i, struct in enumerate(structures)
]

# %%
#
# Another feature of the shape dictionary is that the user can alter the color
# and radius of the bonds.

shape_dict = chemiscope.convert_stk_bonds_as_shapes(
    structures=structures_with_pd_pd_bonds,
    bond_color="#fc5500",
    bond_radius=0.2,
)

# Write the shape string for settings to turn them on automatically.
shape_string = ",".join(shape_dict.keys())


# %%
#
# Now, we show the structure with the new fake bonds overlaid
# (`bonds` and `shape` on).

chemiscope.show(
    structures=structures,
    properties=properties,
    shapes=shape_dict,
    settings=chemiscope.quick_settings(
        x="aspheriticty",
        y="uffenergy",
        color="",
        structure_settings={
            "shape": shape_string,
            "atoms": True,
            "bonds": True,
            "spaceFilling": False,
        },
    ),
)

# %%
#
# Write to json file.

chemiscope.write_input(
    path="shape_example.json.gz",
    structures=structures,
    properties=properties,
    meta=dict(name="Added Pd-Pd bonds overlaid with the stk molecule."),
    settings=chemiscope.quick_settings(
        x="aspheriticty",
        y="uffenergy",
        color="",
        structure_settings={
            "shape": shape_string,
            "atoms": True,
            "bonds": True,
            "spaceFilling": False,
        },
    ),
    shapes=shape_dict,
)
