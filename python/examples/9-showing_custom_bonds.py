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

import ase.io
import stk
import stko
from rdkit.Chem import AllChem as rdkit

import chemiscope


# %%
#
# Generate a list of stk BuildingBlocks (representation of a molecule) with
# properties. This also includes working with rdkit, which comes installed
# with stk.

rdkitmol = rdkit.MolFromSmiles("Cc1ccccc1")
rdkitmol = rdkit.AddHs(rdkitmol)
rdkit.Kekulize(rdkitmol)
params = rdkit.ETKDGv3()
params.randomSeed = 0xF00D
rdkit.EmbedMolecule(rdkitmol, params)

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
# Write their properties using any method, here we show using stko:
# https://stko-docs.readthedocs.io/en/latest/

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
# A chemiscope widget showing the result with standard, stk bonding.

chemiscope.show(
    frames=structures,
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
    frames=structures,
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
# View the same molecule imported with ASE versus stk, where chemiscope will
# automatically detect the bonds because ASE does not contain this information.

structures[0].write("data/stk_cage.xyz")
chemiscope.show(
    frames=[ase.io.read("data/stk_cage.xyz")],
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
# Writing to a json.gz file.

chemiscope.write_input(
    path="vsase_example.json.gz",
    frames=[ase.io.read("data/stk_cage.xyz")],
    properties={i: [properties[i][0]] for i in properties},
    meta=dict(name="Comparing with stk and ase bonding."),
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
# With the custom bond features, we can now overlay connections of interest on
# the structure and existing topology. For example, showing the cage topology
# in metal-organic cages.


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
# Write their properties using any method, here we show using stko:
# https://stko-docs.readthedocs.io/en/latest/

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
# Get the stk bonding information and convert them into shapes.
metal_atoms = [
    [i.get_id() for i in molecule.get_atoms() if i.get_atomic_number() == 46]
    for molecule in structures
]

bonds = [
    [(a1id, a2id) for a1id, a2id in it.combinations(metal_atoms[i], 2)]
    for i, molecule in enumerate(structures)
]

structures_with_pd_pd_bonds = [
    stk.BuildingBlock.init(
        atoms=struct.get_atoms(),
        bonds=tuple(
            stk.Bond(
                atom1=next(struct.get_atoms(a1id)),
                atom2=next(struct.get_atoms(a2id)),
                order=1,
            )
            for a1id, a2id in bonds[i]
        ),
        position_matrix=struct.get_position_matrix(),
    )
    for i, struct in enumerate(structures)
]

shape_dict = chemiscope.convert_stk_bonds_as_shapes(
    frames=structures_with_pd_pd_bonds,
    bond_color="#fc5500",
    bond_radius=0.2,
)

# Write the shape string for settings to turn them on automatically.
shape_string = ",".join(shape_dict.keys())


# %%
#
# Now with added bonding information.

chemiscope.show(
    frames=structures,
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
# Write to json file with added shapes.

chemiscope.write_input(
    path="shape_example.json.gz",
    frames=structures,
    properties=properties,
    meta=dict(name="Added all stk bonds."),
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
