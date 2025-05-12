"""
Showing custom bonds using stk
==============================

This example demonstrates how to add shapes into the chemiscope output such
that custom bonds that would not automatically be assigned can be rendered.

This is done by using `stk <https://stk.readthedocs.io/en/stable/>`_ to
generate and analyse molecules, which comes with topology/bonding information
by default (using the cheminformatic software rdkit).

We use `stko <https://stko-docs.readthedocs.io/en/latest/>`_ to calculate
some rudimentary properties of `stk` molecules. `stko` can be installed with
``pip install stko``.

"""

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

structures = [
    # A building block.
    stk.BuildingBlock(smiles="NCCN"),
    # A mostly optimised cage molecule.
    stk.ConstructedMolecule(
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
    ),
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
    # A host guest molecule.
    stk.ConstructedMolecule(
        topology_graph=stk.host_guest.Complex(
            host=stk.BuildingBlock.init_from_molecule(
                stk.ConstructedMolecule(
                    topology_graph=stk.cage.FourPlusSix(
                        building_blocks=(
                            stk.BuildingBlock(
                                smiles="NC1CCCCC1N",
                                functional_groups=[
                                    stk.PrimaryAminoFactory(),
                                ],
                            ),
                            stk.BuildingBlock(
                                smiles="O=Cc1cc(C=O)cc(C=O)c1",
                                functional_groups=[stk.AldehydeFactory()],
                            ),
                        ),
                        optimizer=stk.MCHammer(),
                    ),
                )
            ),
            guests=stk.host_guest.Guest(
                building_block=stk.BuildingBlock("[Br][Br]"),
            ),
        ),
    ),
    # From rdkit.
    stk.BuildingBlock.init_from_rdkit_mol(rdkitmol),
]


# %%
#
# Write their properties using any method, here we show using stko:
# https://stko-docs.readthedocs.io/en/latest/

energy = stko.UFFEnergy()
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
shape_dict = chemiscope.convert_stk_bonds_as_shapes(
    frames=structures,
    bond_color="#fc5500",
    bond_radius=0.12,
)

# Write the shape string for settings to turn them on automatically.
shape_string = ",".join(shape_dict.keys())


# %%
#
# A chemiscope widget showing the result without the added bonding.

chemiscope.show(frames=structures, properties=properties)


# %%
#
# Writing to a json.gz file, again without added bonding.

chemiscope.write_input(
    path="noshape_example.json.gz",
    frames=structures,
    properties=properties,
    meta=dict(name="Missing bonds by automation."),
    settings=chemiscope.quick_settings(x="aspheriticty", y="uffenergy", color=""),
)


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
            "bonds": False,
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
            "bonds": False,
            "spaceFilling": False,
        },
    ),
    shapes=shape_dict,
)
