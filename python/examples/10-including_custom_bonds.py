"""
Including custom bonds
======================

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

import chemiscope


# %%
#
# Generate a list of stk BuildingBlocks (representation of a molecule) with
# properties. This also includes working with rdkit, which comes installed
# with stk.

structures = [
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
# Now with added bonding information.

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
# Write to json file with added shapes.

chemiscope.write_input(
    path="10-shape_example.json.gz",
    frames=structures,
    properties=properties,
    meta=dict(name="Added all stk bonds."),
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
