"""
Showing custom bonds
====================

This example demonstrates how to add shapes into the chemiscope output such
that custom bonds that would not normally be assigned by ASE can be rendered.

"""

# %%
#
import stk
import stko
import ase.io
import pathlib
import chemiscope

working_path = pathlib.Path(__file__).resolve().parent
data_path = working_path / "data"

########################################################################
# This is just for development, where I build the library of structures.
energy = stko.UFFEnergy()
analyser = stko.molecule_analysis.GeometryAnalyser()
xyz_path = data_path / "stkbuilt.xyz"
string = ""

host = stk.ConstructedMolecule(
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
molecule = stk.ConstructedMolecule(
    topology_graph=stk.host_guest.Complex(
        host=stk.BuildingBlock.init_from_molecule(host),
        guests=stk.host_guest.Guest(
            building_block=stk.BuildingBlock("[Br][Br]"),
        ),
    ),
)
teststring = stk.XyzWriter().to_string(molecule)
teststring = teststring.split("\n")
teststring[1] = (
    f"e='{energy.get_energy(molecule)}' "
    f"rg='{analyser.get_radius_gyration(molecule)}'"
)
string += "\n".join(teststring)

bb1 = stk.BuildingBlock(
    smiles="NCCN",
    functional_groups=[stk.PrimaryAminoFactory()],
)
bb2 = stk.BuildingBlock(
    smiles="O=CC(C=O)C=O",
    functional_groups=[stk.AldehydeFactory()],
)
molecule = stk.ConstructedMolecule(
    topology_graph=stk.cage.FourPlusSix((bb1, bb2)),
)
molecule.write(data_path / "unoptcage.mol")
teststring = stk.XyzWriter().to_string(molecule)
teststring = teststring.split("\n")
teststring[1] = (
    f"e='{energy.get_energy(molecule)}' "
    f"rg='{analyser.get_radius_gyration(molecule)}'"
)
string += "\n".join(teststring)

molecule = stk.ConstructedMolecule(
    topology_graph=stk.cage.FourPlusSix(
        building_blocks=(bb1, bb2),
        optimizer=stk.MCHammer(),
    )
)
teststring = stk.XyzWriter().to_string(molecule)
teststring = teststring.split("\n")
teststring[1] = (
    f"e='{energy.get_energy(molecule)}' "
    f"rg='{analyser.get_radius_gyration(molecule)}'"
)
string += "\n".join(teststring)


palladium_atom = stk.BuildingBlock(
    smiles="[Pd+2]",
    functional_groups=(stk.SingleAtom(stk.Pd(0, charge=2)) for i in range(4)),
    position_matrix=[[0.0, 0.0, 0.0]],
)
bb1 = stk.BuildingBlock(
    smiles=("C1=NC=CC(C2=CC=CC(C3=C" "C=NC=C3)=C2)=C1"),
    functional_groups=[
        stk.SmartsFunctionalGroupFactory(
            smarts="[#6]~[#7X2]~[#6]",
            bonders=(1,),
            deleters=(),
        ),
    ],
)
molecule = stk.ConstructedMolecule(
    stk.cage.M2L4Lantern(
        building_blocks=(palladium_atom, bb1),
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
molecule.write(data_path / "metalcage.mol")
teststring = stk.XyzWriter().to_string(molecule)
teststring = teststring.split("\n")
teststring[1] = (
    f"e='{energy.get_energy(molecule)}' "
    f"rg='{analyser.get_radius_gyration(molecule)}'"
)
string += "\n".join(teststring)


with xyz_path.open("w") as f:
    f.write(string)

########################################################################


print(
    "one change I am interested in, is if you want to use the standard "
    "load xyz > chemiscope workflow, like here. "
    "or, if I use rdkit/stk to provide the graph automatically? "
    "I would not want to add any dependancies in doing so"
)


# %%
#
# Load structures from an extended xyz file.
frames = ase.io.read(xyz_path, ":")


# %%
#
# Write to json file without added shapes and note that the XX-XX bonds are
# missing.

chemiscope.write_input(
    "noshape_example.json.gz",
    frames=frames,
    properties=chemiscope.extract_properties(frames, only=["e", "rg"]),
    meta=dict(name="Missing bonds by automation."),
    settings=chemiscope.quick_settings(x="rg", y="e", color=""),
)


# %%
#
# Compute the shape dicionary from the structures based on XX.
mol1 = stk.BuildingBlock.init_from_file(data_path / "unoptcage.mol")
mol1_bonds = tuple(
    (bond.get_atom1().get_id(), bond.get_atom2().get_id()) for bond in mol1.get_bonds()
)

mol3 = stk.BuildingBlock.init_from_file(data_path / "metalcage.mol")
mol3_bonds = tuple(
    (bond.get_atom1().get_id(), bond.get_atom2().get_id()) for bond in mol3.get_bonds()
)
bonds_per_frames = {
    # Here, we do it by hand to fill in some issues.
    0: ((0, 60), (22, 132), (30, 85), (34, 150), (4, 96), (46, 169)),
    2: ((3, 48), (11, 56), (19, 72)),
    # In these two examples, we are using "fake", unoptimised models and show
    # how to generate these from a mol file.
    # Although these examples currently require a new dependancy, so I can
    # alter them as needed?
    1: mol1_bonds,
    3: mol3_bonds,
}

shape_dict = chemiscope.convert_bonds_as_shapes(
    frames=frames,
    bonds_per_frames=bonds_per_frames,
)
# Write the shape string for settings to turn them on automatically.
shape_string = ",".join(shape_dict.keys())

# %%
#
# Write to json file with added shapes and note that the XX-XX bonds are
# present.
chemiscope.write_input(
    "shape_example.json.gz",
    frames=frames,
    properties=chemiscope.extract_properties(frames, only=["e", "rg"]),
    meta=dict(name="Added some bonds."),
    settings=chemiscope.quick_settings(
        x="rg",
        y="e",
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

raise SystemExit
