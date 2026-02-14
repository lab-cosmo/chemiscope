import warnings
from typing import Sequence

import numpy as np

from ._ase import _remove_invalid_properties


try:
    import chemfiles

    if chemfiles.__version__ < "0.10.0" or chemfiles.__version__ > "0.11.0":
        print(
            "chemiscope requires chemfiles version >=0.10,<0.11; "
            f"but version {chemfiles.__version__} is installed."
        )
        HAVE_CHEMFILES = False
    else:
        HAVE_CHEMFILES = True

    from chemfiles import Frame, MemoryTrajectory, Trajectory

except ImportError:
    HAVE_CHEMFILES = False


def _chemfiles_valid_structures(structures):
    if not HAVE_CHEMFILES:
        return structures, False

    if isinstance(structures, (Trajectory, MemoryTrajectory)):
        return [structures.read_step(i) for i in range(structures.nsteps)], True

    elif isinstance(structures, Frame):
        return [structures], True

    elif (
        isinstance(structures, Sequence)
        and len(structures) > 0
        and isinstance(structures[0], Frame)
    ):
        for structure in structures:
            assert isinstance(structure, Frame)
        return structures, True

    else:
        return structures, False


# fmt: off
ELEMENTS = [
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si",
    "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni",
    "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo",
    "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
    "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
]
# fmt: on

# Standard protein and nucleic acid residues (non-hetatoms)
STANDARD_RESIDUES = {
    # Amino acids
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    # Protonation states / common non-standard names often treated as standard
    "HID",
    "HIE",
    "HIP",
    "CYX",
    "ASH",
    "GLH",
    "LYN",
    # Nucleic acids (DNA/RNA)
    "DA",
    "DC",
    "DG",
    "DT",
    "DI",
    "A",
    "C",
    "G",
    "U",
    "I",
}


def _chemfiles_to_json(frame):
    """Implementation of structures_to_json for chemfiles' ``Frame``."""

    BOND_ORDERS_TO_NUMERIC = {
        chemfiles.BondOrder.Unknown: 1,
        chemfiles.BondOrder.Single: 1,
        chemfiles.BondOrder.Double: 2,
        chemfiles.BondOrder.Triple: 3,
        # 3Dmol seems to cap bond order at 3
        chemfiles.BondOrder.Quadruple: 3,
        chemfiles.BondOrder.Quintuplet: 3,
        chemfiles.BondOrder.Amide: 1,
        chemfiles.BondOrder.Aromatic: 4,
    }

    data = {}
    data["size"] = len(frame.atoms)
    data["names"] = [atom.name for atom in frame.atoms]

    elements = []
    all_have_element = True
    for atom in frame.atoms:
        atomic_number = atom.atomic_number
        if atomic_number == 0:
            all_have_element = False
            break
        else:
            elements.append(ELEMENTS[atomic_number])

    if not all_have_element:
        # If there are unknown elements we assume the element detection failed,
        # and try a more conservative approach assuming the elements is stored
        # in atom type names. Chemfiles uses a similar heuristic, but maps CA
        # and CD to Ca and Cd, while in all likelihood they should be carbon atoms.
        elements = []
        for atom in frame.atoms:
            name = atom.name
            if not name:
                elements.append("X")
                continue
            if name[0] in "HBCNOPS":
                elements.append(name[0])
            elif len(name) >= 2 and name[0:2].capitalize() in ELEMENTS:
                elements.append(name[0:2].capitalize())
            elif name[0] in "FIUKV":
                elements.append(name[0])
            else:
                elements.append("X")

    data["elements"] = elements

    positions = frame.positions
    data["x"] = [float(positions[i][0]) for i in range(data["size"])]
    data["y"] = [float(positions[i][1]) for i in range(data["size"])]
    data["z"] = [float(positions[i][2]) for i in range(data["size"])]

    if frame.cell.shape != chemfiles.CellShape.Infinite:
        data["cell"] = frame.cell.matrix.T.flatten().tolist()

    # bonds
    topology = frame.topology
    if len(topology.bonds) > 0:
        data["bonds"] = [
            [int(bond[0]), int(bond[1]), BOND_ORDERS_TO_NUMERIC[order]]
            for bond, order in zip(topology.bonds, topology.bonds_orders, strict=True)
        ]

    # biomolecule-specific information
    chains = []
    resnames = []
    resids = []
    hetatom = []
    has_biomol_info = False
    for atom_i in range(data["size"]):
        residue = frame.topology.residue_for_atom(atom_i)
        if residue is None:
            resids.append(-1)
            resnames.append("UNK")
            chains.append("")
            hetatom.append(True)
            continue

        has_biomol_info = True
        resids.append(residue.id)
        resnames.append(residue.name)

        residue_properties = residue.list_properties()
        if "chainname" in residue_properties:
            chains.append(residue["chainname"])
        else:
            chains.append("")

        if "is_standard_pdb" in residue_properties:
            hetatom.append(not residue["is_standard_pdb"])
        else:
            # Fallback: check if residue name is in standard list
            # We strip whitespace and uppercase just in case
            is_standard = residue.name.strip().upper() in STANDARD_RESIDUES
            hetatom.append(not is_standard)

    if has_biomol_info:
        data["chains"] = chains
        data["resnames"] = resnames
        data["resids"] = resids
        data["hetatom"] = hetatom

    return data


def _chemfiles_all_atomic_environments(structures, cutoff):
    environments = []
    for frame_i, frame in enumerate(structures):
        for atom_i in range(len(frame.atoms)):
            environments.append((frame_i, atom_i, cutoff))

    return environments


def _chemfiles_get_structure_properties(frames):
    # extract the set of common properties between all frames
    all_properties = {}
    extra = set()

    if len(frames) == 0:
        return all_properties

    for name in frames[0].list_properties():
        all_properties[name] = [frames[0][name]]

    for frame in frames[1:]:
        current_properties = frame.list_properties()
        for name in current_properties:
            if name in all_properties:
                all_properties[name].append(frame[name])
            else:
                extra.add(name)

        for name in list(all_properties.keys()):
            if name not in current_properties:
                all_properties.pop(name, None)
                extra.add(name)

    if len(extra) != 0:
        warnings.warn(
            "the following structure properties are only defined for a subset "
            f"of structures: {list(sorted(extra))}; they will be ignored",
            stacklevel=2,
        )

    # ensures that if a property is a mix of strings and numbers, everything
    # is converted to string (as these should be categorical properties)
    for name in all_properties.keys():
        property = all_properties[name]
        if any(isinstance(x, str) for x in property):
            all_properties[name] = [str(x) for x in property]

    return all_properties


def _chemfiles_get_atom_properties(frames, environments):
    assert environments is not None
    # extract the set of common properties between all atoms in all frames
    all_properties = {}
    extra = set()

    if len(environments) == 0:
        return all_properties

    frame = frames[environments[0][0]]
    atom = frame.atoms[environments[0][1]]
    for name in atom.list_properties():
        all_properties[name] = [atom[name]]

    for frame_i, atom_i, _ in environments[1:]:
        frame = frames[frame_i]
        atom = frame.atoms[atom_i]

        current_properties = atom.list_properties()
        for name in current_properties:
            if name in all_properties:
                all_properties[name].append(atom[name])
            else:
                extra.add(name)

        for name in list(all_properties.keys()):
            if name not in current_properties:
                all_properties.pop(name, None)
                extra.add(name)

    if len(extra) != 0:
        warnings.warn(
            "the following atomic properties are only defined for a subset "
            f"of structures: {list(sorted(extra))}; they will be ignored",
            stacklevel=2,
        )

    return all_properties


def _chemfiles_extract_properties(structures, only=None, environments=None):
    """implementation of ``extract_properties`` for chemfiles"""
    all_properties = _chemfiles_get_structure_properties(structures)
    if only is None:
        selected = all_properties
    else:
        selected = {}
        for name in only:
            if name in all_properties.keys():
                selected[name] = all_properties[name]

    # create property in the format expected by create_input
    properties = {
        name: {"target": "structure", "values": np.asarray(value)}
        for name, value in selected.items()
    }

    _remove_invalid_properties(properties, "chemfiles")

    if environments is None:
        environments = _chemfiles_all_atomic_environments(structures, cutoff=0.0)

    atom_properties = _chemfiles_get_atom_properties(structures, environments)

    if only is None:
        selected = atom_properties
    else:
        selected = {}
        for name in only:
            if name in atom_properties.keys():
                selected[name] = atom_properties[name]

    # create property in the format expected by create_input
    atom_properties = {
        name: {"target": "atom", "values": np.stack(value, axis=0)}
        for name, value in selected.items()
    }
    _remove_invalid_properties(atom_properties, "chemfiles")

    for name, values in atom_properties.items():
        if name in properties:
            warnings.warn(
                f"a property named '{name}' is defined for both atoms and structures, "
                "the atom one will be ignored",
                stacklevel=2,
            )
        else:
            properties[name] = values

    return properties
