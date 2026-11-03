import io
import unittest

import numpy as np

import chemiscope


try:
    import openmm.app
    import openmm.unit
    from openmm import Vec3

    HAVE_OPENMM = True
except ImportError:
    HAVE_OPENMM = False


# positions in nm, OpenMM's native unit
POSITIONS = [[0, 0, 0], [0.1, 0.1, 0.1], [0.2, 0.2, 0.5], [0.3, 0.3, 0.3]]


def create_test_modeller(box=None):
    topology = openmm.app.Topology()
    element = openmm.app.Element.getBySymbol

    chain = topology.addChain("A")
    residue = topology.addResidue("GLY", chain, id="7")
    n = topology.addAtom("N", element("N"), residue)
    ca = topology.addAtom("CA", element("C"), residue)
    topology.addBond(ca, n)

    chain = topology.addChain("B")
    residue = topology.addResidue("LIG", chain, id="101")
    c = topology.addAtom("C1", element("C"), residue)
    o = topology.addAtom("O1", element("O"), residue)
    topology.addBond(c, o)

    if box is not None:
        topology.setPeriodicBoxVectors(box)

    positions = [Vec3(*p) for p in POSITIONS] * openmm.unit.nanometer
    return openmm.app.Modeller(topology, positions)


@unittest.skipUnless(HAVE_OPENMM, "openmm is not installed")
class TestStructures(unittest.TestCase):
    """Conversion of structure data to chemiscope JSON"""

    def test_structures(self):
        data = chemiscope.create_input(create_test_modeller())
        self.assertEqual(len(data["structures"]), 1)

        structure = data["structures"][0]
        self.assertEqual(structure["size"], 4)
        self.assertEqual(structure["names"], ["N", "CA", "C1", "O1"])
        self.assertEqual(structure["elements"], ["N", "C", "C", "O"])
        # nm -> Å
        self.assertTrue(np.allclose(structure["x"], [0, 1, 2, 3]))
        self.assertTrue(np.allclose(structure["y"], [0, 1, 2, 3]))
        self.assertTrue(np.allclose(structure["z"], [0, 1, 5, 3]))
        self.assertEqual(structure["bonds"], [[1, 0, 1], [2, 3, 1]])
        self.assertEqual(structure.get("cell"), None)

        self.assertEqual(structure["resnames"], ["GLY", "GLY", "LIG", "LIG"])
        self.assertEqual(structure["resids"], [7, 7, 101, 101])
        self.assertEqual(structure["chains"], ["A", "A", "B", "B"])
        self.assertEqual(structure["hetatom"], [False, False, True, True])

    def test_cell(self):
        box = [Vec3(2.3, 0, 0), Vec3(0.5, 2.2, 0), Vec3(0, 0.1, 1.1)]
        modeller = create_test_modeller(box * openmm.unit.nanometer)
        data = chemiscope.create_input([modeller, modeller])
        self.assertEqual(len(data["structures"]), 2)
        self.assertTrue(
            np.allclose(data["structures"][0]["cell"], [23, 0, 0, 5, 22, 0, 0, 1, 11])
        )

    def test_pdbfile(self):
        modeller = create_test_modeller()
        expected = chemiscope.create_input(modeller)["structures"][0]

        buffer = io.StringIO()
        openmm.app.PDBFile.writeFile(
            modeller.topology, modeller.positions, buffer, keepIds=True
        )
        buffer.seek(0)
        structure = chemiscope.create_input(openmm.app.PDBFile(buffer))
        structure = structure["structures"][0]

        for key in "xyz":
            self.assertTrue(np.allclose(structure.pop(key), expected.pop(key)))
        self.assertEqual(structure, expected)

    def test_environments(self):
        modeller = create_test_modeller()
        environments = chemiscope.all_atomic_environments([modeller, modeller], 2.0)
        self.assertEqual(len(environments), 8)
        self.assertEqual(environments[5], (1, 1, 2.0))

        with self.assertRaises(RuntimeError):
            chemiscope.extract_properties(modeller)


if __name__ == "__main__":
    unittest.main()
