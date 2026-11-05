import io
import unittest

import numpy as np

from chemiscope.headless import headless


try:
    from PIL import Image

    HAVE_PIL = True
except ImportError:
    HAVE_PIL = False


def polygon(names, center, radius):
    """Atoms on a regular polygon, standing in for real ring geometries"""
    angles = 2 * np.pi * np.arange(len(names)) / len(names)
    return {
        name: np.asarray(center) + radius * np.array([np.cos(a), np.sin(a), 0])
        for name, a in zip(names, angles, strict=True)
    }


def residue_structure(residues):
    """Chemiscope structure from a list of (resname, hetatom, {name: position})"""
    structure = {
        key: []
        for key in ["names", "x", "y", "z", "resnames", "resids", "chains", "hetatom"]
    }
    for resid, (resname, hetatom, atoms) in enumerate(residues):
        for name, position in atoms.items():
            structure["names"].append(name)
            for key, value in zip("xyz", position, strict=True):
                structure[key].append(float(value))
            structure["resnames"].append(resname)
            structure["resids"].append(resid + 1)
            structure["chains"].append("A")
            structure["hetatom"].append(hetatom)
    structure["size"] = len(structure["names"])
    structure["elements"] = [name[0] for name in structure["names"]]
    return structure


def guanosine():
    """Atoms of the base and O5' of a guanosine, as (resname, hetatom, atoms)"""
    atoms = polygon(["N1", "C2", "N3", "C4", "C5", "C6"], [0, 0, 0], 1.39)
    out = (atoms["C4"] + atoms["C5"]) / 2
    out /= np.linalg.norm(out)
    edge = (atoms["C4"] - atoms["C5"]) / np.linalg.norm(atoms["C4"] - atoms["C5"])
    atoms["N9"] = atoms["C4"] + 1.2 * out + 0.4 * edge
    atoms["N7"] = atoms["C5"] + 1.2 * out - 0.4 * edge
    atoms["C8"] = (atoms["C4"] + atoms["C5"]) / 2 + 2.0 * out
    atoms["O5'"] = atoms["N9"] + 4.0 * out
    return ("DG", False, atoms)


def green_pixels(png):
    """Count pixels which are clearly green, which no CPK element color is"""
    image = np.asarray(Image.open(io.BytesIO(png)).convert("RGB")).astype(int)
    r, g, b = image[..., 0], image[..., 1], image[..., 2]
    return int(np.sum((g > 80) & (g > r + 40) & (g > b + 30)))


@unittest.skipIf(not HAVE_PIL, "Pillow is required to check the images")
class TestBiomolCartoons(unittest.TestCase):
    def render(self, structure, settings):
        widget = headless(
            structures=[structure],
            properties={"index": {"target": "structure", "values": [0]}},
            mode="structure",
        )
        try:
            return widget.get_structure_sequence([0] * len(settings), settings=settings)
        finally:
            widget.close()

    def test_snfg_symbols(self):
        # mannose is a green sphere in SNFG
        ring = polygon(["C1", "C2", "C3", "C4", "C5", "O5"], [0, 0, 0], 1.45)
        structure = residue_structure([("MAN", True, ring)])

        snfg, atoms = self.render(structure, [{"snfg": True}, {"snfg": False}])
        self.assertGreater(green_pixels(snfg), 100)
        self.assertEqual(green_pixels(atoms), 0)

    def test_not_a_saccharide(self):
        # GLYCAM-like names are only used when the residue contains a ring
        atoms = {f"C{i + 1}": [1.6 * i, 0, 0] for i in range(5)}
        atoms["O5"] = [8.0, 0, 0]
        structure = residue_structure([("0MA", True, atoms)])

        (image,) = self.render(structure, [{"snfg": True}])
        self.assertEqual(green_pixels(image), 0)

    def test_nucleobases(self):
        # guanine is drawn as a green ring in cartoon mode
        structure = residue_structure([guanosine()])

        cartoon, atoms = self.render(structure, [{"cartoon": True}, {"cartoon": False}])
        self.assertGreater(green_pixels(cartoon), 100)
        self.assertEqual(green_pixels(atoms), 0)


if __name__ == "__main__":
    unittest.main()
