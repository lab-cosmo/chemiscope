import unittest
from typing import Dict, List, Optional

import ase.build
import numpy as np
import torch
from metatensor.torch import Labels, TensorBlock, TensorMap
from metatomic.torch import (
    AtomisticModel,
    ModelCapabilities,
    ModelMetadata,
    ModelOutput,
    System,
)

import chemiscope


STRUCTURES = [ase.build.molecule(name) for name in ["H2O", "NH3", "CH4", "CH3CH2OH"]]
for i, atoms in enumerate(STRUCTURES):
    atoms.positions += i


class SumOfPositions(torch.nn.Module):
    def forward(
        self,
        systems: List[System],
        outputs: Dict[str, ModelOutput],
        selected_atoms: Optional[Labels] = None,
    ) -> Dict[str, TensorMap]:
        values = torch.stack([system.positions.sum(dim=0) for system in systems])

        block = TensorBlock(
            values=values,
            samples=Labels(["system"], torch.arange(len(systems)).reshape(-1, 1)),
            components=[],
            properties=Labels(["xyz"], torch.tensor([[0], [1], [2]])),
        )

        return {"feature": TensorMap(Labels(["_"], torch.tensor([[0]])), [block])}


def sum_of_positions_model():
    capabilities = ModelCapabilities(
        outputs={"feature": ModelOutput(sample_kind="system")},
        atomic_types=[1, 6, 7, 8],
        length_unit="angstrom",
        interaction_range=0.0,
        supported_devices=["cpu"],
        dtype="float64",
    )

    return AtomisticModel(SumOfPositions().eval(), ModelMetadata(), capabilities)


class TestExplore(unittest.TestCase):
    def test_batches(self):
        expected = np.array([atoms.positions.sum(axis=0) for atoms in STRUCTURES])

        for batch_size in [None, 1, 3]:
            featurizer = chemiscope.metatomic_featurizer(
                sum_of_positions_model(), batch_size=batch_size
            )
            features = featurizer(STRUCTURES, None)

            np.testing.assert_allclose(features, expected)

    def test_ecumetric(self):
        try:
            featurizer = chemiscope.get_featurizer("ecumetric-1.0")
        except OSError as e:
            self.skipTest(f"could not download the ecumetric model: {e}")

        features = featurizer(STRUCTURES, None)

        self.assertEqual(features.shape, (len(STRUCTURES), 8))


if __name__ == "__main__":
    unittest.main()
