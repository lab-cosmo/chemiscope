import torch
from torch import nn
from pathlib import Path
import metatomic.torch as mta
import metatensor.torch as mts

from typing import Dict, List, Optional, Union

import warnings

import chemiscope


class MADExplorer(nn.Module):
    def __init__(
        self,
        model: Union[str, Path, mta.AtomisticModel],
        extensions_directory: str = None,
        check_consistency: bool = None,
        input_dim: int = 1024,
        output_dim: int = 2,
        device=None,
    ):
        super().__init__()

        self.input_dim = input_dim
        self.output_dim = output_dim

        if isinstance(model, (str, Path)):
            self.petmad = mta.load_atomistic_model(
                model, extensions_directory=extensions_directory
            )
        else:
            self.petmad = model

        self.check_consistency = check_consistency
        self.device = device

        # MLP layers
        self.fc1 = nn.Linear(input_dim, 512)
        self.fc2 = nn.Linear(512, 256)
        self.fc3 = nn.Linear(256, 128)
        self.output = nn.Linear(128, output_dim)

        capabilities = self.petmad.capabilities()

        self.output_name = "mtt::aux::energy_last_layer_features"
        if self.output_name not in capabilities.outputs:
            raise ValueError(f"this model does not have a '{self.output_name}' output")

        if capabilities.dtype == "float32":
            self.dtype = torch.float32
        else:
            assert capabilities.dtype == "float64"
            self.dtype = torch.float64

    def forward(
        self,
        systems: List[mta.System],
        outputs: Dict[str, mta.ModelOutput],
        selected_atoms: Optional[mts.Labels],
    ) -> Dict[str, mts.TensorMap]:
        if len(selected_atoms) == 0:
            raise ValueError("Provide 'selected_atoms' to see visualisation")

        systems = [s.to(self.dtype, self.device) for s in systems]

        capabilities = self.petmad.capabilities()

        options = mta.ModelEvaluationOptions(
            length_unit=capabilities.length_unit,
            outputs={
                self.output_name: mta.ModelOutput(
                    quantity="length",
                    unit="angstrom",
                    per_atom=selected_atoms is not None,
                )
            },
            selected_atoms=selected_atoms,
        )

        features = self._get_petmad_features(systems, options)
        projections = self._mlp_forward(features)

        block = mts.TensorBlock(
            values=projections,
            samples=mts.Labels("system", torch.arange(len(projections)).reshape(-1, 1)),
            components=[],
            properties=mts.Labels(
                "projection", torch.tensor([[i] for i in range(projections.shape[1])])
            ),
        )

        return {
            "features": mts.TensorMap(
                keys=mts.Labels("_", torch.tensor([[0]])), blocks=[block]
            )
        }

    def _mlp_forward(self, features):
        x = torch.relu(self.fc1(features))
        x = torch.relu(self.fc2(x))
        x = torch.relu(self.fc3(x))
        return self.output(x)

    def get_atomic_types(self):
        return self.petmad.capabilities().atomic_types

    def _get_petmad_features(self, systems, options):
        output = self.petmad(
            systems,
            options,
            check_consistency=self.check_consistency,
        )
        features = output[self.output_name]

        if options.selected_atoms is not None:
            mean = mts.mean_over_samples(features, "atom")
            mean_vals = torch.cat([block.values for block in mean.blocks()], dim=0)

            std = mts.std_over_samples(features, "atom")
            std_vals = torch.cat([block.values for block in std.blocks()], dim=0)

            combined_features = torch.cat([mean_vals, std_vals], dim=1)
            # return combined_features.detach().cpu().numpy()
            return combined_features.detach()
        else:
            return features.block().values.detach()
            # return features.block().values.detach().cpu()


model = MADExplorer("pet-mad-latest.pt")

metadata = mta.ModelMetadata(
    name="mad-explorer",
    description="TODO",
    authors=["TODO"],
    references={
        "architecture": ["https://arxiv.org/abs/2305.19302v3"],
        "model": ["http://arxiv.org/abs/2503.14118"],
    },
)

outputs = {
    "features": mta.ModelOutput(quantity="length", unit="angstrom", per_atom=True),
}

capabilities = mta.ModelCapabilities(
    outputs=outputs,
    length_unit="angstrom",
    supported_devices=["cpu"],
    dtype="float64",
    interaction_range=0.0,
    atomic_types=model.get_atomic_types(),
)

wrapper = mta.AtomisticModel(model.eval(), metadata, capabilities)

featurizer = chemiscope.metatomic_featurizer(wrapper, length_unit="angstrom")
