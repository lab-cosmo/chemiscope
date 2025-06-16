from pathlib import Path
from typing import Dict, List, Optional, Union

import metatensor.torch as mts
import metatomic.torch as mta
import torch
from torch import nn


class MLPProjector(nn.Module):
    """A simple MLP used to project feature vectors to low dimention representations"""

    def __init__(self, input_dim: int = 1024, output_dim: int = 3):
        super().__init__()

        self.input_dim = input_dim
        self.output_dim = output_dim

        self.fc1 = nn.Linear(self.input_dim, 512)
        self.fc2 = nn.Linear(512, 256)
        self.fc3 = nn.Linear(256, 128)
        self.output = nn.Linear(128, self.output_dim)

    def forward(self, features: torch.Tensor) -> torch.Tensor:
        x = torch.relu(self.fc1(features))
        x = torch.relu(self.fc2(x))
        x = torch.relu(self.fc3(x))
        return self.output(x)


class MADExplorer(nn.Module):
    """
    Metatomic wrapper model for that extracts PET-MAD last-layer features and
    projects them into a low-dimensional space using an MLP.

    The model is intended for exploratory analysis and visualization of the
    learned representations.

    :param mtt_model: path to a saved PET-MAD model or an instance of a loaded
        model
    :param extensions_directory: path to model extensions
    :param check_consistency: whether to verify consistency between model and
        system inputs
    :param input_dim : dimensionality of the input PET-MAD features
    :param output_dim: target low dimensionality for the projected features
    :param device: device on which to run the model
    """

    def __init__(
        self,
        mtt_model: Union[str, Path, mta.AtomisticModel],
        mlp_checkpoint: Optional[str] = None,
        extensions_directory: str = None,
        check_consistency: bool = None,
        input_dim: int = 1024,
        output_dim: int = 3,
        device: Optional[Union[str, torch.device]] = "cpu",
    ):
        super().__init__()

        self.check_consistency = check_consistency
        self.device = device

        if isinstance(mtt_model, (str, Path)):
            self.petmad = mta.load_atomistic_model(
                mtt_model, extensions_directory=extensions_directory
            )
        else:
            self.petmad = mtt_model

        capabilities = self.petmad.capabilities()

        self.output_name = "mtt::aux::energy_last_layer_features"
        if self.output_name not in capabilities.outputs:
            raise ValueError(f"this model does not have a '{self.output_name}' output")

        if capabilities.dtype == "float32":
            self.dtype = torch.float32
        else:
            assert capabilities.dtype == "float64"
            self.dtype = torch.float64

        self.projector = MLPProjector(input_dim, output_dim).to(self.device)

        if mlp_checkpoint:
            self.projector_checkpoint = torch.load(mlp_checkpoint, weights_only=False)
            self.projector.load_state_dict(
                self.projector_checkpoint["projector_state_dict"]
            )
        else:
            self.projector_checkpoint = None

    def forward(
        self,
        systems: List[mta.System],
        outputs: Dict[str, mta.ModelOutput],
        selected_atoms: Optional[mts.Labels],
    ) -> Dict[str, mts.TensorMap]:
        if list(outputs.keys()) != ["features"]:
            raise ValueError(
                f"`outputs` keys ({', '.join(outputs.keys())}) contain unsupported "
                "keys. Only 'features' is supported"
            )

        if not outputs["features"].per_atom:
            raise NotImplementedError(
                "Model uses per-atom features to get mean and std features"
            )

        if selected_atoms is None:
            raise ValueError("MADExplorer requires 'selected_atoms' to be provided")

        length_unit = self.petmad.capabilities().length_unit
        options = mta.ModelEvaluationOptions(
            length_unit=length_unit,
            outputs={
                self.output_name: mta.ModelOutput(per_atom=selected_atoms is not None)
            },
            selected_atoms=selected_atoms,
        )

        systems = [s.to(self.dtype, self.device) for s in systems]
        features = self._get_descriptors(systems, options)

        high_scaler = self.projector_checkpoint["scaler_highdim"]
        low_scaler = self.projector_checkpoint["scaler_lowdim"]

        features_np = features.numpy()

        if high_scaler:
            scaled_features = high_scaler.transform(features_np)
            features_tensor = torch.FloatTensor(scaled_features)
        else:
            features_tensor = torch.FloatTensor(features)

        with torch.no_grad():
            X_reduced = self.projector(features_tensor)

        projections = X_reduced.detach().numpy()

        if low_scaler:
            projections = low_scaler.inverse_transform(X_reduced)

        num_atoms, num_projections = projections.shape
        sample_labels = mts.Labels("system", torch.arange(num_atoms).reshape(-1, 1))
        prop_labels = mts.Labels(
            "projection", torch.arange(num_projections).reshape(-1, 1)
        )

        projections_tensor = torch.FloatTensor(projections)

        block = mts.TensorBlock(
            values=projections_tensor,
            samples=sample_labels,
            components=[],
            properties=prop_labels,
        )

        tensor_map = mts.TensorMap(
            keys=mts.Labels("_", torch.tensor([[0]])), blocks=[block]
        )

        return {"features": tensor_map}

    def _get_descriptors(
        self, systems: List[mta.System], options: mta.ModelEvaluationOptions
    ) -> torch.Tensor:
        """
        Compute embeddings for the given systems using the PET-MAD model.

        The method computes per-atom mean and std of features and returns as a
        combined tensor.
        """

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
            return combined_features.detach()
        else:
            return features.block().values.detach()

    def get_atomic_types(self) -> List[int]:
        return self.petmad.capabilities().atomic_types


model = MADExplorer("pet-mad-latest.pt", mlp_checkpoint="mtt_projection_model.pt")

metadata = mta.ModelMetadata(
    name="mad-explorer",
    description="Exploration tool for PET-MAD model features upon SMAP projections",
    authors=["TODO"],
    references={
        "architecture": ["https://arxiv.org/abs/2305.19302v3"],
        "model": ["http://arxiv.org/abs/2503.14118"],
        "implementation": ["https://doi.org/10.1073/pnas.1108486108"],
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

mad_explorer = mta.AtomisticModel(model.eval(), metadata, capabilities)
