import warnings
from pathlib import Path

from ._metatomic import metatomic_featurizer


torch = None


class sketchmap_featurizer:
    def __init__(
        self,
        mtt_model,
        mlp,
        *,
        extensions_directory=None,
        check_consistency=None,
        device=None,
        length_unit="Angstrom",
        high_scaler=None,
        low_scaler=None,
    ):
        self.featurizer = metatomic_featurizer(
            model=mtt_model,
            extensions_directory=extensions_directory,
            check_consistency=check_consistency,
            device=device,
            length_unit=length_unit,
            output="mtt::aux::energy_last_layer_features",
            use_mead_std=True,
        )

        device = self.featurizer.device
        self._init_sketchmap_mlp(mlp, device, high_scaler, low_scaler)

    def __call__(self, frames, environments):
        features = self.featurizer(frames, environments)

        if features.shape[1] != self.smap_model.input_dim:
            warnings.warn(
                "Dimention of the mtt model features does not correspond to "
                "trained mlp checkpoint. No sketchmap projection is done",
                stacklevel=2,
            )
            return features

        # apply sketchmap
        if self.high_scaler:
            scaled_features = self.high_scaler.transform(features)
            features_tensor = torch.FloatTensor(scaled_features)
        else:
            features_tensor = torch.FloatTensor(features)

        with torch.no_grad():
            X_reduced = self.smap_model(features_tensor)

        X_reduced = X_reduced.detach().numpy()

        if self.low_scaler:
            X_reduced = self.low_scaler.inverse_transform(X_reduced)

        return X_reduced

    def _init_sketchmap_mlp(self, mlp, device, high_scaler=None, low_scaler=None):
        if isinstance(mlp, (str, Path)):
            self._load_mlp_from_checkpoint(mlp, device)
        else:
            self.mlp = mlp.to(device).eval()
            self.high_scaler = high_scaler
            self.low_scaler = low_scaler

    def _load_mlp_from_checkpoint(self, checkpoint_path, device):
        # Check if dependencies were installed
        global torch
        if torch is None:
            try:
                import torch
            except ImportError as e:
                raise ImportError(
                    f"Required package not found: {e}. Please install the "
                    "dependencies with `pip install chemiscope[metatomic]`."
                )

        checkpoint = torch.load(checkpoint_path, device, weights_only=False)

        state_dict = checkpoint["model_state_dict"]
        input_dim = state_dict["fc1.weight"].shape[1]
        output_dim = state_dict["output.weight"].shape[0]

        self.smap_model = self._build_mlp(input_dim, output_dim)
        self.smap_model.load_state_dict(state_dict)
        self.smap_model.eval()

        self.high_scaler = checkpoint["scaler_highdim"]
        self.low_scaler = checkpoint["scaler_lowdim"]

    @staticmethod
    def _build_mlp(input_dim, output_dim):
        from torch import nn

        class MLP(nn.Module):
            def __init__(self, input_dim=1024, output_dim=2):
                super(MLP, self).__init__()

                self.input_dim = input_dim
                self.output_dim = output_dim

                self.fc1 = nn.Linear(input_dim, 512)
                self.fc2 = nn.Linear(512, 256)
                self.fc3 = nn.Linear(256, 128)
                self.output = nn.Linear(128, output_dim)

            def forward(self, x):
                x = torch.relu(self.fc1(x))
                x = torch.relu(self.fc2(x))
                x = torch.relu(self.fc3(x))
                x = self.output(x)
                return x

        return MLP(input_dim, output_dim)
