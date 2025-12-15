"""
.. _chemiscope-explore-metatomic:

Using metatomic models for dataset exploration
==============================================

In this example, we demonstrate how to create and use a `metatomic`_ model with
:py:func:`chemiscope.metatomic_featurizer` to extract features from the model, which are
then displayed using a chemiscope widget. To use this function, some additional
dependencies are required. You can install them with the following command:

.. code:: bash

    pip install chemiscope[explore]

.. _metatomic: https://docs.metatensor.org/metatomic/
"""

# %%
#
# Firstly, we import necessary packages and read structures from the dataset.

import ase.io

import chemiscope


frames = ase.io.read("data/explore_c-gap-20u.xyz", ":")

# %%
#
# Using pre-trained models
# ------------------------
#
# Most commonly, you will have an already existing model in metatensor format that
# you'll want to use for dataset exploration. In this case, you'll have to create a
# ``featurizer`` function using :py:func:`chemiscope.metatomic_featurizer`.
#
# ``metatomic_featurizer`` takes an existing model as input. It can be either a
# ``AtomisticModel`` instance or a path to a pre-trained model file (here
# ``"model.pt"``)

featurizer = chemiscope.metatomic_featurizer(model="model.pt")

# %%
#
# From here, you can use :py:func:`chemiscope.explore` to visualize the features
# computed from the structures. For this, we are passing the frames, the ``featurizer``
# function, and — as the model computes per-atom properties — environments.

chemiscope.explore(
    frames=frames,
    featurizer=featurizer,
    environments=chemiscope.all_atomic_environments(frames),
)

# %%
#
# Defining a custom model
# -----------------------
#
# Let's now move on and see how one can define a fully custom model to use through the
# metatomic interface.
#
# Here we will use an atom-centered representation, where each atomic environment is
# represented with the moments of the positions of the neighbors up to a maximal order.
#
# The model computes moments up to a specified maximum order :math:`k_{\text{max}}`,
# computing a representation :math:`F_i^k`
#
# .. math::
#
#     F_i^k = \sum_{j} \frac{r_{ij}}{r_c}^k
#
# where :math:`r_{ij}` is the distance between atom :math:`i` and its neighbor
# :math:`j`, :math:`k` is the moment order and :math:`r_c` is the cutoff radius.
#
# Having computed these moments, the model will take a PCA of their values to
# extract the three most relevant dimensions.

from typing import Dict, List, Optional  # noqa: E402

import torch  # noqa: E402
from metatensor.torch import Labels, TensorBlock, TensorMap  # noqa: E402
from metatomic.torch import (  # noqa: E402
    AtomisticModel,
    ModelCapabilities,
    ModelMetadata,
    ModelOutput,
    NeighborListOptions,
    System,
)


class FeatureModel(torch.nn.Module):
    def __init__(self, cutoff: float, max_k: int):
        super().__init__()
        self.cutoff = cutoff
        self.max_k = max_k

        self._neighbors_options = NeighborListOptions(
            cutoff=cutoff, full_list=True, strict=True
        )

    def requested_neighbor_lists(self) -> List[NeighborListOptions]:
        # our model requires a neighbor list, that will be computed and provided to it
        # automatically.
        return [self._neighbors_options]

    def forward(
        self,
        systems: List[System],
        outputs: Dict[str, ModelOutput],
        selected_atoms: Optional[Labels] = None,
    ) -> Dict[str, TensorMap]:
        results: Dict[str, TensorMap] = {}

        for name, output in outputs.items():
            # feature variant 1: moments + PCA
            if name == "features":
                results["features"] = self._compute_moment_features(systems, output)

            # feature variant 2: cos/sin features of positions
            elif name == "features/cos_sin":
                results["features/cos_sin"] = self._compute_cos_features(
                    systems, output
                )

            else:
                raise ValueError(f"Unknown output '{name}' requested")

        return results

    def _compute_moment_features(
        self, systems: List[System], output: ModelOutput
    ) -> TensorMap:
        if not output.per_atom:
            raise NotImplementedError(
                "per structure features are not implemented for variant 'features'"
            )

        all_features = []
        all_samples = []

        dtype = systems[0].positions.dtype
        device = systems[0].positions.device

        for system_i, system in enumerate(systems):
            n_atoms = len(system.positions)

            # Initialize a tensor to store features for each atom
            features = torch.zeros((n_atoms, self.max_k), dtype=dtype, device=device)

            # get the neighbor list for this system
            neighbors = system.get_neighbor_list(self._neighbors_options)
            i = neighbors.samples.column("first_atom")

            r_ij = torch.linalg.vector_norm(neighbors.values.reshape(-1, 3), dim=1)
            r_ij /= self.cutoff

            for k in range(self.max_k):
                features[i, k] += torch.pow(r_ij, k)

            all_features.append(features)

            # Create labels for each atom in the system
            system_atom_labels = torch.tensor(
                [[system_i, atom_i] for atom_i in range(n_atoms)], device=device
            )
            all_samples.append(system_atom_labels)

        # Concatenate features and labels across all systems
        features_tensor = torch.cat(all_features, dim=0)
        samples_tensor = torch.cat(all_samples, dim=0)

        # Take the PCA of the features
        _, _, V = torch.linalg.svd(features_tensor - features_tensor.mean())
        features_pca = features_tensor @ V[:3].T

        # Add metadata to the output
        block = TensorBlock(
            values=features_pca,
            samples=Labels(names=["system", "atom"], values=samples_tensor),
            components=[],
            properties=Labels(
                names=["feature"],
                values=torch.tensor([[0], [1], [2]], device=device),
            ),
        )
        return TensorMap(
            keys=Labels(names=["_"], values=torch.tensor([[0]], device=device)),
            blocks=[block],
        )

    def _compute_cos_features(
        self, systems: List[System], output: ModelOutput
    ) -> TensorMap:
        all_features = []

        device = systems[0].positions.device

        for system in systems:
            # Compute the norm of each atomic position
            features = system.positions.norm(p=2, dim=1).unsqueeze(1).repeat(1, 2)

            features[:, 0] = torch.cos(features[:, 0])
            features[:, 1] = torch.sin(features[:, 1])

            if not output.per_atom:
                features = features.sum(dim=0, keepdim=True)

            all_features.append(features)

        if output.per_atom:
            samples_list: List[List[int]] = []
            for s, system in enumerate(systems):
                for a in range(len(system)):
                    samples_list.append([s, a])

            samples = Labels(
                ["system", "atom"],
                torch.tensor(samples_list, device=device),
            )
        else:
            samples = Labels(
                ["system"], torch.arange(len(systems), device=device).unsqueeze(1)
            )

        # Add metadata to the output
        block = TensorBlock(
            values=torch.cat(all_features, dim=0),
            samples=samples.to(device),
            components=[],
            properties=Labels(
                ["cos", "sin"], torch.tensor([[1, 0], [0, 1]], device=device)
            ),
        )

        return TensorMap(
            keys=Labels(["_"], torch.tensor([[0]], device=device)),
            blocks=[block],
        )


# %%
#
# With the class defined, we can now create an instance of the model, giving ``cutoff``
# and ``max_k`` as a maximal moment to compute. We don’t need to train this model since
# there are no trainable parameters inside.

model = FeatureModel(cutoff=4.5, max_k=6)


# %%
#
# Next, we set up the model metadata and capabilities:

metadata = ModelMetadata(
    name="Example moment model",
    description=(
        "A model that computes atom-centered features based on the distances of "
        "neighboring atoms"
    ),
)

capabilities = ModelCapabilities(
    outputs={
        "features": ModelOutput(
            per_atom=True, description="PCA of neighbor distance moments"
        ),
        "features/cos_sin": ModelOutput(
            per_atom=True, description="Cosine and sin of atomic position norms"
        ),
    },
    atomic_types=[6],
    length_unit="angstrom",
    interaction_range=0.0,
    supported_devices=["cpu"],
    dtype="float64",
)

mta_model = AtomisticModel(model.eval(), metadata, capabilities)

# %%
#
# For a more detailed example of exporting a model, please check the related
# `tutorial
# <https://docs.metatensor.org/metatomic/latest/examples/1-export-atomistic-model.html>`_
# in metatomic documentation.
#
# Once the model is fully defined, we can use it with
# :py:func:`chemiscope.metatomic_featurizer`:

featurizer = chemiscope.metatomic_featurizer(mta_model, check_consistency=True)
chemiscope.explore(
    frames=frames,
    featurizer=featurizer,
    environments=chemiscope.all_atomic_environments(frames),
)

# %%
#
# The metatomic model can also be easily exported, to be shared with collaborators for
# use in their visualization workflows

mta_model.save("model-exported.pt")

# %%
#
# Using variants of features
# --------------------------
#
# The model we defined provides multiple feature outputs (variants). By default, the
# ``metatomic_featurizer`` uses the main ``"features"`` output (the moments + PCA). You
# can select a variant by passing the ``variant`` argument when creating the featurizer.
#
# For example, to use the ``"features/cos_sin"`` variant:

cos_sin_featurizer = chemiscope.metatomic_featurizer(
    mta_model, variant="cos_sin", check_consistency=True
)

chemiscope.explore(
    frames=frames,
    featurizer=cos_sin_featurizer,
    environments=chemiscope.all_atomic_environments(frames),
)
