"""
.. _chemiscope-explore-metatensor:

Using `metatensor`_ models for dataset exploration
==================================================

In this example, we demonstrate how to create and use a `metatensor`_ model with
:py:func:`chemiscope.metatensor_featurizer` to extract features from the model, which
are then displayed using a chemiscope widget. To use this function, some additional
dependencies are required. You can install them with the following command:

.. code:: bash

    pip install chemiscope[metatensor]

.. _metatensor: https://docs.metatensor.org/latest/index.html
"""

# %%
#
# Firstly, we import necessary packages and read structures from the dataset.

from typing import Dict, List, Optional

import ase.io
import torch
from metatensor.torch import Labels, TensorBlock, TensorMap
from metatensor.torch.atomistic import (
    MetatensorAtomisticModel,
    ModelCapabilities,
    ModelMetadata,
    ModelOutput,
    NeighborListOptions,
    System,
)

import chemiscope

frames = ase.io.read("data/explore_c-gap-20u.xyz", ":")

# %%
#
# Using pre-trained models
# ------------------------
#
# Most commonly, you will have an already existing model in metatensor format that
# you'll want to use for dataset exploration. In this case, you'll have to create a
# ``featurizer`` function using :py:func:`chemiscope.metatensor_featurizer`.
#
# ``metatensor_featurizer`` takes an existing model as input. It can be either a
# ``MetatensorAtomisticModel`` instance or a path to a pre-trained model file (here
# ``"model.pt"``)

featurizer = chemiscope.metatensor_featurizer(model="model.pt")

# %%
#
# From here, you can use :py:func:`chemiscope.explore` to visualize the features
# computed from the structures. For this, we are passing the frames, the ``featurizer``
# function, and — as the model computes per-atom properties — environments.

chemiscope.explore(
    frames=frames,
    featurize=featurizer,
    environments=chemiscope.all_atomic_environments(frames),
)

# %%
#
# Defining a custom model
# -----------------------
#
# Let's now move on and see how one can define a fully custom model to use through the
# metatensor interface.
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
# And then, the model will take a PCA of the above features to extract the three most
# relevant dimensions.


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
        if list(outputs.keys()) != ["features"]:
            raise ValueError(
                "this model can only compute 'features', but outputs contains other "
                f"keys: {', '.join(outputs.keys())}"
            )

        if not outputs["features"].per_atom:
            raise NotImplementedError("per structure features are not implemented")

        all_features = []
        all_samples = []

        for system_i, system in enumerate(systems):
            dtype = system.positions.dtype
            device = system.positions.device
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
                [[system_i, atom_i] for atom_i in range(n_atoms)]
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
                values=torch.tensor([[0], [1], [2]]),
            ),
        )
        return {
            "features": TensorMap(
                keys=Labels(names=["_"], values=torch.tensor([[0]])), blocks=[block]
            )
        }


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
        "features": ModelOutput(per_atom=True),
    },
    atomic_types=[6],
    length_unit="angstrom",
    interaction_range=0.0,
    supported_devices=["cpu"],
    dtype="float64",
)

mta_model = MetatensorAtomisticModel(model.eval(), metadata, capabilities)

# %%
#
# For a more detailed example of exporting a model, please check the related
# documentation `page
# <https://docs.metatensor.org/latest/examples/atomistic/1-export-atomistic-model.html>`_
# in `metatensor`_.
#
# Once the model is fully defined, we can use it with
# :py:func:`chemiscope.metatensor_featurizer`:

featurizer = chemiscope.metatensor_featurizer(mta_model, check_consistency=True)
chemiscope.explore(
    frames=frames,
    featurize=featurizer,
    environments=chemiscope.all_atomic_environments(frames),
)

# %%
# The metatensor model can also be easily exported,
# to be shared with collaborators for use in their
# visualization workflows

mta_model.save("model-exported.pt")
# %%
