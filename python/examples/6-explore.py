"""
.. _explore-example:

Exploring dataset with chemiscope
=================================

The :py:func:`chemiscope.explore` function provides a streamlined way to visualize
datasets as the low dimensional maps. This approach provides a quick and interactive
overview of dataset composition and structure without the need to manually implement and
configure the representation processes. This is particularly useful when the specific
choice of hyperparameters does not significantly impact the resulting low-dimensionality
map.

By passing a list of `ase.Atoms <https://wiki.fysik.dtu.dk/ase/ase/atoms.html>`_ objects
(or similar structures from other libraries) to :py:func:`chemiscope.explore`, it is
possible to generate a chemiscope widget, providing an immediate and intuitive
visualization of the dataset.

By default, the method uses the `PETMADFeaturizer <https://arxiv.org/abs/2506.19674>`_,
which computes representations as the `PET-MAD <https://arxiv.org/abs/2503.14118>`_
features and maps them a low-dimensional MAD dataset latent space.

For more advanced use cases, :py:func:`chemiscope.explore` allows to provide a custom
function for representation and dimensionality reduction.

To use this function, some additional dependencies are required. You can install them
with the following command:

.. code:: bash

    pip install chemiscope[explore]

In this example, we will explore basic and advanced use cases, from simple dataset
visualization to custom featurization.

First, let's import the necessary packages that will be used throughout the examples.
"""

# %%
import ase.io

import chemiscope


# %%
#
# Basic example
# +++++++++++++
#
# This example shows the basic usage of :py:func:`chemiscope.explore`. First, load a
# dataset of structures as `ase.Atoms <https://wiki.fysik.dtu.dk/ase/ase/atoms.html>`_
# objects. Here, we use the samples from the `M3CD dataset
# <https://doi.org/10.24435/materialscloud:rw-t0>`_:

frames = ase.io.read("data/explore_m3cd.xyz", ":")


# %%
#
# Next, pass the frames to :py:func:`chemiscope.explore` to generate an interactive
# Chemiscope. In this basic case, we provide the featurizer version to be used:

chemiscope.explore(frames, featurizer="pet-mad-1.0")

# %%
#
# We can also save the visualization to send it to the colloborators or reopen
# separatelly with :py:func:`chemiscope.read_input`:
chemiscope.explore(frames, featurizer="pet-mad-1.0", write_input="m3cd.chemiscope.json")


# %%
#
# Besides this, it is possible to specify atom-centered environments and properties.
# Environments can be manually defined as a list of tuples in the format
# ``(structure_index, atom_index, cutoff)`` or extracted automatically using
# :py:func:`chemiscope.all_atomic_environments`. We can also configure visualisation
# settings, such as axis and color properties.

properties = chemiscope.extract_properties(frames, only=["energy"])
environments = [(0, 0, 3.5), (1, 0, 3.5), (2, 1, 3.5)]
settings = chemiscope.quick_settings(x="features[1]", y="features[2]", color="energy")
chemiscope.explore(
    frames,
    featurizer="pet-mad-1.0",
    environments=environments,
    properties=properties,
    settings=settings,
)

# %%
#
# Example with custom featurizer
# ++++++++++++++++++++++++++++++
#
# For advanced use cases, you can define a custom featurization function. For example,
# we can describe structures based on their chemical compositions. The function must
# take two arguments: ``frames`` (the input structures) and ``environments`` (optional
# argument for the atom-centered environments). Below, we create a function to calculate
# fractional composition vectors and apply PCA for dimensionality reduction:

import numpy as np  # noqa
from sklearn.decomposition import PCA  # noqa


def fractional_composition_featurize(frames, environments):
    if environments is not None:
        raise ValueError("'environments' are not supported by this featurizer")

    dimentionality = 100

    features = []

    for frame in frames:
        unique, counts = np.unique(frame.numbers, return_counts=True)
        fractions = counts / len(frame.numbers)

        feature_vector = np.zeros(dimentionality)
        for element_number, franction in zip(unique, fractions, strict=True):
            feature_vector[element_number - 1] = franction

        features.append(feature_vector)

    pca = PCA(n_components=3)
    return pca.fit_transform(features)


# %%
#
# Pass the custom featurizer to :py:func:`chemiscope.explore`:

settings = chemiscope.quick_settings(x="features[1]", y="features[2]")
chemiscope.explore(
    frames,
    featurizer=fractional_composition_featurize,
    settings=settings,
)

# %%
#
# For more advanced examples, see the :ref:`next tutorial <advanced-explore-example>`.
