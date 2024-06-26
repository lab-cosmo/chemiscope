"""
Chemiscope.explore example
==========================
This example demonstrates the utilisation of different methods of
dimensionality reduction and its visualisation using `chemiscope.explore`.

First, import the necessary packages:
"""

# %%
import ase.io
import numpy as np
from dscribe.descriptors import SOAP
from mace.calculators import mace_mp, mace_off
from sklearn.decomposition import KernelPCA
from sklearn.manifold import TSNE
from tqdm.auto import tqdm

import chemiscope

# %%
#
# Example with MACE-OFF and t-SNE
# +++++++++++++++++++++++++++++++++++++
#

# %%
#
# Load the dataset of organic molecules.

qm9_frames = ase.io.read("data/explore_qm9.xyz", ":")

# %%
#
# `chemiscope.explore` expects a `featurize` function that retures reduced data
# to be provided. Let's define a function that computes the descriptors using
# MACE-OFF and then uses t-SNE for the dimensionality reduction.


def mace_off_tsne(frames):
    # At first, we initialize a mace_off calculator:
    descriptor_opt = {"model": "small", "device": "cpu", "default_dtype": "float64"}
    calculator = mace_off(**descriptor_opt)

    # After that, let's calculate the features
    descriptors = []
    for frame in tqdm(frames, disable=True):
        structure_avg = np.mean(
            (calculator.get_descriptors(frame, invariants_only=True)),
            axis=0,
        )
        descriptors.append(structure_avg)
    descriptors = np.array(descriptors)

    # Finally, we apply t-SNE
    perplexity = min(30, descriptors.shape[0] - 1)
    reducer = TSNE(n_components=2, perplexity=perplexity)
    return reducer.fit_transform(descriptors)


# %%
#
# Provide the created featurizer function to `chemiscope.explore`.

cs = chemiscope.explore(qm9_frames, featurize=mace_off_tsne)

# %%
#
# Here we display the visualisation of the pre-computed data for 6k structures
# taken from the `QM9 <https://jla-gardner.github.io/load-atoms/index.html>`_ dataset
# which dimensionality reduction was computed using the previously described algorithm.

chemiscope.show_input("data/mace-off-tsne-qm9.json.gz")

# %%
#
# Example with MACE-MP0 and t-SNE
# ++++++++++++++++++++++++++++++++++++++

# %%
#
# Load a reduced M3CD dataset.

m3cd_frames = ase.io.read("data/explore_m3cd.xyz", ":")


# %%
#
# We are defining a function used in `chemiscope.explore` as a featurizer
# that computes the descriptors using MACE-MP0 and then applies t-SNE.


def mace_mp0_tsne(frames):
    # Initialise a mace-mp0 calculator
    descriptor_opt = {"model": "small", "device": "cpu", "default_dtype": "float64"}
    calculator = mace_mp(**descriptor_opt)

    # Calculate the features
    descriptors = []
    for frame in tqdm(frames, disable=True):
        structure_avg = np.mean(
            (calculator.get_descriptors(frame, invariants_only=True)),
            axis=0,
        )
        descriptors.append(structure_avg)
    descriptors = np.array(descriptors)

    # Apply t-SNE
    perplexity = min(30, descriptors.shape[0] - 1)
    reducer = TSNE(n_components=2, perplexity=perplexity)
    return reducer.fit_transform(descriptors)


# %%
#
# Provide the created function to `chemiscope.explore`.

cs = chemiscope.explore(m3cd_frames, featurize=mace_mp0_tsne)


# %%
#
# Loading reduced data for 1k structures pre-computed using
# the described algorithm.

chemiscope.show_input("data/mace-mp-tsne-m3cd.json.gz")

# %%
#
# Example with SOAP and KPCA
# ++++++++++++++++++++++++++++++++++++++

# %%
#
# Load the dataset.

frames = ase.io.read("data/explore_c-gap-20u.xyz", ":")


# %%
# Define the function that calculates SOAP descriptors and then runs KPCA.


def soap_kpca(frames):
    # Initialise soap calculator
    soap = SOAP(
        species=["C"],
        r_cut=4.5,
        n_max=8,
        l_max=6,
        sigma=0.2,
        rbf="gto",
        average="off",
        periodic=True,
        weighting={"function": "pow", "c": 1, "m": 5, "d": 1, "r0": 3.5},
    )

    # Compute SOAP desciptors
    descriptors = []
    for frame in frames:
        descriptors.append(np.mean(soap.create(frame), axis=0))
    descriptors = np.vstack(descriptors)

    # Apply KPCA
    transformer = KernelPCA(n_components=2, gamma=0.05)
    return transformer.fit_transform(descriptors)


# %%
#
# Provide the created function to `chemiscope.explore`.

cs = chemiscope.explore(frames, featurize=soap_kpca)

# %%
#
# Loading pre-computed dimensionality reduction done using the descibed featurizer
# for the C-GAP-20U dataset.

chemiscope.show_input("data/soap_kpca_c-gap-20u.json.gz")
