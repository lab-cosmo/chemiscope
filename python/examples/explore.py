"""
Chemiscope.explore example
==========================
This example illustrates usage of different methods of dimensionality reduction and its
visualisation with `chemiscope.explore`. The primary goal of this function is to enable
quick and effective exploration of a dataset. It follows the following algorithm
internally:

.. image:: /img/chemiscope-explore-algorithm.png
   :align: center

To use this function, some additional dependencies are required. You can install them
with the following command:

.. code:: bash

    pip install chemiscope[explore]

In this example, we will explore several use cases, starting from basic applications to
more customized scenarios.

First, let's import the necessary packages that will be used throughout the examples.
"""

# %%
import ase.io

import chemiscope

# %%
#
# Basic usage
# +++++++++++
#
# This example shows the basic usage of the `chemiscope.explore`. At first, read or
# load the structures from the dataset. Here we use an `ASE package
# <https://wiki.fysik.dtu.dk/ase/>`_ to read the structures from the file and have the
# frames as the `ase.Atoms <https://wiki.fysik.dtu.dk/ase/ase/atoms.html/>`_ objects.

frames = ase.io.read("data/explore_c-gap-20u.xyz", ":")


# %%
#
# Provide the frames to the `chemiscope.explore`. It will generate a Chemiscope
# interactive widget with the reduced dimentionality of data.

chemiscope.explore(frames)

# %%
#
# In this basic case, no featurizer function is provided, so `chemiscope.explore` uses a
# default method that applies `SOAP (Smooth Overlap of Atomic Positions)
# <https://singroup.github.io/dscribe/latest/tutorials/descriptors/soap.html/>`_ to
# compute atomic structure descriptors and then performs `PCA (Principal Component
# Analysis)
# <https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html/>`_
# for dimensionality reduction. The resulting components are then added to the
# properties to be used in visualization.


# %%
#
# Example with SOAP and KPCA
# ++++++++++++++++++++++++++
#
# This part illustrates how to create a custom function for dimensionality reduction
# as an argument (`featurize`) to `chemiscope.explore`. Inside this function, we perform
# descriptor calculation using SOAP and then reduce the dimensionality with `Kernel PCA
# <https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.KernelPCA.html/>`_.
#
# First, let's import the necessary packages.
import os  # noqa

from dscribe.descriptors import SOAP  # noqa
from sklearn.decomposition import KernelPCA  # noqa

# %%
#
# Define the function `soap_kpca` which takes one argument (`frames`). This argument
# contains the structures provided to `chemiscope.explore` and is internally passed
# to the `featurize` function.


def soap_kpca(frames):
    # Initialise soap calculator. The detailed explanation of the provided
    # hyperparameters can be checked in the documentation of the library (`dscribe`).
    soap = SOAP(
        # the dataset used in the example contains only carbon
        species=["C"],
        r_cut=4.5,
        n_max=8,
        l_max=6,
        sigma=0.2,
        rbf="gto",
        average="outer",
        periodic=True,
        weighting={"function": "pow", "c": 1, "m": 5, "d": 1, "r0": 3.5},
    )

    # Get the number of jobs to calculate the features in parallel
    n_jobs = min(len(frames), os.cpu_count())

    # Compute features
    descriptors = soap.create(frames, n_jobs=n_jobs)

    # Apply KPCA. We also use the parallelisation to accelerate the time of computation
    transformer = KernelPCA(n_components=2, gamma=0.05, n_jobs=n_jobs)
    return transformer.fit_transform(descriptors)


# %%
#
# Provide the created function to `chemiscope.explore`.

cs = chemiscope.explore(frames, featurize=soap_kpca)

# %%
#
# Here we are loading pre-computed dimensionality reduction done using the descibed
# `featurize` function for the C-GAP-20U dataset.

chemiscope.show_input("data/soap_kpca_c-gap-20u.json.gz")


# %%
#
# Example with MACE-OFF and t-SNE
# +++++++++++++++++++++++++++++++++++++
#
# In this part, we are going to define another `featurize` function that runs
# calculation of desciptors with MACE-OFF and uses t-SNE for the dimensionality
# reduction.
#
# The dependencies for this example can be installed with the following command:
#
# .. code:: bash
#
#     pip install chemiscope[mace]

# %%
#
# Let's import the necessary libraries.

import numpy as np  # noqa
from mace.calculators import mace_off  # noqa
from sklearn.manifold import TSNE  # noqa

# %%
#
# Load the dataset, in our example we are reading the organic molecules.

qm9_frames = ase.io.read("data/explore_qm9.xyz", ":")

# %%
#
# Now, we are defining a `featurize` function. As on the previous example, it should
# return the reduced data.


def mace_off_tsne(frames):
    # At first, we initialize a mace_off calculator:
    descriptor_opt = {"model": "small", "device": "cpu", "default_dtype": "float64"}
    calculator = mace_off(**descriptor_opt)

    # Calculate MACE features for each frame
    descriptors = []
    for frame in frames:
        structure_avg = np.mean(
            # Only use invariant descriptors (no rotational components)
            (calculator.get_descriptors(frame, invariants_only=True)),
            axis=0,  # Average the descriptors over all atoms in the frame
        )
        descriptors.append(structure_avg)
    descriptors = np.array(descriptors)

    # Get number of jobs for parallelisation
    n_jobs = min(len(frames), os.cpu_count())

    # Apply t-SNE
    perplexity = min(30, descriptors.shape[0] - 1)
    reducer = TSNE(n_components=2, perplexity=perplexity, n_jobs=n_jobs)
    return reducer.fit_transform(descriptors)


# %%
#
# We can also provide the additional properties inside `chemiscope.explore`.
# For example, let's extract diapole moment from the frames using the related helper
# function.

properties = chemiscope.extract_properties(qm9_frames, only=["mu"])

# %%
#
# Provide the created featurizer and the properties to `chemiscope.explore`.

cs = chemiscope.explore(qm9_frames, featurize=mace_off_tsne, properties=properties)

# %%
#
# Here we display the visualisation of the pre-computed data using the described
# function for 6k structures taken from the `QM9
# <https://jla-gardner.github.io/load-atoms/index.html>`_ dataset.

chemiscope.show_input("data/mace-off-tsne-qm9.json.gz")

# %%
#
# Example with MACE-MP0 and t-SNE
# ++++++++++++++++++++++++++++++++++++++
#
# We will define another `featurize` function that uses MACE-MP0
# to calculate the descriptors and t-SNE for the dimensionality
# reduction.
#
# Firstly, import mace library.

from mace.calculators import mace_mp  # noqa

# %%
#
# Load the frames. In this example we are loading the M3CD dataset with the
# reduced number of stuctures.

m3cd_frames = ase.io.read("data/explore_m3cd.xyz", ":")


# %%
#
# We are defining a function used in `chemiscope.explore` as a featurizer
# that computes the descriptors using MACE-MP0 and then applies t-SNE.
# Basically, we repeat the steps done in the previous example but using
# different mace calculator.


def mace_mp0_tsne(frames):
    # Initialise a mace-mp0 calculator
    descriptor_opt = {"model": "small", "device": "cpu", "default_dtype": "float64"}
    calculator = mace_mp(**descriptor_opt)

    # Calculate the features
    descriptors = []
    for frame in frames:
        structure_avg = np.mean(
            (calculator.get_descriptors(frame, invariants_only=True)),
            axis=0,
        )
        descriptors.append(structure_avg)
    descriptors = np.array(descriptors)

    n_jobs = min(len(frames), os.cpu_count())

    # Apply t-SNE
    perplexity = min(30, descriptors.shape[0] - 1)
    reducer = TSNE(n_components=2, perplexity=perplexity, n_jobs=n_jobs)
    return reducer.fit_transform(descriptors)


# %%
#
# Provide the created function to `chemiscope.explore`.

cs = chemiscope.explore(m3cd_frames, featurize=mace_mp0_tsne)

# %%
#
# To show case the result, we are loading pre-computed data using the `mace_mp0_tsne`
# function for 1k structures.

chemiscope.show_input("data/mace-mp-tsne-m3cd.json.gz")
