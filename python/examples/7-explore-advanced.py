"""
Advanced dataset exploration
============================

The :ref:`previous example <explore-example>` introduced :py:func:`chemiscope.explore`
and how to use it for automatic exploration of dataset. In this example, we'll show some
more representations and give additional featurizers you can use in your own code.
"""

# %%
import os

import ase.io
import numpy as np
import requests

import chemiscope


def fetch_dataset(filename, base_url="https://zenodo.org/records/12748925/files/"):
    """Helper function to load the pre-computed examples"""
    local_path = "data/" + filename
    if not os.path.isfile(local_path):
        response = requests.get(base_url + filename)
        with open(local_path, "wb") as file:
            file.write(response.content)


# %%
#
# Example with MACE-OFF and t-SNE
# +++++++++++++++++++++++++++++++++++++
#
# In this part, we are going to define another ``featurize`` function that runs
# calculation of desciptors with `MACE-OFF <https://github.com/ACEsuit/mace>`_ and uses
# `t-SNE
# <https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html>`_ for
# the dimensionality reduction.
#
# The dependencies for this example can be installed with the following command:
#
# .. code:: bash
#
#     pip install mace-torch scikit-learn

# %%
#
# Let's import the necessary libraries.

from mace.calculators import mace_off  # noqa
from sklearn.manifold import TSNE  # noqa

# %%
#
# Load the dataset, in our example we are reading the organic molecules.

qm9_frames = ase.io.read("data/explore_qm9.xyz", ":")

# %%
#
# Now, we are defining a ``featurize`` function. As on the previous example, it should
# return the reduced data.


def mace_off_tsne(frames, environments):
    if environments is not None:
        raise ValueError("'environments' are not supported")
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
# We can also extract the additional properties, for example, dipole moment.

properties = chemiscope.extract_properties(qm9_frames, only=["mu"])

# %%
#
# Provide the created featurizer and the properties to :py:func:`chemiscope.explore`.

cs = chemiscope.explore(qm9_frames, featurize=mace_off_tsne, properties=properties)

# %%
#
# Here we display the visualisation of the pre-computed data using the described
# function for 6k structures taken from the `QM9
# <https://jla-gardner.github.io/load-atoms/index.html>`_ dataset.
# The map is zoomed in to highlight a cluster of zwitterions grouped together
# by application of the previously defined function with `MACE-OFF
# <https://github.com/ACEsuit/mace>`_ and `t-SNE
# <https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html>`_.

fetch_dataset("mace-off-tsne-qm9.json.gz")
chemiscope.show_input("data/mace-off-tsne-qm9.json.gz")

# %%
#
# Example with MACE-MP0 and t-SNE
# ++++++++++++++++++++++++++++++++++++++
#
# We will define another ``featurize`` function that uses MACE-MP0
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
# We are defining a function used in :py:func:`chemiscope.explore` as a featurizer
# that computes the descriptors using MACE-MP0 and then applies t-SNE.
# Basically, we repeat the steps done in the previous example but using
# different mace calculator.


def mace_mp0_tsne(frames, environments):
    if environments is not None:
        raise ValueError("'environments' are not supported")
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
# Provide the created function to :py:func:`chemiscope.explore`.

cs = chemiscope.explore(m3cd_frames, featurize=mace_mp0_tsne)

# %%
#
# To show case the result, we are loading pre-computed data using the ``mace_mp0_tsne``
# function for 1k structures.


fetch_dataset("mace-mp-tsne-m3cd.json.gz")
chemiscope.show_input("data/mace-mp-tsne-m3cd.json.gz")
