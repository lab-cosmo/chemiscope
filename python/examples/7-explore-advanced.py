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

from dscribe.descriptors import SOAP  # noqa
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
        raise ValueError("'environments' are not supported by this featurizer")

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
        raise ValueError("'environments' are not supported by this featurizer")

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

# %%
#
# Example with SOAP, t-SNE and environments
# +++++++++++++++++++++++++++++++++++++++++++++++++
#
# This example demonstrates how to compute descriptors using the SOAP and t-SNE with
# ``environments`` parameter specifying which atoms in the frames are used for
# calculating the descriptors.
#
# We are defining a custom featurizer that takes frames and environments.


def soap_tnse_with_environments(frames, environments):
    if environments is None:
        raise ValueError("'environments' must be provided")

    grouped_envs = {}
    unique_structures = set()

    # Get atom-centered indices from environments
    for [env_index, atom_index, _cutoff] in environments:
        if env_index not in grouped_envs:
            grouped_envs[env_index] = []
        grouped_envs[env_index].append(atom_index)
        unique_structures.add(env_index)
    centers = list(grouped_envs.values())

    # only include frames that are present in the environments
    if len(unique_structures) != len(frames):
        frames = [frames[index] for index in sorted(unique_structures)]

    # Get global species
    species = set()
    for frame in frames:
        species.update(frame.get_chemical_symbols())
    species = list(species)

    # Initialize calculator
    soap = SOAP(
        species=species,
        r_cut=4.5,
        n_max=8,
        l_max=6,
        sigma=0.2,
        rbf="gto",
        average="outer",
        periodic=True,
        weighting={"function": "pow", "c": 1, "m": 5, "d": 1, "r0": 3.5},
        compression={"mode": "mu1nu1"},
    )

    # Calculate descriptors
    feats = soap.create(frames, centers=centers)

    # Compute tsne
    perplexity = min(30, feats.shape[0] - 1)
    reducer = TSNE(n_components=2, perplexity=perplexity)
    return reducer.fit_transform(feats)


# %%
#
# Provide a created function and environments to :py:func:`chemiscope.explore`. The
# environments are manually defined following the format ``[index of structure, index
# of atom, cutoff]``.

chemiscope.explore(
    frames=m3cd_frames,
    featurize=soap_tnse_with_environments,
    environments=[(1, 2, 3.5), (2, 0, 3.5)],
)
