"""
.. _advanced-explore-example:

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


fetch_dataset("mace-off-tsne-qm9.json.gz")

# %%
#
# Example with MACE-OFF and t-SNE
# +++++++++++++++++++++++++++++++
#
# In this part, we are going to define ``featurize`` function that calculates desciptors
# with `MACE-OFF <https://github.com/ACEsuit/mace>`_ and uses `t-SNE
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
# take ``frames`` and ``environments`` as the inputs and return an array of features.


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

cs = chemiscope.explore(qm9_frames, featurizer=mace_off_tsne, properties=properties)

# %%
#
# Here we display the visualization of the pre-computed data using the described
# function for 6k structures taken from the `QM9
# <https://jla-gardner.github.io/load-atoms/index.html>`_ dataset. The map is zoomed in
# to highlight a cluster of zwitterions grouped together by running the previously
# defined ``mace_off_tsne`` featurizer.

chemiscope.show_input("data/mace-off-tsne-qm9.json.gz")

# %%
#
# Example with MACE-MP0, t-SNE and environments
# +++++++++++++++++++++++++++++++++++++++++++++
#
# This example demonstrates how to compute descriptors using the MACE-MP0 and t-SNE with
# ``environments`` parameter specifying which atoms in the frames are used for
# calculating the descriptors.
#
# Firstly, import mace library.

from mace.calculators import mace_mp  # noqa

# %%
#
# Load the frames. In this example we are loading the reduced M3CD dataset.

m3cd_frames = ase.io.read("data/explore_m3cd.xyz", ":")


# %%
#
# We are defining a featurizer function by basically repeating the steps from the
# previous example but using different MACE calculator.


def mace_mp0_tsne(frames, environments):
    # Initialize a mace-mp0 calculator
    descriptor_opt = {"model": "small", "device": "cpu", "default_dtype": "float64"}
    calculator = mace_mp(**descriptor_opt)

    # Calculate the features
    if environments is None:
        descriptors = []
        for frame in frames:
            structure_avg = np.mean(
                (calculator.get_descriptors(frame, invariants_only=True)),
                axis=0,
            )
            descriptors.append(structure_avg)
    else:
        grouped_envs = {}
        unique_structures = set()

        # Group atom indices from environments
        for structure_index, atom_index, _cutoff in environments:
            if structure_index not in grouped_envs:
                grouped_envs[structure_index] = []
            grouped_envs[structure_index].append(atom_index)
            unique_structures.add(structure_index)

        # Compute descriptors per specified atom
        descriptors = []
        for structure_index in sorted(grouped_envs):
            atoms = frames[structure_index]
            atom_indices = grouped_envs[structure_index]

            full_descriptors = calculator.get_descriptors(atoms, invariants_only=True)
            for atom_index in atom_indices:
                descriptors.append(full_descriptors[atom_index])

    descriptors = np.array(descriptors)

    n_jobs = min(len(descriptors), os.cpu_count())

    # Apply t-SNE
    perplexity = min(30, descriptors.shape[0] - 1)
    reducer = TSNE(n_components=2, perplexity=perplexity, n_jobs=n_jobs)
    return reducer.fit_transform(descriptors)


# %%
#
# Provide a created function and environments to :py:func:`chemiscope.explore`. The
# environments are manually defined following the format ``[index of structure, index
# of atom, cutoff]``.

chemiscope.explore(
    m3cd_frames, featurizer=mace_mp0_tsne, environments=[(1, 2, 3.5), (2, 0, 3.5)]
)
