"""
.. _explore-advanced-example:

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
from urllib3.util.retry import Retry

import chemiscope


def fetch_dataset(filename, base_url, local_path=""):
    """Helper function to load the pre-computed examples"""

    local_file = local_path + filename
    if os.path.isfile(local_file):
        return

    # Retry strategy: wait 1s, 2s, 4s, 8s, 16s on 429/5xx errors
    retry_strategy = Retry(
        total=5, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504]
    )
    session = requests.Session()
    session.mount("https://", requests.adapters.HTTPAdapter(max_retries=retry_strategy))

    # Fetch with automatic retry and error raising
    response = session.get(base_url + filename)
    response.raise_for_status()

    with open(local_file, "wb") as file:
        file.write(response.content)


fetch_dataset(
    "mace-off-tsne-qm9.json.gz",
    "https://huggingface.co/datasets/lab-cosmo/chemiscope-visualization/resolve/main/",
    "data/",
)

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

qm9_structures = ase.io.read("data/explore_qm9.xyz", ":")

# %%
#
# Now, we are defining a ``featurize`` function. As on the previous example, it should
# take ``structures`` and ``environments`` as the inputs and return an array of
# features.


def mace_off_tsne(structures, environments):
    if environments is not None:
        raise ValueError("'environments' are not supported by this featurizer")

    # At first, we initialize a mace_off calculator:
    descriptor_opt = {"model": "small", "device": "cpu", "default_dtype": "float64"}
    calculator = mace_off(**descriptor_opt)

    # Calculate MACE features for each structures
    descriptors = []
    for structure in structures:
        structure_avg = np.mean(
            # Only use invariant descriptors (no rotational components)
            (calculator.get_descriptors(structure, invariants_only=True)),
            axis=0,  # Average the descriptors over all atoms in the structure
        )
        descriptors.append(structure_avg)
    descriptors = np.array(descriptors)

    # Get number of jobs for parallelisation
    n_jobs = min(len(structures), os.cpu_count())

    # Apply t-SNE
    perplexity = min(30, descriptors.shape[0] - 1)
    reducer = TSNE(n_components=2, perplexity=perplexity, n_jobs=n_jobs)
    return reducer.fit_transform(descriptors)


# %%
#
# We can also extract the additional properties, for example, dipole moment.

properties = chemiscope.extract_properties(qm9_structures, only=["mu"])

# %%
#
# Provide the created featurizer and the properties to :py:func:`chemiscope.explore`.

cs = chemiscope.explore(qm9_structures, featurizer=mace_off_tsne, properties=properties)

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
# ``environments`` parameter specifying which atoms in the structures are used for
# calculating the descriptors.
#
# Firstly, import mace library.

from mace.calculators import mace_mp  # noqa

# %%
#
# Load the structures. In this example we are loading the reduced MC3D dataset.

mc3d_structures = ase.io.read("data/explore_mc3d.xyz", ":")


# %%
#
# We are defining a featurizer function by basically repeating the steps from the
# previous example but using different MACE calculator.


def mace_mp0_tsne(structures, environments):
    # Initialize a mace-mp0 calculator
    descriptor_opt = {"model": "small", "device": "cpu", "default_dtype": "float64"}
    calculator = mace_mp(**descriptor_opt)

    # Calculate the features
    if environments is None:
        descriptors = []
        for structure in structures:
            structure_avg = np.mean(
                (calculator.get_descriptors(structure, invariants_only=True)),
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
            atoms = structures[structure_index]
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
    mc3d_structures, featurizer=mace_mp0_tsne, environments=[(1, 2, 3.5), (2, 0, 3.5)]
)
