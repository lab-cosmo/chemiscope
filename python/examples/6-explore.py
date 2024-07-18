"""
.. _explore-example:

Exploring dataset with chemiscope
=================================

The :py:func:`chemiscope.explore` function provides a streamlined way to visualize
datasets by automatically computing representation and using dimensionality reduction.
This function simplifies the process of dataset exploration by offering a quick overview
through computed properties and dimensionality reduction, allowing to rapidly gain
insights into the composition and structure of data without need to manually implement
and fine-tune the representation process.

This is particularly useful when the specific choice of hyperparameters does not
significantly impact the resulting 2D map. By passing a list of `ase.Atoms
<https://wiki.fysik.dtu.dk/ase/ase/atoms.html>`_ objects (or similar structures from
other libraries) to :py:func:`chemiscope.explore`, it is possible to generate a
chemiscope widget, providing an immediate and intuitive visualization of the dataset.

Additionally, :py:func:`chemiscope.explore` allows to provide a custom function for
representation and dimensionality reduction, offering flexibility for more advanced
usage.

To use this function, some additional dependencies are required. You can install them
with the following command:

.. code:: bash

    pip install chemiscope[explore]

In this example, we will explore several use cases, starting from basic applications to
more customized scenarios.

First, let's import the necessary packages that will be used throughout the examples.
"""

# %%
import os

import ase.io
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
# Basic example
# +++++++++++++
#
# This example shows the basic usage of the :py:func:`chemiscope.explore`. At first,
# read or load the structures from the dataset. Here we use an `ASE package
# <https://wiki.fysik.dtu.dk/ase>`_ to read the structures from the file and have the
# frames as the `ase.Atoms <https://wiki.fysik.dtu.dk/ase/ase/atoms.html>`_ objects.

frames = ase.io.read("data/explore_c-gap-20u.xyz", ":")


# %%
#
# Provide the frames to the :py:func:`chemiscope.explore`. It will generate a Chemiscope
# interactive widget with the reduced dimensionality of data.

chemiscope.explore(frames)

# %%
#
# In this basic case, no featurizer function is provided, so
# :py:func:`chemiscope.explore` uses a default method that applies
# `SOAP (Smooth Overlap of Atomic Positions)
# <https://singroup.github.io/dscribe/latest/tutorials/descriptors/soap.html>`_ to
# compute atomic structure descriptors and then performs `PCA (Principal Component
# Analysis)
# <https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html>`_
# for dimensionality reduction. The resulting components are then added to the
# properties to be used in visualization.


# %%
#
# Example with custom featurizer and custom properties
# ++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# This part illustrates how to create a custom function for dimensionality reduction
# as an argument (``featurize``) to :py:func:`chemiscope.explore`. Inside this function,
# we perform descriptor calculation using `SOAP
# <https://singroup.github.io/dscribe/latest/tutorials/descriptors/soap.html>`_ and
# then reduce the dimensionality with `Kernel PCA
# <https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.KernelPCA.html>`_.
#
# First, let's import the necessary packages.
from dscribe.descriptors import SOAP  # noqa
from sklearn.decomposition import KernelPCA  # noqa

# %%
#
# Define the function ``soap_kpca`` which takes one argument (``frames``). This argument
# contains the structures provided to :py:func:`chemiscope.explore` and is internally
# passed to the ``featurize`` function.


def soap_kpca(frames):
    # Initialise soap calculator. The detailed explanation of the provided
    # hyperparameters can be checked in the documentation of the library (``dscribe``).
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

    # Compute features
    descriptors = soap.create(frames)

    # Apply KPCA
    transformer = KernelPCA(n_components=2, gamma=0.05)

    # Return a 2D array of reduced features
    return transformer.fit_transform(descriptors)


# %%
#
# Provide the created function to :py:func:`chemiscope.explore`.

cs = chemiscope.explore(frames, featurize=soap_kpca)

# %%
#
# We can also provide the additional properties inside, for example, let's extract
# energy from the frames using :py:func:`chemiscope.extract_properties`.

properties = chemiscope.extract_properties(frames, only=["energy"])
cs = chemiscope.explore(frames, featurize=soap_kpca, properties=properties)

# %%
#
# Note: It is possible to add parallelization when computing the SOAP descriptors and
# performing dimensionality reduction with KernelPCA by providing the ``n_jobs``
# parameter. This allows the computation to utilize multiple CPU cores for faster
# processing. An example of how to include ``n_jobs`` is shown below on this page.
#
# To showcase the results of the ``soap_kpca`` function, we have pre-computed it for
# the 6k structures from the `C-GAP-20U
# <https://jla-gardner.github.io/load-atoms/datasets/C-GAP-20U.html>`_ dataset:

fetch_dataset("soap_kpca_c-gap-20u.json.gz")
chemiscope.show_input("data/soap_kpca_c-gap-20u.json.gz")
