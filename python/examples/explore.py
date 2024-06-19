"""
chemiscope.explore example
==========================

"""

# %%
import time

import ase.io
import numpy as np
from sklearn.decomposition import PCA

import chemiscope

# from mace.calculators import mace_off
# from tqdm.auto import tqdm


# %%
# Load dataset
frames = ase.io.read("data/trajectory.xyz", ":")


# %%
# Function to calculate the features. For example, here we use MACE OFF
def get_descriptors(frames):
    """
    For now I have an issue using mace_off in this project, should be related to its
    installation that I don't do correctly. So, I run the code below in the notebook

    descriptor_opt = {"model": "small", "device": "cpu", "default_dtype": "float64"}
    calculator = mace_off(**descriptor_opt)
    descriptors = []
    for frame in tqdm(frames):
        structure_avg = np.mean(
            (calculator.get_descriptors(frame, invariants_only=True)),
            axis=0,
        )
        descriptors.append(structure_avg)
    np.save("data/trajectory-mace_off_features.npy", np.array(descriptors))
    """

    return np.load("data/trajectory-mace_off_features.npy")


# %%
# Function to provide to the .explore. It computes the dimentionality reduction
# for the frames. Here we use PCA as the technic


def mace_off_pca(frames):
    descriptors = get_descriptors(frames)

    start_time = time.time()
    reducer = PCA(n_components=2)
    PCA_X_reduced = reducer.fit_transform(descriptors)
    execution_time = time.time() - start_time

    print(f"PCA execution time: {execution_time:.3f} seconds")
    return PCA_X_reduced


# %%
# Example usage of chemiscope.explore to visualize data using PCA

chemiscope.explore(frames, reducer=mace_off_pca, mode="default")
