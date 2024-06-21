"""
Chemiscope.explore example
==========================
"""

# %%
import time

import ase.io
import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

import chemiscope

# %%
# Load dataset
frames = ase.io.read("data/trajectory.xyz", ":")


# %%
# Function to calculate features using MACE OFF
def get_mace_off_descriptors(frames):
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

    # Return pre-computed descriptors
    return np.load("data/trajectory-mace_off_features.npy")


# %%
# Function to provide to `chemiscope.explore` for PCA
def mace_off_pca(frames):
    descriptors = get_mace_off_descriptors(frames)

    start_time = time.time()
    reducer = PCA(n_components=2)
    PCA_X_reduced = reducer.fit_transform(descriptors)
    execution_time = time.time() - start_time

    print(f"PCA execution time: {execution_time:.3f} seconds")
    return PCA_X_reduced


# %%
# Function to provide to `chemiscope.explore` for TSNE
def mace_off_tsne(frames):
    descriptors = get_mace_off_descriptors(frames)

    start_time = time.time()
    perplexity = min(30, descriptors.shape[0] - 1)
    reducer = TSNE(n_components=2, perplexity=perplexity)
    X_reduced = reducer.fit_transform(descriptors)
    execution_time = time.time() - start_time

    print(f"TSNE execution time: {execution_time:.3f} seconds")
    return X_reduced


# %%
# Get properties from the frames
properties = chemiscope.extract_properties(frames)


# %%
# Example of `chemiscope.explore` with the default featurizer (SOAP + KPCA)

chemiscope.explore(frames, properties=properties)

# %%
# Example of `chemiscope.explore` with MACE OFF + PCA
chemiscope.explore(frames, featurize=mace_off_pca, properties=properties)

# %%
# Example of `chemiscope.explore` with MACE OFF + TSNE

chemiscope.explore(frames, featurize=mace_off_tsne, properties=properties)
