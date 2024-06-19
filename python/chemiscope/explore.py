import numpy as np
import tqdm
from mace.calculators import mace_off
from sklearn.manifold import TSNE

from .jupyter import show


def explore(frames, properties=None, reducer=None, mode="default"):
    """
    Questions:
    - Which dimensionality reduction algorithm to use by default?

    Use cases:

    1. Basic
       chemiscope.explore(frames)
       Result: default dim. reduction algorithm used; displays only reduced data

    2. With custom reducer function
       chemiscope.explore(frames, reducer=get_mace_pca)
       Result: calls 'get_mace_pca' for reduction and displays the result

    3. With custom properties
       chemiscope.explore(frames, properties) or
       chemiscope.explore(frames, properties, get_mace_pca) etc
       Result: uses default or custom reducer, displays reduced data + user properties

    Flow:
    0. Choose appropriate hypers for computing descriptors (default reducer?)
    1. Compute descriptors
    2. Perform dimensionality reduction on these descriptors
    3. Create the widget
    """
    # Validate inputs
    if reducer is not None and not callable(reducer):
        raise TypeError("'reducer' must be a callable (function)")
    if properties is None:
        properties = {}

    # No reducer -> apply the default dimensionality reduction
    if reducer is None:
        X_reduced = apply_mace_tsne(frames)

    # Reducer provided -> call it
    else:
        X_reduced = reducer(frames)

    # Add dimensionality reduction results to properties
    properties["Component 1"] = X_reduced[:, 0].tolist()
    properties["Component 2"] = X_reduced[:, 1].tolist()

    # Return chemiscope widget
    return show(frames=frames, properties=properties, mode=mode)


def apply_mace_tsne(frames):
    # Initialize MACE OFF calculator
    descriptor_opt = {"model": "small", "device": "cpu", "default_dtype": "float64"}
    calculator_mace_off = mace_off(**descriptor_opt)

    # Compute features
    descriptors = compute_mace_features(frames, calculator_mace_off)

    # Perform dimentionality reduction using TSNE
    return apply_tsne(descriptors)


def compute_mace_features(frames, calculator, invariants_only=True):
    descriptors = []
    for frame in tqdm(frames):
        # Calculate average structure descriptor for the frame
        structure_avg = np.mean(
            (calculator.get_descriptors(frame, invariants_only=invariants_only)),
            axis=0,
        )
        descriptors.append(structure_avg)
    return np.array(descriptors)


def apply_tsne(descriptors):
    # Initialize reducer
    perplexity = min(30, descriptors.shape[0] - 1)
    reducer = TSNE(n_components=2, perplexity=perplexity)

    # Perform dimensionality reduction using t-SNE
    return reducer.fit_transform(descriptors)
