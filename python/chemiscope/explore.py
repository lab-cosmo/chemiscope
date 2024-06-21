import numpy as np
from dscribe.descriptors import SOAP
from sklearn.decomposition import KernelPCA

from .jupyter import show

# import tqdm
# from mace.calculators import mace_off
# from sklearn.manifold import TSNE


def explore(frames, properties=None, featurize=None, mode="default"):
    """
    Questions:
    - Which dimensionality reduction algorithm to use by default?

    Use cases:

    1. Basic
       chemiscope.explore(frames)
       Result: default dim. reduction algorithm used; displays only reduced data

    2. With custom featurize function
       chemiscope.explore(frames, featurize=get_mace_pca)
       Result: calls 'get_mace_pca' for reduction and displays the result

    3. With custom properties
       chemiscope.explore(frames, properties) or
       chemiscope.explore(frames, properties, get_mace_pca) etc
       Result: uses default or custom reducer, displays reduced data + user properties

    Flow:
    0. Choose appropriate hypers for computing descriptors (default featurizer?)
    1. Compute descriptors
    2. Perform dimensionality reduction on these descriptors
    3. Create the widget
    """
    # Validate inputs
    if featurize is not None and not callable(featurize):
        raise TypeError("'featurize' must be a callable (function)")
    if properties is None:
        properties = {}

    # No featurizer -> apply the default dimensionality reduction
    if featurize is None:
        X_reduced = apply_soap_kpca(frames)

    # Featurizer provided -> call it
    else:
        X_reduced = featurize(frames)

    # Add dimensionality reduction results to properties
    properties["Component 1"] = X_reduced[:, 0].tolist()
    properties["Component 2"] = X_reduced[:, 1].tolist()

    # Return chemiscope widget
    return show(frames=frames, properties=properties, mode=mode)


def apply_soap_kpca(frames):
    # Calculate SOAP descriptors
    descriptors = get_soap_descriptors(frames)

    # Reshape SOAP descriptors for KPCA
    n_samples, n_atoms, n_features = descriptors.shape
    descriptors = descriptors.reshape(n_samples, n_atoms * n_features)

    # Apply KPCA
    return get_kpca(descriptors)


def get_soap_descriptors(frames):
    # Get species from frames
    species = set()
    for frame in frames:
        species.update(frame.get_chemical_symbols())
    species = list(species)

    # Create soap calculator
    soap = SOAP(species=species, r_cut=6.0, n_max=8, l_max=6)

    # Compute the desciptors
    descriptors = []
    for frame in frames:
        descriptors.append(soap.create(frame))
    return np.array(descriptors)


def get_kpca(descriptors):
    GAMMA = 0.01
    transformer = KernelPCA(n_components=2, kernel="rbf", gamma=GAMMA)
    return transformer.fit_transform(descriptors)


"""
MACE-OFF + t-SNE

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

"""
