import os

from dscribe.descriptors import SOAP
from sklearn.decomposition import PCA

from .jupyter import show


def explore(frames, featurize=None, properties=None, mode="default"):
    """
    Visualize the dataset with dimensionality reduction.
    If no function is provided as the `featurize` argument the default one is used
    that performs SOAP and PCA.

    :param list frames: list of atomic structures

    :param dict properties: optional, dictionary of additional properties.

    :param str mode: optional, widget mode, either ``default``, ``structure`` or ``map``

    :return: chemiscope widget
    """
    # Validate inputs
    if featurize is not None and not callable(featurize):
        raise TypeError("'featurize' must be a callable (function)")
    if properties is None:
        properties = {}

    # Apply dimensionality reduction from the provided featurizer
    if featurize is not None:
        X_reduced = featurize(frames)

    # Use default featurizer
    else:
        X_reduced = soap_pca(frames)

    # Add dimensionality reduction results to properties
    settings = {"map": {}}
    axis = ["x", "y", "z"]
    for i in range(X_reduced.shape[1]):
        properties[f"Component {i + 1}"] = X_reduced[:, i].tolist()
        if i < len(axis):
            settings["map"][axis[i]] = {"property": f"Component {i + 1}"}

    # Return chemiscope widget
    return show(frames=frames, properties=properties, mode=mode, settings=settings)


def soap_pca(frames):
    """
    Computes SOAP features for a given set of atomic structures and performs
    dimensionality reduction using PCA.

    Note:
    - The SOAP descriptor parameters such as `r_cut`, `n_max`, `l_max`, `sigma`, `rbf`,
    `average`, `periodic`, `weighting`, and `compression` are pre-defined within
    function.
    - The function utilizes all available CPU cores for parallel computation of SOAP
    descriptors.
    """
    # Get global species
    species = set()
    for frame in frames:
        species.update(frame.get_chemical_symbols())
    species = list(species)

    # Check if periodic
    is_periodic = all(all(frame.get_pbc()) for frame in frames)

    # Initialize calculator
    soap = SOAP(
        species=species,
        r_cut=4.5,
        n_max=8,
        l_max=6,
        sigma=0.2,
        rbf="gto",
        average="outer",
        periodic=is_periodic,
        weighting={"function": "pow", "c": 1, "m": 5, "d": 1, "r0": 3.5},
        compression={"mode": "mu1nu1"},
    )

    # Calculate descriptors
    n_jobs = min(len(frames), os.cpu_count())
    feats = soap.create(frames, n_jobs=n_jobs)

    # Compute pca
    pca = PCA(n_components=2)
    return pca.fit_transform(feats)
