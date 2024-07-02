from dscribe.descriptors import SOAP
from sklearn.decomposition import PCA

from .jupyter import show


def explore(frames, featurize, properties=None, mode="default"):
    """
    Visualize the dataset with dimensionality reduction.

    :param list frames: list of atomic structures
    Visualize the dataset with dimensionality reduction.

    :param dict properties: optional, dictionary of additional properties.

    :param str mode: optional, widget mode, either ``default``, ``structure`` or ``map``

    :return: chemiscope widget
    """
    # Validate inputs
    if not callable(featurize):
        raise TypeError("'featurize' must be a callable (function)")
    if properties is None:
        properties = {}

    # Apply dimensionality reduction from the provided featurizer
    if featurize is None:
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
    feats = soap.create(frames, n_jobs=-1, only_physical_cores=True)

    # Compute pca
    pca = PCA(n_components=2)
    return pca.fit_transform(feats)
