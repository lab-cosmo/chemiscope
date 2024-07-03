import os

from dscribe.descriptors import SOAP
from sklearn.decomposition import PCA

from .jupyter import show


def explore(frames, featurize=None, properties=None, mode="default"):
    """
    Visualize the dataset with dimensionality reduction.
    If no function is provided as a `featurize` argument, a default SOAP and PCA based
    featurizer is used. SOAP parameters (e.g., cutoff radius, number of radial and
    angular functions, etc.) are predefined. The SOAP computation uses all available
    CPU cores for parallelization. PCA reduces the dimensionality to two components
    by default.

    :param list frames: list of ASE Atoms objects

    :param callable featurize: optional. Function to compute features and perform
        dimensionality reduction on the `frames`. The function should take `frames` as
        input and return a 2D array of reduced features. If `None`, a default SOAP and
        PCA based featurizer is used.

    :param dict properties: optional. Additional properties to be included in the
        visualization. This dictionary can contain any other relevant data associated
        with the atomic structures. Properties can be extracted from frames with
        :py:func:`extract_properties` or manually defined by the user.

    :param str mode: optional. Visualization mode for the chemiscope widget. Can be one
        of "default", "structure", or "map". The default mode is "default".

    :return: a chemiscope widget for interactive visualization

    Here is an example of usage with and without providing a function to do the
    reprensentation and reduction. The frames are obtained by reading the structures
    from a file that `ase <ase-io_>`_ can read, and performing PCA using `sklearn`_ on
    a descriptor computed with SOAP with `dscribe`_ library.

    .. code-block:: python

        import chemiscope
        from ase.io import read
        from dscribe.descriptors import SOAP
        from sklearn.decomposition import KernelPCA

        # Read the structures from the dataset
        frames = read("trajectory.xyz", ":")

        # Basic usage with the default featurizer (SOAP + PCA)
        chemiscope.explore(frames)

        # Using a custom featurizer
        chemiscope.explore(frames, featurize=soap_kpca)


        def soap_kpca(frames):
            # Compute descriptors
            soap = SOAP(
                species=["C"],
                r_cut=4.5,
                n_max=8,
                l_max=6,
                periodic=True,
            )
            descriptors = soap.create(frames)

            # Apply KPCA
            transformer = KernelPCA(n_components=2, gamma=0.05)

            # Return a 2D array of reduced features
            return transformer.fit_transform(descriptors)

    .. _ase-io: https://wiki.fysik.dtu.dk/ase/ase/io/io.html
    .. _sklearn: https://scikit-learn.org/
    .. _dscribe: https://singroup.github.io/dscribe/latest/
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
