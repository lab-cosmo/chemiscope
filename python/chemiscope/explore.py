import os

from .jupyter import show


def explore(frames, featurize=None, properties=None, environments=None, mode="default"):
    """
    Automatically explore a dataset containing all structures in ``frames``.

    This function computes some low-dimensionality representation of the frames, and
    uses chemiscope to visualize both the resulting embedding and the structures
    simultaneously. ``featurize`` can be used to specify a custom representation and/or
    dimensionality reduction method.

    If no function is provided as a ``featurize`` argument, a default SOAP and PCA based
    featurizer is used. SOAP parameters (e.g., cutoff radius, number of radial and
    angular functions, etc.) are predefined. The SOAP computation uses all available CPU
    cores for parallelization. PCA reduces the dimensionality to two components by
    default.

    :param list frames: list of frames

    :param callable featurize: optional. Function to compute features and perform
        dimensionality reduction on the ``frames``. The function should take ``frames``
        as input and return a 2D array of reduced features. If ``None``, a default SOAP
        and PCA based featurizer is used.

    :param dict properties: optional. Additional properties to be included in the
        visualization. This dictionary can contain any other relevant data associated
        with the atomic structures. Properties can be extracted from frames with
        :py:func:`extract_properties` or manually defined by the user.

    :param environments: optional. List of environments (described as
        ``(structure id, center id, cutoff)``) to include when extracting the
        atomic properties. Can be extracted from frames with
        :py:func:`all_atomic_environments`.
        or manually defined.

    :param str mode: optional. Visualization mode for the chemiscope widget. Can be one
        of "default", "structure", or "map". The default mode is "default".

    :return: a chemiscope widget for interactive visualization

    To use this function, additional dependencies are required, specifically, `dscribe`_
    and `sklearn`_ libraries used for the default dimensionality reduction. They can be
    installed with the following command:

    .. code:: bash

        pip install chemiscope[explore]

    Here is an example using this function with and without a featurizer function. The
    frames are obtained by reading the structures from a file that `ase <ase-io_>`_ can
    read, and performing Kernel PCA using `sklearn`_ on a descriptor computed with SOAP
    with `dscribe`_ library.

    .. code-block:: python

        import chemiscope
        import ase.io
        import dscribe.descriptors
        import sklearn.decomposition

        # Read the structures from the dataset
        frames = ase.io.read("trajectory.xyz", ":")

        # 1) Basic usage with the default featurizer (SOAP + PCA)
        chemiscope.explore(frames)


        # Define a function for dimensionality reduction
        def soap_kpca_featurize(frames, environments):
            if environments is not None:
                raise ValueError("'environments' are not supported by this featurizer")
            # Compute descriptors
            soap = dscribe.descriptors.SOAP(
                species=["C"],
                r_cut=4.5,
                n_max=8,
                l_max=6,
                periodic=True,
            )
            descriptors = soap.create(frames)

            # Apply KPCA
            kpca = sklearn.decomposition.KernelPCA(n_components=2, gamma=0.05)

            # Return a 2D array of reduced features
            return kpca.fit_transform(descriptors)


        # 2) Example with a custom featurizer function
        chemiscope.explore(frames, featurize=soap_kpca_featurize)

    For more examples, see the related `documentation <chemiscope-explore_>`_.

    .. _ase-io: https://wiki.fysik.dtu.dk/ase/ase/io/io.html
    .. _sklearn: https://scikit-learn.org/
    .. _dscribe: https://singroup.github.io/dscribe/latest/
    .. _chemiscope-explore: https://chemiscope.org/docs/examples/6-explore.html
    """

    # Validate inputs
    if featurize is not None and not callable(featurize):
        raise TypeError("'featurize' must be a callable (function)")
    if properties is None:
        properties = {}

    # Apply dimensionality reduction from the provided featurizer
    if featurize is not None:
        X_reduced = featurize(frames, environments)
    else:
        # Run default featurizer
        X_reduced = soap_pca_featurize(frames, environments)

    # Add dimensionality reduction results to properties
    properties["features"] = X_reduced

    # Return chemiscope widget
    return show(
        frames=frames, properties=properties, mode=mode, environments=environments
    )


def soap_pca_featurize(frames, environments=None):
    """
    Computes SOAP features for a given set of atomic structures and performs
    dimensionality reduction using PCA. Custom featurize functions should
    have the same signature.

    Note:
    - The SOAP descriptor parameters are pre-defined.
    - We use all available CPU cores for parallel computation of SOAP descriptors.
    """

    # Check if dependencies were installed
    try:
        from dscribe.descriptors import SOAP
        from sklearn.decomposition import PCA
    except ImportError as e:
        raise ImportError(
            f"Required package not found: {str(e)}. Please install dependency "
            + "using 'pip install chemiscope[explore]'."
        )
    centers = None

    # Get the atom indexes from the environments and pick related frames
    if environments is not None:
        centers = _extract_environment_indices(environments)

    # Pick frames and properties related to the environments if provided
    if environments is not None:
        # Sort environments by structure id and atom id
        environments = sorted(environments, key=lambda x: (x[0], x[1]))

        # Check structure indexes
        unique_structures = list({env[0] for env in environments})
        if any(index >= len(frames) for index in unique_structures):
            raise IndexError(
                "Some structure indices in 'environments' are larger than the number of"
                "frames"
            )

        if len(unique_structures) != len(frames):
            # only include frames that are present in the user-provided environments
            frames = [frames[index] for index in unique_structures]

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
    feats = soap.create(frames, centers=centers, n_jobs=n_jobs)

    # Compute pca
    pca = PCA(n_components=2)
    return pca.fit_transform(feats)


def _extract_environment_indices(envs):
    """
    Convert from chemiscope's environements to DScribe's centers selection

    :param: list envs: each element is a list of [env_index, atom_index, cutoff]
    """
    grouped_envs = {}
    for [env_index, atom_index, _cutoff] in envs:
        if env_index not in grouped_envs:
            grouped_envs[env_index] = []
        grouped_envs[env_index].append(atom_index)
    return list(grouped_envs.values())
