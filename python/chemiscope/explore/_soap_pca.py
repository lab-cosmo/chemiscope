import os


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


def _extract_environment_indices(environments):
    """
    Convert from chemiscope's environments to DScribe's centers selection

    :param: list environments: each element is a list of [env_index, atom_index, cutoff]
    """
    grouped_envs = {}
    for [env_index, atom_index, _cutoff] in environments:
        if env_index not in grouped_envs:
            grouped_envs[env_index] = []
        grouped_envs[env_index].append(atom_index)
    return list(grouped_envs.values())
