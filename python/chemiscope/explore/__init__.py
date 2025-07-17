from ..jupyter import show
from ._metatomic import metatomic_featurizer


__all__ = [
    "explore",
    "metatomic_featurizer",
]


def explore(
    frames,
    featurize=None,
    properties=None,
    environments=None,
    settings=None,
    mode="default",
    device=None,
    batch_size=1,
):
    """
    Automatically generate an interactive Chemiscope visualization of atomic structures.

    This function creates a low-dimensional representation of the input ``frames`` and
    displays them using a Chemiscope widget. It supports automatic featurization with
    `PETMADFeaturizer <https://arxiv.org/abs/2506.19674>`_ or a custom featurization
    function.

    The default ``PETMADFeaturizer`` computes `PET-MAD
    <https://arxiv.org/abs/2503.14118>`_ features from the structures and projects them
    into the 3D MAD latent space.

    :param list frames: list of frames

    :param callable featurize: Optional callable to compute features and perform
        dimensionality reduction on the ``frames``. The callable should take ``frames``
        as the first argument and ``environments`` as the second argument. The return
        value must be a features array of shape ``(n_frames, n_features)`` if
        ``environments`` is ``None``, or ``(n_environments, n_features)`` otherwise. If
        ``None``, a default ``PETMADFeaturizer`` is used.

    :param dict properties: optional. Additional properties to be included in the
        visualization. This dictionary can contain any other relevant data associated
        with the atomic structures. Properties can be extracted from frames with
        :py:func:`extract_properties` or manually defined by the user.

    :param environments: optional. List of environments (described as ``(structure id,
        center id, cutoff)``) to include when extracting the atomic properties. Can be
        extracted from frames with :py:func:`all_atomic_environments` or manually
        defined.

    :param dict settings: optional dictionary of settings to use when displaying the
        data. Possible entries for the ``settings`` dictionary are documented in the
        chemiscope input file reference.

    :param str mode: optional. Visualization mode for the chemiscope widget. Can be one
        of "default", "structure", or "map". The default mode is "default".

    :param str device: torch device to use for the calculation with the default
            ``PETMADFeaturizer``. If `None`, we will try the options in the model's
            ``supported_device`` in order.

    :param int batch_size: optional. Number of structures processed in each batch with
        the default ``PETMADFeaturizer``.

    :return: a chemiscope widget for interactive visualization

    To use this function, additional dependencies are required, specifically, `pet-mad`_
    used for the default dimensionality reduction. They can be installed with the
    following command:

    .. code:: bash

        pip install chemiscope[explore]

    Here is an example using this function with and without a featurizer function. The
    frames are obtained by reading the structures from a file that `ase <ase-io_>`_ can
    read, and performing Kernel PCA using `sklearn`_ on a descriptor computed with SOAP
    using the `dscribe`_ library.

    .. code-block:: python

        import chemiscope
        import ase.io
        import dscribe.descriptors
        import sklearn.decomposition

        # Read the structures from the dataset
        frames = ase.io.read("trajectory.xyz", ":")

        # 1) Basic usage with the default featurizer (PET-MAD featurization + SMAP)
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
    .. _pet_mad: https://github.com/lab-cosmo/pet-mad
    .. _sklearn: https://scikit-learn.org/
    .. _dscribe: https://singroup.github.io/dscribe/latest/
    .. _chemiscope-explore: https://chemiscope.org/docs/examples/6-explore.html
    """
    # Check if dependencies were installed
    try:
        from pet_mad.explore import PETMADFeaturizer
    except ImportError as e:
        raise ImportError(
            f"Required package not found: {e}. Please install the "
            "dependencies with `pip install chemiscope[explore]`."
        )

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
        featurizer = PETMADFeaturizer(
            version="latest", device=device, batch_size=batch_size
        )
        X_reduced = featurizer(frames, environments)

    # Add dimensionality reduction results to properties
    properties["features"] = X_reduced

    return show(
        frames=frames,
        properties=properties,
        shapes=None,
        environments=environments,
        settings=settings,
        mode=mode,
    )
