from ..jupyter import show
from ..structures import extract_properties
from ._metatomic import metatomic_featurizer


__all__ = ["explore", "metatomic_featurizer", "get_featurizer"]

FEATURIZER_VERSION_MAP = {"pet-mad-1.0": "1.0.0"}


def get_featurizer(name, device=None, batch_size=1):
    """
    Get a featurizer by name for feature extraction. Currently available version is:
    "pet-mad-1.0", which returns an instance of `PETMADFeaturizer
    <https://arxiv.org/abs/2506.19674>`_.

    :param str name: name of the featurizer. Must match one of the known versions.
        Currently available is "pet-mad-1.0"

    :param str device: optional. device to run the featurizer on (e.g., "cpu" or "cuda")

    :param int batch_size: optional. number of structures to process per batch during
        featurization, defaults to 1

    .. warning::

        Requires additional dependencies. Install it using:

        .. code:: bash

            pip install chemiscope[explore]
    """
    try:
        from pet_mad.explore import PETMADFeaturizer
    except ImportError as e:
        raise ImportError(
            f"Required package not found: {e}. Please install the "
            "dependencies with `pip install chemiscope[explore]`."
        )

    if not isinstance(name, str):
        raise TypeError(f"Featurizer name must be a string, got {type(name)}")

    versions = FEATURIZER_VERSION_MAP.keys()
    if name not in FEATURIZER_VERSION_MAP:
        raise ValueError(
            f"The featurizer version {name} does not exist. Available options are: "
            f"{', '.join(versions) if len(versions > 1) else versions[0]}"
        )

    version = FEATURIZER_VERSION_MAP.get(name)

    # TEMPORAIRE WORKAROUND TILL FIX IN PET_MAD REPO
    from urllib.request import urlretrieve  # noqa

    url = "https://huggingface.co/lab-cosmo/pet-mad/resolve/v1.1.0/models/pet-mad-v1.1.0.ckpt"
    checkpoint_path = "pet-mad-v1.1.0.ckpt"
    urlretrieve(url, checkpoint_path)

    featurizer = PETMADFeaturizer(
        version=version,
        pet_checkpoint_path=checkpoint_path,
        device=device,
        batch_size=batch_size,
    )
    # END OF WORKAROUND. NEXT LINE TO BE UNCOMMENT ONCE FIXED
    """
    featurizer = PETMADFeaturizer(
        version=version, device=device, batch_size=batch_size
    )
    """

    return featurizer


def explore(
    frames,
    featurizer=None,
    properties=None,
    environments=None,
    settings=None,
    mode="default",
    device=None,
    batch_size=1,
    write_input=None,
    **kwargs,
):
    """
    Automatically generate an interactive Chemiscope visualization of atomic structures.

    This function creates a low-dimensional representation of the input ``frames`` and
    displays them using a Chemiscope widget. It supports automatic featurization with
    `PETMADFeaturizer <https://arxiv.org/abs/2506.19674>`_ or a custom featurization
    function.

    If available, all properties are extracted automatically from the structures.

    If one does not specify a ``featurizer`` (or sets it as a ``None``), only properties
    will be displayed on the map visualizer panel, as long as there are at least two of
    them.

    Overall, the visualization can include: properties extracted from the frames,
    additional user-provided properties, features from either the built-in PET-MAD
    featurizer, with dimensionality reduction, or a custom user-provided featurization
    functions.

    :param list frames: list of frames

    :param featurizer: either string specifying a featurizer version (currently only
        'pet-mad-1.0'), a custom callable function, or None. Used to compute features
        and perform dimensionality reduction on the ``frames``. For automatic default
        option, use ``pet-mad-1.0``. The callable should take ``frames`` as the first
        argument and ``environments`` as the second argument. The return value must be a
        features array of shape ``(n_frames, n_features)`` if ``environments`` is
        ``None``, or ``(n_environments, n_features)`` otherwise. By providing a string
        version, the related ``PETMADFeaturizer`` is used.

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

    :param string write_input: optional. A path to save the chemiscope visualization.
        Afterwards, the file can be loaded using :py:func:`chemiscope.read_input`

    :param kwargs: additional keyword arguments passed to support backward
        compatibility. Currently, only the deprecated ``featurize`` argument is
        supported, which was renamed to ``featurizer``

    :return: a chemiscope widget for interactive visualization

    To use this function, additional dependencies are required, specifically, `pet_mad`_
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
        chemiscope.explore(frames, featurizer="pet-mad-1.0")

        # or
        featurizer = chemiscope.get_featurizer("pet-mad-1.0")
        chemiscope.explore(frames, featurizer=featurizer)


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
        chemiscope.explore(frames, featurizer=soap_kpca_featurize)

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

    if "featurize" in kwargs:
        if featurizer is not None:
            raise ValueError(
                "Both 'featurizer' and deprecated 'featurize' are provided"
            )
        featurizer = kwargs.pop("featurize")
        warnings.warn(
            "'featurize' was deprecated and renamed to 'featurizer'. Explicitly set "
            "this parameter to silence this warning",
            stacklevel=1,
        )

    if properties is None:
        properties = {}

    extracted_properties = extract_properties(frames, environments)
    merged_properties = {**extracted_properties, **properties}

    if isinstance(featurizer, str):
        featurizer_instance = get_featurizer(featurizer, device, batch_size)
        merged_properties["features"] = featurizer_instance(frames, environments)
    elif callable(featurizer):
        merged_properties["features"] = featurizer(frames, environments)
    elif featurizer is not None:
        versions = FEATURIZER_VERSION_MAP.keys()
        raise ValueError(
            "Featurizer must be a version (string) of the default featurizer, "
            "callable (function) or a None. Available default featurizer options are: "
            f"{', '.join(versions) if len(versions) > 1 else versions[0]}"
        )

    if mode != "structure" and len(merged_properties) < 2:
        raise ValueError(
            "Need at least two properties to visualize a map widget. "
            "Either provide additional properties or use a featurizer, e.g. "
            "'pet-mad-1.0'"
        )

    widget = show(
        frames=frames,
        properties=merged_properties,
        shapes=None,
        environments=environments,
        settings=settings,
        mode=mode,
    )

    if write_input is not None:
        widget.save(write_input)

    return widget
