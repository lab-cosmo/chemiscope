from ..jupyter import show
from ._soap_pca import soap_pca_featurize
from ._metatensor import metatensor_featurizer

__all__ = ["explore", "metatensor_featurizer"]


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

    return show(
        frames=frames,
        properties=properties,
        environments=environments,
        mode=mode,
    )
