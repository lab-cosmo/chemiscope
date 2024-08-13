import numpy as np


def metatensor_featurizer(
    model,
    extensions_directory=None,
    check_consistency=False,
    device=None,
):
    """
    Create a featurizer function using a `metatensor`_ model to obtain the features from
    structures. The model must be able to create a ``"feature"`` output.

    :param model: model to use for the calculation. It can be a file path, a Python
        instance of :py:class:`metatensor.torch.atomistic.MetatensorAtomisticModel`, or
        the output of :py:func:`torch.jit.script` on
        :py:class:`metatensor.torch.atomistic.MetatensorAtomisticModel`.
    :param extensions_directory: a directory where model extensions are located
    :param check_consistency: should we check the model for consistency when running,
        defaults to False.
    :param device: a torch device to use for the calculation. If ``None``, the function
        will use the options in model's ``supported_device`` attribute.

    :returns: a function that takes a list of frames and returns the features.

    To use this function, additional dependencies are required. They can be installed
    with the following command:

    .. code:: bash

        pip install chemiscope[metatensor]

    Here is an example using a pre-trained `metatensor`_ model, stored as a ``model.pt``
    file with the compiled extensions stored in the ``extensions/`` directory. To obtain
    the details on how to get it, see metatensor `tutorial
    <https://lab-cosmo.github.io/metatrain/latest/getting-started/usage.html>`_. The
    frames are obtained by reading structures from a file that `ase <ase-io_>`_ can
    read.

    .. code-block:: python

        import chemiscope
        import ase.io

        # Read the structures from the dataset frames =
        ase.io.read("data/explore_c-gap-20u.xyz", ":")

        # Provide model file ("model.pt") to `metatensor_featurizer`
        featurizer = chemiscope.metatensor_featurizer(
            "model.pt", extensions_directory="extensions"
        )

        chemiscope.explore(frames, featurize=featurizer)

    For more examples, see the related :ref:`documentation
    <chemiscope-explore-metatensor>`.

    .. _metatensor: https://docs.metatensor.org/latest/index.html
    .. _chemiscope-explore-metatensor:
        https://chemiscope.org/docs/examples/7-explore-advanced.html#example-with-metatensor-model
    """

    # Check if dependencies were installed
    try:
        from metatensor.torch.atomistic import ModelOutput
        from metatensor.torch.atomistic.ase_calculator import MetatensorCalculator
    except ImportError as e:
        raise ImportError(
            f"Required package not found: {e}. Please install the dependency using "
            "'pip install chemiscope[metatensor]'."
        )

    # Initialize metatensor calculator
    CALCULATOR = MetatensorCalculator(
        model=model,
        extensions_directory=extensions_directory,
        check_consistency=check_consistency,
        device=device,
    )

    def get_features(atoms, environments):
        """Run the model on a single atomic structure and extract the features"""
        outputs = {"features": ModelOutput(per_atom=environments is not None)}
        selected_atoms = _create_selected_atoms(environments)
        output = CALCULATOR.run_model(atoms, outputs, selected_atoms)

        return output["features"].block().values.detach().cpu().numpy()

    def featurize(frames, environments):
        if isinstance(frames, list):
            envs_per_frame = _get_environments_per_frame(environments, len(frames))

            outputs = [
                get_features(frame, envs) for frame, envs in zip(frames, envs_per_frame)
            ]
            return np.vstack(outputs)
        else:
            return get_features(frames, environments)

    return featurize


def _get_environments_per_frame(environments, num_frames):
    """
    Organize the environments for each frame

    :param list environments: a list of atomic environments
    :param int num_frames: total number of frames
    """
    envs_per_frame = [None] * num_frames

    if environments is None:
        return envs_per_frame

    frames_dict = {}

    # Group environments by structure_id
    for env in environments:
        structure_id = env[0]
        if structure_id not in frames_dict:
            frames_dict[structure_id] = []
        frames_dict[structure_id].append(env)

    # Insert environments to the frame indices
    for structure_id, envs in frames_dict.items():
        if structure_id < num_frames:
            envs_per_frame[structure_id] = envs

    return envs_per_frame


def _create_selected_atoms(environments):
    """
    Convert environments into ``Labels`` object, to be used as ``selected_atoms``

    :param environments: a list of atom-centered environments
    """
    import torch
    from metatensor.torch import Labels

    if environments is None:
        return None

    # Extract system and atom indices from environments, overriding the structure id to
    # be 0 (since we will only give a single frame to the calculator at the same time).
    values = torch.tensor([(0, atom_id) for _, atom_id, _ in environments])

    return Labels(names=["system", "atom"], values=values)
