import warnings
from pathlib import Path


mta = None
mts = None
torch = None
vesin_metatomic = None


class metatomic_featurizer:
    """
    Create a featurizer function using a `metatomic`_ model to obtain the features from
    structures. The model must be able to create a ``"features"`` output.

    :param model: model to use for the calculation. It can be a file path, a Python
        instance of :py:class:`metatomic.torch.AtomisticModel`, or the output of
        :py:func:`torch.jit.script` on :py:class:`metatomic.torch.AtomisticModel`.
    :param extensions_directory: a directory where model extensions are located
    :param check_consistency: should we check the model for consistency when running,
        defaults to False.
    :param device: a torch device to use for the calculation. If ``None``, the function
        will use the options in model's ``supported_device`` attribute.
    :param length_unit: Unit of length used in the structures.
    :param variant: selects which feature output variant to use. By default, the main
        ``"features"`` output is used. To choose another variant, provide its name
        (e.g., ``"cos_sin"``), which will select the corresponding
        ``"features/<variant>"`` output.

    :returns: a function that takes a list of structures and returns the features.

    To use this function, additional dependencies are required. They can be installed
    with the following command:

    .. code:: bash

        pip install chemiscope[explore]

    Here is an example using a pre-trained `metatomic`_ model, stored as a ``model.pt``
    file with the compiled extensions stored in the ``extensions/`` directory. The
    structures are obtained by reading structures from a file that `ase <ase-io_>`_ can
    read.

    .. code-block:: python

        import chemiscope
        import ase.io

        # Read the structures from the dataset
        structures = ase.io.read("data/explore_c-gap-20u.xyz", ":")

        # Provide model file ("model.pt") to `metatensor_featurizer`
        featurizer = chemiscope.metatensor_featurizer(
            "model.pt", extensions_directory="extensions"
        )

        chemiscope.explore(structures, featurizer=featurizer)

    For more examples, see the related :ref:`documentation
    <chemiscope-explore-metatomic>`.

    .. _metatomic: https://docs.metatensor.org/metatomic/
    """

    def __init__(
        self,
        model,
        *,
        extensions_directory=None,
        check_consistency=None,
        device=None,
        length_unit="Angstrom",
        variant=None,
    ):
        # Check if dependencies were installed
        global mta, mts, torch, vesin_metatomic
        if mta is None or mts is None or torch is None or vesin_metatomic is None:
            try:
                import metatensor.torch as mts
                import metatomic.torch as mta
                import torch
                import vesin.metatomic as vesin_metatomic
            except ImportError as e:
                raise ImportError(
                    f"Required package not found: {e}. Please install the "
                    "dependencies with `pip install chemiscope[explore]`."
                )

        if isinstance(model, (str, Path)):
            self.model = mta.load_atomistic_model(
                model, extensions_directory=extensions_directory
            )
        else:
            self.model = model

        self.length_unit = length_unit
        self.check_consistency = check_consistency

        capabilities = self.model.capabilities()

        self.feature_output_name = mta.pick_output(
            "features", capabilities.outputs, variant
        )

        if capabilities.dtype == "float32":
            self.dtype = torch.float32
        else:
            assert capabilities.dtype == "float64"
            self.dtype = torch.float64

        self.device = _find_best_device(capabilities.supported_devices, device)

    def __call__(self, structures, environments):
        systems = mta.systems_to_torch(structures)
        vesin_metatomic.compute_requested_neighbors(
            systems,
            self.length_unit,
            self.model,
        )

        systems = [s.to(self.dtype, self.device) for s in systems]

        if environments is not None:
            capabilities = self.model.capabilities()
            if not capabilities.outputs[self.feature_output_name].per_atom:
                raise ValueError(
                    (
                        "this model does not support per-atom features calculation for "
                        f"'{self.feature_output_name}' output"
                    )
                )

            selected_atoms = mts.Labels(
                names=["system", "atom"],
                values=torch.tensor(
                    [(system, atom) for system, atom, _ in environments]
                ),
            )
        else:
            selected_atoms = None

        options = mta.ModelEvaluationOptions(
            length_unit=self.length_unit,
            outputs={
                self.feature_output_name: mta.ModelOutput(
                    per_atom=environments is not None
                )
            },
            selected_atoms=selected_atoms,
        )

        outputs = self.model(
            systems,
            options,
            check_consistency=self.check_consistency,
        )

        # outputs should already be in the correct order because we passed
        # `selected_atoms` to the model
        return outputs[self.feature_output_name].block().values.detach().cpu().numpy()


def _find_best_device(supported_devices, requested_device):
    """
    Find the best device from the list of ``devices`` that is available to the current
    PyTorch installation.
    """

    available = []
    for device in supported_devices:
        if device == "cpu":
            available.append("cpu")
        elif device == "cuda":
            if torch.cuda.is_available():
                available.append("cuda")
        elif device == "mps":
            if (
                hasattr(torch.backends, "mps")
                and torch.backends.mps.is_built()
                and torch.backends.mps.is_available()
            ):
                available.append("mps")
        else:
            warnings.warn(
                f"unknown device in the model's `supported_devices`: '{device}'",
                stacklevel=3,
            )

    if requested_device is None:
        if len(available) == 0:
            warnings.warn(
                "could not find a valid device in the model's `supported_devices`, "
                "falling back to CPU",
                stacklevel=3,
            )
            return "cpu"

        return available[0]

    if requested_device in available:
        return requested_device

    warnings.warn(
        f"the requested device '{requested_device}' is not available on this machine "
        "or not supported by this model, falling back to CPU",
        stacklevel=3,
    )
    return "cpu"
