# -*- coding: utf-8 -*-
import warnings

from ._ase import (  # noqa: F401
    _ase_all_atomic_environments,
    _ase_extract_properties,
    _ase_to_json,
    _ase_valid_structures,
    ase_tensors_to_ellipsoids,
    ase_vectors_to_arrows,
)
from ._mda import (  # noqa: F401
    _mda_all_atomic_environments,
    _mda_extract_properties,
    _mda_to_json,
    _mda_valid_structures,
)
from ._shapes import (  # noqa: F401
    arrow_from_vector,
    ellipsoid_from_tensor,
)
from ._stk import (  # noqa: F401
    _stk_all_atomic_environments,
    _stk_to_json,
    _stk_valid_structures,
)


def _chemiscope_valid_structures(structures):
    """
    Check if the given structures are already in chemiscope format.

    :param structures: iterable over structures
    :return: tuple (structures as list, boolean indicating if structures are valid)
    """

    if not hasattr(structures, "__iter__"):
        return structures, False

    first_structure = structures[0]
    if not isinstance(first_structure, dict):
        return structures, False

    required_keys = {"size", "names", "x", "y", "z"}
    if "size" not in first_structure:
        return structures, False

    for s in structures:
        assert isinstance(s, dict), (
            "inconsistent structure types: "
            f"expected dict but got {s.__class__.__name__}"
        )

        # this is an external structure, so we skip the detailed checks
        if "data" not in s:
            assert required_keys.issubset(set(s.keys())), (
                "invalid chemiscope structure: "
                f"missing keys {required_keys - set(s.keys())}"
            )

    return structures, True


def _guess_adapter(structures):
    """
    Guess which adapter to use for the given structures. This function return the
    structures as a list and a string describing which adapter should be used.
    """

    chemiscope_structures, use_chemiscope = _chemiscope_valid_structures(structures)
    if use_chemiscope:
        return chemiscope_structures, "chemiscope"

    ase_structures, use_ase = _ase_valid_structures(structures)
    if use_ase:
        return ase_structures, "ASE"

    stk_structures, use_stk = _stk_valid_structures(structures)
    if use_stk:
        return stk_structures, "stk"

    mda_structures, use_mda = _mda_valid_structures(structures)
    if use_mda:
        return mda_structures, "mda"

    raise Exception(f"unknown structure type: '{structures[0].__class__.__name__}'")


def structures_to_json(structures):
    """
    Convert the given ``structures`` to the JSON structure used by chemiscope.

    This function is a shim calling specialized implementations for all the
    supported structure types. Currently supported structures types include
    chemiscope-compatible dicts, ``ase.Atoms``, ``stk.BuildingBlocks``, and
    ``MDAnalysis.AtomGroup`` objects.

    :param structures: iterable over structures
    """
    structures, adapter = _guess_adapter(structures)

    if adapter == "chemiscope":
        json_data = structures
    elif adapter == "ASE":
        json_data = [_ase_to_json(s) for s in structures]
    elif adapter == "stk":
        json_data = [_stk_to_json(s) for s in structures]
    elif adapter == "mda":
        # Be careful of the lazy loading of `structures.atoms`, which is updated during
        # the iteration of the trajectory
        json_data = [_mda_to_json(structures) for _ in structures.universe.trajectory]
    else:
        raise Exception("reached unreachable code")

    return json_data


def extract_properties(structures=None, only=None, *, environments=None, frames=None):
    """
    Extract properties defined in the ``structures`` in a chemiscope-compatible format.

    :param structures: iterable over structures
    :param only: optional, list of strings. If not ``None``, only properties with a name
                from this list are included in the output.
    :param environments: optional, list of environments (described as ``(structure id,
        center id, cutoff)``) to include when extracting the atomic properties.
    """
    if frames is not None:
        warnings.warn(
            "`frames` argument is deprecated, use `structures` instead",
            stacklevel=2,
        )
        if structures is not None:
            raise ValueError("cannot use both `structures` and `frames` arguments")

        structures = frames

    structures, adapter = _guess_adapter(structures)

    if adapter == "ASE":
        return _ase_extract_properties(structures, only, environments)

    elif adapter == "mda":
        return _mda_extract_properties(structures, only, environments)

    elif adapter == "stk":
        raise RuntimeError(
            "stk structures do not contain properties, you must manually provide them"
        )

    else:
        raise Exception("reached unreachable code")


def all_atomic_environments(structures=None, cutoff=3.5, *, frames=None):
    """
    Generate a list of environments containing all the atoms in the given
    ``structures``. The optional spherical ``cutoff`` radius is used to display the
    environments in chemiscope.

    :param structures: iterable over structures
    :param float cutoff: spherical cutoff radius used when displaying the
                         environments
    """
    if frames is not None:
        import warnings

        warnings.warn(
            "`frames` argument is deprecated, use `structures` instead",
            stacklevel=2,
        )
        if structures is not None:
            raise ValueError("cannot use both `structures` and `frames` arguments")

        structures = frames

    structures, adapter = _guess_adapter(structures)

    if adapter == "ASE":
        return _ase_all_atomic_environments(structures, cutoff)
    elif adapter == "stk":
        return _stk_all_atomic_environments(structures, cutoff)
    elif adapter == "mda":
        return _mda_all_atomic_environments(structures, cutoff)
    else:
        raise Exception("reached unreachable code")
