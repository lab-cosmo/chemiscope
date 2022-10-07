# -*- coding: utf-8 -*-
from ._ase import (
    _ase_all_atomic_environments,
    _ase_composition_properties,
    _ase_extract_properties,
    _ase_librascal_atomic_environments,
    _ase_list_atom_properties,
    _ase_list_structure_properties,
    _ase_to_json,
    _ase_valid_structures,
)


def _guess_adapter(frames):
    """
    Guess which adapter to use for the given frames. This function return the
    frames as a list and a string describing which adapter should be used.
    """

    ase_frames, use_ase = _ase_valid_structures(frames)
    if use_ase:
        return ase_frames, "ASE"

    raise Exception(f"unknown frame type: '{frames[0].__class__.__name__}'")


def frames_to_json(frames):
    """
    Convert the given ``frames`` to the JSON structure used by chemiscope.

    This function is a shim calling specialized implementations for all the
    supported frame types. Currently only `ase.Atoms` frames are supported.

    :param frames: iterable over structures (typically a list of frames)
    """
    frames, adapter = _guess_adapter(frames)

    if adapter == "ASE":
        return [_ase_to_json(frame) for frame in frames]
    else:
        raise Exception("reached unreachable code")


def _list_atom_properties(frames):
    """
    List existing "atom" properties from the given ``frames``. This is used
    to check if the user might be missing some properties because chemiscope is
    no longer automatically extracting properties
    """
    frames, adapter = _guess_adapter(frames)

    if adapter == "ASE":
        return _ase_list_atom_properties(frames)
    else:
        raise Exception("reached unreachable code")


def _list_structure_properties(frames):
    """
    List existing "structure" properties from the given ``frames``. This is used
    to check if the user might be missing some properties because chemiscope is
    no longer automatically extracting properties
    """
    frames, adapter = _guess_adapter(frames)

    if adapter == "ASE":
        return _ase_list_structure_properties(frames)
    else:
        raise Exception("reached unreachable code")


def extract_properties(frames, only=None, environments=None):
    """
    Extract properties defined in the ``frames`` in a chemiscope-compatible
    format.

    :param frames: iterable over structures (typically a list of frames)
    :param only: optional, list of strings. If not ``None``, only properties
                with a name from this list are included in the output.
    :param environments: optional, list of environments (described as
        ``(structure id, center id, cutoff)``) to include when extracting the
        atomic properties.
    """
    frames, adapter = _guess_adapter(frames)

    if adapter == "ASE":
        return _ase_extract_properties(frames, only, environments)
    else:
        raise Exception("reached unreachable code")


def composition_properties(frames, environments=None):
    """
    Generate properties containing the chemical composition of the given
    ``frames``.

    This create two atomic properties: ``symbol`` (string) and ``number`` (int);
    and multiple structure properties: ``composition`` and ``n_{element}`` for
    each elements in the dataset. The properties are then returned in chemiscope
    format.

    :param frames: iterable over structures (typically a list of frames)
    :param environments: optional, list of environments (described as
        ``(structure id, center id, cutoff)``) to include when generating the
        atomic properties.
    """
    frames, adapter = _guess_adapter(frames)

    if adapter == "ASE":
        return _ase_composition_properties(frames, environments)
    else:
        raise Exception("reached unreachable code")


def all_atomic_environments(frames, cutoff=3.5):
    """
    Generate a list of environments containing all the atoms in the given
    ``frames``. The optional spherical ``cutoff`` radius is used to display the
    environments in chemiscope.

    :param frames: iterable over structures (typically a list of frames)
    :param float cutoff: spherical cutoff radius used when displaying the
                         environments
    """
    frames, adapter = _guess_adapter(frames)

    if adapter == "ASE":
        return _ase_all_atomic_environments(frames, cutoff)
    else:
        raise Exception("reached unreachable code")


def librascal_atomic_environments(frames, cutoff=3.5):
    """
    Generate the list of environments for the given ``frames``, matching the
    behavior used by librascal when computing descriptors for only a subset of
    the atomic centers. The optional spherical ``cutoff`` radius is used to
    display the environments in chemiscope.

    Only ``ase.Atoms`` are supported for the ``frames`` since that's what
    librascal uses.

    :param frames: iterable over ``ase.Atoms``
    :param float cutoff: spherical cutoff radius used when displaying the
                         environments
    """
    frames, adapter = _guess_adapter(frames)

    if adapter != "ASE":
        raise Exception("librascal_atomic_environments only supports ASE frames")

    return _ase_librascal_atomic_environments(frames, cutoff)
