# -*- coding: utf-8 -*-
from ._ase import (
    _ase_all_atomic_environments,
    _ase_atom_properties,
    _ase_librascal_atomic_environments,
    _ase_structure_properties,
    _ase_structures,
    _ase_to_json,
)


def _guess_adapter(frames):
    """
    Guess which adapter to use for the given frames. This function return the
    frames as a list and a string describing which adapter should be used.
    """

    ase_frames, use_ase = _ase_structures(frames)
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


def atom_properties(frames, composition, atoms_mask=None):
    """
    Extract "atom" properties from the given ``frames``, and give them as a
    dictionary compatible with :py:func:`create_input`.

    This function is a shim calling specialized implementations for all the
    supported frame types. Currently only `ase.Atoms` frames are supported.

    :param frames: iterable over structures (typically a list of frames)
    :param composition: whether to also add properties containing information
                        about the chemical composition of the system
    :param atoms_mask: optional list of booleans containing which atoms should
                       be include in the output
    """
    frames, adapter = _guess_adapter(frames)

    if adapter == "ASE":
        return _ase_atom_properties(frames, composition, atoms_mask)
    else:
        raise Exception("reached unreachable code")


def structure_properties(frames, composition):
    """
    Extract "structure" properties from the given ``frames``, and give them as a
    dictionary compatible with :py:func:`create_input`.

    This function is a shim calling specialized implementations for all the
    supported frame types. Currently only `ase.Atoms` frames are supported.

    :param frames: iterable over structures (typically a list of frames)
    :param composition: whether to also add properties containing information
                        about the chemical composition of the system
    """
    frames, adapter = _guess_adapter(frames)

    if adapter == "ASE":
        return _ase_structure_properties(frames, composition)
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
