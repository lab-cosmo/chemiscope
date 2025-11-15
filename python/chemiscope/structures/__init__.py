# -*- coding: utf-8 -*-
from ._ase import (  # noqa: F401
    _ase_all_atomic_environments,
    _ase_extract_properties,
    _ase_to_json,
    _ase_valid_structures,
    ase_merge_pi_frames,
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
    center_shape,
    ellipsoid_from_tensor,
)
from ._stk import (  # noqa: F401
    _stk_all_atomic_environments,
    _stk_to_json,
    _stk_valid_structures,
    convert_stk_bonds_as_shapes,
)


def _chemiscope_valid_structures(frames):
    """
    Check if the given frames are already in chemiscope format.

    :param frames: iterable over structures (typically a list of frames)
    :return: tuple (frames as list, boolean indicating if frames are valid)
    """

    if not hasattr(frames, "__iter__"):
        return frames, False

    first_frame = frames[0]
    if not isinstance(first_frame, dict):
        return frames, False

    required_keys = {"size", "names", "x", "y", "z"}
    if "size" not in first_frame:
        return frames, False

    for frame in frames:
        assert isinstance(frame, dict), (
            "inconsistent frame types: "
            f"expected dict but got {frame.__class__.__name__}"
        )

        # this is an external frame, so we skip the detailed checks
        if "data" not in frame:
            assert required_keys.issubset(set(frame.keys())), (
                "invalid chemiscope frame: "
                f"missing keys {required_keys - set(frame.keys())}"
            )

    return frames, True


def _guess_adapter(frames):
    """
    Guess which adapter to use for the given frames. This function return the
    frames as a list and a string describing which adapter should be used.
    """

    chemiscope_frames, use_chemiscope = _chemiscope_valid_structures(frames)
    if use_chemiscope:
        return chemiscope_frames, "chemiscope"

    ase_frames, use_ase = _ase_valid_structures(frames)
    if use_ase:
        return ase_frames, "ASE"

    stk_frames, use_stk = _stk_valid_structures(frames)
    if use_stk:
        return stk_frames, "stk"

    mda_frames, use_mda = _mda_valid_structures(frames)
    if use_mda:
        return mda_frames, "mda"

    raise Exception(f"unknown frame type: '{frames[0].__class__.__name__}'")


def frames_to_json(frames):
    """
    Convert the given ``frames`` to the JSON structure used by chemiscope.

    This function is a shim calling specialized implementations for all the
    supported frame types. Currently supported frames types include
    chemiscope-compatible dicts, `ase.Atoms`, `stk.BuildingBlocks`_, and
    `MDAnalysis.AtomGroup`_ objects.

    :param frames: iterable over structures (typically a list of frames)
    """
    frames, adapter = _guess_adapter(frames)

    json_frames = []
    if adapter == "chemiscope":
        json_frames = frames
    elif adapter == "ASE":
        json_frames = [_ase_to_json(frame) for frame in frames]
    elif adapter == "stk":
        json_frames = [_stk_to_json(frame) for frame in frames]
    elif adapter == "mda":
        # Be careful of the lazy loading of `frames.atoms`, which is updated during the
        # iteration of the trajectory
        json_frames = [_mda_to_json(frames) for _ in frames.universe.trajectory]
    else:
        raise Exception("reached unreachable code")

    return json_frames


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

    elif adapter == "mda":
        return _mda_extract_properties(frames, only, environments)

    elif adapter == "stk":
        raise RuntimeError(
            "stk structures do not contain properties, you must manually provide them"
        )

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
    elif adapter == "stk":
        return _stk_all_atomic_environments(frames, cutoff)
    elif adapter == "mda":
        return _mda_all_atomic_environments(frames, cutoff)
    else:
        raise Exception("reached unreachable code")
