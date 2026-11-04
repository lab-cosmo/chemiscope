from __future__ import annotations

from collections.abc import Sequence, Sized
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from ._chemfiles import ELEMENTS


if TYPE_CHECKING:
    from metatomic import System


Environment = tuple[int, int, float]


def _as_numpy(array: Any) -> NDArray[Any]:
    if hasattr(array, "detach"):
        return array.detach().cpu().numpy()
    if hasattr(array, "__array__"):
        return np.asarray(array)
    return np.from_dlpack(array)


def _metatomic_valid_structures(structures: object) -> tuple[object, bool]:
    try:
        from metatomic import System
    except ImportError:
        return structures, False

    if isinstance(structures, System):
        return [structures], True
    if isinstance(structures, Sequence) and structures:
        if isinstance(structures[0], System):
            if not all(isinstance(s, System) for s in structures):
                raise TypeError("all structures must be metatomic.System objects")
            return structures, True
    return structures, False


def _metatomic_to_json(system: System) -> dict[str, Any]:
    from metatomic.utils import unit_conversion_factor

    if not system.length_unit:
        raise ValueError("a metatomic System must have a length unit for visualization")
    factor = unit_conversion_factor(system.length_unit, "angstrom")
    positions = _as_numpy(system.positions) * factor
    types = _as_numpy(system.types)
    data = {
        "size": len(system),
        "names": [ELEMENTS[t] if 0 <= t < len(ELEMENTS) else str(t) for t in types],
        "elements": [ELEMENTS[t] if 0 <= t < len(ELEMENTS) else "X" for t in types],
        "x": positions[:, 0].tolist(),
        "y": positions[:, 1].tolist(),
        "z": positions[:, 2].tolist(),
    }
    cell = _as_numpy(system.cell) * factor
    if np.all(np.linalg.norm(cell, axis=1) != 0):
        data["cell"] = cell.reshape(-1).tolist()
    return data


def _metatomic_all_atomic_environments(
    structures: Sequence[Sized], cutoff: float
) -> list[Environment]:
    return [(i, a, cutoff) for i, s in enumerate(structures) for a in range(len(s))]
