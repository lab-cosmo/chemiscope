from typing import Optional

import numpy as np


def arrow_from_vector(
    position,
    *,
    scale: Optional[float] = 1.0,
    radius: Optional[float] = 0.1,
    head_radius_scale: Optional[float] = 1.75,
    head_length_scale: Optional[float] = 2.0,
):
    """
    Draws an arrow from the origin to the specified 3D ``position``. Returns a custom
    shape in the form required by the chemiscope input. Use ``None`` for the arrow shape
    parameters to leave them undefined (so that they can be specified in the global
    parameters).

    :param scale: conversion from the units of the vector to the units of the atomic
        positions (usually Å)
    :param radius: radius of the stem of the arrow (same units as the atomic positions,
        typically Å)
    :param head_radius_scale: radius of the arrow tip, relative to the stem radius
    :param head_length_scale: length of the arrow tip, relative to the stem radius
    """

    data = {"vector": [v * scale for v in position]}
    if radius is not None:
        data["baseRadius"] = radius
    if head_radius_scale is not None:
        data["headRadius"] = radius * head_radius_scale
    if head_length_scale is not None:
        data["headLength"] = radius * head_length_scale

    return data


def ellipsoid_from_tensor(
    tensor,
    *,
    scale=1.0,
    force_positive=False,
):
    """
    Returns an ellipsoid (semiaxes + quaternion) representation of a positive definite
    tensor (e.g. a polarizability), in the form required by the chemiscope input.

    :param tensor: a positive-definite tensor (3x3 or a 6-array [xx,yy,zz,xy,xz,yz])
    :param scale: conversion from the units of the tensor to the units of the atomic
        positions (usually Å)
    :param force_positive: takes the absolute value of eigenvalues, to handle
        non-positive tensors
    """

    tensor = np.asarray(tensor) * (scale**2)
    if len(tensor.flatten()) == 9:
        matrix = tensor.reshape((3, 3))
    elif len(tensor) == 6:
        matrix = np.array(
            [
                [tensor[0], tensor[3], tensor[4]],
                [tensor[3], tensor[1], tensor[5]],
                [tensor[4], tensor[5], tensor[2]],
            ]
        )
    # Compute the eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)

    if force_positive:
        eigenvalues = np.abs(eigenvalues)

    # The lengths of the principal axes are the square roots of the eigenvalues
    old_settings = np.seterr(invalid="raise")
    try:
        ax = np.sqrt(eigenvalues[0])
        ay = np.sqrt(eigenvalues[1])
        az = np.sqrt(eigenvalues[2])
    except FloatingPointError as e:
        raise ValueError(
            f"Non-positive definite tensor found with eigenvalues {eigenvalues}.\n"
            "If this is acceptable, set `force_positive=True` to take the "
            "absolute values.",
        ) from e
    np.seterr(**old_settings)

    # makes sure the rotation is proper
    if np.linalg.det(eigenvectors) < 0:
        eigenvectors[:, 0] *= -1

    # Form the rotation matrix from eigenvectors
    rotation = eigenvectors.T

    # converts to quaternion
    try:
        from scipy.spatial.transform import Rotation
    except ImportError as e:
        raise RuntimeError(
            "scipy is required to construct ellipsoids from tensors"
        ) from e
    quaternion = Rotation.from_matrix(rotation.T).as_quat()

    return dict(
        semiaxes=[ax, ay, az],
        orientation=list(quaternion),
    )
