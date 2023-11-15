import numpy as np


def center_shape(shape):
    """
    Takes a dictionary that describes a custom shape, and centers it, subtracting the
    mean of the vertex positions. Returns a shallow copy, with a new list for the
    vertices.

    :param shape: A dictionary describing a custom shape
    """

    if shape["kind"] != "custom":
        raise ValueError("Only `custom` shapes can be centered")

    new_shape = {k: shape[k] for k in shape}
    points = np.array(shape["vertices"])
    points -= points.mean(axis=0)
    new_shape["vertices"] = points.tolist()

    return new_shape


def _oriented_circle(radius, vec, n_points=20):
    """
    Generates a circle in 3D centered at the origin ands oriented according to the
    given vector
    """
    # makes sure the normal is unit
    nvec = vec / np.linalg.norm(vec)

    # Generate an arbitrary vector not collinear with n
    if nvec[0] or nvec[1]:
        vx = np.array([0, 0, 1])
    else:
        vx = np.array([0, 1, 0])

    # generate orthogonal vectors in the plane defined by nvec
    u = vx - np.dot(vx, nvec) * nvec
    u = u / np.linalg.norm(u)
    v = np.cross(nvec, u)

    # generate n_points in the plane defined by nvec, centered at vec
    angles = np.linspace(0, 2 * np.pi, n_points, endpoint=False)
    circle_points = radius * np.outer(np.cos(angles), u) + np.outer(np.sin(angles), v)

    return circle_points


def arrow_from_vector(
    vec, scale=1.0, radius=0.1, head_radius_scale=1.75, head_length_scale=2.0
):
    """
    Draws an arrow from the origin to the specified 3D position. Returns a custom shape
    in the form required by the chemiscope input. Use `None` for the arrow shape
    parameters to leave them undefined (so that they can be specified in the global
    parameters).

    :param scale: conversion from the units of the vector to the units of the atomic
        positions (usually Å)
    :param radius: radius of the stem of the arrow (same units as the atomic positions,
        typically Å)
    :param head_radius_scale: radius of the arrow tip, relative to the stem radius
    :param head_length_scale: length of the arrow tip, relative to the stem radius
    """

    data = {"vector": [v * scale for v in vec]}
    if radius is not None:
        data["baseRadius"] = radius
    if head_radius_scale is not None:
        data["headRadius"] = radius * head_radius_scale
    if head_length_scale is not None:
        data["headLength"] = radius * head_length_scale

    return data


def ellipsoid_from_tensor(tensor, scale=1.0, force_positive=False):
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
