import warnings

SHAPE_PARAMS = {
    "ellipsoid": {"required": ["semiaxes"], "optional": [], "computable": []},
    "sphere": {"required": ["radius"], "optional": [], "computable": []},
    "custom": {
        "required": ["vertices"],
        "optional": [],
        "computable": ["simplices"],
    },
}


def _compute_simplices(dictionary, key):
    """
    Computes the indices of the facets for the
    point cloud stored in dictionary['vertices']
    using convex triangulation.
    """
    from scipy.spatial import ConvexHull

    return list(ConvexHull(dictionary[f"{key}_vertices"]).simplices)


SHAPE_FUNCS = {"simplices": _compute_simplices}

for k, v in SHAPE_PARAMS.items():
    for p in v["computable"]:
        assert p in SHAPE_FUNCS


def _get_shape_params(prefix, shape_key, dictionary):
    shape_dict = {"kind": shape_key}
    for p in SHAPE_PARAMS[shape_key]["required"]:
        if f"{prefix}_{p}" not in dictionary:
            raise KeyError("Missing required parameter {}.".format(f"{prefix}_{p}"))
        else:
            shape_dict[p] = dictionary[f"{prefix}_{p}"]

    for p in SHAPE_PARAMS[shape_key]["optional"]:
        if f"{prefix}_{p}" not in dictionary:
            warnings.warn("Missing optional parameter {}.".format(f"{prefix}_{p}"))
        else:
            shape_dict[p] = dictionary[f"{prefix}_{p}"]

    for p in SHAPE_PARAMS[shape_key]["computable"]:
        if f"{prefix}_{p}" not in dictionary:
            shape_dict[p] = SHAPE_FUNCS[p](dictionary, prefix)
        else:
            shape_dict[p] = dictionary[f"{prefix}_{p}"]

    return shape_dict


def _get_shape_params_i(prefix, shape_key, dictionary, i):
    shape_dict = {"kind": shape_key}
    for p in SHAPE_PARAMS[shape_key]["required"]:
        if f"{prefix}_{p}" not in dictionary:
            raise KeyError("Missing required parameter {}.".format(p))
        else:
            shape_dict[p] = dictionary[f"{prefix}_{p}"][i]

    for p in SHAPE_PARAMS[shape_key]["optional"]:
        if f"{prefix}_{p}" not in dictionary:
            warnings.warn("Missing optional parameter {}.".format(p))
        else:
            shape_dict[p] = dictionary[f"{prefix}_{p}"][i]

    for p in SHAPE_PARAMS[shape_key]["computable"]:
        if f"{prefix}_{p}" not in dictionary:
            shape_dict[p] = SHAPE_FUNCS[p](
                {k: v[i] for k, v in dictionary.items()}, prefix
            )
        else:
            shape_dict[p] = dictionary[f"{prefix}_{p}"]

    return shape_dict


def _extract_shapes(frame, key="shape"):
    if key in frame.info:
        if frame.info[key] not in SHAPE_PARAMS:
            raise KeyError(
                "The currently-supported shape types are {}, received {}.".format(
                    ", ".join(SHAPE_PARAMS.keys()), frame.info[key]
                )
            )
        shape = _get_shape_params(key, frame.info[key], frame.info)
        if "orientation" in frame.arrays:
            return {
                key: [
                    {**shape, "orientation": list(o)}
                    for o in frame.arrays["orientation"]
                ]
            }
        else:
            return {key: [shape for _ in frame]}

    elif key in frame.arrays:
        shapes = []
        for si, s in enumerate(frame.arrays[key]):
            if s not in SHAPE_PARAMS:
                raise KeyError(
                    "The currently-supported shape types are {}, received {}.".format(
                        ", ".join(SHAPE_PARAMS.keys()), s
                    )
                )
            shape = _get_shape_params_i(
                key,
                s,
                frame.arrays,
                si,
            )
            if "orientation" in frame.arrays:
                shape["orientation"] = list(frame.arrays["orientation"][si])
            shapes.append(shape)
        return {key: shapes}
