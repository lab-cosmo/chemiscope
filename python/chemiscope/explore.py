from .jupyter import show


def explore(frames, featurize, properties=None, mode="default"):
    """
    Visualize the dataset with dimensionality reduction.

    :param list frames: list of atomic structures

    :param dict properties: optional, dictionary of additional properties.

    :param str mode: optional, widget mode, either ``default``, ``structure`` or ``map``
    """
    # Validate inputs
    if not callable(featurize):
        raise TypeError("'featurize' must be a callable (function)")
    if properties is None:
        properties = {}

    # Apply dimensionality reduction from the provided featurizer
    X_reduced = featurize(frames)

    # Add dimensionality reduction results to properties
    settings = {"map": {}}
    axis = ["x", "y", "z"]
    for i in range(X_reduced.shape[1]):
        properties[f"Component {i + 1}"] = X_reduced[:, i].tolist()
        if i < len(axis):
            settings["map"][axis[i]] = {"property": f"Component {i + 1}"}

    # Return chemiscope widget
    return show(frames=frames, properties=properties, mode=mode, settings=settings)
