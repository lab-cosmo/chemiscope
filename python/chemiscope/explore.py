from .jupyter import show


def explore(frames, properties=None, featurize=None, mode="default"):
    """
    TODO: add a short description

    Use cases:

    1. Basic
        chemiscope.explore(frames, featurize=get_mace_pca)
        Result: calls 'get_mace_pca' for reduction and displays the result

    2. With custom properties
        chemiscope.explore(frames, get_mace_pca, properties)
        Result: uses reducer, displays reduced data + user properties
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
