from .jupyter import show


def explore(frames, properties=None, reducer=None, mode="default"):
    """
    Questions:
    - Which red. dimentionality algo to use par default?
    - What is user wants to display custom properties? => we take it from the input

    Use cases:

    1. Basic
    chemiscope.explore(frames)
    Result: default algo, display only reduced data

    2. With reducer
    chemiscope.explore(frames, reducer=get_mace_pca)
    Result: call 'get_mace_pca' and display the result

    3. With properties
    chemiscope.explore(frames, properties)
    or chemiscope.explore(frames, properties, reducer=get_mace_pca)
    Result: call default algo / reducer, display reduced data + user properties

    Flow general:
    0. Choose appropriate hyperparameters for the computation of descriptors (?)
    1. Compute descriptors
    2. Perform a dimensionality reduction on those features
    3. Create the widget
    """
    # Validate inputs
    if reducer is not None and not callable(reducer):
        raise TypeError("'reducer' must be a callable (function)")
    if properties is None:
        properties = {}

    # No reducer -> apply the default dimentionality reduction
    if reducer is None:
        X_reduced = apply_mace_tsne(frames)

    # Reducer provided -> call it
    else:
        X_reduced = reducer(frames)

    # Add properties with the dimensionality reduction
    properties["Component 1"] = X_reduced[:, 0].tolist()
    properties["Component 2"] = X_reduced[:, 1].tolist()

    # Return chemiscope widget
    return show(frames=frames, properties=properties, mode=mode)


def apply_mace_tsne(frames):
    pass
