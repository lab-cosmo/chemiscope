### Chemiscope Directive Documentation

**.. chemiscope**

The objective of this directive is to enable embedding the chemiscope widget into HTML pages, providing interactive visualization.

### Configuration Steps

1. **Add the `chemiscope.sphinx` extension to your `docs/conf.py` file:**
    ```python
    extensions = [
        ...
        "chemiscope.sphinx",
    ]
    ```

2. **Use the `chemiscope` directive in your `.rst` file:**
    This directive requires two arguments:
    - `filepath`: The absolute path to the JSON file containing the dataset.
    - `mode`: The type of widget to use. Options are `default`, `structure`, or `map`.

    Example usage:
    ```rst
    .. chemiscope::
        :filepath: /home/username/chemiscope-directive/docs/src/datasets/fig_base_001.json.gz
        :mode: default
    ```

Don't forget to change the filepath to your actual value!

3. **Build the HTML files from the `.rst` files:**
    In this repository, you can build the HTML files using `tox`:
    ```sh
    tox -e docs
    ```

4. **Run a local HTTP server to view the HTML files:**
    The `file://` protocol is not supported for this directive, so you need to serve the files over HTTP. Navigate to the `docs/build/html` directory and run:
    ```sh
    cd docs/build/html
    python3 -m http.server
    ```

Happy chemiscoping!

---
