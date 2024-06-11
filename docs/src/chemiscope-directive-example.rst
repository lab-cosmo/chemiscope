.. _example-chemiscope:

==========================
Chemiscope Widget Example
==========================

This example demonstrates how to embed the chemiscope widget into an HTML page using the chemiscope directive.

Chemiscope Visualization
========================

The following directive embeds a chemiscope widget into this page. Make sure to specify the absolute path to your JSON dataset and choose the appropriate mode (`default`, `structure`, or `map`).

.. chemiscope::
    :filepath: datasets/fig_base_001.json.gz 
    :mode: structure 

Building the Documentation
==========================

To build the HTML documentation, follow these steps:

1. Add the `chemiscope.sphinx` extension to your `docs/conf.py` file:
   ```python
   extensions = ["chemiscope.sphinx"]
   ```

2. **Build the HTML files using `tox`:**
    ```sh
    tox -e docs
    ```

3. **Run a local HTTP server to view the HTML files:**
    ```sh
    cd docs/build/html
    python3 -m http.server
    ```
