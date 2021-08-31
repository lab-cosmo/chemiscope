Using the standalone visualizer
===============================

The default chemiscope interface lives online, at https://chemiscope.org/. But
there are some cases where you do not want to use an online tool for your own
dataset, such as scientific article supplementation information. For these use
cases, a standalone, mostly offline visualizer exists that uses the same input
file format as the default interface. You can download the latest version of the
standalone viewer at
:download:`https://chemiscope.org/chemiscope_standalone.html`.

This file contains all the required HTML and JavaScript code for chemiscope.
You can then add your own dataset by adding the corresponding JSON file at the
end of the ``chemiscope_standalone.html`` file.

.. code-block:: bash

    cat chemiscope_standalone.html my-dataset.json > my-dataset.html

To re-build the ``chemiscope_standalone.html`` file from sources, please follow
the steps below:

.. code-block:: bash

    git clone https://github.com/cosmo-epfl/chemiscope
    cd chemiscope
    npm install
    npm run build
    python3 ./utils/generate_standalone.py
