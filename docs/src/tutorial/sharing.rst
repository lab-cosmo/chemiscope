Sharing datasets with collaborators
===================================

Once you have converted your data in the :ref:`format used by chemiscope
<input-format>`, you might want to share it with collaborators. There are
multiple ways to do this, we'll go over them in this section.

Online visualizer at chemiscope.org
-----------------------------------

Uploading datasets
^^^^^^^^^^^^^^^^^^

The simplest way to share chemiscope dataset is to send the corresponding file
(``my-dataset.json`` or ``my-dataset.json.gz``) to your collaborators.

They can then load this file from the main chemiscope website:
https://chemiscope.org/. Here they can use the *Load/Save* menu to upload the
dataset to the website.

Loading from an URL
^^^^^^^^^^^^^^^^^^^

You can also host your dataset in any publicly accessible web server such as
your home page. Given a file hosted at
``https://university.edu/~myself/dataset.json``, you can then load this file
directly in chemiscope by going to
``https://chemiscope.org/?load=https://university.edu/~myself/dataset.json``.
This allow you to create links to directly open a given dataset in the main
chemiscope website, loading such dataset from your own webpage.

In general, you can set the ``load`` GET parameter on ``https://chemiscope.org``
to any url-encoded URL, and chemiscope will try to load the dataset from this
URL.

Saving visualization state
^^^^^^^^^^^^^^^^^^^^^^^^^^

You can save the current visualization settings from the website using the
*Load/Save* menu. This will allow you to download a JSON file, which you can use
later to reset visualization state.

If you are loading a file from an URL by creating a link that looks like
``https://chemiscope.org/?load=https://university.edu/dataset.json.gz``, you can
use the ``settings`` GET parameter to specify the URL of saved settings:
``https://chemiscope.org/?load=https://university.edu/dataset.json.gz&settings=https://university.edu/settings.json``

Standalone offline visualizer
-----------------------------

There are some cases where you do not want to depend on an online tool when
sharing your dataset, such as scientific article supplementation information.
For these use cases, a standalone, offline visualizer exists that uses the same
input file format as the default interface. You can download the latest version
of the standalone viewer at
:download:`https://chemiscope.org/chemiscope_standalone.html`.

This file contains all the required HTML and JavaScript code for chemiscope. You
can then add your own dataset by adding the corresponding JSON file at the end
of the ``chemiscope_standalone.html`` file.

**WARNING:** Only JSON, not compressed JSON (``.json.gz``) files are supported
with the standalone HTML visualizer.

.. code-block:: bash

    cat chemiscope_standalone.html my-dataset.json > my-dataset.html

You can then share the ``my-dataset.html`` file with others, who can open it in
any web browser.

To re-build the ``chemiscope_standalone.html`` file from sources, please follow
the steps below:

.. code-block:: bash

    git clone https://github.com/cosmo-epfl/chemiscope
    cd chemiscope
    npm install
    npm run build
    python3 ./utils/generate_standalone.py
