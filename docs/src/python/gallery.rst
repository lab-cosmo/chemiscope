.. _gallery:

Use with `sphinx-gallery`
=========================

The chemiscope widgets double up in the documentation of python packages
using :py:mod:`sphinx_gallery` - a sphinx extension that executes python 
scripts and formats their output as static pages inside the documentation. 
You can look at :ref:`the examples <examples>` to see how the widgets get 
displayed inside a documentation, and read the documentation of 
``sphinx_gallery`` if you never used it.

Here we just provide a short explanation of how to modify an existing 
``sphinx-gallery`` setup to include chemiscope widgets in the examples.


Using widgets in the examples
-----------------------------

Showing a widget in an example file is as simple as using 
:py:func:`chemiscope.show` in the example file. The widgets will
be shown in a similar style as that shown in the 
:ref:`jupyter interface <jupyter>`, including the possibility
of showing only the structure or map viewers. 
In practice, this relies on the same custom RST directive used
for the :ref:`sphinx extension <sphinx>`. 


Setting up ``sphinx-gallery``
-----------------------------

In order to enable the display of chemiscope widgets in the documentation,
assuming you already have a ``sphinx_gallery`` setup, you need to activate
the custom scraper from the chemiscope module. 

The ``conf.py`` file in the documentation should include the appropriate 
extension (usually you'll be adding this to a list of other extensions)

.. code-block:: python

    extensions = [ "chemiscope.sphinx" ]

Furthermore, you need to specify that ``sphinx_gallery`` should use the
custom chemiscope scraper (again, you'll typically add this to a tuple
of other scrapers)

.. code-block:: python

    sphinx_gallery_conf = {
        "image_scrapers": ( ChemiscopeScraper(), ),
    }



