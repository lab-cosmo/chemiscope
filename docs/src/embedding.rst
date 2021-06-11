Chemiscope as a library
=======================

It is possible to use chemiscope as a software library when writing your own
web-based interface. This page document how to get the library and give a few
usage examples. You may also want to look at the `API documentation
<api/index.html>`_ for all classes, interfaces and functions in chemiscope.

Dependencies
^^^^^^^^^^^^

Chemiscope relies on different external dependencies that you should load in all
the HTML pages using it. You can serve these from your own web server, or use a
CDN to deliver them.

- `Bootstrap <https://getbootstrap.com/>`_ for HTML styling and basic UI;
- Bootstrap relies on the ubiquitous `JQuery <https://jquery.com/>`_ and
  `JQueryUI <https://jqueryui.com/>`_

Getting and building the code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pre-built version
-----------------

The easiest way to do so is to download the latest release from `the release
page <https://github.com/cosmo-epfl/chemiscope/releases>`_ on GitHub. The main
file is ``chemiscope.min.js``, containing the code required to create the
default visualizer. This file exports a single global object ``Chemiscope``,
which contains references to

- `Chemiscope.DefaultVisualizer <DefaultVisualizer_>`_
- `Chemiscope.addWarningHandler <addWarningHandler_>`_
- `Chemiscope.MetadataPanel <MetadataPanel_>`_
- `Chemiscope.ViewersGrid <ViewersGrid_>`_
- `Chemiscope.PropertiesMap <PropertiesMap_>`_
- `Chemiscope.EnvironmentInfo <EnvironmentInfo_>`_
- `Chemiscope.EnvironmentIndexer <EnvironmentIndexer_>`_

Partial builds are also available, in particular ``molecule-viewer.min.js``
which only contains code related to `Chemiscope.MoleculeViewer
<MoleculeViewer_>`_, making the minified JavaScript file much smaller. Other
partial builds containing only part of chemiscope can be added upon request.

.. _DefaultVisualizer: api/classes/main.defaultvisualizer.html
.. _addWarningHandler: api/modules/utils.html#addwarninghandler-1
.. _ViewersGrid: api/classes/structure.viewersgrid.html
.. _PropertiesMap: api/classes/map.propertiesmap-1.html
.. _EnvironmentInfo: api/classes/info.environmentinfo-1.html
.. _MetadataPanel: api/classes/main.metadatapanel.html
.. _EnvironmentIndexer: api/classes/utils.environmentindexer.html
.. _MoleculeViewer: api/classes/structure.moleculeviewer.html

Build from sources
------------------

Chemiscope is written in TypeScript, a statically typed language which compiles
to JavaScript. It uses the standard JavaScript ecosystem tools for dependency
management, ``npm``. To build chemiscope from sources, you will first need to
get the sources, either as an archive from `the release page
<https://github.com/cosmo-epfl/chemiscope/releases>`_, or using git

.. code-block:: bash

    git clone https://github.com/cosmo-epfl/chemiscope

You will also need `node.js <https://nodejs.org/en/>`_ and `npm
<https://docs.npmjs.com/cli/npm>`_, which you can install with your favorite
package manager.

.. code-block:: bash

    cd chemiscope
    npm install
    npm run build

This should create a ``dist`` directory, containing all the minified JavaScript
libraries.

Usage inside a JavaScript project
---------------------------------

If you already have a JavaScript project using ``npm`` or ``yarn``, you can use
the version of Chemiscope available on npm. Add it to your project with

.. code-block:: bash

    npm install chemiscope

And import the library inside your own code using

.. code-block:: js

    const Chemiscope = require("chemiscope");

If your are using TypeScript, definition files are also provided with the npm
package, and should give you auto-completion, inline documentation and interface
checking.

Usage example
^^^^^^^^^^^^^

Below is the minimal HTML and JavaScript code needed to load and use the default
chemiscope interface. Additional code would be required to load the dataset,
using JSON files or directly.

.. literalinclude:: embedded-example.html
    :language: html
