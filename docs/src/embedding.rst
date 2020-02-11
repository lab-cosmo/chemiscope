Chemiscope as a library
=======================

It is possible to use chemiscope as a software library when writing your own
web-based interface. This page document how to get the library and give a few
usage examples. You may also want to look at the `API documentation
<api/index.html>`_ for all classes, interfaces and functions in chemiscope.

Dependencies
^^^^^^^^^^^^

Chemiscope relies on different external dependencies that you should load in all
the HTML pages using it. You can serve these from your own webserver, or use a
CDN to deliver them.

- `JSmol <http://jmol.org/>`_ for structure visualization;
- `Bootstrap <https://getbootstrap.com/>`_ for HTML styling and basic UI;
- Both JSmol and Bootstrap rely on the ubiquitous `JQuery
  <https://jquery.com/>`_ and `JQueryUI <https://jqueryui.com/>`_

Getting and building the code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pre-built version
-----------------

The easiest way to do so is to download the latest release from `the release
page <https://github.com/cosmo-epfl/chemiscope/releases>`_ on GitHub. The main
file is ``chemiscope.min.js``, containing the needed JavaScript code to create
the default visualizer. This file exports a single global object ``Chemiscope``,
which contains references to

- `Chemiscope.DefaultVizualizer <DefaultVizualizer_>`_
- `Chemiscope.addWarningHandler <addWarningHandler_>`_
- `Chemiscope.StructureViewer <StructureViewer_>`_
- `Chemiscope.PropertiesMap <PropertiesMap_>`_
- `Chemiscope.EnvironmentInfo <EnvironmentInfo_>`_
- `Chemiscope.EnvironmentIndexer <EnvironmentIndexer_>`_

Partial builds are also available, in particular ``jsmol-widget.min.js`` which
only contains code related to `Chemiscope.JSmolWidget <JSmolWidget_>`_, making
the minified JavaScript file much smaller. Other partial builds containing only
part of chemiscope can be added upon request.

.. _DefaultVizualizer: api/classes/main.defaultvizualizer.html
.. _addWarningHandler: api/modules/utils.html#addwarninghandler
.. _StructureViewer: api/classes/structure.structureviewer.html
.. _PropertiesMap: api/classes/map.propertiesmap.html
.. _EnvironmentInfo: api/classes/info.environmentinfo.html
.. _EnvironmentIndexer: api/classes/utils.environmentindexer.html
.. _JSmolWidget: api/classes/structure.jsmolwidget.html

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

.. code-block:: html

    <!DOCTYPE html>
    <html>
        <head>
            <meta charset="utf-8">
            <title>Chemiscope basic example</title>

            <!-- Load all dependencies -->
            <!-- jquery -->
            <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.4.1/jquery.min.js" integrity="sha256-CSXorXvZcTkaix6Yvo6HppcZGetbYMGWSFlBw8HfCJo=" crossorigin="anonymous"></script>
            <script src="https://cdnjs.cloudflare.com/ajax/libs/jqueryui/1.12.1/jquery-ui.min.js" integrity="sha256-KM512VNnjElC30ehFwehXjx1YCHPiQkOPmqnrWtpccM=" crossorigin="anonymous"></script>

            <!-- bootstrap -->
            <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/4.4.1/css/bootstrap.min.css" integrity="sha256-L/W5Wfqfa0sdBNIKN9cG6QA5F2qx4qICmU2VgLruv9Y=" crossorigin="anonymous" />
            <script src="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/4.4.1/js/bootstrap.min.js" integrity="sha256-WqU1JavFxSAMcLP2WIOI+GB2zWmShMI82mTpLDcqFUg=" crossorigin="anonymous"></script>

            <!-- JSmol -->
            <script type="text/javascript" src="https://chemapps.stolaf.edu/jmol/jsmol-2019-10-30/JSmol.min.nojq.js"></script>

            <!-- Chemiscope code and default viewer code -->
            <script defer type="text/javascript" src="chemiscope.min.js"></script>
        </head>


        <body>
            <!-- Create basic DOM element for the different panels to live in -->
            <div id=map></div>
            <div id=structure></div>
            <div id=info></div>
        </body>


        <script type="text/javascript">
            // load data from anywhere
            const dataset = {
                meta: // to be loaded
                structures: // to be loaded
                properties: // to be loaded
                environments: // to be loaded (optional)
            }

            const config = {
                // id of the different elements
                map:       'map',
                info:      'info',
                structure: 'structure',
                // path to load J2S files for JSmol
                j2sPath:   'https://chemapps.stolaf.edu/jmol/jsmol-2019-10-30/j2s/',
            };

            Chemiscope.DefaultVizualizer.load(config, dataset);
        </script>
    </html>
