User tutorial
=============

This tutorial will present how to use the `default chemiscope visualizer
<chemiscope_>`_ with your own database: the different panels and related
settings; as well as how to create an input file for it.

.. _chemiscope: https://chemiscope.org

Introduction to structural properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before we get started, we will introduce a few concepts used in chemiscope, and
in particular how structural properties can be used to create a *"map"* of
structures in 2D or 3D, in which similar structures are grouped together.
Everything starts with  **environments**. Chemiscope can work with either
full-structures environments or atom-centered environments. These environments
are fully defined by the positions of all atoms in the structure or around a
central atom. To be able to compare different environments, we use descriptors
based for example on `atom density representation <soap>`_ or `Behler-Parrinello
symmetry functions <Behler-Parrinello>`_. These descriptors are usually
high-dimensional vectors, hard to visualize and interpret. The next step is to
use a dimensionality reduction algorithm, such as `PCA`_, `sketch-map`_, *etc.*
The interpretation of the resulting map will differ depending on both the
descriptor used to represent the environments and the dimensionality reduction
algorithm.

.. figure:: img/mol-to-map.*
    :width: 65 %

    Illustration of the process used to create structural properties from a
    molecule.

Chemiscope is completly agnostic with respect to how structural properties are
generated, and do not provide any facilities to generate such structural
properties. In the rest of this document, we will refer to properties describing
the structure of an environment as *structural properties*; and other properties
associated with the environment (such as energy, density, ...) as *physical
properties*.

.. _soap: https://doi.org/10.1063/1.5090481
.. _Behler-Parrinello: https://doi.org/10.1103/physrevlett.98.146401
.. _PCA: https://en.wikipedia.org/wiki/Principal_component_analysis
.. _sketch-map: https://doi.org/10.1073/pnas.1108486108

Different panels and settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The default chemiscope visualizer is organized in three main panels: the map,
the structure viewer and the environment information display. Additionally,
clicking on the dataset title (on top of the map) will display some metadata
about the dataset (description, authors, references). This section will
present each one, as well as the main settings accessible to customize the
display.

The map is a 2D or 3D scatter plot showing properties for all the environments
in the dataset. You can set which properties (structural or physical) should be
used a the x, y, and potentially z axis; as well as for color and size of the
points. Additionally, properties which have string values (an not numeric
values) can be used as category data to set the symbols used for the points. To
open the settings modal window, click on the hamburger menu (the ☰ symbol) on
the left of the dataset title.

.. figure:: img/map.png
    :width: 80 %

    The map panel in 2D mode and the related settings

The structure panel is a 3D molecular viewer based on `Jmol`_. The settings are
accessible through the hamburger menu (☰) on the right of the viewer. The
settings are grouped into **representation** (how is the molecule rendered);
**supercell** (how many copies of the unit cell to display); **environments**
(how atom-centered environments are displayed); **camera** (reset the camera in
along one of the given axis); and **trajectory** (playback related settings).

.. figure:: img/structure.png
    :width: 80 %

    The structure panel and related settings

Finally, the environments information panel features sliders and text input to
allow for an easy selection of the environment of interest. The play button on
the left of the sliders activates the trajectory playback, looping over the
structures in the datasets or the atoms in a structure. By clicking on the
labels at the top (*structure XXX* and*atom XXX*), one can hide or show the
full property tables. These tables show all properties in the dataset for the
currently selected environment.

.. figure:: img/info.png
    :width: 40 %

    The environment information panel fully expanded

.. _Jmol: http://jmol.org

Input file format for chemiscope
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When using the default chemiscope interface, all the structures and properties
in a dataset are loaded from a single JSON file. This sections describe how to
generate such JSON file, either using a pre-existing python script that does
most of the work for you, or by writing the JSON file directly. Since the
resulting JSON file can be quite large and thus harder to share with
collaborators, the default chemiscope interface also allows to load JSON files
compressed with gzip.

Creating an input file
----------------------

The easiest way to create a JSON input file is to use the `chemiscope_input`_
Python 3 script that lives inside chemiscope's `github`_ repository. Download
the script and place it somewhere it can be imported by Python. Then, in your
own script, run the ``write_chemiscope_input`` function to generate the JSON
file. This script assumes you use the `ase`_ Python module to read the
structures.

.. autofunction:: chemiscope_input::write_chemiscope_input

.. _chemiscope_input: https://github.com/cosmo-epfl/chemiscope/blob/master/utils/chemiscope_input.py
.. _github: https://github.com/cosmo-epfl/chemiscope
.. _ase: https://wiki.fysik.dtu.dk/ase/index.html

Input file structure
--------------------

If you can not or do not want to use the script mentionned above, you can also
directly write the JSON file conforming to the schema described here. The input
file follows closely the `Dataset`_ typescript interface used in the library.
Using a pseudo-JSON format, the file should contains the following fields and
values:

.. code-block:: javascript

    {
        // metadata of the dataset. `description`, `authors` and `references`
        // will be rendered as markdown.
        "meta": {
            // the name of the dataset
            "name": "this is my name"
            // description of the dataset, OPTIONAL
            "description": "This contains data from ..."
            // authors of the dataset, OPTIONAL
            "authors": ["John Doe", "Mr Green, green@example.com"],
            // references for the dataset, OPTIONAL
            "references": [
                "'A new molecular construction', Journal of Random Words 19 (1923) pp 3333, DOI: 10.0000/0001100",
                "'nice website' http://example.com",
            ],
        },

        // list of properties in this dataset
        "properties": {
            // each property have a name, a target and some values
            <name>: {
                // the property target: is it defined per atom or for the full
                // structures
                "target": "atom" | "structure",
                // values of the properties can either be numbers or strings.
                // string properties are assumed to represent categories of
                // data.
                "values": [1, 2, 3, ...] | ["first", "second", "first", ...]
            }
        }

        // list of structures in this dataset
        "structures": [
            {
                // names of the atoms in the structure
                "names": ["H", "O", "C", "C", ...],
                // x cartesian coordinate of all the atoms, in Angstroms
                "x": [0, 1.5, 5.2, ...],
                // y cartesian coordinate of all the atoms, in Angstroms
                "y": [5.7, 7, -2.4, ...],
                // z cartesian coordinate of all the atoms, in Angstroms
                "z": [8.1, 2.9, -1.3, ...],
                // OPTIONAL: unit cell of the system, if any.
                //
                // This should be given as [ax ay az bx by bz cx cy cz], where
                // a, b, and c are the unit cell vectors. All values are
                // expressed in Angstroms.
                "cell": [10, 0, 0, 0, 10, 0, 0, 0, 10],
            },
            // other structures as needed
            ...
        ],

        // OPTIONAL: atom-centered environments descrptions
        //
        // If present, there should be one environment for each atom in each
        // structure.
        "environments": [
            {
                // index of the structure in the above structures list
                "structure": 0,
                // index of the central atom in structures
                "center": 8,
                // spherical cutoff radius, expressed in Angstroms
                "cutoff": 3.5,
            },
            // more environments
            ...
        ]
    }

.. _Dataset: api/interfaces/main.dataset.html


Using the standalone visualizer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The default chemiscope interface lives online, at https://chemiscope.org/. But
there are some cases where you do not want to use an online tool for your own
dataset, such as scientific article supplementation information. For these use
cases, a standalone, mostly offline visualizer exists that uses the same input
file format as the default interface.

To create a standalone visualizer with your own dataset, please follow the steps
below:

.. code-block:: bash

    git clone https://github.com/cosmo-epfl/chemiscope
    cd chemiscope
    npm install
    npm run build
    python3 ./utils/generate_standalone.py

This will create a ``standalone.html`` file containing all the required HTML and
javascript. You can then add your own dataset by adding the corresponding JSON
file at the end of the ``standalone.html`` file.

.. code-block:: bash

    cat standalone.html my-dataset.json > my-dataset.html
