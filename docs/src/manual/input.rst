.. _input:

Chemiscope input files
===============================

When using the default chemiscope interface, all the structures and properties
in a dataset are loaded from a single JSON file. These sections describe how to
generate such JSON file, either using a pre-existing python script that does
most of the work for you, or by writing the JSON file directly. Since the
resulting JSON file can be quite large and thus harder to share with
collaborators, the default chemiscope interface also allows to load JSON files
compressed with gzip.

tl;dr if you would like to generate a simple chemiscope for your dataset, we
have a `Google Colab notebook <https://colab.research.google.com/drive/1NU0gjtaHcB5Oc3NbFZiQYtctY2190hDu>`_
that can help!

Tools that can create chemiscope inputs
---------------------------------------

Note that chemiscope does not compute structural representations or
dimensionality reduction, and the command line interface works iff
there are mappable quantities in the file. You can generate such representations
or reductions with packages such as `ASAP`_ or `scikit-matter`_.
The `ASAP`_ structural analysis package is another tool that can directly
generate an output in chemiscope format.

The easiest way to create a JSON input file is to use the ``chemiscope``
:ref:`Python module <python-module>`.
Install the package with ``pip install chemiscope``, and use
:py:func:`chemiscope.write_input` or :py:func:`chemiscope.create_input` in your
own script to generate the JSON file.

If all the properties you want to include into chemiscope are already stored in
a file `ase`_ can read, the ``chemiscope`` python package also installs a
:ref:`chemiscope-input <chemiscope-input-cli>` command line script.


Input file reference
--------------------

If you can not or do not want to use the ``chemiscope`` python package to create
your input files, you can also directly write the JSON file conforming to the
schema described here. The input file follows closely the `Dataset`_ typescript
interface used in the library. Using a pseudo-JSON format, the file should
contains the following fields and values:

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
            // Each property have at least a name, a target and some values.
            // Optional entries for the units and descriptions can also be added.
            <name>: {
                // the property target: is it defined per atom or for the full
                // structures
                "target": "atom" | "structure",
                // values of the properties can either be numbers, strings
                // or array of numbers.
                //
                // string properties are assumed to represent categories of
                // data.
                //
                // the first dimension of the multidimensional property corresponds
                // to the number atoms/structures, the second dimension corresponds
                // to the size of the array of values per atom/structure.
                "values": [1, 2, 3, ...] | ["first", "second", "first", ...] | [[1, 3, 5], [2, 4, 6], ...],

                // OPTIONAL: units of the property' value
                "units": "A/fs^2",
                // OPTIONAL: free-form description of the property as a string
                "description": "acceleration of the atoms in the structure ...",
                // OPTIONAL: an array containing a single string of the name of
                // the parameter (from the `parameters`` object below). This is
                // required multidimensional properties
                "parameter": ["parameter_name"]
            }
        }
        // OPTIONAL: list of parameters to be used with multidimensional properties
        "parameters": {
            // each Parameter must have a name, and an array of values that should match
            // the second dimension of the associated multidimensional properties
            <name>: {
                // an array of numbers containing the values of the parameter
                // the size should correspond to the second dimension of the
                // corresponding multidimensional property
                "values": [0, 0.1, 0.2]

                // OPTIONAL free-form description of the parameter as a string
                "name": "a short description of this parameter"
                // OPTIONAL units of the values in the values array
                "units": "eV"

            }
        }

        // list of structures in this dataset
        "structures": [
            {
                // number of atoms in the structure
                "size": 42,
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
                // OPTIONAL: bonds in the system, if any.
                //
                // This should be given as [i, j, order], where i, j, and
                // order are integers. i and j are the indices of
                // the atoms bonded together, and order is the bond order,
                // which can be 1 (for single bonds) or 2 (for double bonds)
                // or 3 (for triple bonds).
                "bonds": [[0, 1, 1], [1, 2, 2],...],
            },
            // other structures as needed
            ...
        ],

        // OPTIONAL: shapes to display in the structure viewer.
        // Multiple shapes groups with different names are supported.
        //
        // Each shape is defined by parameters that can be specified globally,
        // at the structure level, or for individual atoms
        "shapes": {
            <name>: {
                "kind" : <"sphere", "ellipsoid", "cylinder", "arrow", "custom">,
                "parameters" : {
                    "global" : { <global_parameters_dictionary> },
                    "structure" : [ <list_of_structure_parameter_dictionaries> ],
                    "atom" : [ <list_of_atom_parameter_dictionaries> ]
                }
            },
            // Sphere shapes, with the given `radius`
            <other_name>: {
                "kind" : "sphere"
                "parameters" : {
                    "global" : { "radius" : 0.2 }
                }
            },
            // Ellipsoid shapes, with the given `[ax, ay, az]` semi-axes
            <other_name>: {
                "kind" : "ellipsoid"
                "parameters" : {
                    "global" : {},
                    "structure" : [ {"semiaxes": [1, 1, 2]}, ... ]
                }
            },
            // Cylinder, with the given radiys, and `vector` direction
            <other_name>: {
                "kind" : "cylinder"
                "parameters" : {
                    "global" : { "radius" : 0.2 },
                    "atom" : [ {"vector" : [0,0,1]}, {"vector": [0,1,1]}, ... ]
                }
            },
            // Arrow, with the given shape parameters, and `vector` direction
            <other_name>: {
                "kind" : "arrow"
                "parameters" : {
                    "global" : { "baseRadius" : 0.2, 'headRadius': 0.3, 'headLength' : 0.4 },
                    "atom" : [ {"vector" : [0,0,1]}, {"vector": [0,1,1]}, ... ]
                }
            },
            // Custom shapes. Must provide list of vertices, and the vertex
            // indices associated with simplices (the latter are autocalculated)
            // if omitted
            <yet_another> : {
                "kind" : "custom",
                "parameters" : {
                    "global" : { "vertices" : [[x1,y1,z1], [x2,y2,z2], ...],
                                 "simplices" : [[0,1,2], [1,3,4], [0,5,5]] },
                    "atom" : [ {"scale" : 1}, {"scale" : 0.5}, ... ]
                }
            }

        }

        // OPTIONAL: atom-centered environments descriptions
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

        // OPTIONAL: setting for each panel
        //
        // Adding these values allow to setup how a given dataset should be
        // visualized in chemiscope.
        //
        // Each value inside the settings group is optional
        "settings": {
            // Visualization display target, either per atom-centered environments or per structure.
            // Supported in default and structure visualizers; the atom visualizer uses the "atom"
            // target by default. To use "atom" target, make sure to provide a list of environments.
            "target": "atom" | "structure",
            // settings related to the map
            "map": {
                // x axis settings
                "x": {
                    // name of the property to use for this axis, this must be
                    // one of the key from the root `properties` table.
                    "property": "<name>",
                    // should the axis use linear or logarithmic scaling
                    "scale": "linear" | "log",
                    // lower bound of the axis
                    "min": -0.23,
                    // upper bound of the axis
                    "max": 1.42,
                },
                // y axis setting, using the the same keys as x axis setting
                "y": {
                    // ...
                },
                // z axis setting, using the the same keys as x axis setting
                "z": {
                    // property can be set to an empty string to get a 2D map
                    "property": "",
                    // ...
                },
                // name of the property to use for markers symbols, this must be
                // one of the key from the root `properties` table. The
                // associated property should have string values
                "symbol": "<name>",
                // point color setting, using the the same keys as x axis setting
                "color": {
                    // property can be set to an empty string for uniform color
                    "property": "",
                    // ...
                },
                // Color palette to use, default to 'inferno'
                "palette": "cividis",
                // settings related to the markers sizes
                "size": {
                    // scaling factor for the axis, between 1 and 100
                    "factor": 55,
                    // mode to scale the markers with respect to the properties
                      // `constant`: all markers are same size, scaled by `factor`
                      // `linear`: markers are directly proportional to the property
                      // `log`: markers are proportional to the logarithm of the property
                      // `sqrt`: markers are proportional to the square root of the property
                      // `inverse`: markers are inversely proportional to the property
                    "mode": "constant" | "linear" | "log" | "sqrt | "inverse"",
                    // name of the property to use for the markers size, this
                    // must be one of the key from the root `properties` table.
                    "property": "<name>",
                    // if false, markers scale from smallest to largest property value
                    // if true, marker scale from largest to smallest property value
                    // in the case of `inverse` scaling, this is reversed.
                    "reverse": false | true,
                },
                // whether points in the map should have a thin black outline
                "markerOutline": true | false,
                // whether the points in the map should be linked by a thin black trace
                "joinPoints": false | true,
            },
            // Settings related to the structure viewers grid. This is an array
            // containing the settings for each separate viewer
            "structure": [
                {
                    // show bonds between atoms
                    "bonds": true,
                    // show the atoms
                    "atoms": true,
                    // use space filling representation
                    "spaceFilling": false,
                    // show atoms labels
                    "atomLabels": false,
                    // show unit cell information and lines
                    "unitCell": false,
                    // number of repetitions in the `a/b/c` direction for the supercell
                    "supercell": [2, 2, 3],
                    // make the molecule spin
                    "rotation": false,
                    // which axis system to use
                    "axes": "none" | "xyz" | "abc",
                    // keep the orientation constant when loading a new structure
                    "keepOrientation": false,
                    // options related to atom-centered environments
                    "environments": {
                        // should we display environments & environments options
                        "activated": true,
                        // automatically center the environment when loading it
                        "center": false,
                        // the cutoff value for spherical environments
                        "cutoff": 3.5
                        // which style for atoms not in the environment
                        "bgStyle": "licorice" | "ball-stick" | "hide",
                        // which colors for atoms not in the environment
                        // it is possible to color those atoms by the property
                        // currently selected
                        "bgColor": "grey" | "CPK" | "property",
                    };
                    // options related to the coloring of the atoms
                    "color": {
                        // name of the property to use for coloring, this must be
                        // one of the key from the root `properties` table.
                        // the default value is "element"
                        "property": "element" | "<name>",
                        // if the atoms should not be colored by element,
                        // this is the transformation to apply to the property
                        // the default value is "linear"
                        // if the value of the selected property value of an atom
                        // is missing, the atom will be colored in light grey
                        // if the value is not a real number or infinite,
                        // the atom will be colored in dark grey
                        "transform": "linear" | "log10" | "sqrt" | "inverse",
                        // minimum property value to use for the color scale
                        // the color corresponding to this value will be used
                        // for atoms with a smaller value
                        // the min value should not be bigger than the max value
                        "min": "<number>",
                        // maximum property value to use for the color scale
                        // the color corresponding to this value will be used
                        // for atoms with a bigger value
                        // the max value should not be bigger than the min value
                        "max": "<number>",
                        // color palette used to color the atoms, default to 'bwr'
                        // coloring atoms from blue to white to red according to
                        // the property value.
                        "palette": "bwr",
                    };
                },
                // ...
            ]
            // List of environments to display (up to 9). These environments
            // will be shown in the structure viewer grid and indicated on
            // the map.
            //
            // This list should containg 0-based indexes of the environment in
            // the root "environments" object; or of the structure in the root
            // "environments" if no environments are present.
            //
            // If both this list and the "structure" settings list above are
            // present, they should have the same size and will be used together
            // (first element of "structure" setting used for the first "pinned"
            // value; and so on).
            //
            // This defaults to [0], i.e. showing only the first
            // environment/structure.
            "pinned": [
                33, 67, 12, 0,
            ]
        }
    }

.. _Dataset: api/interfaces/main.dataset.html

.. _ase: https://wiki.fysik.dtu.dk/ase/index.html
.. _ASAP: https://github.com/BingqingCheng/ASAP
.. _scikit-matter: https://scikit-matter.readthedocs.io/en/latest/
