.. _input-reference:

Input file reference
====================

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
                // string properties are assumed to represent categories of
                // data.
                // the first dimension of the multidimensional property corresponds
                // to the number atoms/structures, the second dimension corresponds
                // to the size of the array of values per atom/structure. 
                "values": [1, 2, 3, ...] | ["first", "second", "first", ...] | [[1, 3, 5], [2, 4, 6], ...],

                // OPTIONAL: units of the property' value
                "units": "A/fs^2",
                // OPTIONAL: free-form description of the property as a string
                "description": "acceleration of the atoms in the structure ...",
                // OPTIONAL: an array containing a single string of the name of the parameter if the property is mutidimensional
                "parameter": ["name_parameter"]
            }
        }
        // OPTIONAL: list of parameters to be used with multidimensional properties
        "parameters": {
            // each Parameter must have a name, and an array of values that should match
            // the second dimension of the associated multidimensional properties
            <name>: {
                // an array of numbers containing the values of the parameter
                // the size should correspond to the second dimension of the corresponding
                // multidimensionl property
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
            },
            // other structures as needed
            ...
        ],

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
            },
            // Settings related to the structure viewers grid. This is an array
            // containing the settings for each separate viewer
            "structure": [
                {
                    // show bonds between atoms
                    "bonds": true,
                    //use space filling representation
                    "spaceFilling": false,
                    // show atoms labels
                    "atomLabels": false,
                    // show unit cell information and lines
                    "unitCell": false,
                    // displayed unit cell as a packed cell
                    "packedCell": false,
                    // number of repetitions in the `a/b/c` direction for the supercell
                    "supercell": [2, 2, 3],
                    // make the molecule spin
                    "rotation": false,
                    // which axis system to use
                    "axes": "none" | "xyz" | "abc",
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
                        "bgColor": "grey" | "CPK",
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
