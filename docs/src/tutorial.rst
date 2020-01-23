User tutorial
=============

The different panels and settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Input file format for chemiscope
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Creating an input file
----------------------

.. autofunction:: chemiscope_input::write_chemiscope_input


Input file structure
--------------------

The input file structure follows closely the `Dataset`_ interface from the code.
Here is another representation of what the JSON file should contain, in
pseudo-JSON format.

.. code-block:: javascript

    {
        // metadata of the dataset
        "meta": {
            // the name of the dataset
            "name": "this is my name"
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
