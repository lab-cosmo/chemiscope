.. _input:

Chemiscope input files
======================

Chemiscope loads datasets from a single JSON file containing structures, properties, and
optional metadata. For large files, use gzip compression (e.g., ``output.json.gz``) to
ease sharing and loading.

This section covers creation methods, from provided Python tools or by manual JSON
writing.

For a simple way to generate an input file, we have a `Google Colab notebook
<https://colab.research.google.com/drive/1NU0gjtaHcB5Oc3NbFZiQYtctY2190hDu>`_.

The easiest programmatic way is the ``chemiscope`` Python package. Install it via:

.. code-block:: bash

   pip install chemiscope

You can then either use :py:func:`chemiscope.write_input` to save a JSON file and
display the widget using :py:func:`chemiscope.show_input` or see it through
chemiscope.org or use :py:func:`chemiscope.show` with the same arguments as
:py:func:`chemiscope.write_input` to just see it in your Jupyter notebook. This method
is best when you have structures (e.g., as ASE Atoms list) and properties in Python
scripts. See the :ref:`Python module documentation <python-module>` for full API
details.


Python quick example
~~~~~~~~~~~~~~~~~~~~
.. code-block:: python

   from chemiscope import write_input

   write_input(
       "output.json",
       structures=[
        ase.Atoms('H2O', positions=[(0, 0, 0), (0.76, 0.59, 0), (-0.76, 0.59, 0)])
       ],
       properties={"PCA[1]": [...], "PCA[2]": [...]},  # Your properties
   )


Input API
~~~~~~~~~

The chemiscope JSON file consists of these top-level entries: 

.. list-table::
   :header-rows: 1
   :widths: 20 50 15

   * - Key
     - Description
     - Required?
   * - ``meta``
     - Dataset metadata (name, authors, etc.)
     - No
   * - ``structures``
     - List of atomic structures (coordinates, cell, bonds)
     - Yes
   * - ``properties``
     - Properties mapped to atoms/structures (energy, forces, PCA, etc)
     - Yes
   * - ``shapes``
     - Custom shapes (spheres, arrows) to visualize in structures
     - No
   * - ``settings``
     - Default visualization settings (map axes, colors, etc.)
     - No
   * - ``environments``
     - Atom-centered neighborhoods
     - No
   * - ``parameters``
     - Parameters for multidimensional properties (e.g., time series)
     - No


Below is the detailed description of the values types and examples for each entry.

Metadata (``meta``)
~~~~~~~~~~~~~~~~~~~

Optional. Contains description of your dataset. The fields will be rendered as markdown.

.. list-table:: Metadata fields
   :header-rows: 1
   :widths: 15 15 30 15 15
   :class: tight-table

   * - Field
     - Type
     - Description
     - Required
     - Example
   * - ``name``
     - string
     - Short title of your dataset
     - Yes
     - ``"XXX dataset"``
   * - ``description``
     - string
     - Detailed explanation of content/origin
     - No
     - ``"A dataset from ..."``
   * - ``authors``
     - string[]
     - List of dataset authors
     - No
     - ``["Author 1", "Author 2"]``
   * - ``references``
     - string[]
     - Citations/links (one per array item)
     - No
     - ``["DOI:10.1234/abc", "'A new molecular construction', Journal of Random Words 19
       (1923) pp 3333, DOI: 10.0000/0001100", "'nice website' http://example.com"]``


Example:
++++++++

.. code-block:: json

   "meta": {
     "name": "MAD PCA",
     "description": "1000 validation structures from the *MAD dataset* under PCA",
     "authors": ["Author 1"],
     "references": ["https://arxiv.org/abs/2506.19674"]
   }


Properties (``properties``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Required. Defines all dataset properties to be available for a display inside the
widget. The field is a dictionary where each key is the property name (e.g., `"energy"`,
`"charge"`) and each value is the property definition, either in *full form* (explicit
dictionary with target and metadata) or *short form* (values only, target inferred from
array length).

Full definition
+++++++++++++++

Each property is a dictionary with the following fields:

.. list-table:: Property fields
   :header-rows: 1
   :widths: 15 15 30 15 15
   :class: tight-table

   * - Field
     - Type
     - Description
     - Required
     - Example
   * - ``target``
     - string
     - Scope of the property: ``"atom"`` or ``"structure"``
     - Yes
     - ``"atom"``
   * - ``values``
     - number[] | string[] | number[][]
     - Values (scalars or arrays; see shapes below)
     - Yes
     - ``[1.2, 3.4]``
   * - ``units``
     - string
     - Physical units of the values
     - No
     - ``"eV/Å"``
   * - ``description``
     - string
     - Property description
     - No
     - ``"DFT-calculated forces"``
   * - ``parameter``
     - string[]
     - For multidimensional: single parameter name (links to ``parameters`` top-level)
     - If multidimensional properties
     - ``["time"]``

Example:

.. code-block:: python

   properties = {
       "forces": {
           "target": "atom",
           "values": [[0.1, 0.0, -0.1], [0.2, -0.1, 0.0]],
           "units": "eV/Å",
           "description": "DFT-calculated forces"
       }
   }


Shortened definition
++++++++++++++++++++

You can omit the dictionary and provide only the array of values. The target (``"atom"``
or ``"structure"``) is automatically inferred by comparing the array length to
``n_atoms`` and ``n_structures``.

Example:

.. code-block:: python

   properties = {
       "energy": [-5.0, -4.8],  # inferred as "structure" if 2 structures
   }

If the number of atoms equals the number of structures, Chemiscope assumes
``"structure"`` by default and issues a warning.

Value types and shapes
++++++++++++++++++++++

.. list-table:: Supported value formats
   :header-rows: 1
   :widths: 15 25 40 20
   :class: tight-table

   * - Type
     - Format
     - Details
     - Example
   * - **Scalar**
     - 1D list or ndarray
     - Length must match ``n_atoms`` (for ``target="atom"``) or ``n_structures`` (for
       ``target="structure"``).
     - ``[0.1, 0.2, 0.3]``
   * - **Multi-dimensional**
     - 2D ndarray (``[n_samples, n_components]``)
     - Generates multiple properties (e.g., ``PCA[1]``, ``PCA[2]``)
     - ``[[0.1, 0.2], [0.3, 0.4]]``
   * - **Categorical**
     - List of strings
     - Used for symbols setting
     - ``["H", "O", "O"]``


Settings (``settings``)
~~~~~~~~~~~~~~~~~~~~~~~

Optional. Preconfigures the default visualization settings for the dataset, such as map
axes, colors, and structure viewer options. The ``settings`` field is a dictionary with
top-level keys for different panels.

Top-level settings
++++++++++++++++++

.. list-table:: Top-level settings
   :header-rows: 1
   :widths: 15 15 30 15 15
   :class: tight-table

   * - Field
     - Type
     - Description
     - Required
     - Example
   * - ``target``
     - string
     - Default view: ``"atom"`` or ``"structure"`` (requires environments for
       ``"atom"``)
     - No
     - ``"atom"``
   * - ``map``
     - dict
     - Map panel settings (axes, colors, etc.)
     - No
     - See below
   * - ``structure``
     - array[dict]
     - Structure viewer settings (one dict per viewer; up to 9)
     - No
     - See below
   * - ``pinned``
     - integer[]
     - Indices of environments/structures to pin in viewers (up to 9; defaults to [0])
     - No
     - ``[0, 5, 10]``

Map settings (``map``)
++++++++++++++++++++++

Configures the scatter plot. Sub-keys for axes (x/y/z), color, size, etc.

.. list-table:: Map settings
   :header-rows: 1
   :widths: 15 15 30 15 25
   :class: tight-table

   * - Field
     - Type
     - Description
     - Required
     - Example
   * - ``x`` / ``y`` / ``z``
     - dict
     - Axis config: ``{"property": "<name>", "scale": "linear"|"log", "min": number,
       "max": number}``
     - No
     - ``{"property": "PCA[1]", "scale": "linear"}``
   * - ``color``
     - dict
     - Color config: like axis, plus ``"palette"`` option for colormap selection.  
       Supported palettes:  
       ``"inferno"``, ``"magma"``, ``"plasma"``, ``"viridis"``, ``"cividis"``,  
       ``"hsv"``, ``"twilight"``, ``"twilight_shifted"``.
     - No
     - ``{"property": "energy", "palette": "viridis"}``
   * - ``size``
     - dict
     - Size config: ``{"property": "<name>", "mode":
       "constant"|"linear"|"log"|"sqrt"|"inverse", "factor": 1–100, "reverse": bool}``
     - No
     - ``{"property": "volume", "mode": "sqrt", "factor": 50}``
   * - ``symbol``
     - string
     - Property name for categorical symbols (string values)
     - No
     - ``"phase"``
   * - ``markerOutline``
     - bool
     - Thin black outline on markers
     - No
     - ``true``
   * - ``joinPoints``
     - bool
     - Connect points with a thin black line (e.g., for trajectories)
     - No
     - ``false``


Structure settings (``structure``)
++++++++++++++++++++++++++++++++++

Array of dicts (one per viewer).

.. list-table:: Structure settings
   :header-rows: 1
   :widths: 15 15 30 15 15
   :class: tight-table

   * - Field
     - Type
     - Description
     - Required
     - Example
   * - ``bonds``
     - bool
     - Show bonds
     - No
     - true
   * - ``atoms``
     - bool
     - Show atoms
     - No
     - true
   * - ``spaceFilling``
     - bool
     - Use space-filling representation
     - No
     - false
   * - ``atomLabels``
     - bool
     - Show atom labels
     - No
     - false
   * - ``unitCell``
     - bool
     - Show unit cell
     - No
     - true
   * - ``supercell``
     - integer[3]
     - Repetitions in a/b/c directions
     - No
     - ``[2, 2, 2]``
   * - ``rotation``
     - bool
     - Auto-rotate molecule
     - No
     - false
   * - ``axes``
     - string
     - Axis system: ``"none"|"xyz"|"abc"``
     - No
     - ``"xyz"``
   * - ``keepOrientation``
     - bool
     - Maintain orientation when changing structures
     - No
     - true
   * - ``environments``
     - dict
     - Environment options: ``{"activated": bool, "center": bool, "cutoff": number,
       "bgStyle": "licorice"|"ball-stick"|"hide", "bgColor": "grey"|"CPK"|"property"}``
     - No
     - ``{"activated": true, "cutoff": 3.5}``
   * - ``color``
     - dict
     - Atom coloring: ``{"property": "element"|"<name>", "transform":
       "linear"|"log10"|"sqrt"|"inverse", "min": number, "max": number, "palette":
       "bwr"}``
     - No
     - ``{"property": "charge", "palette": "bwr"}``


Examples
++++++++

1. Basic 2D map with color:

.. code-block:: json

   "settings": {
     "map": {
       "x": {"property": "PCA[1]"},
       "y": {"property": "PCA[2]"},
       "z": {"property": ""},
       "color": {"property": "energy"}
     }
   }

2. Trajectory mode with pinned structures:

.. code-block:: json

   "settings": {
     "target": "structure",
     "map": {"joinPoints": true},
     "structure": [{
       "keepOrientation": true,
       "playbackDelay": 50
     }],
     "pinned": [0, 10, 20]
   }

Environments (``environments``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Optional. Specifies atom-centered environments for datasets where properties are
associated with specific atoms rather than entire structures. Omitting this field
defaults to structure-level visualization only.

The ``environments`` field is an array of objects, each defining an environment. To
generate, use helper function :py:func:`chemiscope.all_atomic_environments`.

Environment definition
++++++++++++++++++++++

Each environment is a dictionary with:

.. list-table:: Environment fields
   :header-rows: 1
   :widths: 15 15 30 15 15
   :class: tight-table

   * - Field
     - Type
     - Description
     - Required
     - Example
   * - ``structure``
     - integer
     - 0-based index of the structure in ``structures``
     - Yes
     - 0
   * - ``center``
     - integer
     - 0-based index of the central atom in the structure
     - Yes
     - 8
   * - ``cutoff``
     - number
     - Spherical cutoff radius (in Å; must be > 0)
     - Yes
     - 3.5


Examples
++++++++

.. code-block:: json

   "environments": [
     {"structure": 0, "center": 0, "cutoff": 3.5},
     {"structure": 0, "center": 1, "cutoff": 3.5},
     {"structure": 0, "center": 2, "cutoff": 3.5}
   ]

Shapes (``shapes``)
~~~~~~~~~~~~~~~~~~~

Optional. Defines custom 3D shapes (e.g., spheres, arrows) to overlay in the structure
viewer. This is useful for visualizing additional features like atomic forces (as
arrows), uncertainty ellipsoids, or custom meshes.

The ``shapes`` field is a dictionary where each key is a shape group name (e.g.,
"forces"), and the value is the shape group definition. Multiple groups can be defined
for different visualizations.

Shape group definition
++++++++++++++++++++++

Each shape group is a dictionary with the following fields:

.. list-table:: Shape group fields
   :header-rows: 1
   :widths: 15 15 30 15 15
   :class: tight-table

   * - Field
     - Type
     - Description
     - Required
     - Example
   * - ``kind``
     - string
     - Shape type: ``"sphere"``, ``"ellipsoid"``, ``"cylinder"``, ``"arrow"``, or
       ``"custom"``
     - Yes
     - ``"arrow"``
   * - ``parameters``
     - dict
     - Parameters at different levels: ``"global"`` (all shapes), ``"structure"`` (per
       structure), ``"atom"`` (per atom/environment)
     - Yes
     - See below

Parameter levels
++++++++++++++++

Parameters are specified in a nested dictionary and merged from general to specific:

- ``global``: Dictionary of default parameters applied to all shapes in the group.
- ``structure``: List of dictionaries, one per structure. Overrides ``global`` for that structure.
- ``atom``: Flat list of dictionaries, one per atom. Overrides ``structure`` and ``global`` for that atom.

Common parameters
+++++++++++++++++

These apply to any shape kind:

.. list-table:: Common parameters
   :header-rows: 1
   :widths: 15 15 30 15

   * - Parameter
     - Type
     - Description
     - Required
   * - ``position``
     - number[3]
     - Center position (defaults to structure origin or atom position)
     - No
   * - ``scale``
     - number
     - Scaling factor for size
     - No
   * - ``orientation``
     - number[4]
     - Rotation quaternion (x, y, z, w); ignored for ``cylinder`` and ``arrow``
     - No
   * - ``color``
     - string | hex
     - Color (e.g., ``"red"`` or ``0xFF0000``)
     - No

Kind-specific parameters
++++++++++++++++++++++++

Each shape kind supports additional parameters (in addition to the common ones above).

**Sphere**

.. list-table::
   :header-rows: 1
   :widths: 20 20 45 15
   :class: tight-table

   * - Parameter
     - Type
     - Description
     - Required
   * - ``radius``
     - number
     - Sphere radius.
     - Yes

**Ellipsoid**

.. list-table::
   :header-rows: 1
   :widths: 20 20 45 15
   :class: tight-table

   * - Parameter
     - Type
     - Description
     - Required
   * - ``semiaxes``
     - number[3]
     - Semi-axes lengths along x, y, z.
     - Yes

**Cylinder**

.. list-table::
   :header-rows: 1
   :widths: 20 20 45 15
   :class: tight-table

   * - Parameter
     - Type
     - Description
     - Required
   * - ``vector``
     - number[3]
     - Direction and length; tip at vector end.
     - Yes
   * - ``radius``
     - number
     - Cylinder radius.
     - Yes

**Arrow**

.. list-table::
   :header-rows: 1
   :widths: 20 20 45 15
   :class: tight-table

   * - Parameter
     - Type
     - Description
     - Required
   * - ``vector``
     - number[3]
     - Direction and length; tip at vector end.
     - Yes
   * - ``baseRadius``
     - number
     - Shaft radius.
     - Yes
   * - ``headRadius``
     - number
     - Arrowhead radius.
     - Yes
   * - ``headLength``
     - number
     - Arrowhead length; may extend beyond base if the vector is short.
     - Yes

**Custom**

.. list-table::
   :header-rows: 1
   :widths: 20 20 45 15
   :class: tight-table

   * - Parameter
     - Type
     - Description
     - Required
   * - ``vertices``
     - number[N][3]
     - List of vertex positions.
     - Yes
   * - ``simplices``
     - integer[M][3]
     - Mesh triangulation indices.
     - No

Examples
++++++++

1. Global spheres (same for all structures/atoms):

.. code-block:: json

   "shapes": {
     "highlight_spheres": {
       "kind": "sphere",
       "parameters": {
         "global": {"radius": 0.5, "color": "red"}
       }
     }
   }

2. Per-atom arrows (e.g., for forces):

.. code-block:: json

   "shapes": {
     "forces": {
       "kind": "arrow",
       "parameters": {
         "global": {"baseRadius": 0.1, "headRadius": 0.2, "headLength": 0.3},
         "atom": [  // length = total atoms
           {"vector": [1.0, 0.0, 0.0]},
           {"vector": [0.0, 1.0, 0.0]},
           /* ... */
         ]
       }
     }
   }

3. Custom mesh with per-structure scaling:

.. code-block:: json

   "shapes": {
     "tetrahedron": {
       "kind": "custom",
       "parameters": {
         "global": {
           "vertices": [[0,0,0], [1,0,0], [0,1,0], [0,0,1]],
           "simplices": [[0,1,2], [0,1,3], [0,2,3], [1,2,3]]
         },
         "structure": [  // length = number of structures
           {"scale": 1.0},
           {"scale": 2.0},
           /* ... */
         ]
       }
     }
   }



Structures (``structures``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Required. Contains all atomic configurations in your dataset. In the most common use
case, it is automatically converted internally from the list of ``ASE.Atoms`` object. In
details, each structure is defined as an object with the following fields:

.. list-table:: Structure fields
   :header-rows: 1
   :widths: 15 15 30 15 15
   :class: tight-table

   * - Field
     - Type
     - Description
     - Required
     - Example
   * - ``size``
     - integer
     - Number of atoms in the structure
     - Yes
     - ``5``
   * - ``names``
     - string[]
     - Chemical symbols (length must match ``size``)
     - Yes
     - ``["H", "O", "H"]``
   * - ``x``
     - number[]
     - X coordinates (Å)
     - Yes
     - ``[0.0, 1.5]``
   * - ``y``
     - number[]
     - Y coordinates (Å)
     - Yes
     - ``[0.0, 0.0]``
   * - ``z``
     - number[]
     - Z coordinates (Å)
     - Yes
     - ``[0.0, -1.5]``
   * - ``cell``
     - number[9]
     - Unit cell vectors as ``[ax,ay,az,bx,by,bz,cx,cy,cz]`` (Å)
     - No
     - ``[10,0,0,0,10,0,0,0,10]``
   * - ``bonds``
     - integer[][3]
     - Bonds as ``[[i,j,order],...]`` (0-based indices)
     - No
     - ``[[0,1,1]]``


Example:
++++++++

.. code-block:: json

   {
     "size": 3,
     "names": ["O", "H", "H"],
     "x": [0.0, 0.76, -0.76],
     "y": [0.0, 0.59, 0.59],
     "z": [0.0, 0.0, 0.0]
   }


Parameters (``parameters``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Optional. Defines the variables for multidimensional properties. This is required when
properties use the ``parameters`` field in the ``properties``.

The ``parameters`` field is a dictionary where each key is a name (referenced in
property ``parameters``), and the value is the parameter definition.

Parameter definition
++++++++++++++++++++

Each parameter is a dictionary with the following fields:

.. list-table:: Parameter fields
   :header-rows: 1
   :widths: 15 15 30 15 15
   :class: tight-table

   * - Field
     - Type
     - Description
     - Required
     - Example
   * - ``values``
     - number[]
     - Array of values (length must match the second dimension of linked properties)
     - Yes
     - ``[0, 10, 20]``
   * - ``name``
     - string
     - Display name
     - No
     - ``"Time steps"``
   * - ``units``
     - string
     - Units
     - No
     - ``"fs"``

Example
+++++++

.. code-block:: json

    "parameters": {
      "time": {
        "values": [0, 10, 20, 30],  // matches length in linked property
        "name": "Simulation time",
        "units": "fs"
      }
    }

    // linked property (from properties section)
    "properties": {
      "energy": {
        "target": "structure",
        "values": [[-1.0, -1.1, -1.2, -1.3], /* ... */],
        "parameters": ["time"],
        "units": "eV"
      }
    }

.. _ase: https://wiki.fysik.dtu.dk/ase/index.html
