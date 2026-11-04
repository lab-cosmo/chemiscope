Metatomic Core
==============

Core ``metatomic.System`` objects work with ``show``, ``create_input``,
``write_input``, ``write_external_structures`` and ``all_atomic_environments``.
ASE is not needed for this path. Install ``chemiscope[metatomic]``.

Structures
----------

.. code-block:: python

    import numpy as np
    import chemiscope
    from metatomic import System

    system = System(
        length_unit="angstrom",
        types=np.array([8, 1, 1], dtype=np.int32),
        positions=np.array([[0., 0., 0.], [0.76, 0.59, 0.], [-0.76, 0.59, 0.]]),
        cell=np.zeros((3, 3)),
        pbc=np.zeros(3, dtype=bool),
    )
    chemiscope.write_input("water.json", structures=[system])

Conversion conventions
----------------------

* Positions and cells are converted to Angstroms; an empty length unit is rejected.
* Types from 0 to 118 use the atomic-number convention (0 is ``X``). Other integer
  types keep their integer names and use element ``X``. For a model using a different
  type convention, supply explicit ``names`` and ``elements`` through Chemiscope
  structure dictionaries instead.
* Cells are included when all three cell vectors are nonzero, as for ASE inputs.
  Otherwise only atomic positions are displayed. Per-axis periodicity is not stored.
* NumPy and host-convertible Torch/JAX arrays are supported. Conversion does not
  change the input System's array backend. Raw DLPack arrays must be host-accessible.
