"""
Simple chemiscope input
=======================

This example demonstrates the basic usage of the ``chemiscope`` package, reading
structures and properties from structures stored as
`ASE <https://wiki.fysik.dtu.dk/ase>`_  objects, and
preparing a chemiscope file to visualize them using :py:func:`chemiscope.show`. First,
import dependencies:
"""

# %%
#

import glob
import os

import ase.io

import chemiscope


# %%
#
# Load structures from an extended xyz file:

frames = ase.io.read("data/showcase.xyz", ":")

# %%
#
# A chemiscope widget can be used to visualize structures and properties.
# This generates a Chemiscope object that is rendered to an interactive
# widget when executed in a Jupyter notebook.

chemiscope.show(
    frames=frames,
    # quickly extract properties from the ASE frames. nb: if you're doing this for
    # sharing, don't forget to also include metadata such as units and description
    properties=chemiscope.extract_properties(frames, only=["dipole_ccsd", "ccsd_pol"]),
    # it's always good to set some metadata to explain what the dataset - title is bare
    # minimum
    meta=dict(name="Dipole and polarizability"),
    # it is possible to set all visualization parameters with a dictionary format.
    # this is a shortcut for the most basic ones
    settings=chemiscope.quick_settings(
        x="ccsd_pol[1]", y="ccsd_pol[2]", color="dipole_ccsd[1]"
    ),
)


# %%
#
# For sharing with collaborators, or when one does not want to use an interactive
# notebook, one can also write a JSON (or compressed JSON) file that contains all
# information about structures and properties, and can be viewed at
# `chemiscope.org <https://chemiscope.org>`_

# Save as a file that can be viewed at chemiscope.org

chemiscope.write_input(
    "showcase.json.gz",
    frames=frames,
    properties=chemiscope.extract_properties(frames, only=["dipole_ccsd", "ccsd_pol"]),
    meta=dict(name="Dipole and polarizability"),
    settings=chemiscope.quick_settings(
        x="ccsd_pol[1]", y="ccsd_pol[2]", color="dipole_ccsd[1]"
    ),
)

# %%
#
# In a notebook it is also possible to load a `.json` file and create an interactive
# widget from it. This is another way to share datasets with collaborators.

chemiscope.show_input("showcase.json.gz")

# %%
#
# When working with large datasets, it is also possible to save separately the structure
# data, and have chemiscope load them on demand. This can be done with a couple of
# utility functions. Note that you will need to share the structure files alongside the
# main dataset file, and that it will not be possible to use the standalone viewer
# at `chemiscope.org`, as it requires all data to be included in the JSON file.

# This will write the external structures as separate files `structure-*.json`
external_frames = chemiscope.write_external_structures(frames, prefix="structure")

# We also use this to demonstrate the 'structure' mode of chemiscope
chemiscope.show(
    frames=external_frames,
    mode="structure",
)

# The dataset file is smaller and will take up less browser memory when loaded
chemiscope.write_input(
    "showcase-nostructures.json.gz",
    frames=external_frames,
)

print("\nCompressed dataset files:")
for f in sorted(glob.glob("showcase*.json.gz")):
    size = os.path.getsize(f)
    print(f"  {f}  ({size} bytes)")

print("External structure files:")
for f in sorted(glob.glob("structure-*.json.gz")):
    print("  ", f)
