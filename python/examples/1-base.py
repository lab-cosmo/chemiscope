"""
Simple chemiscope input
=======================

This example demonstrates the basic usage of the ``chemiscope`` package, reading
structures and properties from an `ASE package <https://wiki.fysik.dtu.dk/ase>`_ and
preparing a chemiscope file to visualize them using :py:func:`chemiscope.show`. First,
import dependencies:
"""

# %%
#

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
# To share the visualization with collaborators or view it on `<chemiscope.org>`_, save
# the dataset as a JSON (or compressed JSON) file. This file contain all information
# about structures, properties, and visualization settings.

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
# In a Jupyter notebook, you can load a previously saved ``.json`` or ``.json.gz`` file
# to recreate the interactive widget.

chemiscope.show_input("showcase.json.gz")
