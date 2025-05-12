"""
Simple chemiscope input
=======================

This example demonstrates the basic usage of the chemiscope package, reading structures
and properties from an ASE package and preparing a chemiscope file to visualize them.
Use `chemiscope.show` in a Jupyter notebook for interactive visualization

"""

# %%
#

import ase.io

import chemiscope


# %%
#
# Load structures from an extended xyz file

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
    # it is possible to set _all_ visualization parameters with a dictionary format.
    # this is a shortcut for the most basic ones
    settings=chemiscope.quick_settings(
        x="ccsd_pol[1]", y="ccsd_pol[2]", color="dipole_ccsd[1]"
    ),
)


# %%
#
# For sharing with collaborators, or when one does not want to use an interactive
# notebook, one can also write a JSON (or compressed JSON) file that contains all
# information about structures and properties, and can be viewed at chemiscope.org
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
# widget from it.

chemiscope.show_input("showcase.json.gz")
