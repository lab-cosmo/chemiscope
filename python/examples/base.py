"""
Simple chemiscope input
=======================

This example demonstrates the basic usage of the chemiscope package, reading structures
and properties from an ASE package and preparing a chemiscope file to visualize them.
Use `chemiscope.show` in a Jupyter notebook for interactive visualization

"""

import ase.io

import chemiscope

# load structures
frames = ase.io.read("data/showcase.xyz", ":")

# write the chemiscope input file. use `chemiscope.show` to view directly in a
# Jupyter environment
chemiscope.write_input(
    path="showcase.json.gz",
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

chemiscope.show(frames, properties=chemiscope.extract_properties(frames), mode="structure")
# cs = chemiscope.show(frames, properties=chemiscope.extract_properties(frames), mode="structure") 
# print(cs._repr_html_())
