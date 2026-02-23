Chemiscope: interactive structure/property explorer for materials and molecules
===============================================================================

Welcome to the documentation of the `chemiscope`_ visualization tool, an
interactive structure/property explorer for materials and molecules. The goal of
chemiscope is to provide interactive exploration of large databases of materials
and molecules and help researchers find structure-properties correlations
inside such databases. The screenshot below shows an example of such database
being visualized with chemiscope. 
The :ref:`first part of this documentation <getting-started>`
describes the default interface of chemiscope and how to use it with your own database,
both using the stand-alone viewer or through the :ref:`python module <python-module>`
and/or the :ref:`jupyter widget <jupyter>`.

.. figure:: img/screenshot.png
    :align: center

    Screenshot of the `Qm7b`_ database visualized in the default chemiscope viewer


Chemiscope is built around reusable components that can be arranged in
different manners to create visualization adapted to different kinds of data. The
:ref:`second part of this documentation <dev-manual>` explains how to build the
code and use it in your own website to create new interfaces.


Features and capabilities
-------------------------

Chemiscope is designed to visualize chemical structures (using a viewer
based on `3dmol.js <https://3dmol.csb.pitt.edu/>`_) together with an interactive
scatter plot (using `plotly.js <https://plotly.com/javascript/>`_) that visualizes
the properties associated with each configuration. Properties can be associated
with structures (e.g. total energy, band gap, etc.) or with atoms (e.g. atomic
charges, coordination numbers, etc.). 

The scatter plot can be configured to
show different properties on each axis, in two or three dimensions, and further
properties can be used to color the points and size them. The user can
interactively select points in the scatter plot to visualize the corresponding
structures in the viewer, and vice versa. The structure viewer also supports 
coloring atoms according to atomic properties, and to display additional data 
in terms of "shapes" (e.g. vectors for forces, ellipsoid for tensorial quantities, 
etc.). 

The code is designed to handle large databases (up to hundreds of thousands
of structures) by using optimized data structures and rendering techniques.
Structures can be loaded dynamically from external files to reduce the
initial loading time, and the scatter plot rendering is optimized to handle
large numbers of points, dynamically resampling the data to provide fast interaction
at different levels of detail.

Chemiscope can be used as a web-based application, as a 
:ref:`jupyter widget <jupyter>`, as 
a visualization component embedded in 
:ref:`sphinx <sphinx>` or :ref:`sphinx-gallery <gallery>` documentation,
as well as a :ref:`streamlit component <streamlit>`. 


Getting and citing chemiscope
-----------------------------

Chemiscope is distributed under an open-source license, and you are welcome to
use it and incorporate it into your own research and software projects.
You can get the source from the
`GitHub repository <https://github.com/lab-cosmo/chemiscope>`_.
If you find it useful, we would appreciate a citation to the chemiscope
`paper`_:

G. Fraux, R. K. Cersonsky, M. Ceriotti, *Chemiscope: Interactive
Structure-Property Explorer for Materials and Molecules.* **Journal of Open
Source Software** 5 (51), 2117 (2020)

If you incorporate chemiscope components into a software project, a link
back to the `chemiscope`_ homepage is the preferred form of acknowledgment.


What's in this documentation?
-----------------------------

.. toctree::
    :maxdepth: 2

    Chemiscope <self>
    getting-started/index
    examples/index
    python/index
    embedding


.. _chemiscope: https://chemiscope.org
.. _paper: https://doi.org/10.21105/joss.02117
.. _Qm7b: https://doi.org/10.1088/1367-2630/15/9/095003
