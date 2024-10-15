
Introduction to structural properties
=====================================

Before we get started, we will introduce a few concepts that underlie the
concept and the usage of chemiscope. Chemiscope is designed to help navigating
*structure-property maps*, i.e. 2D or 3D representations of a set of
atomic scale entities that reflect how structure influences materials
properties.

Chemiscope can work with two kinds of entities: full structures, or
atom-centred environments. A structure consists in a set of atoms, possibly
representing the periodic repeat unit of an infinite structure. An
environment consists in a set of neighbors that surround a central atom,
within a set cutoff radius.
In both cases, these entities are fully defined by the position and nature
of the atoms present in the structure, or in the neighborhood of the
environment center.

For each structure or environment, one may have computed *properties*,
e.g. the cohesive energy of a molecule, or the NMR chemical shielding of
a nucleus, or *structural representations*, i.e. functions of the
spatial arrangement of the atoms that achieve a description of the structure 
that is as complete as possible, yet concise. 
Examples of such representations are for instance
`atom density representationis <soap>`_ or `Behler-Parrinello
symmetry functions <Behler-Parrinello>`_. These representations are usually
high-dimensional vectors, hard to visualize and interpret. For this reason, one
usually applies a dimensionality reduction algorithm, such as `PCA`_, `sketch-map`_,
`PCovR`_, *etc.*  The interpretation of the resulting descriptor will differ depending on
both the descriptor used to represent the structures or environments and the
dimensionality reduction algorithm applied.

Chemiscope simplifies visualizing the correlations between structural
representations and properties associated with structures and environments,
by representing in an interactive fashion these atomic-scale entities as points
on a map, and by associating these points with an explicit, 3D visualization
of the structure of the material or molecule.

.. figure:: ../img/mol-to-map.*
    :width: 65 %

    Illustration of the process used to create structural properties from a
    molecule.

Chemiscope is completly agnostic with respect to how properties and structural
representations are generated, and do not provide any facilities to generate them.
In the rest of this document, we will refer to properties describing
the structure of an environment or structure as *structural properties*
and other associated properties associated (such as energy, density, ...) as
*physical properties*.

.. _soap: https://doi.org/10.1063/1.5090481
.. _Behler-Parrinello: https://doi.org/10.1103/physrevlett.98.146401
.. _PCA: https://en.wikipedia.org/wiki/Principal_component_analysis
.. _sketch-map: https://doi.org/10.1073/pnas.1108486108
.. _PCovR: https://doi.org/10.1088/2632-2153/aba9ef
