
Introduction to structural properties
=====================================

This section provides an overview of the concepts underlying the usage of chemiscope.
The tool is designed to help navigating *structure-property maps*, which are 2d- or
3d-dimensional embeddings that reflect how variations in atomic structure influence
materials properties.

Chemiscope works with two kinds of entities: full structures and atom-centered
environments. A structure contains all atoms in a configuration, e.g. the unit
cell of a crystal or the geometry of a molecule. An environment consists of a set of
neighbors that surround a central atom within a chosen cutoff radius. In both cases,
these entities are fully defined by the positions of the atoms and their chemical
elements.

Each structure or environment can be associated with *properties*. Some of these, such
as energy, forces, density, or NMR shielding, represent *physical properties* that
describe the system. Others are *structural representations* that encode the geometry
itself. These include a wide family of descriptors, ranging from traditional
hand-crafted functions, such as `Behler-Parrinello symmetry functions`_, to more recent
approaches such as `SOAP`_, Atomic Cluster Expansion (`ACE`_), or latent embeddings
extracted from graph neural network potentials (e.g. `PET-MAD`_). These
representations are usually high-dimensional vectors, hard to visualize and interpret.
For this reason, it is common to apply dimensionality reduction algorithm, such as
`PCA`_, `sketch-map`_, `PCovR`_, *etc.* These methods compress the descriptor into a
small number of components that can be visualized on a map. The interpretability of the
low-dimensional embedding depends both on the representation that was originally used
and on the algorithm used to reduce its dimensionality.

Chemiscope provides an interface that links each structure or environment to a point in
a reduced-dimensionality map and to a 3d-dimensional visualization of its atomic
configuration.

.. figure:: ../img/mol-to-map.*
    :width: 65 %

    Illustration of the process used to create structural properties from a
    molecule.

Chemiscope is completely agnostic with respect to how properties and structural
representations are generated.

There are two main ways to prepare data for chemiscope. One option is to compute
structural descriptors or low-dimensional features using external tools, such as
`scikit-matter`_ or `ASAP`_, and then load them into chemiscope. The other option is to
let chemiscope compute these representations automatically, using
:py:func:`chemiscope.explore`.

.. _ACE: https://doi.org/10.1103/PhysRevB.99.014104
.. _ASAP: https://github.com/BingqingCheng/ASAP
.. _Behler-Parrinello symmetry functions: https://doi.org/10.1103/physrevlett.98.146401
.. _PCovR: https://doi.org/10.1088/2632-2153/aba9ef
.. _PCA: https://en.wikipedia.org/wiki/Principal_component_analysis
.. _PET-MAD: https://arxiv.org/abs/2503.14118
.. _sketch-map: https://doi.org/10.1073/pnas.1108486108
.. _scikit-matter: https://scikit-matter.readthedocs.io/en/latest/
.. _SOAP: https://doi.org/10.1103/PhysRevB.87.184115
