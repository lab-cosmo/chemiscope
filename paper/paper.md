---
title: 'Chemiscope: interactive structure-property explorer for materials and molecules'
tags:
  - TypeScript
  - JavaScript
  - chemistry
  - material science
  - machine learning
authors:
  - name: Guillaume Fraux
    orcid: 0000-0003-4824-6512
    affiliation: 1
  - name: Rose K. Cersonsky
    orcid: 0000-0003-4515-3441
    affiliation: 1
  - name: Michele Ceriotti
    orcid: 0000-0003-2571-2832
    affiliation: 1
affiliations:
 - name: Laboratory of Computational Science and Modeling, IMX, École Polytechnique Fédérale de Lausanne, 1015 Lausanne, Switzerland
   index: 1
date: 30 January 2020
bibliography: paper.bib
---

# Summary

The number of materials or molecules that can be created by combining different
chemical elements in various proportions and spatial arrangements is enormous.
Computational chemistry can be used to generate databases containing billions of
potential structures [@Ruddigkeit2012], and predict some of the associated
properties [@Montavon2013; @Ramakrishnan2014]. Unfortunately, the very large
number of structures makes exploring such database — to understand
structure-property relations or find the *best* structure for a given
application — a daunting task. In recent years, multiple molecular
*representations* [@Behler2007; @Bartok2013; @Willatt2019] have been developed
to compute structural similarities between materials or molecules, incorporating
physically-relevant information and symmetries. The features associated with
these representations can be used for unsupervised machine learning
applications, such as clustering or classification of the different structures,
and high-throughput screening of database for specific properties [@Maier2007;
@De2017; @Hautier2019]. Unfortunately, the dimensionality of these features (as
well as most of other descriptors used in chemical and materials informatics) is
very high, which makes the resulting classifications, clustering or mapping very
hard to visualize. Dimensionality reduction algorithms [@Schlkopf1998;
@Ceriotti2011; @McInnes2018] can reduce the number of relevant dimensions to a
handful, creating 2D or 3D maps of the full database.

![The Qm7b database [@Montavon2013] visualized with chemiscope](screenshot.png)

Chemiscope is a graphical tool for the interactive exploration of materials and
molecular databases, correlating local and global structural descriptors with
the physical properties of the different systems. The interface consists of
two panels. The left panel displays a 2D or 3D scatter plot, in which each
point corresponds to a chemical entity. The axes, color, and style of each point
can be set to represent a property or a structural descriptor to visualize
structure-property relations directly. Structural descriptors are not computed
directly by chemiscope, but must be obtained from one of the many codes
implementing general-purpose atomic representation [@librascal; @QUIP] or more specialized descriptors. Since the most common
descriptors can be very high dimensional, it can be convenient to apply a
dimensionality reduction algorithm that maps them to a lower-dimensional space
for easier visualization. For example the sketch-map algorithm [@Ceriotti2011]
was used with the Smooth Overlap of Atomic Positions representation [@Bartok2013] to
generate the visualization in Figure 1. The right panel displays the
three-dimensional structure of the chemical entities, possibly including
periodic repetition for crystals. Visualizing the chemical structure can help
in finding an intuitive rationalization of the layout of the dataset and the
structure-property relations.

Whereas similar tools [@Gong2013; @Gutlein2014; @Probst2017; @ISV] only allow
visualizing maps and structures in which each data point corresponds to a
molecule, or a crystal  structure, a distinctive feature of chemiscope is the
possibility of visualizing maps in which points correspond to atom-centred
environments. This is useful, for instance, to rationalize the relationship
between structure and atomic properties such as nuclear chemical shieldings
(Figure 2). This is also useful as a diagnostic tool for the many
machine-learning schemes that decompose properties into atom-centred
contributions [@Behler2007; @Bartok2010].


![Database of chemical shieldings [@Paruzzo2018] in chemiscope demonstrating the use of a 3D plot and highlighting of atomic environments](./screenshot-3d.png)

Chemiscope took strong inspiration from a previous similar graphical software,
the interactive sketch-map visualizer [@ISV]. This previous software was used in
multiple research publication, related to the exploration of large-scale
databases, and the mapping of structure-property relationships [@De2016;
@De2017; @Musil2018].

# Implementation

Chemiscope is implemented using the web platform: HTML5, CSS and WebGL to
display graphical elements, and TypeScript (compiled to JavaScript) for
interactivity. It uses [Plotly.js](https://plot.ly/javascript/) to render and
animate 2D and 3D plots; and the JavaScript version of [Jmol](http://jmol.org/)
to display atomic structures. The visualization is fast enough to be used with
datasets containing up to a million points, reacting to user input within a few
hundred milliseconds in the default 2D mode. More elaborate visualizations are
slower, while still handling 100k points easily.

The use of web technologies makes chemiscope usable from different operating
systems without the need to develop, maintain and package the code for each
operating system. It also means that we can provide an online service at
http://chemiscope.org that allows users to visualize their own dataset without any
local installation. Chemiscope is implemented as a library of re-usable
components linked together via callbacks. This makes it easy to modify the
default interface to generate more elaborate visualizations, for example,
displaying multiple maps generated with different parameters of a dimensionality
reduction algorithm. Chemiscope can also be distributed in a standalone mode,
where the code and a predefined dataset are merged together as a single HTML
file. This standalone mode is useful for archival purposes, for example as
supplementary information for a published article and for use in corporate
environments with sensitive datasets.

# Acknowledgements

The development of chemiscope have been funded by the [NCCR
MARVEL](http://nccr-marvel.ch/), the [MAX](http://max-centre.eu/) European
centre of excellence, and the European Research Council (Horizon 2020 grant
agreement no. 677013-HBMAP).

# References
