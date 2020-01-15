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

The number of different materials or molecules that can created by combining
different chemical elements in various proportions and spatial arrangements is
enormous. We can use computational chemistry to generate databases containing up
to millions of potential structures, and predict some of the associated
properties. Unfortunately, the very large number of structures makes exploring
such database — to understand structure-property relations or find the *best*
structure for a given application — a daunting task. In the recent years,
multiple molecular *descriptors* [@Willatt2019; @Behler2007] have been developed
to compute structural similarities between materials or molecules, incorporating
physically-relevant information and symmetries. These descriptors can be used
for unsupervised machine learning applications, such as clustering or
classification of the different structures, and high-throughput screening of
database for specific properties [@Hautier2019]. Unfortunately, the
dimensionality of most descriptors is very high, which makes the resulting
classifications, clustering or mapping very hard to visualize. Additional
dimensionality reduction algorithm can reduce the number of relevant dimensions
to a handful, creating 2D or 3D maps of the full database.

![The Qm7b database [@Montavon2013] visualized with chemiscope](screenshot.png)

Chemiscope is an graphical tool for the interactive exploration of materials and
molecular databases, correlating local and global structural descriptors with
the physical properties of the different systems. Structural properties are
represented by by a descriptor mapped onto a smaller sub-space using a
dimensionality reduction algorithm, and the structure themselves are shown in
relation with their position in this sub-space. For example the sketchmap
algorithm [@Ceriotti2011] was used with the Smooth Overlap of Atomic Positions
descriptor [@Bartok2013] to generate the visualization above.

Compared to similar tools (@Gong2013; @Gutlein2014; @Probst2017; @ISV)
chemiscope supports displaying three dimensional molecular systems and spherical
atomic environments (see the second example below) in the molecular viewer. Such
spherical atomic environments are at the basis of a lot of different structural
descriptors. Chemiscope is also able to display crystalline super-cells and
non-orthorhombic unit cells, which are especially important when working with
databases of materials.

![Database of chemical shielding [@Paruzzo2018] in chemiscope showing the use of a 3D plot and atomic environments highlighting](./screenshot-3d.png)

Chemiscope took strong inspiration from a previous similar graphical software,
the interactive sketchmap visualizer [@ISV]. This previous software was used in
multiple research publication, related to the exploration of large-scale
databases, and the mapping of structure-property relationships (@De2016;
@De2017; @Musil2018).

# Implementation

Chemiscope is implemented using the web platform: HTML5, CSS and WebGL to
display graphical elements, and TypeScript (compiled to JavaScript) for the
interactivity. It uses [Plotly.js](https://plot.ly/javascript/) to render and
animate 2D and 3D plots; and the JavaScript version of [Jmol](http://jmol.org/)
to display atomic structures. The web platform makes chemiscope usable from a
lot of different operating systems without requiring porting, packaging and
distributing for each separate operating system. It also make it available to
users without installation, since everyone can visualize their own datasets
online, at http://chemiscope.org.

Chemiscope is implemented as library of re-usable components linked together via
callbacks. This makes creating new visualizations easy: for example using
parametric dimensionality reduction algorithms could re-use the molecular viewer
and have multiple maps. The visualization is fast enough to be used with
datasets containing up to a million points, reacting to user input within a few
hundred milliseconds in the default 2D mode. More elaborate visualization are
slower, while still handling 100k points easily.

Chemiscope also provides an offline mode, where it uses a single, standalone
HTML file to display and visualize the dataset. This standalone mode is useful
for archival, for example in as supplementary information for a published
article; and for use in corporate environments with sensitive datasets.

# Acknowledgements

The development of chemiscope have been funded by the [NCCR
MARVEL](http://nccr-marvel.ch/) and the [MAX](http://max-centre.eu/) European
center of excellence.

# References
