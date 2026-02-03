---
title: 'Chemiscope 1.0: a multi-platform atomistic data explorer'
tags:
    - TypeScript
    - JavaScript
    - Python
    - chemistry
    - material science
    - machine learning
    - visualization
authors:
    - name: Sofiia Chorna
      orcid: 0009-0008-7426-0856
      affiliation: 1
    - name: Jakub Lala
      orcid: 0000-0002-5424-5260
      affiliation: "1, 2"
    - name: Qianjun Xu
      orcid: 0000-0003-0778-7208
      affiliation: 1
    - name: Rose K. Cersonsky
      orcid: 0000-0003-4515-3441
      affiliation: 3
      # TODO: confirm affiliation
    - name: Guillaume Fraux
      orcid: 0000-0003-4824-6512
      affiliation: 1
    - name: Michele Ceriotti
      orcid: 0000-0003-2571-2832
      affiliation: 1
affiliations:
    - name: Laboratory of Computational Science and Modeling, IMX, École Polytechnique Fédérale de Lausanne, 1015 Lausanne, Switzerland
      index: 1
    - name: Department of Materials, Imperial College London, London SW7 2AZ, United Kingdom
      index: 2
    - name: Department of Chemical and Biological Engineering, University of Wisconsin-Madison, Madison, WI 53705, United States
      index: 3
date: 31 January 2026
bibliography: paper.bib
---

# Summary
Chemiscope is an interactive visualization tool for exploring structure–property
relationships in molecular and materials datasets [@Fraux2020]. It links a map view,
e.g., a low-dimensional embedding or property–property scatter plot, to an interactive
3D structure viewer, enabling rapid inspection of clusters and outliers by moving
between points in feature space and the corresponding atomic configurations.

Chemiscope 1.0 turns the original browser-based visualizer into a multi-platform tool
that fits into Python-centric workflows. The same visualization can be rendered as a
standalone web viewer, embedded as a Jupyter widget [@Jupyter; @IPython], included in
Streamlit web-applications, or integrated into Sphinx-built documentation and
sphinx-gallery examples for reproducible teaching and software manuals [@sphinx].
Chemiscope 1.0 also provides support for visualizing datasets directly from widely used
atomistic Python toolkits, including ASE [@ase-paper], MDAnalysis [@MDAnalysis], and stk
[@STK].

![Overview of chemiscope v1.0 cross-platform architecture. The Python API accepts structures from ASE, MDAnalysis, and stk, along with user-defined properties and visualization settings. These inputs can be rendered as an interactive Jupyter widget, embedded in Streamlit applications, integrated into Sphinx documentation, or exported for the standalone web application at chemiscope.org.](chemiscope-v1.0.0.png)

Computational chemistry and materials science increasingly rely on the large datasets that combine atomic structures, local environments, computed properties and often visualized through descriptors (TODO: citations). Identifying correlations, outliers, and structure–property relationships in such data often requires simultaneous inspection of reduced-dimensional representations, scalar or vector properties, and the underlying geometry. While the original chemiscope addressed the need for such visualization [@Fraux2020], modern computational workflows increasingly rely on Python environments for data analysis.

Recent chemiscope 1.0 addresses this gap through several advancements. Native Jupyter integration provides interactive widgets, supporting bidirectional communication between the visualization and the notebook kernel. Streamlit component enable building interactive web applications, while Sphinx extensions allow direct embedding of chemiscope viewers in documentation with sphinx-gallery integration for executable examples [@sphinx]. A unified Python API provides a consistent interface for creating visualizations across ASE [@ase-paper], MDAnalysis [@MDAnalysis], and stk structures [@STK]. Finally, built-in support for machine learning featurizers, including integration with metatomic models [@metatensor], enables automated dimensionality reduction and exploration of chemical space without manual descriptor engineering.

# Implementation
**Very drafty**

Chemiscope 1.0 maintains its core implementation in TypeScript and Python.

To handle datasets with hundreds of thousands of data points while maintaining interactive frame rates, chemiscope 1.0 implements adaptive Level of Detail (LOD) rendering.

The Jupyter widget implementation uses the ipywidgets framework to establish bidirectional communication between Python and JavaScript [@IPython]. Users can create visualizations with a single function call:

```python
import chemiscope
import ase.io

structures = ase.io.read("trajectory.xyz", ":")
chemiscope.show(structures, properties={"energy": energies})
```

The widget supports programmatic control of the visualization, including selection synchronization, settings modification, and export of map snapshots.

Chemiscope 1.0 introduces support for visualizing atom-centered shapes, enabling the
display of tensorial and vectorial properties directly on atomic structures. For
example, users can render ellipsoids to represent polarizability tensors or atomic
displacement parameters, arrows to visualize forces or dipole moments, and custom shapes defined by arbitrary triangular meshes.

For biomolecular systems, chemiscope now supports cartoon representations ...

The structure panel grid layout now enables simultaneous comparison of up to nine different structures or atomic environments. 

Finally, the `chemiscope.explore()` function provides automatic featurization and
dimensionality reduction based on last-layer features of PET-MAD universal interatomic potential [@Mazitov2025] projected into reduced latent feature space upon MAD dataset [@MAD],  without requiring user-specified descriptors or dimensionality reduction parameters.

```python
chemiscope.explore(structures, featurizer="pet-mad-1.0")
```

Chemiscope datasets can be exported as standalone artifacts that remain viewable in a
browser, supporting dissemination alongside publications and long-term access to
interactive figures. Integration with the Materials Cloud is an additional route to
archive computational data and make interactive visualization available directly from
the record page [@Talirz_2020]. (Sofiia: should we cite some papers that do this?)

Chemiscope is distributed as an open-source package that can be installed from PyPI, and
the default standalone viewer is available online at https://chemiscope.org for quick
inspection of datasets without local setup.

# Acknowledgements

The development of chemiscope 1.0 has been funded primarily by 
the [NCCR MARVEL](http://nccr-marvel.ch/). 

# References
