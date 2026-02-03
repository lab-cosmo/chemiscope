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
Chemiscope is an interactive visualization viewer for exploring structure–property
relationships in molecular and materials datasets [@Fraux2020]. The original project
provided a lightweight, browser-based interface to link embedding maps with 3D
renderings of structures and atomic environments. With the _1.0.0_ release the project matures into a cross-platform visualizer that embeds the interactive viewer into Jupyter notebooks [@Jupyter], web application frameworks (Streamlit), and Sphinx-generated documentation [@sphinx], while preserving the standalone web viewer and format. The release also extends support for common atomistic Python toolkits [@ase-paper; @MDAnalysis; @STK] and scales interactive exploration to datasets containing hundreds of thousands of points.

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

Chemiscope package can be installed via pip (`pip install chemiscope`), with optional dependencies for the explore functionality (`pip install chemiscope[explore]`) and Streamlit support (`pip install chemiscope[streamlit]`). The standalone web application is accessible at https://chemiscope.org.

TODO: say about materials cloud intergration?

# Acknowledgements

The development of chemiscope 1.0 has been funded primarily by 
the [NCCR MARVEL](http://nccr-marvel.ch/). 

# References
