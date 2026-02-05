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
3D structure viewer, streamlining inspection of clusters and outliers by moving between
points in feature space and the corresponding atomic configurations.

Chemiscope 1.0 turns the original browser-based visualizer into a multi-platform tool
that fits into Python-centric workflows. The same visualization can be rendered as a
standalone web viewer, embedded as a Jupyter widget [@Jupyter; @IPython], included in
Streamlit web-applications, or integrated into Sphinx-built documentation and
sphinx-gallery examples for reproducible software manuals [@sphinx]. Chemiscope 1.0 also
provides support for visualizing datasets directly from widely used atomistic Python
toolkits, including ASE [@ase-paper], MDAnalysis [@MDAnalysis], stk [@STK], and
Chemfiles [@chemfiles].

![Overview of Chemiscope 1.0 cross-platform architecture. The Python API accepts
structures from ASE, MDAnalysis, stk, and Chemfiles, along with user-defined properties and
visualization settings. These inputs can be rendered as an interactive Jupyter widget,
embedded in Streamlit applications, integrated into Sphinx documentation, or exported
for the standalone web application at chemiscope.org.](chemiscope-v1.0.svg)

Atomistic modeling workflows produce collections of molecular and materials structures
together with associated quantities, including energies, forces, charges, and other
scalar or tensorial properties. These datasets are commonly explored using
post-processing analysis, including property–property correlations [@Huang2020;
@Wurger2021] and low-dimensional projections [@Helfrecht2020; @Jorgensen2026; @orlov2025; @Tamura2022;
@HernandezLeon2024], to relate abstract representations to the underlying atomic
configurations [@Chapman2022; @Nicholas2020]. Interactive visualization provides a practical
means to interpret structure–property relationships [@Wurger2021], verify computational
results, identify unexpected patterns [@xie2018], and explore learned representations
[@Walsh2025mapping; @Cheng2020; @De2016].

For this purpose, Chemiscope has been adopted across multiple atomistic modeling and coarse-grained studies, with interactive viewers shared alongside publications and archived datasets on
platforms such as Materials Cloud [@Talirz_2020]. While complementary visualization
tools exist, from desktop applications such as VMD and OVITO [@Humphrey1996;
@Stukowski2010] to WebGL-based molecular viewers such as 3Dmol.js and NGLview [@Rego2015;
@Nguyen2017], Chemiscope distinguishes itself by providing a single dataset
representation and rendering stack that can be reused across multiple contexts. This is
especially important in Python-based workflows, where the same visualization is often
needed in a Jupyter notebook for analysis, a web view for sharing, and documentation for
reproducibility and teaching [@JupyterNotebook; @Goscinski2025scicodewidgets; @Du2024].

# Implementation

Chemiscope 1.0 is implemented as a TypeScript visualization library with the Python
package providing platform-specific integrations. The Python API builds the chemiscope
dataset from atomic structures, associated properties, and visualization settings, and
exports it in the JSON schema consumed by the JavaScript renderer. The interface is
organized into linked map, structure, and information panels. The map panel uses Plotly.js to
render 2D and 3D scatter plots [@plotlyjs], the structure panel uses 3Dmol.js for
molecular rendering and supports both atomistic styles and biomolecular cartoons.

The map rendering is a primary performance bottleneck for large datasets. Chemiscope 1.0
introduces adaptive Level of Detail (LOD) rendering for scatter views, which downsample
large datasets based on screen-space density, i.e., how many points would overlap in the
current view. As users zoom or change view parameters, the displayed subset is updated
to preserve both responsiveness and visual structure. In practice, this handles maps with
more than 500,000 points on commodity hardware, without requiring users to pre-filter or
manually decimate their data.

Chemiscope 1.0 can render atom-centered shapes to represent vectorial and tensorial
properties, including arrows (e.g., dipoles or forces), ellipsoids (e.g.,
polarizabilities), and user-defined triangular meshes. For biomolecular systems, it
supports cartoon representations based on residue and chain information. The structure
viewer handles a grid layout for side-by-side comparison of multiple structures or local
environments.

In Jupyter notebooks, the viewer is exposed as a widget with bidirectional communication
between Python and the JavaScript runtime, implemented via traitlets [@Jupyter;
@IPython]. The widget supports programmatic control of the visualization, including
selection synchronization, settings modification, and export of map snapshots. Users can
create a visualization by preparing structures and associated properties and calling
`chemiscope.show`:

```python
import ase.io
import chemiscope

structures = ase.io.read("trajectory.xyz", ":")

# Extract properties present in the trajectory (e.g., energy, forces)
properties = chemiscope.extract_properties(structures)

# Set default settings for multi-frame trajectories
settings = chemiscope.quick_settings(trajectory=True)

# Display the viewer
chemiscope.show(structures=structures, properties=properties, settings=settings)
```

For web applications built with Streamlit, Chemiscope component renders a viewer from an
in-memory dataset and propagates user interactions (e.g., selection and settings changes)
back to Python, coupling to other Streamlit widgets. For reproducible documentation,
Chemiscope includes a Sphinx extension that embeds interactive viewers alongside
narrative text and executable examples [@sphinx].

[LOD plot upon PET-MAD featurizer here ?]

Finally, the package includes `explore` workflow that generates interactive
visualizations starting from structures alone. It integrates metatomic models
[@metatensor], particularly, the PET-MAD model [@Mazitov2025] which is used by default,
to derive informative representations and produce map coordinates without requiring
manual descriptor engineering or an explicit dimensionality reduction step [@MAD]:

```python
chemiscope.explore(structures, featurizer="pet-mad-1.0")
```

Chemiscope is distributed as an open-source package that can be installed from PyPI, and
the default standalone viewer is available online at https://chemiscope.org for quick
inspection of datasets without local setup. Optional features can be installed via
extras: `pip install 'chemiscope[streamlit]'` and `pip install 'chemiscope[explore]'`.

# Acknowledgements

The development of Chemiscope 1.0 has been funded primarily by the [NCCR
MARVEL](http://nccr-marvel.ch/).

# References
