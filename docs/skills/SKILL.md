---
name: chemiscope-python
description: >-
  Author chemiscope visualizations in Python — build interactive structure/property
  explorers for molecules and materials. Use when creating chemiscope.show() /
  write_input() datasets, visualizing forces or vectors as arrows, atom-centered
  (per-environment) properties, custom shapes (spheres, ellipsoids, arrows, meshes),
  configuring map/structure display settings, or embedding chemiscope in Sphinx docs or
  a Streamlit app. Covers ASE / chemfiles / MDAnalysis / stk inputs and the JSON format.
---

# Chemiscope in Python

[Chemiscope](https://chemiscope.org) is an interactive explorer for collections of
molecules and materials. A dataset pairs **structures** (atomic geometries) with
**properties** (scalars/vectors per structure or per atom). Three linked panels stay in
sync: a 2D **property map** (Plotly scatter), a 3D **structure viewer** (3Dmol), and an
**info** panel.

> Requires `pip install chemiscope`. Reading structures needs a backend such as `ase`
> (`pip install ase`), or `chemfiles` / `MDAnalysis` / `stk`.

## Reference files (load on demand)

- **`settings-reference.md`** — full `settings` dict schema (every `map` / `structure`
  field, defaults, accepted values, palettes). Read it when hand-authoring settings
  beyond `quick_settings`.
- **`shapes-reference.md`** — full shapes spec (all kinds, parameters, helper
  signatures). Read it when building custom shapes beyond the cookbook recipes below.

## Choosing the right output

Pick the delivery mechanism before writing code:

| Goal | Use |
|---|---|
| Interactive view in a Jupyter / marimo / Colab notebook | `chemiscope.show(...)` (or `show_input(file)`) |
| Shareable file, open online | `chemiscope.write_input("data.json.gz", ...)` → load at chemiscope.org |
| Single offline HTML file, no server | `write_input("data.json", ...)` (uncompressed!) then `cat chemiscope_standalone.html data.json > out.html` |
| Embed in Sphinx documentation | `chemiscope.sphinx` extension + `.. chemiscope::` directive |
| Embed in a Streamlit app | `chemiscope.streamlit.viewer(dataset, ...)` |
| Scripted PNG images / fully headless | widget `save_structure_image` / `save_map_image`, or `chemiscope.headless` |
| Build a map from learned ML features automatically | `chemiscope.explore(...)` |

There is **no** one-call Python API that emits a self-contained HTML; the standalone
route is "JSON + prebuilt `chemiscope_standalone.html`" (concatenation, uncompressed JSON
only). See "Other delivery mechanisms" at the bottom for Sphinx/Streamlit/headless detail.

`show` / `show_input` / `streamlit.viewer` / the `.. chemiscope::` directive all accept
`mode=` `"default"` (map + structure) | `"structure"` (3D only) | `"map"` (map only).

## Minimal recipe

```python
import ase.io
import chemiscope

structures = ase.io.read("trajectory.xyz", ":")  # ":" reads ALL frames

chemiscope.show(
    structures=structures,
    properties=chemiscope.extract_properties(structures),  # pull props from the files
    metadata={"name": "My dataset"},
    settings=chemiscope.quick_settings(x="energy", y="volume", map_color="energy"),
)
```

Swap `show(...)` for `write_input("out.json.gz", ...)` (same args) to save instead.
`create_input(...)` returns the dataset `dict` without writing — useful for inspection.

## Input structures

`structures` is a list of any supported object (don't mix types):

| Object | Backend | Notes |
|---|---|---|
| `ase.Atoms` | `ase` | Most common. `ase.io.read(path, ":")` reads extxyz, CIF, PDB, POSCAR, LAMMPS dumps, … |
| `MDAnalysis.AtomGroup` | `MDAnalysis` | Each trajectory frame → one structure; good for biomolecules. |
| chemfiles `Frame` / `Trajectory` | `chemfiles` (>=0.10,<0.11) | PDB, trajectories, biomolecule metadata. |
| `stk.Molecule` | `stk` | Bonds included; **no cell, no stored properties**. |
| chemiscope structure `dict` | — | Hand-built (schema below). |

ASE emits a `cell` only for periodic systems (non-zero lattice).

### Hand-built structure dict

```python
{
    "size": 3,                       # required: number of atoms
    "names": ["O", "H", "H"],        # required: per-atom symbol/name (len == size)
    "x": [0.0, 0.76, -0.76],         # required: Cartesian coords (len == size)
    "y": [0.0, 0.59, 0.59],          # required
    "z": [0.0, 0.0, 0.0],            # required
    "cell": [a1,a2,a3, b1,b2,b3, c1,c2,c3],  # optional: flattened 3x3, row-major
    # optional: "bonds": [[i,j,order],...], "elements", "chains", "resnames",
    #           "resids", "hetatom"
}
```

### Large datasets: external structures

```python
external = chemiscope.write_external_structures(structures, prefix="structure")
chemiscope.write_input("dataset.json.gz", structures=external, properties=...)
# Structure files must be shipped alongside the JSON; such a dataset can't be opened
# on chemiscope.org (which needs everything inline).
```

## Properties

A dict keyed by name. Full form (preferred for sharing — units/description matter) or
shorthand (target auto-deduced by length):

```python
properties = {
    "energy": {                      # full form
        "target": "structure",       # "structure" (one per frame) or "atom" (one per atom)
        "values": [1.2, 0.9, 1.5],   # list, 1D np.ndarray, or 2D np.ndarray
        "units": "eV", "description": "DFT energy",
    },
    "volume": [12.0, 12.1, 11.9],    # shorthand
}
```

- A **2D numpy array** (N×M, M>1) becomes M properties `name[1] … name[M]` (PCA
  components, eigenvalues, tensor entries…). Reference them as `x="pca[1]"`.
- `target="atom"` values length = **total atoms across all structures** (or # selected
  environments). `NaN` hides that point on the map.

### Extract from the structures

```python
props = chemiscope.extract_properties(structures)                  # everything
props = chemiscope.extract_properties(structures, only=["energy", "forces"])
props = chemiscope.extract_properties(structures, environments=envs)  # atom props for a subset
```

For ASE: structure props come from `atoms.info`, atom props from `atoms.arrays`;
attached-calculator `energy`/`forces` are injected automatically (so **forces** show up
as an atom property named `"forces"`). Only props present in *every* frame are kept.

### Parameterized (multidimensional) properties

```python
properties = {"spectrum": {"target": "structure", "values": spectra_2d,
                           "parameters": ["freq"]}}
parameters = {"freq": {"values": freq_grid_1d, "units": "cm^-1", "name": "frequency"}}
chemiscope.write_input("out.json.gz", structures=s, properties=properties,
                       parameters=parameters)
```

## Environments (atom-centered properties)

To color/select per atom, declare which atoms carry data and their display cutoff, as
`(structure_id, center_id, cutoff)` tuples:

```python
environments = chemiscope.all_atomic_environments(structures, cutoff=3.5)  # all atoms
environments = [(0, 4, 3.5), (0, 7, 3.5), (2, 0, 4.0)]                     # manual subset
```

Pass `environments=...` **and** use `target="atom"` in the settings for per-atom views.

## Settings (common case)

```python
settings = chemiscope.quick_settings(
    x="energy", y="volume", z=None,        # map axes (property names; use name[i] for components)
    map_color="energy", size=None, symbol=None,
    structure_color=None,                  # color atoms by a property
    target="structure",                    # "atom" for per-environment view
    trajectory=False,                      # True: join map points + keep 3D orientation
    periodic=False,                        # True: show unit cell + 3x3x3 supercell
)
```

For full control (every `map`/`structure` field, palettes, etc.), read
**`settings-reference.md`** and pass a settings dict directly.

## Shapes (common case)

Pass `shapes={name: definition}`, then activate via the structure setting
`"shape": "name"`. Helpers build definitions for you:

```python
chemiscope.ase_vectors_to_arrows(structures, "forces", scale=1, radius=0.15)
chemiscope.ase_tensors_to_ellipsoids(structures, "alpha", force_positive=True, scale=0.2)
```

For custom meshes, batched spheres/cylinders, combined shapes, layered per-atom
parameters, and the full kind list, read **`shapes-reference.md`**.

## Cookbook

### Forces / vectors as arrows

```python
chemiscope.write_input(
    "trajectory.json.gz",
    structures=structures,
    properties=chemiscope.extract_properties(structures),
    shapes={"forces": chemiscope.ase_vectors_to_arrows(structures, "forces", scale=1, radius=0.15)},
    settings={
        "structure": [{"keepOrientation": True, "playbackDelay": 100, "shape": "forces"}],
        "map": {"joinPoints": True},
    },
)
```

Works for any N×3 atom array (velocities, dipoles, displacements…).

### Atom-centered property: color atoms + labels

```python
chemiscope.show(
    structures=structures,
    properties={"charge": {"target": "atom", "values": charges_flat, "units": "e"}},
    environments=chemiscope.all_atomic_environments(structures, cutoff=3.5),
    settings={
        "target": "atom",
        "structure": [{"color": {"property": "charge", "palette": "bwr"},
                       "atomLabels": True, "labelsProperty": "charge"}],
    },
)
```

### Tensors as ellipsoids (e.g. polarizability)

```python
chemiscope.write_input(
    "tensors.json.gz", structures=structures,
    properties=chemiscope.extract_properties(structures, only=["alpha"]),
    shapes={"alpha": chemiscope.ase_tensors_to_ellipsoids(structures, "alpha",
                                                          force_positive=True, scale=0.2)},
    settings={"structure": [{"atoms": False, "shape": "alpha"}]},
    environments=chemiscope.all_atomic_environments(structures),
)
```

### Custom mesh shape, per-atom scaling

```python
cube = {
    "kind": "custom",
    "parameters": {
        "global": {"vertices": [[-1,-1,-1],[1,-1,-1],[-1,1,-1],[1,1,-1],
                                 [-1,-1,1],[1,-1,1],[-1,1,1],[1,1,1]]},  # hull auto-computed
        "atom": [{"scale": 0.3, "color": 0x44AAFF} for _ in range(total_atoms)],
    },
}
chemiscope.write_input("mesh.json.gz", structures=structures, shapes={"cube": cube},
                       settings={"structure": [{"shape": "cube"}]})
```

### Trajectory / periodic / biomolecules

```python
settings = chemiscope.quick_settings(x="time", y="energy", trajectory=True)  # path + fixed view
settings = chemiscope.quick_settings(periodic=True, target="atom")           # cell + supercell

import MDAnalysis as mda
u = mda.Universe("protein.pdb")
chemiscope.show(structures=u.atoms, mode="structure",
                settings=chemiscope.quick_settings(structure_settings={"cartoon": True}))
```

## Metadata

```python
metadata = {"name": "Dataset title",          # required (bare minimum)
            "description": "...", "authors": ["Ada Lovelace"],
            "references": ["doi:10.xxxx/yyyy"]}
```

## Automatic dimensionality reduction

```python
chemiscope.explore(structures, featurizer="pet-mad-1.0",
                   properties=chemiscope.extract_properties(structures),
                   environments=chemiscope.all_atomic_environments(structures))
# write straight to a file: explore(..., write_input="out.json.gz")
# custom model: featurizer = chemiscope.metatomic_featurizer(model="model.pt")
```

## Other delivery mechanisms

### Self-contained offline HTML
```python
chemiscope.write_input("data.json", structures=s, properties=p)  # uncompressed JSON!
```
Then download `chemiscope_standalone.html` from
`https://chemiscope.org/chemiscope_standalone.html` and concatenate:
```bash
cat chemiscope_standalone.html data.json > out.html   # opens in any browser, no server
```
Only uncompressed `.json` (not `.json.gz`) works with the standalone HTML.

### Sphinx documentation
Add the extension in `conf.py`: `extensions = ["chemiscope.sphinx"]`. Then in an `.rst`:
```rst
.. chemiscope:: ../datasets/showcase.json.gz
    :mode: structure       # default | structure | map (optional, default "default")
    :warning_timeout: -1   # ms; 0 = persistent, -1 = disable (optional)
```
The path is to a `.json`/`.json.gz` file relative to the `.rst`. The widget loads data
dynamically, so **serve the built docs over HTTP** (e.g. `python3 -m http.server` in
`build/html`); opening the raw file won't work. (For sphinx-gallery examples, assigning a
`show()` widget to a variable named `___` auto-generates the directive via the bundled
`ChemiscopeScraper`.)

### Streamlit app (`pip install chemiscope[streamlit]`)
```python
import chemiscope, ase.io, streamlit as st
structures = ase.io.read("structures.xyz", ":")
dataset = chemiscope.create_input(structures)
chemiscope.streamlit.viewer(dataset, mode="structure")   # run with: streamlit run app.py
```
`viewer(dataset, *, settings=None, mode="default", no_info_panel=False, width="stretch",
height=550, key="chemiscope_viewer", selected_index=None, on_select=None,
on_settings_change=None)`. `on_select` receives the selected structure index (or `None`);
`selected_index` (e.g. bound to `st.session_state`) drives the selection externally.

### Headless / scripted images
```python
w = chemiscope.show(structures=structures, properties=props, settings=settings)
await w.save_structure_image("frame.png")   # await the futures in a notebook
await w.save_map_image("map.png")
await w.save_structure_sequence([0, 1, 2], ["a.png", "b.png", "c.png"])
# fully headless (no notebook), Playwright-based:  from chemiscope import headless
```

## Gotchas & validation cheat-sheet

- `settings` must be a JSON-serializable **dict** (round-tripped through JSON).
- `target="atom"` / atom-coloring requires `environments`.
- `mode="default"`/`"map"` needs **≥2 properties** for the map.
- Properties-only dataset (no structures) requires every property `target="structure"`.
- Atom-property length = total atoms; structure-property length = number of structures.
- Change a widget's view by reassigning `.settings` with a **complete** dict (in-place
  nested mutation is ignored).
- `frames=` / `meta=` are deprecated aliases for `structures=` / `metadata=`.
- `extract_properties` drops properties not present in *every* structure (with a warning).

## Using this skill elsewhere

This skill is self-contained, portable, and harness-agnostic Markdown. In the chemiscope
repository the canonical copy lives at `docs/skills/` (this `SKILL.md` plus the two
reference files); Claude Code also picks it up via the symlink
`.claude/skills/chemiscope-python -> ../../docs/skills`.

To reuse it in another project or agent:
- **Claude Code** — copy (or symlink) the three files into a named folder under
  `~/.claude/skills/<name>/` (available in every project) or
  `<project>/.claude/skills/<name>/` (scoped to one project).
- **Other agents** — point the agent at these files, e.g. from an `AGENTS.md` pointer or
  by referencing the path directly. The YAML frontmatter is Claude-specific metadata and
  is harmless to other tools.

Either way, `pip install chemiscope` (plus a backend like `ase`) in the working
environment.
