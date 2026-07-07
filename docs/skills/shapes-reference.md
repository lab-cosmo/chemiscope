# Chemiscope shapes reference

Read this when hand-authoring custom shapes or when you need the full list of shape
kinds and their parameters. Common cases (forces as arrows, tensors as ellipsoids,
a custom mesh) are shown in SKILL.md's cookbook; this file is the exhaustive spec.

Shapes are extra geometry drawn in the 3D viewer. Pass `shapes={name: definition}` to
`show` / `write_input` / `create_input`, then **activate** a shape via the structure
setting `"shape": "name"` (comma-separated for several, e.g. `"forces,dipole"`).

## Layered parameters

```python
shapes = {
    "my_shape": {
        "kind": "arrow",
        "parameters": {
            "global": {...},      # applies to every instance
            "structure": [ {...}, ... ],  # optional: one dict per structure
            "atom": [ {...}, ... ],       # optional: flat list over ALL atoms
        },
    }
}
```

Merge order for atom *k* in structure *j*: `global` ← `structure[j]` ← `atom[k]`
(more specific overrides). Lengths must match exactly: `structure` == number of
structures; `atom` == total number of atoms across all structures.

## General parameters (any kind)

- `position`: `[x, y, z]` — defaults to the origin (structure shapes) or the atom
  position (atom shapes).
- `scale`: float — uniform size scaling.
- `orientation`: quaternion `[x, y, z, w]`.
- `color`: hex int `0xRRGGBB` or string `"#RRGGBBAA"`.

## Kinds and their specific parameters

| kind | specific params | notes |
|---|---|---|
| `sphere` | `radius` | no orientation |
| `ellipsoid` | `semiaxes` `[ax, ay, az]` | |
| `cylinder` | `vector` `[x,y,z]`, `radius` | `orientation` ignored (vector sets it) |
| `arrow` | `vector` `[x,y,z]`, `baseRadius`, `headRadius`, `headLength` | `orientation` ignored |
| `spheres` | `centers` (N×3), `radii` (scalar or N, default 1.0), `colors` (one or N) | batched |
| `cylinders` | `vectors` (N×3), `bases` (N×3, optional), `radii` (scalar or N, default 0.1), `colors` | batched |
| `custom` | `vertices` (N×3), `simplices` (M×3 int triangles, optional) | if `simplices` omitted, the convex hull is computed (needs `scipy`) |
| `combined` | `shapes`: list of sub-shape definitions | sub-shapes use the standard `parameters` format; **no nested `combined`** |

`combined` example:

```python
{"kind": "combined", "shapes": [
    {"kind": "cylinders", "parameters": {...}},
    {"kind": "spheres",   "parameters": {...}},
]}
```

## Helper functions

Build single parameter dicts (use `None` for a param to leave it for `global`):

```python
chemiscope.arrow_from_vector([1, 0, 0], scale=1.0, radius=0.1,
                             head_radius_scale=1.75, head_length_scale=2.0)
# -> {"position": [...], "baseRadius": ..., "headRadius": ..., "headLength": ...}

chemiscope.ellipsoid_from_tensor(tensor, scale=1.0, force_positive=False)
# tensor: 3x3 or 6-array [xx, yy, zz, xy, xz, yz]; needs scipy
# -> {"semiaxes": [...], "orientation": [x, y, z, w]}
```

Build a whole shape definition straight from an ASE property:

```python
chemiscope.ase_vectors_to_arrows(structures, key="forces", target=None,  # auto-detect
                                 scale=1.0, radius=0.1,
                                 head_radius_scale=1.75, head_length_scale=2.0)
# key must be an N×3 atom (or structure) property

chemiscope.ase_tensors_to_ellipsoids(structures, key="alpha", target=None,
                                     scale=1.0, force_positive=False)
# key must be a 6- or 9-component tensor property
```

`scale` converts from the property's units to the position units (usually Å). For
non-positive-definite tensors pass `force_positive=True`.

## Validation reminders

- `shapes` is a dict of `{name: {"kind": ..., "parameters": ...}}`.
- For `combined`, the only allowed top-level keys are `kind` and `shapes`; for all other
  kinds they are `kind` and `parameters`.
- Unknown kinds/keys/params raise `ValueError`/`TypeError`.
```
