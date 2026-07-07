# Chemiscope settings reference

Read this when hand-authoring a full `settings` dict (beyond what `quick_settings`
covers) or when you need exact field names, defaults, and accepted values.

`settings` must be a JSON-serializable dict. It has four top-level keys:

| key | type | meaning |
|---|---|---|
| `map` | object | the 2D property-map panel (single object) |
| `structure` | **list** of objects | one object per 3D viewer panel |
| `pinned` | list[int] | environment indices pinned into viewers; default `[0]` |
| `target` | `"structure"` \| `"atom"` | what selections refer to; `"atom"` requires `environments` |

## `map`

```python
"map": {
    # axes — each of x / y / z is an object
    "x": {"property": "pca[1]", "scale": "linear",   # scale: "linear" | "log"
          "min": float("nan"), "max": float("nan")}, # NaN = autoscale
    "y": {"property": "pca[2]"},
    "z": {"property": ""},                # "" = no z axis → 2D map (default)

    "color": {
        "property": "energy",             # "" = fixed color
        "mode": "linear",                 # "linear" | "log" | "sqrt" | "inverse"
        "min": float("nan"), "max": float("nan"),
        "palette": "inferno",
        "opacity": 100,                   # 1..100
    },
    "size": {
        "property": "",                   # "" = fixed size
        "factor": 50,                     # 1..100
        "mode": "linear",                 # linear|log|sqrt|inverse|flip-linear|proportional
    },
    "symbol": "",                         # name of a categorical property, or ""
    "joinPoints": False,                  # draw a line through points (trajectories)
    "markerOutline": True,
    "useLOD": True,                       # level-of-detail rendering for >50k points
    # "camera": {"eye": ..., "center": ..., "up": ..., "zoom": >0}  # optional
}
```

## `structure` (list — one dict per viewer)

```python
"structure": [{
    "bonds": True,
    "atoms": True,
    "spaceFilling": False,
    "atomLabels": False,
    "labelsProperty": "element",          # which atom property to label with
    "unitCell": False,
    "supercell": [1, 1, 1],               # positive ints, repeats along a/b/c
    "axes": "off",                        # "off" | "abc" | "xyz"
    "keepOrientation": False,             # keep camera across frames (trajectories)
    "playbackDelay": 700,                 # ms between frames
    "rotation": False,                    # auto-spin
    "cartoon": False,                     # cartoon representation (biomolecules)
    "shape": "",                          # shape name(s) to display, comma-separated
    "color": {                            # color atoms by a property
        "property": "element",
        "transform": "linear",            # "linear" | "log" | "sqrt" | "inverse"
        "min": 0, "max": 0,
        "palette": "bwr",
    },
    "environments": {
        "activated": True,
        "center": False,                  # auto-center the environment
        "cutoff": 4.0,                    # >= 0
        "bgStyle": "ball-stick",          # "ball-stick" | "licorice" | "cartoon" | "hide"
        "bgColor": "grey",                # "grey" | "CPK" | "property"
    },
    # "camera": {"eye": ..., "center": ..., "up": ..., "zoom": >0}  # optional
}]
```

## Palettes

Shared by `map.color.palette` and `structure[].color.palette`:

`inferno`, `magma`, `plasma`, `viridis`, `cividis`, `seismic`, `brg`, `bwr`, `rwg`,
`twilight (periodic)`, `twilight dark (periodic)`, `hsv (periodic)`, `tab10`, `tab20`,
`tab20b` (the `tab*` ones are categorical).

## Notes & legacy quirks

- The TypeScript validators are the source of truth for accepted values. Some older docs
  show values the current code rewrites: `color.scale` → `color.mode`; a top-level
  `settings.palette` → `color.palette`; size `mode: "constant"` → `property: ""`;
  `select.mode` `"range"`/`"category"` → `"range-gray"`/`"category-gray"`.
- A viewer reassigns its whole `settings`; to change a widget's view at runtime assign a
  **complete** dict to `widget.settings` (in-place mutation of a nested key is ignored).
- `target="atom"` requires `environments` in the dataset.
```
