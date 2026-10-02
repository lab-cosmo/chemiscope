"""
Chemiscope inside marimo.

Loads the showcase structures shipped with the Python examples and displays them
with :func:`chemiscope.marimo.viewer`. A second cell reads ``selected_ids``, so it
re-runs when a point or structure is selected.

Run from the repository root:

    marimo run python/marimo/example/app.py
"""

import marimo


__generated_with = "0.25.1"
app = marimo.App(width="full")


@app.cell
def _():
    from pathlib import Path

    import ase.io
    import marimo as mo

    import chemiscope

    data_path = (
        Path(__file__).resolve().parents[2] / "examples" / "data" / "showcase.xyz"
    )
    structures = ase.io.read(data_path, ":")
    properties = chemiscope.extract_properties(
        structures, only=["dipole_ccsd", "ccsd_pol"]
    )
    return chemiscope, mo, properties, structures


@app.cell
def _(chemiscope, mo, properties, structures):
    viewer = chemiscope.marimo.viewer(
        structures,
        properties=properties,
        metadata={"name": "Dipole and polarizability"},
        settings=chemiscope.quick_settings(
            x="ccsd_pol[1]",
            y="ccsd_pol[2]",
            map_color="dipole_ccsd[1]",
        ),
    )
    mo.vstack([mo.md("# Chemiscope inside marimo"), viewer])
    return (viewer,)


@app.cell
def _(mo, viewer):
    selected = viewer.selected_ids
    mo.md(
        f"""
        **Selection** — click a point on the map or a structure in the viewer.

        `{selected}`
        """
    )
    return


if __name__ == "__main__":
    app.run()
