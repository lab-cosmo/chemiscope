"""
Metatrain training dashboard in marimo.

Open it in edit mode so the cells themselves are the Python interpreter:

    marimo edit python/marimo/example/app.py

The cells wire inputs, the viewer, and the section layout. Checkpoint
loading, inference, and the matplotlib figures live in ``dashboard.py``
next to this file.

The sidebar takes a training-run directory, the structure file, and an optional
predictions file. Checkpoints in that directory can be reloaded from the
dropdown. The viewer is the default two-panel layout: a parity plot (reference
energy against the checkpoint prediction) on the left, the structure on the
right. Axis and structure settings sit above it. Training curves come from
``train.csv`` / ``train.log``. Weight norms and the parameter distribution
come from the selected ``model_*.ckpt``.

The defaults point at the MAD subset and the experimental.lorem run under
``test-lorem``.
"""

import marimo


__generated_with = "0.25.1"
app = marimo.App(width="full")


@app.cell
def _():
    import io
    import traceback
    from pathlib import Path

    import dashboard as dash
    import marimo as mo

    import chemiscope

    # Shared bag so the scratch interpreter can see the latest data without
    # subscribing the scratch cell to every widget change.
    state = {}

    def infer(atoms, checkpoint_path=None):
        """Evaluate one ASE structure with the checkpoint selected above."""
        path = checkpoint_path or state.get("checkpoint_path")
        if not path:
            raise ValueError("Select a checkpoint before running inference.")
        return dash.infer(atoms, path)

    state["infer"] = infer
    return Path, chemiscope, dash, infer, io, mo, state, traceback


@app.cell
def _(mo):
    run_dir = mo.ui.text(
        value="/Users/ericboittier/metawork/test-lorem/outputs/2026-10-02/08-05-45",
        label="Training run directory",
        full_width=True,
    )
    structures_path = mo.ui.text(
        value="/Users/ericboittier/metawork/test-lorem/mad_subset.xyz",
        label="Structures (.xyz)",
        full_width=True,
    )
    predictions_path = mo.ui.text(
        value="/Users/ericboittier/metawork/test-lorem/lorem_predictions.xyz",
        label="Predictions (.xyz, optional)",
        full_width=True,
    )
    n_show = mo.ui.slider(
        start=100,
        stop=4709,
        step=100,
        value=400,
        label="Structures on the map",
        show_value=True,
        debounce=True,
        full_width=True,
    )
    mo.sidebar(
        mo.vstack(
            [
                mo.md(
                    """
                    ## Metatrain run

                    Changing a path re-runs the cells that read it.
                    """
                ),
                run_dir,
                structures_path,
                predictions_path,
                n_show,
            ]
        )
    )
    return n_show, predictions_path, run_dir, structures_path


@app.cell
def _(Path, mo, run_dir):
    run = Path(run_dir.value).expanduser()

    def _epoch(path: Path) -> int:
        suffix = path.stem.split("_")[-1]
        return int(suffix) if suffix.isdigit() else 0

    checkpoints = sorted(run.glob("model_*.ckpt"), key=_epoch)
    checkpoint_options = {path.name: str(path) for path in checkpoints} or {
        "(none)": ""
    }
    names = list(checkpoint_options)
    checkpoint = mo.ui.dropdown(checkpoint_options, value=names[-1], label="Checkpoint")
    compare = mo.ui.dropdown(
        checkpoint_options,
        value=names[-2] if len(checkpoints) > 1 else None,
        allow_select_none=True,
        label="Compare weights with",
    )
    mo.hstack([checkpoint, compare], widths="equal", gap=1)
    return checkpoint, compare, run


@app.cell
def _(dash, n_show, predictions_path, state, structures_path):
    _loaded = dash.build_map(
        structures_path.value,
        predictions_path.value,
        int(n_show.value),
    )
    state["frames"] = _loaded["shown"]
    state["properties"] = _loaded["properties"]
    state["indices"] = _loaded["indices"]
    state["prediction_note"] = _loaded["prediction_note"]
    state["n_frames"] = _loaded["n_frames"]
    prediction_note = _loaded["prediction_note"]
    properties = _loaded["properties"]
    property_names = _loaded["property_names"]
    shown = _loaded["shown"]
    return prediction_note, properties, property_names, shown


@app.cell
def _(dash, mo, property_names):
    _numeric, _x_default, _y_default, _color_default = dash.axis_defaults(
        property_names
    )
    x_axis = mo.ui.dropdown(_numeric, value=_x_default, label="Parity x")
    y_axis = mo.ui.dropdown(_numeric, value=_y_default, label="Parity y")
    color_by = mo.ui.dropdown(_numeric, value=_color_default, label="Color")
    palette = mo.ui.dropdown(dash.PALETTES, value="bwr", label="Palette")
    marker_size = mo.ui.slider(
        10,
        80,
        value=32,
        label="Marker size",
        show_value=True,
        debounce=True,
    )
    opacity = mo.ui.slider(
        10,
        100,
        value=90,
        label="Opacity",
        show_value=True,
        debounce=True,
    )
    bonds = mo.ui.checkbox(value=True, label="Bonds")
    space_filling = mo.ui.checkbox(value=False, label="Space filling")
    atom_labels = mo.ui.checkbox(value=False, label="Atom labels")
    controls = mo.vstack(
        [
            mo.hstack(
                [x_axis, y_axis, color_by, palette],
                wrap=True,
                gap=1,
                align="end",
            ),
            mo.hstack(
                [marker_size, opacity, bonds, space_filling, atom_labels],
                wrap=True,
                gap=1,
                align="center",
            ),
        ]
    )
    return (
        atom_labels,
        bonds,
        color_by,
        controls,
        marker_size,
        opacity,
        palette,
        space_filling,
        x_axis,
        y_axis,
    )


@app.cell
def _(
    atom_labels,
    bonds,
    chemiscope,
    color_by,
    controls,
    marker_size,
    mo,
    opacity,
    palette,
    prediction_note,
    properties,
    shown,
    space_filling,
    x_axis,
    y_axis,
):
    settings = chemiscope.quick_settings(
        x=x_axis.value,
        y=y_axis.value,
        map_color=color_by.value,
        structure_settings={
            "bonds": bool(bonds.value),
            "spaceFilling": bool(space_filling.value),
            "atomLabels": bool(atom_labels.value),
            "atoms": True,
        },
        map_settings={
            "color": {
                "property": color_by.value,
                "palette": palette.value,
                "opacity": int(opacity.value),
            },
            "size": {"property": "", "factor": int(marker_size.value)},
        },
    )
    viewer = chemiscope.marimo.viewer(
        shown,
        properties=properties,
        metadata={
            "name": "MAD subset — experimental.lorem",
            "description": prediction_note,
        },
        settings=settings,
        mode="default",
        warning_timeout=-1,
    )
    mo.vstack(
        [
            controls,
            mo.md(f"_{prediction_note}_"),
            viewer,
        ]
    )
    return (viewer,)


@app.cell
def _(mo, state, viewer):
    selected = viewer.selected_ids or {}
    index = selected.get("structure")
    state["selected_index"] = index
    if index is None or "frames" not in state:
        detail = mo.md("Click a point to inspect that structure.")
    else:
        frame = state["frames"][index]
        props = {
            name: spec["values"][index] for name, spec in state["properties"].items()
        }
        source_index = state["indices"][index]
        lines = [
            f"**{frame.get_chemical_formula()}** — map index `{index}`, "
            f"file index `{source_index}`",
            "",
        ]
        for name, value in props.items():
            if isinstance(value, float):
                lines.append(f"- `{name}`: {value:.6g}")
            else:
                lines.append(f"- `{name}`: {value}")
        detail = mo.md("\n".join(lines))
    detail
    return


@app.cell
def _(checkpoint, dash, infer, mo, properties, state, traceback, viewer):
    _pieces = [mo.md("## Inference")]
    _prediction_fig, _summary = dash.prediction_figure(properties)
    if _summary is None:
        _pieces.append(
            mo.md(
                "No predictions file is aligned with the structures, so the "
                "map plots are empty. Click a point to run the checkpoint."
            )
        )
    else:
        _pieces.extend(
            [
                mo.hstack(
                    [
                        mo.stat(
                            f"{_summary['rmse']:.3g}",
                            label="energy RMSE (eV)",
                            bordered=True,
                        ),
                        mo.stat(
                            f"{_summary['mae']:.3g}",
                            label="energy MAE (eV)",
                            bordered=True,
                        ),
                        mo.stat(
                            f"{_summary['rmse_atom']:.0f}",
                            label="RMSE (meV/atom)",
                            bordered=True,
                        ),
                        mo.stat(
                            f"{_summary['bias']:+.3g}",
                            label="mean error (eV)",
                            bordered=True,
                        ),
                    ],
                    gap=1,
                ),
                mo.md(
                    f"{_summary['n']} structures from the predictions file. "
                    "The dashed line on the histogram is the mean error."
                ),
                _prediction_fig,
            ]
        )

    _map_index = (viewer.selected_ids or {}).get("structure")
    if not checkpoint.value:
        _pieces.append(mo.md("Select a checkpoint to evaluate a structure."))
    elif _map_index is None or "frames" not in state:
        _pieces.append(
            mo.md(
                "Click a point for the short-range, long-range, and per-atom "
                "breakdown. `infer(atoms)` runs the same checkpoint from the "
                "Python cell."
            )
        )
    else:
        try:
            _result = infer(state["frames"][_map_index], checkpoint.value)
        except Exception:
            _pieces.append(mo.md(f"```\n{traceback.format_exc()}\n```"))
        else:
            state["inference"] = _result
            _signed = (
                None
                if _result["reference"] is None
                else _result["energy"] - _result["reference"]
            )
            _cards = [
                mo.stat(
                    f"{_result['energy']:.6g}",
                    label="final energy (eV)",
                    caption=f"epoch {_result['epoch']}",
                    bordered=True,
                ),
                mo.stat(
                    f"{_result['sr_energy']:.6g}",
                    label="short-range (eV)",
                    bordered=True,
                ),
                mo.stat(
                    f"{_result['lr_energy']:.6g}",
                    label="long-range (eV)",
                    bordered=True,
                ),
            ]
            if _signed is not None:
                _cards.append(
                    mo.stat(
                        f"{_signed:+.6g}",
                        label="final − reference (eV)",
                        bordered=True,
                    )
                )
            _pieces.extend(
                [
                    mo.md(
                        f"### {_result['formula']}\n\n"
                        f"Checkpoint epoch {_result['epoch']}. The red dashed "
                        "line is the reference energy. Final energy includes "
                        "the composition baseline and the output scale."
                    ),
                    mo.hstack(_cards, gap=1),
                    dash.atom_figure(_result),
                    mo.ui.table(_result["intermediates"], selection=None),
                    mo.ui.table(_result["per_atom"], page_size=8, selection=None),
                ]
            )
    inference_view = mo.vstack(_pieces)
    inference_view
    return


@app.cell
def _(mo):
    matrix_aspect = mo.ui.dropdown(
        {
            "fit": "auto",
            "square cells": "equal",
            "wide cells": 0.5,
            "tall cells": 2.0,
        },
        value="fit",
        label="Aspect",
    )
    matrix_transform = mo.ui.dropdown(
        {
            "linear": "linear",
            "absolute value": "abs",
            "log10(1+|x|)": "log10",
            "signed log": "signed_log",
            "asinh": "asinh",
            "signed square root": "sqrt",
        },
        value="linear",
        label="Transform",
    )
    matrix_controls = mo.hstack(
        [matrix_aspect, matrix_transform],
        wrap=True,
        gap=1,
        align="end",
    )
    return matrix_aspect, matrix_controls, matrix_transform


@app.cell
def _(
    checkpoint,
    dash,
    infer,
    matrix_aspect,
    matrix_controls,
    matrix_transform,
    mo,
    state,
    traceback,
    viewer,
):
    if not checkpoint.value:
        shapes_view = mo.md("## Shapes\n\nSelect a checkpoint.")
    else:
        _loaded, _, _ = dash.read_checkpoint(checkpoint.value)
        _hypers = _loaded["model_data"]["model_hypers"]
        _panels = dash.parameter_panels(_loaded["model_state_dict"], _hypers)
        _transform = matrix_transform.value
        _aspect = matrix_aspect.value
        _shape_fig = dash.parameter_figure(
            _panels, transform=_transform, aspect=_aspect
        )
        _adam_fig, _adam_note = dash.adam_figure(
            _panels,
            _loaded["model_state_dict"],
            _loaded.get("optimizer_state_dict"),
            transform=_transform,
            aspect=_aspect,
        )
        _result = None
        _hidden_fig = None
        _hidden_note = (
            "Click a point to draw the hidden state of that structure. "
            "A training batch uses the same tensors with `n_atoms` equal to "
            "every atom in the batch concatenated, not a padded "
            "`(batch, max_atoms, features)` array."
        )
        _map_index = (viewer.selected_ids or {}).get("structure")
        if _map_index is not None and "frames" in state:
            try:
                _result = infer(state["frames"][_map_index], checkpoint.value)
            except Exception:
                _hidden_note = f"```\n{traceback.format_exc()}\n```"
            else:
                state["inference"] = _result
                _hidden_fig = dash.hidden_figure(
                    _result["hidden"], transform=_transform, aspect=_aspect
                )
                _hidden_note = dash.hidden_caption(_result)
        shapes_view = mo.vstack(
            [
                mo.md(dash.shapes_intro(_hypers)),
                matrix_controls,
                mo.ui.table(dash.shape_ledger(_hypers, _result), selection=None),
                mo.md("### Parameters"),
                _shape_fig,
                mo.md("### Adam\n\n" + _adam_note),
                *([_adam_fig] if _adam_fig is not None else []),
                mo.md("### Hidden state\n\n" + _hidden_note),
                *([_hidden_fig] if _hidden_fig is not None else []),
            ]
        )
    shapes_view
    return


@app.cell
def _(dash, mo, run):
    metrics = dash.read_metrics(run)
    if not metrics:
        training = mo.md(f"No `train.csv` in `{run}`.")
    else:
        _cards = dash.training_headline(metrics)
        training = mo.vstack(
            [
                mo.md("## Training"),
                mo.hstack(
                    [
                        mo.stat(str(_cards["epoch"]), label="epoch", bordered=True),
                        mo.stat(
                            f"{_cards['val_loss']:.3g}",
                            label="val loss",
                            caption=f"train {_cards['train_loss']:.3g}",
                            direction=_cards["loss_direction"],
                            target_direction="decrease",
                            bordered=True,
                        ),
                        mo.stat(
                            f"{_cards['val_rmse']:.0f}",
                            label="val RMSE meV/atom",
                            caption=f"{_cards['delta']:+.0f} vs previous log",
                            direction=_cards["rmse_direction"],
                            target_direction="decrease",
                            bordered=True,
                        ),
                        mo.stat(
                            f"{_cards['learning_rate']:.3e}",
                            label="learning rate",
                            bordered=True,
                        ),
                    ],
                    gap=1,
                ),
                mo.ui.table(metrics),
                dash.training_figure(metrics),
            ]
        )
    training
    return


@app.cell
def _(checkpoint, compare, dash, mo, state):
    if not checkpoint.value:
        weights = mo.md("No checkpoint in this run directory.")
    else:
        _loaded, _weight_rows, _flats = dash.read_checkpoint(checkpoint.value)
        _weight_rows = [dict(row) for row in _weight_rows]
        if compare.value and compare.value != checkpoint.value:
            _, _other_rows, _ = dash.read_checkpoint(compare.value)
            _weight_rows = dash.with_l2_delta(_weight_rows, _other_rows)
        _best = _loaded.get("best_metric")
        state["checkpoint"] = _loaded
        state["checkpoint_path"] = checkpoint.value
        state["weight_rows"] = _weight_rows
        weights = mo.vstack(
            [
                mo.md("## Weights"),
                mo.hstack(
                    [
                        mo.stat(
                            f"{sum(row['numel'] for row in _weight_rows):,}",
                            label="parameters",
                            bordered=True,
                        ),
                        mo.stat(
                            str(_loaded.get("epoch")),
                            label="checkpoint epoch",
                            bordered=True,
                        ),
                        mo.stat(
                            f"{_best:.4g}" if isinstance(_best, float) else "—",
                            label=f"best metric (epoch {_loaded.get('best_epoch')})",
                            bordered=True,
                        ),
                        mo.stat(
                            str(_loaded.get("architecture_name", "—")),
                            label="architecture",
                            bordered=True,
                        ),
                    ],
                    gap=1,
                ),
                dash.weight_figure(_weight_rows, _flats),
                mo.ui.table(_weight_rows, page_size=12, selection=None),
            ]
        )
    weights
    return


@app.cell
def _(dash, mo, run):
    _log_path = run / "train.log"
    _log_text = _log_path.read_text() if _log_path.is_file() else ""
    _blocks = [
        mo.stat(
            row["count"],
            label=f"{row['split'].lower()} structures",
            caption=f"energy {row['mean']:.2f} ± {row['std']:.1f} eV",
            bordered=True,
        )
        for row in dash.dataset_summaries(_log_text)
    ]
    _options_path = run / "options_restart.yaml"
    _options_text = _options_path.read_text() if _options_path.is_file() else ""
    _summary = (
        mo.hstack(_blocks, gap=1)
        if _blocks
        else mo.md("No dataset summary in the log.")
    )
    mo.vstack(
        [
            mo.md("## Run log"),
            _summary,
            mo.accordion(
                {
                    "train.log": mo.md(f"```\n{_log_text[-4000:]}\n```"),
                    "options_restart.yaml": mo.md(f"```yaml\n{_options_text}\n```"),
                }
            ),
        ]
    )
    return


@app.cell
def _(mo):
    scratch = mo.ui.code_editor(
        value=(
            "index = state.get('selected_index')\n"
            "frame = state['frames'][0 if index is None else index]\n"
            "result = infer(frame)  # or infer(some_atoms, checkpoint_path)\n"
            "print(result['formula'], result['n_atoms'], 'atoms')\n"
            "print('final energy', result['energy'])\n"
            "for row in result['intermediates']:\n"
            "    print(f\"{row['name']}: {row['value']}\")\n"
        ),
        language="python",
        min_height=180,
        debounce=True,
        label="Python",
    )
    run_scratch = mo.ui.run_button(label="Run this Python")
    return run_scratch, scratch


@app.cell
def _(dash, infer, io, mo, run_scratch, scratch, state, traceback):
    import contextlib

    output = mo.md(
        "Every cell in this notebook is Python. `infer(atoms)` evaluates the "
        "selected checkpoint. `dash` is the helper module (figures, checkpoint "
        "loading, the map sample). `state` holds `frames`, `properties`, "
        "`indices`, `checkpoint`, `checkpoint_path`, `weight_rows`, and "
        "`inference`."
    )
    if run_scratch.value:
        buffer = io.StringIO()
        namespace = {"state": state, "mo": mo, "infer": infer, "dash": dash}
        try:
            with contextlib.redirect_stdout(buffer):
                exec(compile(scratch.value, "<scratch>", "exec"), namespace)
            printed = buffer.getvalue() or "(no output)"
            output = mo.md(f"```\n{printed}\n```")
        except Exception:
            output = mo.md(f"```\n{traceback.format_exc()}\n```")
    mo.vstack([mo.md("## Python"), scratch, run_scratch, output])
    return


if __name__ == "__main__":
    app.run()
