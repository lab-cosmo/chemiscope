"""
Metatrain training dashboard in marimo.

Open it in edit mode so the cells themselves are the Python interpreter:

    marimo edit python/marimo/example/app.py

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
    import functools
    import io
    import re
    import traceback
    from pathlib import Path

    import ase.io
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np

    import chemiscope

    # Shared bag so the scratch interpreter can see the latest data without
    # subscribing the scratch cell to every widget change.
    state = {}

    @functools.lru_cache(maxsize=2)
    def read_frames(path: str):
        return list(ase.io.read(path, index=":"))

    def comment_values(path: str, key: str) -> list[str]:
        values = []
        token = key + "="
        pattern = re.compile(rf"(?:^|\s){re.escape(key)}=(\S+)")
        with open(path) as handle:
            for line in handle:
                if token not in line or "Properties=" not in line:
                    continue
                match = pattern.search(line)
                values.append(match.group(1) if match else "")
        return values

    @functools.lru_cache(maxsize=4)
    def read_checkpoint(path: str):
        import metatomic.torch  # noqa: F401
        import torch

        checkpoint = torch.load(path, map_location="cpu", weights_only=False)
        rows = []
        flats = []
        for name, tensor in checkpoint["model_state_dict"].items():
            flat = tensor.detach().float().reshape(-1).cpu()
            rows.append(
                {
                    "name": name,
                    "shape": "×".join(str(size) for size in tensor.shape) or "scalar",
                    "numel": int(flat.numel()),
                    "mean": float(flat.mean()),
                    "std": float(flat.std(unbiased=False)),
                    "min": float(flat.min()),
                    "max": float(flat.max()),
                    "l2": float(torch.linalg.vector_norm(flat)),
                }
            )
            flats.append(flat.numpy())
        return checkpoint, rows, flats

    def read_metrics(path: Path):
        csv_path = path / "train.csv"
        if not csv_path.is_file():
            return []
        parsed = []
        for line in csv_path.read_text().splitlines()[1:]:
            if not line.strip() or line.startswith(","):
                continue
            epoch, lr, train_loss, train_rmse, val_loss, val_rmse = line.split(",")
            parsed.append(
                {
                    "epoch": int(epoch),
                    "learning_rate": float(lr),
                    "train_loss": float(train_loss),
                    "train_rmse_meV": float(train_rmse),
                    "val_loss": float(val_loss),
                    "val_rmse_meV": float(val_rmse),
                }
            )
        return parsed

    return (
        Path,
        chemiscope,
        comment_values,
        io,
        mo,
        np,
        plt,
        re,
        read_checkpoint,
        read_frames,
        read_metrics,
        state,
        traceback,
    )


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
def _(
    comment_values,
    n_show,
    np,
    predictions_path,
    read_frames,
    state,
    structures_path,
):
    frames = read_frames(structures_path.value)
    reference = [
        float(value)
        for value in comment_values(structures_path.value, "ecumetric_energy")
    ]
    groups = comment_values(structures_path.value, "dataset_group")
    predicted = []
    prediction_note = "No predictions file."
    pred_path = predictions_path.value.strip()
    if pred_path:
        predicted_text = comment_values(pred_path, "energy")
        if len(predicted_text) == len(frames):
            predicted = [float(value) for value in predicted_text]
            prediction_note = f"Predictions aligned with {len(predicted)} structures."
        else:
            prediction_note = (
                f"Predictions file has {len(predicted_text)} frames, "
                f"structures file has {len(frames)}. Error column omitted."
            )

    take = min(int(n_show.value), len(frames))
    chosen = np.linspace(0, len(frames) - 1, take, dtype=int)
    shown = [frames[int(i)] for i in chosen]
    energies = [reference[int(i)] for i in chosen]
    per_atom = [
        energy / max(len(frame), 1)
        for energy, frame in zip(energies, shown, strict=True)
    ]
    n_atoms = [len(frame) for frame in shown]
    dataset_group = [groups[int(i)] if groups else "" for i in chosen]

    properties = {
        "energy": {
            "target": "structure",
            "values": energies,
            "units": "eV",
            "description": "MAD ecumetric_energy",
        },
        "energy_per_atom": {
            "target": "structure",
            "values": per_atom,
            "units": "eV/atom",
            "description": "ecumetric_energy divided by the number of atoms",
        },
        "n_atoms": {
            "target": "structure",
            "values": n_atoms,
            "description": "Number of atoms",
        },
    }
    if groups:
        properties["dataset_group"] = {
            "target": "structure",
            "values": dataset_group,
            "description": "MAD dataset_group",
        }
    if predicted:
        pred = [predicted[int(i)] for i in chosen]
        properties["predicted_energy"] = {
            "target": "structure",
            "values": pred,
            "units": "eV",
            "description": "Energy written in the predictions file",
        }
        properties["error"] = {
            "target": "structure",
            "values": [
                left - right for left, right in zip(pred, energies, strict=True)
            ],
            "units": "eV",
            "description": "predicted_energy - ecumetric_energy",
        }

    property_names = list(properties)
    state["frames"] = shown
    state["properties"] = properties
    state["indices"] = [int(i) for i in chosen]
    state["prediction_note"] = prediction_note
    state["n_frames"] = len(frames)
    return prediction_note, properties, property_names, shown


@app.cell
def _(mo, property_names):
    numeric = [name for name in property_names if name != "dataset_group"]
    x_default = "energy" if "energy" in numeric else numeric[0]
    if "predicted_energy" in numeric:
        y_default = "predicted_energy"
    elif "error" in numeric:
        y_default = "error"
    else:
        y_default = numeric[min(1, len(numeric) - 1)]
    color_default = "error" if "error" in numeric else x_default
    x_axis = mo.ui.dropdown(numeric, value=x_default, label="Parity x")
    y_axis = mo.ui.dropdown(numeric, value=y_default, label="Parity y")
    color_by = mo.ui.dropdown(numeric, value=color_default, label="Color")
    palette = mo.ui.dropdown(
        [
            "bwr",
            "seismic",
            "inferno",
            "magma",
            "plasma",
            "viridis",
            "cividis",
            "twilight (periodic)",
            "tab10",
        ],
        value="bwr",
        label="Palette",
    )
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
def _():
    # evaluate the model on the selected structure, and print intermediates, and final energy, etc

    return


@app.cell
def _(mo, plt, read_metrics, run):
    metrics = read_metrics(run)
    if not metrics:
        training = mo.md(f"No `train.csv` in `{run}`.")
    else:
        epochs = [row["epoch"] for row in metrics]
        loss_fig, loss_axes = plt.subplots(1, 3, figsize=(11, 3.1))
        train_loss = [row["train_loss"] for row in metrics]
        val_loss = [row["val_loss"] for row in metrics]
        loss_axes[0].plot(epochs, train_loss, "o-", label="train")
        loss_axes[0].plot(epochs, val_loss, "o-", label="val")
        loss_axes[0].set_title("Loss")
        loss_axes[0].set_xlabel("epoch")
        loss_axes[0].legend()
        train_rmse = [row["train_rmse_meV"] for row in metrics]
        val_rmse = [row["val_rmse_meV"] for row in metrics]
        loss_axes[1].plot(epochs, train_rmse, "o-", label="train")
        loss_axes[1].plot(epochs, val_rmse, "o-", label="val")
        loss_axes[1].set_title("Energy RMSE (meV/atom)")
        loss_axes[1].set_xlabel("epoch")
        loss_axes[1].legend()
        loss_axes[2].plot(epochs, [row["learning_rate"] for row in metrics], "o-")
        loss_axes[2].set_title("Learning rate")
        loss_axes[2].set_xlabel("epoch")
        loss_fig.tight_layout()
        latest = metrics[-1]
        previous = metrics[-2] if len(metrics) > 1 else latest
        delta = latest["val_rmse_meV"] - previous["val_rmse_meV"]
        loss_direction = (
            "decrease" if latest["val_loss"] <= previous["val_loss"] else "increase"
        )
        cards = mo.hstack(
            [
                mo.stat(str(latest["epoch"]), label="epoch", bordered=True),
                mo.stat(
                    f"{latest['val_loss']:.3g}",
                    label="val loss",
                    caption=f"train {latest['train_loss']:.3g}",
                    direction=loss_direction,
                    target_direction="decrease",
                    bordered=True,
                ),
                mo.stat(
                    f"{latest['val_rmse_meV']:.0f}",
                    label="val RMSE meV/atom",
                    caption=f"{delta:+.0f} vs previous log",
                    direction="decrease" if delta <= 0 else "increase",
                    target_direction="decrease",
                    bordered=True,
                ),
                mo.stat(
                    f"{latest['learning_rate']:.3e}",
                    label="learning rate",
                    bordered=True,
                ),
            ],
            gap=1,
        )
        training = mo.vstack(
            [mo.md("## Training"), cards, mo.ui.table(metrics), loss_fig]
        )
    training
    return


@app.cell
def _(checkpoint, compare, mo, np, plt, read_checkpoint, state):
    if not checkpoint.value:
        weights = mo.md("No checkpoint in this run directory.")
    else:
        loaded, weight_rows, flats = read_checkpoint(checkpoint.value)
        weight_rows = [dict(row) for row in weight_rows]
        if compare.value and compare.value != checkpoint.value:
            _, other_rows, _ = read_checkpoint(compare.value)
            other = {row["name"]: row["l2"] for row in other_rows}
            for row in weight_rows:
                baseline = other.get(row["name"])
                row["l2_other"] = baseline
                row["l2_delta"] = None if baseline is None else row["l2"] - baseline
        values = np.concatenate(flats)
        low, high = np.percentile(values, [1, 99])
        clipped = values[(values >= low) & (values <= high)]
        weight_fig, weight_axes = plt.subplots(1, 2, figsize=(11, 3.2))
        order = np.argsort([row["l2"] for row in weight_rows])[::-1][:18]
        weight_axes[0].barh(
            [weight_rows[i]["name"] for i in order][::-1],
            [weight_rows[i]["l2"] for i in order][::-1],
        )
        weight_axes[0].set_xscale("log")
        weight_axes[0].set_title("Largest tensor L2 norms")
        weight_axes[1].hist(clipped, bins=80, color="#4c1d95")
        weight_axes[1].set_title("Parameter values (1st–99th percentile)")
        weight_fig.tight_layout()
        n_params = sum(row["numel"] for row in weight_rows)
        best_metric = loaded.get("best_metric")
        state["checkpoint"] = loaded
        state["weight_rows"] = weight_rows
        weights = mo.vstack(
            [
                mo.md("## Weights"),
                mo.hstack(
                    [
                        mo.stat(f"{n_params:,}", label="parameters", bordered=True),
                        mo.stat(
                            str(loaded.get("epoch")),
                            label="checkpoint epoch",
                            bordered=True,
                        ),
                        mo.stat(
                            (
                                f"{best_metric:.4g}"
                                if isinstance(best_metric, float)
                                else "—"
                            ),
                            label=f"best metric (epoch {loaded.get('best_epoch')})",
                            bordered=True,
                        ),
                        mo.stat(
                            str(loaded.get("architecture_name", "—")),
                            label="architecture",
                            bordered=True,
                        ),
                    ],
                    gap=1,
                ),
                weight_fig,
                mo.ui.table(weight_rows, page_size=12, selection=None),
            ]
        )
    weights
    return


@app.cell
def _(mo, re, run):
    log_path = run / "train.log"
    log_text = log_path.read_text() if log_path.is_file() else ""
    blocks = []
    pattern = (
        r"(Training|Validation|Test) dataset:\s+Dataset containing (\d+) structures"
        r".*?mean\s+([-\d.]+)\s+eV\s+- std\s+([\d.]+)\s+eV"
    )
    for match in re.finditer(pattern, log_text, flags=re.S):
        split, count, mean, std = match.groups()
        blocks.append(
            mo.stat(
                count,
                label=f"{split.lower()} structures",
                caption=f"energy {float(mean):.2f} ± {float(std):.1f} eV",
                bordered=True,
            )
        )
    options_path = run / "options_restart.yaml"
    options_text = options_path.read_text() if options_path.is_file() else ""
    summary = (
        mo.hstack(blocks, gap=1) if blocks else mo.md("No dataset summary in the log.")
    )
    mo.vstack(
        [
            mo.md("## Run log"),
            summary,
            mo.accordion(
                {
                    "train.log": mo.md(f"```\n{log_text[-4000:]}\n```"),
                    "options_restart.yaml": mo.md(f"```yaml\n{options_text}\n```"),
                }
            ),
        ]
    )
    return


@app.cell
def _(mo):
    scratch = mo.ui.code_editor(
        value=(
            "import numpy as np\n"
            "energies = state['properties']['energy']['values']\n"
            "print(len(state['frames']), 'structures on the map')\n"
            "print('energy mean', float(np.mean(energies)))\n"
            "print('checkpoint epoch', state['checkpoint'].get('epoch'))\n"
            "print(state['weight_rows'][0]['name'], state['weight_rows'][0]['l2'])\n"
        ),
        language="python",
        min_height=180,
        debounce=True,
        label="Python",
    )
    run_scratch = mo.ui.run_button(label="Run this Python")
    return run_scratch, scratch


@app.cell
def _(io, mo, run_scratch, scratch, state, traceback):
    output = mo.md(
        "Every cell in this notebook is Python. Edit them directly in "
        "`marimo edit`. The box below runs extra code against `state` "
        "(`frames`, `properties`, `indices`, `checkpoint`, `weight_rows`)."
    )
    if run_scratch.value:
        buffer = io.StringIO()
        namespace = {"state": state, "mo": mo}
        try:
            with mo.redirect_stdout(buffer):
                exec(compile(scratch.value, "<scratch>", "exec"), namespace)
            printed = buffer.getvalue() or "(no output)"
            output = mo.md(f"```\n{printed}\n```")
        except Exception:
            output = mo.md(f"```\n{traceback.format_exc()}\n```")
    mo.vstack([mo.md("## Python"), scratch, run_scratch, output])
    return


@app.cell
def _():
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
