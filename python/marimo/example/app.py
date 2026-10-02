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

    _lorem_models = {}

    def _load_lorem(path: str):
        cached = _lorem_models.get(path)
        if cached is not None:
            return cached
        loaded, _, _ = read_checkpoint(path)
        from metatrain.experimental.lorem.model import LOREM

        model = LOREM.load_checkpoint(loaded, context="restart")
        model.eval()
        model._lorem_epoch = loaded.get("epoch")
        _lorem_models[path] = model
        return model

    def infer(atoms, checkpoint_path=None):
        """Run the selected checkpoint on one ASE structure.

        Returns the final energy plus the short-range and long-range
        contributions, feature norms, and per-atom values. ``checkpoint_path``
        defaults to the checkpoint chosen in the dashboard.
        """
        import torch
        from metatomic.torch import System
        from metatrain.utils.neighbor_lists import get_system_with_neighbor_lists

        path = checkpoint_path or state.get("checkpoint_path")
        if not path:
            raise ValueError("Select a checkpoint before running inference.")
        model = _load_lorem(str(path))
        dtype = next(
            parameter.dtype
            for parameter in model.parameters()
            if parameter.is_floating_point()
        )
        system = System(
            positions=torch.tensor(atoms.positions, dtype=dtype),
            types=torch.tensor(atoms.numbers, dtype=torch.int32),
            cell=torch.tensor(atoms.cell.array, dtype=dtype),
            pbc=torch.tensor(list(atoms.pbc)),
        )
        system = get_system_with_neighbor_lists(
            system, model.requested_neighbor_lists()
        )
        with torch.no_grad():
            (
                nodes_scalar,
                distances,
                nodes_spherical,
                sr_energy,
                snapshots,
            ) = model.sr([system])
            lr_energy, spherical_updates = model.lr(
                [system], nodes_scalar, distances, nodes_spherical
            )
            predicted = model([system], model.outputs)
            charges = model.lr.charges(nodes_scalar, nodes_spherical)

        target = next(iter(predicted))
        atomic = predicted[target].block().values.detach().reshape(-1).cpu()
        sr = sr_energy.detach().reshape(-1).cpu()
        lr = lr_energy.detach().reshape(-1).cpu()
        charge = charges[:, 0].detach().cpu()
        scalar_norm = torch.linalg.vector_norm(nodes_scalar.detach(), dim=-1).cpu()
        spherical_norm = torch.linalg.vector_norm(
            nodes_spherical.detach().flatten(1), dim=-1
        ).cpu()
        final = float(atomic.sum())
        sr_total = float(sr.sum())
        lr_total = float(lr.sum())
        reference = atoms.info.get("ecumetric_energy")
        reference = None if reference is None else float(reference)
        intermediates = [
            {"name": "short-range energy (eV)", "value": sr_total},
            {"name": "long-range energy (eV)", "value": lr_total},
            {"name": "raw energy, sr + lr (eV)", "value": sr_total + lr_total},
            {"name": "final energy (eV)", "value": final},
        ]
        if reference is not None:
            intermediates.append(
                {"name": "reference energy (eV)", "value": reference}
            )
            intermediates.append(
                {"name": "final − reference (eV)", "value": final - reference}
            )
        intermediates.extend(
            [
                {"name": "atoms", "value": int(len(atoms))},
                {"name": "neighbor pairs", "value": int(distances.shape[0])},
                {
                    "name": "message-passing steps",
                    "value": int(snapshots.shape[0] - 1),
                },
                {
                    "name": "mean scalar-feature L2",
                    "value": float(scalar_norm.mean()),
                },
                {
                    "name": "mean spherical-feature L2",
                    "value": float(spherical_norm.mean()),
                },
                {"name": "scalar charge sum", "value": float(charge.sum())},
            ]
        )
        symbols = list(atoms.symbols)
        pair_samples = system.get_neighbor_list(model.requested_nl).samples.values
        pair_center = pair_samples[:, 0].detach().cpu().long()
        pair_neighbor = pair_samples[:, 1].detach().cpu().long()
        pair_distance = distances.detach().reshape(-1).cpu()
        if pair_distance.shape[0] != pair_center.shape[0]:
            pair_center = pair_center[:0]
        per_atom = []
        for i in range(len(atoms)):
            mask = pair_center == i
            row = {
                "atom": i,
                "symbol": symbols[i],
                "sr_energy": float(sr[i]),
                "lr_energy": float(lr[i]),
                "final_energy": float(atomic[i]),
                "scalar_charge": float(charge[i]),
                "scalar_feature_l2": float(scalar_norm[i]),
                "min_pair": None,
                "min_distance": None,
                "max_pair": None,
                "max_distance": None,
            }
            if int(mask.sum()) > 0:
                chosen = pair_distance[mask]
                partners = pair_neighbor[mask]
                near = int(torch.argmin(chosen))
                far = int(torch.argmax(chosen))
                near_atom = int(partners[near])
                far_atom = int(partners[far])
                row["min_pair"] = f"{symbols[near_atom]} {near_atom}"
                row["min_distance"] = float(chosen[near])
                row["max_pair"] = f"{symbols[far_atom]} {far_atom}"
                row["max_distance"] = float(chosen[far])
            per_atom.append(row)
        return {
            "formula": atoms.get_chemical_formula(),
            "n_atoms": len(atoms),
            "n_pairs": int(distances.shape[0]),
            "energy": final,
            "sr_energy": sr_total,
            "lr_energy": lr_total,
            "reference": reference,
            "target": target,
            "epoch": getattr(model, "_lorem_epoch", None),
            "intermediates": intermediates,
            "per_atom": per_atom,
            "hidden": {
                "nodes_scalar": nodes_scalar.detach().cpu().numpy(),
                "nodes_spherical": nodes_spherical.detach().cpu().numpy(),
                "charges": charges.detach().cpu().numpy(),
                "spherical_updates": spherical_updates.detach().cpu().numpy(),
                "snapshots": snapshots.detach().cpu().numpy(),
            },
        }

    state["infer"] = infer

    return (
        Path,
        chemiscope,
        comment_values,
        infer,
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
def _(checkpoint, infer, mo, np, plt, properties, state, traceback, viewer):
    _reference = np.asarray(properties["energy"]["values"], dtype=float)
    _counts = np.asarray(properties["n_atoms"]["values"], dtype=float)
    _predicted_spec = properties.get("predicted_energy")
    _pieces = [mo.md("## Inference")]
    if _predicted_spec is None:
        _pieces.append(
            mo.md(
                "No predictions file is aligned with the structures, so the "
                "map plots are empty. Click a point to run the checkpoint."
            )
        )
    else:
        _predicted = np.asarray(_predicted_spec["values"], dtype=float)
        _error = _predicted - _reference
        _per_atom_mev = _error / np.maximum(_counts, 1.0) * 1000.0
        _rmse = float(np.sqrt(np.mean(_error**2)))
        _mae = float(np.mean(np.abs(_error)))
        _rmse_atom = float(np.sqrt(np.mean(_per_atom_mev**2)))
        _bias = float(np.mean(_error))
        output_fig, output_axes = plt.subplots(
            2, 2, figsize=(11.2, 6.6), layout="constrained"
        )
        _span = float(
            max(
                abs(_reference.min()),
                abs(_reference.max()),
                abs(_predicted.min()),
                abs(_predicted.max()),
            )
        )
        _limit = max(_span, 1.0)
        _color_limit = float(np.percentile(np.abs(_per_atom_mev), 95))
        _color_limit = max(_color_limit, 1.0)
        _scatter = output_axes[0, 0].scatter(
            _reference,
            _predicted,
            c=_per_atom_mev,
            cmap="coolwarm",
            vmin=-_color_limit,
            vmax=_color_limit,
            s=16,
            alpha=0.85,
        )
        output_axes[0, 0].plot(
            [-_limit, _limit], [-_limit, _limit], color="0.35", lw=1
        )
        output_axes[0, 0].set_xlim(-_limit, _limit)
        output_axes[0, 0].set_ylim(-_limit, _limit)
        output_axes[0, 0].set_xlabel("reference energy (eV)")
        output_axes[0, 0].set_ylabel("predicted energy (eV)")
        output_axes[0, 0].set_title("Parity")
        output_fig.colorbar(
            _scatter, ax=output_axes[0, 0], label="error (meV/atom)", fraction=0.046
        )
        output_axes[0, 1].hist(_error, bins=40, color="#4c1d95", alpha=0.9)
        output_axes[0, 1].axvline(0.0, color="0.3", lw=1)
        output_axes[0, 1].axvline(_bias, color="#b45309", lw=1, ls="--")
        output_axes[0, 1].set_xlabel("predicted − reference (eV)")
        output_axes[0, 1].set_ylabel("structures")
        output_axes[0, 1].set_title("Energy error")
        output_axes[1, 0].scatter(_reference, _error, s=16, alpha=0.8, c="#1d4ed8")
        output_axes[1, 0].axhline(0.0, color="0.3", lw=1)
        output_axes[1, 0].set_xlabel("reference energy (eV)")
        output_axes[1, 0].set_ylabel("error (eV)")
        output_axes[1, 0].set_title("Error against reference")
        _groups = properties.get("dataset_group")
        if _groups:
            _labels = np.asarray(_groups["values"])
            _order = []
            for _name in _labels:
                if _name not in _order:
                    _order.append(_name)
            _order.sort(key=lambda name: int(np.sum(_labels == name)), reverse=True)
            _keep = _order[:8]
            output_axes[1, 1].boxplot(
                [_per_atom_mev[_labels == name] for name in _keep],
                showfliers=False,
            )
            output_axes[1, 1].set_xticks(
                range(1, len(_keep) + 1),
                [name.replace("_", " ")[:18] for name in _keep],
                rotation=30,
                ha="right",
            )
            output_axes[1, 1].axhline(0.0, color="0.3", lw=1)
            output_axes[1, 1].set_ylabel("error (meV/atom)")
            output_axes[1, 1].set_title("Error by dataset group")
        else:
            output_axes[1, 1].scatter(
                _counts, _per_atom_mev, s=16, alpha=0.8, c="#b45309"
            )
            output_axes[1, 1].axhline(0.0, color="0.3", lw=1)
            output_axes[1, 1].set_xlabel("atoms")
            output_axes[1, 1].set_ylabel("error (meV/atom)")
            output_axes[1, 1].set_title("Per-atom error")
        _pieces.extend(
            [
                mo.hstack(
                    [
                        mo.stat(
                            f"{_rmse:.3g}",
                            label="energy RMSE (eV)",
                            bordered=True,
                        ),
                        mo.stat(
                            f"{_mae:.3g}",
                            label="energy MAE (eV)",
                            bordered=True,
                        ),
                        mo.stat(
                            f"{_rmse_atom:.0f}",
                            label="RMSE (meV/atom)",
                            bordered=True,
                        ),
                        mo.stat(
                            f"{_bias:+.3g}",
                            label="mean error (eV)",
                            bordered=True,
                        ),
                    ],
                    gap=1,
                ),
                mo.md(
                    f"{len(_reference)} structures from the predictions file. "
                    "The dashed line on the histogram is the mean error."
                ),
                output_fig,
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
            _atoms = np.arange(_result["n_atoms"])
            _sr = np.array([row["sr_energy"] for row in _result["per_atom"]])
            _lr = np.array([row["lr_energy"] for row in _result["per_atom"]])
            _final = np.array([row["final_energy"] for row in _result["per_atom"]])
            _charge = np.array([row["scalar_charge"] for row in _result["per_atom"]])
            _features = np.array(
                [row["scalar_feature_l2"] for row in _result["per_atom"]]
            )
            atom_fig, atom_axes = plt.subplots(
                2, 2, figsize=(11.2, 6.4), layout="constrained"
            )
            atom_axes[0, 0].plot(_atoms, _sr, "o-", ms=3, label="short-range")
            atom_axes[0, 0].plot(_atoms, _lr, "o-", ms=3, label="long-range")
            atom_axes[0, 0].plot(_atoms, _final, "o-", ms=3, label="final")
            atom_axes[0, 0].set_xlabel("atom")
            atom_axes[0, 0].set_ylabel("energy (eV)")
            atom_axes[0, 0].set_title("Per-atom energy")
            atom_axes[0, 0].legend()
            atom_axes[0, 1].axhline(0.0, color="0.75", lw=1)
            atom_axes[0, 1].bar(
                ["short-range", "long-range", "raw", "final"],
                [
                    _result["sr_energy"],
                    _result["lr_energy"],
                    _result["sr_energy"] + _result["lr_energy"],
                    _result["energy"],
                ],
                color=["#1d4ed8", "#b45309", "#6d28d9", "#111827"],
            )
            if _result["reference"] is not None:
                atom_axes[0, 1].axhline(
                    _result["reference"], color="#dc2626", lw=1, ls="--"
                )
            atom_axes[0, 1].set_ylabel("energy (eV)")
            atom_axes[0, 1].set_title("Structure total")
            atom_axes[1, 0].plot(_atoms, _charge, "o-", ms=3, color="#0f766e")
            atom_axes[1, 0].axhline(0.0, color="0.75", lw=1)
            atom_axes[1, 0].set_xlabel("atom")
            atom_axes[1, 0].set_ylabel("scalar charge")
            atom_axes[1, 0].set_title("Long-range charges")
            atom_axes[1, 1].plot(_atoms, _features, "o-", ms=3, color="#4c1d95")
            atom_axes[1, 1].set_xlabel("atom")
            atom_axes[1, 1].set_ylabel("L2")
            atom_axes[1, 1].set_title("Scalar feature norm")
            _reference_energy = _result["reference"]
            _signed = (
                None
                if _reference_energy is None
                else _result["energy"] - _reference_energy
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
                    atom_fig,
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
    infer,
    matrix_aspect,
    matrix_controls,
    matrix_transform,
    mo,
    np,
    plt,
    read_checkpoint,
    state,
    traceback,
    viewer,
):
    def _as_tiles(values):
        array = np.asarray(values, dtype=float)
        if array.ndim <= 1:
            return array.reshape(1, -1)
        if array.ndim == 2:
            return array
        blocks = array.reshape(-1, *array.shape[-2:])
        count, height, width = blocks.shape
        columns = int(np.ceil(np.sqrt(count)))
        rows = int(np.ceil(count / columns))
        canvas = np.full(
            (rows * (height + 1) - 1, columns * (width + 1) - 1),
            np.nan,
        )
        for index, block in enumerate(blocks):
            row, column = divmod(index, columns)
            top = row * (height + 1)
            left = column * (width + 1)
            canvas[top : top + height, left : left + width] = block
        return canvas

    def _show_matrix(ax, values, title):
        image = _as_tiles(values)
        kind = matrix_transform.value
        if kind == "abs":
            image = np.abs(image)
            diverging = False
        elif kind == "log10":
            image = np.log10(1.0 + np.abs(image))
            diverging = False
        elif kind == "signed_log":
            image = np.sign(image) * np.log10(1.0 + np.abs(image))
            diverging = True
        elif kind == "asinh":
            image = np.arcsinh(image)
            diverging = True
        elif kind == "sqrt":
            image = np.sign(image) * np.sqrt(np.abs(image))
            diverging = True
        else:
            diverging = True
        finite = image[np.isfinite(image)]
        if finite.size == 0:
            low, high = -1.0, 1.0
        elif diverging:
            limit = max(float(np.percentile(np.abs(finite), 99)), 1e-8)
            low, high = -limit, limit
        else:
            low = float(np.percentile(finite, 1))
            high = float(np.percentile(finite, 99))
            if high <= low:
                high = low + 1e-8
        cmap = plt.get_cmap("coolwarm" if diverging else "magma").copy()
        cmap.set_bad("#f4f4f5")
        ax.imshow(
            np.ma.masked_invalid(image),
            cmap=cmap,
            vmin=low,
            vmax=high,
            aspect=matrix_aspect.value,
            interpolation="nearest",
        )
        ax.set_title(title, fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])

    if not checkpoint.value:
        shapes_view = mo.md("## Shapes\n\nSelect a checkpoint.")
    else:
        _loaded, _, _ = read_checkpoint(checkpoint.value)
        _hypers = _loaded["model_data"]["model_hypers"]
        _features = int(_hypers["num_features"])
        _spherical = int(_hypers["num_spherical_features"])
        _species = int(_hypers["num_species"])
        _radial = int(_hypers["num_radial"])
        _degree = int(_hypers["max_degree"])
        _degree_lr = int(_hypers["max_degree_lr"])
        _components = (_degree + 1) ** 2
        _components_lr = (_degree_lr + 1) ** 2
        _stages = int(_hypers["num_message_passing"]) + 1
        _weights = _loaded["model_state_dict"]
        _panels = []
        _radial_out = _weights.get("sr.radial_coefficients.2.weight")
        if _radial_out is not None and tuple(_radial_out.shape) == (
            _radial * _features,
            _features,
        ):
            _radial_out = _radial_out.reshape(_radial, _features, _features)
        for _title, _key, _values in (
            (
                f"embedding ({_species} species)",
                "sr.chemical_embedding.weight",
                _weights.get("sr.chemical_embedding.weight"),
            ),
            (
                "scalar map, species to features",
                "sr.dense0.0.weight",
                _weights.get("sr.dense0.0.weight"),
            ),
            (
                "scalar message, features x features",
                "sr.dense1.weight",
                _weights.get("sr.dense1.weight"),
            ),
            (
                "features to spherical channels",
                "sr.dense2.weight",
                _weights.get("sr.dense2.weight"),
            ),
            (
                "radial MLP, in",
                "sr.radial_coefficients.0.weight",
                _weights.get("sr.radial_coefficients.0.weight"),
            ),
            (
                f"radial MLP, {_radial} feature blocks",
                "sr.radial_coefficients.2.weight",
                _radial_out,
            ),
            (
                "spherical tensor product",
                "sr.tensor_dense.tensor_weight",
                _weights.get("sr.tensor_dense.tensor_weight"),
            ),
            (
                "short-range energy readout",
                "sr.energy_mlp.4.weight",
                _weights.get("sr.energy_mlp.4.weight"),
            ),
            (
                "scalar charge MLP",
                "lr.scalar_charge_mlp.0.weight",
                _weights.get("lr.scalar_charge_mlp.0.weight"),
            ),
            (
                "spherical charges",
                "lr.spherical_charge_dense.dense.weight",
                _weights.get("lr.spherical_charge_dense.dense.weight"),
            ),
            (
                "potential to features",
                "lr.potential_to_features.weight",
                _weights.get("lr.potential_to_features.weight"),
            ),
            (
                "long-range energy readout",
                "lr.energy_mlp.4.weight",
                _weights.get("lr.energy_mlp.4.weight"),
            ),
        ):
            if _values is not None:
                _panels.append((_title, _key, _values.detach().cpu().numpy()))
        _columns = 4
        _rows = int(np.ceil(len(_panels) / _columns))
        shape_fig, shape_axes = plt.subplots(
            _rows,
            _columns,
            figsize=(11.4, 2.15 * _rows),
            layout="constrained",
        )
        _flat_axes = np.atleast_1d(shape_axes).ravel()
        for _ax, _panel in zip(_flat_axes, _panels, strict=False):
            _show_matrix(_ax, _panel[2], _panel[0])
        for _ax in _flat_axes[len(_panels) :]:
            _ax.axis("off")
        _optimizer = _loaded.get("optimizer_state_dict") or {}
        _adam_state = _optimizer.get("state") or {}
        _buffer_marks = (
            "buffer",
            "bernstein_coeff",
            "type_to_index",
            "smearing",
            "prefactor",
        )
        _param_names = [
            name
            for name in _weights
            if not any(mark in name for mark in _buffer_marks)
        ]
        _adam_by_name = {}
        if len(_param_names) == len(_adam_state):
            for _index, _name in enumerate(_param_names):
                _adam_by_name[_name] = _adam_state[_index]
        _group = (_optimizer.get("param_groups") or [{}])[0]
        _betas = _group.get("betas", (0.9, 0.999))
        _adam_lr = _group.get("lr")
        _adam_step = None
        if _adam_state:
            _adam_step = int(next(iter(_adam_state.values()))["step"])
        adam_fig = None
        if _adam_by_name:
            adam_fig, adam_axes = plt.subplots(
                len(_panels),
                3,
                figsize=(11.2, 1.2 * len(_panels)),
                layout="constrained",
            )
            adam_axes = np.atleast_2d(adam_axes)
            for _row, (_title, _key, _view) in enumerate(_panels):
                _moment = _adam_by_name.get(_key)
                _show_matrix(
                    adam_axes[_row, 0],
                    _view,
                    "parameter" if _row == 0 else "",
                )
                adam_axes[_row, 0].set_ylabel(_title, fontsize=7)
                if _moment is None:
                    adam_axes[_row, 1].axis("off")
                    adam_axes[_row, 2].axis("off")
                    continue
                _first = _moment["exp_avg"].detach().cpu().numpy()
                _second = _moment["exp_avg_sq"].detach().cpu().numpy()
                if _first.size == _view.size:
                    _first = _first.reshape(_view.shape)
                    _second = _second.reshape(_view.shape)
                _show_matrix(
                    adam_axes[_row, 1],
                    _first,
                    "first moment m" if _row == 0 else "",
                )
                _show_matrix(
                    adam_axes[_row, 2],
                    np.sqrt(np.maximum(_second, 0.0)),
                    "root second moment" if _row == 0 else "",
                )
            _lr_text = f"{_adam_lr:.4g}" if isinstance(_adam_lr, float) else "—"
            _adam_note = (
                f"Adam keeps `exp_avg` and `exp_avg_sq` for each of "
                f"{len(_adam_by_name)} parameters, at step {_adam_step}. "
                f"betas are {_betas[0]} and {_betas[1]}, learning rate {_lr_text}. "
                "Each row is the parameter, its first moment, and the root of "
                "its second moment, tiled the same way. "
                "`best_optimizer_state_dict` stores the moments from the best "
                "epoch; these are the moments saved with the weights above."
            )
        else:
            _adam_note = "This checkpoint has no Adam state."
        _ledger = [
            {
                "tensor": "nodes_scalar",
                "shape": f"(n_atoms, {_features})",
                "role": "invariant hidden state",
            },
            {
                "tensor": "nodes_spherical",
                "shape": f"(n_atoms, {_components}, {_spherical})",
                "role": f"(max_degree+1)^2 = {_components} equivariant channels",
            },
            {
                "tensor": "snapshots",
                "shape": f"({_stages}, n_atoms, {_components}, {_spherical})",
                "role": "spherical state after each message-passing stage",
            },
            {
                "tensor": "radial basis",
                "shape": f"(n_pairs, {_radial})",
                "role": f"Bernstein basis inside the {_hypers['cutoff']} A cutoff",
            },
            {
                "tensor": "charges",
                "shape": f"(n_atoms, {1 + _components_lr})",
                "role": "scalar charge plus one long-range (l, m) channel",
            },
            {
                "tensor": "atomic energy",
                "shape": "(n_atoms,)",
                "role": "short-range plus long-range, before the structure sum",
            },
            {
                "tensor": "structure energy",
                "shape": "(n_structures,)",
                "role": "sum of atomic energies; this is the training target",
            },
            {
                "tensor": "forces",
                "shape": "(n_atoms, 3)",
                "role": "gradient of the structure energy w.r.t. positions",
            },
            {
                "tensor": "stress",
                "shape": "(n_structures, 3, 3)",
                "role": "gradient of the structure energy w.r.t. strain",
            },
        ]
        _map_index = (viewer.selected_ids or {}).get("structure")
        _hidden_note = (
            "Click a point to draw the hidden state of that structure. "
            "A training batch uses the same tensors with `n_atoms` equal to "
            "every atom in the batch concatenated, not a padded "
            "`(batch, max_atoms, features)` array."
        )
        _hidden_fig = None
        if _map_index is not None and "frames" in state:
            try:
                _result = infer(state["frames"][_map_index], checkpoint.value)
            except Exception:
                _hidden_note = f"```\n{traceback.format_exc()}\n```"
            else:
                state["inference"] = _result
                _hidden = _result["hidden"]
                _scalar = _hidden["nodes_scalar"]
                _equivariant = _hidden["nodes_spherical"]
                _charges = _hidden["charges"]
                _updates = _hidden["spherical_updates"]
                for _row, _array in (
                    (_ledger[0], _scalar),
                    (_ledger[1], _equivariant),
                    (_ledger[2], _hidden["snapshots"]),
                    (_ledger[4], _charges),
                ):
                    _row["this structure"] = "×".join(
                        str(size) for size in _array.shape
                    )
                _ledger[3]["this structure"] = str(_result["n_pairs"])
                _ledger[5]["this structure"] = str(_result["n_atoms"])
                hidden_fig, hidden_axes = plt.subplots(
                    2, 2, figsize=(11.2, 6.2), layout="constrained"
                )
                _show_matrix(
                    hidden_axes[0, 0],
                    _scalar,
                    f"scalar features {_scalar.shape}",
                )
                _show_matrix(
                    hidden_axes[0, 1],
                    np.transpose(_equivariant, (1, 0, 2)),
                    f"spherical features, {_components} (l, m) blocks",
                )
                _show_matrix(
                    hidden_axes[1, 0],
                    _charges,
                    f"charges {_charges.shape}",
                )
                _show_matrix(
                    hidden_axes[1, 1],
                    np.transpose(_updates, (1, 0, 2)),
                    f"long-range spherical update {_updates.shape}",
                )
                _hidden_fig = hidden_fig
                _hidden_note = (
                    f"**{_result['formula']}**, {_result['n_atoms']} atoms, "
                    f"{_result['n_pairs']} neighbor pairs. "
                    "Spherical panels are stacked by (l, m). "
                    "Training sums the atomic energies and backpropagates "
                    "forces and stress; eval mode also adds the composition "
                    "baseline and the output scale inside `forward`."
                )
        shapes_view = mo.vstack(
            [
                mo.md(
                    "## Shapes\n\n"
                    "Each atom carries one scalar vector and one spherical "
                    f"tensor. Here that is `{_features}` features and "
                    f"`{_components}` × `{_spherical}` equivariant channels "
                    f"(`max_degree={_degree}`). "
                    f"`num_message_passing` is "
                    f"{_hypers['num_message_passing']}, so the spherical "
                    f"state has `{_stages}` stage. "
                    "Color scales are independent and clipped to the 99th "
                    "percentile after the transform. Aspect and transform "
                    "apply to the parameter, Adam, and hidden-state matrices."
                ),
                matrix_controls,
                mo.ui.table(_ledger, selection=None),
                mo.md("### Parameters"),
                shape_fig,
                mo.md("### Adam\n\n" + _adam_note),
                *([adam_fig] if adam_fig is not None else []),
                mo.md("### Hidden state\n\n" + _hidden_note),
                *([_hidden_fig] if _hidden_fig is not None else []),
            ]
        )
    shapes_view
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
        state["checkpoint_path"] = checkpoint.value
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
def _(infer, io, mo, run_scratch, scratch, state, traceback):
    import contextlib

    output = mo.md(
        "Every cell in this notebook is Python. `infer(atoms)` evaluates the "
        "selected checkpoint. `state` holds `frames`, `properties`, `indices`, "
        "`checkpoint`, `checkpoint_path`, `weight_rows`, and `inference`."
    )
    if run_scratch.value:
        buffer = io.StringIO()
        namespace = {"state": state, "mo": mo, "infer": infer}
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
