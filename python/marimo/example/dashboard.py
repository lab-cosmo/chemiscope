"""Reusable pieces for the metatrain marimo dashboard.

``app.py`` only builds the widgets. Loading a run, evaluating one structure,
and drawing the figures happen here, so the notebook cells stay short and
these functions can be called from the scratch cell or another script.
"""

import functools
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


PALETTES = [
    "bwr",
    "seismic",
    "inferno",
    "magma",
    "plasma",
    "viridis",
    "cividis",
    "twilight (periodic)",
    "tab10",
]

# Names that sit in model_state_dict but are buffers, not Adam parameters.
_BUFFER_MARKS = (
    "buffer",
    "bernstein_coeff",
    "type_to_index",
    "smearing",
    "prefactor",
)


@functools.lru_cache(maxsize=2)
def read_frames(path: str):
    import ase.io

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


def dataset_summaries(log_text: str) -> list[dict]:
    pattern = (
        r"(Training|Validation|Test) dataset:\s+Dataset containing (\d+) structures"
        r".*?mean\s+([-\d.]+)\s+eV\s+- std\s+([\d.]+)\s+eV"
    )
    rows = []
    for match in re.finditer(pattern, log_text, flags=re.S):
        split, count, mean, std = match.groups()
        rows.append(
            {
                "split": split,
                "count": count,
                "mean": float(mean),
                "std": float(std),
            }
        )
    return rows


def split_paths(value) -> list[str]:
    """One path per line, or a sequence of paths. Blank lines are ignored."""
    if value is None:
        return []
    if isinstance(value, (str, Path)):
        return [line.strip() for line in str(value).splitlines() if line.strip()]
    return [str(item).strip() for item in value if str(item).strip()]


def _unique_label(path: Path, paths: list[Path], *, qualify: bool = False) -> str:
    name = path.name
    if not qualify and sum(item.name == name for item in paths) == 1:
        return name
    short = f"{path.parent.name}/{name}"
    if sum(f"{item.parent.name}/{item.name}" == short for item in paths) == 1:
        return short
    return f"{path.parent.parent.name}/{short}"


def _checkpoint_epoch(path: Path) -> int:
    suffix = path.stem.split("_")[-1]
    return int(suffix) if suffix.isdigit() else 0


def discover_checkpoints(run_dirs, extra_paths=None) -> dict:
    """Checkpoints from every run directory, plus any extra ``.ckpt`` paths.

    ``model_*.ckpt`` in each directory is included. Labels stay unique when
    several runs contain ``model_60.ckpt``. ``runs`` is the directories that
    exist, in the order they were given.
    """
    directories = [Path(item).expanduser() for item in split_paths(run_dirs)]
    extras = [Path(item).expanduser() for item in split_paths(extra_paths)]
    missing = []
    runs = []
    found: list[Path] = []
    for directory in directories:
        if not directory.is_dir():
            missing.append(str(directory))
            continue
        runs.append(directory)
        found.extend(sorted(directory.glob("model_*.ckpt"), key=_checkpoint_epoch))
    for path in extras:
        if path.is_file():
            found.append(path)
        else:
            missing.append(str(path))
    def _identity(path: Path) -> str:
        return str(path.resolve()) if path.exists() else str(path)

    unique: list[Path] = []
    seen = set()
    for path in found:
        key = _identity(path)
        if key not in seen:
            seen.add(key)
            unique.append(path)
    rank = {_identity(path): index for index, path in enumerate(unique)}
    unique.sort(key=lambda path: (str(path.parent), _checkpoint_epoch(path), path.name))
    qualify = len({path.parent for path in unique}) > 1
    labels = [_unique_label(path, unique, qualify=qualify) for path in unique]
    options = {label: str(path) for label, path in zip(labels, unique, strict=True)}
    if not options:
        return {
            "options": {"(none)": ""},
            "runs": runs,
            "missing": missing,
            "default": "(none)",
            "compare": None,
        }
    best = max(
        unique,
        key=lambda path: (
            _checkpoint_epoch(path),
            rank[_identity(path)],
        ),
    )
    earlier = [
        path
        for path in unique
        if path.parent == best.parent
        and _checkpoint_epoch(path) < _checkpoint_epoch(best)
    ]
    compare = (
        _unique_label(max(earlier, key=_checkpoint_epoch), unique, qualify=qualify)
        if earlier
        else None
    )
    return {
        "options": options,
        "runs": runs,
        "missing": missing,
        "default": _unique_label(best, unique, qualify=qualify),
        "compare": compare,
    }


def build_map(structures, predictions, n_show: int) -> dict:
    """Sample structures from one or more files for the parity map.

    ``structures`` and ``predictions`` are a path, a newline-separated list,
    or a sequence of paths. Prediction files pair with structure files in
    order. One structure file keeps the previous single-file behavior.
    """
    structure_paths = [Path(item).expanduser() for item in split_paths(structures)]
    prediction_paths = [Path(item).expanduser() for item in split_paths(predictions)]
    if not structure_paths:
        raise ValueError("Add at least one structure file, one path per line.")
    notes = []
    paired: list[Path | None] = [None] * len(structure_paths)
    if prediction_paths and len(prediction_paths) != len(structure_paths):
        notes.append(
            f"{len(prediction_paths)} prediction files for "
            f"{len(structure_paths)} structure files. "
            "List them in the same order, one path per line."
        )
    elif prediction_paths:
        paired = list(prediction_paths)

    labels = [_unique_label(path, structure_paths) for path in structure_paths]
    records = []
    for label, path, pred_path in zip(labels, structure_paths, paired, strict=True):
        if not path.is_file():
            notes.append(f"Missing structure file: {path}")
            continue
        frames = read_frames(str(path))
        reference = comment_values(str(path), "ecumetric_energy")
        groups = comment_values(str(path), "dataset_group")
        predicted: list[float | None] = [None] * len(frames)
        if pred_path is not None and not pred_path.is_file():
            notes.append(f"Missing predictions file: {pred_path}")
        elif pred_path is not None:
            text = comment_values(str(pred_path), "energy")
            if len(text) == len(frames):
                predicted = [float(value) for value in text]
                if len(structure_paths) == 1:
                    notes.append(
                        f"Predictions aligned with {len(predicted)} structures."
                    )
                else:
                    notes.append(f"{label}: predictions aligned ({len(frames)}).")
            else:
                notes.append(
                    f"{label}: predictions have {len(text)} frames, "
                    f"structures have {len(frames)}. Error omitted for this file."
                )
        if len(reference) != len(frames):
            notes.append(
                f"{label}: ecumetric_energy has {len(reference)} values "
                f"for {len(frames)} frames."
            )
        for index, frame in enumerate(frames):
            energy = float(reference[index]) if index < len(reference) else float("nan")
            records.append(
                {
                    "frame": frame,
                    "energy": energy,
                    "group": groups[index] if index < len(groups) else "",
                    "predicted": predicted[index],
                    "source": label,
                    "index": index,
                }
            )
    if not records:
        detail = "\n".join(notes) or "No structures could be read."
        raise ValueError(detail)

    take = min(int(n_show), len(records))
    chosen = [
        records[int(i)] for i in np.linspace(0, len(records) - 1, take, dtype=int)
    ]
    energies = [row["energy"] for row in chosen]
    shown = [row["frame"] for row in chosen]
    per_atom = [
        energy / max(len(frame), 1)
        for energy, frame in zip(energies, shown, strict=True)
    ]
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
            "values": [len(frame) for frame in shown],
            "description": "Number of atoms",
        },
    }
    sources = [row["source"] for row in chosen]
    if len(set(sources)) > 1:
        properties["input_file"] = {
            "target": "structure",
            "values": sources,
            "description": "Structure file this frame was read from",
        }
    if any(row["group"] for row in chosen):
        properties["dataset_group"] = {
            "target": "structure",
            "values": [row["group"] for row in chosen],
            "description": "MAD dataset_group",
        }
    if chosen and all(row["predicted"] is not None for row in chosen):
        pred = [float(row["predicted"]) for row in chosen]
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
    if not notes:
        notes.append("No predictions file.")
    return {
        "shown": shown,
        "properties": properties,
        "property_names": list(properties),
        "prediction_note": " ".join(notes),
        "indices": [row["index"] for row in chosen],
        "sources": sources,
        "n_frames": len(records),
    }


_STRING_PROPERTIES = {"dataset_group", "input_file"}


def axis_defaults(property_names: list[str]):
    numeric = [name for name in property_names if name not in _STRING_PROPERTIES]
    x_default = "energy" if "energy" in numeric else numeric[0]
    if "predicted_energy" in numeric:
        y_default = "predicted_energy"
    elif "error" in numeric:
        y_default = "error"
    else:
        y_default = numeric[min(1, len(numeric) - 1)]
    color_default = "error" if "error" in numeric else x_default
    return numeric, x_default, y_default, color_default


_lorem_models = {}


def load_lorem(path: str):
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


def infer(atoms, checkpoint_path: str):
    """Run one checkpoint on one ASE structure.

    Returns the final energy, short-range and long-range contributions,
    per-atom values (including the closest and farthest neighbor), and the
    hidden tensors as numpy arrays.
    """
    import torch
    from metatomic.torch import System
    from metatrain.utils.neighbor_lists import get_system_with_neighbor_lists

    model = load_lorem(str(checkpoint_path))
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
    system = get_system_with_neighbor_lists(system, model.requested_neighbor_lists())
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
        intermediates.append({"name": "reference energy (eV)", "value": reference})
        intermediates.append(
            {"name": "final − reference (eV)", "value": final - reference}
        )
    intermediates.extend(
        [
            {"name": "atoms", "value": int(len(atoms))},
            {"name": "neighbor pairs", "value": int(distances.shape[0])},
            {"name": "message-passing steps", "value": int(snapshots.shape[0] - 1)},
            {"name": "mean scalar-feature L2", "value": float(scalar_norm.mean())},
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


def prediction_figure(properties: dict):
    """Parity, error histogram, and group boxplot for a predictions file.

    Returns ``(figure, summary)``, or ``(None, None)`` when that file is absent.
    """
    predicted_spec = properties.get("predicted_energy")
    if predicted_spec is None:
        return None, None
    reference = np.asarray(properties["energy"]["values"], dtype=float)
    counts = np.asarray(properties["n_atoms"]["values"], dtype=float)
    predicted = np.asarray(predicted_spec["values"], dtype=float)
    error = predicted - reference
    per_atom_mev = error / np.maximum(counts, 1.0) * 1000.0
    summary = {
        "n": int(reference.shape[0]),
        "rmse": float(np.sqrt(np.mean(error**2))),
        "mae": float(np.mean(np.abs(error))),
        "rmse_atom": float(np.sqrt(np.mean(per_atom_mev**2))),
        "bias": float(np.mean(error)),
    }
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 6.6), layout="constrained")
    span = float(
        max(
            abs(reference.min()),
            abs(reference.max()),
            abs(predicted.min()),
            abs(predicted.max()),
        )
    )
    limit = max(span, 1.0)
    color_limit = max(float(np.percentile(np.abs(per_atom_mev), 95)), 1.0)
    scatter = axes[0, 0].scatter(
        reference,
        predicted,
        c=per_atom_mev,
        cmap="coolwarm",
        vmin=-color_limit,
        vmax=color_limit,
        s=16,
        alpha=0.85,
    )
    axes[0, 0].plot([-limit, limit], [-limit, limit], color="0.35", lw=1)
    axes[0, 0].set_xlim(-limit, limit)
    axes[0, 0].set_ylim(-limit, limit)
    axes[0, 0].set_xlabel("reference energy (eV)")
    axes[0, 0].set_ylabel("predicted energy (eV)")
    axes[0, 0].set_title("Parity")
    fig.colorbar(scatter, ax=axes[0, 0], label="error (meV/atom)", fraction=0.046)
    axes[0, 1].hist(error, bins=40, color="#4c1d95", alpha=0.9)
    axes[0, 1].axvline(0.0, color="0.3", lw=1)
    axes[0, 1].axvline(summary["bias"], color="#b45309", lw=1, ls="--")
    axes[0, 1].set_xlabel("predicted − reference (eV)")
    axes[0, 1].set_ylabel("structures")
    axes[0, 1].set_title("Energy error")
    axes[1, 0].scatter(reference, error, s=16, alpha=0.8, c="#1d4ed8")
    axes[1, 0].axhline(0.0, color="0.3", lw=1)
    axes[1, 0].set_xlabel("reference energy (eV)")
    axes[1, 0].set_ylabel("error (eV)")
    axes[1, 0].set_title("Error against reference")
    groups = properties.get("dataset_group")
    if groups:
        labels = np.asarray(groups["values"])
        order = []
        for name in labels:
            if name not in order:
                order.append(name)
        order.sort(key=lambda name: int(np.sum(labels == name)), reverse=True)
        keep = order[:8]
        axes[1, 1].boxplot(
            [per_atom_mev[labels == name] for name in keep],
            showfliers=False,
        )
        axes[1, 1].set_xticks(
            range(1, len(keep) + 1),
            [name.replace("_", " ")[:18] for name in keep],
            rotation=30,
            ha="right",
        )
        axes[1, 1].axhline(0.0, color="0.3", lw=1)
        axes[1, 1].set_ylabel("error (meV/atom)")
        axes[1, 1].set_title("Error by dataset group")
    else:
        axes[1, 1].scatter(counts, per_atom_mev, s=16, alpha=0.8, c="#b45309")
        axes[1, 1].axhline(0.0, color="0.3", lw=1)
        axes[1, 1].set_xlabel("atoms")
        axes[1, 1].set_ylabel("error (meV/atom)")
        axes[1, 1].set_title("Per-atom error")
    return fig, summary


def atom_figure(result: dict):
    """Per-atom energies, the structure total, charges, and feature norms."""
    atoms = np.arange(result["n_atoms"])
    rows = result["per_atom"]
    sr = np.array([row["sr_energy"] for row in rows])
    lr = np.array([row["lr_energy"] for row in rows])
    final = np.array([row["final_energy"] for row in rows])
    charge = np.array([row["scalar_charge"] for row in rows])
    features = np.array([row["scalar_feature_l2"] for row in rows])
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 6.4), layout="constrained")
    axes[0, 0].plot(atoms, sr, "o-", ms=3, label="short-range")
    axes[0, 0].plot(atoms, lr, "o-", ms=3, label="long-range")
    axes[0, 0].plot(atoms, final, "o-", ms=3, label="final")
    axes[0, 0].set_xlabel("atom")
    axes[0, 0].set_ylabel("energy (eV)")
    axes[0, 0].set_title("Per-atom energy")
    axes[0, 0].legend()
    axes[0, 1].axhline(0.0, color="0.75", lw=1)
    axes[0, 1].bar(
        ["short-range", "long-range", "raw", "final"],
        [
            result["sr_energy"],
            result["lr_energy"],
            result["sr_energy"] + result["lr_energy"],
            result["energy"],
        ],
        color=["#1d4ed8", "#b45309", "#6d28d9", "#111827"],
    )
    if result["reference"] is not None:
        axes[0, 1].axhline(result["reference"], color="#dc2626", lw=1, ls="--")
    axes[0, 1].set_ylabel("energy (eV)")
    axes[0, 1].set_title("Structure total")
    axes[1, 0].plot(atoms, charge, "o-", ms=3, color="#0f766e")
    axes[1, 0].axhline(0.0, color="0.75", lw=1)
    axes[1, 0].set_xlabel("atom")
    axes[1, 0].set_ylabel("scalar charge")
    axes[1, 0].set_title("Long-range charges")
    axes[1, 1].plot(atoms, features, "o-", ms=3, color="#4c1d95")
    axes[1, 1].set_xlabel("atom")
    axes[1, 1].set_ylabel("L2")
    axes[1, 1].set_title("Scalar feature norm")
    return fig


def _as_tiles(values):
    """Lay every leading axis out as a mosaic, with a one-pixel gap."""
    array = np.asarray(values, dtype=float)
    if array.ndim <= 1:
        return array.reshape(1, -1)
    if array.ndim == 2:
        return array
    blocks = array.reshape(-1, *array.shape[-2:])
    count, height, width = blocks.shape
    columns = int(np.ceil(np.sqrt(count)))
    rows = int(np.ceil(count / columns))
    canvas = np.full((rows * (height + 1) - 1, columns * (width + 1) - 1), np.nan)
    for index, block in enumerate(blocks):
        row, column = divmod(index, columns)
        top = row * (height + 1)
        left = column * (width + 1)
        canvas[top : top + height, left : left + width] = block
    return canvas


def _show_matrix(ax, values, title, transform, aspect):
    image = _as_tiles(values)
    if transform == "abs":
        image = np.abs(image)
        diverging = False
    elif transform == "log10":
        image = np.log10(1.0 + np.abs(image))
        diverging = False
    elif transform == "signed_log":
        image = np.sign(image) * np.log10(1.0 + np.abs(image))
        diverging = True
    elif transform == "asinh":
        image = np.arcsinh(image)
        diverging = True
    elif transform == "sqrt":
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
        aspect=aspect,
        interpolation="nearest",
    )
    ax.set_title(title, fontsize=8)
    ax.set_xticks([])
    ax.set_yticks([])


def _layout(hypers: dict) -> dict:
    degree = int(hypers["max_degree"])
    degree_lr = int(hypers["max_degree_lr"])
    return {
        "features": int(hypers["num_features"]),
        "spherical": int(hypers["num_spherical_features"]),
        "species": int(hypers["num_species"]),
        "radial": int(hypers["num_radial"]),
        "degree": degree,
        "components": (degree + 1) ** 2,
        "components_lr": (degree_lr + 1) ** 2,
        "stages": int(hypers["num_message_passing"]) + 1,
        "cutoff": hypers["cutoff"],
        "message_passing": hypers["num_message_passing"],
    }


def shapes_intro(hypers: dict) -> str:
    layout = _layout(hypers)
    return (
        "## Shapes\n\n"
        "Each atom carries one scalar vector and one spherical "
        f"tensor. Here that is `{layout['features']}` features and "
        f"`{layout['components']}` × `{layout['spherical']}` equivariant channels "
        f"(`max_degree={layout['degree']}`). "
        f"`num_message_passing` is {layout['message_passing']}, so the spherical "
        f"state has `{layout['stages']}` stage. "
        "Color scales are independent and clipped to the 99th "
        "percentile after the transform. Aspect and transform "
        "apply to the parameter, Adam, and hidden-state matrices."
    )


def parameter_panels(weights, hypers: dict):
    """Named weight matrices, with the radial output split into feature blocks."""
    layout = _layout(hypers)
    features = layout["features"]
    radial = layout["radial"]
    radial_out = weights.get("sr.radial_coefficients.2.weight")
    if radial_out is not None and tuple(radial_out.shape) == (
        radial * features,
        features,
    ):
        radial_out = radial_out.reshape(radial, features, features)
    specs = (
        (
            f"embedding ({layout['species']} species)",
            "sr.chemical_embedding.weight",
            weights.get("sr.chemical_embedding.weight"),
        ),
        (
            "scalar map, species to features",
            "sr.dense0.0.weight",
            weights.get("sr.dense0.0.weight"),
        ),
        (
            "scalar message, features x features",
            "sr.dense1.weight",
            weights.get("sr.dense1.weight"),
        ),
        (
            "features to spherical channels",
            "sr.dense2.weight",
            weights.get("sr.dense2.weight"),
        ),
        (
            "radial MLP, in",
            "sr.radial_coefficients.0.weight",
            weights.get("sr.radial_coefficients.0.weight"),
        ),
        (
            f"radial MLP, {radial} feature blocks",
            "sr.radial_coefficients.2.weight",
            radial_out,
        ),
        (
            "spherical tensor product",
            "sr.tensor_dense.tensor_weight",
            weights.get("sr.tensor_dense.tensor_weight"),
        ),
        (
            "short-range energy readout",
            "sr.energy_mlp.4.weight",
            weights.get("sr.energy_mlp.4.weight"),
        ),
        (
            "scalar charge MLP",
            "lr.scalar_charge_mlp.0.weight",
            weights.get("lr.scalar_charge_mlp.0.weight"),
        ),
        (
            "spherical charges",
            "lr.spherical_charge_dense.dense.weight",
            weights.get("lr.spherical_charge_dense.dense.weight"),
        ),
        (
            "potential to features",
            "lr.potential_to_features.weight",
            weights.get("lr.potential_to_features.weight"),
        ),
        (
            "long-range energy readout",
            "lr.energy_mlp.4.weight",
            weights.get("lr.energy_mlp.4.weight"),
        ),
    )
    panels = []
    for title, key, values in specs:
        if values is not None:
            panels.append((title, key, values.detach().cpu().numpy()))
    return panels


def parameter_figure(panels, *, transform, aspect):
    columns = 4
    rows = int(np.ceil(len(panels) / columns))
    fig, axes = plt.subplots(
        rows,
        columns,
        figsize=(11.4, 2.15 * rows),
        layout="constrained",
    )
    flat = np.atleast_1d(axes).ravel()
    for ax, panel in zip(flat, panels, strict=False):
        _show_matrix(ax, panel[2], panel[0], transform, aspect)
    for ax in flat[len(panels) :]:
        ax.axis("off")
    return fig


def _adam_by_name(weights, optimizer) -> dict:
    state = (optimizer or {}).get("state") or {}
    names = [
        name for name in weights if not any(mark in name for mark in _BUFFER_MARKS)
    ]
    if len(names) != len(state):
        return {}
    return {name: state[index] for index, name in enumerate(names)}


def adam_figure(panels, weights, optimizer, *, transform, aspect):
    """Parameter, first moment, and root second moment for each panel.

    Returns ``(figure, note)``. The figure is ``None`` when this checkpoint
    has no matching Adam state.
    """
    optimizer = optimizer or {}
    by_name = _adam_by_name(weights, optimizer)
    group = (optimizer.get("param_groups") or [{}])[0]
    if not by_name:
        return None, "This checkpoint has no Adam state."
    state = optimizer.get("state") or {}
    betas = group.get("betas", (0.9, 0.999))
    learning_rate = group.get("lr")
    step = int(next(iter(state.values()))["step"]) if state else None
    fig, axes = plt.subplots(
        len(panels),
        3,
        figsize=(11.2, 1.2 * len(panels)),
        layout="constrained",
    )
    axes = np.atleast_2d(axes)
    for row, (title, key, view) in enumerate(panels):
        moment = by_name.get(key)
        _show_matrix(
            axes[row, 0],
            view,
            "parameter" if row == 0 else "",
            transform,
            aspect,
        )
        axes[row, 0].set_ylabel(title, fontsize=7)
        if moment is None:
            axes[row, 1].axis("off")
            axes[row, 2].axis("off")
            continue
        first = moment["exp_avg"].detach().cpu().numpy()
        second = moment["exp_avg_sq"].detach().cpu().numpy()
        if first.size == view.size:
            first = first.reshape(view.shape)
            second = second.reshape(view.shape)
        _show_matrix(
            axes[row, 1],
            first,
            "first moment m" if row == 0 else "",
            transform,
            aspect,
        )
        _show_matrix(
            axes[row, 2],
            np.sqrt(np.maximum(second, 0.0)),
            "root second moment" if row == 0 else "",
            transform,
            aspect,
        )
    rate = f"{learning_rate:.4g}" if isinstance(learning_rate, float) else "—"
    note = (
        f"Adam keeps `exp_avg` and `exp_avg_sq` for each of "
        f"{len(by_name)} parameters, at step {step}. "
        f"betas are {betas[0]} and {betas[1]}, learning rate {rate}. "
        "Each row is the parameter, its first moment, and the root of "
        "its second moment, tiled the same way. "
        "`best_optimizer_state_dict` stores the moments from the best "
        "epoch; these are the moments saved with the weights above."
    )
    return fig, note


def shape_ledger(hypers: dict, result: dict | None = None) -> list[dict]:
    layout = _layout(hypers)
    ledger = [
        {
            "tensor": "nodes_scalar",
            "shape": f"(n_atoms, {layout['features']})",
            "role": "invariant hidden state",
        },
        {
            "tensor": "nodes_spherical",
            "shape": f"(n_atoms, {layout['components']}, {layout['spherical']})",
            "role": (f"(max_degree+1)^2 = {layout['components']} equivariant channels"),
        },
        {
            "tensor": "snapshots",
            "shape": (
                f"({layout['stages']}, n_atoms, "
                f"{layout['components']}, {layout['spherical']})"
            ),
            "role": "spherical state after each message-passing stage",
        },
        {
            "tensor": "radial basis",
            "shape": f"(n_pairs, {layout['radial']})",
            "role": f"Bernstein basis inside the {layout['cutoff']} A cutoff",
        },
        {
            "tensor": "charges",
            "shape": f"(n_atoms, {1 + layout['components_lr']})",
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
    if result is None:
        return ledger
    hidden = result["hidden"]
    for row, array in (
        (ledger[0], hidden["nodes_scalar"]),
        (ledger[1], hidden["nodes_spherical"]),
        (ledger[2], hidden["snapshots"]),
        (ledger[4], hidden["charges"]),
    ):
        row["this structure"] = "×".join(str(size) for size in array.shape)
    ledger[3]["this structure"] = str(result["n_pairs"])
    ledger[5]["this structure"] = str(result["n_atoms"])
    return ledger


def hidden_figure(hidden: dict, *, transform, aspect):
    scalar = hidden["nodes_scalar"]
    equivariant = hidden["nodes_spherical"]
    charges = hidden["charges"]
    updates = hidden["spherical_updates"]
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 6.2), layout="constrained")
    _show_matrix(
        axes[0, 0],
        scalar,
        f"scalar features {scalar.shape}",
        transform,
        aspect,
    )
    _show_matrix(
        axes[0, 1],
        np.transpose(equivariant, (1, 0, 2)),
        f"spherical features, {equivariant.shape[1]} (l, m) blocks",
        transform,
        aspect,
    )
    _show_matrix(
        axes[1, 0],
        charges,
        f"charges {charges.shape}",
        transform,
        aspect,
    )
    _show_matrix(
        axes[1, 1],
        np.transpose(updates, (1, 0, 2)),
        f"long-range spherical update {updates.shape}",
        transform,
        aspect,
    )
    return fig


def hidden_caption(result: dict) -> str:
    return (
        f"**{result['formula']}**, {result['n_atoms']} atoms, "
        f"{result['n_pairs']} neighbor pairs. "
        "Spherical panels are stacked by (l, m). "
        "Training sums the atomic energies and backpropagates "
        "forces and stress; eval mode also adds the composition "
        "baseline and the output scale inside `forward`."
    )


def training_series(runs: list[Path]) -> list[tuple[str, list[dict]]]:
    labels = [_unique_label(run, list(runs)) for run in runs]
    series = []
    for label, run in zip(labels, runs, strict=True):
        metrics = read_metrics(run)
        if metrics:
            series.append((label, metrics))
    return series


def training_overlay(series: list[tuple[str, list[dict]]]):
    """Validation curves, one line per run. A single run keeps train and val."""
    if len(series) == 1:
        return training_figure(series[0][1])
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.1))
    for name, metrics in series:
        epochs = [row["epoch"] for row in metrics]
        axes[0].plot(epochs, [row["val_loss"] for row in metrics], "o-", label=name)
        axes[1].plot(epochs, [row["val_rmse_meV"] for row in metrics], "o-", label=name)
        axes[2].plot(
            epochs, [row["learning_rate"] for row in metrics], "o-", label=name
        )
    axes[0].set_title("Validation loss")
    axes[1].set_title("Validation RMSE (meV/atom)")
    axes[2].set_title("Learning rate")
    for axis in axes:
        axis.set_xlabel("epoch")
        axis.legend(fontsize=8)
    fig.tight_layout()
    return fig


def metrics_rows(series: list[tuple[str, list[dict]]]) -> list[dict]:
    rows = []
    for name, metrics in series:
        for row in metrics:
            item = {"run": name, **row}
            rows.append(item)
    return rows


def run_logs(runs: list[Path]) -> list[dict]:
    sections = []
    labels = [_unique_label(run, list(runs)) for run in runs]
    for label, run in zip(labels, runs, strict=True):
        log_path = run / "train.log"
        options_path = run / "options_restart.yaml"
        log_text = log_path.read_text() if log_path.is_file() else ""
        sections.append(
            {
                "name": label,
                "log": log_text,
                "options": (options_path.read_text() if options_path.is_file() else ""),
                "summaries": dataset_summaries(log_text),
            }
        )
    return sections


def training_figure(metrics: list[dict]):
    epochs = [row["epoch"] for row in metrics]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.1))
    axes[0].plot(epochs, [row["train_loss"] for row in metrics], "o-", label="train")
    axes[0].plot(epochs, [row["val_loss"] for row in metrics], "o-", label="val")
    axes[0].set_title("Loss")
    axes[0].set_xlabel("epoch")
    axes[0].legend()
    axes[1].plot(
        epochs, [row["train_rmse_meV"] for row in metrics], "o-", label="train"
    )
    axes[1].plot(epochs, [row["val_rmse_meV"] for row in metrics], "o-", label="val")
    axes[1].set_title("Energy RMSE (meV/atom)")
    axes[1].set_xlabel("epoch")
    axes[1].legend()
    axes[2].plot(epochs, [row["learning_rate"] for row in metrics], "o-")
    axes[2].set_title("Learning rate")
    axes[2].set_xlabel("epoch")
    fig.tight_layout()
    return fig


def training_headline(metrics: list[dict]) -> dict:
    latest = metrics[-1]
    previous = metrics[-2] if len(metrics) > 1 else latest
    delta = latest["val_rmse_meV"] - previous["val_rmse_meV"]
    return {
        "epoch": latest["epoch"],
        "val_loss": latest["val_loss"],
        "train_loss": latest["train_loss"],
        "loss_direction": (
            "decrease" if latest["val_loss"] <= previous["val_loss"] else "increase"
        ),
        "val_rmse": latest["val_rmse_meV"],
        "delta": delta,
        "rmse_direction": "decrease" if delta <= 0 else "increase",
        "learning_rate": latest["learning_rate"],
    }


def with_l2_delta(weight_rows, other_rows):
    """Copy ``weight_rows`` and add the other checkpoint's L2 and the delta."""
    rows = [dict(row) for row in weight_rows]
    other = {row["name"]: row["l2"] for row in other_rows}
    for row in rows:
        baseline = other.get(row["name"])
        row["l2_other"] = baseline
        row["l2_delta"] = None if baseline is None else row["l2"] - baseline
    return rows


def weight_figure(weight_rows, flats):
    values = np.concatenate(flats)
    low, high = np.percentile(values, [1, 99])
    clipped = values[(values >= low) & (values <= high)]
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.2))
    order = np.argsort([row["l2"] for row in weight_rows])[::-1][:18]
    axes[0].barh(
        [weight_rows[i]["name"] for i in order][::-1],
        [weight_rows[i]["l2"] for i in order][::-1],
    )
    axes[0].set_xscale("log")
    axes[0].set_title("Largest tensor L2 norms")
    axes[1].hist(clipped, bins=80, color="#4c1d95")
    axes[1].set_title("Parameter values (1st–99th percentile)")
    fig.tight_layout()
    return fig
