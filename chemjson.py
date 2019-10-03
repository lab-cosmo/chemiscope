# -*- coding: utf-8 -*-
import json
from chemfiles import Trajectory

system = "Arginine-Dipeptide"

structures = []
with Trajectory("../isv/example/traj-{}.xyz".format(system)) as fd:
    for _ in range(fd.nsteps):
        frame = fd.read()
        positions = frame.positions
        atoms = []
        for i in range(len(frame.atoms)):
            atoms.append({
                "l": frame.atoms[i].name,
                "x": positions[i, 0],
                "y": positions[i, 1],
                "z": positions[i, 2],
            })
        structures.append({"mol": {"a": atoms}})

data = []
with open("../isv/example/{}.dat".format(system)) as fd:
    header = fd.readline().split()
    for col in header:
        data.append([col])

    for line in fd:
        for i, val in enumerate(line.split()):
            data[i].append(val)

sketchmap = {}
for col in data:
    sketchmap[col[0]] = col[1:]

print(json.dumps({
    "sketchmap": sketchmap,
    "structures": structures,
}))
