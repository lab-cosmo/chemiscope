import numpy as np
import ase.io as aseio
import chemiscope


alphaml = aseio.read("alphaml.xyz", ":")

for a in alphaml:
    a.arrays["alpha"] = np.array(
        [
            [axx, ayy, azz, axy, axz, ayz]
            for (axx, ayy, azz, axy, axz, ayz) in zip(
                a.arrays["axx"],
                a.arrays["ayy"],
                a.arrays["azz"],
                a.arrays["axy"],
                a.arrays["axz"],
                a.arrays["ayz"],
            )
        ]
    )
    a.arrays["alpha-vec"] = np.array(
        [
            [axx, ayy, azz]
            for (axx, ayy, azz) in zip(
                a.arrays["axx"], a.arrays["ayy"], a.arrays["azz"]
            )
        ]
    )

cs = chemiscope.write_input(
    "shapes.json.gz",
    frames=alphaml,
    properties=chemiscope.extract_properties(alphaml, only=["alpha", "alpha-vec"]),
    shapes={
        "cube_smooth": dict(
            kind="custom",
            parameters={
                "global": {
                    "vertices": [
                        [0, 0, 0],
                        [1, 0, 0],
                        [0, 1, 0],
                        [1, 1, 0],
                        [0, 0, 1],
                        [1, 0, 1],
                        [0, 1, 1],
                        [1, 1, 1],
                    ],
                    "simplices": [
                        [0, 2, 1],
                        [1, 2, 3],
                        [4, 6, 5],
                        [5, 6, 7],
                        [0, 1, 4],
                        [0, 4, 2],
                        [1, 3, 5],
                        [2, 6, 3],
                        [1, 5, 4],
                        [2, 4, 6],
                        [3, 7, 5],
                        [3, 6, 7],
                    ],
                },
                "structure": [
                    {"position": [0, 4, 0], "color": 0xFF0000},
                    {"position": [3, 2, 1], "color": 0x00FF00},
                ],
            },
        ),
        "cube": dict(
            kind="custom",
            parameters={
                "global": {
                    "vertices": [
                        [0, 0, 0],
                        [0, 1, 0],
                        [1, 1, 0],
                        [1, 0, 0],
                        [0, 0, 0],
                        [0, 0, 1],
                        [0, 1, 1],
                        [0, 1, 0],
                        [0, 1, 0],
                        [0, 1, 1],
                        [1, 1, 1],
                        [1, 1, 0],
                        [1, 1, 0],
                        [1, 1, 1],
                        [1, 0, 1],
                        [1, 0, 0],
                        [1, 0, 0],
                        [1, 0, 1],
                        [0, 0, 1],
                        [0, 0, 0],
                        [0, 0, 1],
                        [1, 0, 1],
                        [1, 1, 1],
                        [0, 1, 1],
                    ],
                    "simplices": [
                        [0, 1, 2],
                        [2, 3, 0],
                        [4, 5, 6],
                        [6, 7, 4],
                        [8, 9, 10],
                        [10, 11, 8],
                        [12, 13, 14],
                        [14, 15, 12],
                        [16, 17, 18],
                        [18, 19, 16],
                        [20, 21, 22],
                        [22, 23, 20],
                    ],
                },
                "structure": [
                    {"position": [0, 4, 0], "color": 0xFF0000},
                    {"position": [3, 2, 1], "color": 0x00FF00},
                ],
            },
        ),
        "structure_shape": dict(
            kind="sphere",
            parameters={
                "global": {"radius": 1.0},
                "structure": [
                    {"position": [0, 0, 0], "color": 0xFF0000},
                    {"position": [1, 2, 3], "color": 0x00FF00},
                ],
            },
        ),
        "atom_shape": dict(
            kind="sphere",
            parameters={
                "global": {"radius": 1.0},
                "structure": [
                    {"radius": 1.5},
                    {"radius": 0.7},
                ],
                "atom": [{"radius": 1 + np.random.uniform()} for i in range(39)]
                + [{"color": 0xC08000}],
            },
        ),
        "alpha": chemiscope.extract_tensors_from_ase(
            alphaml, "alpha", force_positive=True, scale=0.2
        ),
    },
    settings={
        "structure": [
            {
                "spaceFilling": False,
                "atomLabels": False,
                "shape": "cube",
                "axes": "off",
                "keepOrientation": False,
                "playbackDelay": 700,
            }
        ]
    },
    environments=chemiscope.all_atomic_environments(alphaml),
)
