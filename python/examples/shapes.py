import numpy as np
import ase
import ase.io as aseio
import chemiscope


alphaml = aseio.read("alphaml.xyz", ":")

for a in alphaml:
    a.arrays["alpha"] = np.array([[ axx, ayy, azz, axy, axz, ayz ] for 
                               (axx, ayy, azz, axy, axz, ayz) in 
                         zip(a.arrays["axx"],a.arrays["ayy"],a.arrays["azz"],
                             a.arrays["axy"],a.arrays["axz"],a.arrays["ayz"],) ] )
    a.arrays["alpha-vec"] = np.array([[ axx, ayy, azz] for (axx, ayy, azz) in 
                         zip(a.arrays["axx"],a.arrays["ayy"],a.arrays["azz"]) ] )
    
cs = chemiscope.write_input("shapes.json.gz",
    frames=alphaml,
    properties=chemiscope.extract_properties(alphaml, only=["alpha", "alpha-vec"]),
    shapes={"structure_shape" : dict(
            kind = "sphere",
            parameters = { 
              "global" : {"radius":1.0},
              "structure" : [ 
                  {"position": [0,0,0], "color": 0xff0000},
                  {"position": [1,2,3], "color": 0x00ff00}
              ]
            }
          ), 
          "atom_shape" : dict(
           kind = "sphere",
          parameters = {
              "global" : { "radius": 1.0 },
              "structure": [
                  { "radius": 1.5 },
                  { "radius": 0.7 },
                  ],
              "atom" : [ { "radius": 1+np.random.uniform() } for i in range(40) ],
              }
        )
        }, 
    settings={'structure': [{
   'spaceFilling': False,
   'atomLabels': False,
   'shape': 'structure_shape',
   'axes': 'off',
   'keepOrientation': False,
   'playbackDelay': 700}]},
    environments=chemiscope.all_atomic_environments(alphaml),
    )     
