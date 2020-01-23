# -*- coding: utf-8 -*-
'''
Generate JSON input files for the default visualizer using ASE to read
structures.
'''
import numpy as np
import json

import ase
from ase import neighborlist


def _typetransform(data):
    assert isinstance(data, list) and len(data) > 0
    if isinstance(data[0], str):
        return list(map(str, data))
    elif isinstance(data[0], bytes):
        return list(map(lambda u: u.decode('utf8'), data))
    else:
        try:
            return [float(value) for value in data]
        except ValueError:
            raise Exception('unsupported type in property')


def _linearize(name, property):
    '''
    Transform 2D arrays in multiple 1D arrays, converting types to fit json as
    needed.
    '''
    data = {}
    if isinstance(property['values'], list):
        data[name] = {
            'target': property['target'],
            'values': _typetransform(property['values']),
        }
    elif isinstance(property['values'], np.ndarray):
        if len(property['values'].shape) == 1:
            data[name] = {
                'target': property['target'],
                'values': _typetransform(list(property['values'])),
            }
        elif len(property['values'].shape) == 2:
            for i in range(property['values'].shape[1]):
                data[f'{name}[{i + 1}]'] = {
                    'target': property['target'],
                    'values': _typetransform(list(property['values'][:, i])),
                }
        else:
            raise Exception('unsupported ndarray property')
    else:
        raise Exception(f'unknown type for property {name}')

    return data


def _frame_to_json(frame):
    data = {}
    data['names'] = list(frame.symbols)
    data['x'] = [float(value) for value in frame.positions[:, 0]]
    data['y'] = [float(value) for value in frame.positions[:, 1]]
    data['z'] = [float(value) for value in frame.positions[:, 2]]

    if (frame.cell.lengths() != [0.0, 0.0, 0.0]).all():
        data['cell'] = list(np.concatenate(frame.cell))

    return data


def _generate_environments(frames, cutoff):
    environments = []
    for frame_id, frame in enumerate(frames):
        for center in range(len(frame)):
            environments.append({
                'structure': frame_id,
                'center': center,
                'cutoff': cutoff,
            })
    return environments


def write_chemiscope_input(filename, name, frames, extra, cutoff=None):
    '''
    Write the json file expected by the sketchmap vizualizer at ``filename``.

    ---Arguments---
    filename: name of the file to use to save the json data
    name: name of the dataset
    frames: ase.Atoms list containing all the structures
    extra: dictionary of extra properties to display.
        For example:

            extra = {
                'CS': {
                    'target': 'atom',
                    'values': [2, 3, 4, 5, 6, 7]
                }
            }

        'target' can be 'atom' or 'structure', and 'values' can contain either
        numbers or strings
    cutoff: real, optional. If present, will be used to generate atom-centered
            environments, containing all other atoms inside the spherical cutoff
    '''

    data = {
        'meta': {
            'name': name,
        }
    }

    properties = {}
    for name, value in extra.items():
        properties.update(_linearize(name, value))

    from_frames = {}

    # target: structure
    from_frames.update({
        name: {
            'target': 'structure',
            'values': []
        }
        for name in frames[0].info.keys()
    })
    for frame in frames:
        for name, value in frame.info.items():
            from_frames[name]['values'].append(value)

    # target: atom
    from_frames.update({
        name: {
            'target': 'atom',
            'values': value
        }
        for name, value in frames[0].arrays.items() if name != 'positions'
    })
    for frame in frames[1:]:
        for name, value in frame.arrays.items():
            if name == 'positions':
                continue
            from_frames[name]['values'] = np.concatenate([from_frames[name]['values'], value])

    for name, value in from_frames.items():
        properties.update(_linearize(name, value))

    data['properties'] = properties
    data['structures'] = [_frame_to_json(frame) for frame in frames]

    if cutoff is not None:
        data['environments'] = _generate_environments(frames, cutoff)

    with open(filename, 'w') as file:
        json.dump(data, file)
