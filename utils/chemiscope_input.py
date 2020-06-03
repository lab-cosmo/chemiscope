#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Generate JSON input files for the default chemiscope visualizer.
'''
import warnings
import numpy as np
import json
import gzip

IGNORED_ASE_ARRAYS = ['positions', 'numbers']


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
    data['size'] = len(frame)
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


def write_chemiscope_input(filename, frames, meta=None, extra=None, cutoff=None):
    '''
    Write the json file expected by the default chemiscope visualizer at
    ``filename``.

    :param str filename: name of the file to use to save the json data. If it
                         ends with '.gz', a gzip compressed file will be written
    :param list frames: list of `ase.Atoms`_ objects containing all the
                        structures
    :param dict meta: optional metadata of the dataset, see below
    :param dict extra: optional dictionary of additional properties, see below
    :param float cutoff: optional. If present, will be used to generate
                         atom-centered environments

    The dataset metadata should be given in the ``meta`` dictionary, the
    possible keys are:

    .. code-block:: python

        meta = {
            'name': '...',         # str, dataset name
            'description': '...',  # str, dataset description
            'authors': [           # list of str, dataset authors, OPTIONAL
                '...',
            ],
            'references': [        # list of str, references for this dataset,
                '...',             # OPTIONAL
            ],
        }

    The written JSON file will contain all the properties defined on the
    `ase.Atoms`_ objects. Values in ``ase.Atoms.arrays`` are mapped to
    ``target = "atom"`` properties; while values in ``ase.Atoms.info`` are
    mapped to ``target = "structure"`` properties. The only exception is
    ``ase.Atoms.arrays["numbers"]``, which is always ignored. If you want to
    have the atomic numbers as a property, you should add it to ``extra``
    manually.

    Additional properties can be added with the ``extra`` parameter. This
    parameter should be a dictionary containing one entry for each property.
    Each entry contains a ``target`` attribute (``'atom'`` or ``'structure'``)
    and a set of values. ``values`` can be a Python list of float or string; a
    1D numpy array of numeric values; or a 2D numpy array of numeric values. In
    the later case, multiple properties will be generated along the second axis.
    For example, passing

    .. code-block:: python

        extra = {
            'cheese': {
                'target': 'atom',
                'values': np.zeros((300, 4))
            }
        }

    will generate four properties named ``cheese[1]``, ``cheese[2]``,
    ``cheese[3]``,  and ``cheese[4]``, each containing 300 values.

    .. _`ase.Atoms`: https://wiki.fysik.dtu.dk/ase/ase/atoms.html
    '''

    if not (filename.endswith('.json') or filename.endswith('.json.gz')):
        raise Exception('filename should end with .json or .json.gz')

    data = {'meta': {}}

    if meta is not None:
        if 'name' in meta:
            data['meta']['name'] = str(meta['name'])

        if 'description' in meta:
            data['meta']['description'] = str(meta['description'])

        if 'authors' in meta:
            data['meta']['authors'] = list(map(str, meta['authors']))

        if 'references' in meta:
            data['meta']['references'] = list(map(str, meta['references']))

        for key in meta.keys():
            if key not in ['name', 'description', 'authors', 'references']:
                warnings.warn('ignoring unexpected metadata: {}'.format(key))

    if 'name' not in data['meta'] or not data['meta']['name']:
        data['meta']['name'] = filename

    properties = {}
    if extra is not None:
        for name, value in extra.items():
            properties.update(_linearize(name, value))

    # Read properties coming from the ase.Atoms objects
    from_frames = {}

    # target: structure properties
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

    # target: atom properties
    from_frames.update({
        name: {
            'target': 'atom',
            'values': value
        }
        for name, value in frames[0].arrays.items() if name not in IGNORED_ASE_ARRAYS
    })
    for frame in frames[1:]:
        for name, value in frame.arrays.items():
            if name in IGNORED_ASE_ARRAYS:
                continue
            from_frames[name]['values'] = np.concatenate([from_frames[name]['values'], value])

    for name, value in from_frames.items():
        properties.update(_linearize(name, value))

    data['properties'] = properties
    data['structures'] = [_frame_to_json(frame) for frame in frames]

    if cutoff is not None:
        data['environments'] = _generate_environments(frames, cutoff)

    if filename.endswith(".gz"):
        with gzip.open(filename, 'w', 9) as file:
            file.write(json.dumps(data).encode("utf8"))
    else:
        with open(filename, 'w') as file:
            json.dump(data, file)


def main():
    """
    Command-line utility to generate an input for chemiscope — the interactive
    structure-property explorer. Parses an input file containing atomic
    structures using the ASE I/O module, and converts it into a JSON file that
    can be loaded in chemiscope. Frame and environment properties must be
    written in the same file containing atomic structures: we recommend the
    extended xyz format, which is flexible and simple. In all cases, this
    utility will simply write to the JSON file anything that is readable by
    ASE.
    """

    # command-line execution. requires ASE IO module
    import argparse
    try:
        import ase.io as ase_io
    except ImportError:
        raise ImportError("chemiscope_input needs ASE modules to parse structure inputs")

    # Tweak the autogenerated help output to look nicer
    formatter = lambda prog: argparse.HelpFormatter(prog, max_help_position=22)
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=formatter)
    parser.add_argument('input', type=str,
                        help='input file containing the structures and properties')
    parser.add_argument('-o', '--output', type=str,
                        help='chemiscope output file in JSON format')
    parser.add_argument('-c', '--cutoff', type=float,
                        help='generate atom-centred environments with the given cutoff')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--only-atoms', action="store_true",
                       help='only use per-atom properties from the input file')
    group.add_argument('--only-structures', action="store_true",
                       help='only use per-structure properties from the input file')
    parser.add_argument('--name', default="", type=str,
                        help='name of the dataset')
    parser.add_argument('--description', default="", type=str,
                        help='description of the dataset')
    parser.add_argument('--authors', nargs='*', type=str, default=[],
                        help='list of dataset authors')
    parser.add_argument('--references', nargs='*', type=str, default=[],
                        help='list of references for the dataset')
    args = parser.parse_args()

    if args.only_atoms and args.cutoff is None:
        raise Exception("--only-atoms requires to give --cutoff")
    if args.only_structures and args.cutoff is not None:
        raise Exception("--only-structure can not be given with --cutoff")

    # read file with ASE and remove extraneous properties
    frames = ase_io.read(args.input, ':')
    if args.only_structures:
        for frame in frames:
            for key in list(frame.arrays.keys()):
                if key not in IGNORED_ASE_ARRAYS:
                    del frame.arrays[key]
    elif args.only_atoms:
        for frame in frames:
            frame.info = {}

    # determine output file name automatically if missing
    output = args.output or args.input + "_chemiscope.json.gz"

    write_chemiscope_input(
        filename=output,
        frames=frames,
        meta={
            "name": args.name,
            "description": args.description,
            "authors": args.authors,
            "references": args.references
        },
        extra=None,
        cutoff=args.cutoff
    )


if __name__ == '__main__':
    main()
