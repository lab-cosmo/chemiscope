#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import json

from .input import write_input
from .structures import all_atomic_environments, extract_properties


def _chemiscope_input_parser():
    import argparse

    # Tweak the autogenerated help output to look nicer
    def formatter(prog):
        return argparse.HelpFormatter(prog, max_help_position=22)

    parser = argparse.ArgumentParser(
        description=main.__doc__, formatter_class=formatter
    )
    parser.add_argument(
        "input", type=str, help="input file containing the structures and properties"
    )
    parser.add_argument(
        "-o", "--output", type=str, help="chemiscope output file in JSON format"
    )
    parser.add_argument(
        "--properties",
        type=str,
        default="",
        help="comma-separated list of properties that should be extracted. defaults to all",
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--only-atoms",
        action="store_true",
        help="only use per-atom properties from the input file",
    )
    group.add_argument(
        "--only-structures",
        action="store_true",
        help="only use per-structure properties from the input file (default)",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=3.5,
        help="spherical cutoff radius that should be visualized around environments",
    )
    parser.add_argument("--name", default="", type=str, help="name of the dataset")
    parser.add_argument(
        "--description", default="", type=str, help="description of the dataset"
    )
    parser.add_argument(
        "--authors", nargs="*", type=str, default=[], help="list of dataset authors"
    )
    parser.add_argument(
        "--references",
        nargs="*",
        type=str,
        default=[],
        help="list of references for the dataset",
    )
    parser.add_argument(
        "--settings",
        type=str,
        default="",
        help="visualization settings, as a JSON string, following the chemiscope format",
    )
    return parser


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

    try:
        # command-line execution. requires ASE IO module
        import ase.io as ase_io
    except ImportError:
        raise ImportError(
            "chemiscope_input needs ASE modules to parse structure inputs"
        )

    parser = _chemiscope_input_parser()
    args = parser.parse_args()

    if args.only_atoms and args.cutoff is None:
        raise Exception("--only-atoms requires to give --cutoff")
    if args.only_structures and args.cutoff is not None:
        raise Exception("--only-structure can not be given with --cutoff")

    # read file with ASE and remove extraneous properties
    frames = ase_io.read(args.input, ":")
    if args.only_structures:
        for frame in frames:
            for key in list(frame.arrays.keys()):
                if key not in ["positions", "numbers", "momenta"]:
                    del frame.arrays[key]
        environments = None
    elif args.only_atoms:
        for frame in frames:
            frame.info = {}
        environments = all_atomic_environments(frames, cutoff=args.cutoff)

    if args.properties == "":
        filter_properties = None
    else:
        filter_properties = args.properties.split(",")

    # determine output file name automatically if missing
    output = args.output or args.input + "_chemiscope.json.gz"
    properties = extract_properties(frames, only=filter_properties)
    has_environment_property = False
    for prop in properties.values():
        if prop["target"] == "atom":
            has_environment_property = True
            break
    if has_environment_property:
        # assumes properties are associated with all atomic environments
        environments = all_atomic_environments(frames, cutoff=args.cutoff)
    else:
        environments = None

    if args.settings == "":
        settings = None
    else:
        # reads settings as a dictionary, following chemiscope format and schema
        settings = json.loads(args.settings)

    write_input(
        path=output,
        frames=frames,
        properties=properties,
        environments=environments,
        meta={
            "name": args.name,
            "description": args.description,
            "authors": args.authors,
            "references": args.references,
        },
        settings=settings,
    )


if __name__ == "__main__":
    main()
