# Python helpers for chemiscope

This package contains Python code to help generate input files for the
[chemiscope](https://chemiscope.org) default visualizer, and integrate
chemiscope with jupyter notebooks.

## Installation

You should use pip to install this package:

```bash
pip install chemiscope
```

This installs both a `chemiscope-input` command line tool, and the `chemiscope`
package.

## Usage

To create a new chemiscope input file:

```python
import chemiscope
import ase.io

# read frames using ase
frames = ase.io.read("structures.xyz", ":")

# add additional properties to display
properties = {
    "<property name>": {
        target: "atom",
        values: [3, 4, 2, 8, 9, 10],
    }
}

chemiscope.write_input("my-input.json.gz", frames=frames, properties=properties)
```

To display a chemiscope widget inside a jupyter notebook:

```python
import chemiscope
import ase.io

# read frames using ase
frames = ase.io.read("structures.xyz", ":")

# add additional properties to display
properties = {
    "<property name>": [3, 4, 2, 8, 9, 10],
}

chemiscope.show(frames=frames, properties=properties)
```

To create a new chemiscope input file using `stk`_:

.. _`stk`: https://github.com/lukasturcani/stk

`stk` is an optional dependancy:

```bash
pip install stk
```

```python
import chemiscope
import stk

# read frames using ase
frames = [stk.BuildingBlock(smiles="NCCN")]

# add additional properties to display
properties = {
    "<property name>": {
        target: "atom",
        values: [3, 4, 2, 8, 9, 10],
    }
}

chemiscope.write_input(
    "my-input.json.gz",
    frames=frames,
    properties=properties,
)
```