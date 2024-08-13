# Chemiscope: interactive structure-property explorer for materials and molecules

[![Tests](https://github.com/lab-cosmo/chemiscope/actions/workflows/tests.yml/badge.svg)](https://github.com/lab-cosmo/chemiscope/actions/workflows/tests.yml)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02117/status.svg)](https://doi.org/10.21105/joss.02117)

Chemiscope is a graphical tool for the interactive exploration of materials and
molecular databases, correlating local and global structural descriptors with
the physical properties of the different systems; as well as a library of
re-usable components useful to create new interfaces.

![Default interface of chemiscope](docs/src/img/screenshot.png)

## Citing chemiscope

Chemiscope is distributed under an open-source license, and you are welcome to
use it and incorporate it into your own research and software projects.
If you find it useful, we would appreciate a citation to the chemiscope
[paper](https://doi.org/10.21105/joss.02117):

> G. Fraux, R. K. Cersonsky, M. Ceriotti, _Chemiscope: Interactive
> Structure-Property Explorer for Materials and Molecules._ **Journal of Open
> Source Software** 5 (51), 2117 (2020)

If you incorporate chemiscope components into a software project, a link back to
the chemiscope homepage (https://chemiscope.org) is the preferred form of
acknowledgement.

## [Documentation](https://chemiscope.org/docs/)

You may be interested in particular about how to [create a visualization of your
own dataset](https://chemiscope.org/docs/manual/input.html).

If you would like to generate a simple chemiscope for your dataset, we
have a [Google Colab notebook](https://colab.research.google.com/drive/1NU0gjtaHcB5Oc3NbFZiQYtctY2190hDu)
that can help!

## Getting help for using chemiscope

If you want to get help when using chemiscope either as a JavaScript/TypeScript
library inside your own project; or for creating input files for the default
visualizer at https://chemiscope.org, you can open a [Github
issue](https://github.com/lab-cosmo/chemiscope/issues/new) with your question;
or send an email to the developers (you can find these emails on the lab
webpage: https://www.epfl.ch/labs/cosmo/people/)

## Getting the python package and using chemiscope in Jupyter notebooks

Using chemiscope in a Jupyter notebook should be as easy as

```bash
pip install chemiscope
```

This also allows to generate chemiscope JSON files that can be viewed on
http://chemiscope.org

If you need to build and install a development version, you should have all the
npm stack installed, and then just run

```bash
git clone https://github.com/lab-cosmo/chemiscope
cd chemiscope
pip install .
```

## Getting and running the web app locally

```bash
git clone https://github.com/lab-cosmo/chemiscope
cd chemiscope
npm install
npm start

# navigate to localhost:8080
```

## Building the code to use it in other projects

```bash
git clone https://github.com/lab-cosmo/chemiscope
cd chemiscope
npm install
npm run build

# Include dist/chemiscope.min.js or dist/molecule-viewer.min.js
# in your own web page
```

See [app/] or the [documentation](https://chemiscope.org/docs/embedding.html)
for an example of how to create a webpage using chemiscope.

## `chemiscope.explore` option

Chemiscope provides a way to automatically explore datasets, using machine
learning representations and dimensionality reduction. For examples and detailed
usage, refer to the related
[documentation](https://chemiscope.org/docs/examples/explore.html).

To use the explore functionality, you'll need to install the required
dependencies with:

```bash
pip install chemiscope[explore]
```

To use `chemiscope.metatensor_featurizer` for providing your trained model
to get the features for `chemiscope.explore`, install the dependencies with:
```bash
pip install chemiscope[metatensor]
```

## sphinx and sphinx-gallery integration

Chemiscope provides also extensions for `sphinx` and `sphinx-gallery` to
include chemiscope viewers within the documentation of a Python package.
See the [documentation](https://chemiscope.org/docs/python/sphinx.html)
for a discussion of the setup and a few examples.

## License and contributions

If you are interested in contributing to chemiscope, please have a look at our
[contribution guidelines](Contributing.md)

Chemiscope itself is distributed under the 3-Clauses BSD license. By
contributing to this repository, you agree to distribute your contributions
under the same license.
