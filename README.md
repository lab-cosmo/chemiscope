# Chemiscope: interactive structure-property explorer for materials and molecules

![tests](https://github.com/lab-cosmo/chemiscope/workflows/Tests%20&%20Lints/badge.svg)
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
own dataset](https://chemiscope.org/docs/tutorial/input.html).

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

## Getting and running the code

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
for a examples of how to create a webpage using chemiscope.

## License and contributions

If you are interested in contributing to chemiscope, please have a look at our
[contribution guidelines](Contributing.md)

Chemiscope itself is distributed under the 3-Clauses BSD license. By
contributing to this repository, you agree to distribute your contributions
under the same license.
