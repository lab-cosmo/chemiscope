# Chemiscope: interactive structure-property explorer for materials and molecules

Chemiscope is an graphical tool for the interactive exploration of materials and
molecular databases, correlating local and global structural descriptors with
the physical properties of the different systems; as well as a library of
re-usable components useful to create new interfaces.

![Default interface of chemiscope](paper/screenshot.png)

## [Documentation](https://chemiscope.org/docs/)

You may be interested in particular about how to [create a visualization of your
own dataset](https://chemiscope.org/docs/tutorial.html#input-file-format-for-chemiscope).

## Getting and running the code

```bash
git clone https://github.com/cosmo-epfl/chemiscope
cd chemiscope
npm install
npm run download-example-input
npm start

# navigate to localhost:8080
```

## Building the code to use it in other projects

```bash
git clone https://github.com/cosmo-epfl/chemiscope
cd chemiscope
npm install
npm run build

# Include dist/chemiscope.min.js in your own web page
```

See [app/] or the [documentation](https://chemiscope.org/docs/embedding.html)
for a examples of how to create a webpage using chemiscope.

## License and contributions

Chemiscope itself is distributed under the 3-Clauses BSD license. By
contributing to this repository, you agree to distribute your contributions
under the same license.

Chemiscope depends on JSmol, which is distributed under the LGPLv2.1 license.
