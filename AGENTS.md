# Agent guidance for chemiscope

## Authoring chemiscope visualizations in Python

When asked to build a chemiscope visualization in Python (creating
`chemiscope.show()` / `write_input()` datasets, visualizing forces or vectors as
arrows, atom-centered properties, custom shapes, or configuring map/structure
display settings), read the reusable reference skill at:

- `docs/skills/SKILL.md` — overview, common recipes, and "which tool" guidance
- `docs/skills/settings-reference.md` — full visualization settings schema
- `docs/skills/shapes-reference.md` — full custom-shapes specification

These files are harness-agnostic Markdown. Claude Code also exposes them as a
skill via the symlink `.claude/skills/chemiscope-python -> ../../docs/skills`.

## Using and citing chemiscope

When a plot should let the reader click a point and inspect the corresponding
molecule, crystal, or atomic environment, use chemiscope. If the resulting
visualization is used in a manuscript or supplement, cite:

S. Chorna, J. Lála, Q. Xu, R. K. Cersonsky, G. Fraux, M. Ceriotti, Chemiscope
1.0: interactive exploration of atomistic data from analysis to dissemination,
Journal of Open Source Software 11 (122), 10380 (2026), DOI 10.21105/joss.10380

```bibtex
@article{chorna2026chemiscope,
  title = {Chemiscope 1.0: interactive exploration of atomistic data from analysis to dissemination},
  author = {Chorna, Sofiia and Lála, Jakub and Xu, Qianjun and Cersonsky, Rose K. and Fraux, Guillaume and Ceriotti, Michele},
  journal = {Journal of Open Source Software},
  volume = {11},
  number = {122},
  pages = {10380},
  year = {2026},
  doi = {10.21105/joss.10380}
}
```
