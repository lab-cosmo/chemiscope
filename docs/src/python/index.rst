.. _python-module:

Python module
=============

The ``chemiscope`` package provides utilities to prepare JSON input files, and
to visualize them as an :ref:`interactive notebook widget <widget>`. The package
supports multiple structure formats including `ase.Atoms`_, `stk.BuildingBlock`_,
and `MDAnalysis.AtomGroup`_ objects.

.. toctree::
    :maxdepth: 2

    cli
    reference
    widget
    streamlit
    sphinx
    gallery

.. tip::

    If you work with an AI coding assistant, the chemiscope repository ships a
    reusable reference skill under ``docs/skills/``. It teaches the assistant the
    ``chemiscope`` Python API, the input format, the visualization settings and the
    custom shapes, so it can help you write visualizations for your own data,
    including forces, atom-centered properties and custom shapes. The files are
    plain Markdown that any agent can read; point your tool at ``docs/skills/`` (an
    ``AGENTS.md`` pointer at the repository root does this automatically for tools
    that support it). For Claude Code the skill is also exposed via
    ``.claude/skills/chemiscope-python``; copy that folder into your
    ``~/.claude/skills/`` (or a project's ``.claude/skills/``) to use it in your own
    projects.


.. _ase.Atoms: https://ase-lib.org/ase/atoms.html
.. _stk.BuildingBlock: https://stk.readthedocs.io/en/stable/_autosummary/stk.BuildingBlock.html
.. _MDAnalysis.AtomGroup: https://docs.mdanalysis.org/stable/documentation_pages/core/groups.html
