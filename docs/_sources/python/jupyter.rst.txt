.. _jupyter:

Jupyter notebooks
=================

Chemiscope can be used as a widget in Jupyter notebooks, that should work in both Jupyter classic and JupyterLab. 
The widget can be created in `default` mode (showing both a structure and a map panel), or used to display only
structures or only properties. 

Once created, it is possible to interact with the widget using a traitlet interface, modeled after 
`Jupyter widgets <http://ipywidgets.readthedocs.io>`_.

Creating a chemiscope widget
----------------------------

.. autofunction:: chemiscope.show
.. autofunction:: chemiscope.show_input
