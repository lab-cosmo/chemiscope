.. _examples:

Examples
========

This page lists a few examples of visualizations based on chemiscope. 
Rather than providing a demonstration of a scientifc application, as for the 
examples on `the main page <http://chemiscope.org>`_, these examples
consist in short Python scripts that demonstrate how to make a chemiscope
viewer file that uses some particular features (e.g. atom coloring or shape
visualization). 

.. raw:: html

    <h2>Basic usage</h2>

.. figure:: img/example_base.png
    :align: center

    Simple visualization of molecular data from an ASE input file that
    contains both structure and properties
    
.. raw:: html

   <details>
   <summary style="cursor:default"><b> ▶ Show source code...</b></summary>
   
.. literalinclude:: ../../python/examples/base.py
    :language: python

.. raw:: html

   </details>
   

.. raw:: html

    <h2>Trajectory data</h2>

.. figure:: img/example_md.png
    :align: center

    Trajectory data from a MD simulation of allyl alcohol, including 
    visualization of force data.

.. raw:: html

   <details>
   <summary style="cursor:default"><b> ▶ Show source code...</b></summary>
       
.. literalinclude:: ../../python/examples/trajectory.py
    :language: python

.. raw:: html

   </details>
       

.. raw:: html

    <h2>Structure-property maps</h2>

.. figure:: img/example_pca.png
    :align: center

    PCA map of the atomic environments in a dataset of allyl alcohol configurations. 
    The three group of points correspond to O environments, to the OH-bearing carbon,
    and to the aliphatic carbon environments.

.. raw:: html

   <details>
   <summary style="cursor:default"><b> ▶ Show source code...</b></summary>
   
.. literalinclude:: ../../python/examples/structure_map.py
    :language: python

.. raw:: html

   </details>
   
       
.. raw:: html

    <h2>Color atoms by property</h2>       

.. figure:: img/example_colors.png
    :align: center

    Atom coloring based on the values of a local property

.. raw:: html

   <details>
   <summary style="cursor:default"><b> ▶ Show source code...</b></summary>
   
.. literalinclude:: ../../python/examples/colors.py
    :language: python

.. raw:: html

   </details>
   
   
.. raw:: html

    <h2>Add shapes to visualize data</h2>

.. figure:: img/example_shapes.png
    :align: center

    Visualization of atomic and molecular data using vectors, ellipsoids ... and more

.. raw:: html

   <details>
   <summary style="cursor:default"><b> ▶ Show source code...</b></summary>
   
.. literalinclude:: ../../python/examples/shapes.py
    :language: python

.. raw:: html

   </details>
