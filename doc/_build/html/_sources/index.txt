.. image:: orbkit.png
   :height: 200px
   :width: 200 px
   :alt: orbkit Logo
   :align: right

orbkit
======

orbkit is a parallel Python program package for post-processing 
wave function data from output files of quantum chemical programs.

:orbkit: Copyright (C) 2015, Gunter Hermann, Vincent Pohl, and Axel Schild.
:Website: http://sourceforge.net/p/orbkit

The computational capabilities of orbkit range from grid-based quantities, e.g., molecular orbitals or 
electron density, to non grid-based quantities for instance Mulliken population charges or
analytical overlap integrals between molecular orbitals. 
There are several options and features to control the respective calculations like grid types and parameters. 
The required data can be extracted from MOLPRO_ (Molden File Format), 
TURBOMOLE_ (AOMix file format), GAMESS-US_, and Gaussian_ (.log File and Formatted Checkpoint File)
output files.
A complete list of all input and output formats as well as quantities and features of orbkit is given in 
:doc:`./general`.

This documentation should serve as an overview of how to operate orbkit as a standalone program 
and how to use it in your own Python programs. The latter is facilitated by its modular design.
For a quick introduction to orbkit, see the :doc:`./quick` or stick to the example files in 
the orbkit program packet.
Tutorials for the calculation of grid-based and non grid-based quantities are given in section 
:doc:`./gridbased/index` and :doc:`./nongridbased/index`.
Supplementary, a detailed description of the existing options, main variables, and functions 
is given in :doc:`./refs/index`.
Finally, some additional applications of orbkit are presented in :doc:`./adtutorials/index`.


Distribution
------------

The source code and multiple example files of orbkit can be freely downloaded from the web page
 
:Website: http://sourceforge.net/p/orbkit

For installation instructions, please check :doc:`./install`.

Citation
--------

If you use orbkit in your work, please cite it as follows:

Gunter Hermann, Vincent Pohl, Axel Schild, "orbkit: A Toolbox for Post-Processing 
Quantum Chemical Wavefunction Data." available via http://sourceforge.net/p/orbkit (2015).

Contact
-------

The orbkit support team, Axel, Gunter, and Vincent, welcomes every new
user and will be available to answer your questions. For any change
requests, do not hesitate to contact the orbkit support team via

:Support: http://sourceforge.net/p/orbkit/support

Contents:
---------

.. toctree::
   :maxdepth: 2

   install
   general
   quick
   gridbased/index
   nongridbased/index   
   adtutorials/index
   refs/index
   licence
   impressum

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _MOLPRO: https://www.molpro.net/
.. _TURBOMOLE: http://www.turbomole.com/
.. _Gamess-US: http://www.msg.chem.iastate.edu/gamess/
.. _Gaussian: http://www.gaussian.com/
.. _ZIBAmira: http://amira.zib.de/

