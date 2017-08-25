.. image:: orbkit.png
   :height: 200 px
   :width: 200 px
   :alt: ORBKIT Logo
   :align: right

ORBKIT
======

ORBKIT is a parallel Python program package for post-processing 
wave function data from output files of quantum chemical programs.

:ORBKIT: Copyright (C) 2017, Gunter Hermann, Vincent Pohl, Axel Schild, and Lukas Eugen Marsoner Steinkasserer.
:Website: https://github.com/orbkit/orbkit

The computational capabilities of ORBKIT range from grid-based quantities, e.g., molecular orbitals or 
electron density, to non grid-based quantities for instance Mulliken population charges or
analytical overlap integrals between molecular orbitals. 
There are several options and features to control the respective calculations like grid types and parameters. 
The required data can be extracted from MOLPRO_ (Molden File Format), 
TURBOMOLE_ (AOMix file format), GAMESS-US_, PROAIMS/AIMPAC (wfn/wfx file format), and Gaussian_ (.log File and Formatted Checkpoint File)
output files. Futhermore, an interface to cclib, a parser for quantum chemical logfiles, is provided.
A complete list of all input and output formats as well as quantities and features of ORBKIT is given in 
:doc:`./general`.

This documentation should serve as an overview of how to operate ORBKIT as a standalone program 
and how to use it in your own Python programs. The latter is facilitated by its modular design.
For a quick introduction to ORBKIT, see the :doc:`./quick` or stick to the example files in 
the ORBKIT program packet.
Tutorials for the calculation of grid-based and non grid-based quantities are given in section 
:doc:`./gridbased/index` and :doc:`./nongridbased/index`.
Supplementary, a detailed description of the existing options, main variables, and functions 
is given in :doc:`./refs/index`.
Finally, some additional applications of ORBKIT are presented in :doc:`./adtutorials/index`. 

**Refactoring** : A major refactoring of the code has recently been undertaken. The most importand changes include:

- Reorganization of AO and MO data into dedicated classes (see  :doc:`./refs/cvars`). 
- Refactoring of readers which amongst other remove sthe need to specify the ``itype`` of a file and allows to read data from compressed (.tar, .tar.gz, tar.bz2) files. 
- The QCinfo class now has a save and read function which allows to save QCinfo instances to file restart from files directly.

**NEW:** :ref:`detCI\@ORBKIT` extends ORBKIT's functionality to multi-determinantal wave functions.

Distribution
------------

The source code and multiple example files of ORBKIT can be freely downloaded from the web page
 
:Website:  https://github.com/orbkit/orbkit

For installation instructions, please check :doc:`./install`.

Citation
--------

If you use ORBKIT in your work, please cite it as follows:

Gunter Hermann, Vincent Pohl, Jean Christophe Tremblay, Beate Paulus, Hans-Christian Hege, and Axel Schild,
"ORBKIT: A Modular Python Toolbox for Cross-Platform Postprocessing of Quantum Chemical Wavefunction Data", 
*J. Comput. Chem.* **2016**, `DOI: 10.1002/jcc.24358`__.

.. note::

  The paper is also freely available on `arXiv <https://arxiv.org/abs/1601.03069>`_.

BibTex (:download:`orbkit.bib`):

.. literalinclude:: orbkit.bib

__ http://dx.doi.org/10.1002/jcc.24358


Contact
-------

The ORBKIT support team, Axel, Gunter, Vincent, and Lukas, welcomes every new
user and will be available to answer your questions. For any change
requests, do not hesitate to contact the ORBKIT support team via

https://github.com/orbkit/orbkit/issues

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
   detci/index
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

