General Aspects
===============

This section lists all available input formats from the several quantum chemistry programs and the 
requirements for their proper processing with ORBKIT.

Alongside to the available program output files, there is an `interface to cclib`_.
This platform can extract the data from additional computational chemistry packages.

At the end of this section, all existing quantities, features and output formats of ORBKIT are introduced.

.. note:: 

    Although, quantum chemistry programs often support multiple output file formats,
    all files do not necessarily have the same quality.
    In order to prevent frustration, here is a list of recommended file formats for different quantum chemistry programs:  

    +------------------------------------+-----------------------------+
    |                                    | **Recommanded File Format** |
    +------------------------------------+-----------------------------+
    | **MOLPRO**                         | `Molden File`_              |
    +------------------------------------+-----------------------------+
    | **GAUSSIAN**                       | `GAUSSIAN .log File`_       |
    +------------------------------------+-----------------------------+
    | **GAMESS-US**                      | `GAMESS-US Output File`_    |
    +------------------------------------+-----------------------------+
    | **TURBOMOLE**                      | `AOMix File`_               |
    +------------------------------------+-----------------------------+
    | **Psi4**                           | `Molden File`_              |
    +------------------------------------+-----------------------------+
    | **ORCA**                           | `Molden File`_              |
    +------------------------------------+-----------------------------+


.. contents:: Table of Contents:
  :local:
  :depth: 2

Supported Input File Formats
----------------------------

Subsequently, you find a brief overview of all available input file formats 
and some advices for the input file preparation.

.. _`Molden File`:

Molden File Format
..................

  * Starts with :literal:`[Molden Format]`
  * Contains the sections :literal:`[Atoms]`, :literal:`[GTO]`, :literal:`[MO]`
  * If more than one :literal:`[Molden Format]` keyword is present, ORBKIT 
    provides an interactive selection.
  * Contains `Cartesian Harmonic Gaussian`_ basis set by default
  * For `Real-Valued (Pure) Spherical Harmonic`_ basis set the following keywords are present: 
    :literal:`[5D]`, :literal:`[7F]`, :literal:`[5D7F]`, :literal:`[9G]`
  * For more information see 
  * How to create Molden files:

    * **MOLPRO**: Here is an example input file :download:`molpro.inp <downloads/molpro.inp>` (http://www.molpro.net/info/2012.1/doc/manual/node102.html)
    * **PSI4**: Here is an example input file :download:`psi4.in <downloads/psi4.in>` (http://www.psicode.org/psi4manual/master/molden.html)
    * **ORCA**: (https://sites.google.com/site/orcainputlibrary/printing-and-visualization)

      .. code-block:: bash

          $ orca_2mkl gbw-basename -molden

.. attention::
  
  For TURBOMOLE, please use the `AOMix File`_ format and **not** the molden file format,
  since here, the norm of the atomic orbitals and the order of molecular orbital coefficients are not consistent.
  
    .. * TURBOMOLE: tm2molden

.. _`AOMix File`:

AOMix File Format
.................

  * Very similar to the molden file format
  * Starts with :literal:`[AOMix Format]`
  * Contains the sections :literal:`[Atoms]`, :literal:`[GTO]`, :literal:`[MO]`
  * If more than one :literal:`[AOMix Format]` keyword is present, ORBKIT 
    provides an interactive selection.
  * Contains `Cartesian Harmonic Gaussian`_ basis set by default
  * How to create AOMix files:

    * **TURBOMOLE**: Run ``t2aomix`` in the TURBOMOLE working directory

.. _`GAUSSIAN .log File`:

GAUSSIAN .log File
..................

  * Use the following parameters in your root section

    .. code-block:: bash

        gfinput IOP(6/7=3)

  * `Real-Valued (Pure) Spherical Harmonic`_ basis set is chosen by default
  * You may switch manually to `Cartesian Harmonic Gaussian`_ basis set using :literal:`6D 10F`
  * If more than one *"linked"* file/geometry/atomic orbitals/molecular orbitals 
    section is present in the .log file, ORBKIT provides an interactive selection.

.. _`GAUSSIAN Formatted Checkpoint File`:

GAUSSIAN Formatted Checkpoint File
..................................

  * Contains `Cartesian Harmonic Gaussian`_ basis set by default
  * Not applicable for natural orbitals => occupation numbers are not printed
  * Labels of the molecular orbitals are also not printed
  * How to create FChk files:

    * **GAUSSIAN**:
      
      * Add :literal:`%Chk=chkpt-file` to your Gaussian input file
      * Use :literal:`formchk` to convert the chk file:

        .. code-block:: bash

            $ formchk chkpt-file formatted-file
    * **PSI4**: Here is an example input file :download:`psi4.in <downloads/psi4.in>` (http://psicode.org/psi4manual/master/fchk.html)


.. _`GAMESS-US Output File`:

GAMESS-US Output File
.....................

  * Please use `Cartesian Harmonic Gaussian`_ basis set (**default**)
  * Hint: GAMESS-US uses a non-standard order of basis functions. Thus, the 
    "exp_list" is explicitly defined in ``qc.ao_spec`` 
    (cf. :ref:`Central Variables <qc.ao_spec>` for details)

.. _`wfn/wfx File`:

wfn/wfx Files
.............

  * Contains `Cartesian Harmonic Gaussian`_ basis set by default
  * How to create wfn/wfx files:

    * **ORCA**: (https://sites.google.com/site/orcainputlibrary/orbital-and-density-analysis)

      .. code-block:: bash

          $ orca_2aim gbw-basename

.. _`interface to cclib`:
          
Interface to cclib Library
..........................

The cclib_ library is an open source Python package which allows for the parsing and interpretation of data stemming from different quantum chemistry packages.
It is well checked for multiple versions of different programs.
The interface for cclib_ that we have implemented converts data extracted with cclib into the data structure of ORBKIT.
A tutorial for the usage of this interface is given in :doc:`./adtutorials/cclib`.

.. hint: 

  The cclib interface is only tested for a few examples. If anything does not work properly, please contact us.
  
Capabilities of ORBKIT
----------------------

ORBKIT is designed with a modular structure. This allows to use it not only 
as a standalone program but also to combine its individual  modules or functions 
into new, user-written Python programs. 

Each module consists of different functions accomplishing specific tasks. 
Thus, there are three ways to use ORBKIT:

1. As a standalone program via the Terminal (:doc:`./gridbased/terminal`)

2. With a Python script setting options and calling the main function of ORBKIT (:doc:`./gridbased/highlevel`) 

3. With a user-written Python program using the built-in functions of ORBKIT (:doc:`./gridbased/lowlevel`)

Detailed tutorials for the three variants are given in the respective sections. 
All grid-based quantities and most of the options can be applied in each of these variants while
non grid-based quantities are solely available via the low-level interface.

The complete list of all quantities, options, and output formats can be seen below.

.. |cm| unicode:: U+2714 .. Check Mark

.. |bm| unicode:: U+2718 .. Ballot Mark




Computable Quantities
.....................

+------------------------------------+------------------------+--------------------------+-------------------------+
| **Quantity**                       | **Usage via Terminal** | **High-Level Interface** | **Low-Level Interface** |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **Electron Density**               |          |cm|          |           |cm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **Atomic and Molecular Orbitals**  |          |cm|          |           |cm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **Orbital Derivatives**            |          |cm|          |           |cm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **Spin Density**                   |          |cm|          |           |cm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **Gross Atomic Density**           |          |cm|          |           |cm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **Total Dipole Moment**            |          |bm|          |           |bm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **Nuclear Dipole Moment**          |          |bm|          |           |bm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **Mulliken and Lowdin Charges**    |          |bm|          |           |bm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **Center of Charge and Mass**      |          |bm|          |           |bm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+

Options and Features
....................

+------------------------------------+------------------------+--------------------------+-------------------------+
| **Grid Options**                   | **Usage via Terminal** | **High-Level Interface** | **Low-Level Interface** |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **Cartesian Equidistant Grid**     |          |cm|          |           |cm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **Spherical Equidistant Grid**     |          |cm|          |           |cm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **Arbitrary Vector Grid**          |          |cm|          |           |cm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **Random Grid**                    |          |bm|          |           |bm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **Adaption to Molecule Structure** |          |cm|          |           |cm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **Center Grid around Nuclei**      |          |cm|          |           |cm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **Symmetry Operations on Grid**    |          |bm|          |           |bm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+

+------------------------------------+------------------------+--------------------------+-------------------------+
| **Special Features**               | **Usage via Terminal** | **High-Level Interface** | **Low-Level Interface** |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **Multiple File Handling**         |          |bm|          |           |bm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **Ordering of Molecular Orbitals** |          |bm|          |           |bm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **Interpolation**                  |          |bm|          |           |bm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **MO Transition Flux Density**     |          |bm|          |           |bm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+

Output Formats
..............

+------------------------------------+------------------------+--------------------------+-------------------------+
| **Output Formats**                 | **Usage via Terminal** | **High-Level Interface** | **Low-Level Interface** |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **HDF5 Files**                     |          |cm|          |           |cm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **Gaussian Cube Files**            |          |cm|          |           |cm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **VMD Script Files**               |          |cm|          |           |cm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **ZIBAmira Mesh Files**            |          |cm|          |           |cm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **ZIBAmira Network Files**         |          |cm|          |           |cm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **Mayavi Visualization**           |          |cm|          |           |cm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+
| **XYZ and PBE Files**              |          |bm|          |           |bm|           |           |cm|          |
+------------------------------------+------------------------+--------------------------+-------------------------+

Notes on Gaussian Basis Sets
----------------------------

In modern quantum chemistry for finite systems, there are two widely used basis set types:
Cartesian harmonic Gaussian basis sets and real-valued (pure) spherical harmonic Gaussian basis sets.
While ORBKIT internally uses the former type, it is able to handle the latter using a transformation.

.. _`Cartesian Harmonic Gaussian`:

Cartesian Harmonic Gaussian Basis Sets
......................................

Internally, ORBKIT works with cartesian harmonic Gaussian basis sets. Unless the exponents are explicitly defined using "exp_list" 
(cf. :ref:`Central Variables <qc.ao_spec>` for details), it assumes the standard Molden basis function order for the exponents 

:math:`(l_x,l_y,l_z)`:

  * s:
    (0,0,0)
  * p:
    (1,0,0), (0,1,0), (0,0,1)
  * d: 
    (2,0,0), (0,2,0), (0,0,2),
    (1,1,0), (1,0,1), (0,1,1)
  * f: 
    (3,0,0), (0,3,0), (0,0,3), 
    (1,2,0), (2,1,0), (2,0,1), 
    (1,0,2), (0,1,2), (0,2,1), 
    (1,1,1)
  * g: 
    (4,0,0), (0,4,0), (0,0,4), 
    (3,1,0), (3,0,1), (1,3,0), 
    (0,3,1), (1,0,3), (0,1,3), 
    (2,2,0), (2,0,2), (0,2,2), 
    (2,1,1), (1,2,1), (1,1,2)

.. hint:: 

  Also note that ORBKIT is restricted to s, p, d, f, and g atomic orbitals (Molden file 
  limitation).

.. _`Real-Valued (Pure) Spherical Harmonic`:

Real-Valued (Pure) Spherical Harmonic Gaussian basis sets
.........................................................

ORBKIT supports Spherical Harmonic Gaussian basis sets currently up to g atomic orbitals. 
After computing the Cartesian Gaussian basis set,
it converts the atomic orbitals to a Spherical Harmonic Gaussian basis. 
The conversion procedure is adapted from

  H.B. Schlegel and M.J. Frisch, *International Journal of Quantum Chemistry*, 
  **54**, 83 (1995).

.. _cclib: https://github.com/cclib/cclib
