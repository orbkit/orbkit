General Aspects
===============

This section lists all available input formats from the several quantum chemistry programs and the 
requirements for their proper processing with orbkit.
Besides, we have implemented an interface to cclib_.
This platform can extract the data from additional computational chemistry packages.
At the end of this section, all existing quantities, features and output formats of orbkit are introduced.

.. contents:: Table of Contents:
  :local:
  :depth: 2

Requirements for the Input Files
--------------------------------

In order to operate without problems, the input files of orbkit have to fulfill
some requirements.

Cartesian Harmonic Gaussian basis sets
......................................

Internally, orbkit works with Cartesian Harmonic Gaussian basis sets. Unless 
otherwise stated (cf. :ref:`Central Variables <qc.ao_spec>` for details), it 
assumes the standard Molden basis function order for the exponents 
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

  Unless the exponents are not defined explicitly using "exp_list" in 
  ``qc.ao_spec`` (cf. :ref:`Central Variables <qc.ao_spec>` for details), 
  orbkit is restricted to s, p, d, f, and g atomic orbitals (Molden file 
  limitation).

Real-Valued (Pure) Spherical Harmonic Gaussian basis sets
.........................................................

orbkit supports Spherical Harmonic Gaussian basis sets currently for 
Gaussian .log files (up to g atomic orbitals). After computing the Cartesian Gaussian basis set,
it converts the atomic orbitals to a Spherical Harmonic Gaussian basis. 
The conversion procedure is adapted from

  H.B. Schlegel and M.J. Frisch, *International Journal of Quantum Chemistry*, 
  **54**, 83 (1995).

Notes on the Input File Formats
...............................

Subsequently, you can find a brief overview of all available input file formats 
and some advices for the input file preparation.

**Molden File Format:**

  * Contains Cartesian Harmonics by default
  * Starts with :literal:`[Molden Format]`
  * Contains the sections :literal:`[Atoms]`, :literal:`[GTO]`, :literal:`[MO]`
  * If more than one :literal:`[Molden Format]` keyword is present, orbkit 
    provides an interactive selection.
  * How to create Molden files:

    * MOLPRO: http://www.molpro.net/info/2012.1/doc/manual/node102.html

    .. * TURBOMOLE: tm2molden

**AOMix File Format:**

  * Contains Cartesian Harmonics by default
  * Starts with :literal:`[AOMix Format]`
  * Contains the sections :literal:`[Atoms]`, :literal:`[GTO]`, :literal:`[MO]`
  * If more than one :literal:`[AOMix Format]` keyword is present, orbkit 
    provides an interactive selection.
  * How to create AOMix files:

    * TURBOMOLE: t2aomix


**GAMESS-US Output File:**

  * Please use Cartesian Harmonics (**default**)
  * Hint: GAMESS-US uses a non-standard order of basis functions. Thus, the 
    "exp_list" is explicitly defined in ``qc.ao_spec`` 
    (cf. :ref:`Central Variables <qc.ao_spec>` for details)

**GAUSSIAN .log File:**

  * Spherical Harmonics are chosen by default
  * Use the following parameters in your root section

    .. code-block:: bash

        gfinput IOP(6/7=3)

  * You may switch manually to Cartesian Harmonics using `6D 10F`
  * If more than one "linked" file/geometry/atomic orbitals/molecular orbitals 
    section is present in the .log file, orbkit provides an interactive selection.

**GAUSSIAN Formatted Checkpoint File:**

  * Contains Cartesian Harmonics by default
  * Not applicable for natural orbitals => occupation numbers are not printed
  * How to create FChk files:

    * Add :literal:`%Chk=chkpt-file` to your Gaussian input file
    * Use :literal:`formchk` to convert the chk file:

      .. code-block:: bash

          $ formchk chkpt-file formatted-file

**wfn/wfx Files:**

  * Contains Cartesian Harmonics by default
  * How to create wfn/wfx files:

    * For ORCA: 

      .. code-block:: bash

          $ orca_2aim gbw-basename

Interface to cclib Library
..........................

The cclib_ library is an open source Python package which allows for the parsing and interpreting data of quantum chemistry packages.
It is well checked for multiple versions of different programs.
The interface for cclib_ that we have implemented converts data extracted with cclib into the data structure of orbkit.
A tutorial for the usage of this interface is given in :doc:`./adtutorials/cclib`.


Capabilities
------------

orbkit is designed with a modular structure. This allows to use it not only 
as a standalone version but also to combine its individual  modules or functions 
in user-written Python programs. Each module consists of different functions 
accomplishing specific tasks. 
Thus, there are three ways to use orbkit:

1. As a standalone program via the Terminal (:doc:`./gridbased/terminal`)

2. With a Python script setting options and calling the main function of orbkit (:doc:`./gridbased/highlevel`) 

3. With a user-written Python program using the built-in functions of orbkit (:doc:`./gridbased/lowlevel`)

Detailed tutorials for the three variants are given in the respective sections. 
All grid-based quantities and most of the options can be applied in each of these variants.
The non grid-based quantities are solely available via the low-level interface.
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
| **Mulliken and Löwdin Charges**    |          |bm|          |           |bm|           |           |cm|          |
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

.. _cclib: https://github.com/cclib/cclib
