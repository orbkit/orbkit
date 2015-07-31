Requirements for the Input Files
================================

In order to operate without problems, the input files of orbkit have to fulfill
some requirements.

.. contents:: Table of Contents:
  :local:
  :depth: 1

Cartesian Harmonic Gaussian basis sets
--------------------------------------

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
---------------------------------------------------------

orbkit supports Spherical Harmonic Gaussian basis sets currently for 
Gaussian .log files (up to g atomic orbitals). After computing the Cartesian Gaussian basis set,
it converts the atomic orbitals to a Spherical Harmonic Gaussian basis. 
The conversion procedure is adapted from

  H.B. Schlegel and M.J. Frisch, *International Journal of Quantum Chemistry*, 
  **54**, 83 (1995).

Notes on the Input File Formats
-------------------------------

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
