Requirements for the Input Files
================================

In order to operate without problems, the input files of orbkit have to fulfill
some requirements.

orbkit requires Cartesian Harmonic Gaussian basis sets. Unless otherwise 
stated (cf. :ref:`Central Variables <qc.ao_spec>` for details), it assumes the 
standard Molden basis function order for the exponents :math:`(l_x,l_y,l_z)`:

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

Subsequently, you can find a brief overview of all available input file formats 
and some advices for the input file preparation.

**Molden File Format:**

  * Contains Cartesian Harmonics by default
  * Starts with :literal:`[Molden Format]`
  * Contains the sections :literal:`[Atoms]`, :literal:`[GTO]`, :literal:`[MO]`
  * If more than one :literal:`[Molden Format]` keyword is present, orbkit 
    provides an interactive selection.
  * How to create Molden Files:

    * MOLPRO: http://www.molpro.net/info/2012.1/doc/manual/node102.html
    * TURBOMOLE: tm2molden


**GAMESS-US Output File:**

  * Please use Cartesian Harmonics (**default**)
  * Hint: GAMESS-US uses a non-standard order of basis functions. Thus, the 
    "exp_list" is explicitly defined in ``qc.ao_spec`` 
    (cf. :ref:`Central Variables <qc.ao_spec>` for details)

**GAUSSIAN Formatted Checkpoint File:**

  * Contains Cartesian Harmonics by default
  * Not applicable for natural orbitals => occupation numbers are not printed
  * How to create FChk files:

    * Add :literal:`%Chk=chkpt-file` to your Gaussian input file
    * Use :literal:`formchk` to convert the chk file:

      .. code-block:: bash

          $ formchk chkpt-file formatted-file

**GAUSSIAN .log File:**

  * Spherical Harmonics are chosen by default! You have to switch **manually** 
    to Cartesian Harmonics (:literal:`6D 10F`)!
  * Use the following parameters in your root section

    .. code-block:: bash

        6D 10F gfinput IOP(6/7=3)

  * If more than one geometry/atomic orbitals/molecular orbitals section is 
    present in the .log file, orbkit provides an interactive selection.

