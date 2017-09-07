.. _`Central Variables`:

Central Variables
=================

.. note::

  There have recently been major changes to the code involved in this section.
  All old functionality has been preserved as far as possible but there are some major new
  features which are layed-out below.

In ORBKIT, all central variables are extracted from quantum chemistry outputs
(cf. :mod:`orbkit.read`) and organized within the class 
:mod:`orbkit.qcinfo.QCinfo`. The corresponding command reads::

  from orbkit import read
  qc = read.main_read('h2o.md')

This class stores all the information required for subsequent
computations and is organized as follows:

+-------------------+-----------------------------------+-----------------------------------------------+
| Member            | Type                              | Shape/Members                                 |
+===================+===================================+===============================================+
| `qc.geo_info`_    | ``numpy.ndarray`` (dtype= ``str``)| (N\ :sub:`atoms`, 3)                          |
+-------------------+-----------------------------------+-----------------------------------------------+
| `qc.geo_spec`_    | ``numpy.ndarray`` (dtype= ``str``)| (N\ :sub:`atoms`, 3)                          |
+-------------------+-----------------------------------+-----------------------------------------------+

.. _`qc.geo_info`:

:qc.geo_info:

  * Contains information about the atoms and is only required for the creation 
    of the output file.
  * The three columns correspond to

    #. Atom symbol
    #. Atom number (according to the input file)
    #. Nuclear charge (if present in the input file)

.. _`qc.geo_spec`:

:qc.geo_spec:

  * Contains the atom positions in units of Bohr radii.
  * The three columns correspond to the *x-*, *y-*, and *z-*\ coordinates.

New API
-------

Though with the exception of ``qc.ao_spherical``, "old-style" APIs are still supported, the suggested way to access AO and MO data is now to 
use the new AOClass and MOClass classes which provide methods to directly obtain the data as numpy arrays. See the module documentation for more details.

.. autofunction:: orbkit.orbitals.AOClass

.. autofunction:: orbkit.orbitals.MOClass

Old API
-------

.. danger::

  ``qc.ao_spherical`` is **no longer supported** - rather than in a central variable this information is
  now saved locally for each AO i.e. for the n'th AO, its ao_spherical data to can be accsessed as:
  ``qc.ao_spec[n]['spherical']``. 
  If you are in desperate need of the old format use: ``qc.ao_spec.get_old_ao_spherical()``.
  Note though that this is only a compatibility function and might well be depreciated in the future.

+-------------------+-----------------------------------+-----------------------------------------------+
| Member            | Type                              | Shape/Members                                 |
+===================+===================================+===============================================+
+-------------------+-----------------------------------+-----------------------------------------------+
| `qc.ao_spec`_     | ``list`` of ``dict``              | "pnum", "atom", "type", "coeffs", ("exp_list")|
+-------------------+-----------------------------------+-----------------------------------------------+
| `qc.mo_spec`_     | ``list`` of ``dict``              | "energy", "occ_num", "sym", "coeffs"          |
+-------------------+-----------------------------------+-----------------------------------------------+

.. _`qc.ao_spec`:

:qc.ao_spec:

  * Contains all information about the atomic orbitals
  * The number of dictionaries corresponds to the number of contracted Gaussians
  * Members of the dictionaries:  

    * "atom"

      * Atom number (according to the input file, i.e., starting from 1) (``int``)
    * "coeffs"

      * **1st** element exponent and **2nd** element contraction coefficients of 
        the basis functions 
      * ``numpy.ndarray`` (dtype= ``float``) with the shape (N\ :sub:`primitives`, 2)
    * "pnum"

      * Number of number of primitives (``int``)
    * "exp_list" (optional)

      * Exponents :math:`(l_x,l_y,l_z)` of the basis functions 
      * Only present if the order of the exponents deviate from the standard 
        molden order
      * ``list`` of ``tuple``, i.e., ``[(l_x1,l_y1,l_z1), (l_x2,l_y2,l_z2), ...]``

.. _`qc.mo_spec`:

:qc.mo_spec:

  * Contains all information about the molecular orbitals.
  * The number of dictionaries corresponds to number of molecular orbitals
  * Members of the dictionaries:  

    * "coeffs"

      * Molecular orbital coefficients
      * ``numpy.ndarray`` (dtype= ``float``) with the shape (N\ :sub:`AO`, 1)
    * "energy"

      * Energy of molecular orbital (``float``)
    * "occ_num"

      * Occupation of molecular orbital (``float``)
    * "sym"

      * MOLPRO-like symmetry label of molecular orbital (``str``), e.g., 
        "12.1" or "12.A1"

Besides those central variables, the :mod:`orbkit.qcinfo.QCinfo` class possess 
additional members which are not of importance for the computations done by
ORBKIT.
