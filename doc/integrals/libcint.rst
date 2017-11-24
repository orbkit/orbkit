Libcint Interface
=================

The ``libcint_interface`` module provides a class ``AOIntegrals``, that gives access to various integrals calculated from the information available in a QCinfo object. Make sure libcint is installed  and added to the PATH environment variable (see :ref:`Install libcint`). First create an instance of ``AOIntegrals`` from a QCinto object, e.g.:

.. code-block:: python

    >> from orbkit import read, libcint_interface
    >> qc = read.main_read('h2o.molden', all_mo=True)
    >> ao = integrals.AOIntegrals(qc)

By default, the Integrals are calculated by using cartesian Gaussians. Switching to real spherical Gaussians is possible by disabling the ``cartesian`` flag:

.. code-block:: python

    >> ao = libcint_interface.AOIntegrals(qc, cartesian=False)

Integrals can then be calculated by calling ``ao.int1e(operator)`` or ``ao.int2e(operator)`` for 1- and 2-electron integrals respectively, and passing the desired ``operator`` base name. A list of available operators can be found in the `libcint documentation`_.

.. _libcint documentation: https://github.com/sunqm/libcint/tree/master/doc

There is a small list of shortcuts available:

.. code-block:: python

    >> S = ao.overlap()     # atomic orbital overlap;         = ao.int1e('ovlp')
    >> T = ao.kinetic()     # electron kinetic energy;        = ao.int1e('kin')
    >> Vne = ao.Vne()       # electron-nuclear interaction;   = ao.int1e('nuc')
    >> Hcore = ao.Hcore()   # = T + Vne;                      = ao.int1e('kin') + ao.int1e('nuc')
    >> Vee = ao.Vee()       # 2-electron repulsion integrals; = ao.int2e('')

By default, the integrals are calculated in the basis of atomic orbitals (AO). Transformation into the molecular orbital (MO) basis is requested by passing ``asMO=True``, e.g.:

.. code-block:: python

    >> S = ao.overlap(asMO=True)     # = numpy.eye(ao.Norb)

Requesting subsets of orbitals
-------------------------------

If the basis set is large, calculating all integrals takes time and requires a lot of memory to store them. Especially when requesting integrals in the MO basis, which requires to calculate all integrals in the AO basis first.

If not all integrals are required, it is possible to specify the subset needed. *Before* starting the integral calculation, you can define multiple blocks of orbitals:

.. code-block:: python

    >> ao.add_AO_block_1e(AOrangei=AOrangei, AOrangej=AOrangej)
    >> ao.add_AO_block_2e(AOrangei=AOrangei, AOrangej=AOrangej, AOrangek=AOrangek, AOrangel=AOrangel)
    >> ao.add_MO_block_1e(MOrangei=MOrangei, MOrangej=MOrangej)
    >> ao.add_MO_block_2e(MOrangei=MOrangei, MOrangej=MOrangej, MOrangek=MOrangek, MOrangel=MOrangel)

where the different range parameters define the orbitals for each index. Those can be a Python range object (e.g. ``range(3, 6)``), a list of orbital indices (e.g. ``[1, 3, 7]``), or ``None`` (default) for all orbitals.
The integral calculation will then return a list of matrices, if more than one block is requested.

An example:

.. code-block:: python

    >> ao.add_MO_block_1e(MOrangei=range(10), MOrangej=range(10))
    >> ao.add_MO_block_1e(MOrangei=range(10), MOrangej=range(10, 20))
    >> S = ao.overlap(asMO=True)
    >> S[0] == numpy.eye(10)
    >> S[1] == numpy.zeros((10,10))

For MO integrals, the code automatically selects only those AO integrals which will contribute to the requested MO integrals (based on the MO coefficients).
In connection with symmetry this may save save some time and memory.

To reset the blocks use one of the following methods:

.. code-block:: python

    >> ao.clear_all_blocks()
    >> ao.clear_AO_blocks_1e()
    >> ao.clear_AO_blocks_2e()
    >> ao.clear_MO_blocks_1e()
    >> ao.clear_MO_blocks_2e()

AO slicing when calculating MO integrals
----------------------------------------

For calculating MO integrals all relevant AO integrals are required first. If the basis set is large, the amount of AO integrals to be stored may exceed memory limitations of the hardware. This can be bypassed by calculating the AO integrals in slices (of shells) along the first index (:math:`i`). Each slice is calculated, transformed in MO basis and added to the final result separately, which means only one slice at a time needs to be kept in memory. However, matrix properties (hermitian and exchange of electronic coordinates) can only be exploited within a slice. This means a decent amount of integrals needs to be recalculated. Thus: **Smaller slices require less memory, but take more time to calculate.**

Slicing is controlled by the ``max_dims`` parameter, and primarily usefull for 2-electron integrals with large basis sets and limited memory:

.. code-block:: python

    >> Vee = ao.Vee(asMO=True, max_dims=10)

The default value for ``max_dims`` is 0 and disables slicing. Values larger than that specify the maximum number of basis functions in a slice (along index :math:`i`).

For 2-electron integrals, you may also pass the parameter ``max_mem``, which defines a rough memory limit (in MB) for the AO matrix and ao2mo transformation. The appropiate value for ``max_dims`` is then determined automatically.

Alternatively, if you can exploit symmetry, consider setting only a limited number of MO blocks at a time, to reduce the number of AO integrals to be kept in memory:

.. code-block:: python

    >> ao.add_MO_block_2e(AOMOrangei=irrep1, MOrangej=irrep1, MOrangek=irrep1, MOrangel=irrep1)
    >> Vee[numpy.ix_(irrep1, irrep1, irrep1, irrep1)] = ao.Vee(asMO=True)
    >> ao.clear_MO_blocks_2e()
    >> ao.add_MO_block_2e(AOMOrangei=irrep1, MOrangej=irrep1, MOrangek=irrep2, MOrangel=irrep2)
    >> Vee[numpy.ix_(irrep1, irrep1, irrep2, irrep2)] = ao.Vee(asMO=True)
