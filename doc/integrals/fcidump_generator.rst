FCIDUMP Generator
=================

The FCIDUMP is an ASCII file containing all one- and two-electron integrals in molecular orbital basis required to set up a CI Hamiltonian. The file format is adopted from `Molpro`_.

.. _Molpro: http://www.molpro.net/info/current/doc/manual/node475.html

Creating a FCIDUMP file
-----------------------

The file can be generated from a QCinfo object. Assuming we start from some molden file, run in a Python Shell:

.. code-block:: python

    >> from orbkit import read, libcint_interface
    >> qc = read.main_read('h2o.molden', all_mo=True)
    >> libcint_interface.generate_fcidump(qc)

This will create a file named ``FCIDUMP``. The filename can be controlled by passing the ``filename`` parameter, an empty string disables writing to file. The return value of the function is a ``FCIDUMP`` object (see below).

Furthermore you can pass a point group to exploit symmetry, provided the molecular orbitals in the QCinfo object include an IRREP label.

.. code-block:: python

    >> libcint_interface.generate_fcidump(qc, filename='FCIDUMP_H2O_C2v', sym='c2v')

The active space can be controlled by the optional ``occ`` and ``closed`` parameters, which work similar to the corresponding cards in Molpro: They accept a list of integers, one for each IRREP. A given integer :math:`n` then sets the first :math:`n` orbitals in this IRREP to be *occupied* or *core* respectively.

.. code-block:: python

    >> libcint_interface.generate_fcidump(qc, filename='FCIDUMP_H2O_C2v', sym='c2v', core=[1,0,0,0], occ=[2,1,1,0])


Accessing integrals from the FCIDUMP
------------------------------------

The ``FCIDUMP`` class may also be used to access the MO integrals. First, either calculate them or read them from some FCIDUMP file:

.. code-block:: python

    >> fcidump = libcint_interface.generate_fcidump(qc, filename='')
    >> fcidump = libcint_interface.load_fcidump('FCIDUMP')

The ``fcidump`` object then gives access to the different integrals. Integrals are internally stored in spatial orbital basis, and only unique integrals are kept. For example ``fcidump.get_H()`` returns a upper triangular matrix with the 1-electron integrals, while the lower triangular part contains only zeros. If the complete matrix is needed, pass ``full=True``. For unrestricted orbitals the parameter ``spin='alpha'`` (default) or ``spin='beta'`` controlls the required subset. Integrals in spin orbital basis are obtained by ``fcidump.get_h()``.

.. code-block:: python

    >> H = fcidump.get_H(full=True)                          # 1-electron integrals for spatial orbitals
                                                             # equals fcidump.get_H(spin='alpha') for unrestricted case
    >> H_beta = fcidump.get_H(full=True, spin='beta')        # integrals for beta orbitals for unrestricted case
    >> H_spin = fcidump.get_h(full=True)                     # integrals over spin orbitals

In a very similar way 2-electron integrals can be accessed by calling ``fcidump.get_G()`` and ``fcidump.get_g()``. The former one accepts for the optional ``spin`` parameter the values ``alpha`` (default), ``beta`` and ``alphabeta``.

.. code-block:: python

    >> G = fcidump.get_G(full=True)                          # 2-electron integrals for spatial orbitals
    >> G_beta = fcidump.get_G(full=True, spin='beta')        # integrals for beta orbitals for unrestricted case
    >> G_beta = fcidump.get_G(full=True, spin='alphabeta')   # integrals for beta orbitals for unrestricted case
    >> G_spin = fcidump.get_g(full=True)                     # integrals over spin orbitals

The integral indices are ordered according to the chemists' notation:

**For spatial orbitals**:

    :math:`(ij|kl) = \int\int \mathrm{d}r_1 \mathrm{d}r_2 \phi_i^*(r_1) \phi_j(r_1) r_{12}^{-1} \phi_k^*(r_2) \phi_l(r_2)`

**For spin orbitals**:

    :math:`[ij|kl] = \int\int \mathrm{d}x_1 \mathrm{d}x_2 \phi_i^*(x_1) \phi_j(x_1) r_{12}^{-1} \phi_k^*(x_2) \phi_l(x_2)`

Furthermore the FCIDUMP holds some system parameters:

.. code-block:: python

    >> fcidump.nuclear_repulsion    # Coulomb repulsion energy of the nuclei
    >> fcidump.Nelec                # Number of electrons
    >> fcidump.Norb                 # Number of orbitals
    >> fcidump.spin                 # total spin of the electrons (=2S)
    >> fcidump.OrbSym               # IRREP of each orbital as a list of integers

It is possible to change the order of the orbitals with ``fcidump.change_order(order)``, where ``order`` is a list of indices. The active space can be reduces with ``fcidump.reduce_active_space(core, occ)``, where the parameters ``core`` and ``occ`` work as explained above.

These changes can be save to a (new) file by ``fcidump.store(filename)``.
