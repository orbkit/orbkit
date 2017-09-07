Introduction
============

This chapter explains in general how to use detCI\@ORBKIT to process multi-determinantal wave functions
to compute expectation values of different one-electron operators of two Configuration Interaction (CI) wavefunctions 
:math:`\left|\Phi_a\right\rangle` and :math:`\left|\Phi_b\right\rangle`:

.. math::

  \left\langle \Phi_a \middle| \hat{F} \middle| \Phi_b \right\rangle 
  = \sum_{p,q} C_p^{a}C_q^{b} \left\langle \phi_p \middle| \hat{F} \middle| \phi_q\right\rangle 

where :math:`\left|\phi_p\right\rangle` are the Slater determinants and :math:`C_p^{a}` are the respective CI coefficients.
The Slater determinats themselves are built up from the molecular orbitals :math:`\left|\varphi_m\right\rangle`.

For specific examples, please refer to the ORBKIT `examples folder`__.
There, you can find two tutorials: one for a complex superposition states of :math:`{\rm H}_3^+` and 
one for a superposition state of :math:`{\rm LiH}`.
They describe how to reproduce the results shown in

__ https://github.com/orbkit/orbkit/tree/cython/examples/detci

  Vincent Pohl, Gunter Hermann, and Jean Christophe Tremblay,
  "An Open-Source Framework of Analyzing *N*-electron Dynamics: I. Multi-Determinantal Wave Functions", 
  `arXiv:1701.06885`__ (2017).

__ http://arxiv.org/abs/1701.06885

detCI\@ORBKIT currently supports the following output file formats:

+----------------+---------------------+------------------------------+----------------+------------------------------------------------------------+ 
| **QC Program** | **Level of Theory** | **Ground State QC-Data**     | **CI-Data**    | **Example File**                                           | 
+----------------+---------------------+------------------------------+----------------+------------------------------------------------------------+ 
| PSI4           | All CI calculations | :ref:`Molden File`           | Output file    | :download:`psi4_fci.in <downloads/psi4_fci.in>`            | 
+----------------+---------------------+------------------------------+----------------+------------------------------------------------------------+  
| GAMESS-US      | CIS                 | :ref:`GAMESS-US Output File` | Output file    |                                                            | 
+----------------+---------------------+------------------------------+----------------+------------------------------------------------------------+  
| Turbomole      | TD-DFT              | :ref:`AOMix File`            | "sing_a" file  |                                                            | 
+----------------+---------------------+------------------------------+----------------+------------------------------------------------------------+  
| MOLPRO         | MCSCF               | :ref:`Molden File`           | Output file    | :download:`molpro_casscf.inp <downloads/molpro_casscf.inp>`| 
+----------------+---------------------+------------------------------+----------------+------------------------------------------------------------+  

.. Note::
  
  The underlying ORBKIT modules are called using :ref:`Low-Level Interface`.

.. --------------------+
..  **itype**          |
.. --------------------+
..  **'psi4_detci'**   |
.. --------------------+
..  **'gamess_cis'**   |
.. --------------------+
..  **'tmol_tddft'**   |
.. --------------------+
..  **'molpro_mcscf'** |
.. --------------------+

How to Read the QC Output
=========================

After importing the ``orbkit.read`` module and the detCI\@ORBKIT module, i.e.::
  
  from orbkit import read, detci
  
we can read the ground state quantum chemistry data, i.e., 
the molecular geometry and the atomic and molecular orbital data::

  qc = read.main_read(fid_molden,all_mo=True) 

Here, the quantum chemistry output is parsed into an instance of the ``QCinfo`` class 
(cf. :ref:`Central Variables`). Please note that we have to read the occupied 
*and* the virtual molecular orbitals (``all_mo=True``).

In a similar manner, we can read the Configuration Interaction (CI) output::

  qc,ci = detci.ci_read.main_ci_read(qc,fid_psi4,itype='psi4_detci',threshold=0.0)

where ``threshold`` specifies a read threshold for the CI coefficients,
which can considerably reduce computional time.
The output variable is a list of instances of the ``CIinfo`` class, which contains

  - ``ci[a].info``: A dictionary with several information on the electronic state such as the state energy, name, multiplicity, etc. 
  - ``ci[a].coeffs``: An array with CI coefficients
  - ``ci[a].occ``: An array with the corresponding occupation patterns 

.. Attention::
 
  The function ``main_ci_read`` changes its first argument (the ``QCinfo`` instance).
  Within this class the molecular orbitals are reordered according to their symmetry.
  This is required for all subsequent calculations. 

How to Prepare All Subsequent Calculations
==========================================

The starting point for all grid-based detCI\@ORBKIT calculations
is the computation of the molecular orbitals :math:`\varphi_m({\bf r})` and the derivatives thereof 
(:math:`\vec{\nabla}\varphi_m({\bf r})` and :math:`\vec{\nabla}^2\varphi_m({\bf r})`)::

  from orbkit import grid,core
  # Set up the grid
  grid.adjust_to_geo(qc,extend=5.0,step=0.1)
  grid.grid_init()    
  print(grid.get_grid())

  # Compute the molecular orbitals and their derviatives
  molist = core.rho_compute(qc,
                            calc_mo=True,
                            slice_length=1e4,           # Length of grid slice
                            drv=[None,                  # No derivative
                                'x','y','z',            # First derviatives
                                'xx','yy','zz'],        # Second derivatives
                            numproc=4)                  # Number of subprocesses
  molistdrv = molist[1:4]                               # \vec{\nabla} of MOs
  molistdrv2 = molist[-3:]                              # \vec{\nabla}^2 of MOs
  molist = molist[0]                                    # MOs

Non grid-based calculations, i.e., the electron number, dipole moments in length and velocity gauge,
require several expectation values, i.e., :math:`\langle\varphi_m|\varphi_n\rangle`, 
:math:`\langle\varphi_m|{\bf r}|\varphi_n\rangle`, and :math:`\langle\varphi_m|\vec{\nabla}|\varphi_n\rangle` ::


  from orbkit.analytical_integrals import get_ao_overlap,get_mo_overlap_matrix
  from orbkit.analytical_integrals import get_ao_dipole_matrix
  aoom = get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec,
                        ao_spherical=qc.ao_spherical,
                        drv=[None,                      # No derivative
                            'x','y','z'])               # First derviatives
  dm_aoom = get_ao_dipole_matrix(qc,component=['x','y','z'])

  coeff = qc.mo_spec.get_coeff()
  moom = get_mo_overlap_matrix(coeff,coeff,aoom[0],
                               numproc=4)               # <m|n>
  omr = numpy.zeros((3,) + moom.shape)                  # <m|r|n>
  omv = numpy.zeros((3,) + moom.shape)                  # <m|\vec{\nabla}|n>
  for i in range(3):
    omr[i] = get_mo_overlap_matrix(coeff,coeff,dm_aoom[i],numproc=4)
    omv[i] = get_mo_overlap_matrix(coeff,coeff,aoom[i+1] ,numproc=4)


How to Compare CI States
========================

Consider two electronic states :math:`|\Phi_a\rangle` and 
:math:`|\Phi_b\rangle`.
According to the Slater-Condon rules, an expectation value of a one-electron operator 

.. math::

  \left\langle \Phi_a \middle| \hat{F} \middle| \Phi_b \right\rangle 
  = \sum_{p,q} C_p^{a}C_q^{b} \left\langle \phi_p \middle| \hat{F} \middle| \phi_q\right\rangle 

where :math:`\left|\phi_p\right\rangle` are the Slater determinants and :math:`C_p^{a}` are the respective CI coefficients,
is defined by the following contributions 

  1. Two identical Slater determinats (``zero``):
    
    .. math::

      \left\langle \phi_p \middle| \hat{F} \middle| \phi_p\right\rangle  
      = \sum_m \left\langle \varphi_m \middle| \hat{F} \middle| \varphi_m\right\rangle

  2. Two Slater determinants differing by a single orbital :math:`m\rightarrow n` (``sing``)
  
    .. math::

      \left\langle \phi_p \middle| \hat{F} \middle| \phi_q\right\rangle  
      = \left\langle \varphi_m \middle| \hat{F} \middle| \varphi_n\right\rangle

  3. Two Slater determinants differing by two or more orbitals
  
    .. math::

      \left\langle \phi_p \middle| \hat{F} \middle| \phi_q\right\rangle  = 0

Beforehand, these determinants have to brought into maximum coincidence.  

All these non-zero contributions can be identified by calling::

  zero,sing = detci.occ_check.compare(ci[a],ci[b],numproc=4)

.. _`detCI:Electron Density`:

Electron Density
================

The electron (transition) density is defined by

.. math::

  \rho_{ab}({\bf r}) &= \left\langle \Phi_a \middle| \hat{\rho}({\bf r}) \middle| \Phi_b \right\rangle  \\
                   &= \Phi_a({\bf r})\Phi_b({\bf r})

It can be computed with::

  rho_ab = detci.ci_core.rho(zero,sing,molist,slice_length=1e4,numproc=4)

The integral over :math:`{\bf r}`, which may be evaluated analytically by::

  elnum = detci.ci_core.enum(zero,sing,moom)

yields the number of electrons for identical states and zero otherwise.

.. _`detCI:Electronic Flux Density`:

Electronic Flux Density and Its Divergence
==========================================

The electronic flux density (or Current Density) is defined as

.. math::

  j_{ab}({\bf r}) &= \left\langle \Phi_a \middle| \hat{j}({\bf r}) \middle| \Phi_b \right\rangle  \\
                  &= -\frac{i \hbar}{2m_e} \left(\Phi_a({\bf r})\vec{\nabla}\Phi_b({\bf r})-\Phi_b({\bf r})\vec{\nabla}\Phi_a({\bf r})\right)

where :math:`\vec{\nabla}` symbolizes the spatial derviatives with respect to the electronic coordinates. 
For the real-valued electronic eigenstates this quantity is purely imaginary. 
Thus, detCI\@ORBKIT exclusively computes the non-vanishing imaginary part of
the electronic flux density :math:`\mbox{Im}\left[j_{ab}\right]({\bf r})`::

  j_ab = detci.ci_core.jab(zero,sing,molist,molistdrv,slice_length=1e4,numproc=4)

Please note that the electronic flux density vanishes for identical states.
The divergence of the electronic flux density :math:`\vec{\nabla}\cdot j_{ab}`, which can be computed simply by::

  nabla_j_ab = detci.ci_core.jab(zero,sing,molist,molistdrv,slice_length=1e4,numproc=4)

is directly related to the time-derivative of the electron density via the 
electonic continuity equation

.. math::

  \frac{\partial \rho_{ab}({\bf r})}{\partial t} = - \vec{\nabla} \cdot j_{ab}({\bf r})

This relation may be reformulated into a convergence test:

.. math::

  - \vec{\nabla} \cdot \mbox{Im}\left[j_{ab}\right]({\bf r}) = -\frac{(E_b-E_a)}{\hbar} \rho_{ab}

Please note that the missing cusp relation for Gaussian basis sets leads to a  poor convergence around the nuclei, 
as shown in `arXiv:1701.06885`__.

__ http://arxiv.org/abs/1701.06885

.. _`detCI:Electronic Dipole Moment`:

Electronic Dipole Moment
========================

The electronic dipole moment can be defined in length gauge

.. math::

  \left\langle\mu_r\right\rangle_{ab} &= \left\langle \Phi_a \middle| \hat{\mu_r}({\bf r}) \middle| \Phi_b \right\rangle  \\
                  &= -e\int {\rm d}{\bf r} \left({\bf r}\rho_a({\bf r})\right)

and in velocity gauge 

.. math::

  \left\langle\mu_v\right\rangle_{ab} &= \left\langle \Phi_a \middle| \hat{\mu_v}({\bf r}) \middle| \Phi_b \right\rangle  \\
                  &= -e \int {\rm d}{\bf r}\ j_{ab}({\bf r})

Both quantities can be directly related to each other

.. math::

  \left\langle(\mu_v)_r\right\rangle_{ab} &= \frac{i\hbar}{(E_b-E_a)} 
		      \left\langle \Phi_a \middle| \hat{\mu_v}({\bf r}) \middle| \Phi_b \right\rangle \\
		  &= -\frac{e\hbar}{(E_b-E_a)}  \int {\rm d}{\bf r}\ \mbox{Im}\left[j_{ab}\right]({\bf r}) 

and may be used as another convergence test to directly evaluate the quality of the electronic flux density.

All three quantities can be calculated via::

  mu_ab = detci.ci_core.mu(ci[a],ci[b],qc,zero,sing,omr,omv)

Here, ``mu_ab[0]`` is :math:`\left\langle\mu_r\right\rangle_{ab}`, ``mu_ab[1]`` corresponds to :math:`\left\langle\mu_v\right\rangle_{ab}`, and 
``mu_ab[2]`` is :math:`\left\langle(\mu_v)_r\right\rangle_{ab}`.

For identical states, the latter is not defined., and besides, 
the permanent dipole moment of the molecule, i.e., including the nuclear contributions, is computed in that case.
