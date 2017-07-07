# -*- coding: iso-8859-1 -*-
r'''
This file is part of orbkit. See the main program or documentation 
for information on the license.

It shows two examples of how to use orbkit for evaluating 
analytical integrals.

First Part:

  Compute the electron number analytically using the standard density formula.

Second Part:

  Compute the dipole moments analytically. This expectation value is more complex.
  As an explanation, consider the following expectation value between two primitive 
  $p_z$ orbitals:

  .. math::

    <\phi_k|z|\phi_l> = \int dxdydz (z - Z_k)^{k_z} e^{\alpha_k r_k^2} (z - Z_l)^{l_z+1} e^{\alpha_l r_l^2} 
                        + Z_l \int dxdydz (z - Z_k)^{k_z} e^{\alpha_k r_k^2} (z - Z_l)^{l_z} e^{\alpha_l r_l^2} 
                      = ao_part_1 + ao_part_2

with :math:`r_l = \sqrt{(x-X_l)^2 + (y-Y_l)^2 + (z - Z_l)^2)}`
'''
from orbkit.read import main_read
from orbkit.tools import *
from orbkit.analytical_integrals import get_ao_overlap,get_mo_overlap,print2D

in_fid = 'h2o.molden'
# Read the input file
qc = main_read(in_fid,itype='molden',all_mo=False)

# Compute atomic orbital overlap matrix 
ao_overlap_matrix = get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec)

# Compute the overlap of the molecular orbitals and weight it with the occupation number
electron_number = 0.
for i_mo in qc.mo_spec:
  electron_number += i_mo['occ_num'] * get_mo_overlap(i_mo['coeffs'],
                                                      i_mo['coeffs'],
                                                      ao_overlap_matrix)

print('The total number of electrons is %.8f' % electron_number)
# Compute the x-, y-, and z-component of the dipole moment
dipole_moment = []
for component in range(3):
  
  # Compute the first part of the expectation value:
  # Get the the exponents lx, ly, lz for the primitive Cartesian Gaussians of
  # the `Ket` basis set, and increase lz by one.
  lxlylz_b = qc.ao_spec.get_lxlylz()
  lxlylz_b[:,component] += 1

  ao_part_1 = get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec,
                             lxlylz_b=lxlylz_b)

  # Compute the second part of the expectation value:
  ao_part_2 = get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec)

  i = 0
  for sel_ao in range(len(qc.ao_spec)):
    l = l_deg(l=qc.ao_spec[sel_ao]['type'].lower(),
              cartesian_basis=(not qc.ao_spec.spherical))
    for ll in range(l):
      ao_part_2[:,i] *= qc.geo_spec[qc.ao_spec[sel_ao]['atom'],component]
      i += 1
  
  ao_dipole_matrix = (ao_part_1+ao_part_2)
  
  # Print the atomic orbital dipole matrix
  if 0:
    print2D(ao_dipole_matrix) 
  
  # Compute the electronic part of the dipole moment
  dm = 0.
  for i,i_mo in enumerate(qc.mo_spec):
    dm -= i_mo['occ_num'] * get_mo_overlap(i_mo['coeffs'],i_mo['coeffs'],ao_dipole_matrix) 
  
  # Add the nuclear part
  for i_nuc in range(len(qc.geo_spec)):
    dm += float(qc.geo_info[i_nuc,2])*qc.geo_spec[i_nuc,component]
  
  dipole_moment.append(dm)

print('The dipole moment is \mu = (%.8f, %.8f, %.8f) ea_0' % tuple(dipole_moment))
print('\tRead from the molden_file h2o.md: \mu = (0,0,0.8174022121802662) ea_0')
