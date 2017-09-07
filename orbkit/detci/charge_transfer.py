'''
Module for the analysis of the charge transfer character of a molecular system.
'''

# Import general modules
import numpy

# Import orbkit modules
from .density_matrix import DM
from orbkit.output import molden_writer


def get_nto(qc,ci,md_filename=None):
  ''' 
  Function to calculate natural transition orbitals (NTO) for CIS-type wavefunctions 
  according to R.L. Martin, J. Chem. Phys. 2003, 118(11), 4775.

  **Parameters:**
  ci : CIinfo list of ci instances 
       See :ref:`Central Variables` for details.
  qc : class QCinfo
       See :ref:`Central Variables` for details.
    
  **Returns:**
    qc_nto is a list of qc instances containing attributes for geo_spec, geo_info, 
    ao_spec, mo_spec for all natural transition orbitals:
       See :ref:`Central Variables` for details of the QC Class.
      
  '''
  # Creates matrix for occupied and virtual orbitals
  LUMO = qc.mo_spec.get_lumo()
  mo_coeffs = qc.mo_spec.get_coeff()
  occ_mo = mo_coeffs[:LUMO]
  virt_mo = mo_coeffs[LUMO:]
  
  #Creation of symmetry and spin labels for NTOs
  sym = numpy.array([['NTO_h']*len(qc.mo_spec)][0], dtype=str)
  sym[LUMO:] = 'NTO_p'
  
  #Initialize QC Class for NTOs
  qc_nto = []
  
  # Calculate NTOs (assuming that ci[0] is the ground state)
  for i in range(len(ci)):
    
    # Initialize NTOs
    ntos = {}
    ntos['vec'] = numpy.zeros((mo_coeffs.shape))
    ntos['val'] = numpy.zeros(len(qc.mo_spec))
    dmat = DM(ci[0],ci[i],qc)
    rdm = (1/numpy.sqrt(2))*dmat.Tij[:LUMO,LUMO:]
    # Calculate reduced density matrix
    T = rdm
    U = numpy.dot(T,T.T)
    V = numpy.dot(T.T,T)

    u_val, u_vec = numpy.linalg.eigh(U)
    v_val, v_vec = numpy.linalg.eigh(V)
    v_val = v_val[::-1]

    occ_nto = numpy.dot(occ_mo.transpose(),u_vec)
    occ_nto = (occ_nto).transpose()
    ntos['vec'][:LUMO] = occ_nto
    virt_nto = numpy.dot(virt_mo.T,v_vec)
    virt_nto = virt_nto.T
    virt_nto = virt_nto[::-1]
    ntos['vec'][LUMO:] = virt_nto
    
    dmat = DM(ci[0],ci[i],qc)
    rdm = (1/numpy.sqrt(2))*dmat.Tij[:LUMO,LUMO:]
    
    # Eigenvalue equation
    (u_vec, sqrtlmbd, v_vec) = numpy.linalg.svd(rdm)
    lmbd = sqrtlmbd * sqrtlmbd
    
    # Eigenvalues for occupied and virtual orbitals
    ntos['val'][:LUMO] = -lmbd[::-1]
    ntos['val'][LUMO:(LUMO+len(lmbd))] = lmbd
    
    # Eigenvectors for occupied and virtual orbitals
    occ_nto = (numpy.dot(occ_mo.T,u_vec)).T
    ntos['vec'][:LUMO] = occ_nto[::-1]
    
    virt_nto = numpy.dot(virt_mo.T,v_vec)
    virt_nto = (virt_nto.T)
    ntos['vec'][LUMO:] = virt_nto
    
    # Write new molden-Files
    qc_nto.append(qc.copy())
    qc_nto[-1].mo_spec.set_coeff(ntos['vec'])
    qc_nto[-1].mo_spec.set_occ(ntos['val'])
    qc_nto[-1].mo_spec.set_sym(sym)
    if md_filename:
      molden_writer(qc,filename='nto_%s_%s' % (md_filename,i))

  return qc_nto
