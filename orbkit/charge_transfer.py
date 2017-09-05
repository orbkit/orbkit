'''
Module for the analysis of the charge transfer character of a molecular system.
'''

# Import general modules
import numpy

# Import orbkit modules
from detci.density_matrix import DM




def get_nto(qc,ci):
  ''' 
  Function to calculate natural transition orbitals for CIS-type wavefunctions 
  according to 
  '''
  # Creates matrix for occupied and virtual orbitals
  LUMO = qc.mo_spec.get_lumo()
  mo_coeffs = qc.mo_spec.get_coeff()
  occ_mo = mo_coeffs[:LUMO]
  virt_mo = mo_coeffs[LUMO:]
  
  # Initialize NTOs
  ntos = {}
  ntos['vec'] = numpy.zeros(((len(ci),) + mo_coeffs.shape))
  ntos['val'] = numpy.zeros((len(ci),len(qc.mo_spec)))
  
  # Calculate NTOs (assuming that ci[0] is the ground state)
  for i in range(2):#len(ci)):
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
    ntos['vec'][i,:LUMO] = occ_nto
    virt_nto = numpy.dot(virt_mo.T,v_vec)
    virt_nto = virt_nto.T
    virt_nto = virt_nto[::-1]
    ntos['vec'][i,LUMO:] = virt_nto
    
    #dmat = DM(ci[0],ci[i],qc)
    #rdm = (1/numpy.sqrt(2))*dmat.Tij[:LUMO,LUMO:]
    
    ## Eigenvalue equation
    #(u_vec, sqrtlmbd, v_vec) = numpy.linalg.svd(rdm)
    #lmbd = sqrtlmbd * sqrtlmbd
    
    ## Eigenvalues for occupied and virtual orbitals
    #ntos['val'][i,:LUMO] = -lmbd[::-1]
    #ntos['val'][i,LUMO:(LUMO+len(lmbd))] = lmbd
    
    ## Eigenvectors for occupied and virtual orbitals
    #occ_nto = (numpy.dot(occ_mo.T,u_vec)).T
    #ntos['vec'][i,:LUMO] = occ_nto[::-1]
    
    #virt_nto = numpy.dot(virt_mo.T,v_vec)
    #virt_nto = (virt_nto.T)
    #ntos['vec'][i,LUMO:] = virt_nto
  
  return ntos,virt_nto
