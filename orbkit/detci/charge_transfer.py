'''
Module for the analysis of the charge transfer character of a molecular system.
'''

# Import general modules
import numpy,h5py

# Import orbkit modules
from .density_matrix import DM
from orbkit.output import molden_writer
from orbkit import grid,core,output
from orbkit.display import display, init_display


def compute_nto(qc,ci=[],Tij=[],md_filename=None):
  ''' 
  Function to calculate natural transition orbitals (NTO) for CIS-type wavefunctions 
  according to R.L. Martin, J. Chem. Phys. 2003, 118(11), 4775.

  **Parameters:**
  ci : CIinfo list of ci instances 
       See :ref:`Central Variables` for details.
  qc : class QCinfo
       See :ref:`Central Variables` for details.
  Tij : Transition density matrix.
       
    
  **Returns:**
    qc_nto is a list of qc instances containing attributes for geo_spec, geo_info, 
    ao_spec, mo_spec for all natural transition orbitals:
       See :ref:`Central Variables` for details of the QC Class.
      
  '''
  # Determine length of states
  if ci != []:
    stlen = len(ci)
  elif Tij != []:
    stlen = len(Tij)
  else:
    raise IOError('Variables ci and Tij are both empty.\n')
  
  display('Calculate natural transition orbitals.')
  # Creates matrix for occupied and virtual orbitals
  LUMO = qc.mo_spec.get_lumo()
  mo_coeffs = qc.mo_spec.get_coeffs()
  occ_mo = mo_coeffs[:LUMO]
  virt_mo = mo_coeffs[LUMO:]
  
  #Creation of symmetry and spin labels for NTOs
  sym = ['%d.NTO_h' % (i) for i in range(LUMO,0,-1)]
  sym += ['%d.NTO_p' % (i-LUMO+1) for i in range(LUMO,len(qc.mo_spec))]
  
  #Initialize QC Class for NTOs
  qc_nto = []
  
  # Calculate NTOs (assuming that ci[0] is the ground state)
  for i in range(stlen):

    # Initialize NTOs
    ntos = {}
    ntos['vec'] = numpy.zeros((mo_coeffs.shape))
    ntos['val'] = numpy.zeros(len(qc.mo_spec))
    if Tij == []:
      dmat = DM(ci[0],ci[i],qc)
      rdm = (1/numpy.sqrt(2))*dmat.Tij[:LUMO,LUMO:]
    else:
      rdm = Tij[i]
      
    # Eigenvalue equation
    (u_vec, sqrtlmbd, v_vec) = numpy.linalg.svd(rdm)
    lmbd = sqrtlmbd * sqrtlmbd

    # Eigenvalues for occupied and virtual orbitals
    ntos['val'][:LUMO] = -lmbd[::-1]
    ntos['val'][LUMO:(LUMO+len(lmbd))] = lmbd
    
    # Eigenvectors for occupied and virtual orbitals
    occ_nto = (numpy.dot(occ_mo.T,u_vec)).T
    ntos['vec'][:LUMO] = occ_nto[::-1]
    
    virt_nto = numpy.dot(v_vec,virt_mo)
    ntos['vec'][LUMO:] = virt_nto
    
    # Write new molden-Files
    qc_nto.append(qc.copy())
    qc_nto[-1].mo_spec.set_coeffs(ntos['vec'])
    qc_nto[-1].mo_spec.set_occ(ntos['val'])
    qc_nto[-1].mo_spec.set_sym(numpy.array(sym))
    if md_filename:
      molden_writer(qc_nto[-1],filename='nto_%s_%s' % (md_filename,i))

  return qc_nto

def compute_p_h_rho(qc_nto,min_val=1e-6,numproc=1,slice_length=1e4,hdf5_fid=None):
  
  display('Calculate particle and hole densities.')
  # Initialize grid
  if not grid.is_initialized:
    grid.adjust_to_geo(qc_nto[0],extend=5.0,step=0.4)
    display('\nSetting up the grid...')
    grid.grid_init(is_vector=False)
    display(grid.get_grid())   # Display the grid
    slice_length = grid.N_[1]*grid.N_[2]/2

  rho_h = []
  rho_p = []
  for i in range(len(qc_nto)):
    
    # MO selection
    qc = qc_nto[i].copy()
    if numpy.sum(abs(qc.mo_spec.get_occ())) != 0.0:
      bmo = abs(qc.mo_spec.get_occ()) >= min_val
      qc.mo_spec = qc.mo_spec[bmo]
      
      # Get indices for HOMO and LUMO
      sym = qc.mo_spec.get_sym()
      LUMO = numpy.sum([1 if s[-1]=='h' else 0 for s in sym])
      HOMO = LUMO - 1
      
      # Calculate MOs
      display('\nCalculating molecular orbitals...\n')
      molist = core.rho_compute(qc,calc_mo=True,drv=None,numproc=numproc)

      # Calculate hole and particle density
      rh = numpy.zeros(molist.shape[1:])
      rp = numpy.zeros(molist.shape[1:])
      for j in range(LUMO):
        rh += abs(qc.mo_spec[j]['occ_num'])*molist[j]*molist[j]
      for j in range(LUMO,len(qc.mo_spec)):
        rp += abs(qc.mo_spec[j]['occ_num'])*molist[j]*molist[j]
      rho_p.append(rp)
      rho_h.append(rh)
  
  rho_h = numpy.array(rho_h)
  rho_p = numpy.array(rho_p)
  
  if hdf5_fid:
    display('Save densities in HDF5 file...\n')
    # Initialize HDF5-File
    fid = hdf5_fid if hdf5_fid.endswith('.h5') else '%s.h5' % hdf5_fid
    output.hdf5_write(fid,mode='w',gname='general_info',x=grid.x,y=grid.y,z=grid.z,N=grid.N_,
            geo_spec=qc.geo_spec,geo_info=qc.geo_info,
            grid_info=numpy.array(grid.is_vector,dtype=int))
    f = h5py.File(fid, 'a')
    f['rho_p'] = rho_p
    f['rho_h'] = rho_h
    f.close()

  return rho_p,rho_h