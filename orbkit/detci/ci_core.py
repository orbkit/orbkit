from . import cy_ci
from ..analytical_integrals import get_nuclear_dipole_moment
from ..display import display
from ..tools import require
from .. import omp_functions
import numpy

def mu(cia,cib,qc,zero,sing,omr,omv):
  '''Computes the analytic transition dipole moments in length and velocity gauge
  using analytical integrals.
  
  mu_r = <psi_a|-er|psi_b>
  mu_v = <psi_a|-e\hat{j}|psi_b> 
  mur_v = (i hbar)/(E_b-E_a) mu_v

  **Parameters:**
  
    cia : CIinfo class instance 
      See :ref:`Central Variables` for details.
    cib : CIinfo class instance 
      See :ref:`Central Variables` for details.
    qc : class QCinfo
      See :ref:`Central Variables` for details.
    zero : list of two lists
      Identical Slater-determinants. 
      | First member: Prefactor (occupation * CI coefficient)
      | Second member: Indices of occupied orbitals
    sing : list of two lists
      Effective single excitation.
      |  First member: Product of CI coefficients
      |  Second member: Indices of the two molecular orbitals 
    omr : numpy.ndarray, shape=(3,NMO,NMO)
      Contains the molecular orbital dipole expectation value matrix.
    omv : numpy.ndarray, shape=(3,NMO,NMO)
      Contains the molecular orbital derivative expectation value matrix.
  
  **Returns:**
  
    mur : numpy.ndarray, shape=(3,)
      Contains the dipole moment in length gauge.
    muv : numpy.ndarray, shape=(3,)
      Contains the imaginary part of the dipole moment in velocity gauge. (There is no real part!)
    mur_v : numpy.ndarray, shape=(3,)
      Contains the dipole moment in length gauge converted from muv.
  '''
  mur,muv = cy_ci.get_mu(zero,sing,omr,omv) 

  if cia == cib:
    for component in range(3):
      mur[component] = mur[component] + get_nuclear_dipole_moment(qc,component=component)
  
  delta_v = float(cib.info['energy'])-float(cia.info['energy'])
  mur_v = numpy.zeros((3,))
  if delta_v != 0.0:
    mur_v = -muv/delta_v
  
  return numpy.array([mur, muv, mur_v])


def enum(zero,sing,moom):
  '''Computes the norm of the wavefunction using analytical integrals.
  
  **Parameters:**
  
    zero : list of two lists
      Identical Slater-determinants. 
      | First member: Prefactor (occupation * CI coefficient)
      | Second member: Indices of occupied orbitals
    sing : list of two lists
      Effective single excitation.
      |  First member: Product of CI coefficients
      |  Second member: Indices of the two molecular orbitals 
    moom : numpy.ndarray, shape=(NMO,NMO)
      Contains the molecular orbital overlap matrix.
  
  **Returns:**
  
    enum : float
      Contains the number of electrons.

  '''
  return cy_ci.get_enum(zero,sing,moom)


def slice_rho(ij):  
  '''Computes the electron (transition) density on a grid for a single slice.
  '''
  return cy_ci.get_rho(ij[0],ij[1],multici['zero'],multici['sing'],multici['molist'])

def rho(zero,sing,molist,slice_length=1e4,numproc=1):
  '''Computes the electron (transition) density on a grid.
  
  **Parameters:**
  
    zero : list of two lists
      Identical Slater-determinants. 
      | First member: Prefactor (occupation * CI coefficient)
      | Second member: Indices of occupied orbitals
    sing : list of two lists
      Effective single excitation.
      |  First member: Product of CI coefficients
      |  Second member: Indices of the two molecular orbitals 
    molist : numpy.ndarray, shape=(NMO,N)
      Contains the NMO=len(qc.mo_spec) molecular orbitals on a grid.
    slice_length : int, optional
      Specifies the number of points per subprocess.
    numproc : int, optional
      Specifies the number of subprocesses for multiprocessing.
  
  **Returns:**
  
    rho : numpy.ndarray, shape=(N,)
      Contains the density on a grid.

  '''
  global multici
  # Reshape the input array if required
  molist = require(molist,dtype='f')
  shape = molist.shape
  molist.shape = (shape[0],-1)
  
  # Set the global array
  multici = {'zero': zero, 'sing': sing, 'molist':molist}
  
  slice_length = min(molist.shape[1],slice_length)
  ij = numpy.arange(0,molist.shape[1]+1,abs(int(slice_length)),dtype=numpy.intc)
  ij = list(zip(ij[:-1],ij[1:]))
  
  data = numpy.zeros(molist.shape[1])
  return_val = omp_functions.run(slice_rho,x=ij,numproc=min(len(ij),numproc),display=display)
  for k,(i,j) in enumerate(ij):
    data[i:j] = return_val[k]
  
  # Restore the original shape
  molist.shape = shape
  return data.reshape(shape[1:],order='C')

def slice_jab(ij):   
  '''Computes the electronic (transition) flux density on a grid for a single slice.
  '''  
  return cy_ci.get_jab(ij[0],ij[1],multici['zero'],multici['sing'],multici['molist'],multici['molistdrv'])


def jab(zero,sing,molist,molistdrv,slice_length=1e4,numproc=1):
  r'''Computes the imaginary part of the electronic transition flux density on a grid.

  .. math::
  
    j_{a,b} = -(i hbar)/(2me) * (\psi_a \nabla \psi_b - \psi_b \nabla \psi_a)    
            = 0 + i Im( j_{a,b} )
  
  **Parameters:**
  
    zero : list of two lists
      Identical Slater-determinants. 
      | First member: Prefactor (occupation * CI coefficient)
      | Second member: Indices of occupied orbitals
    sing : list of two lists
      Effective single excitation.
      |  First member: Product of CI coefficients
      |  Second member: Indices of the two molecular orbitals 
    molist : numpy.ndarray, shape=(NMO,N)
      Contains the NMO=len(qc.mo_spec) molecular orbitals on a grid (N points).
    molistdrv : numpy.ndarray, shape=(3,NMO,N)
      Contains the three (x,y,z) derivatives of the molecular orbitals on a grid.
    slice_length : int, optional
      Specifies the number of points per subprocess.
    numproc : int, optional
      Specifies the number of subprocesses for multiprocessing.
  
  **Returns:**
  
    rho : numpy.ndarray, shape=(3,N)
      Contains the imaginary part of the electronic (transition) flux density on a grid.
  
  .. note:: 
    
    Since the electronic eigenfunctions are real-valued the corresponding 
    electronic transition flux density is purely imaginary. 
  '''  
  global multici
  # Reshape the input array if required
  molist = require(molist,dtype='f')
  molistdrv = require(molistdrv,dtype='f')
  shape = molist.shape
  molist.shape = (shape[0],-1)
  molistdrv.shape = (3,shape[0],-1)
  
  # Set the global array
  multici = {'zero': zero, 'sing': sing, 'molist':molist, 'molistdrv':molistdrv}
    
  slice_length = min(molist.shape[1],slice_length)
  ij = numpy.arange(0,molist.shape[1]+1,abs(int(slice_length)),dtype=numpy.intc)
  ij = list(zip(ij[:-1],ij[1:]))
  
  data = numpy.zeros((3,molist.shape[1]))
  return_val = omp_functions.run(slice_jab,x=ij,numproc=min(len(ij),numproc),display=display)
  for k,(i,j) in enumerate(ij):
    data[:,i:j] = return_val[k]

  # Restore the original shape
  molist.shape = shape
  molistdrv.shape = (3,) + shape
  return data.reshape((3,) + shape[1:],order='C')

def slice_a_nabla_b(ij): 
  '''Computes the \psi_a \nabla \psi_b on a grid for a single slice.
  '''    
  return cy_ci.get_a_nabla_b(ij[0],ij[1],multici['zero'],multici['sing'],multici['molist'],multici['molistdrv'])


def a_nabla_b(zero,sing,molist,molistdrv,slice_length=1e4,numproc=1):
  r'''Computes the following quantity (`a_nabla_b`) on a grid:

  .. math::
  
    \psi_a \nabla \psi_b
    
  **Parameters:**
  
    zero : list of two lists
      Identical Slater-determinants. 
      | First member: Prefactor (occupation * CI coefficient)
      | Second member: Indices of occupied orbitals
    sing : list of two lists
      Effective single excitation.
      |  First member: Product of CI coefficients
      |  Second member: Indices of the two molecular orbitals 
    molist : numpy.ndarray, shape=(NMO,N)
      Contains the NMO=len(qc.mo_spec) molecular orbitals on a grid (N points).
    molistdrv : numpy.ndarray, shape=(3,NMO,N)
      Contains the three (x,y,z) derivatives of the molecular orbitals on a grid.
    slice_length : int, optional
      Specifies the number of points per subprocess.
    numproc : int, optional
      Specifies the number of subprocesses for multiprocessing.
  
  **Returns:**
  
    a_nabla_b : numpy.ndarray, shape=(3,N)
      Contains the data on a grid.

  '''  
  global multici
  # Reshape the input array if required
  molist = require(molist,dtype='f')
  molistdrv = require(molistdrv,dtype='f')
  shape = molist.shape
  molist.shape = (shape[0],-1)
  molistdrv.shape = (3,shape[0],-1)
  
  # Set the global array
  multici = {'zero': zero, 'sing': sing, 'molist':molist, 'molistdrv':molistdrv}
    
  slice_length = min(molist.shape[1],slice_length)
  ij = numpy.arange(0,molist.shape[1]+1,abs(int(slice_length)),dtype=numpy.intc)
  ij = list(zip(ij[:-1],ij[1:]))
  
  data = numpy.zeros((3,molist.shape[1]))
  return_val = omp_functions.run(slice_a_nabla_b,x=ij,numproc=min(len(ij),numproc),display=display)
  for k,(i,j) in enumerate(ij):
    data[:,i:j] = return_val[k]

  # Restore the original shape
  molist.shape = shape
  molistdrv.shape = (3,) + shape
  return data.reshape((3,) + shape[1:],order='C')
