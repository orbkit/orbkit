from . import cy_ci
from ..analytical_integrals import get_nuclear_dipole_moment
from ..display import display
from ..core import require
from .. import omp_functions
import numpy, scipy

'''Implementation of ci_core using the reduced density matrix'''

def slice_rho(ij):  
  '''Computes the electron (transition) density on a grid for a single slice.
  '''
  return cy_ci.get_rho(ij[0],ij[1],multici['Ai'],multici['Vi'],multici['molist'])

def rho(dm,molist,slice_length=1e4,numproc=1):
  '''Computes the electron (transition) density on a grid.
  
  **Parameters:**
  
    dm : DM instance.
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
  assert isinstance(dm, orbkit.detci.reduced_density.DM):
  Ai = dm.get_eigenvalues()
  Vi = dm.get_eigenvalues()

  global multici
  # Reshape the input array if required
  molist = require(molist,dtype='f')
  shape = molist.shape
  molist.shape = (shape[0],-1)
  
  # Set the global array
  multici = {'Ai': Ai, 'Vi': Vi, 'molist':molist}
  
  slice_length = min(molist.shape[1],slice_length)
  ij = numpy.arange(0,molist.shape[1]+1,abs(int(slice_length)),dtype=numpy.intc)
  ij = zip(ij[:-1],ij[1:])
  
  data = numpy.zeros(molist.shape[1])
  return_val = omp_functions.run(slice_rho,x=ij,numproc=min(len(ij),numproc),display=display)
  for k,(i,j) in enumerate(ij):
    data[i:j] = return_val[k]
  
  # Restore the original shape
  molist.shape = shape
  return data.reshape(shape[1:],order='C')

