from . import cy_occ_check
from .. import omp_functions
from ..display import display
import numpy

def slice_occ(ij):
  '''Compares a slice of occupation patterns.
  '''
  if any([multici['cia'].method == i for i in ['mcscf','detci','fci','ci']]):
    ab_sorting = multici['cia'].method == 'mcscf' # Alternating alpha beta orbitals
    zero,sing = cy_occ_check.mcscf_ab(ij[0],ij[1],
                                      multici['cia'].coeffs,multici['cia'].occ,
                                      multici['cib'].coeffs,multici['cib'].occ,
                                      multici['moocc'],
                                      ab_sorting)
  elif any([multici['cia'].method == i for i in ['cis','tddft']]):
    zero,sing = cy_occ_check.cis_ab(ij[0],ij[1],
                                    multici['cia'].coeffs,multici['cia'].occ,
                                    multici['cib'].coeffs,multici['cib'].occ,
                                    multici['moocc'])
  return zero,sing


def compare(cia,cib,moocc=None,numproc=1):
  '''Compares occupation patterns of two CI vectors, and extracts identical 
  determinants and formal single excitations. 
  
  This is prerequisite for all subsequent detCI@orbkit calculations.
  
  **Parameters:**
  
    cia : CIinfo class instance 
      See :ref:`Central Variables` for details.
    cib : CIinfo class instance 
      See :ref:`Central Variables` for details.
    moocc : None or numpy.array, dtype=numpy.intc, optional
      Specifies the closed molecular orbitals. If None, it has to be a member
      of `cia` AND `cib`.
    numproc : int
      Specifies number of subprocesses for multiprocessing.
      
  **Returns:**
  
    zero : list of two lists
      Identical Slater-determinants. 
      | First member: Prefactor (occupation * CI coefficient)
      | Second member: Indices of occupied orbitals
    sing : list of two lists
      Effective single excitation.
      |  First member: Product of CI coefficients
      |  Second member: Indices of the two molecular orbitals 
  '''
  global multici
  if moocc is None:
    moocc = cia.moocc
    assert moocc is not None, '`moocc` is not given'
    assert numpy.array_equal(moocc,cib.moocc), '`cia.moocc` has to be the same as `cib.moocc`'
  
  multici = {'cia': cia, 'cib': cib, 'moocc': moocc}
  
  numproc = min(len(cia.coeffs),max(1,numproc))
  ij = numpy.array(numpy.linspace(0, len(cia.coeffs), num=numproc+1, endpoint=True),  
                   dtype=numpy.intc) # backward compatibility
  ij = list(zip(ij[:-1],ij[1:]))
  
  display('\nComparing the occupation patterns \nof the determinants of the two states...')
  return_value = omp_functions.run(slice_occ,x=ij,numproc=numproc,display=display)

  zero = [[],[]] 
  sing = [[],[]]
  for z,s in return_value:
    zero[0].extend(z[0])
    zero[1].extend(z[1])
    sing[0].extend(s[0])
    sing[1].extend(s[1])
  return zero,sing
