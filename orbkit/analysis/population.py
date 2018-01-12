'''Module for atomic population analysis.
'''
import itertools
import numpy
from orbkit.analytical_integrals import get_atom2mo,get_lc
from orbkit.analytical_integrals import get_ao_overlap

def calc_P(qc):
  '''Calculates density matrix.'''
  C = qc.mo_spec.get_occ()[None,:]*qc.mo_spec.get_coeffs().T
  P = numpy.dot(C, qc.mo_spec.get_coeffs())
  return P

def mulliken(qc):
  '''Calculates the Mulliken populations (gross atomic populations) and Mulliken
  charges for each atom in the system.

  **Parameters:**

    qc.geo_spec, qc.geo_info, qc.ao_spec, qc.mo_spec :
      See :ref:`Central Variables` for details.

  **Returns:**
    mulliken_pop : dict
      Contains information of Mulliken charge analysis and has following members:
        :population: - Mulliken population for each atom.
        :charges: - Mulliken charges for each atom.

  '''
  ao_dim = qc.mo_spec.get_coeffs().shape[1]

  # Calculate AO-overlap matrix
  S = get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec)

  # Calculate density matrix
  P = calc_P(qc)

  # Calculate Mulliken population
  N = 0
  for m in itertools.product(range(ao_dim)):
    N += (P[:,m] * S[m,:])

  GP = N.diagonal()
  atom2mo = get_atom2mo(qc)
  GP_A = numpy.zeros(len(qc.geo_spec))
  for a in range(len(qc.geo_spec)):
    GP_A[a] = GP[get_lc(a,atom2mo)].sum()

  # Save Mulliken population and charges to dictionary
  mulliken_pop = {'population': GP_A,
                  'charge': numpy.array(qc.geo_info[:,2],dtype=float)-GP_A}

  return mulliken_pop

def lowdin(qc):
  '''Calculates the Lowdin populations and charges for each atom
      in the system.

  **Parameters:**

    qc.geo_spec, qc.geo_info, qc.ao_spec, qc.mo_spec :
      See :ref:`Central Variables` for details.

  **Returns:**
    lowdin_pop : dict
      Contains information of Lowdin charge analysis and has following members:
        :population: - Lowdin population for each atom.
        :charges: - Lowdin charges for each atom.

  '''
  ao_dim = qc.mo_spec.get_coeffs().shape[1]

  # Calculate AO-overlap matrix
  S = get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec)

  # Orthogonalize basis set
  s_eval, S_evec = numpy.linalg.eigh(S)
  S_eval = numpy.eye(len(s_eval))*s_eval**(0.5)
  S12 = numpy.dot(S_evec,numpy.dot(S_eval,S_evec.T))

  # Calculate density matrix
  P = calc_P(qc)

  # Calculate Lowdin population
  PS = 0
  for m in itertools.product(range(ao_dim)):
    PS += (P[:,m] * S12[m,:])

  N = 0
  for m in itertools.product(range(ao_dim)):
    N += (S12[:,m] * PS[m,:])

  GP = N.diagonal()
  atom2mo = get_atom2mo(qc)

  GP_A = numpy.zeros(len(qc.geo_spec))
  for a in range(len(qc.geo_spec)):
    GP_A[a] = GP[get_lc(a,atom2mo)].sum()

  lowdin_pop = {'population': GP_A,
              'charge': numpy.array(qc.geo_info[:,2],dtype=float)-GP_A}

  return lowdin_pop
