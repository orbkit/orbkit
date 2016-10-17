from . import cy_occ_check
from .. import omp_functions
import numpy

def slice_occ(ij):  
  identical_states = multici['cia'] == multici['cib']
  if any([multici['cia'].method == i for i in ['mcscf','detci','fci']]):
    zero,sing = cy_occ_check.mcscf_ab(ij[0],ij[1],
                                      multici['cia'].coeffs,multici['cia'].occ,
                                      multici['cib'].coeffs,multici['cib'].occ,
                                      multici['moocc'],
                                      identical_states)
  elif any([multici['cia'].method == i for i in ['cis','tddft']]):
    if identical_states:
      zero,sing = cy_occ_check.cis_aa(ij[0],ij[1],
                                      multici['cia'].coeffs,multici['cia'].occ,
                                      multici['moocc'])
    else:
      zero,sing = cy_occ_check.cis_ab(ij[0],ij[1],
                                      multici['cia'].coeffs,multici['cia'].occ,
                                      multici['cib'].coeffs,multici['cib'].occ,
                                      multici['moocc'])
  return zero,sing


def compare(cia,cib,moocc,numproc=1):
  
  global multici
  multici = {'cia': cia, 'cib': cib,'moocc': moocc}
  
  numproc = min(len(cia.coeffs),max(1,numproc))
  ij = numpy.linspace(0, len(cia.coeffs), num=numproc+1, endpoint=True,  dtype=numpy.intc)
  ij = zip(ij[:-1],ij[1:])
  
  return_value = omp_functions.run(slice_occ,x=ij,numproc=numproc,display=lambda x: None)
  
  zero = [[],[]] 
  sing = [[],[]]
  for z,s in return_value:
    zero[0].extend(z[0])
    zero[1].extend(z[1])
    sing[0].extend(s[0])
    sing[1].extend(s[1])
  return zero,sing
