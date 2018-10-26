#Tools for detci readers
import numpy
from orbkit.orbitals import MOClass

def point_groups():
  point_groups = {'c1': 1, 'cs': 2, 'c2': 2, 'ci': 2, 
                  'c2v': 4, 'c2h': 4, 'd2': 4, 'd2h': 8}
  return point_groups

def multiplicity():
  multiplicity = ['Unknown','Singlet','Doublet','Triplet','Quartet','Quintet']
  return multiplicity

def molpro_mo_order_ci(occ_info,mo_spec,irreps=None,nIRREP=None,order_sym=False):
  '''Orders the molecular orbitals according to the occupation information 
  of a CI class. (For processing MCSCF MOLPRO output)  
  
  **Parameters:**
  
    occ_info : ci.info["occ_info"]
      Contains a dict with the number of orbitals in the different 
      IRREPS for "core", "closed", "active", and "external" orbitals.
    mo_spec : list of dictionaries 
      See :ref:`Central Variables` for details.
  
  **Returns:**
  
    closed,active,external : list of mo_spec
      Contains mo_spec assigned to the "closed", "active", and "external" 
      orbitals.
  '''
  if nIRREP is None: nIRREP = len(occ_info['closed']) 
  if irreps is None: 
    irreps = dict(zip(map(str,range(1,nIRREP+1)),range(nIRREP)))
  else:
    nIRREP = len(irreps)
    irreps = dict(zip(irreps,range(nIRREP)))
  
  mo_index = [[] for i in range(nIRREP)]
  for i,i_mo in enumerate(mo_spec):
    info = i_mo['sym'].split('.')
    mo_index[irreps[info[1]]].append((i,int(info[0])))
  
  mo_sym = [[] for i in range(nIRREP)]
  for i in range(nIRREP):
    if mo_index[i] == []:
      continue
    mo_index[i] = numpy.array(mo_index[i])
    j = numpy.argsort(mo_index[i],axis=0)[:,1]
    mo_index[i] = mo_index[i][j]
    mo_sym[i] = MOClass([mo_spec[j] for j,k in mo_index[i]])
    #mo_sym[i].update()
  if order_sym:
    return mo_sym
  
  closed = MOClass()
  active = MOClass()
  external = MOClass()
  for i in range(nIRREP):
    a = occ_info['core'][i] + occ_info['closed'][i]
    closed.extend(mo_sym[i][:a])
    b = a+occ_info['active'][i]
    active.extend(mo_sym[i][a:b])
    external.extend(mo_sym[i][b:])
  #closed.update()
  #active.update()
  #external.update()
  return closed,active,external

def orthonorm(ci,reorth=False,**kwargs):
  '''Orthonomalize CI coefficients after Grahm-Schmidt 
  
  **Parameters:**
  
    ci: 
    reorth : bool (optional)
      if True, ci eigenstates will be reorthonormalized.
      else ci eigenstates will be renormalized.
  
  **Returns:**
  
    ci: 
      Contains the orthonormalized CI coefficients for each Slater-determinant
  '''
  # Grahm-Schmidt orthonormalization
  for ii in range(len(ci)):
    scal = numpy.dot(ci[ii].coeffs,ci[ii].coeffs)
    ci[ii].coeffs /= numpy.sqrt(scal)
    if reorth:
      for jj in range(ii+1,len(ci)):
        scal = numpy.dot(ci[ii].coeffs,ci[jj].coeffs)
        ci[jj].coeffs -= ci[ii].coeffs*scal

  return ci

