# -*- coding: iso-8859-1 -*-
'''Module performing analytical integrals between atomic and molecular orbitals.

Code for the computation of the overlap between primitive atomic basis functions
adapted from 

  M. Hô, J. M. Hernandez-Perez: "Evaluation of Gaussian Molecular Integrals", DOI:10.3888/tmj.14-3
'''
'''
orbkit
Gunter Hermann, Vincent Pohl, Lukas Eugen Marsoner Steinkasserer, Axel Schild, and Jean Christophe Tremblay

Institut fuer Chemie und Biochemie, Freie Universitaet Berlin, 14195 Berlin, Germany

This file is part of orbkit.

orbkit is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or any later version.

orbkit is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with orbkit.  If not, see <http://www.gnu.org/licenses/>.
'''
import numpy
from multiprocessing import Pool

from . import cy_overlap
from .tools import *
from .orbitals import AOClass, MOClass

def get_ao_overlap(coord_a, coord_b, ao_spec, lxlylz_b=None,
                   drv=None):
  '''Computes the overlap matrix of a basis set, where `Bra` basis set
  corresponds to the geometry :literal:`coord_a` and `Ket` basis set corresponds 
  to the geometry :literal:`coord_b`.
  
  In order to enable the computation of analytical expectation values, 
  the exponents lx, ly, lz for the primitive Cartesian Gaussians of the `Ket`
  basis set can be set manually with :literal:`lxlylz_b`.
  Please note that for the normalization of the primitive Cartesian Gaussians 
  the exponents from :literal:`ao_spec` are used.
  
  **Parameters:**
  
  coord_a : geo_spec
    Specifies the geometry of the `Bra` basis set. 
    See :ref:`Central Variables` in the manual for details.  
  coord_b : geo_spec
    Specifies the geometry of the `Ket` basis set. 
    See :ref:`Central Variables` in the manual for details.  
  ao_spec : 
    See :ref:`Central Variables` in the manual for details.   
  lxlylz_b : numpy.ndarray, dtype=numpy.int64, shape = (NAO,3), optional
    Contains the expontents lx, ly, lz for the primitive Cartesian Gaussians of
    the `Ket` basis set. 
  
  **Returns:**
  
  ao_overlap_matrix : numpy.ndarray, shape = (NAO,NAO)
    Contains the overlap matrix.
  '''
  if not isinstance(ao_spec, AOClass):
    raise TypeError('ao_spec must be an instance of the AOClass')
    
  if isinstance(drv, list) or (isinstance(drv, str) and len(drv) > 1):
    aoom = []
    for ii_d in drv:
      aoom.append(get_ao_overlap(coord_a,coord_b,ao_spec,
                                 lxlylz_b=lxlylz_b,
                                 drv=ii_d))
    return aoom
  lxlylz_a = ao_spec.get_lxlylz()

  if lxlylz_b is None:
    lxlylz_b = numpy.array(lxlylz_a,copy=True)
  else:
    try:
      lxlylz_b = numpy.array(lxlylz_b,dtype=numpy.intc)
    except ValueError:
      raise ValueError('The keyword argument `lxlylz` has to be convertable ' + 
                       'into a numpy integer array.')
    if lxlylz_a.shape != lxlylz_b.shape:
      raise ValueError('The exponents lxlylz for basis set a and basis set b ' +
                      'have to have the same shape.')

  # Derivative Calculation requested?  
  drv = validate_drv(drv)
  if drv > 3:
    raise ValueError('Only first derivatives are currently supported for ' +
                     'analytical integrals.')
  
  ao_overlap_matrix = cy_overlap.aooverlap(coord_a,
                                           coord_b,
                                           lxlylz_a,
                                           lxlylz_b,
                                           ao_spec.get_nlxlylz_per_cont(),
                                           ao_spec.get_prim_coeffs(),
                                           ao_spec.get_nprim_per_cont(),
                                           ao_spec.get_assign_cont_to_atoms(),
                                           drv,
                                           ao_spec.get_normalized())

  if 'N' in ao_spec[0]:
    for i in range(len(ao_overlap_matrix)):
      ao_overlap_matrix[i,:] *= ao_spec[0]['N'][i]*ao_spec[0]['N'][:,0]

  if ao_spec.spherical:
    # Convert the overlap matrix to the real-valued spherical harmonic basis.
    ao_overlap_matrix = cartesian2spherical_aoom(ao_overlap_matrix,ao_spec)
  return ao_overlap_matrix

def cartesian2spherical_aoom(ao_overlap_matrix,ao_spec):
  '''Transforms the atomic orbitals overlap matrix from a Cartesian Gaussian 
  basis set to a (real) pure spherical harmonic Gaussian basis set.
  
  Adapted from H.B. Schlegel and M.J. Frisch,
  International Journal of Quantum Chemistry, Vol. 54, 83-87 (1995).
  
  **Parameters:**
  
  ao_overlap_matrix : numpy.ndarray, shape = (NAO,NAO)
    Contains the overlap matrix of the Cartesian basis set. 
  ao_spec,ao_spherical :
    See :ref:`Central Variables` in the manual for details.
  
  **Returns:**
  
  aoom_sph : numpy.ndarray, shape = (NAO,NAO)
    Contains the overlap matrix of the spherical harmonic Gaussian basis set.
  
  ..hint: 
  
    Only supported up to g atomic orbitals and only for contracted 
    atomic orbitals.
  '''

  # Get the exponents of the Cartesian basis functions
  lxlylz = ao_spec.get_lxlylz()
  assign = ao_spec.get_assign_lxlylz_to_cont()
  ao_spherical  = ao_spec.get_old_ao_spherical()
  if ao_overlap_matrix.shape != (len(lxlylz),len(lxlylz)):
    raise IOError('No contraction is currently not supported for a '+ 
                  'spherical harmonics. Please come back'+
                  ' manually after calling `contract_ao_overlap_matrix()`.')

  l = [[] for i in ao_spec]
  for i,j in enumerate(assign):
    l[j].append(i) 

  indices = []
  c = 0
  for i0,(j0,k0) in enumerate(ao_spherical):
    sph0 = get_cart2sph(*k0)    
    for c0 in range(len(sph0[0])):
      for i,j in enumerate(l[j0]):
        if tuple(lxlylz[j]) == sph0[0][c0]:
          indices.append(i + l[j0][0])
      c += 1
  if len(indices) == 0:
    print('Here\n',ao_spec.up2date)
    exit()
  c = 0
  aoom_sph = numpy.zeros((len(ao_spherical),len(ao_spherical)))
  for i0,(j0,k0) in enumerate(ao_spherical):
    sph0 = get_cart2sph(*k0)    
    for c0 in range(len(sph0[0])):
      d = 0 
      for i1,(j1,k1) in enumerate(ao_spherical):
        sph1 = get_cart2sph(*k1)
        for c1 in range(len(sph1[0])):
          aoom_sph[i0,i1] += (sph0[1][c0]*sph0[2] * sph1[1][c1]*sph1[2]
                              * ao_overlap_matrix[indices[c],indices[d]])
          d += 1
      c+=1

  return aoom_sph

def get_mo_overlap(mo_a,mo_b,ao_overlap_matrix):
  '''Computes the overlap of two molecular orbitals.
  
  **Parameters:**
  
  mo_a : numpy.ndarray with shape = (,NAO)
     Contains the molecular orbital coefficients of all `Bra` orbitals.
  mo_b : numpy.ndarray with shape = (,NAO)
     Contains the molecular orbital coefficients of all `Ket` orbitals.
  ao_overlap_matrix : numpy.ndarray, shape = (NAO,NAO)
    Contains the overlap matrix of the basis set.
  
  **Returns:**
  
  mo_overlap : float
    Contains the overlap of the two input molecular orbitals.
  '''
  shape = numpy.shape(ao_overlap_matrix)
  if isinstance(mo_a,dict):
    mo_a = numpy.array(mo_a['coeffs'])
  if mo_a.ndim != 1 or len(mo_a) != shape[0]:
    raise ValueError('The coefficients of mo_a have to be a vector of the ' + 
                     'length of the ao_overlap_matrix.')
  if isinstance(mo_b,dict):
    mo_b = numpy.array(mo_b['coeffs'])
  if mo_b.ndim != 1 or len(mo_b) != shape[1]:
    raise ValueError('The coefficients of mo_b have to be a vector of the ' + 
                     'length of the ao_overlap_matrix.')

  return cy_overlap.mooverlap(require(mo_a,dtype='f'),require(mo_b,dtype='f'),
                              ao_overlap_matrix)

def initializer(gargs):
  global global_args
  global_args = gargs

def get_slice(x):  
  return cy_overlap.mooverlapmatrix(global_args['mo_a'],global_args['mo_b'],
                             global_args['ao_overlap_matrix'],x[0],x[1])

def get_mo_overlap_matrix(mo_a,mo_b,ao_overlap_matrix,numproc=1):
  '''Computes the overlap of two sets of molecular orbitals.
  
  **Parameters:**
  
  mo_a : numpy.ndarray with shape = (NMO,NAO), dict, or MOClass instance
     Contains the molecular orbital coefficients of all `Bra` orbitals.
  mo_b : numpy.ndarray with shape = (NMO,NAO), dict, or MOClass instance
     Contains the molecular orbital coefficients of all `Ket` orbitals.
  ao_overlap_matrix : numpy.ndarray, shape = (NAO,NAO)
    Contains the overlap matrix of the basis set.
  numproc : int
    Specifies number of subprocesses for multiprocessing.
  
  **Returns:**
  
  mo_overlap_matrix : numpy.ndarray, shape = (NMO,NMO)
    Contains the overlap matrix between the two sets of input molecular orbitals.
  '''
  if isinstance(mo_a, MOClass):
    mo_a = mo_a.get_coeffs()
  elif isinstance(mo_a,dict):
    mo_a = numpy.array(mo_a['coeffs'])

  if isinstance(mo_b, MOClass):
    mo_b = mo_b.get_coeffs()
  elif isinstance(mo_b,dict):
    mo_b = numpy.array(mo_b['coeffs'])

  global_args = {'mo_a': mo_a,
                 'mo_b': mo_b,
                 'ao_overlap_matrix': ao_overlap_matrix}

  if ((global_args['mo_a'].shape[1] != ao_overlap_matrix.shape[0]) or
      (global_args['mo_b'].shape[1] != ao_overlap_matrix.shape[1])):
    raise ValueError('mo_a and mo_b have to correspond to the same basis set, '+
                     'i.e., shape_a[1] != shape_b[1]')
  
  numproc = min(len(global_args['mo_a']),max(1,numproc))
  ij = numpy.array(numpy.linspace(0, len(global_args['mo_a']), 
                      num=numproc+1, endpoint=True),  dtype=numpy.intc)
  ij = list(zip(ij[:-1],ij[1:]))

  # Start the worker processes
  if numproc > 1:
    pool = Pool(processes=numproc, initializer=initializer, initargs=(global_args,))
    it = pool.imap(get_slice, ij)
  else:
    initializer(global_args)
  
  mo_overlap_matrix = numpy.zeros((len(mo_a),len(mo_b)),dtype=numpy.float64)
  
  #--- Send each task to single processor
  for l,[m,n] in enumerate(ij):
    #--- Call function to compute one-electron density
    mo_overlap_matrix[m:n,:] = it.next() if numproc > 1 else get_slice(ij[l])
  
  #--- Close the worker processes
  if numproc > 1:  
    pool.close()
    pool.join()
  
  #cy_overlap.mooverlapmatrix(moom,mo_a,mo_b,ao_overlap_matrix,0,len(moom))
  return mo_overlap_matrix

def get_moom_atoms(atoms,qc,mo_a,mo_b,ao_overlap_matrix,numproc=1):
  '''Computes the molecular orbital overlap matrix for selected atoms.
    
  **Parameters:**
  
  atoms : int or list of int
    Contains the indices of the selected atoms.
  mo_a : numpy.ndarray with shape = (NMO,NAO) or mo_spec (cf. :ref:`Central Variables`)
     Contains the molecular orbital coefficients of all `Bra` orbitals.  
  mo_b : numpy.ndarray with shape = (NMO,NAO) or mo_spec (cf. :ref:`Central Variables`)
     Contains the molecular orbital coefficients of all `Ket` orbitals.  
  ao_overlap_matrix : numpy.ndarray, shape = (NAO,NAO)
    Contains the overlap matrix of the basis set.
  numproc : int
    Specifies number of subprocesses for multiprocessing.
  
  **Returns:**
  
  mo_overlap_matrix : numpy.ndarray, shape = (NMO,NMO)
    Contains the overlap matrix between the two sets of input molecular orbitals.
  '''
  mo_a = mo_a.get_coeffs()
  mo_b = mo_b.get_coeffs()
  indices = get_lc(atoms,get_atom2mo(qc),strict=True)
  ao_overlap_matrix = numpy.ascontiguousarray(ao_overlap_matrix[:,indices])
  return get_mo_overlap_matrix(numpy.ascontiguousarray(mo_a),
                               numpy.ascontiguousarray(mo_b[:,indices]),
                               ao_overlap_matrix,numproc=numproc)

def check_mo_norm(qc):
  geo_spec = qc.geo_spec
  ao_spec = qc.ao_spec

  aoom = get_ao_overlap(geo_spec, geo_spec, ao_spec)
  moom = get_mo_overlap_matrix(mo_spec, mo_spec, aoom, numproc=options.numproc)

  deviation = numpy.linalg.norm(moom - numpy.eye(len(moom)))
  return deviation

def get_dipole_moment(qc,component=['x','y','z']):
  '''Computes the dipole moment analytically.
  
  **Parameters:**
  
  qc : class
    QCinfo class. (See :ref:`Central Variables` for details.)
  component : string or list of strings, {'x','y', or 'z'}
    Specifies the compontent(s) of the dipole moment which shall be computed.
  
  **Returns:**
  
  dipole_moment : 1D numpy.array, shape[0]=len(component)
    Contains the dipole moment.
  '''

  try:
    component = list(component)
  except TypeError: 
    component = [component]
  
  dipole_moment = numpy.zeros((len(component),))
  coeffs = qc.mo_spec.get_coeffs()
  occ = qc.mo_spec.get_occ()
  for i_d,c in enumerate(component):
    ao_dipole_matrix = get_ao_dipole_matrix(qc,component=c)
    for i_mo in range(len(qc.mo_spec)):
      dipole_moment[i_d] -= occ[i_mo] * get_mo_overlap(coeffs[i_mo,:],
                                                       coeffs[i_mo,:],
                                                       ao_dipole_matrix)

    # Add the nuclear part
    dipole_moment[i_d] += get_nuclear_dipole_moment(qc,component=c)
  return dipole_moment

def get_ao_dipole_matrix(qc,component='x'):
  '''Computes the expectation value of the dipole moment operator between 
  all atomic orbitals.
  
  **Parameters:**
  
  qc : class
    QCinfo class. (See :ref:`Central Variables` for details.)
  component : int or string, {'x','y', 'z'}
    Specifies the compontent of the dipole moment operator which shall be applied.
  
  **Returns:**
  
  ao_dipole_matrix : numpy.ndarray, shape=(NAO,NAO)
    Contains the expectation value matrix.
  '''
  
  if isinstance(component, list):
    aoom = []
    for ii_d in component:
      aoom.append(get_ao_dipole_matrix(qc,component=ii_d))
    return aoom
  if not isinstance(component, int):
    component = 'xyz'.find(component)
  if component == -1: # Was the selection valid?
    raise ValueError("The selection of the component was not valid!" +
              " (component = 'x' or 'y' or 'z')")

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
    if 'lxlylz' in qc.ao_spec[sel_ao].keys():
      l = len(qc.ao_spec[sel_ao]['lxlylz'])
    else:
      l = l_deg(l=qc.ao_spec[sel_ao]['type'].lower(),
              cartesian_basis=(not qc.ao_spec.spherical))
    for ll in range(l):
      ao_part_2[:,i] *= qc.geo_spec[qc.ao_spec[sel_ao]['atom'],component]
      i += 1

  # the atomic orbital overlap matrix  
  return (ao_part_1+ao_part_2) 

def get_nuclear_dipole_moment(qc,component='x'):
  '''Computes the nuclear part of the dipole moment.
  
  **Parameters:**
  
  qc : class
    QCinfo class. (See :ref:`Central Variables` for details.)
  component : string, {'x','y', 'z'}
    Specifies the compontent of the dipole moment operator which shall be applied.
  
  **Returns:**
  
  nuclear_dipole_moment : float
    Contains the nuclear dipole moment.
  '''
  if not isinstance(component, int):
    component = 'xyz'.find(component)
  if component == -1: # Was the selection valid?
    raise ValueError("The selection of the component was not valid!" +
              " (component = 'x' or 'y' or 'z')")
  
  nuclear_dipole_moment = 0.
  for i_nuc in range(len(qc.geo_spec)):
    nuclear_dipole_moment += float(qc.geo_info[i_nuc,2])*qc.geo_spec[i_nuc,component]
  return nuclear_dipole_moment
  
def get_atom2mo(qc):
  '''Assigns atom indices to molecular orbital coefficients.
  
  **Parameters:**
  
  qc.ao_spec :
      See :ref:`Central Variables` for details.
  
  **Returns:**
  
  atom2mo : numpy.ndarray, shape = (NAO,)
    Contains indices of atoms assigned to the molecular orbital coefficients.
  '''
  atom2mo = []
  a2mo_type = []
  b = 0  
  for sel_ao in range(len(qc.ao_spec)):
    a = qc.ao_spec[sel_ao]['atom']
    if 'lxlylz' in qc.ao_spec[sel_ao].keys():
      l = len(qc.ao_spec[sel_ao]['lxlylz'])
    else:
      l = l_deg(l=qc.ao_spec[sel_ao]['type'].lower(),
              cartesian_basis=(not qc.ao_spec.spherical))
    for i in range(l):
      atom2mo.append(a)
      b += 1
  
  return numpy.array(atom2mo,dtype=int)
  
def get_lc(atoms,atom2mo,strict=False):
  '''Returns indices of molecular orbital coefficients corresponding 
  to the selected atoms.
  
  **Parameters:**
  
  atoms : int or list of int
    Contains the indices of the selected atoms.
  atom2mo : numpy.ndarray, shape = (NAO,)
    Contains indices of atoms assigned to the molecular orbital coefficients.
    >> atom2mo = get_atom2mo(qc)
  
  **Returns:**
  
  lc : numpy.ndarray, shape = (NAO_atom,)
    Contains the NAO_atom indices molecular orbital coefficients corresponding 
    to the selected atoms.
  
  **Example:**
  
    >> atom2mo = get_atom2mo(qc)
    >> get_lc([0,3],atom2mo)
  '''
  if isinstance(atoms,int):
    atoms = [atoms]
  if not strict:
    lc = numpy.zeros(len(atom2mo),dtype=bool)
    for i in atoms:
      lc = numpy.logical_or(atom2mo==i,lc)

    return numpy.nonzero(lc)[0]
  else:
    lc = []
    for i in atoms:
      for k,j in enumerate(atom2mo):
        if i==j: lc.append(k)
    return numpy.array(lc,dtype=numpy.intc)

