# -*- coding: iso-8859-1 -*-
'''Module performing analytical integrals between atomic and molecular orbitals.

Code for the computation of the overlap between primitive atomic basis functions
adapted from 

  M. Hô, J. M. Hernandez-Perez: "Evaluation of Gaussian Molecular Integrals", DOI:10.3888/tmj.14-3
'''
'''
orbkit
Gunter Hermann, Vincent Pohl, and Axel Schild

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
# test how to import weave
try:
    from scipy import weave
except:    
    import weave

from orbkit import cSupportCode
from orbkit.core import exp,lquant,slicer,get_lxlylz,get_cart2sph,l_deg

def get_ao_overlap(coord_a,coord_b,ao_spec,lxlylz_b=None,contraction=True,
                   drv=None,ao_spherical=None):
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
  
  contraction : bool, optional
    If True, the basis set will be contracted after the computation of the 
    overlap matrix.
  
  ao_spherical : optional
    Specifies if the input is given in a spherical harmonic Gaussian basis set.
    See :ref:`Central Variables` in the manual for details.
  
  **Returns:**
  
  ao_overlap_matrix : numpy.ndarray, shape = (NAO,NAO)
    Contains the overlap matrix.
  '''
  lxlylz_a = get_lxlylz(ao_spec)

  if lxlylz_b is None:
    lxlylz_b =  numpy.array(lxlylz_a,copy=True)
  else:
    try:
      lxlylz_b = numpy.array(lxlylz_b,dtype=numpy.int64)
    except ValueError:
      raise ValueError('The keyword argument `lxlylz` has to be convertable ' + 
                       'into a numpy integer array.')
    if lxlylz_a.shape != lxlylz_b.shape:
      raise ValueError('The exponents lxlylz for basis set a and basis set b ' +
                      'have to have the same shape.')
  
  # Derivative Calculation requested? 
  if drv is not None:
    # With respect to which variable the derivative shall be computed?
    if not isinstance(drv, (int, long)):
      drv = 'xyz'.find(drv)
    if drv == -1: # Was the selection valid? If not drv='x'
      drv = 0
      display("The selection of the derivative variable was not valid!" +
                " (drv = 'x' or 'y' or 'z')")
      display("Calculating the derivative with respect to x...")

  ra = numpy.array([])
  ra.shape = (0,3)
  rb = numpy.array(ra, copy=True)

  coeff = numpy.array([])
  coeff.shape = (0,2)

  index = []
  b = 0
  for sel_ao in range(len(ao_spec)):
    if 'exp_list' in ao_spec[sel_ao].keys():
      l = ao_spec[sel_ao]['exp_list']
    else:
      l = exp[lquant[ao_spec[sel_ao]['type']]]
    
    ra = numpy.append(ra,coord_a[ao_spec[sel_ao]['atom']][numpy.newaxis,:],axis=0)
    rb = numpy.append(rb,coord_b[ao_spec[sel_ao]['atom']][numpy.newaxis,:],axis=0)
    
    for i in l:
      coeff = numpy.append(coeff,ao_spec[sel_ao]['coeffs'],axis=0)
      for j in ao_spec[sel_ao]['coeffs']:
        index.append([sel_ao,b])
      b += 1

    
  index = numpy.array(index,dtype=numpy.int64)
  if contraction:
    ao_overlap_matrix = numpy.zeros((len(lxlylz_a),len(lxlylz_b)))
  else:
    ao_overlap_matrix = numpy.zeros((len(index),len(index)))

  # Ask for respective the C++ code
  code = ao_overlap_code(is_drv=(drv is not None))
  # A list of Python variable names that should be transferred from
  # Python into the C/C++ code. 
  arg_names = ['ra', 'rb', 'coeff', 'index', 
               'lxlylz_a', 'lxlylz_b', 'ao_overlap_matrix','drv']
  support_code = cSupportCode.overlap + cSupportCode.norm
  
  weave.inline(code, arg_names = arg_names, 
               support_code = support_code,
               headers=["<vector>"],verbose = 1)
  
  if not (ao_spherical is None or ao_spherical == []):
    if not contraction:
      raise IOError('No contraction is currently not supported for a '+ 
                    'spherical harmonics. Please apply `cartesian2spherical_aoom()`'+
                    ' manually after the contraction.')
    # Convert the overlap matrix to the real-valued spherical harmonic basis.
    ao_overlap_matrix = cartesian2spherical_aoom(ao_overlap_matrix,ao_spec,ao_spherical)
  
  return ao_overlap_matrix

def contract_ao_overlap_matrix(ao_uncontracted,ao_spec):
  '''Converts an overlap matrix built up from primitive (uncontracted) Cartesian 
  Gaussians to an overlap matrix of built up from contracted Cartesian Gaussians.
  
  **Parameters:**
  
  ao_uncontracted : numpy.ndarray, shape = (NAO,NAO)
    Contains the uncontracted overlap matrix.
  
  ao_spec : 
     See :ref:`Central Variables` in the manual for details.
  
  **Returns:**
  
  ao_overlap_matrix : numpy.ndarray, shape = (NAO,NAO)
    Contains the contracted overlap matrix.
  '''
  lxlylz_a = get_lxlylz(ao_spec)
  index = []
  b = 0
  for sel_ao in range(len(ao_spec)):
    if 'exp_list' in ao_spec[sel_ao].keys():
      l = ao_spec[sel_ao]['exp_list']
    else:
      l = exp[lquant[ao_spec[sel_ao]['type']]]
    for i in l:
      for j in ao_spec[sel_ao]['coeffs']:
        index.append([sel_ao,b])
      b += 1
  
  index = numpy.array(index,dtype=numpy.int64)
  ao_overlap_matrix = numpy.zeros((len(lxlylz_a),len(lxlylz_a)))
  
  code = '''
  int i_l, j_l;
  
  for (int i=0; i<Nindex[0]; i++)
  {
    for (int j=0; j<Nindex[0]; j++)
    {
      i_l = INDEX2(i,1);
      j_l = INDEX2(j,1);
      
      AO_OVERLAP_MATRIX2(i_l,j_l) += AO_UNCONTRACTED2(i,j);
    }
  }
  '''
  arg_names = ['index','ao_uncontracted', 'ao_overlap_matrix']

  weave.inline(code, arg_names = arg_names,verbose = 1)
    
  return ao_overlap_matrix

def cartesian2spherical_aoom(ao_overlap_matrix,ao_spec,ao_spherical):
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
  exp_list,assign = get_lxlylz(ao_spec,get_assign=True)
  
  if ao_overlap_matrix.shape != (len(exp_list),len(exp_list)):
    raise IOError('No contraction is currently not supported for a '+ 
                  'spherical harmonics. Please come back'+
                  ' manually after calling `contract_ao_overlap_matrix()`.')
  
  l = [[] for i in ao_spec]
  for i,j in enumerate(assign):
    l[j].append(i) 
  
  aoom_sph = numpy.zeros((len(ao_spherical),len(ao_spherical)))
  for i0,(j0,k0) in enumerate(ao_spherical):
    sph0 = get_cart2sph(*k0)
    for i1,(j1,k1) in enumerate(ao_spherical):
      sph1 = get_cart2sph(*k1)
      for c0 in range(len(sph0[0])):
        for i,j in enumerate(l[j0]):
          if tuple(exp_list[j]) == sph0[0][c0]:
            index0 = i + l[j0][0]
        for c1 in range(len(sph1[0])):
          for i,j in enumerate(l[j1]):
            if tuple(exp_list[j]) == sph1[0][c1]:
              index1 = i + l[j1][0]
          
          aoom_sph[i0,i1] += (sph0[1][c0]*sph0[2]*sph1[1][c1]*sph1[2]
                              *ao_overlap_matrix[index0,index1])
  return aoom_sph
  
def get_mo_overlap(mo_a,mo_b,ao_overlap_matrix):
  '''Computes the overlap of two molecular orbitals.
  
  **Parameters:**
  
  mo_a : numpy.ndarray, shape = (NAO,)
     Contains the molecular orbital coefficients of the `Bra` orbital.
  
  mo_b : numpy.ndarray, shape = (NAO,)
     Contains the molecular orbital coefficients of the `Ket` orbital.
  
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
  
  code = '''
  double overlap = 0.;
  for (int k=0; k<Nao_overlap_matrix[0]; k++)
  {
    for (int l=0; l<Nao_overlap_matrix[1]; l++)
    {
      overlap += MO_A1(k) * MO_B1(l) * AO_OVERLAP_MATRIX2(k,l);
    }
  }
  return_val = overlap;
  '''
  arg_names = ['mo_a', 'mo_b', 'ao_overlap_matrix']
  
  return weave.inline(code, arg_names = arg_names,verbose = 1)
  
def run_slice(sector):
  mo_overlap_matrix = numpy.zeros((sector[1]-sector[0],multiov['shape'][1]))
  mo_a = multiov['mo_a']
  mo_b = multiov['mo_b']
  ao_overlap_matrix = multiov['ao_overlap_matrix']
  arg_names = ['mo_a', 'mo_b', 'ao_overlap_matrix', 'mo_overlap_matrix','sector']
  weave.inline(moom_code, arg_names = arg_names, 
          support_code = cSupportCode.overlap + cSupportCode.norm,verbose = 1)
  return mo_overlap_matrix

def get_mo_overlap_matrix(mo_a,mo_b,ao_overlap_matrix,numproc=1):
  '''Computes the overlap of two sets of molecular orbitals.
  
  **Parameters:**
  
  mo_a : numpy.ndarray with shape = (NMO,NAO) or mo_spec (cf. :ref:`Central Variables`)
     Contains the molecular orbital coefficients of all `Bra` orbitals.
  
  mo_b : numpy.ndarray with shape = (NMO,NAO) or mo_spec (cf. :ref:`Central Variables`)
     Contains the molecular orbital coefficients of all `Ket` orbitals.
  
  ao_overlap_matrix : numpy.ndarray, shape = (NAO,NAO)
    Contains the overlap matrix of the basis set.
  
  **Returns:**
  
  mo_overlap_matrix : numpy.ndarray, shape = (NMO,NMO)
    Contains the overlap matrix between the two sets of input molecular orbitals.
  '''
  global multiov
  multiov = {'ao_overlap_matrix': ao_overlap_matrix,
             'mo_a': create_mo_coeff(mo_a,name='mo_a'),
             'mo_b': create_mo_coeff(mo_b,name='mo_b'),
             }
  
  shape_a = numpy.shape(multiov['mo_a'])
  shape_b = numpy.shape(multiov['mo_b'])
  multiov['shape'] = (shape_a[0],shape_b[0])
  
  if shape_a[1] != shape_b[1]:
    raise ValueError('mo_a and mo_b have to correspond to the same basis set, '+
                     'i.e., shape_a[1] != shape_b[1]')
  
  mo_overlap_matrix = numpy.zeros(multiov['shape']) 
  
  xx = slicer(N=len(multiov['mo_a']),
              vector=round(len(multiov['mo_a'])/float(numproc)),
              numproc=numproc)

  #--- Start the worker processes
  if numproc > 1:
    pool = Pool(processes=numproc)
    it = pool.imap(run_slice, xx)
  
  #--- Send each task to single processor
  for l,[m,n] in enumerate(xx):
    #--- Call function to compute one-electron density
    mo_overlap_matrix[m:n,:] = it.next() if numproc > 1 else run_slice(xx[l])
    
  #--- Close the worker processes
  if numproc > 1:  
    pool.close()
    pool.join()
  
  del multiov
  
  return mo_overlap_matrix

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
  for i,c in enumerate(component):
    ao_dipole_matrix = get_ao_dipole_matrix(qc,component=c)
    for i_mo in qc.mo_spec:
      dipole_moment[i] -= i_mo['occ_num'] * get_mo_overlap(i_mo['coeffs'],
                                                           i_mo['coeffs'],
                                                           ao_dipole_matrix)
    
    # Add the nuclear part
    dipole_moment[i] += get_nuclear_dipole_moment(qc,component=c)
  
  return dipole_moment

def get_ao_dipole_matrix(qc,component='x'):
  '''Computes the expectation value of the dipole moment operator between 
  all atomic orbitals.
  
  **Parameters:**
  
  qc : class
    QCinfo class. (See :ref:`Central Variables` for details.)
  component : int or string, {'x','y', 'z', 0, 1, 2}
    Specifies the compontent of the dipole moment operator which shall be applied.
  
  **Returns:**
  
  ao_dipole_matrix : numpy.ndarray, shape=(NAO,NAO)
    Contains the expectation value matrix.
  '''
  if not isinstance(component, (int, long)):
    component = 'xyz'.find(component)
  if component == -1: # Was the selection valid?
    raise ValueError("The selection of the component was not valid!" +
              " (component = 'x' or 'y' or 'z')")

  # Get the the exponents lx, ly, lz for the primitive Cartesian Gaussians of
  # the `Ket` basis set, and increase lz by one.
  lxlylz_b = get_lxlylz(qc.ao_spec)
  lxlylz_b[:,component] += 1

  ao_part_1 = get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec,
                             lxlylz_b=lxlylz_b,ao_spherical=qc.ao_spherical)

  # Compute the second part of the expectation value:
  ao_part_2 = get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec,
                             ao_spherical=qc.ao_spherical) 

  i = 0
  for sel_ao in range(len(qc.ao_spec)):
    if 'exp_list' in qc.ao_spec[sel_ao].keys():
      l = len(qc.ao_spec[sel_ao]['exp_list'])
    else:
      l = l_deg(l=qc.ao_spec[sel_ao]['type'].lower(),
              cartesian_basis=(qc.ao_spherical is None or qc.ao_spherical == []))
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
  component : int or string, {'x','y', 'z', 0, 1, 2}
    Specifies the compontent of the dipole moment operator which shall be applied.
  
  **Returns:**
  
  nuclear_dipole_moment : float
    Contains the nuclear dipole moment.
  '''
  if not isinstance(component, (int, long)):
    component = 'xyz'.find(component)
  if component == -1: # Was the selection valid?
    raise ValueError("The selection of the component was not valid!" +
              " (component = 'x' or 'y' or 'z')")
  
  nuclear_dipole_moment = 0.
  for i_nuc in range(len(qc.geo_spec)):
    nuclear_dipole_moment += float(qc.geo_info[i_nuc,2])*qc.geo_spec[i_nuc,component]
  return nuclear_dipole_moment

def is_mo_spec(mo):
  '''Checks if :literal:`mo` is of :literal:`mo_spec` type. 
  (See :ref:`Central Variables` for details.)'''
  if not isinstance(mo,list):
    return False
  return_val = True
  for i in mo:
    try:
      return_val = return_val and 'coeffs' in i.keys()
    except:
      return_val = False
  
  return return_val

def create_mo_coeff(mo,name='mo'):
  '''Converts the input variable to an :literal:`mo_coeff` numpy.ndarray.
  
  **Parameters:**
  
  mo : list, numpy.ndarray, or mo_spec (cf. :ref:`Central Variables`)
    Contains the molecular orbital coefficients of all orbitals.
  name : string, optional
    Contains a string describing the input variable. 
  
  **Returns:**
  
  mo : numpy.ndarray, shape = (NMO,NAO)
    Contains the molecular orbital coefficients of all orbitals.
  '''
  if (not is_mo_spec(mo)):
    if (not isinstance(mo,(list,numpy.ndarray))):
      raise ValueError('%s has to be mo_spec or an numpy coefficient array.'%s)
  else:
    tmp = []
    for i in mo:
      tmp.append(i['coeffs'])
    mo = tmp
  mo = numpy.array(mo, dtype=numpy.float64)  
  if mo.ndim != 2:
    raise ValueError('%s has to be 2-dimensional.'%s)  
  return mo
  
def get_atom2mo(qc):
  '''Assigns atoms to molecular orbital coefficients.
  
  **Parameters:**
  
  qc.ao_spec :
      See :ref:`Central Variables` for details.
  
  **Returns:**
  
  atom2mo : numpy.ndarray, shape = (NAO)
    Contains indices of atoms assigned to the molecular orbital coefficients.
  '''
  atom2mo = []
  a2mo_type = []
  b = 0
  for sel_ao in range(len(qc.ao_spec)):
    a = qc.ao_spec[sel_ao]['atom']
    l = exp[lquant[qc.ao_spec[sel_ao]['type']]]
    for i in l:
      atom2mo.append(a)
      a2mo_type.append([a,qc.ao_spec[sel_ao]['type'],l])
      b += 1
  
  return numpy.array(atom2mo,dtype=int)
  
def get_lc(indices,atom2mo):
  '''Converts the input variable to an :literal:`mo_coeff` numpy.ndarray.
  
  **Parameters:**
  
  mo : list, numpy.ndarray, or mo_spec (cf. :ref:`Central Variables`)
    Contains the molecular orbital coefficients of all orbitals.
  name : string, optional
    Contains a string describing the input variable. 
  
  **Returns:**
  
  mo : numpy.ndarray, shape = (NMO,NAO)
    Contains the molecular orbital coefficients of all orbitals.
  '''
  if isinstance(indices,int):
    indices = [indices]
  lc = numpy.zeros(len(atom2mo),dtype=bool)
  for i in indices:
    lc = numpy.logical_or(atom2mo==i,lc)

  return numpy.nonzero(lc)[0]

def print2D(x,format='%+.2f ',start='\t'):
  '''Prints a 2D matrix.
  
  **Parameters:**
  
  x : numpy.ndarray, shape = (n,m)
    Contains a 2D matrix.
  
  format : str
    Specifies the output format.
  '''
  shape = numpy.shape(x)
  for i in range(shape[0]):
    s = start
    for j in range(shape[1]):
      s += format % x[i,j]
    print s

def ao_overlap_code(is_drv=False):
  '''Returns the requested C++ code.
  
  **Parameters:**
  
  is_drv : bool
    If True, returns the code for the computation of the
    derivative of an atomic orbital.
  '''
  
  if not is_drv:
    code = '''
    // Hint: ao_spec coeff -> [exponent, coefficient]
    S_Primitive A, B;
    int i_ao, j_ao;
    int i_l, j_l;
    int index_i,index_j;
    bool contraction = (Nindex[0] != Nao_overlap_matrix[0]);
    
    std::vector<double> norm(Nindex[0]);
    
    for (int i=0; i<Nindex[0]; i++)
    {
      i_l = INDEX2(i,1);
      norm.at(i) = ao_norm(LXLYLZ_A2(i_l,0), LXLYLZ_A2(i_l,1), LXLYLZ_A2(i_l,2),&COEFF2(i,0));
    }
    
    
    for (int i=0; i<Nindex[0]; i++)
    {
      for (int j=0; j<Nindex[0]; j++)
      {
        i_ao = INDEX2(i,0);
        j_ao = INDEX2(j,0);
        i_l = INDEX2(i,1);
        j_l = INDEX2(j,1);
        if (contraction)
        {
          index_i = INDEX2(i,1);
          index_j = INDEX2(j,1);
        }
        else
        {
          index_i = i;
          index_j = j;
        }
        
        for (int rr=0; rr<3; rr++)
        {
          A.R[rr] = RA2(i_ao,rr);
          B.R[rr] = RB2(j_ao,rr);
          A.l[rr] = LXLYLZ_A2(i_l,rr);
          B.l[rr] = LXLYLZ_B2(j_l,rr);
        }
        
        A.alpha = COEFF2(i,0);
        B.alpha = COEFF2(j,0);
        
        AO_OVERLAP_MATRIX2(index_i,index_j) += COEFF2(i,1) * COEFF2(j,1) *
            norm.at(i) * norm.at(j) * get_overlap(&A, &B);
      }
    }
    
    '''
  else:
    code = '''
    // Hint: ao_spec coeff -> [exponent, coefficient]
    S_Primitive A, B;
    int i_ao, j_ao;
    int i_l, j_l;
    int index_i,index_j;
    bool contraction = (Nindex[0] != Nao_overlap_matrix[0]);
    
    std::vector<double> norm(Nindex[0]);
    
    for (int i=0; i<Nindex[0]; i++)
    {
      i_l = INDEX2(i,1);
      norm.at(i) = ao_norm(LXLYLZ_A2(i_l,0), LXLYLZ_A2(i_l,1), LXLYLZ_A2(i_l,2),&COEFF2(i,0));
    }
    
    
    for (int i=0; i<Nindex[0]; i++)
    {
      for (int j=0; j<Nindex[0]; j++)
      {
        i_ao = INDEX2(i,0);
        j_ao = INDEX2(j,0);
        i_l = INDEX2(i,1);
        j_l = INDEX2(j,1);
        if (contraction)
        {
          index_i = INDEX2(i,1);
          index_j = INDEX2(j,1);
        }
        else
        {
          index_i = i;
          index_j = j;
        }
        
        for (int rr=0; rr<3; rr++)
        {
          A.R[rr] = RA2(i_ao,rr);
          B.R[rr] = RB2(j_ao,rr);
          A.l[rr] = LXLYLZ_A2(i_l,rr);
          B.l[rr] = LXLYLZ_B2(j_l,rr);
        }
        
        A.alpha = COEFF2(i,0);
        B.alpha = COEFF2(j,0);
        
        if (B.l[drv] == 0)
        {
          B.l[drv] = LXLYLZ_B2(j_l,drv)+1;
          AO_OVERLAP_MATRIX2(index_i,index_j) += (-2 * B.alpha) * COEFF2(i,1) * 
                      COEFF2(j,1) * norm.at(i) * norm.at(j) * get_overlap(&A, &B);
        }
        else
        {
          B.l[drv] = LXLYLZ_B2(j_l,drv)-1;
          AO_OVERLAP_MATRIX2(index_i,index_j) += LXLYLZ_B2(j_l,drv) * COEFF2(i,1) * 
                        COEFF2(j,1) * norm.at(i) * norm.at(j) * get_overlap(&A, &B);
          B.l[drv] = LXLYLZ_B2(j_l,drv)+1;
          AO_OVERLAP_MATRIX2(index_i,index_j) += (-2. * B.alpha) * COEFF2(i,1) * 
                      COEFF2(j,1) * norm.at(i) * norm.at(j) * get_overlap(&A, &B);
        }
      }
    }
    
    '''
  return code

moom_code = '''
int i0 = int(sector[0]);
for (int i=i0; i<int(sector[1]); i++)
{
  for (int j=0; j<Nmo_b[0]; j++)
  {
    for (int k=0; k<Nao_overlap_matrix[0]; k++)
    {
      for (int l=0; l<Nao_overlap_matrix[1]; l++)
      {
        MO_OVERLAP_MATRIX2(i-i0,j)  += MO_A2(i,k) * MO_B2(j,l) * AO_OVERLAP_MATRIX2(k,l);
      }
    }
  }
}
'''
