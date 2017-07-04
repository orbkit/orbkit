import numpy
import string

# Assign the quantum number l to every AO symbol (s,p,d,etc.) 
orbit = 'spd' + string.ascii_lowercase[5:].replace('s','').replace('p','')
lquant = dict([(j, i) for i,j in enumerate(orbit)])

def l_deg(l=0,ao=None,cartesian_basis=True):
  '''Calculates the degeneracy of a given atomic orbitals.
  
  **Options:**
  
  Works with the molpro output nomenclature for Cartesian Harmonics:
    s->'s', p->['x','y','z'], d-> ['xx','yy', etc.], etc.
    e.g., l_deg(ao='xxy')
  
  Works with quantum number l for the Cartesian Harmonic:
    e.g., l_deg(l=1)
  
  Works with name of the Cartesian Harmonic:
    e.g., l_deg(l='p')
  ''' 
  if ao != None:
    if ao == 's':
      return 1
    else:
      l = len(ao)
  elif isinstance(l,str):
    l = lquant[l]
  return int((l+1)*(l+2)/2) if cartesian_basis else int(2*l+1)
  # l_deg 

# Molden AO order 
exp = []
exp.append([(0,0,0)])                   # s orbitals

exp.append([(1,0,0), (0,1,0), (0,0,1)]) # p orbitals

exp.append([(2,0,0), (0,2,0), (0,0,2),
            (1,1,0), (1,0,1), (0,1,1)]) # d orbitals

exp.append([(3,0,0), (0,3,0), (0,0,3),
            (1,2,0), (2,1,0), (2,0,1),
            (1,0,2), (0,1,2), (0,2,1),
            (1,1,1)])                   # f orbitals
      
exp.append([(4,0,0), (0,4,0), (0,0,4),
            (3,1,0), (3,0,1), (1,3,0),
            (0,3,1), (1,0,3), (0,1,3),
            (2,2,0), (2,0,2), (0,2,2),
            (2,1,1), (1,2,1), (1,1,2)]) # g orbitals

# wfn order of exponents 
exp_wfn = exp[:3]                           # s,p,d orbitals 

exp_wfn.append([(3,0,0), (0,3,0), (0,0,3),                   
                (2,1,0), (2,0,1),(0,2,1),                    
                (1,2,0), (1,0,2), (0,1,2),                   
                (1,1,1)])                   # f orbitals     

exp_wfn.append(exp[4]) # g orbitals     


'''                                                                             
Transformation Between Cartesian and (Real) Pure Spherical Harmonic Gaussians   

adapted from H.B. Schlegel and M.J. Frisch 
International Journal of Quantum Chemistry, Vol. 54, 83-87 (1995).
'''
sqrt = numpy.sqrt
cart2sph = [ #: Transformation Between Cartesian and (Real) Pure Spherical Harmonic Gaussians
  [
  [[(0,0,0)], [1.], 1.]
  ],                                    # s orbitals
  [
  [[(0,1,0)], [1.], 1.],
  [[(0,0,1)], [1.], 1.],
  [[(1,0,0)], [1.], 1.],
  ],                                    # p orbitals
  [
  [[(1,1,0)], [1.], 1.],
  [[(0,1,1)], [1.], 1.],
  [[(0,0,2),(2,0,0),(0,2,0)], [1., -1/2., -1/2.], 1.],
  [[(1,0,1)], [1.], 1.],
  [[(2,0,0),(0,2,0)], [1.,-1.], sqrt(3)/2.],
  ],                                    # d orbitals
  [
  [[(0,3,0),(2,1,0)], [-sqrt(5), 3.], 1/(2.*sqrt(2))],
  [[(1,1,1)], [1.], 1.],
  [[(0,1,2),(0,3,0),(2,1,0)], [sqrt(3/5.), -sqrt(3)/4., -sqrt(3)/(4.*sqrt(5))], sqrt(2)] ,
  [[(0,0,3),(2,0,1),(0,2,1)], [1.,-3/(2*sqrt(5)),-3/(2*sqrt(5))], 1.],
  [[(1,0,2),(3,0,0),(1,2,0)], [sqrt(3/5.), -sqrt(3)/4., -sqrt(3)/(4.*sqrt(5))], sqrt(2)],
  [[(2,0,1),(0,2,1)], [1.,-1.], sqrt(3)/2.],
  [[(3,0,0),(1,2,0)], [sqrt(5), -3.], 1/(2.*sqrt(2))],
  ],                                    # f orbitals
  [
  [[(3,1,0), (1,3,0)], [1.,-1.], sqrt(2) * sqrt(5/8.)],
  [[(0,3,1), (2,1,1)], [-sqrt(5)/4.,3/4.], sqrt(2)],
  [[(1,1,2), (3,1,0), (1,3,0)], [3/sqrt(14), -sqrt(5)/(2*sqrt(14)), -sqrt(5)/(2*sqrt(14))], sqrt(2)],
  [[(0,3,1), (0,3,1), (2,1,1)], [sqrt(5/7.), -3*sqrt(5)/(4.*sqrt(7)), -3/(4.*sqrt(7))], sqrt(2)],
  [[(0,0,4), (4,0,0), (0,4,0), (2,0,2), (0,2,2), (2,2,0)], [1., 3/8., 3/8., -3*sqrt(3)/sqrt(35), -3*sqrt(3)/sqrt(35), -1/4.], sqrt(2)],
  [[(1,0,3), (3,0,1), (1,2,1)], [sqrt(5/7.), -3*sqrt(5)/(4.*sqrt(7)), -3/(4.*sqrt(7))], sqrt(2)],
  [[(2,0,2), (0,2,2), (4,0,0), (0,4,0)], [3*sqrt(3)/(2.*sqrt(14)), -3*sqrt(3)/(2.*sqrt(14)), -sqrt(5)/(4.*sqrt(2)), sqrt(5)/(4.*sqrt(2))], sqrt(2)],
  [[(3,0,1), (1,2,1)], [sqrt(5)/4., -3/4.], sqrt(2)],
  [[(4,0,0), (0,4,0), (2,2,0)], [sqrt(35)/(8.*sqrt(2)), sqrt(35)/(8.*sqrt(2)), -3*sqrt(3)/(4.*sqrt(2))], sqrt(2)],
  ],                                    # g orbitals
]

def get_cart2sph(l,m):
  '''Returns the linear combination required for the transformation Between 
  the Cartesian and (Real) Pure Spherical Harmonic Gaussian basis.
  
  Adapted from H.B. Schlegel and M.J. Frisch,
  International Journal of Quantum Chemistry, Vol. 54, 83-87 (1995).
  
  **Parameters:**
  
  l : int
    Angular momentum quantum number.
  m : int
    Magnetic quantum number.
  
  **Returns:**
  
  cart2sph[l][l+m] : list
    Contains the conversion instructions with three elements
      
      1. Exponents of Cartesian basis functions (cf. `core.exp`): list of tuples
      2. The corresponding expansion coefficients: list of floats 
      3. Global factor  
  
  ..hint: 
  
    The conversion is currently only supported up to g atomic orbitals.
  '''
  return cart2sph[l][l+m]

def validate_drv(drv):
  if drv is None or drv == 'None' or drv == '': return 0
  elif drv == 'x': return 1
  elif drv == 'y': return 2
  elif drv == 'z': return 3
  elif drv == 'xx' or drv == 'x2': return 4
  elif drv == 'yy' or drv == 'y2': return 5
  elif drv == 'zz' or drv == 'z2': return 6
  elif drv == 'xy' or drv == 'yx': return 7
  elif drv == 'xz' or drv == 'zx': return 8
  elif drv == 'yz' or drv == 'zy': return 9
  elif not (isinstance(drv,int) and 0 <= drv <= 9):
    raise ValueError("The selection `drv=%s` is not valid!"  % drv)
  else:
    return drv

def each_ao_is_normalized(ao_spec):
  is_normalized = []
  for sel_ao in range(len(ao_spec)):
    is_normalized.append((ao_spec[sel_ao]['pnum'] < 0))
  
  if all(is_normalized) != any(is_normalized):
    raise ValueError('Either all or none of the atomic orbitals have to be normalized!')
  return all(is_normalized)

def prepare_ao_calc(ao_spec):    
  pnum_list = []
  atom_indices = []
  ao_coeffs = numpy.zeros((0,2))  
  for sel_ao in range(len(ao_spec)):
    atom_indices.append(ao_spec[sel_ao]['atom'])
    c = ao_spec[sel_ao]['coeffs']
    ao_coeffs = numpy.append(ao_coeffs,c,axis=0)
    pnum_list.append(len(c))
      
  pnum_list = require(pnum_list, dtype='i')
  atom_indices = require(atom_indices, dtype='i')
  ao_coeffs = require(ao_coeffs, dtype='f')
  return ao_coeffs,pnum_list,atom_indices

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
      raise ValueError('%s has to be mo_spec or an numpy coefficient array.'%name)
  else:
    tmp = []
    for i in mo:
      tmp.append(i['coeffs'])
    mo = tmp
  mo = numpy.array(mo, dtype=numpy.float64)  
  if mo.ndim != 2:
    raise ValueError('%s has to be 2-dimensional.'%name)  
  return mo

def is_mo_spec(mo):
  '''Checks if :literal:`mo` is of :literal:`mo_spec` type. 
  (See :ref:`Central Variables` for details.)'''
  #Avoids circular inports:
  from orbkit.orbitals import MOClass
  if not isinstance(mo, MOClass):
    return False
  return_val = True
  for i in mo:
    try:
      return_val = return_val and 'coeffs' in i.keys()
    except:
      return_val = False
  
  return return_val

def require(data,dtype='f',requirements='CA'):
  if dtype == 'f':
    dtype = numpy.float64
  elif dtype == 'i':
    dtype = numpy.intc
  return numpy.require(data, dtype=dtype, requirements='CA')
