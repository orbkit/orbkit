import string
import numpy
from os import path

from orbkit.units import u_to_me

nist_mass = None
# Standard atomic masses as "Linearized ASCII Output", see http://physics.nist.gov
nist_file = path.join(path.dirname(path.realpath(__file__)),
                      'supporting_data/Atomic_Weights_NIST.html')
# see http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii2&isotype=some


def read_nist():
  '''Reads and converts the atomic masses from the "Linearized ASCII Output", 
  see http://physics.nist.gov.
  '''
  global nist_mass

  f = open(nist_file,'r')
  flines = f.readlines()
  f.close()
  
  nist_mass = []
  index = None
  new = True
  
  def rm_brackets(text,rm=['(',')','[',']']):
    for i in rm:
      text = text.replace(i,'')
    return text
  
  for line in flines:
    thisline = line.split()
    if 'Atomic Number =' in line:
      i = int(thisline[-1]) - 1
      new = (i != index)
      if new:
        nist_mass.append(['',0])
      index = i
    elif 'Atomic Symbol =' in line and new:
      nist_mass[index][0] = thisline[-1]
    elif 'Standard Atomic Weight =' in line and new:
      nist_mass[index][1] = float(rm_brackets(thisline[-1]))

def standard_mass(atom):
  '''Returns the standard atomic mass of a given atom.
    
  **Parameters:**
  
  atom : int or str
    Contains the name or atomic number of the atom.
  
  **Returns:**
  
  mass : float
    Contains the atomic mass in atomic units.
  '''
  if nist_mass is None:
    read_nist()  
  try:
    atom = int(atom) - 1
    return nist_mass[atom][1] * u_to_me
  except ValueError:
    return dict(nist_mass)[atom.title()] * u_to_me
    
def get_atom_symbol(atom):
  '''Returns the atomic symbol of a given atom.
    
  **Parameters:**
  
  atom : int or str
    Contains the atomic number of the atom.
  
  **Returns:**
  
  symbol : str
    Contains the atomic symbol.
  '''
  if nist_mass is None:
    read_nist()  
  try:
    atom = int(atom) - 1
    return nist_mass[atom][0]
  except ValueError:    
    return atom.title()

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
  elif not isinstance(l,int):
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

def remove_from_list(inlist, value):
  return [item for item in inlist if item != value]

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

# Compatability function for old numpy versions
def moveaxis(array, source, target):
  transpose = array.transpose
  order = [n for n in range(array.ndim) if n not in source]
  for dest, src in sorted(zip(target, source)):
    order.insert(dest, src)
  result = transpose(order)
  return result

def require(data,dtype='f',requirements='CA'):
  if dtype == 'f':
    dtype = numpy.float64
  elif dtype == 'i':
    dtype = numpy.intc
  return numpy.require(data, dtype=dtype, requirements='CA')

def convert(data,was_vector,N):
  data = numpy.array(data,order='C')
  if not was_vector:
    data = data.reshape(data.shape[:-1] + N,order='C')
  return data

def zeros(shape,name,hdf5_file=None,chunks=True):
  if hdf5_file is None:
    return numpy.zeros(shape)
  else:
    return hdf5_file.create_dataset(name,shape,dtype=numpy.float64,chunks=chunks)

def reshape(data,shape,save_hdf5):
  if not save_hdf5:
    return data.reshape(shape)
  else:
    data.attrs['shape'] = shape
    return data[...].reshape(shape)
  
def print2D(x,format='%+.2f ',start='\t',end=''):
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
    print(s + end)

def pmat(matrix,vmax=lambda x: numpy.max(numpy.abs(x))):
  import matplotlib.pyplot as plt
  fig = plt.figure()
  if matrix.dtype == complex:
    print('plotting real part of matrix')
    matrix = matrix.real
  vm = vmax(numpy.abs(matrix)) if callable(vmax) else vmax
  plt.imshow(matrix,interpolation=None,vmin=-vm,vmax=vm,cmap='seismic_r')
  plt.colorbar()
  return fig
