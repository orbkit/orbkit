# -*- coding: iso-8859-1 -*-
'''Module for creating and manipulating the grid on which all computations 
are performed.'''
'''
orbkit
Gunter Hermann, Vincent Pohl, Lukas Eugen Marsoner Steinkasserer, Axel Schild, and Jean Christophe Tremblay

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

# Import general modules
import sys
import numpy

# Import orbkit modules
from orbkit import cy_grid

def grid_init(is_vector=False, force=False):
  '''Sets up the regular x-, y-, z-grid 
  specified by the global lists: 

    :min\_: List of minimum grid values
    :max\_: List of maximum grid values
    :N\_: List of number of grid points

  **Parameters:**
  
    is_vector : bool, optional
      If True, converts the regular grid to a vector grid.

  '''
  
  # All grid related variables should be globals 
  global x, y, z, d3r, min_, max_, N_, delta_, grid, is_initialized, is_regular
  
  if is_initialized and not force:
    return 0
  
  # Initialize a list for the grid 
  grid = [[],[],[]]
    
  # Loop over the three dimensions 
  for ii in range(3):
    if max_[ii] == min_[ii]:
      # If min-value is equal to max-value, write only min-value to grid  
      grid[ii]   = numpy.array([min_[ii]],dtype=numpy.float64)
      delta_[ii] = 1
      N_[ii]     = 1
    else:
      # Calculate the grid using the input parameters 
      if delta_[ii]:
        grid[ii] = numpy.arange(min_[ii],max_[ii]+delta_[ii],delta_[ii],dtype=numpy.float64)
        N_[ii] = len(grid[ii])
      else:
        grid[ii] = numpy.array(numpy.linspace(min_[ii],max_[ii],N_[ii]),dtype=numpy.float64) # backward compatibility
        delta_[ii] = grid[ii][1]-grid[ii][0]
  
  # Write grid 
  x = grid[0]  
  y = grid[1]  
  z = grid[2]
  d3r = numpy.product(delta_)
  
  is_initialized = True
  is_regular = True
  
  if is_vector:
    grid2vector()
  else:
    setattr(sys.modules[__name__],'is_vector',False)  
  # grid_init 

# Synonyms
init = grid_init
init_grid = grid_init

def get_grid(start='\t'):
  '''Returns a string describing the current x-, y-, z-grid.
  '''
  coord = ['x', 'y', 'z']
  grid = [x, y, z]
  display = ''
  for ii in range(3):
    display += ('%(s)s%(c)s[0] = %(min).2f %(c)s[-1] = %(max).2f N%(c)s = %(N)d ' % 
      {'s': start, 'c': coord[ii], 'min': grid[ii][0], 'max': grid[ii][-1], 
       'N': len(grid[ii])})
    if max_[ii] != min_[ii] and delta_[ii] != 0.:
      # Print the delta values only if min-value is not equal to max-value 
      display += 'd%(c)s = %(d).3f' % {'c': coord[ii], 'd': delta_[ii]}
    display += '\n'
  
  return display
  # get_grid 

def tolist():
  '''Returns a list containing the current x-, y-, z-grid.
  '''
  return [numpy.copy(x), numpy.copy(y), numpy.copy(z)]
  
def todict():
  '''Returns a dictionary containing the current x-, y-, z-grid.
  '''
  return {'x': x, 'y': y, 'z': z}

def get_shape():
 '''Returns the shape of the grid.
 '''
 if not is_initialized:
   raise ValueError('`grid.get_shape` requires the grid to be initialized.')
 
 return (len(x),) if is_vector else tuple(N_)

def set_grid(xnew,ynew,znew,is_vector):
  '''Sets the x-, y-, z-grid.
  '''
  global x, y, z, is_initialized, is_regular
  coord = ['x', 'y', 'z']  
  NumberTypes = (int, float) #: Contains the supported types.
  
  reset_grid() # Reset the current grid parameters

  length = []
  grid = [xnew,ynew,znew]
  for i,c in enumerate(grid):
    # Check the type of the grid
    if isinstance(c,NumberTypes):
      c = numpy.array([c],dtype=numpy.float64)      
    elif isinstance(c,(list,tuple)): 
      c = numpy.array(c,dtype=numpy.float64)    
    elif not isinstance(c,numpy.ndarray):
      raise TypeError('%s (dimension %d) is of inappropriate type. (%s)' %(coord[i],i,type(c)))
    # Reshape if necessary
    if c.ndim != 1:
      c = c.reshape((-1,))
    # Save new grid
    grid[i] = c
    length.append(len(c))
  
  # Produce some information about the grid.
  info_string = 'Grid has been set up...'
  info_string += ('\n\tIf the input coordinates will be used for a regular grid,' +
                  '\n\tit will contain %dx%dx%d=%d data points.' % 
                  (tuple(length) + (numpy.product(length),)) 
                  )
  if length[0] == length[1] == length[2]:
    info_string += ('\n\n\tIf the input coordinates will be used for a vector grid,' +
                    '\n\tit will contain %d data points.' % length[0]
                    )
  else:
    info_string += ('\n\n\tAttention: Due to their different length, the grid variables' +
                    '\n\tcannot be used for a computation using a vector grid!'
                    )
  
  # Write grid 
  x = grid[0]  
  y = grid[1]  
  z = grid[2]
  
  is_initialized = True
  if isinstance(is_vector,bool):
    setattr(sys.modules[__name__],'is_vector',is_vector)
    info_string += ('\n\nThe variable `grid.is_vector` has been set to %s.' % 
                    is_vector)

  is_regular = is_vector==False
  
  set_boundaries(is_regular)
  return info_string
  # set_grid 

def set_boundaries(is_regular,Nx=None,Ny=None,Nz=None):
  global min_, max_, delta_, N_
  min_ = [x.min(),y.min(),z.min()]
  max_ = [x.max(),y.max(),z.max()]
  N_ = [len(x),len(y),len(z)]
  if is_regular:
    f = lambda x: 1.0 if len(x) <= 1 else x[1]-x[0]
    delta_ = [f(x), f(y), f(z)]
  elif all([Nx,Ny,Nz]):
    N_ = [Nx,Ny,Nz]
    grid = numpy.array([x,y,z]).reshape(3,Nx,Ny,Nz)
    delta_ = [grid[1,0,0]-grid[0,0,0],
              grid[0,1,0]-grid[0,0,0],
              grid[0,0,1]-grid[0,0,0]]

def get_bbox():
  bbox = numpy.zeros(6)
  bbox[::2] = min_
  bbox[1::2] = max_
  return bbox

def grid2vector():
  '''Converts the regular grid characterized by x-, y-, z-vectors
  to a (3, (Nx*Ny*Nz)) grid matrix (vector grid). 
  Reverse operation: :mod:`orbkit.grid.vector2grid` 
  '''
  # All grid related variables should be globals 
  global x, y, z, is_vector, is_regular
  
  if not is_initialized:
    raise ValueError('You have to initialize a grid before calling '+
                 ' `grid.grid2vector`, i.e., `grid.is_initialized=True`.')
  
  x,y,z = cy_grid.grid2vector(x,y,z)    
  is_vector = True
  is_regular = True
  
def vector2grid(Nx,Ny,Nz):
  '''Converts the (3, (Nx*Ny*Nz)) grid matrix (vector grid) back to the regular grid 
  characterized by the x-, y-, z-vectors.
  Reverse operation: :mod:`orbkit.grid.grid2vector`
  '''
  # All grid related variables should be globals 
  global x, y, z, is_vector
  
  if not is_initialized:
    raise ValueError('You have to initialize a grid before calling '+
                     ' `grid.vector2grid`. '+
                     '(`grid.is_initialized == True`)')
  if not is_regular:
    raise ValueError('The grid has to regular. '+
                     '(`grid.is_regular == True`)')
  if not (len(x) == len(y) == len(z)):
    raise ValueError('Not a valid vector grid, i.e., dimensions of '+
                     'x-, y-, and z- coordinate differ.') 
  if (Nx*Ny*Nz) != len(x):
    raise ValueError('It has to hold that `len(x) = (N_x * N_y * N_z)`')
  
  x,y,z = cy_grid.vector2grid(x,y,z,Nx,Ny,Nz)  
  is_vector = False
  
def matrix_grid2vector(matrix): 
  '''Converts the (Nx,Ny,Nz) data matrix back to the regular grid (Nx,Nz,Ny)
  '''
  matrix = numpy.asarray(matrix,dtype=float)
  
  if matrix.ndim != 3:
    raise ValueError('`matrix` has to be 3d matrix.')
  
  return numpy.reshape(matrix,(-1,))
  
def matrix_vector2grid(matrix,Nx=None,Ny=None,Nz=None): 
  '''Converts the (Nx*Ny*Nz) data matrix back to the (Nx,Nz,Ny)
  '''
  
  matrix = numpy.asarray(matrix,dtype=float)
  
  if matrix.ndim != 1 or (Nx*Ny*Nz) != len(matrix):
    raise ValueError('`matrix` has to be one dimensional with the length '+
                     '`len(matrix) = N_x * N_y * N_z`.'+
                     'For Nd matrices use the function `grid.mv2g`.')
  
  return numpy.reshape(matrix,(Nx,Ny,Nz))

def mv2g(**kwargs):
  '''Converts all `numpy.ndarrays` given as the keyword arguments 
  (`**kwargs`) from a vector grid of `shape=(..., Nx*Ny*Nz, ...,)` to a regular 
  grid of `shape=(..., Nx, Ny, Nz, ...,)`, and, if more than one `**kwargs` is 
  given, returns it as a dictionary.
  
  Hint: The global values for the grid dimensionality, i.e., :mod:`grid.N_`,
  are used for reshaping.
  '''
  import itertools
  return_val = {}
  for i,j in kwargs.items():
    j = numpy.asarray(j,dtype=float)
    shape = numpy.shape(j)
    where = numpy.argwhere(shape==numpy.product(N_))[0,0]
    return_val[i] = numpy.zeros(shape[:where]+tuple(N_)+shape[where+1:])
    for key in itertools.product(*[range(k) for k in (shape[:where] + shape[where+1:])]):
      obj = [slice(k,k+1) for k in key]
      for r in range(3):
        obj.insert(where,slice(None,None))
      return_val[i][obj] = matrix_vector2grid(j[obj[:where]+obj[where+2:]].reshape((-1,)), 
                                          **dict(zip(['Nx','Ny','Nz'],N_)))

  return list(return_val.values())[0] if len(return_val.values()) == 1 else return_val

def grid_sym_op(symop):
  '''Executes given symmetry operation on vector grid 
  '''
  # All grid related variables should be globals 
  global x, y, z, is_vector, is_regular
  
  symop = numpy.asarray(symop,dtype=numpy.float64)  
  if symop.shape != (3,3):
    raise ValueError('`symop` needs to be a numpy array with shape=(3,3)')
  
  if not is_initialized:
    raise ValueError('You have to initialize a grid before executing a '+
                 ' symmetry operation on it. (`grid.is_initialized == True`)')
  
  if not is_vector:
    grid2vector()
  
  x, y, z = numpy.dot(symop,numpy.array([x,y,z]))
  
  # The grid is not regular anymore.
  is_regular = False
  

def grid_translate(dx,dy,dz):
  '''Translates the grid by (dx,dy,dz).
  '''
  global x,y,z
  x += dx
  y += dy
  z += dz

def rot(ang,axis):
 '''Creates matrix representation for arbitrary rotations
 Angle has to be defined in radians, e.g., numpy.pi/2.0
 Axis has to be specified as follows: 
 x-axis -> axis=0,
 y-axis -> axis=1,
 z-axis -> axis=2,
 '''
 # Initialize cosine, sinus, and additional numpy functions
 cos = numpy.cos
 sin = numpy.sin
 array = numpy.array
 insert = numpy.insert
 
 # Create rotation matrix around defined rotations axis
 rotmatrix = array([[ cos(ang), sin(ang)],
                 [-sin(ang), cos(ang)]])
 rotmatrix = insert(insert(rotmatrix,axis,0,axis=0),axis,0,axis=1)
 rotmatrix[axis,axis] = 1
 
 return rotmatrix
 # rot

def reflect(plane):
 '''Creates matrix representation for reflection
 Plane has to be specified as follows:
 xy-plane -> plane= numpy.array([0,1])
 xz-plane -> plane= numpy.array([0,2])
 yz-plane -> plane= numpy.array([1,2])
 '''
 
 # Create reflection matrix for defined plane
 sigma = numpy.array([[1,0,0],[0,1,0],[0,0,1]],dtype=float)
 axis = 3-numpy.sum(plane)
 sigma[axis,axis] *= -1.0

 return sigma
 # reflect

def inversion():
 '''Transfer matrix representation for inversion
 '''

 # Inversion matrix
 inv = numpy.array([[-1,0,0],[0,-1,0],[0,0,-1]],dtype=float) # inversion

 return inv
 # inversion

def sph2cart_vector(r,theta,phi):
  '''Converts a spherical regular grid matrix (r, theta, phi)
  to a Cartesian grid matrix (vector grid) with the shape (3, (Nr*Ntheta*Nphi)).

  **Parameters:**
  
    r : numpy.ndarray, shape=(Nr,)
      Specifies radial distance.
    theta : numpy.ndarray, shape=(Ntheta,)
      Specifies polar angle. 
    phi : numpy.ndarray, shape=(Nphi,)
      Specifies azimuth angle. 
  '''
  # All grid related variables should be globals 
  global x, y, z, is_initialized, is_vector, is_regular
  x,y,z = cy_grid.sph2cart(numpy.asarray(r,dtype=numpy.float64),
                      numpy.asarray(theta,dtype=numpy.float64),
                      numpy.asarray(phi,dtype=numpy.float64))
  
  is_initialized = True
  is_vector = True
  is_regular = False
  # sph2cart_vector 

def cyl2cart_vector(r,phi,zed):
  '''Converts a cylindrical regular grid matrix (r, phi, zed)
  to a Cartesian grid matrix (vector grid) with the shape (3, (Nr*Nphi*Nzed)).

  **Parameters:**
  
    r : numpy.ndarray, shape=(Nr,)
      Specifies radial distance.
    phi : numpy.ndarray, shape=(Nphi,)
      Specifies azimuth angle. 
    zed : numpy.ndarray, shape=(Nz,)
      Specifies z distance.
  '''
  # All grid related variables should be globals 
  global x, y, z, is_initialized, is_vector, is_regular
  
  x,y,z = cy_grid.cyl2cart(numpy.asarray(r,dtype=numpy.float64),
                      numpy.asarray(phi,dtype=numpy.float64),
                      numpy.asarray(zed,dtype=numpy.float64))
  
  is_initialized = True
  is_vector = True
  is_regular = False
  # cyl2cart_vector 

def random_grid(geo_spec,N=1e6,scale=0.5):
  '''Creates a normally distributed grid around the atom postions (geo_spec).

  **Parameters:**

    geo_spec : 
        See :ref:`Central Variables` for details.
    N : int
        Number of points distributed around each atom
    scale : float
        Width of normal distribution
  '''
  # All grid related variables should be globals 
  global x, y, z, is_initialized, is_vector, is_regular
  geo_spec = numpy.array(geo_spec)
  at_num = len(geo_spec)
  # Initialize a list for the grid 
  grid = numpy.zeros((3,at_num,N))
  
  # Loop over the three dimensions 
  for ii_d in range(3):
    for ii_a in range(at_num):
      grid[ii_d,ii_a,:] = numpy.random.normal(loc=geo_spec[ii_a,ii_d],scale=0.5,size=N)
  
  grid = numpy.reshape(grid,(3,N*at_num))

  # Write grid 
  x = grid[0]  
  y = grid[1]  
  z = grid[2]
  
  is_initialized = True
  is_vector = True
  is_regular = False  
  # random_grid 

def read(filename, comment='#'):
  '''Reads a grid from a plain text file.
  
  **Parameters:**
  
    fid : str
      Specifies the filename of the grid file. 
  
  **Returns:**
  
    is_vector : bool
      If True, a vector grid is used for the computations.

  **Supported Formats:**
  
    Regular Grid::
    
      # Regular Grid Example File
      # Format:
      # x xmin xmax Nx
      # y ymin ymax Ny
      # z zmin zmax Nz
      
      x -5 5 101
      y -2 2  51
      z  0 0   1 
    
    Vector Grid::
    
      # Vectorized Grid Example File
      # The header 'x y z' is mandatory!

      x       y       z
      -5      -5      0
      -4      -5      0
      -3      -5      0
      0       0       0
      2       -1e-1   9.78
  
  **Hint:** If a line starts with '#', it will be skipped. Please, do not use '#' at the end of a line!
  '''
  # All grid related variables should be globals 
  global x, y, z, min_, max_, N_, delta_, is_initialized, is_regular, is_vector
  
  def check(i, is_vector):
    if (len(i) == 3) and (is_vector is None or is_vector == True):
      return True
    elif (len(i) == 4) and (is_vector is None or is_vector == False):
      return False
    else:
      raise IOError('Inconsistency in Grid File in "%s"' % i) 

  # Go through the file line by line 
  is_vector = None
  is_dx = [False for i in range(3)] # Last column grid spacing not number of points

  grid = [[] for i in range(3)]
  dim = 'xyz'
  index = [[] for i in range(3)]

  with open(filename) as fileobject:
    for l,line in enumerate(fileobject):
      cl = line.split()                 # The Current Line split into segments
      
      if not (cl == [] or cl[0] == comment): 
        is_vector = check(cl, is_vector)
        if is_vector:
          for i,j in enumerate(cl):
            if index[i] == []: 
              index[i] = dim.find(j)
            else:              
              grid[index[i]].append(j)
        else:                  
          grid[dim.find(cl[0].lower())] = cl[1:]
          is_dx[dim.find(cl[0].lower())] = '.' in cl[-1]

  # Convert the variables 
  grid = numpy.array(grid,dtype=numpy.float64)
  if is_vector:
    x = grid[0,:]
    y = grid[1,:]
    z = grid[2,:]
    is_initialized = True # The grid will be seen as initialized
    is_regular = False # The grid is assumed to be non-regular
  else:
    min_ = grid[:,0]
    max_ = grid[:,1]
    for i in range(3):
      if is_dx[i]:
        delta_[i] = grid[i,2]
      else:
        N_[i] = int(grid[i,2])
  
  return is_vector   

def adjust_to_geo(qc,extend=5.0,step=0.1):
  '''Adjusts the grid boundaries to the molecular geometry.
  
  **Parameters:**
  
  qc : QCinfo class
    See :ref:`Central Variables` for details.
  extend : float
    Specifies the value by which the grid boundaries are extended in each 
    direction.
  step : float
    Specifies the grid spacing.
  
  '''
  global min_, max_, N_, delta_, is_vector, is_initialized
  
  for i in range(3):
    min_[i] = min(qc.geo_spec[:,i]) - (extend)
    max_[i] = max(qc.geo_spec[:,i]) + (extend)
    dist = (max_[i] - min_[i])
    N_[i] = int(numpy.ceil(dist/step))+1
    rest = (N_[i]-1)*step - dist
    # Correct minimum and maximum value, if necessary
    min_[i] -= rest/2.
    max_[i] += rest/2.
    delta_[i] = step
  
  is_vector = False
  is_initialized = False

def check_atom_select(atom,geo_info,geo_spec,interactive=True,
                      display=sys.stdout.write):
  if not((isinstance(atom, int)) and (0 < atom <= len(geo_spec))):
    display('Not a Valid atom number for centering the grid')
    display('Coose a valid index:')
    for i,j in enumerate(geo_info): 
      display('\t%d\t%s\t%s' % (i+1,j[0],geo_spec[i]))
    if interactive:
      while not((isinstance(atom, int)) and (0 < atom <= len(geo_spec))):
        try:
          atom = int(input('Please insert a correct index: '))
        except ValueError:
          pass
    else: 
      raise IOError('Insert a correct filename for the MO list!')
  return atom

def center_grid(ac,display=sys.stdout.write):
  '''Centers the grid to the point ac and to the origin (0,0,0).
  '''
  # All grid related variables should be globals 
  global x, y, z, d3r, min_, max_, N_, delta_
  
  if not is_regular:
    raise ValueError('Center grid is only supported for regular grids.')
  
  was_vector = is_vector
  if was_vector:
    vector2grid(*N_)
  
  P=[numpy.zeros((3,1)), numpy.reshape(ac,(3,1))]
  
  d_tilde = numpy.abs(P[0] - P[1])
  N_tilde = numpy.round(numpy.abs(d_tilde / delta_))
  
  for ii in range(3): 
    if N_tilde[ii] != 0:
      delta_[ii] = d_tilde[ii] / N_tilde[ii]
  
  grid = [x, y, z]
  
  for ii in range(3):
    if len(grid[ii]) != 1:
      position = numpy.nonzero(ac[ii] <= grid[ii])[0][0]
      g = numpy.abs(grid[ii][position] - ac[ii]);
      c = 1/2.*delta_[ii] - g;
      grid[ii] += c;
  
  x = grid[0]  
  y = grid[1]  
  z = grid[2]
  d3r = numpy.product(delta_)
  
  min_ = [min(grid[0]), min(grid[1]), min(grid[2])]
  max_ = [max(grid[0]), max(grid[1]), max(grid[2])]
  N_   = [len(grid[0]), len(grid[1]), len(grid[2])]
  
  display('Centered Grid to (%.2f %.2f %.2f): \n' % (ac[0], ac[1], ac[2]))
  display(get_grid())
  
  for ii in range(3):
    if len(numpy.nonzero(0. == numpy.round(grid[ii]*10000))[0])!= 0: 
      display('Warning!\n\tAt least one grid point is equal to zero.\n')
  
  if was_vector:
    grid2vector()
  # center_grid 

def reset_grid():
  '''Resets the grid parameters.'''
  global is_initialized, is_vector, is_regular, min_, max_, N_
  is_initialized = False
  is_vector = True
  is_regular = False
  min_ = [-8.0, -8.0, -8.0]
  max_ = [ 8.0,  8.0,  8.0]
  N_   = [ 101,  101,  101]
  delta_ = numpy.zeros((3,1))
  # reset_grid 

# Default values for the grid parameters 
min_ = [-8.0, -8.0, -8.0]   #: Specifies minimum grid values (regular grid).
max_ = [ 8.0,  8.0,  8.0]   #: Specifies maximum grid values (regular grid).
N_   = [ 101,  101,  101]   #: Specifies the number of grid points (regular grid).

# Initialize some lists 
x = numpy.array([0.0])        #: Contains the x-coordinates. 
y = numpy.array([0.0])        #: Contains the y-coordinates. 
z = numpy.array([0.0])        #: Contains the z-coordinates. 
delta_ = numpy.zeros((3,1))   #: Contains the grid spacing.
d3r = 0.0                     #: A volume element

is_initialized = False        #: If True, the grid is assumed to be initialized.
is_vector = True              #: If True, the grid is assumed to be vector grid.
is_regular = False            #: If True, the grid is assumed to be regular, i.e., a conversion of a vector grid to a regular grid is possible, if ``N_`` is set.
