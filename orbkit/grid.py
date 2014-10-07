# -*- coding: iso-8859-1 -*-
'''Module for creating and manipulating the grid on which all computations 
are performed.'''
'''
orbkit
Gunter Hermann, Vincent Pohl, and Axel Schild

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

# test how to import weave
try:
    from scipy import weave
except:    
    import weave

# Import orbkit modules
from orbkit import cSupportCode

def grid_init(is_vector=False):
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
  global x, y, z, d3r, min_, max_, N_, delta_, grid, is_initialized
  
  if is_initialized:
    return 0
  
  # Initialize a list for the grid 
  grid = [[],[],[]]
  delta_ = numpy.zeros((3,1))
  
  # Loop over the three dimensions 
  for ii in range(3):
    if max_[ii] == min_[ii]:
      # If min-value is equal to max-value, write only min-value to grid  
      grid[ii]   = numpy.array([min_[ii]])
      delta_[ii] = 1
    else:
      # Calculate the grid using the input parameters 
      delta_[ii] = (max_[ii]-min_[ii]) / float(N_[ii] - 1)
      grid[ii] = min_[ii] + numpy.arange(N_[ii]) * delta_[ii]
  
  # Write grid 
  x = grid[0]  
  y = grid[1]  
  z = grid[2]
  d3r = numpy.product(delta_)
  
  if is_vector:
    grid2vector()
  
  is_initialized = True
  
  return 0
  # grid_init 

def get_grid(start='\t'):
  '''Returns a string describing the current x-, y-, z-grid.
  '''
  coord = ['x', 'y', 'z']
  grid = [x,y,z]
  display = ''
  for ii in range(3):
    display += ('%(s)s%(c)smin = %(min).2f %(c)smax = %(max).2f N%(c)s = %(N)d ' % 
      {'s': start, 'c': coord[ii], 'min': min(grid[ii]), 'max': max(grid[ii]), 
       'N': len(grid[ii])})
      #{'s': start, 'c': coord[ii], 'min': min_[ii], 'max': max_[ii], 'N': N_[ii]})
    if max_[ii] != min_[ii] and delta_[ii] != 0.:
      # Print the delta values only if min-value is not equal to max-value 
      display += 'd%(c)s = %(d).3f' % {'c': coord[ii], 'd': delta_[ii]}
    display += '\n'
  
  return display
  # get_grid 

def set_grid(grid):
  '''Returns a string describing the current x-, y-, z-grid.
  '''
  global x, y, z
  coord = ['x', 'y', 'z']
  delta_ = numpy.zeros((3,1)) #: Contains the grid spacing.
  
  NumberTypes = (int, long, float) #: Contains the supported types.
  
  # Check the input variable
  correct_type = (isinstance(grid,list) or isinstance(grid,numpy.ndarray))
  if not correct_type or len(grid) != 3:
    raise TypeError('The `grid` variable has to be a list or a numpy array with ' + 
                    'three dimensions.')
  
  length = []
  for i,c in enumerate(grid):
    # Check the type of the grid
    if isinstance(c,NumberTypes):
      c = numpy.array([c],dtype=float)      
    elif isinstance(c,(list,tuple)): 
      c = numpy.array(c,dtype=float)    
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
  
  return info_string
  # set_grid 

def reset_grid():
  '''Resets the grid parameters.'''
  global is_initialized, min_, max_, N_
  is_initialized = False
  min_ = [-8.0, -8.0, -8.0]
  max_ = [ 8.0,  8.0,  8.0]
  N_   = [ 101,  101,  101]
  # reset_grid 

def grid2vector():
  '''Converts the regular grid characterized by x-, y-, z-vectors
  to a (3, (Nx*Ny*Nz)) grid matrix (vector grid). 
  Reverse operation: :mod:`orbkit.grid.vector2grid` 
  '''
  # All grid related variables should be globals 
  global x, y, z, is_vector
  
  # Initialize a list for the grid 
  grid = numpy.zeros((3,numpy.product(N_)))
  
  grid_code = """
  int count=0;
  
  for (int i=0; i<Nx[0]; i++)
  {
    for (int j=0; j<Ny[0]; j++)
    {
      for (int k=0; k<Nz[0]; k++)
      {
        GRID2(0,count) = x[i];
        GRID2(1,count) = y[j];
        GRID2(2,count) = z[k];
        count += 1;
      }
    }
  }
  """
  weave.inline(grid_code, ['x','y','z','grid'], verbose = 1, support_code = cSupportCode.math)
  
  # Write grid 
  x = grid[0,:]  
  y = grid[1,:]  
  z = grid[2,:]
  
  is_vector = True
  return 0
  # grid2vector 
  
def vector2grid(Nx=None,Ny=None,Nz=None):
  '''Converts the (3, (Nx*Ny*Nz)) grid matrix (vector grid) back to the regular grid 
  characterized by the x-, y-, z-vectors.
  Reverse operation: :mod:`orbkit.grid.grid2vector`
  '''
  # All grid related variables should be globals 
  global x, y, z, is_vector
  
  # Initialize a list for the grid
  grid = [[],[],[]]
  grid[0] = x
  grid[1] = y
  grid[2] = z
  grid = numpy.array(grid,dtype=float)
  
  # Initialize x, y, and z for C++-code
  x = numpy.zeros(Nx)
  y = numpy.zeros(Ny)
  z = numpy.zeros(Nz)
  
  grid_code = """
  int count=0;
  for (int i=0; i<Nz[0]; i++)
  {
    Z1(i) = GRID2(2,i);
    for (int j=0; j<Ny[0]; j++)
    {
      Y1(j) = GRID2(1,j*Nz[0]);
      for (int k=0; k<Nx[0]; k++)
      {
        X1(k) = GRID2(0,(k*Nz[0]*Ny[0]));
      }  
    }
  }
  """
  weave.inline(grid_code, ['x','y','z','grid'], verbose = 1, support_code = cSupportCode.math)
  
  is_vector = False
  
  return 0
  # vector2grid 

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
  global x, y, z, is_initialized, is_vector
  
  grid = numpy.zeros((3,numpy.product([len(r),len(theta),len(phi)])))
  grid_code = """
  int count=0;

  for (int i=0; i<Nr[0]; i++)
  {
    for (int j=0; j<Ntheta[0]; j++)
    {
      for (int k=0; k<Nphi[0]; k++)
      {
        GRID2(0,count) = r[i] * sin(theta[j]) * cos(phi[k]);
        GRID2(1,count) = r[i] * sin(theta[j]) * sin(phi[k]);
        GRID2(2,count) = r[i] * cos(theta[j]);
        count += 1;
      }
    }
  }
  """
  weave.inline(grid_code, ['r','theta','phi','grid'], verbose = 1, support_code = cSupportCode.math)

  # Write grid 
  x = grid[0,:]  
  y = grid[1,:]  
  z = grid[2,:]
  
  is_initialized = True
  is_vector = True
  
  return 0
  # sph2cart_vector 

def random_grid(geo_spec,N=1e6,scale=0.5):
  '''Creates a normally distributed grid around the atom postions (geo_spec).

  **Parameters:**

    geo_spec : 
        See `Central Variables`_ for details.
    N : int
        Number of points distributed around each atom
    scale : float
        Width of normal distribution
  '''
  # All grid related variables should be globals 
  global x, y, z, is_initialized, is_vector
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
  
  return 0
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
  
    Regular Grid:
    
      .. literalinclude:: ../examples/grid_reg.txt
            :language: bash
    
    Vector Grid:
    
      .. literalinclude:: ../examples/grid_vec.txt
            :language: bash
  
  **Hint:** If a line starts with '#', it will be skipped. Please, do not use '#' at the end of a line!
  '''
  # All grid related variables should be globals 
  global x, y, z, min_, max_, N_, is_initialized
  
  def check(i, is_vector):
    if (len(i) == 3) and (is_vector is None or is_vector == True):
      return True
    elif (len(i) == 4) and (is_vector is None or is_vector == False):
      return False
    else:
      raise IOError('Inconsistency in Grid File in "%s"' % i) 

  # Go through the file line by line 
  is_vector = None

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

  # Convert the variables 
  grid = numpy.array(grid,dtype=numpy.float64)
  if is_vector:
    x = grid[0,:]
    y = grid[1,:]
    z = grid[2,:]
    is_initialized = True # The grid will be seen as initialized
  else:
    min_ = grid[:,0]
    max_ = grid[:,1]
    N_   = numpy.array(grid[:,2],dtype=int)
  
  return is_vector   

def center_grid(ac,display=sys.stdout.write):
  '''Centers the grid to the point ac and to the origin (0,0,0).
  '''
  # All grid related variables should be globals 
  global x, y, z, d3r, min_, max_, N_, delta_
  
  P=[numpy.zeros((3,1)), numpy.reshape(ac,(3,1))]
  
  d_tilde = numpy.abs(P[0] - P[1])
  N_tilde = numpy.round(numpy.abs(d_tilde / delta_))
  
  for ii in range(3): 
    if N_tilde[ii] != 0:
      delta_[ii] = d_tilde[ii] / N_tilde[ii]
  
  grid = [x, y, z]
  
  for ii in range(3):
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
  
  return 0
  # center_grid 

# Default values for the grid parameters 
min_ = [-8.0, -8.0, -8.0]   #: Specifies minimum grid values (regular grid).
max_ = [ 8.0,  8.0,  8.0]   #: Specifies maximum grid values (regular grid).
N_   = [ 101,  101,  101]   #: Specifies the number of grid points (regular grid).

# Initialize some lists 
x = numpy.array([0])        #: Contains the x-coordinates. 
y = numpy.array([0])        #: Contains the y-coordinates. 
z = numpy.array([0])        #: Contains the z-coordinates. 
delta_ = numpy.zeros((3,1)) #: Contains the grid spacing.

is_initialized = False      #: If True, the grid is assumed to be initialized.
is_vector = False           #: If True, the grid is assumed to be vectorized.