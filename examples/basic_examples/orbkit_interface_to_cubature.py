# -*- coding: iso-8859-1 -*-
'''Orbkit interface to cubature:
http://ab-initio.mit.edu/wiki/index.php/Cubature
https://github.com/saullocastro/cubature -> python interface
'''
import numpy
from orbkit import read,core,grid,options
from orbkit.analytical_integrals import get_ao_overlap,get_mo_overlap

from time import time
import resource

create_plot = False     #: Warning! Do not create plots when using a high precision.
save_data_values = False #: Saves the grid and output values computed. This requires more RAM and time.

# create_plot needs the data values
if create_plot and not save_data_values:
  save_data_values = True

numproc = 4             #: Specifies number of subprocesses.
slice_length = 1e4      #: Specifies number of points per subprocess.

in_fid  = 'h2o.molden'      #: Specifies input file name.

# Open molden file and read parameters
qc = read.main_read(in_fid,itype='molden',all_mo=False)

# Compute atomic orbital overlap matrix 
ao_overlap_matrix = get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec)

# Compute the overlap of the molecular orbitals and weight it with the occupation number
analytical_integral = 0.
for i_mo in qc.mo_spec:
  analytical_integral += i_mo['occ_num'] * get_mo_overlap(i_mo['coeffs'],
                                                      i_mo['coeffs'],
                                                      ao_overlap_matrix)

print('Analytical Integral: %.12f' % analytical_integral)

# Disable orbkit terminal output for each run
options.quiet = True
options.no_log = True


# Initialize some variables
t = [time()]
if save_data_values:
  vec = numpy.zeros((0,3))
  rho = numpy.zeros((0,))
  vec_mo = numpy.zeros((0,3))
  mos = numpy.zeros((len(qc.mo_spec),0))

def func(x_array,*args):
  '''Calls orbkit.
  
  **Parameters:**
    
  x_array : 1-D numpy.ndarray, shape[0]=ndim*npt if vectorized=True else ndim
    Contains the grid.
  args : tuple or list
    |args[0] : int
    |  Contains the number of points (npt). Variable does only exist if 
    |  vectorized=True.
    |args[-2] : bool
    |  If True, assumes a vectorized grid. (vectorized=True)
    |args[-1] : bool
    |  If True, computes the molecular orbitals, and integrates their squared
    |  values iteratively.
  
  ** Returns:**
  
  out : 1-D numpy.ndarray, shape[0]=ndim*npt if vectorized=True else ndim
    Contains the output data.  
  '''
  global count_calls,vec,rho,vec_mo,mos
  
  if args[-2]: # vectorized
    x_array = x_array.reshape((-1,3))
    grid.x = numpy.array(x_array[:,0],copy=True)
    grid.y = numpy.array(x_array[:,1],copy=True)
    grid.z = numpy.array(x_array[:,2],copy=True)
    
    count_calls += args[0] # Number of points
  else:
    x,y,z = x_array
    x_array = x_array[numpy.newaxis,:]
    
    grid.x = numpy.array([x])
    grid.y = numpy.array([y])
    grid.z = numpy.array([z])
    
    count_calls += 1
  
  # We have already initialized a grid for orbkit!
  grid.is_initialized = True
  
  # Run orbkit
  out = core.rho_compute(qc,
                         calc_mo=args[-1],
                         slice_length=slice_length,
                         drv=None,
                         numproc=numproc)
  
  if save_data_values:
    if args[-1]:
      vec_mo = numpy.append(vec_mo,x_array,axis=0)
      mos = numpy.append(mos,out,axis=1)
    else:
      vec = numpy.append(vec,x_array,axis=0)
      rho = numpy.append(rho,out,axis=0)
  
  if args[-1]: # calc_mo
    out **= 2. 
  
  return out.transpose()

'''
Cubature
========
'''
from cubature import cubature
ndim = 3                                         #: Specifies the number of dimensions being integrated
fdim = 1                                         #: Specifies the length of the output vector of func
xmin = numpy.array([-20.,-20.,-20.],dtype=float) #: Specifies the minimum integration limit for each variable
xmax = numpy.array([ 20., 20., 20.],dtype=float) #: Specifies the maximum integration limit for each variable
abserr = 1e-3                                    #: Specifies the absolute error: |error| < abserr requested (If zero, it will be ignored.)
relerr = 1e-3                                    #: Specifies the relative error: |error| < relerr*|integral| requested (If zero, it will be ignored.)
vectorized = True                                #: If True, uses a vector of points in cubature, instead of a single point calculation. (Much faster!)

'''
Example One
===========
Run the cubature routine together with orbkit and integrate the density.
'''
count_calls = 0
calc_mo = False

# Call the cubature routine together with orbkit.
try:
  integral,error = cubature(func, ndim, fdim, xmin, xmax, 
                            args=[vectorized,calc_mo], 
                            adaptive='h', abserr=abserr, relerr=relerr, 
                            norm=0, maxEval=0,
                            vectorized=vectorized)
except TypeError:
  raise RuntimeError('Calculation failed. Should work with cubature 0.13.1 which has a different syntax than earlier versions.')

print('After %d function calls the integral is %.14f. (Error: %.4e)' % 
      (count_calls,integral,error))

'''
Example Two
===========
Run the cubature routine together with orbkit  and integrate all MOs (squared) 
at once.

Note: The function saves the molecular orbital values. (Not the squared values!)
'''
count_calls = 0
calc_mo = True
fdim = len(qc.mo_spec)           #: now func will return the number of orbitals

# Call the cubature routine together with orbkit.
try:
  integral_mo,error_mo = cubature(func, ndim, fdim, xmin, xmax, 
                                  args=[vectorized,calc_mo], 
                                  adaptive='h', abserr=abserr, relerr=relerr, 
                                  norm=0, maxEval=0, vectorized=vectorized)
except TypeError:
  raise RuntimeError('Calculation failed. Should work with cubature 0.13.1 which has a different syntax than earlier versions.')

print('After %d function calls the integral of...' % count_calls)
for i,(inte,err) in enumerate(zip(integral_mo,error_mo)):
  print('\tMO %s is %.14f. (Error: %.4e)' % (qc.mo_spec[i]['sym'],inte,err))

# Print a final comment on the required time and RAM usage.
t.append(time())
print('The calculation took %.3fs' % (t[-1]-t[0]))
ram_requirement = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
print('and required %.2f MB of RAM' % (ram_requirement/1000.))

if create_plot and save_data_values:
  # Create a 3d plot of the grid points including a colorbar for the data values
  from mpl_toolkits.mplot3d import Axes3D
  import matplotlib.pyplot as plt
  from matplotlib import colors
  
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  
  ax.scatter(vec[:,0],vec[:,1], vec[:,2], c=numpy.log10(rho), cmap='cool',
             marker='.',linewidths=0,alpha=0.6)
  
  # Plot the geometry
  c = qc.geo_spec
  ax.scatter(c[:,0],c[:,1], c[:,2], c='k', marker='x',s=100.,lw=1.5)
  
  ax.set_xlabel('x')
  ax.set_ylabel('y')
  ax.set_zlabel('z')
  
  plt.show()
