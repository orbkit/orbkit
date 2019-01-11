# -*- coding: iso-8859-1 -*-
'''Angular Electron Density of (C3H3)+

Example for the orbkit low-level interface in the interplay with cubature.
(see http://orbkit.sourceforge.net/adtutorials/cubature.html)

This example is based on the cubature algorithm:
http://ab-initio.mit.edu/wiki/index.php/Cubature
https://github.com/saullocastro/cubature -> python interface

Due to the small system size we perform multithreading over the angle segments. 
'''
import numpy
from orbkit import read,core,grid,options,omp_functions

from time import time
import resource
from cubature import cubature

# Disable orbkit terminal output for each run
options.quiet = True
options.no_log = True
t = time()

def func(x_array,*args):
  '''Calls orbkit.
  
  **Parameters:**
    
  x_array : 1-D numpy.ndarray, shape[0]=ndim*npt 
    Vector of grid points for r, phi, z-values (in this order)
    This example function is only available for cubature parameter
    vectorized=True.
  args : tuple or list
    |args[0] : int
    |  Contains the number of points (npt). Variable does only exist if 
    |  vectorized=True.
    |args[-1] : bool
    |  If True, computes the molecular orbitals, and integrates their squared
    |  values iteratively.
  
  '''
  grid.is_initialized = True
  x_array = x_array.reshape((-1,3))
  grid.x = numpy.multiply(x_array[:,0],numpy.cos(x_array[:,1]))
  grid.y = numpy.multiply(x_array[:,0],numpy.sin(x_array[:,1]))
  grid.z = numpy.array(x_array[:,2],copy=True)
  
  out = numpy.array(core.rho_compute_no_slice(qc,
                                              calc_mo=args[-1],
                                              drv=None))
  
  if args[-1]: # calc_mo
    out **= 2.
  
  out = numpy.transpose(out*x_array[:,0])
  if args[-2]:
    # Old syntax of cubature
    out = numpy.reshape(out,(-1,))
  return out

def get_segment(iphi):
  '''Input function for `omp_functions.run`. 
  Computes one angular segment.
  
  **Parameters:**
  
  iphi : int
    Index of lower integration bound in phi.
    
  **Returns:**
  intVal : 1-D numpy.ndarray, shape[0]=ndim
    Contains the ndim integral values 
  error : 1-D numpy.ndarray, shape[0]=ndim
    Contains the ndim values of the integration error.
  
  '''
  
  # Specify the integration limits for each variable
  xmin = numpy.array([ 0.,phi[iphi]  +phi0, -8],dtype=float)
  xmax = numpy.array([15.,phi[iphi+1]+phi0,  8],dtype=float)
  # use cubature and function <okCylinder> to integrate
  try:    # old cubature syntax    
    intVal,error = cubature(3,func,xmin,xmax,args=[True,calc_mo],norm=0,\
                            adaptive='h',abserr=abserr,relerr=relerr,maxEval=0,\
                            vectorized=True)
  except: # new cubature syntax
    # number of return values of func
    if calc_mo: fdim = len(qc.mo_spec)
    else: fdim = 1
    intVal,error = cubature(func,3,fdim,xmin,xmax,args=[False,calc_mo],norm=0,\
                            adaptive='h',abserr=abserr,relerr=relerr,maxEval=0,\
                            vectorized=True)
  return intVal,error

def display(x):
  '''Naive print function.'''
  print(x)

numproc = 4
abserr = 1e-8 #: Specifies the absolute error: |error| < abserr (If zero, it will be ignored.)
relerr = 1e-8 #: Specifies the relative error: |error| < relerr*|integral| (If zero, it will be ignored.)

in_fid  = 'c3h3+.molden'

# Open molden file and read parameters
qc = read.main_read(in_fid,itype='molden',all_mo=False)

# set nuclear center of mass to zero 
qc.geo_spec = qc.geo_spec - qc.get_com()

# grid in angular direction
phi = numpy.linspace(0,2*numpy.pi,361)
iterphi = numpy.arange(len(phi)-1).reshape((-1,1))
dphi = phi[1]-phi[0]

# Starting point of the integration
phi0 = -numpy.pi/2.

# Define the indices for the angular segments
x = numpy.arange(len(phi)-1).reshape((-1,1))

# Integrate the density for each angular segment
calc_mo = False
return_val = numpy.array(omp_functions.run(get_segment,x=x,numproc=numproc,
                                           display=display))
integral = return_val[:,0] 
error = return_val[:,1]
print('Sum of angular integrals is %.8f. (Error: %.4e)' % 
      (integral.sum(),error.sum()))

# Integrate the norm of the molecular orbitals for each angular segment
calc_mo = True
return_val = numpy.array(omp_functions.run(get_segment,x=x,numproc=numproc
                                           ,display=display))
integral_mo = return_val[:,0] 
error_mo = return_val[:,1]

print('Sum of angular integrals is %.8f. (Error: %.4e)' % 
      (integral.sum(),error.sum()))

print('The sum of the angular integrals of...' )
for i,(inte,err) in enumerate(zip(integral_mo.T,error_mo.T)):
  print('\tMO %s is %.8f. (Error: %.4e)' % (qc.mo_spec[i]['sym'],inte.sum(),err.sum()))

print('The calculation took %.3fs' % (time()-t))
ram_requirement = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
print('and required %.2f MB of RAM' % (ram_requirement/1000.))

# Create a simple plot of the angular density
import pylab as plt
fig, axs = plt.subplots(2, 1,sharex='all')
for i in integral.T:
  axs[0].plot(numpy.rad2deg(phi[:-1]),i)

axs[0].set_ylabel(r'$\rho_{\rm ang}$')

for i in integral_mo.T:
  axs[1].plot(numpy.rad2deg(phi[:-1]),i)

axs[1].set_xlabel(r'$\phi$ ($^\circ$)')
axs[1].set_xlim(0,360)
axs[1].set_ylabel(r'$\rho^{\rm MO}_{\rm ang}$')

plt.show()
