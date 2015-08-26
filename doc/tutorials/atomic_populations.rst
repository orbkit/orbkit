Tutorial for Atomic Population Analysis
=======================================

This short tutorial shows, how to perform a Mulliken and Löwdin population
analysis with orbkit analytically, and how to write the output to an ``xyz`` or
a ``PDB`` file.

Computation of the Population Analysis
--------------------------------------

First, we have to import some modules and set some of orbkits options::
  from orbkit import read, atomic_populations

Then, we have to read the input file::

  # Open molden file and read parameters
  qc = orbkit.main_read('h2o.md',itype='molden',all_mo=True)

Now, we have to call the respective functions:
  
For **Mulliken Population Analysis**::
  
  pop = atomic_populations.mulliken(qc)

and for **Löwdin Population Analysis**::
  
  pop = atomic_populations.lowdin(qc)

The return value ``pop`` is a dictionary, which contains information of 
population analysis and has following members:

+-----------------+----------------------------------------+
| Variable        | Contents                               |
+-----------------+----------------------------------------+
|``population``   | Contains the population for each atom. |
+-----------------+----------------------------------------+
|``charges``      | Contains the charges for each atom.    |
+-----------------+----------------------------------------+

Creation of an Output File
--------------------------

orbkit provides two possible output file formats for pop

Cubature requires a function with a special structure. It provides a 1d array
for the grid which we have bring ilowdinnto orbkit's shape. Then, we have to call 
orbkit's main computational function ``orbkit.rho_compute`` with the respective
options, i.e., we compute the norm of all occupied molecular orbitals in this 
example. Finally, we have to reshape orbkit's output::

  def func(x_array,*args):
    '''Calls orbkit.
    
    **Parameters:**
      
    x_array : 1-D numpy.ndarray, shape[0]=ndim*npt 
      Contains the grid.
    args : tuple or list
      |args[0] : int
      |  Contains the number of points (npt). 
    
    '''
    
    # Set the grid
    x_array = x_array.reshape((-1,3))
    orbkit.grid.x = numpy.array(x_array[:,0],copy=True)
    orbkit.grid.y = numpy.array(x_array[:,1],copy=True)
    orbkit.grid.z = numpy.array(x_array[:,2],copy=True)
    
    # Compute the squared molecular orbitals with orbkit
    out = orbkit.rho_compute(qc,
			     calc_mo=True,
			     vector=vector,
			     drv=None,
			     laplacian=False, 
			     numproc=numproc)
    out **= 2.
    
    # Return the results
    return numpy.reshape(numpy.transpose(out),(-1,))

.. note::
  
  We only consider the "vectorized" cubature which is much faster than the 
  single point cubature variant.

Run Cubature
------------

We integrate in each dimension from -5 to 5 :math:`{\rm a}_0` until we reach a 
relative error of 1e-3 without considering the absolute error of the integral::

  ndim = 3                                      #: Number of dimensions being integrated
  xmin = numpy.array([-5.,-5.,-5.],dtype=float) #: Minimum integration limit for each variable
  xmax = numpy.array([ 5., 5., 5.],dtype=float) #: Maximum integration limit for each variable
  abserr = 0    #: Absolute error: |error| < abserr (If zero, it will be ignored.)
  relerr = 1e-3 #: Relative error: |error| < relerr*|integral| (If zero, it will be ignored.)

  # Call the cubature routine together with orbkit.
  integral_mo,error_mo = cubature(ndim, func, xmin, xmax, 
				  args=[], 
				  adaptive='h', abserr=abserr, relerr=relerr, 
				  norm=0, maxEval=0, vectorized=True)

Finally, we can print the output::

  print('The integral of...')
  for i,(inte,err) in enumerate(zip(integral_mo,error_mo)):
    print('\tMO %s is %.14f. (Error: %.4e)' % (qc.mo_spec[i]['sym'],inte,err))
