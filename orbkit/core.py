# -*- coding: iso-8859-1 -*-
'''Module performing all computational tasks.'''

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

# Import general modules
import string
import time

import numpy

from multiprocessing import Pool

# Import orbkit modules
from orbkit import grid,cSupportCode
from orbkit.display import display

# test how to import weave
try:
    from scipy import weave
except:    
    import weave

def l_creator(geo_spec,ao_spec,sel_ao,exp_list=None,coeff_list=None,
              at_pos=None,is_vector=False,drv=None,
              x=None,y=None,z=None,N=None):
  '''Calculates the contracted atomic orbitals of quantum number l or its
  derivative with respect to a specific variable (e.g. drv = 'x' or drv = 0)
  for the atomic orbitals: ao_spec[sel_ao].
  
  **Parameters:**
  
  geo_spec,ao_spec :
    See `Central Variables`_ for details.
  sel_ao : int
    Index of the requested atomic orbitals
  exp_list : array_like, shape=(NDEG, 3), optional
    If not None, list of xyz-exponents of the NDEG 
    degenerate atomic orbitals ,i.e., NDEG=len(exp_list),
    else the standard molden exponents (exp) for quantum number l will be 
    used.
  coeff_list : array_like, shape=(PNUM, 2), optional
    If not None, list of the PNUM primitive atomic orbital
    exponents [:,0] and coefficients [:,1],
    else the coefficients from ao_spec[sel_ao] will be used.
  at_pos : array_like, shape=(3,), optional
    If not None, xyz-coordinates where the atomic orbital is centered,
    else the position geo_spec[ao_spec[sel_ao]['atom']] will be used.
  is_vector : bool, optional
    If True, a vectorized grid will be applied.
  drv : int or string, {None, 'x', 'y', 'z', 0, 1, 2}, optional
    If not None, a derivative calculation of the atomic orbitals 
    is requested.
  compile_only : bool, optional
    If True, compiles only the C++ code.
  x,y,z : list of floats, optional
    If not None, provides a list of Cartesian coordinates, 
    else the respective coordinates of the module :mod:`orbkit.grid` will 
    be used.
  N : tuple
    If not None, provides the shape of the grid.
  
  **Returns:**
  
  ao_list : numpy.ndarray, shape=((NDEG,) + N)
      Contains the computed NGED=len(exp_list) atomic orbitals on a grid.
  
  **Information:**
  
  We use scipy.weave.inline to run C++ Code within the python environment.
  http://docs.scipy.org/doc/numpy/user/c-info.python-as-glue.html#inline-c-code
  '''
  # Create the grid
  if x is None: x = grid.x
  if y is None: y = grid.y
  if z is None: z = grid.z
  
  if not is_vector:
    if N is None: N = numpy.array(grid.N_)
  else:
    if len(x) != len(y) or len(x) != len(z):
      display("Dimensions of x-, y-, and z- coordinate differ!")
      return 0
    else:
      N = (len(x),)
  
  # Build up the numpy arrays for the AO compuation
  if exp_list is None:
    l = lquant[ao_spec[sel_ao]['type']]
    exp_list = numpy.array(exp[l])
  if coeff_list is None:
    coeff_list = numpy.array(ao_spec[sel_ao]['coeffs'])
  if at_pos is None:
    at_pos = numpy.array(geo_spec[ao_spec[sel_ao]['atom']])
  ao_list = numpy.zeros(((len(exp_list),) + tuple(N)))
  
  # number of primitive atomic oribtals
  ao_num = numpy.shape(coeff_list)[0]
  
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
  
  # Ask for respective the C++ code
  code = ao_code(is_vector=is_vector,is_drv=(drv is not None))  
  # A list of Python variable names that should be transferred from
  # Python into the C/C++ code. 
  arg_names = ['x','y','z','ao_num','exp_list',
               'coeff_list','at_pos','ao_list','drv']
  # A string of valid C++ code declaring extra code
  support_code = cSupportCode.norm + cSupportCode.xyz
  
  try:
    # Compute the atomic orbitals
    weave.inline(code, arg_names = arg_names, 
            support_code = support_code,verbose = 1)
  except (weave.build_tools.CompileError, ImportError):
    display(('-'*80) + '''
    You have tried to compile the C++ inline code simulatniously using mulitiple
    processes... Please ignore the above messages!
    Waiting 2s for a second attempt...\n''' + ('-'*80))
    time.sleep(2)
    
    # Compute the atomic orbitals
    weave.inline(code, arg_names = arg_names, 
                 support_code = support_code,verbose = 1)
    

  return ao_list
  # l_creator 

def ao_creator(geo_spec,ao_spec,exp_list=False,
               is_vector=False,drv=None,
               x=None,y=None,z=None,N=None):
  '''Calculates all contracted atomic orbitals or its
  derivatives with respect to a specific variable (e.g. drv = 'x' or drv = 0).
  
  **Parameters:**
  
  geo_spec,ao_spec :
    See `Central Variables`_ in the manual for details.
  sel_ao : int
    Index of the requested atomic orbital
  exp_list : bool, optional
    If True, takes the xyz-exponents from ao_spec[i]['Exponents'],
    else the standard molden exponents (exp) for quantum number l will be 
    used.
  is_vector : bool, optional
    If True, a vectorized grid will be applied
  drv : int or string, {None, 'x', 'y', 'z', 0, 1, 2}, optional
    If not None, an analytical  calculation of the derivatives for 
    the atomic orbitals with respect to DRV is requested.
  x,y,z : None or list of floats, optional
    If not None, provides a list of Cartesian coordinates, 
    else the respective coordinates of grid. will be used
  N : None or tuple, optional
    If not None, provides the shape of the grid.
  
  **Returns:**
  
  ao_list : numpy.ndarray, shape=((NAO,) + N)
    Contains the computed NAO atomic orbitals on a grid.
  '''
  # Create the grid
  if x is None: x = grid.x
  if y is None: y = grid.y
  if z is None: z = grid.z

  if not is_vector:
    if N is None: N = numpy.array(grid.N_)
  else:
    if len(x) != len(y) or len(x) != len(z):
      display("Dimensions of x-, y-, and z- coordinate differ!")
      return 0
    else:
      N = (len(x),)
  
  if not exp_list:
    ii_exp = None # Use the standard molden xyz-exponents (exp)
  
  # Generalized AO creator
  for ii in range(len(ao_spec)):
    if exp_list or 'type' not in ao_spec[ii].keys():
      # Read the user-defined xyz-exponents
      ii_exp = numpy.array([ao_spec[ii]['exp_list']]) # Exponents
    
    # Compute the atomic orbitals
    if not ii:
      ao_list = l_creator(geo_spec,ao_spec,ii,drv=drv,
                x=x,y=y,z=z,N=N,is_vector=is_vector,
                exp_list=ii_exp)
    else:
      ao_list = numpy.append(ao_list,l_creator(geo_spec,ao_spec,ii,drv=drv,
                    x=x,y=y,z=z,N=N,is_vector=is_vector,
                    exp_list=ii_exp), axis = 0)
    
  return ao_list
  # ao_creator 

def calc_single_mo(xx):
  '''Computes a single molecular orbital. 
  
  This function is called by the multiprocessing module in the :mod:`orbkit.core.mo_creator`.
  
  **Parameters:**
  
  xx : int
    Specifies which molecular orbital shall be computed, i.e., mo_spec[xx].
  Spec : dict, global
    Dictionary containing all required varibles:
      :ao_list: : List of atomic orbitals on a grid.
      :mo_spec: : List of dictionaries (see `Central Variables`_ for details).
      :N: : Tuple containing the shape of the grid.

  **Returns:**
  
  mo : numpy.ndarray, shape=(N)
    Contains the molecular orbitals on a grid.
  '''
  try:
    ao_list = Spec['ao_list']
    mo_spec = Spec['mo_spec']
    N = Spec['N']
    mo = numpy.zeros(N)
    for jj in range(len(ao_list)):
      mo += mo_spec[xx]['coeffs'][jj] * ao_list[jj]
    return mo
  except KeyboardInterrupt:
    # Catch keybord interrupt signal to prevent a hangup of the worker processes 
    return 0

def mo_creator(ao_list,mo_spec,is_vector=False,
            x=None,y=None,z=None,N=None,mo_coeff=None,
            HDF5_save=False,h5py=False,
            numproc=1,s=0):
  '''Calculates the molecular orbitals.
  
  If a string (filename) is given for the argument :literal:`HDF5_save`, the slice s of each
  molecular orbital will be saved to the disk. This module is used by
  :mod:`orbkit.extras.save_mo_hdf5()`, where the HDF5 file is initialized.
  
  **Parameters:**
  
  ao_list : numpy.ndarray, shape=((NAO,) + N)
    Contains the NAO atomic orbitals on a grid.
  mo_spec : List of dictionaries
    See `Central Variables`_ for details.
  is_vector : bool, optional
    If True, a vectorized grid will be applied.
  x,y,z : None or list of floats, optional
    If not None, provide a list of Cartesian coordinates, 
    else the respective coordinates of grid will be used.
  N : None or tuple, optional
    If not None, provides the shape of the grid.
  HDF5_save : False or string, optional
    If not False, filename of HDF5 file for storing the molecular orbitals.
    (Requires Parameters: h5py and s)
  h5py : python module
    required if HDF5_save is not False
  numproc : int, required if HDF5_save is not False
    Specifies number of subprocesses for multiprocessing.
  s : int, required if HDF5_save is not False
    Specifies which slice of the molecular orbital has to be computed.

  **Returns:**
  
  mo_list : numpy.ndarray, shape=((NMO,) + N)
    Contains the NMO=len(mo_spec) molecular orbitals on a grid.
  '''
  
  if x is None: x = grid.x
  if y is None: y = grid.y
  if z is None: z = grid.z
  
  if not is_vector:
    if N is None: N = tuple(grid.N_)
  else:
    if len(x) != len(y) or len(x) != len(z):
      display("Dimensions of x-, y-, and z- coordinate differ!")
      return 0
    else:
      N = (len(x),)
  
  if not HDF5_save:
    # Standard mo_creator 
    mo_list = []
    if mo_coeff is None:
      for ii in range(len(mo_spec)):
        mo_list.append(numpy.zeros(N))
        for jj in range(len(ao_list)):
          mo_list[ii] += mo_spec[ii]['coeffs'][jj] * ao_list[jj]
    else:
      for ii in range(len(mo_coeff)):
        mo_list.append(numpy.zeros(N))
        for jj in range(len(ao_list)):
          mo_list[ii] += mo_coeff[ii][jj] * ao_list[jj]
    return mo_list
  else:
    # Save the MOs directly to an HDF5 file 
    f = h5py.File(HDF5_save, 'a')
    
    global Spec    
    Spec = {'ao_list': ao_list, 'mo_spec': mo_spec, 'N': N}
    
    if numproc > len(mo_spec): numproc = len(mo_spec)
    
    # Start the worker processes --
    pool = Pool(processes=numproc)
    
    # Write the slices in x to an array xx 
    xx=[]
    for ii in range(len(mo_spec)):
      xx.append(ii)
    
    it = pool.imap(calc_single_mo, xx)
    
    for ii in range(len(mo_spec)):
      dID = 'MO:%s' % mo_spec[ii]['sym']
      a = it.next()
      f[dID][s,:,:] = a
    
    # Close the worker processes --
    
    pool.close()
    pool.join()
    
    f.close()
    return 0
  # mo_creator 

def slice_rho(xx):
  '''Calculates the density, the molecular orbitals, or the derivatives thereof
  with respect to Spec['Derivative'] for one slice (xx)
  
  This function is called by the multiprocessing module in the :mod:`orbkit.core.rho_compute`.
  
  **Parameters:**
  
  xx : [float] or [int, int]
    Specifies which slice in x-direction shall be computed.
    
      | **If not is_vector:** One slice at x=xx will be computed.
      | **Else:**  One slice from index xx[0] to xx[1] will be calculated.
      
  Spec : dict, global
    Dictionary containing all required varibles:
      :geo_spec: List of floats, shape=(NATOMS, 3) (see `Central Variables`_ for details).
      :ao_spec: List of dictionaries (see `Central Variables`_ for details).
      :mo_spec: List of dictionaries (see `Central Variables`_ for details).
      :calc_mo: Bool if only the molecular orbitals are requested.
      :is_vector: Bool if a vectorized grid is used.
      :Derivative: List of strings, choices={'x','y', or 'z'}. 
                   If not None, derivative calculation will be carried out.
  grid : module or class, global
    Contains the grid, i.e., grid.x, grid.y, and grid.z.

  **Returns:**
  
  :if calc_mo and drv is None: 
    - mo_list
  :if calc_mo and drv is not None:
    - delta_mo_list
  :if not calc_mo and drv is None: 
    - rho, mo_norm
  :if not calc_mo and drv is not None: 
    - rho, mo_norm, delta_rho
  
  mo_list : numpy.ndarray, shape=((NMO,) + N)
    Contains the NMO=len(mo_spec) molecular orbitals on a grid.
  delta_mo_list : numpy.ndarray, shape=((NDRV,NMO) + N)
    Contains the derivatives with respect to drv (NDRV=len(drv)) of the 
    NMO=len(mo_spec) molecular orbitals on a grid.
  rho : numpy.ndarray, shape=(N)
    Contains the density on a grid.
  delta_rho : numpy.ndarray, shape=((NDRV,) + N)
    Contains the derivatives with respect to drv (NDRV=len(drv)) of 
    the density on a grid.
  '''
  try:
    # All desired information is stored in the Global variable Spec 
    geo_spec = Spec['geo_spec']
    ao_spec = Spec['ao_spec']
    mo_spec = Spec['mo_spec']
    is_vector = Spec['is_vector']
    drv = Spec['Derivative']
    calc_mo = Spec['calc_mo']
    
    if is_vector:
      # Set up Grid
      x = grid.x[xx[0]:xx[1]]
      y = grid.y[xx[0]:xx[1]]
      z = grid.z[xx[0]:xx[1]]
      N = (len(x),)
    else:
      x = xx
      y = grid.y
      z = grid.z
      N = (len(x),len(y),len(z))
    
    if drv is not None and calc_mo:
      delta_mo_list = []
      for ii_d in drv:
        # Calculate the derivatives of the AOs and MOs for this slice 
        delta_ao_list = ao_creator(geo_spec,ao_spec,drv=ii_d,
                    x=x,y=y,z=z,N=N,is_vector=is_vector)
        delta_mo_list.append(mo_creator(delta_ao_list,mo_spec,
                    x=x,y=y,z=z,N=N,is_vector=is_vector))
      return numpy.array(delta_mo_list)
    
    # Calculate the MOs and AOs for this slice 
    ao_list = ao_creator(geo_spec,ao_spec,x=x,y=y,z=z,N=N,is_vector=is_vector)
    mo_list = mo_creator(ao_list,mo_spec,x=x,y=y,z=z,N=N,is_vector=is_vector)
    
    if calc_mo:
      return numpy.array(mo_list)
    
    # Initialize a numpy array for the density 
    rho = numpy.zeros(N)
    
    # Initialize a numpy array for the norm of the MOs 
    mo_norm = numpy.zeros((len(mo_list)))
    
    # Calculate the density and the norm 
    for ii_mo in range(len(mo_list)): 
      mo_norm[ii_mo] = numpy.sum(numpy.square(mo_list[ii_mo]))
      rho += mo_spec[ii_mo]['occ_num'] * numpy.square(numpy.abs(mo_list[ii_mo]))
      
    if drv is None:
      # Return the density and the norm 
      return rho, mo_norm
    else:
      # Initialize a numpy array for the derivative of the density 
      delta_rho = numpy.zeros((len(drv),) + N)
      for i,ii_d in enumerate(drv):
        # Calculate the derivatives of the AOs and MOs for this slice 
        delta_ao_list = ao_creator(geo_spec,ao_spec,drv=ii_d,
                    x=x,y=y,z=z,N=N,is_vector=is_vector)
        delta_mo_list = mo_creator(delta_ao_list,mo_spec,
                    x=x,y=y,z=z,N=N,is_vector=is_vector)
        
        # Calculate the derivative of the density
        for ii_mo in range(len(mo_list)): 
          delta_rho[i] += (mo_spec[ii_mo]['occ_num'] * 
                    2 * delta_mo_list[ii_mo]*mo_list[ii_mo])
      # Return the derivative of the density 
      return rho, mo_norm, delta_rho
  except KeyboardInterrupt:
    # Catch keybord interrupt signal to prevent a hangup of the worker processes 
    return 0
  # slice_rho 

def rho_compute(qc,calc_mo=False,vector=None,drv=None,numproc=1):
  '''Calculate the density, the molecular orbitals, or the derivatives thereof.
  
  orbkit divides 3-dimensional regular grids into 2-dimensional slices and 
  1-dimensional vector grids into 1-dimensional slices of equal length. By default,
  3-dimensional grids are used (:literal:`vector=None`).
  The computational tasks are distributed to the worker processes.
  
  **Parameters:**
  
  qc : class or dict
    QCinfo class or dictionary containing the following attributes/keys.
    See `Central Variables`_ for details.
  qc.geo_spec : array_like, shape=(3,NATOMS) 
    See `Central Variables`_ for details.
  qc.ao_spec : List of dictionaries
    See `Central Variables`_ for details.
  qc.mo_spec : List of dictionaries
    See `Central Variables`_ for details.
  calc_mo : bool, optional
    If True, the computation of  the molecular orbitals requested is only
    carried out.
  vector : None or int, optional
    If not None, performs the computations on a vectorized grid, i.e., 
    with x, y, and z as vectors.
  drv : string or list of strings {None,'x','y', or 'z'}, optional
    If not None, computes the analytical derivative of the requested 
    quantities with respect to DRV.
  numproc : int
    Specifies number of subprocesses for multiprocessing.
  grid : module or class, global
    Contains the grid, i.e., grid.x, grid.y, and grid.z.

  **Returns:**
  
  :if calc_mo and drv is None: 
    - mo_list
  :if calc_mo and drv is not None:  
    - delta_mo_list
  :if not calc_mo and drv is None: 
    - rho
  :if not calc_mo and drv is not None: 
    - rho, delta_rho
  
  mo_list : numpy.ndarray, shape=((NMO,) + N)
    Contains the NMO=len(qc.mo_spec) molecular orbitals on a grid.
  delta_mo_list : numpy.ndarray, shape=((NDRV,NMO) + N)
    Contains the derivatives with respect to drv (NDRV=len(drv)) of the 
    NMO=len(qc.mo_spec) molecular orbitals on a grid.
  mo_norm : numpy.ndarray, shape=(NMO,)
    Contains the numerical norms of the molecular orbitals.
  rho : numpy.ndarray, shape=(N)
    Contains the density on a grid.
  delta_rho : numpy.ndarray, shape=((NDRV,) + N)
    Contains derivatives with respect to drv (NDRV=len(drv)) of 
    the density on a grid.
  '''
  
  if drv is not None:
    is_drv = True
    try:
      drv = list(drv)
    except TypeError: 
      drv = [drv]
  else:
    is_drv = False
  
  # Specify the global variable containing all desired information needed 
  # by the function slice_rho 
  global Spec
  
  if isinstance(qc,dict):
    Spec = qc
  else:
    Spec = qc.todict()
  Spec['calc_mo'] = calc_mo
  Spec['Derivative'] = drv
  Spec['is_vector'] = (vector is not None)
  
  mo_num = len(Spec['mo_spec'])
  
  if vector is None:
    is_vector = False
    N = tuple(grid.N_)
    sDim = 0
    sNum = N[sDim]
  else:
    is_vector = True
    if len(grid.x) != len(grid.y) or len(grid.x) != len(grid.z):
      raise ValueError("Dimensions of x-, y-, and z- coordinate differ!")
    N = (len(grid.x),)
    sNum = int(numpy.floor(N[0]/vector)+1)
  
  # The number of worker processes is capped to the number of 
  # grid points in x-direction.  
  if numproc > sNum: numproc = sNum
  
  # Print information regarding the density calculation 
  if not calc_mo:
    display("\nStarting the density calculation...")
  display("The grid has been separated into %d %sslices and the" % 
                (sNum, '2d-' if vector is None else ''))
  if numproc == 1:
    display("calculation will be carried out with 1 subprocess.\n" + 
    "\n\tThe number of subprocesses can be changed with -p\n")
  else:
    display("calculation will be carried out with %d subprocesses." 
            % numproc)
  display("\nThere are %d contracted AOs and %d MOs to be calculated."
            % (len(Spec['mo_spec'][0]['coeffs']), mo_num))
  
  # Initialize some additional user information 
  status_old = 0
  s_old = 0
  t = [time.time()]
  
  # Make slices 
  # Initialize an array to store the results 
  mo_norm = numpy.zeros((mo_num,))
  if calc_mo:
    mo_list = numpy.zeros(((mo_num,) if drv is None 
            else (len(drv),mo_num)) + tuple(N))
  else:
    rho = numpy.zeros(N)
    if is_drv:
      delta_rho = numpy.zeros((len(drv),) + N)

  # Start the worker processes
  if numproc > 0:
    pool = Pool(processes=numproc)
  
  # Write the slices in x to an array xx 
  xx = []
  i = 0
  for s in range(sNum):
    if not is_vector:
      xx.append((numpy.array([grid.x[s]])))
    else:
      if i == N[0]:
        sNum -= 1
        break
      elif (i + vector) >= N[0]:
        xx.append((numpy.array([i,N[0]],dtype=int)))      
      else:
        xx.append((numpy.array([i,i + vector],dtype=int)))
      i += vector 
  
  # Compute the density slice by slice   
  for s in range(sNum):
    # Which slice do we compute 
    if not is_vector:
      i = s
      j = s+1
    else:
      i = xx[s][0]
      j = xx[s][1]
    
    if (s == 0):
      # Check if the code for the atomic orbitals is compiled already
      if not is_compiled(ao_code(is_vector=is_vector,is_drv=is_drv)):
        display(('-'*80) + '''
        The C++ atomic orbital code seems not to be compiled yet.
        Running the first slice on one CPU...\n''' + ('-'*80))
        result = slice_rho(xx[s])
        if numproc > 0:
          it = pool.imap(slice_rho, xx[1:])
        display(('-'*80))
      else:
        if numproc > 0:
          it = pool.imap(slice_rho, xx)
        result = it.next() if numproc > 0 else slice_rho(xx[s])
    else:
      # Perform the compution for the current slice 
      result = it.next() if numproc > 0 else slice_rho(xx[s])
    # What output do we expect 
    if calc_mo:
      if not is_drv:
        mo_list[:,i:j] = result[:,:]
      else:
        for ii_d in range(len(drv)):
          mo_list[ii_d,:,i:j] = result[ii_d,:,:,]
    else:
      rho[i:j] = result[0]
      mo_norm += result[1]
      if is_drv:
        for ii_d in range(len(drv)):
          delta_rho[ii_d,i:j] = result[2][ii_d,:]
    
    
    # Print out the progress of the computation 
    status = round(s*100/float(sNum))
    if not status % 20 and status != status_old: 
      t.append(time.time())
      display("\tFinished %(f)d%% (%(s)d slices in %(t).3fs)" 
                % {'f': status,
                's': s + 1 - s_old,
                't': t[-1]-t[-2]})
      status_old = status
      s_old = s + 1
  
  # Close the worker processes
  if numproc > 0:
    pool.close()
    pool.join()
  
  if not is_vector:
    # Print the norm of the MOs 
    display("\nNorm of the MOs:")
    for ii_mo in range(len(mo_norm)):
      if calc_mo:
        norm = numpy.sum(numpy.square(mo_list[ii_mo]))*grid.d3r
      else:
        norm = mo_norm[ii_mo]*grid.d3r
      display("\t%(m).6f\tMO %(n)s" 
                % {'m':norm, 'n':Spec['mo_spec'][ii_mo]['sym']})
  
  if calc_mo:
    return mo_list
  
  if not is_vector:
    # Print the number of electrons 
    display("We have " + str(numpy.sum(rho)*grid.d3r) + " electrons.")
  if not is_drv:
    return rho
  else:
    return rho, delta_rho
  # rho_compute 

def rho_compute_no_slice(qc,calc_mo=False,is_vector=False,drv=None,
                         return_components=False,x=None,y=None,z=None):
  '''Calculates the density, the molecular orbitals, or the derivatives thereof
  without slicing the grid.
  
  **Parameters:**
  
  qc : class or dict
    QCinfo class or dictionary containing the following attributes/keys.
    See `Central Variables`_ for details.
  qc.geo_spec : array_like, shape=(3,NATOMS) 
    See `Central Variables`_ for details.
  qc.ao_spec : List of dictionaries
    See `Central Variables`_ for details.
  qc.mo_spec : List of dictionaries
    See `Central Variables`_ for details.
  calc_mo : bool, optional
    If True, the computation of  the molecular orbitals requested is only
    carried out.
  is_vector : bool, optional
    If True, performs the computations for a vectorized grid, i.e., 
    with x, y, and z as vectors.
  drv : string or list of strings {None,'x','y', or 'z'}, optional
    If not None, computes the analytical derivative of the requested 
    quantities with respect to DRV.
  return_components : bool, optional
    If True, returns the atomic and molecular orbitals, and the density, 
    and if requested, the derivatives thereof as well.
  grid : module or class, global
    Contains the grid, i.e., grid.x, grid.y, and grid.z.

  **Returns:**
  
  :if not return_components:
  
    :if calc_mo and drv is None: 
      - mo_list
    :if calc_mo and drv is not None:
      - delta_mo_list
    :if not calc_mo and drv is None:
      - rho
    :if not calc_mo and drv is not None: 
      - rho, delta_rho
  
  :else:
    | 
    
    :if calc_mo and drv is None:
      - ao_list,mo_list
    :if calc_mo and drv is not None: 
      - delta_ao_list,delta_mo_list
    :if not calc_mo and drv is None: 
      - ao_list,mo_list,rho
    :if not calc_mo and drv is not None: 
      - ao_list, mo_list, rho, delta_ao_list, delta_mo_list, delta_rho
  
  ao_list : numpy.ndarray, shape=((NAO,) + N)
    Contains the NAO=len(ao_spec) atomic orbitals on a grid.
  delta_ao_list : numpy.ndarray, shape=((NDRV,NAO) + N)
    Contains the derivatives with respect to drv (NDRV=len(drv)) of the 
    NAO=len(ao_spec) atomic orbitals on a grid.
  mo_list : numpy.ndarray, shape=((NMO,) + N)
    Contains the NMO=len(qc.mo_spec) molecular orbitals on a grid.
  delta_mo_list : numpy.ndarray, shape=((NDRV,NMO) + N)
    Contains the derivatives with respect to drv (NDRV=len(drv)) of the 
    NMO=len(qc.mo_spec) molecular orbitals on a grid.
  mo_norm : numpy.ndarray, shape=(NMO,)
    Contains the numerical norms of the molecular orbitals.
  rho : numpy.ndarray, shape=(N)
    Contains the density on a grid.
  delta_rho : numpy.ndarray, shape=((NDRV,) + N)
    Contains the derivatives with respect to drv (NDRV=len(drv)) of 
    the density on a grid.
  '''
  # Create the grid
  if x is None: x = grid.x
  if y is None: y = grid.y
  if z is None: z = grid.z
  
  if not isinstance(qc,dict):
    qc = qc.todict()
  
  #FIXME inaccurate implementation
  geo_spec = qc['geo_spec']
  ao_spec = qc['ao_spec']
  mo_spec = qc['mo_spec']
    
  if not is_vector:
    N = (len(x),len(y),len(z))
  else:
    if len(x) != len(y) or len(x) != len(z):
      display("Dimensions of x-, y-, and z- coordinate differ!")
      return 0
    else:
      N = (len(x),)
  
  if drv is not None:
    try:
      drv = list(drv)
    except TypeError: 
      drv = [drv]
    delta_mo_list = []
    for ii_d in drv:
      # Calculate the derivatives of the AOs and MOs for this slice 
      delta_ao_list = ao_creator(geo_spec,ao_spec,drv=ii_d,
                                 is_vector=is_vector,
                                 x=x,y=y,z=z)
      delta_mo_list.append(mo_creator(delta_ao_list,mo_spec,
                                      is_vector=is_vector,
                                      x=x,y=y,z=z))
    delta_mo_list = numpy.array(delta_mo_list)
    if calc_mo:
      return ((delta_ao_list, delta_mo_list) if return_components 
        else delta_mo_list)
  
  # Calculate the AOs and MOs 
  ao_list = ao_creator(geo_spec,ao_spec,is_vector=is_vector,x=x,y=y,z=z)
  mo_list = mo_creator(ao_list,mo_spec,is_vector=is_vector,x=x,y=y,z=z)
  
  if not is_vector:
    d3r = numpy.product([x[1]-x[0],y[1]-y[0],z[1]-z[0]])
    # Print the norm of the MOs 
    display("Norm of the MOs:")
    for ii_mo in range(len(mo_list)): 
      display("\t%(m).6f\tMO %(n)s" 
       % {'m':numpy.sum(mo_list[ii_mo]**2)*d3r, 'n':mo_spec[ii_mo]['sym']})
  
  if calc_mo:
    return ((ao_list, mo_list) if return_components 
      else mo_list)
  
  # Initialize a numpy array for the density 
  rho = numpy.zeros(N)
  
  # Calculate the density 
  for ii_mo in range(len(mo_list)): 
    rho += numpy.square(numpy.abs(mo_list[ii_mo])) * mo_spec[ii_mo]['occ_num']
  
  if not is_vector:
    d3r = numpy.product([x[1]-x[0],y[1]-y[0],z[1]-z[0]])    
    # Print the number of electrons 
    display("We have " + str(numpy.sum(rho)*d3r) + " electrons.")
  
  if drv is None:
    return ((ao_list, mo_list, rho) if return_components else rho)
  
  # Print information 
  display('\nCalculating the derivative of the density...')
  delta_rho = numpy.zeros((len(drv),) + N)
  
  # Loop over spatial directions 
  for ii_d in range(len(drv)):
    display('\t...with respect to %s' % drv[ii_d])
    # Calculate the derivative of the density
    for ii_mo in range(len(mo_list)): 
      delta_rho[ii_d] += (mo_spec[ii_mo]['occ_num'] * 
            2 * delta_mo_list[ii_d,ii_mo]*mo_list[ii_mo])
  
  return ((ao_list, mo_list, rho, delta_ao_list, delta_mo_list, delta_rho) 
        if return_components else (rho, delta_rho))
  # rho_compute_no_slice 


#--- Support Code ---# 

# Assign the quantum number l to every AO symbol (s,p,d,etc.) 
orbit = 'spd' + string.lowercase[5:].replace('s','').replace('p','')
lquant = dict([(j, i) for i,j in enumerate(orbit)])
del i,j

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

def l_deg(l=0,ao=None):
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
  return (l+1)*(l+2)/2
  # l_deg 

def integration(matrix,x=None,y=None,z=None):
  from scipy import integrate
  
  if x is None: x = grid.x
  if y is None: y = grid.y
  if z is None: z = grid.z  
  
  if matrix.squeeze().ndim == 3:
    integral = integrate.simps(matrix, x, axis=0, even='avg')
    integral = integrate.simps(integral, y, axis=0, even='avg')
    integral = integrate.simps(integral, z, axis=0, even='avg')
  elif matrix.squeeze().ndim == 2:
    if len(x) == 1:
      r = y
      matrix = matrix[0,:,:]
    elif len(y) == 1:
      r = x
      matrix = matrix[:,0,:]
    else:
      print('dim(z) = 1! No cylindrical coordinates...')
      return 255
    [Z,R] = numpy.meshgrid(z,r)
    integral = 2*numpy.pi*integrate.simps(R*matrix, r, axis=0, even='avg')
    integral = integrate.simps(integral, z, axis=0, even='avg')
  else: 
    #print('Error in integration! ndim is not 3 or 2...')
    return numpy.sum(matrix)
  return integral

# C++ Code 
def ao_code(is_vector=False,is_drv=False):
  '''Returns the requested C++ code.
  
  **Parameters:**
  
  is_vector : bool
    If True, returns the code for the computation of an atomic
    orbital on a vecotrized grid_file.
  is_drv : bool
    If True, returns the code for the computation of the
    derivative of an atomic orbital.
  '''
  if not is_vector and not is_drv:
    code =  '''
    double Norm[ao_num][Nao_list[0]];
    double X, Y, Z;
    int lx[Nao_list[0]], ly[Nao_list[0]], lz[Nao_list[0]];
    double rr, ao_l0[ao_num], ao_xyz[Nao_list[0]];
    

    for (int il=0; il<Nao_list[0]; il++)
    {
      lx[il] = EXP_LIST2(il,0);
      ly[il] = EXP_LIST2(il,1);
      lz[il] = EXP_LIST2(il,2);
      
      for (int ii=0; ii<ao_num; ii++)
      {
        Norm[ii][il] = ao_norm(lx[il],ly[il],lz[il],&COEFF_LIST2(ii,0));
      }
    }
    
    for (int i=0; i<Nx[0]; i++)
    {
      for (int j=0; j<Ny[0]; j++)
      {
        for (int k=0; k<Nz[0]; k++)
        {
        X = x[i]-at_pos[0];
        Y = y[j]-at_pos[1];
        Z = z[k]-at_pos[2];
        
        rr = pow(X,2)+pow(Y,2)+pow(Z,2);
        
        for (int il=0; il<Nao_list[0]; il++)
        {
            ao_xyz[il] = xyz(X, Y, Z, lx[il], ly[il], lz[il]);
        }
        
        for (int ii=0; ii<ao_num; ii++)
        {
            ao_l0[ii] = COEFF_LIST2(ii,1) * exp(-COEFF_LIST2(ii,0) * rr);
            
            for (int il=0; il<Nao_list[0]; il++)
            {
            AO_LIST4(il,i,j,k) += Norm[ii][il] * ao_xyz[il] * ao_l0[ii];
            }
          }
        }
      }
    }
    '''
  elif is_vector and not is_drv:
    code =  '''
    double Norm[ao_num][Nao_list[0]];
    double X, Y, Z;
    int lx[Nao_list[0]], ly[Nao_list[0]], lz[Nao_list[0]];
    double rr, ao_l0[ao_num], ao_xyz[Nao_list[0]];
    

    for (int il=0; il<Nao_list[0]; il++)
    {
      lx[il] = EXP_LIST2(il,0);
      ly[il] = EXP_LIST2(il,1);
      lz[il] = EXP_LIST2(il,2);
      
      for (int ii=0; ii<ao_num; ii++)
      {
        Norm[ii][il] = ao_norm(lx[il],ly[il],lz[il],&COEFF_LIST2(ii,0));
      }
    }
    
    for (int i=0; i<Nx[0]; i++)
    {
      X = x[i]-at_pos[0];
      Y = y[i]-at_pos[1];
      Z = z[i]-at_pos[2];
      
      rr = pow(X,2)+pow(Y,2)+pow(Z,2);
      
      for (int il=0; il<Nao_list[0]; il++)
      { 
        ao_xyz[il] = xyz(X, Y, Z, lx[il], ly[il], lz[il]);
      }
      
      for (int ii=0; ii<ao_num; ii++)
      {
        ao_l0[ii] = COEFF_LIST2(ii,1) * exp(-COEFF_LIST2(ii,0) * rr);
        
        for (int il=0; il<Nao_list[0]; il++)
        {
        AO_LIST2(il,i) += Norm[ii][il] * ao_xyz[il] * ao_l0[ii];
        }
      }
    }
    '''
  elif not is_vector and is_drv:
    code = '''
    double Norm[ao_num][Nao_list[0]];
    double X, Y, Z;
    int lx[Nao_list[0]], ly[Nao_list[0]], lz[Nao_list[0]];
    double rr, ao_l0[ao_num], ao_xyz;
    

    for (int il=0; il<Nao_list[0]; il++)
    {
      lx[il] = EXP_LIST2(il,0);
      ly[il] = EXP_LIST2(il,1);
      lz[il] = EXP_LIST2(il,2);
      
      for (int ii=0; ii<ao_num; ii++)
      {
        Norm[ii][il] = ao_norm(lx[il],ly[il],lz[il],&COEFF_LIST2(ii,0));
      }
    }
    
    for (int i=0; i<Nao_list[1]; i++)
    {
      for (int j=0; j<Nao_list[2]; j++)
      {
        for (int k=0; k<Nao_list[3]; k++)
        {
        X = x[i]-at_pos[0];
        Y = y[j]-at_pos[1];
        Z = z[k]-at_pos[2];
        
        rr = pow(X,2)+pow(Y,2)+pow(Z,2);
        
        for (int ii=0; ii<ao_num; ii++)
        {
            ao_l0[ii] = COEFF_LIST2(ii,1) * exp(-COEFF_LIST2(ii,0) * rr);
            
            for (int il=0; il<Nao_list[0]; il++)
            {
              switch(drv)
              {
              case 0:
              {
                if (lx[il] == 0)
                {
                    ao_xyz = - 2 * COEFF_LIST2(ii,0) * 
                        xyz(X, Y, Z, lx[il]+1, ly[il], lz[il]);
                }
                else
                {
                    ao_xyz = lx[il] * xyz(X, Y, Z, lx[il]-1, ly[il], lz[il]) - 
                        2 * COEFF_LIST2(ii,0) * 
                        xyz(X, Y, Z, lx[il]+1, ly[il], lz[il]);    
                }
              } break;
              
              case 1:
              {
              if (ly[il] == 0)
                {
                    ao_xyz = - 2 * COEFF_LIST2(ii,0) * 
                        xyz(X, Y, Z, lx[il], ly[il]+1, lz[il]);
                }
                else
                {
                    ao_xyz = ly[il] * xyz(X, Y, Z, lx[il], ly[il]-1, lz[il]) 
                        - 2 * COEFF_LIST2(ii,0) * 
                        xyz(X, Y, Z, lx[il], ly[il]+1, lz[il]);   
                }
              } break;
              
              case 2:
              {
                if (lz[il] == 0)
                {
                    ao_xyz = - 2 * COEFF_LIST2(ii,0) *
                        xyz(X, Y, Z, lx[il], ly[il], lz[il]+1);
                }
                else
                {
                    ao_xyz = lz[il] * xyz(X, Y, Z, lx[il], ly[il], lz[il]-1) 
                            - 2 * COEFF_LIST2(ii,0) * 
                            xyz(X, Y, Z, lx[il], ly[il], lz[il]+1);    
                }
              } break;
              
              default:
              {
                std::cout << "False statement for derivative variable!" 
                          << std::endl;
              }
            }
            AO_LIST4(il,i,j,k) += Norm[ii][il] * ao_xyz * ao_l0[ii];
            }
          }
        }
      }
    }
    
    '''
  else:
    code = '''
    double Norm[ao_num][Nao_list[0]];
    double X, Y, Z;
    int lx[Nao_list[0]], ly[Nao_list[0]], lz[Nao_list[0]];
    double rr, ao_l0[ao_num], ao_xyz;
    

    for (int il=0; il<Nao_list[0]; il++)
    {
      lx[il] = EXP_LIST2(il,0);
      ly[il] = EXP_LIST2(il,1);
      lz[il] = EXP_LIST2(il,2);
      
      for (int ii=0; ii<ao_num; ii++)
      {
        Norm[ii][il] = ao_norm(lx[il],ly[il],lz[il],&COEFF_LIST2(ii,0));
      }
    }
    
    for (int i=0; i<Nx[0]; i++)
    {
      X = x[i]-at_pos[0];
      Y = y[i]-at_pos[1];
      Z = z[i]-at_pos[2];
      
      rr = pow(X,2)+pow(Y,2)+pow(Z,2);
      
      for (int ii=0; ii<ao_num; ii++)
      {
        ao_l0[ii] = COEFF_LIST2(ii,1) * exp(-COEFF_LIST2(ii,0) * rr);
        
        for (int il=0; il<Nao_list[0]; il++)
        {
        switch(drv)
        {
            case 0:
            {
            if (lx[il] == 0)
            {
            ao_xyz = - 2 * COEFF_LIST2(ii,0) * 
                xyz(X, Y, Z, lx[il]+1, ly[il], lz[il]);
            }
            else
            {
            ao_xyz = lx[il] * xyz(X, Y, Z, lx[il]-1, ly[il], lz[il]) 
                    - 2 * COEFF_LIST2(ii,0) * 
                    xyz(X, Y, Z, lx[il]+1, ly[il], lz[il]);     
          }
        } break;
        
            case 1:
            {
            if (ly[il] == 0)
            {
            ao_xyz = - 2 * COEFF_LIST2(ii,0) * 
                    xyz(X, Y, Z, lx[il], ly[il]+1, lz[il]);
            }
            else
            {
            ao_xyz = ly[il] * xyz(X, Y, Z, lx[il], ly[il]-1, lz[il]) 
                    - 2 * COEFF_LIST2(ii,0) * 
                    xyz(X, Y, Z, lx[il], ly[il]+1, lz[il]);    
            }
            } break;
            
            case 2:
            {
            if (lz[il] == 0)
            {
            ao_xyz = - 2 * COEFF_LIST2(ii,0) * 
                    xyz(X, Y, Z, lx[il], ly[il], lz[il]+1);
            }
            else
            {
            ao_xyz = lz[il] * xyz(X, Y, Z, lx[il], ly[il], lz[il]-1) 
                    - 2 * COEFF_LIST2(ii,0) * 
                    xyz(X, Y, Z, lx[il], ly[il], lz[il]+1);    
            }
            } break;
            
            default:
            {
            std::cout << "False statement for derivative variable!" 
                        << std::endl;
            }
        }
        AO_LIST2(il,i) += Norm[ii][il] * ao_xyz * ao_l0[ii];
        }
      }
    }
    '''
  return code

def is_compiled(code):
  '''Checks if the C++ code is already compiled.
  
  Adaped from :func:`weave.inline_tools`.'''
  # 1. try local cache
  try:
    weave.inline_tools.function_cache[code]
    return True
  except KeyError:
    pass
  # 2. try catalog cache
  if weave.inline_tools.function_catalog.get_functions_fast(code) != []:
    return True
  # 3. try persistent catalog
  if weave.inline_tools.function_catalog.get_functions(code) != []:
    return True
  # The funciton was not found
  return False
    
def slicer(N,vector=1e4,numproc=1):
  i = 0
  sNum = int((N/(vector))+1) if int(vector) > 0 else int(N)
  xx = []
  if numproc > 1:
    for s in range(sNum):
      if i == N:
        N -= 1
        break
      elif (i + vector) >= N:
        xx.append((numpy.array([i,N],dtype=int)))      
      else:
        xx.append((numpy.array([i,i + vector],dtype=int)))
      i += vector
  else:
    xx.append((numpy.array([0,N],dtype=int))) 
  return xx
