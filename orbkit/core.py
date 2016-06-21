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
from orbkit import grid,cy_grid,cy_core
from orbkit.display import display

def ao_creator(geo_spec,ao_spec,ao_spherical=None,drv=None,
               x=None,y=None,z=None,is_vector=None):
  '''Calculates all contracted atomic orbitals or its
  derivatives with respect to a specific variable (e.g. drv = 'x' or drv = 0).
  
  **Parameters:**
  
  geo_spec,ao_spec :
    See :ref:`Central Variables` in the manual for details.
  sel_ao : int
    Index of the requested atomic orbital
  drv : int or string, {None, 'x', 'y', 'z', 0, 1, 2}, optional
    If not None, an analytical  calculation of the derivatives for 
    the atomic orbitals with respect to DRV is requested.
  x,y,z : None or list of floats, optional
    If not None, provides a list of Cartesian coordinates, 
    else the respective coordinates of grid. will be used
  is_vector : bool, optional
    If True, a vector grid will be applied
  
  **Returns:**
  
  ao_list : numpy.ndarray, shape=((NAO,) + N)
    Contains the computed NAO atomic orbitals on a grid.
  '''
  # Create the grid
  if all(v is None for v in [x,y,z,is_vector]) and not grid.is_initialized:
    display('\nSetting up the grid...')
    grid.grid_init(is_vector=True)
    display(grid.get_grid())   # Display the grid    
  if x is None: x = grid.x
  if y is None: y = grid.y
  if z is None: z = grid.z
  if is_vector is None: is_vector = grid.is_vector
  
  was_vector = is_vector  
  if not is_vector:
    N = (len(x),len(y),len(z))
    # Convert regular grid to vector grid
    x,y,z = cy_grid.grid2vector(x.copy(),y.copy(),z.copy())
  else:
    if len(x) != len(y) or len(x) != len(z):
      raise ValueError("Dimensions of x-, y-, and z- coordinate differ!")
    N = (len(x),)

  x = require(x, dtype='f')
  y = require(y, dtype='f')
  z = require(z, dtype='f') 
  
  geo_spec = require(geo_spec, dtype='f')
  lxlylz,assign = get_lxlylz(ao_spec,get_assign=True,bincount=True)
  ao_coeffs,pnum_list,atom_indices = prepare_ao_calc(ao_spec)
  
  is_normalized = each_ao_is_normalized(ao_spec)
  drv = validate_drv(drv)
  
  lxlylz = require(lxlylz,dtype='i')
  assign = require(assign,dtype='i')
  ao_list = cy_core.aocreator(lxlylz,assign,ao_coeffs,pnum_list,geo_spec,
                              atom_indices,x,y,z,drv,is_normalized)
  
  if not (ao_spherical is None or ao_spherical == []):
    ao_list = cartesian2spherical(ao_list,ao_spec,ao_spherical)
  
  if not was_vector: return ao_list.reshape((len(ao_list),) + N,order='C')
  return ao_list
 
def mo_creator(ao_list,mo_spec):
  '''Calculates the molecular orbitals.
  
  **Parameters:**
  
  ao_list : numpy.ndarray, shape=((NAO,) + N)
    Contains the NAO atomic orbitals on a grid.
  mo_spec : List of dictionaries
    See :ref:`Central Variables` for details.
  mo_coeff : numpy.ndarray, shape = (NMO,NAO)
    Contains the molecular orbital coefficients of all orbitals.
    
  **Returns:**
  
  mo_list : numpy.ndarray, shape=((NMO,) + N)
    Contains the NMO=len(mo_spec) molecular orbitals on a grid.
  '''
  ao_list = require(ao_list,dtype='f')
  shape = ao_list.shape
  ao_list.shape = (shape[0],-1)
  mo_coeff = create_mo_coeff(mo_spec,name='The argument `mo_spec`')
  mo_list = cy_core.mocreator(ao_list,mo_coeff).reshape(
                                                 ((len(mo_coeff),) + shape[1:]),
                                                 order='C')
  ao_list.shape = shape
  return mo_list


def cartesian2spherical(ao_list,ao_spec,ao_spherical):
  '''Transforms the atomic orbitals from a Cartesian Gaussian basis to a 
  (real) pure spherical harmonic Gaussian basis set.
  
  Adapted from H.B. Schlegel and M.J. Frisch,
  International Journal of Quantum Chemistry, Vol. 54, 83-87 (1995).
  
  **Parameters:**
  
  ao_list : numpy.ndarray, shape=((NAO,) + N)
    Contains the NAO atomic orbitals on a grid.  
  ao_spec,ao_spherical :
    See :ref:`Central Variables` in the manual for details.
  
  **Returns:**
  
  ao_list_sph : numpy.ndarray, shape=((NAO,) + N)
    Contains the NAO spherical atomic orbitals on a grid. 
  
  ..hint: 
  
    The conversion is currently only supported up to g atomic orbitals.
  '''
  lxlylz,assign = get_lxlylz(ao_spec,get_assign=True)

  l = [[] for i in ao_spec]
  for i,j in enumerate(assign):
    l[j].append(i) 
  
  shape = list(ao_list.shape)
  shape[0] = len(ao_spherical)
  ao_list_sph = numpy.zeros(shape)
  for i0,(j0,k0) in enumerate(ao_spherical):
    sph0 = get_cart2sph(*k0)
    for c0 in range(len(sph0[0])):
      for i,j in enumerate(l[j0]):
        if tuple(lxlylz[j]) == sph0[0][c0]:
          index0 = i + l[j0][0]
      ao_list_sph[i0,:] += sph0[1][c0]*sph0[2]*ao_list[index0,:]
  
  return ao_list_sph

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
      :geo_spec: List of floats, shape=(NATOMS, 3) (see :ref:`Central Variables` for details).
      :ao_spec: List of dictionaries (see :ref:`Central Variables` for details).
      :mo_spec: List of dictionaries (see :ref:`Central Variables` for details).
      :calc_mo: Bool if only the molecular orbitals are requested.
      :is_vector: Bool if a vector grid is used.
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
    ao_spherical = Spec['ao_spherical']
    mo_spec = Spec['mo_spec']
    drv = Spec['Derivative']
    calc_mo = Spec['calc_mo']
    
    # Set up Grid
    x = grid.x[xx[0]:xx[1]]
    y = grid.y[xx[0]:xx[1]]
    z = grid.z[xx[0]:xx[1]]
    N = (len(x),)
    
    if drv is not None and calc_mo:
      delta_mo_list = []
      for ii_d in drv:
        # Calculate the derivatives of the AOs and MOs for this slice 
        delta_ao_list = ao_creator(geo_spec,ao_spec,ao_spherical=ao_spherical,drv=ii_d,
                    x=x,y=y,z=z,is_vector=True)
        delta_mo_list.append(mo_creator(delta_ao_list,mo_spec))
      return numpy.array(delta_mo_list)
    # Calculate the MOs and AOs for this slice 
    ao_list = ao_creator(geo_spec,ao_spec,ao_spherical=ao_spherical,x=x,y=y,z=z,is_vector=True)
    mo_list = mo_creator(ao_list,mo_spec)
    
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
        delta_ao_list = ao_creator(geo_spec,ao_spec,ao_spherical=ao_spherical,drv=ii_d,
                                   x=x,y=y,z=z,is_vector=True)
        delta_mo_list = mo_creator(delta_ao_list,mo_spec)        
        if len(ii_d) == 2:
          ao_0 = ao_creator(geo_spec,ao_spec,ao_spherical=ao_spherical,drv=ii_d[0],
                                  x=x,y=y,z=z)
          if '2' in ii_d or ii_d[0] == ii_d[1]:          
            delta2_mo_list = numpy.array(mo_creator(ao_0,mo_spec))**2
          else:
            ao_1 = ao_creator(geo_spec,ao_spec,ao_spherical=ao_spherical,drv=ii_d[1],
                                  x=x,y=y,z=z,is_vector=True)
            delta2_mo_list = (numpy.array(mo_creator(ao_0,mo_spec)) *
                              numpy.array(mo_creator(ao_1,mo_spec)))
        # Calculate the derivative of the density
        for ii_mo in range(len(mo_list)): 
          delta_rho[i] += (mo_spec[ii_mo]['occ_num'] * 
                    2 * delta_mo_list[ii_mo]*mo_list[ii_mo])
          if len(ii_d) == 2:
            delta_rho[i] += mo_spec[ii_mo]['occ_num'] * 2 * delta2_mo_list[ii_mo]
      # Return the derivative of the density 
      return rho, mo_norm, delta_rho
  except KeyboardInterrupt:
    # Catch keybord interrupt signal to prevent a hangup of the worker processes 
    return 0
  # slice_rho 

def initializer(global_args):
  global Spec
  Spec = global_args
  
def rho_compute(qc,calc_mo=False,drv=None,laplacian=False,
                numproc=1,slice_length=1e4,vector=None,save_hdf5=False,
                **kwargs):
  r'''Calculate the density, the molecular orbitals, or the derivatives thereof.
  
  orbkit divides 3-dimensional regular grids into 2-dimensional slices and 
  1-dimensional vector grids into 1-dimensional slices of equal length. By default,
  3-dimensional grids are used (:literal:`vector=None`).
  The computational tasks are distributed to the worker processes.
  
  **Parameters:**
  
  qc : class or dict
    QCinfo class or dictionary containing the following attributes/keys.
    See :ref:`Central Variables` for details.
  qc.geo_spec : numpy.ndarray, shape=(3,NATOMS) 
    See :ref:`Central Variables` for details.
  qc.ao_spec : List of dictionaries
    See :ref:`Central Variables` for details.
  qc.mo_spec : List of dictionaries
    See :ref:`Central Variables` for details.
  calc_mo : bool, optional
    If True, the computation of  the molecular orbitals requested is only
    carried out.
  slice_length : int, optional
    If not None, performs the computations on a vector grid, i.e., 
    with x, y, and z as vectors.
  drv : string or list of strings {None,'x','y', 'z', 'xx', 'xy', ...}, optional
    If not None, computes the analytical derivative of the requested 
    quantities with respect to DRV.
  laplacian : bool, optional
    If True, computes the laplacian of the density.
  numproc : int
    Specifies number of subprocesses for multiprocessing.
  grid : module or class, global
    Contains the grid, i.e., grid.x, grid.y, and grid.z. If grid.is_initialized
    is not True, functions runs grid.grid_init().

  **Returns:**
  
  :if calc_mo and drv is None: 
    - mo_list
  :if calc_mo and drv is not None:  
    - delta_mo_list
  :if not calc_mo and drv is None: 
    - rho
  :if not calc_mo and drv is not None: 
    - rho, delta_rho
  :if not calc_mo and laplacian:
    - rho, delta_rho, laplacian_rho      
  
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
  laplacian_rho : numpy.ndarray, shape=(N)
    Contains the laplacian of the density on a grid, i.e. 
    :math:`\nabla^2 \rho = \nabla^2_x \rho + \nabla^2_y \rho + \nabla^2_z \rho`.
  '''  
  slice_length = slice_length if not vector else vector
  if slice_length <= 0:
    return rho_compute_no_slice(qc,calc_mo=calc_mo,drv=drv,
                                laplacian=laplacian,**kwargs)
  if laplacian:
    if not (drv is None or drv == ['xx','yy','zz'] or drv == ['x2','y2','z2']):
      display('Note: You have set the option `laplacian` and specified values\n' +
              'for `drv`. Both options are not compatible.\n' +
              'The option `drv` has been changed to `drv=["xx","yy","zz"]`.')
    drv = ['xx','yy','zz']
  
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
  if isinstance(qc,dict):
    Spec = qc
  else:
    Spec = qc.todict()
  Spec['calc_mo'] = calc_mo
  Spec['Derivative'] = drv
  mo_num = len(Spec['mo_spec'])
  
  if not grid.is_initialized:
    display('\nSetting up the grid...')
    grid.grid_init(is_vector=True)
    display(grid.get_grid())   # Display the grid
  
  was_vector = grid.is_vector
  N = (len(grid.x),) if was_vector else (len(grid.x),len(grid.y),len(grid.z))
  if not was_vector:
    grid.grid2vector()
    display('Converting the regular grid to a vector grid containing ' +
            '%.2e grid points...' % len(grid.x))
  
  # Define the slice length
  npts = len(grid.x)
  sNum = int(numpy.floor(npts/slice_length)+1)
  
  # The number of worker processes is capped to the number of 
  # grid points in x-direction.  
  if numproc > sNum: numproc = sNum
  
  # Print information regarding the density calculation 
  display("\nStarting the calculation of the %s..." % 
         ("molecular orbitals" if calc_mo else "density"))
  display("The grid has been separated into %d slices each having %.2e grid points." % 
                (sNum, slice_length))
  if numproc <= 1:
    display("The calculation will be carried out using only one process.\n" + 
    "\n\tThe number of subprocesses can be changed with -p\n")
  else:
    display("The calculation will be carried out with %d subprocesses." 
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
  
  def zeros(shape,name,save_hdf5):
    if not save_hdf5:
      return numpy.zeros(shape)
    else:
      return f.create_dataset(name,shape,dtype=numpy.float64,chunks=shape[:-1] + (slice_length,))
  def reshape(data,shape):
    if not save_hdf5:
      return data.reshape(shape)
    else:
      data.attrs['shape'] = shape
      return data[...].reshape(shape)
  
  if save_hdf5:
    import h5py
    f = h5py.File(str(save_hdf5), 'w')
    f['x'] = grid.x
    f['y'] = grid.y
    f['z'] = grid.z
  
  if calc_mo:
    mo_list = zeros((mo_num,npts) if drv is None 
                    else (len(drv),mo_num,npts),"mo_list",save_hdf5)
  else:
    rho = zeros(npts,"rho",save_hdf5)
    if is_drv:
      delta_rho = zeros((len(drv),npts),"rho",save_hdf5)
  
  # Write the slices in x to an array xx 
  xx = []
  i = 0
  for s in range(sNum):
    if i == npts:
      sNum -= 1
      break
    elif (i + slice_length) >= npts:
      xx.append((numpy.array([i,npts],dtype=int)))      
    else:
      xx.append((numpy.array([i,i + slice_length],dtype=int)))
    i += slice_length 
  
  # Start the worker processes
  if numproc > 1:
    pool = Pool(processes=numproc, initializer=initializer, initargs=(Spec,))
    it = pool.imap(slice_rho, xx)
  else:
    initializer(Spec)
  
  # Compute the density slice by slice   
  for s in range(sNum):
    # Which slice do we compute 
    i = xx[s][0]
    j = xx[s][1]    
    # Perform the compution for the current slice 
    result = it.next() if numproc > 1 else slice_rho(xx[s])
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
    status = numpy.floor(s*10/float(sNum))*10
    if not status % 10 and status != status_old:
      t.append(time.time())
      display("\tFinished %(f)d %% (%(s)d slices in %(t).3f s)" 
                % {'f': status,
                's': s + 1 - s_old,
                't': t[-1]-t[-2]})
      status_old = status
      s_old = s + 1
  
  # Close the worker processes
  if numproc > 1:
    pool.close()
    pool.join()
    
  if not was_vector:
    grid.vector2grid(*N)
    display('Converting the output from a vector grid to a regular grid...')
  
  if not was_vector and drv is None:
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
    #if not was_vector: 
    mo_list = reshape(mo_list,((mo_num,) if drv is None 
                                         else (len(drv),mo_num,)) + N)
    if save_hdf5: f.close()
    return mo_list
  
  if not was_vector:
    # Print the number of electrons
    display("We have " + str(numpy.sum(rho)*grid.d3r) + " electrons.")
  
  #if not was_vector: 
  rho = reshape(rho,N)
  if not is_drv:
    if save_hdf5: f.close()
    return rho
  else:
    #if not was_vector: 
    delta_rho = reshape(delta_rho,(len(drv),) + N)
    if save_hdf5: f.close()
    if laplacian: return rho, delta_rho, delta_rho.sum(axis=0)
    return rho, delta_rho
  # rho_compute 

def rho_compute_no_slice(qc,calc_mo=False,drv=None,
                         laplacian=False,return_components=False,
                         x=None,y=None,z=None,is_vector=None,**kwargs):
  r'''Calculates the density, the molecular orbitals, or the derivatives thereof
  without slicing the grid.
  
  **Parameters:**
  
  qc : class or dict
    QCinfo class or dictionary containing the following attributes/keys.
    See :ref:`Central Variables` for details.
  qc.geo_spec : numpy.ndarray, shape=(3,NATOMS) 
    See :ref:`Central Variables` for details.
  qc.ao_spec : List of dictionaries
    See :ref:`Central Variables` for details.
  qc.mo_spec : List of dictionaries
    See :ref:`Central Variables` for details.
  calc_mo : bool, optional
    If True, the computation of  the molecular orbitals requested is only
    carried out.
  is_vector : bool, optional
    If True, performs the computations for a vector grid, i.e., 
    with x, y, and z as vectors.
  drv : string or list of strings {None,'x','y', or 'z'}, optional
    If not None, computes the analytical derivative of the requested 
    quantities with respect to DRV.
  laplacian : bool, optional
    If True, computes the laplacian of the density.
  return_components : bool, optional
    If True, returns the atomic and molecular orbitals, and the density, 
    and if requested, the derivatives thereof as well.
  x,y,z : numpy.ndarray, optional
    If not None, provides a list of Cartesian coordinates, 
    else the respective coordinates of the module :mod:`orbkit.grid` will 
    be used.

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
    :if not calc_mo and laplacian:
      - rho, delta_rho, laplacian_rho      
  
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
    :if not calc_mo and laplacian: 
      - ao_list, mo_list, rho, delta_ao_list, delta_mo_list, delta_rho, laplacian_rho
  
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
  laplacian_rho : numpy.ndarray, shape=(N)
    Contains the laplacian of the density on a grid, i.e. 
    :math:`\nabla^2 \rho = \nabla^2_x \rho + \nabla^2_y \rho + \nabla^2_z \rho`.
  '''  
  
  # Create the grid
  if all(v is None for v in [x,y,z,is_vector]) and not grid.is_initialized:
    display('\nSetting up the grid...')
    grid.grid_init(is_vector=True)
    display(grid.get_grid())   # Display the grid    
  if x is None: x = grid.x
  if y is None: y = grid.y
  if z is None: z = grid.z
  if is_vector is None: is_vector = grid.is_vector
  
  def convert(data,was_vector,N):
    data = numpy.array(data,order='C')
    if not was_vector:
      data = data.reshape(data.shape[:-1] + N,order='C')
    return data
  
  was_vector = is_vector  
  if not is_vector:
    N = (len(x),len(y),len(z))
    d3r = numpy.product([x[1]-x[0],y[1]-y[0],z[1]-z[0]])
    # Convert regular grid to vector grid
    x,y,z = cy_grid.grid2vector(x.copy(),y.copy(),z.copy())
    is_vector = True
    display('Converting the regular grid to a vector grid containing ' +
            '%.2e grid points...' % len(grid.x))
  else:
    if len(x) != len(y) or len(x) != len(z):
      raise ValueError("Dimensions of x-, y-, and z- coordinate differ!")
    N = (len(x),)
  
  if not isinstance(qc,dict):
    qc = qc.todict()
  
  geo_spec = qc['geo_spec']
  ao_spec = qc['ao_spec']
  ao_spherical = qc['ao_spherical']
  mo_spec = qc['mo_spec']
    
  if laplacian:
    if not (drv is None or drv == ['xx','yy','zz'] or drv == ['x2','y2','z2']):
      display('Note: You have set the option `laplacian` and specified values\n' +
              'for `drv`. Both options are not compatible.\n' +
              'The option `drv` has been changed to `drv=["xx","yy","zz"]`.')
    drv = ['xx','yy','zz']
  
  display('\nStarting the calculation without slicing the grid...')
  
  if drv is not None:
    try:
      drv = list(drv)
    except TypeError: 
      drv = [drv]
    
    display('\nCalculating the derivatives of the atomic and molecular orbitals...')
    delta_ao_list = [[] for ii_d in drv]
    delta_mo_list = [[] for ii_d in drv]
    for i,ii_d in enumerate(drv):
      display('\t...with respect to %s' % ii_d)
      # Calculate the derivatives of the AOs and MOs
      delta_ao_list[i] = ao_creator(geo_spec,ao_spec,ao_spherical=ao_spherical,
                                      drv=ii_d,
                                      is_vector=True,x=x,y=y,z=z)
      delta_mo_list[i] =  mo_creator(delta_ao_list[i],mo_spec)                     
    delta_ao_list = convert(delta_ao_list,was_vector,N)
    delta_mo_list = convert(delta_mo_list,was_vector,N)
    if calc_mo:  
      return ((delta_ao_list,delta_mo_list) if return_components 
              else delta_mo_list)
    
    delta2_mo_list = [None for ii_d in drv]
    for i,ii_d in enumerate(drv):
      if len(ii_d) == 2:
        display('\t...with respect to %s' % ii_d[0])
        ao_0 = ao_creator(geo_spec,ao_spec,ao_spherical=ao_spherical,
                          drv=ii_d[0],
                          is_vector=True,x=x,y=y,z=z)
        if '2' in ii_d or ii_d == 'xx'  or ii_d == 'yy' or ii_d == 'zz':
          delta2_mo_list[i] = mo_creator(ao_0,mo_spec)**2
        else: 
          display('\t...with respect to %s' % ii_d[1])
          ao_1 = ao_creator(geo_spec,ao_spec,ao_spherical=ao_spherical,
                            drv=ii_d[1],
                            is_vector=True,x=x,y=y,z=z)
          delta2_mo_list[i] = (mo_creator(ao_0,mo_spec) *
                               mo_creator(ao_1,mo_spec))    
        delta2_mo_list[i] = convert(delta2_mo_list[i],was_vector,N)
  
  display('\nCalculating the atomic and molecular orbitals...')
  # Calculate the AOs and MOs 
  ao_list = ao_creator(geo_spec,ao_spec,ao_spherical=ao_spherical,
                       is_vector=True,x=x,y=y,z=z)
  mo_list = convert(mo_creator(ao_list,mo_spec),was_vector,N)
  ao_list = convert(ao_list,was_vector,N)
  if not was_vector:
    # Print the norm of the MOs 
    display("\nNorm of the MOs:")
    for ii_mo in range(len(mo_list)): 
      display("\t%(m).6f\tMO %(n)s" 
       % {'m':numpy.sum(mo_list[ii_mo]**2)*d3r, 'n':mo_spec[ii_mo]['sym']})
  
  if calc_mo:
    return ((ao_list, mo_list) if return_components 
            else mo_list)
  
  # Initialize a numpy array for the density 
  rho = numpy.zeros(N)
  
  display('\nCalculating the density...') 
  for ii_mo in range(len(mo_list)): 
    rho += numpy.square(numpy.abs(mo_list[ii_mo])) * mo_spec[ii_mo]['occ_num']
  
  if not was_vector:
    # Print the number of electrons 
    display("We have " + str(numpy.sum(rho)*d3r) + " electrons.")
  
  if drv is None:
    return ((ao_list, mo_list, rho) if return_components else rho)
  
  # Print information 
  display('\nCalculating the derivative of the density...')
  delta_rho = numpy.zeros((len(drv),) + N)
  # Loop over spatial directions 
  for i,ii_d in enumerate(drv):
    display('\t...with respect to %s' % ii_d)
    # Calculate the derivative of the density
    for ii_mo in range(len(mo_list)): 
      delta_rho[i] += (mo_spec[ii_mo]['occ_num'] * 
            2 * delta_mo_list[i,ii_mo]*mo_list[ii_mo])
      if len(ii_d) == 2:
        delta_rho[i] += mo_spec[ii_mo]['occ_num'] * 2 * delta2_mo_list[i][ii_mo]
  
  delta = (delta_rho,delta_rho.sum(axis=0)) if laplacian else (delta_rho,)
  
  return ((ao_list, mo_list, rho, delta_ao_list, delta_mo_list,) + delta 
        if return_components else (rho,) + delta)
  # rho_compute_no_slice 


#--- Support Code ---# 
# Information on atomic orbitals

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

def get_lxlylz(ao_spec,get_assign=False,bincount=False):
  '''Extracts the exponents lx, ly, lz for the Cartesian Gaussians.
  
  **Parameters:**
  
  ao_spec : 
    See :ref:`Central Variables` in the manual for details.
  get_assign : bool, optional
    Specifies, if the index of the atomic orbital shall be returned as well.
  
  **Returns:**
  
  lxlylz : numpy.ndarray, dtype=numpy.intc, shape = (NAO,3)
    Contains the expontents lx, ly, lz for the Cartesian Gaussians.
  assign : list of int, optional
    Contains the index of the atomic orbital in ao_spec.
  '''
  lxlylz = []
  assign = []
  for sel_ao in range(len(ao_spec)):
    if 'exp_list' in ao_spec[sel_ao].keys():
      l = ao_spec[sel_ao]['exp_list']
    else:
      l = exp[lquant[ao_spec[sel_ao]['type']]]
    lxlylz.extend(l)
    assign.extend([sel_ao]*len(l))
  if get_assign:
    if bincount:
      assign = numpy.bincount(assign)
    return (numpy.array(lxlylz,dtype=numpy.intc,order='C'), 
            numpy.array(assign,dtype=numpy.intc,order='C'))
  
  return numpy.array(lxlylz,dtype=numpy.intc,order='C') 

def validate_drv(drv):
  if drv is None or drv == '': return 0
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
  if not isinstance(mo,list):
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
    return numpy.sum(matrix)
  return integral

def slicer(N,vector=1e4,numproc=1):
  i = 0
  vector = 1 if int(vector) <= 0.0 else int(vector)
  sNum = int((N/(vector))+1)
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
