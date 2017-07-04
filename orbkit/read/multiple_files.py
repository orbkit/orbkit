#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
'''Module for processing multiple output files of quantum chemical software.

**Capabilities**

  - Read a list of files
  - Order the molecular orbital coefficients by using analytical integrals,
    by extrapolation, or by hand (e.g. for the preparation of an interpolation)
  - Save and read the information obtained to and from an HDF5-file
  - Depiction of molecular orbitals

Example for the Execution::

  # Create a list of input filenames
  fid_list = []
  for i in range(101):
    fid_list.append('~/molden_files/%03d.molden' % i)

  init_display(name = 'MO_ordering')  # Specify a filename for the oklog file 
                                      # (We want a oklog file but we have no 
                                      # options.outputname here.)

  # Start reading and ordering routine.
  order_using_analytical_overlap(fid_list)
'''

from copy import deepcopy
import numpy
from time import time

from orbkit.display import display,init_display
from orbkit.analytical_integrals import get_ao_overlap, get_mo_overlap_matrix

from .tar import is_tar_file, get_all_files_from_tar
from . import main_read

geo_spec_all  = [] #: Contains all molecular geometries, i.e., :literal:`geo_spec`. (See :ref:`Central Variables` for details.)
geo_info      = [] #: See :ref:`Central Variables` for details.
ao_spec       = [] #: See :ref:`Central Variables` for details.
mo_coeff_all  = [] #: Contains all molecular orbital coefficients. List of numpy.ndarray
mo_energy_all = [] #: Contains all molecular orbital energies. List of numpy.ndarray
mo_occ_all    = [] #: Contains all molecular orbital occupations. List of numpy.ndarray
sym           = [] #: Python dictionary containing the molecular orbital symmetries and the corresponding position in mo_coeff_all, mo_energy_all, and mo_occ_all, respectively.
index_list    = [] #: After the execution of the ordering routine, it contains the new indices of the molecular orbitals. If index < 0, the molecular orbital changes its sign. shape=(Nfiles,NMO)

geo_spec_tck  = []
mo_coeff_tck  = []
mo_energy_tck = []
mo_occ_tck    = []

def read(fid_list,itype=None,all_mo=True,nosym=False, sort=True, **kwargs):
  '''Reads a list of input files.
  
  **Parameters:**
  
    fid_list : list of str
      List of input file names.
    itype : str, choices={None, 'tar', 'molden', 'gamess', 'gaussian.log', 'gaussian.fchk'}
      Specifies the type of the input files.
    sort: bool
      Sort input files by name.  
  
  **Global Variables:**
  
  geo_spec_all, geo_info, ao_spec, mo_coeff_all, mo_energy_all, mo_occ_all, sym
  '''
  global geo_spec_all, geo_info, ao_spec, mo_coeff_all, mo_energy_all, mo_occ_all, sym
  
  geo_spec_all = []
  MO_Spec = []
  mo_coeff_all = []
  mo_energy_all = []
  mo_occ_all = []
  # geo_info and ao_info have to stay unchanged
  geo_old = []
  ao_old = []
  
  sym_list = {}
  n_ao = {}

  #Check if fname poits to a tar archive and
  #read all files from archive if that is the case 
  if is_tar_file(fid_list):
    fid_list, itypes = get_all_files_from_tar(fid_list, sort=sort)

  n_r = len(fid_list)
  
  for i,fname in enumerate(fid_list):
    qc = main_read(fname, itype=itypes[i], all_mo=all_mo, **kwargs)
    # Geo Section
    if i > 0 and (geo_old != qc.geo_info).sum():
      raise IOError('qc.geo_info has changed!')
    else:
      geo_old = deepcopy(qc.geo_info)
    geo_spec_all.append(qc.geo_spec)
    # AO Section
    if (i > 0 and not
        numpy.alltrue([numpy.allclose(ao_old[j]['coeffs'],qc.ao_spec[j]['coeffs'])
                      for j in range(len(ao_old))]
                      )):
        raise IOError('qc.ao_spec has changed!')
    else:
        ao_old = deepcopy(qc.ao_spec)
    # MO Section
    sym = {}    
    MO_Spec.append(qc.mo_spec)
    
    for i,mo in enumerate(qc.mo_spec):
      if nosym:
        qc.mo_spec[i]['sym'] = '%d.1' % (i+1)
      key = mo['sym'].split('.')
      if key[1] not in sym.keys():
        sym[key[1]] = 0
        n_ao[key[1]] = len(qc.mo_spec[0]['coeffs'])
      sym[key[1]] += 1
    for k,it in sym.items():
      if k in sym_list:
        sym_list[k] = max(sym_list[k],it)
      else:
        sym_list[k] = it
  
  geo_spec_all = numpy.array(geo_spec_all)
  geo_info = qc.geo_info
  ao_spec = qc.ao_spec
  # Presorting of the MOs according to their symmetry
 
  sym = []
  for k,it in sym_list.items():
    sym.append((k,len(sym)))
    mo_coeff_all.append(numpy.zeros((n_r,it,n_ao[k])))
    mo_energy_all.append(numpy.zeros((n_r,it)))
    mo_occ_all.append(numpy.zeros((n_r,it)))
 
  sym = dict(sym)
  
  for i,spec in enumerate(MO_Spec):
    for j,mo in enumerate(spec):
      index,k = mo['sym'].split('.')
      
      index = int(index)-1
      
      mo_coeff_all[sym[k]][i,index,:] = mo['coeffs']
      mo_energy_all[sym[k]][i,index] = mo['energy']
      mo_occ_all[sym[k]][i,index] = mo['occ_num']

  return

def get_extrapolation(r1,r2,mo_coeff,deg=1,grid1d=None):
  '''Extrapolates the molecular orbital coefficients :literal:`mo_coeff` 
  using a polynomial of degree :literal:`deg`.
  
  **Paramerters:**
  
  r1 : int
    Specifies the index of the last known molecular orbital.
  r2 : int
    Specifies the index to which the molecular orbital coefficients are 
    extrapolated.
  deg : int
    Specifies the degree of the extrapolation polynomial.
  grid1d : list or numpy.1darray, optional
    Specifies the grid for the extrapolation.
  
  **Returns:**
  
  epol : numpy.ndarray, shape=(NMO,NAO))
    Contains the extrapolated molecular orbital coefficients.  
  '''
  if grid1d is None:
    grid1d = range(r2+1)
  
  if deg < 2:
    m = (mo_coeff[r1-1,:,:] - mo_coeff[r1,:,:])/float(grid1d[r1-1] - grid1d[r1])
    epol = (m * (grid1d[r2] - grid1d[r1]) + mo_coeff[r1,:,:])
  else:
    shape = mo_coeff.shape
    epol = numpy.zeros(shape[1:])
    for i in range(shape[1]):
      for j in range(shape[2]):
        x = grid1d[:r2]
        y = mo_coeff[:r2,i,j]
        z = numpy.polyfit(x, y, deg)
        epol[i,j] = numpy.poly1d(z)(grid1d[r2])
  return epol

def order_using_analytical_overlap(fid_list,itype=None,deg=0,numproc=1,
                                   **kwargs):
  '''Performs an ordering routine using analytical overlap integrals between 
  molecular orbitals. Set fid_list to None to omit the reading of input files.
  
  If :literal:`deg` is set to a value larger than zero, the molecular orbital 
  coefficients are extrapolated with a polynomial of degree :literal:`deg`,
  before computing the molecular orbital overlap matrix.
  
  **Paramerters:**
  
  fid_list : list of str or None
    If not None, it contains the list of input file names.
  itype : str, choices={None, 'tar', 'molden', 'gamess', 'gaussian.log', 'gaussian.fchk'}
    Specifies the type of the input files.
  deg : int, optional
    If greater than zero, specifies the degree of the extrapolation polynomial
    for the molecular orbital coefficients.  
  
  **Returns:**
  
  index_list : numpy.ndarray, shape=(Nfiles,NMO)
    Contains the new indices of the molecular orbitals. If index < 0, the 
    molecular orbital changes its sign.
  mo_overlap : numpy.ndarray, shape=((Nfiles - 1),NMO,NMO)
    Contains the overlap matrix between the molecular orbitals of two neighboring
    geometries, i.e., mo_overlap[i,j,k] corresponds to overlap between the 
    jth molecular orbital at geometry i to the kth molecular orbital at 
    geometry (i+1).
  
  **Global Variables:**
  
  geo_spec_all, geo_info, ao_spec, mo_coeff_all, mo_energy_all, mo_occ_all, sym, index_list
  '''
  global geo_spec_all, geo_info, ao_spec, mo_coeff_all, mo_energy_all, mo_occ_all, sym, index_list
  
  if fid_list is not None:
    read(fid_list,itype=itype,**kwargs)
    
  display('\nStarting the ordering routine using the molecular orbital overlap...')

  iterate= list(range(1,len(geo_spec_all)))
  
  if deg > 0:
    display('\tThe molecular orbital coefficients will be extrapolated')
    display('\tusing a least squares polynomial fit of degree %d.' % deg)
    std = numpy.array([numpy.std(i-geo_spec_all[0]) for i in geo_spec_all])
  
  mo_overlap = [[] for i in sym.keys()]
  index_list = [[] for i in sym.keys()]
  for s in sym.values():
    shape = numpy.shape(mo_coeff_all[s])
    index_list[s] = numpy.ones((shape[0],shape[1]),dtype=int)
    index_list[s] *= numpy.arange(shape[1],dtype=int)
  
  c = 0
  t = time()
  for rr in iterate:
    r1 = rr-1
    r2 = rr

    if (deg is None) or (deg > 0 and r1 >= deg):
      ao_overlap = get_ao_overlap(geo_spec_all[r2],geo_spec_all[r2],ao_spec)
    else:
      ao_overlap = get_ao_overlap(geo_spec_all[r1],geo_spec_all[r2],ao_spec)
    cs = 0
    for s in sym.values():
      mo_coeff = mo_coeff_all[s]
      shape = numpy.shape(mo_coeff)
      if deg is not None and deg > 0 and r1 >= deg:
        mo_r1 = get_extrapolation(r1,r2,mo_coeff,grid1d=std,deg=deg)
      else:
        mo_r1 = mo_coeff[r1]
      overlap = get_mo_overlap_matrix(mo_r1,mo_coeff[r2],ao_overlap,
                                      numproc=numproc)
      for i in range(shape[1]):
        # Iterate the rows of the overlap matrix
        line_max = None # variable for maximum value in the current row
        line_sort = numpy.argsort(numpy.abs(overlap[i,:]))[::-1] # sort the row
        for k in line_sort[::-1]:
          # Is this value the maximum in the current column?
          col_max = numpy.argmax(numpy.abs(overlap[:,k])) 
          if i == col_max:
            line_max = k
            break
        if line_max is not None:
          # Interchange the coefficients
          mo_coeff[r2,[i,line_max],:] = mo_coeff[r2,[line_max,i],:]
          overlap[:,[i,line_max]] = overlap[:,[line_max,i]]
          index_list[s][r2,[i,line_max]] = index_list[s][r2,[line_max,i]]
      
      for i in range(shape[1]):
        # Change the signs
        mo_coeff[r2,i,:] *= numpy.sign(overlap[i,i])
        overlap[:,i] *= numpy.sign(overlap[i,i])
        index_list[s][r2,i] *= numpy.sign(overlap[i,i])
      
      mo_overlap[cs].append(overlap)
      cs += 1
      
      mo_coeff_all[s] = mo_coeff
      
      index = numpy.abs(index_list[s])[r2,:]
      mo_energy_all[s][r2,:] = mo_energy_all[s][r2,index]
      mo_occ_all[s][r2,:] = mo_occ_all[s][r2,index]
    
    c += 1
    #if not c % int(numpy.ceil(len(iterate)/10.)):
    display('\tFinished %d of %d geometries (%.1f s)' % (c, len(iterate), time()-t))
    t = time()
  
  tmp = []
  for i in mo_overlap:
    tmp.append(numpy.array(i))

  mo_overlap = tmp
  
  return index_list, mo_overlap

def order_using_extrapolation(fid_list,itype=None,deg=1,
                              use_mo_values=False,matrix=None,**kwargs):
  '''Performs an ordering routine using extrapolation of quantities related to 
  the molecular orbitals. Set fid_list to None to omit the reading of input 
  files.
  
  The molecular orbital coefficients (If use_mo_values is False) are 
  extrapolated with a polynomial of degree :literal:`deg` and ordered by 
  minimizing a selected norm (default: Euclidian norm).
  
  **Paramerters:**
  
  fid_list : list of str or None
    If not None, it contains the list of input file names.
  itype : str, choices={None, 'tar', 'molden', 'gamess', 'gaussian.log', 'gaussian.fchk'}
    Specifies the type of the input files.
  deg : int
    Specifies the degree of the extrapolation polynomial.
  use_mo_values : bool, optional
    If True, some molecular orbital values and their derivatives are computed
    at the nuclear positions. The ordering routine is applied for those values
    instead.
  matrix : None or numpy.ndarray with shape=(Nfiles,N,M)
    If not None, contains the data to be ordered.
  
  **Returns:**
  
  :if matrix is None:
    - index_list
  :else:
    - matrix, index_list
    
  index_list : numpy.ndarray, shape=(Nfiles,NMO)
    Contains the new indices of the molecular orbitals. If index < 0, the 
    molecular orbital changes its sign.
  matrix : numpy.ndarray, shape=(Nfiles,N,M)
    Contains the ordered matrix.
  
  **Global Variables:**
  
  geo_spec_all, geo_info, ao_spec, mo_coeff_all, mo_energy_all, mo_occ_all, sym, index_list
  '''
  global geo_spec_all, geo_info, ao_spec, mo_coeff_all, mo_energy_all, mo_occ_all, sym, index_list
  
  # Read all input files
  if fid_list is not None:
    read(fid_list,itype=itype,**kwargs)

  radius = range(len(geo_spec_all)) #: We assume an equally spaced grid
  
  if deg < 2:
    function = order_mo
  else:
    function = order_mo_higher_deg
  
  if matrix is not None:
    display('\tOdering backward')
    matrix, index_list = function(matrix,index_list=index_list[ii_s],backward=True,mu=mu,deg=deg)
    display('\tOdering forward')
    matrix, index_list = function(matrix,index_list=index_list[ii_s],backward=False,mu=mu,deg=deg)
    return matrix, index_list
  
  index_list = [None for i in sym.keys()]
  
  for s,ii_s in sym.items():
    display('Starting ordering of MOs of symmetry %s' % s)
    
    shape = numpy.shape(mo_coeff_all[ii_s])
    
    mu = 5e-2
    
    matrix = mo_coeff_all[ii_s]
    if use_mo_values:
      display('\tComputing molecular orbitals at the nuclear positions')
      matrix = compute_mo_list(geo_spec_all,ao_spec,matrix,
                               iter_drv=[None, 'x', 'y', 'z'])
    
    
    display('\tOdering backward')
    matrix, index_list[ii_s] = function(matrix,index_list=index_list[ii_s],backward=True,mu=mu,deg=deg)
    display('\tOdering forward')
    matrix, index_list[ii_s] = function(matrix,index_list=index_list[ii_s],backward=False,mu=mu,deg=deg)
    
    for rr in range(shape[0]):
      index = numpy.abs(index_list[ii_s])[rr,:]
      sign = (-1)**(index_list[ii_s][rr] < 0)
      mo_energy_all[ii_s][rr,:] = mo_energy_all[ii_s][rr,index]
      mo_occ_all[ii_s][rr,:] = mo_occ_all[ii_s][rr,index]
      mo_coeff_all[ii_s][rr,:,:] = sign[:,numpy.newaxis]*mo_coeff_all[ii_s][rr,index,:] # numpy.array(matrix,copy=True)
    
  return index_list

def order_manually(matrix,i_0,i_1,r_range,using_sign=True):
  '''Performs the ordering manually.
  '''
  
  def sign(x): return -1 if x < 0 and using_sign else 1
  
  for rr in r_range:
    matrix[rr,[abs(i_0),abs(i_1)]] = sign(i_1)*matrix[rr,[abs(i_1),abs(i_0)]]
  
  return matrix

def order_mo(mo,index_list=None,backward=True,mu=1e-1,use_factor=False,**kwargs):
  '''Orders a 3d-matrix (shape=(Nfiles,NMO,NAO)) by interchanging the axis=1, 
  i.e., NMO, applying linear extrapolation.'''
  
  shape = numpy.shape(mo)
  
  if index_list == None:
    index_list = numpy.ones((shape[0],shape[1]),dtype=int)
    index_list *= numpy.arange(shape[1],dtype=int)
  
  if 'criterion' in kwargs:
    if kwargs['criterion'] == '1-norm':
      test = lambda x,y: numpy.sum(numpy.abs(x)) < numpy.sum(numpy.abs(y))
    if kwargs['criterion'] == '2-norm':
      test = lambda x,y: ((x**2).sum()) < ((y**2).sum())
    elif kwargs['criterion'] == 'infty-norm':
      test = lambda x,y: abs(x).max() < abs(y).max()
    elif kwargs['criterion'] == 'perc': 
      test = lambda x,y: ((x**2 < y**2).sum()/float(shape[2])) > 1./2.
    else:
      raise ValueError('creterion %s is not defined!' % kwargs['criterion'])
  else:
    # Take 2-norm by default
    test = lambda x,y: ((x**2).sum()) < ((y**2).sum())
  
  if backward:
    st = [-1, 1,-1]
  else:
    st = [ 0,-2, 1]
  
  x = 2.*st[2] # Extrapolate linearly to the next point    
  
  for i,i_0 in enumerate(range(shape[1])[:-1]):
    for rr in range(shape[0])[st[0]:st[1]:st[2]]:       
      f = numpy.ones(shape[2])
      if use_factor:
        is_larger = abs(mo[rr,i_0,:]) > mu
        f[is_larger] = 1/abs(mo[rr,i_0,is_larger])
      
      for ii_s,sign in enumerate([1,-1]):
        m = (sign*mo[rr+st[2],i_0,:] - mo[rr,i_0,:])/float(st[2])
        epol = (m * x + mo[rr,i_0,:])
        
        cp = ((f[:]*mo[rr+2*st[2],i_0,:] - f[:]*epol[:]))
        cm = ((-1*f[:]*mo[rr+2*st[2],i_0,:] - f[:]*epol[:]))
        is_smaller = test(cm,cp)
        current = cm if is_smaller else cp
        if ii_s == 0:
          diff = current
          i_1 = i_0
          new_signs = [sign,(-1)**is_smaller]
        elif test(current,diff):
          diff = current
          i_1 = i_0
          new_signs = [sign,(-1)**is_smaller]
        # Check other molecular orbitals
        for ik_index,ik in enumerate(range(shape[1])[i+1:]):
          cp = ((f[:]*mo[rr+2*st[2],ik,:] - f[:]*epol[:]))
          cm = ((-1*f[:]*mo[rr+2*st[2],ik,:] - f[:]*epol[:]))
          is_smaller = test(cm,cp)
          current = cm if is_smaller else cp
          
          if test(current,diff):
            diff = current
            i_1 = ik
            new_signs = [sign,(-1)**is_smaller]
      if i_0 != i_1:
        mo[rr+2*st[2],[i_0,i_1],:] = mo[rr+2*st[2],[i_1,i_0],:]
        index_list[rr+2*st[2],[i_0,i_1]] = index_list[rr+2*st[2],[i_1,i_0]]
      mo[rr+st[2],i_0,:] *= new_signs[0]
      mo[rr+2*st[2],i_0,:] *= new_signs[1]
      index_list[rr+st[2],i_0] *= new_signs[0]
      index_list[rr+2*st[2],i_0] *= new_signs[1]
  
  return mo, index_list

def order_mo_higher_deg(mo,index_list=None,backward=True,mu=1e-1,deg=2,**kwargs):
  '''Orders a 3d-matrix (shape=(Nfiles,NMO,NAO)) by interchanging the axis=1, 
  i.e., NMO, applying an extrapolation a polynomial fit with a Vandermonde matrix
  as implemented in numpy.'''
  
  shape = numpy.shape(mo)
  
  # Check if degree is correctly set
  if not isinstance(deg, int) or deg < 1 or deg > (shape[0]-1):
    raise IOError('Wrong choice for degree of the fitting polynomial!')
  
  display('\tusing a least squares polynomial fit of degree %d.' % deg)
  
  if index_list == None:
    index_list = numpy.ones((shape[0],shape[1]),dtype=int)
    index_list *= numpy.arange(shape[1],dtype=int)
  
  if 'criterion' in kwargs:
    if kwargs['criterion'] == '1-norm':
      test = lambda x,y: numpy.sum(numpy.abs(x)) < numpy.sum(numpy.abs(y))
    if kwargs['criterion'] == '2-norm':
      test = lambda x,y: ((x**2).sum()) < ((y**2).sum())
    elif kwargs['criterion'] == 'infty-norm':
      test = lambda x,y: abs(x).max() < abs(y).max()
    elif kwargs['criterion'] == 'perc': 
      test = lambda x,y: ((x**2 < y**2).sum()/float(shape[2])) > 1./2.
    else:
      raise ValueError('creterion %s is not defined!' % kwargs['criterion'])
  else:
    # Take 2-norm by default
    test = lambda x,y: ((x**2).sum()) < ((y**2).sum())
  
  if backward:
    st = [-(deg + 1), 0,-1]  
    x = numpy.arange(0,deg+1)
  else:
    st = [deg, -1,1]
    x = numpy.arange(-deg,1)
  
  for i,i_0 in enumerate(range(shape[1])[:-1]):
    for rr in range(shape[0])[st[0]:st[1]:st[2]]:
      epol = numpy.zeros(shape[2])
      for k in range(shape[2]):
        if mo[rr,i_0,k] != 0.:
          xnew = rr+st[2]
          y = mo[rr+x,i_0,k]
          z = numpy.polyfit(rr+x, y, deg)
          epol[k] = numpy.poly1d(z)(xnew)
      
      cp = ((mo[rr+st[2],i_0,:] - epol[:])**2)
      cm = ((-1*mo[rr+st[2],i_0,:] - epol[:])**2)
      is_smaller = test(cm,cp)
      current = cm if is_smaller else cp
      
      diff = current
      i_1 = i_0
      new_signs = [1,(-1)**is_smaller]
      # Check other molecular orbitals
      for ik_index,ik in enumerate(range(shape[1])[i+1:]):
        cp = ((mo[rr+st[2],ik,:] - epol[:])**2)
        cm = ((-1*mo[rr+st[2],ik,:] - epol[:])**2)
        is_smaller = test(cm,cp)
        current = cm if is_smaller else cp
        
        if test(current,diff):
          diff = current
          i_1 = ik
          new_signs = [1,(-1)**is_smaller]
      if i_0 != i_1:
        mo[rr+st[2],[i_0,i_1],:] = mo[rr+st[2],[i_1,i_0],:]
        index_list[rr+st[2],[i_0,i_1]] = index_list[rr+st[2],[i_1,i_0]]
      mo[rr,i_0,:] *= new_signs[0]
      mo[rr+st[2],i_0,:] *= new_signs[1]
      index_list[rr,i_0] *= new_signs[0]
      index_list[rr+st[2],i_0] *= new_signs[1]
  
  return mo, index_list

def order_pm(x,y,backward=True,mu=1e-1,use_factor=False):
  '''Outdated function to order exclusively the sign of a data set.
  '''
  if backward:
    st = [-2,1,-1]
  else:
    st = [1,-2,1]
  
  if numpy.ndim(y) == 1:
    diff = numpy.zeros(2)
    for rr in range(len(y))[st[0]:st[1]:st[2]]:
      m = (y[rr+st[2]]-y[rr])/(x[rr+st[2]]-x[rr])
      epol = m * (x[rr+2*st[2]]-x[rr]) + y[rr]
      for ii_d in range(2):
        diff[ii_d] = (((-1)**ii_d * y[rr+2*st[2]])-epol)**2
      if numpy.argmin(numpy.abs(diff)) == 1:
        y[rr+2*st[2]] = -y[rr+2*st[2]]
  elif numpy.ndim(y) == 2:
    y = numpy.array(y)
    shape = numpy.shape(y)
    for rr in range(shape[0])[st[0]:st[1]:st[2]]:
      for ii_s,sign in [(0,-1),(1,+1)]:
        f = numpy.ones(shape[1])
        if use_factor:
          f[numpy.abs(y[rr,:]) > mu] = 1/numpy.abs(y[rr,numpy.abs(y[rr,:]) > mu])
        m = (y[rr+st[2],:]-y[rr,:])/(x[rr+st[2]]-x[rr])
        epol = f[:]*(m * (x[rr+2*st[2]]-x[rr]) + y[rr,:])
        # Euclidean norm (2 norm)
        current = numpy.sum((f[:]*sign*y[rr+2*st[2],:] - epol[:])**2)
        # Current value
        if ii_s == 0:
          diff = current
          new_sign = sign
        elif current < diff:
          new_sign = sign
      y[rr+2*st[2],:] *= new_sign
        
  else:
    display('Function order_pm only works for vectors and 2D matrices')
  return y

def save_hdf5(fid,variables=['geo_info',
                             'geo_spec_all',
                             'ao_spec',
                             'mo_coeff_all', 
                             'mo_energy_all', 
                             'mo_occ_all', 
                             'sym',
                             'index_list'],hdf5mode='w',**kwargs):
  '''Writes all global variables specified in :literal:`variables` and all
  additional :literal:`kwargs` to an HDF5 file. 
  '''
  
  from orbkit.output import hdf5_open,hdf5_append
  
  # Save HDF5 File
  display('Saving Hierarchical Data Format file (HDF5) to %s...' % fid)
  data_stored = []
  for HDF5_file in hdf5_open(fid,mode=hdf5mode):
    for i in variables:
      if i in globals():
        data = globals()[i]
        if not (data == [] or data is None):          
          if i == 'sym':
            data = numpy.array([[k,l] for k,l in data.items()])
          hdf5_append(data,HDF5_file,name=i)
          data_stored.append(i)
      elif i not in kwargs:
        raise ValueError('Variable `%s` is not in globals() or in **kwargs' % i)
    
    for j in kwargs.keys():
      hdf5_append(kwargs[j],HDF5_file,name=j)
      data_stored.append(j)
    
  display('\tContent: ' + ', '.join(data_stored))

def read_hdf5(fid,variables=['geo_info',
                             'geo_spec_all',
                             'ao_spec',
                             'mo_coeff_all', 
                             'mo_energy_all', 
                             'mo_occ_all', 
                             'sym',
                             'index_list']):
  '''Reads all variables specified in :literal:`variables` from an HDF5 file
  created with :literal:`write_hdf5` and appends this data to the globals() of 
  this module. 
  '''
  
  from orbkit.output import hdf5_open,hdf52dict
  
  # Read HDF5 File
  display('Reading Hierarchical Data Format file (HDF5) File from %s' % fid)
  data_stored = []
  for HDF5_file in hdf5_open(fid,mode='r'):
    for i in variables:
      try:
        globals()[i] = hdf52dict(i,HDF5_file)
        data_stored.append(i)
        if i == 'sym':
          s = dict(globals()[i])
          globals()[i] = {}
          for k,l in s.items():
            globals()[i][k] = int(l)
      except KeyError:
        pass
  
  if not data_stored:
    raise IOError('Could not find any data in `%s` for the selected `variables`.' 
                  % fid)
  
  display('\tFound: ' + ', '.join(data_stored))

def construct_qc():
  '''Converts all global variables to a list of `QCinfo` classes.
  '''
  from orbkit.qcinfo import QCinfo
  QC = []
  for rr in range(len(geo_spec_all)):
    QC.append(QCinfo())
    QC[rr].geo_spec = geo_spec_all[rr]
    QC[rr].geo_info = geo_info
    QC[rr].ao_spec = ao_spec
    QC[rr].mo_spec = []
    for s,ii_s in sym.items():
      for i,coeffs in enumerate(mo_coeff_all[ii_s][rr]):
        QC[rr].mo_spec.append({'coeffs': coeffs,
                               'energy' : mo_energy_all[ii_s][rr,i],
                               'occ_num' : mo_occ_all[ii_s][rr,i],
                               'sym': '%d.%s' % (i+1,s)})
  
  return QC

def compute_mo_list(geo_spec_all,ao_spec,mo_matrix,
                    iter_drv=[None, 'x', 'y', 'z']):
  '''Computes the values of the molecular orbitals and, if requested, their 
  derivatives at the nuclear positions for a complete 
  mo_matrix (shape=(Nfiles,NMO,NAO)).'''
  from orbkit.core import ao_creator
  
  shape = numpy.shape(mo_matrix)
  mo_list = numpy.zeros((shape[0],shape[1],4*numpy.shape(geo_spec_all)[1]))

  for rr in range(shape[0]):
    geo_spec = geo_spec_all[rr]
    x = geo_spec[:,0]
    y = geo_spec[:,1]
    z = geo_spec[:,2]
    N = len(x)
    for i,drv in enumerate(iter_drv):
      ao_list = ao_creator(geo_spec,ao_spec,
                           exp_list=False,
                           is_vector=True,
                           drv=drv,
                           x=x,y=y,z=z)
      for i_mo in range(shape[1]):
        for i_ao in range(shape[2]):
          mo_list[rr,i_mo,N*i+numpy.arange(N)] += mo_matrix[rr,i_mo,i_ao] * ao_list[i_ao,:]
    
    return mo_list

def data_interp(x,y,xnew,k=3,der=0,s=0,**kwargs):
  '''Interpolates a dataset y(x) to y(xnew) using B-Splines of order k.'''
  from scipy import interpolate
  tck = interpolate.splrep(x,y,s=s,k=k)
  ynew = interpolate.splev(xnew,tck,der=der)
  
  return ynew

def splrep_all(x,k=3,**kwargs):
  from scipy import interpolate
  global geo_spec_tck, mo_coeff_tck, mo_energy_tck, mo_occ_tck
  
  geo_spec_tck  = []
  mo_coeff_tck  = []
  mo_energy_tck = []
  mo_occ_tck    = []
  
  shape = geo_spec_all.shape
  for i in range(shape[1]):
    geo_spec_tck.append([])
    for j in range(shape[2]):
      geo_spec_tck[-1].append(interpolate.splrep(x,geo_spec_all[:,i,j],
                              k=k,**kwargs))
  
  for i_mo in range(len(mo_coeff_all)):
    mo_coeff_tck.append([])
    mo_energy_tck.append([])
    mo_occ_tck.append([])
    shape = mo_coeff_all[i_mo].shape    
    for i in range(shape[1]):
      mo_coeff_tck[-1].append([])
      mo_energy_tck[-1].append(interpolate.splrep(x,mo_energy_all[i_mo][:,i],
                               k=k,**kwargs))
      mo_occ_tck[-1].append(interpolate.splrep(x,mo_occ_all[i_mo][:,i],
                            k=k,**kwargs))
      for j in range(shape[2]):
        mo_coeff_tck[-1][-1].append(interpolate.splrep(x,
                                    mo_coeff_all[i_mo][:,i,j],
                                    k=k,**kwargs))

def interpolate_all(x,xnew,k=3,**kwargs):
  '''Interpolates a dataset y(x) to y(xnew) using B-Splines of order k.'''
  from scipy import interpolate
  global geo_spec_all, mo_coeff_all, mo_energy_all, mo_occ_all 
  
  shape = list(geo_spec_all.shape)
  shape[0] = len(xnew)
  tmp = numpy.zeros(shape)  
  for i in range(shape[1]):
    for j in range(shape[2]):
      tmp[:,i,j] = data_interp(x,geo_spec_all[:,i,j],xnew,k=k,**kwargs)
  geo_spec_all = numpy.copy(tmp)
  
  for i_mo in range(len(mo_coeff_all)):
    shape = list(mo_coeff_all[i_mo].shape)
    shape[0] = len(xnew)
    tmp = numpy.zeros(shape)  
    for i in range(shape[1]):
      for j in range(shape[2]):
        tmp[:,i,j] = data_interp(x,mo_coeff_all[i_mo][:,i,j],xnew,k=k,**kwargs)
    mo_coeff_all[i_mo] = numpy.copy(tmp)
    
    shape = list(mo_energy_all[i_mo].shape)
    shape[0] = len(xnew)
    tmp = numpy.zeros(shape)  
    for i in range(shape[1]):
      tmp[:,i] = data_interp(x,mo_energy_all[i_mo][:,i],xnew,k=k,**kwargs)
    mo_energy_all[i_mo] = numpy.copy(tmp)
    
    shape = list(mo_occ_all[i_mo].shape)
    shape[0] = len(xnew)
    tmp = numpy.zeros(shape)  
    for i in range(shape[1]):
      tmp[:,i] = data_interp(x,mo_occ_all[i_mo][:,i],xnew,k=k,**kwargs)
    mo_occ_all[i_mo] = numpy.copy(tmp)

def plot(mo_matrix,symmetry='1',title='All',x_label='index',
         y_label='MO coefficients',output_format='png',
         plt_dir='Plots',ylim=None,thresh=0.1,x0=0,grid=True,x_grid=None,**kwargs):
  '''Plots all molecular orbital coefficients of one symmetry.'''
  import pylab as plt
  from matplotlib.ticker import MultipleLocator
  import os
  
  display('Plotting data of symmetry %s to %s/' % (symmetry,plt_dir))
  if not os.path.exists(plt_dir):
    os.makedirs(plt_dir)
  
  if numpy.ndim(mo_matrix) == 2:
    mo_matrix = mo_matrix[:,numpy.newaxis,:]
  
  shape = numpy.shape(mo_matrix)
  
  colors = 'bgrcmyk'
  
  def plot_mo(i):
    fig=plt.figure()
    plt.rc('xtick', labelsize=16) 
    plt.rc('ytick', labelsize=16)
    ax = plt.subplot(111)
    curves=[]
    for ij in range(shape[2]):
      Y = mo_matrix[:,i,ij]
      if x_grid is None:
        X = numpy.arange(len(Y))+x0
      else:
        X = x_grid
      if max(numpy.abs(Y)) > thresh:
        curves.append(ax.plot(X,Y, colors[ij%len(colors)]+'.-' ,linewidth=1.5))
    
    
    plt.xlabel(x_label, fontsize=16);
    plt.ylabel(y_label, fontsize=16);
    plt.title('%s: %d.%s'%  (title,i+1,symmetry))
    plt.ylim(ylim)
    #ax.xaxis.set_minor_locator(MultipleLocator(1))
    #ax.xaxis.grid(grid, which='minor')
    #ax.grid(grid, which='both')
    return fig
  
  if output_format == 'pdf':
    from matplotlib.backends.backend_pdf import PdfPages
    output_fid = '%s.%s.pdf'% (title,symmetry.replace(' ','_'))
    display('\t%s' % output_fid)
    with PdfPages(os.path.join(plt_dir,output_fid)) as pdf:
      for i in range(shape[1]):
        fig = plot_mo(i)
        pdf.savefig(fig,**kwargs)
        plt.close()
  elif output_format == 'png':
    for i in range(shape[1]):
      fig = plot_mo(i)
      output_fid = '%d.%s.png' % (i+1,symmetry.replace(' ','_'))
      display('\t%s' % output_fid)
      fig.savefig(os.path.join(plt_dir, output_fid),format='png',**kwargs)
      plt.close()
  else:
    raise ValueError('output_format `%s` is not supported' % output_format)

def show_selected_mos(selected_mos,r0=0,steps=1,select_slice='xz',where=0.0,
                      npts=[26,51],minpts=[-3,-6],maxpts=[3,6],nuclear_pos='x'):
  '''Uses orbkit to compute selected molecular orbitals and plots it with
  :func:`contour_mult_mo`.'''
  from orbkit import grid
  from orbkit.core import ao_creator,mo_creator
  r = range(r0,r0+steps)
  
  grid.N_ = [1,1,1]
  grid.min_ = [0,0,0]
  grid.max_ = [0,0,0]
  if select_slice == 'xy':
    k = [0,1]
    grid.min_[2] += where
    grid.max_[2] += where
  elif select_slice == 'yz':
    k = [1,2]
    grid.min_[0] += where
    grid.max_[0] += where
  elif select_slice == 'xz':
    k = [0,2]
    grid.min_[1] += where
    grid.max_[1] += where
  else:
    raise ValueError('`show_selected_mos` currently only' + 
                     'supports slices parallel to the following planes:' +
                     'select_slice = `xy`, `yz`, or `xz`')
  for i,j in enumerate(k):
    grid.N_[j] = npts[i]
    grid.min_[j] = minpts[i]
    grid.max_[j] = maxpts[i]
  
  # Initialize grid
  grid.is_initialized = False
  grid.grid_init(force=True)
  xyz = grid.x,grid.y,grid.z
  for mo_sel in selected_mos:
    i,j = mo_sel.split('.')
    mo = []
    for rr in r:
      ao_list = ao_creator(geo_spec_all[rr],ao_spec)
      mo.append(mo_creator(ao_list,mo_coeff_all[sym[j]][rr,int(i)-1,numpy.newaxis])[0].reshape(tuple(npts)))
    
    f, pics = contour_mult_mo(xyz[k[0]],xyz[k[1]],mo,
                    xlabel=select_slice[0],ylabel=select_slice[1],
                    title='MO:%s' % mo_sel,r0=r0)
    for i,pic in enumerate(pics):
      pic.plot(geo_spec_all[rr,:,k[1]],geo_spec_all[rr,:,k[0]],nuclear_pos,
               markersize=10,markeredgewidth=2)

def contour_mult_mo(x,y,mo,xlabel='x',ylabel='y',title='',r0=0):
  '''Uses matplotlib to show slices of a molecular orbitals.'''
  import matplotlib.pyplot as plt
  
  # Plot slices
  f, pics = \
              plt.subplots(len(mo),1,sharex=True,sharey=True,figsize=(6,2+4*len(mo)))
  plt.suptitle(title)
  vmax = numpy.max(numpy.abs(mo))
  for i,pic in enumerate(pics):
    pic.contour(y,x,mo[i],50,linewidths=0.5,colors='k')
    pic.contourf(\
        y,x,mo[i],50,cmap=plt.cm.rainbow,vmax=vmax,vmin=-vmax)
    pic.set_ylabel(xlabel)  
    pic.set_xlabel(ylabel)  
    pic.set_title('Data Point %d' % (r0+i))
  
  f.subplots_adjust(left=0.15,bottom=0.05,top=0.95,right=0.95)
  f.show()
  return f,pics
