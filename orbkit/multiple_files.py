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
  order_using_analytical_overlap(fid_list,itype='molden')
'''

from copy import deepcopy
import numpy

from orbkit.read import main_read
from orbkit.display import display,init_display
from orbkit.analytical_integrals import get_ao_overlap, get_mo_overlap_matrix

geo_spec_all  = [] #: Contains all molecular geometries, i.e., :literal:`geo_spec`. (See `Central Variables`_ for details.)
geo_info      = [] #: See `Central Variables`_ for details.
ao_spec       = [] #: See `Central Variables`_ for details.
mo_coeff_all  = [] #: Contains all molecular orbital coefficients. List of numpy.ndarray
mo_energy_all = [] #: Contains all molecular orbital energies. List of numpy.ndarray
mo_occ_all    = [] #: Contains all molecular orbital occupations. List of numpy.ndarray
sym           = [] #: Python dictionary containing the molecular orbital symmetries and the corresponding position in mo_coeff_all, mo_energy_all, and mo_occ_all, respectively.
index_list    = [] #: After the execution of the ordering routine, it contains the new indices of the molecular orbitals. If index < 0, the molecular orbital changes its sign. shape=(Nfiles,NMO)

def read(fid_list,itype='molden',all_mo=True,nosym=False,**kwargs):
  '''Reads a list of input files.
  
  **Parameters:**
  
    fid_list : list of str
      List of input file names.
    itype : str, choices={'molden', 'gamess', 'gaussian.log', 'gaussian.fchk'}
        Specifies the type of the input files.
  
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
  n_r = len(fid_list)

  for i,filename in enumerate(fid_list):
    qc = main_read(filename, itype=itype, all_mo=all_mo,**kwargs)
    # Geo Section
    if i > 0 and (geo_old != qc.geo_info).sum():
      raise IOError('geo_info has changed!')
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
    
    for k,it in sym.iteritems():
      if k in sym_list:
        sym_list[k] = max(sym_list[k],it)
      else:
        sym_list[k] = it
      
        
    
  geo_spec_all = numpy.array(geo_spec_all)
  geo_info = qc.geo_info
  ao_spec = qc.ao_spec
  # Presorting of the MOs according to their symmetry

  sym = []
  for k,it in sym_list.iteritems():
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

def order_using_analytical_overlap(fid_list,itype='molden',**kwargs):
  '''Ordering routine using analytical overlap integrals between molecular
  orbitals. Set fid_list to None to omit the reading of input files.
  
  **Paramerters:**
  
  fid_list : list of str or None
    If not None, it contains the list of input file names.
  itype : str, choices={'molden', 'gamess', 'gaussian.log', 'gaussian.fchk'}
    Specifies the type of the input files.
  
  **Returns:**
  
  index_list : numpy.ndarray, shape=(Nfiles,NMO)
    Contains the new indices of the molecular orbitals. If index < 0, the 
    molecular orbital changes its sign.
  mo_overlap : numpy.ndarray, shape=((Nfiles - 1),NMO,NMO)
    Contains the overlap matrix between the molecular orbitals of two subsequent
    geometries, i.e., mo_overlap[i,j,k] corresponds to overlap between the 
    jth molecular orbital at geometry i to the kth molecular orbital at 
    geometry (i+1).
  
  **Global Variables:**
  
  geo_spec_all, geo_info, ao_spec, mo_coeff_all, mo_energy_all, mo_occ_all, sym, index_list
  '''
  global geo_spec_all, geo_info, ao_spec, mo_coeff_all, mo_energy_all, mo_occ_all, sym, index_list
  
  if fid_list is not None:
    read(fid_list,itype=itype,**kwargs)

  iterate= range(len(geo_spec_all)-1)

  mo_overlap = [[] for i in sym.iterkeys()]
  index_list = [[] for i in sym.iterkeys()]

  for s in sym.itervalues():
    shape = numpy.shape(mo_coeff_all[s])
    index_list[s] = numpy.ones((shape[0],shape[1]),dtype=int)
    index_list[s] *= numpy.arange(shape[1],dtype=int)
  
  display('Starting the ordering routine using the MO overlap...')
  c = 0
  for rr in iterate:
    c += 1
    if not c % (len(iterate)/10):
      display('\tFinished %d of %d geometries' % (c, len(iterate)))
    r1 = rr
    r2 = rr+1
    ao_overlap = get_ao_overlap(geo_spec_all[r1],geo_spec_all[r2],ao_spec)
    
    cs = 0
    for s in sym.itervalues():
      mo_coeff = mo_coeff_all[s]
      shape = numpy.shape(mo_coeff)
      overlap = get_mo_overlap_matrix(mo_coeff[r1],mo_coeff[r2],ao_overlap)
      
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
          mo_coeff[rr+1,[i,line_max],:] = mo_coeff[rr+1,[line_max,i],:]
          overlap[:,[i,line_max]] = overlap[:,[line_max,i]]
          index_list[s][rr+1,[i,line_max]] = index_list[s][rr+1,[line_max,i]]
      
      for i in range(shape[1]):
        # Change the signs
        mo_coeff[rr+1,i,:] *= numpy.sign(overlap[i,i])
        overlap[:,i] *= numpy.sign(overlap[i,i])
        index_list[s][rr+1,i] *= numpy.sign(overlap[i,i])
      
      mo_overlap[cs].append(overlap)
      cs += 1
      
      mo_coeff_all[s] = mo_coeff
      
      index = numpy.abs(index_list[s])[rr,:]
      mo_energy_all[s][rr,:] = mo_energy_all[s][rr,index]
      mo_occ_all[s][rr,:] = mo_occ_all[s][rr,index]
  
  tmp = []
  for i in mo_overlap:
    tmp.append(numpy.array(i))

  mo_overlap = tmp
  
  return index_list, mo_overlap

def order_using_extrapolation(fid_list,itype='molden',deg=1,
                              use_mo_values=False,matrix=None,**kwargs):
  '''Ordering routine using extrapolation of quantities related to the 
  molecular orbitals. Set fid_list to None to omit the reading of input files.
  
  The molecular orbital coefficients (If use_mo_values is False) are 
  extrapolated with a polynomial of degree :literal:`deg` and ordered by 
  minimizing a selected norm (default: Euclidian norm).
  
  **Paramerters:**
  
  fid_list : list of str or None
    If not None, it contains the list of input file names.
  itype : str, choices={'molden', 'gamess', 'gaussian.log', 'gaussian.fchk'}
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
  
  index_list = [None for i in sym.iterkeys()]
  
  for s,ii_s in sym.iteritems():
    display('Starting ordering of MOs of symmetry %s' % s)
    
    shape = numpy.shape(mo_coeff_all[ii_s])
    
    mu = 5e-2
    
    matrix = mo_coeff_all[ii_s]
    if use_mo_values:
      display('\tComputing molecular orbitals at the nuclear positions')
      matrix = compute_mo_list(geo_spec_all,ao_spec,matrix,iter_drv=[None, 'x', 'y', 'z'])
    
    
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

def order_manually(mo,i_0,i_1,r_range,index_list=None):
  '''Performs the ordering manually.
  '''
  shape = numpy.shape(mo)  
  
  if index_list == None:
    index_list = numpy.ones((shape[0],shape[1]),dtype=int)
    index_list *= numpy.arange(shape[1],dtype=int)
  
  def sign(x): return -1 if x < 0 else 1
  
  for rr in r_range:
    mo[rr,[abs(i_0),abs(i_1)],:] = sign(i_1)*mo[rr,[abs(i_1),abs(i_0)],:]
    index_list[rr,[abs(i_0),abs(i_1)]] = index_list[rr,[abs(i_1),abs(i_0)]]
  
  return mo, index_list


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
                             'index_list'],**kwargs):
  '''Writes all global variables specified in :literal:`variables` and all
  additional :literal:`kwargs` to an HDF5 file. 
  '''
  
  from orbkit.output import hdf5_open,hdf5_append
  
  # Save HDF5 File
  display('Saving HDF5 File to %s' % fid)
  for HDF5_file in hdf5_open(fid,mode='w'):
    for i in variables:
      if i == 'sym':
        s = globals()[i]
        sym = []
        for k,l in s.iteritems():
          sym.append([k,l])
        hdf5_append(numpy.array(sym),HDF5_file,name=i)
      elif i in globals():
        hdf5_append(globals()[i],HDF5_file,name=i)
      elif i not in kwargs:
        raise ValueError('Variable `%s` is not in globals() or in **kwargs' % i)
    
    for j in kwargs.iterkeys():
      hdf5_append(kwargs[j],HDF5_file,name=j)
    
  

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
  display('Reading HDF5 File to %s' % fid)
  for HDF5_file in hdf5_open(fid,mode='r'):
    for i in variables:
      globals()[i] = hdf52dict(i,HDF5_file)
      if i == 'sym':
        s = dict(globals()[i])
        globals()[i] = {}
        for k,l in s.iteritems():
          globals()[i][k] = int(l)
  

def construct_qc():
  '''Converts all global variables to a list of `QCinfo` classses.
  '''
  from orbkit.qcinfo import QCinfo
  QC = []
  for rr in range(len(geo_spec_all)):
    QC.append(QCinfo())
    QC[rr].geo_spec = geo_spec_all[rr]
    QC[rr].geo_info = geo_info
    QC[rr].ao_spec = ao_spec
    QC[rr].mo_spec = []
    for s,ii_s in sym.iteritems():
      for i,coeffs in enumerate(mo_coeff_all[ii_s][rr]):
        QC[rr].mo_spec.append({'coeffs': coeffs,
                               'energy' : mo_energy_all[ii_s][rr,i],
                               'occ_num' : mo_occ_all[ii_s][rr,i],
                               'sym': '%d.%s' % (i+1,s)})
  
  return QC

def compute_mo_list(geo_spec_all,ao_spec,mo_matrix,iter_drv=[None, 'x', 'y', 'z']):
  '''Computes the values of the molecular orbitals and, if requested, their 
  derivatives at the nuclear positions for a complete 
  mo_matrix (shape=(Nfiles,NMO,NAO)).'''
  from orbkit.core import ao_creator
  
  shape = numpy.shape(mo_matrix)
  mo_list = numpy.zeros((shape[0],shape[1],4*numpy.shape(geo_spec_all)[1]))

  for rr in range(shape[0]):
    geo_spec = geo_spec_all[rr]
    ao_spec = ao_spec
    x = geo_spec[:,0]
    y = geo_spec[:,1]
    z = geo_spec[:,2]
    N = len(x)
    for i,drv in enumerate(iter_drv):
      ao_list = ao_creator(geo_spec,ao_spec,exp_list=False,
                    is_vector=True,drv=drv,
                    x=x,y=y,z=z,N=len(x))
      for i_mo in range(shape[1]):
        for i_ao in range(shape[2]):
          mo_list[rr,i_mo,N*i+numpy.arange(N)] += mo_matrix[rr,i_mo,i_ao] * ao_list[i_ao,:]
    
    return mo_list

def data_interp(x,y,xnew,k=3):
  '''Interpolates a dataset y(x) to y(xnew) using B-Splines of order k.'''
  from scipy import interpolate
  tck = interpolate.splrep(x,y,s=0,k=k)
  ynew = interpolate.splev(xnew,tck,der=0)
  
  return ynew

def plot(mo_matrix,ydata2=None,symmetry='1',title='All',output_format='png',
         plt_dir='Plots',ylim=[-5, 5],thresh=0.1,x0=0):
  '''Plots all molecular orbital coefficients of one symmetry.'''
  import pylab as plt
  from matplotlib.ticker import MultipleLocator
  import os
  
  display('Plotting data of symmetry %s to %s' % (symmetry,plt_dir))
  if not os.path.exists(plt_dir):
    os.makedirs(plt_dir)
  
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
      X = numpy.arange(len(Y))+x0
      if max(numpy.abs(Y)) > thresh:
        curves.append(ax.plot(Y, colors[ij%len(colors)]+'o-' ,linewidth=1.5))
    
    x_label='${\\sf index}$';y_label='MO coefficients';
    plt.xlabel(x_label, fontsize=16);
    plt.ylabel(y_label, fontsize=16);
    plt.title('%s: %d.%s'%  (title,i+1,symmetry))
    plt.ylim(ylim)
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.grid(True, which='minor')
    ax.grid(True, which='both')
    return fig
  
  if output_format == 'pdf':
    from matplotlib.backends.backend_pdf import PdfPages
    output_fid = '%s.%s.pdf'% (title,symmetry.replace(' ','_'))
    display('\t%s' % output_fid)
    with PdfPages(os.path.join(plt_dir,output_fid)) as pdf:
      for i in range(shape[1]):
        fig = plot_mo(i)
        pdf.savefig(fig)
        plt.close()
  elif output_format == 'png':
    for i in range(shape[1]):
      fig = plot_mo(i)
      output_fid = '%d.%s.png' % (i+1,symmetry.replace(' ','_'))
      display('\t%s' % output_fid)
      fig.savefig(os.path.join(plt_dir, output_fid) ,format='png')
      plt.close()
  else:
    raise ValueError('output_format `%s` is not supported' % output_format)

def show_selected_mos(selected_mos,r0=0,steps=1,
                      N_=None,max_=None,min_=None):
  '''Uses orbkit to compute selected molecular orbitals and plots it with
  :func:`contour_mult_mo`.'''
  from orbkit import grid
  from orbkit.core import ao_creator,mo_creator
  r = range(r0,r0+steps)
  # Set grid parameters 
  grid.N_   = [  26,   1,  51] if N_   is None else N_
  grid.max_ = [ 3.0,   0,   6] if max_ is None else max_
  grid.min_ = [-3.0,   0,  -6] if min_ is None else min_
  
  if grid.N_[1] != 1:
    raise ValueError('show_selected_mos currently only supports an xz slice.' +
                     'N_[1] has to be 1. Please use contour_mult_mo directly,' + 
                     'if you want a different 2D-slice.')
  
  # Initialize grid
  grid.is_initialized = False
  grid.grid_init()
  for mo_sel in selected_mos:
    i,j = mo_sel.split('.')
    mo = []
    for rr in r:
      ao_list = ao_creator(geo_spec_all[rr],ao_spec)
      mo.append(mo_creator(ao_list,None,
                mo_coeff=mo_coeff_all[sym[j]][rr,int(i)-1,numpy.newaxis])[0][:,0,:])
    
    contour_mult_mo(grid.x,grid.z,mo,xlabel='x',ylabel='z',
                    title='MO:%s' % mo_sel,r0=r0)

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
    pic.set_title('Data Point %d' % (r0+i))
  
  f.subplots_adjust(left=0.15,bottom=0.05,top=0.95,right=0.95)
  f.show()
