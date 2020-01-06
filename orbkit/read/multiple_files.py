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
from collections import OrderedDict as odict

from orbkit.qcinfo import QCinfo
from orbkit.orbitals import AOClass, MOClass
from orbkit.display import display,init_display
from orbkit.analytical_integrals import get_ao_overlap, get_mo_overlap_matrix

from .tar import is_tar_file, get_all_files_from_tar
from . import main_read

class Multi():
  def __init__(self, data=None):
    self.geo_spec_all  = [] #: Contains all molecular geometries, i.e., :literal:`geo_spec`. (See :ref:`Central Variables` for details.)
    self.geo_info      = [] #: See :ref:`Central Variables` for details.
    self.ao_spec       = [] #: See :ref:`Central Variables` for details.
    self.mo_coeff_all  = [] #: Contains all molecular orbital coefficients. List of numpy.ndarray
    self.mo_energy_all = [] #: Contains all molecular orbital energies. List of numpy.ndarray
    self.mo_occ_all    = [] #: Contains all molecular orbital occupations. List of numpy.ndarray
    self.sym           = [] #: Python dictionary containing the molecular orbital self.symmetries and the corresponding position in self.mo_coeff_all, self.mo_energy_all, and self.mo_occ_all, respectively.
    self.index_list    = [] #: After the execution of the ordering routine, it contains the new indices of the molecular orbitals. If index < 0, the molecular orbital changes its sign. shape=(Nfiles,NMO)

    self.geo_spec_tck  = []
    self.mo_coeff_tck  = []
    self.mo_energy_tck = []
    self.mo_occ_tck    = []
    self.MO_Spec = []
    self.QC = []

    if data:
      if isinstance(data['ao_spec'], numpy.ndarray):
        ao_spec = data['ao_spec'][numpy.newaxis][0]
      else:
        ao_spec = data['ao_spec']
      self.ao_spec = AOClass(restart=ao_spec)
      self.geo_info = data['geo_info']
      self.geo_spec_all = data['geo_spec_all']
      if isinstance(data['mo_data'], numpy.ndarray):
        mo_data = data['mo_data'][numpy.newaxis][0]
      else:
        mo_data = data['mo_data']
      self.mo_list_parsing(mo_data)
      if isinstance(data['sym'], numpy.ndarray):
        self.sym = data['sym'][numpy.newaxis][0]
      else:
        self.sym = data['sym']
      self.index_list = data['index_list']

  def read(self,fid_list,itype='auto',all_mo=True,nosym=False, sort=True, **kwargs_all):
    '''Reads a list of input files.
    
    **Parameters:**
    
      fid_list : list of str
        List of input file names.
      itype : str, choices={'auto', 'tar', 'molden', 'gamess', 'gaussian.log', 'gaussian.fchk'}
        Specifies the type of the input files.
      sort: bool
        Sort input files by name.
    '''
    # self.geo_info and ao_info have to stay unchanged
    geo_old = []
    ao_old = []
    
    sym_list = {}
    n_ao = {}

    #Check if fname poits to a tar archive and
    #read all files from archive if that is the case 
    if is_tar_file(fid_list):
      fid_list, itypes = get_all_files_from_tar(fid_list, sort=sort)
    else:
      itypes = [[itype]*len(fid_list)][0]

    for i,fname in enumerate(fid_list):

      kwargs = kwargs_all['kwargs'][i] if 'kwargs' in kwargs_all.keys() else kwargs_all
      qc = main_read(fname, itype=itypes[i], all_mo=all_mo, **kwargs)
      # Geo Section
      if i > 0 and (geo_old != qc.geo_info).sum():
        raise IOError('qc.geo_info has changed!')
      else:
        geo_old = deepcopy(qc.geo_info)
      self.geo_spec_all.append(qc.geo_spec)
      # AO Section
      if (i > 0 and not
          numpy.alltrue([numpy.allclose(ao_old[j]['coeffs'],qc.ao_spec[j]['coeffs'])
                        for j in range(len(ao_old))]
                        )):
          raise IOError('qc.ao_spec has changed!')
      else:
          ao_old = deepcopy(qc.ao_spec)

      # MO Section
      sym_tmp = {}    
      self.MO_Spec.append(qc.mo_spec)
      for i,mo in enumerate(qc.mo_spec):
        if nosym:
          qc.mo_spec[i]['sym'] = '%d.1' % (i+1)
        key = mo['sym'].split('.')
        if key[1] not in sym_tmp.keys():
          sym_tmp[key[1]] = 0
          n_ao[key[1]] = len(qc.mo_spec[0]['coeffs'])
        sym_tmp[key[1]] += 1

      for k,it in sym_tmp.items():
        if k in sym_list:
          sym_list[k] = max(sym_list[k],it)
        else:
          sym_list[k] = it
    
    self.geo_spec_all = numpy.array(self.geo_spec_all)
    self.geo_info = qc.geo_info
    self.ao_spec = qc.ao_spec

    # Presorting of the MOs according to their self.symmetry 
    n_r = len(fid_list)
    self.sym = []
    for k in sorted(sym_list.keys()):
      it = sym_list[k]
      self.sym.append((k,len(self.sym)))
      self.mo_coeff_all.append(numpy.zeros((n_r,it,n_ao[k])))
      self.mo_energy_all.append(numpy.zeros((n_r,it)))
      self.mo_occ_all.append(numpy.zeros((n_r,it)))

    self.sym = odict(self.sym)

    for i,spec in enumerate(self.MO_Spec):
      for j,mo in enumerate(spec):
        index,k = mo['sym'].split('.')
        
        index = int(index)-1
        
        self.mo_coeff_all[self.sym[k]][i,index,:] = mo['coeffs']
        self.mo_energy_all[self.sym[k]][i,index] = mo['energy']
        self.mo_occ_all[self.sym[k]][i,index] = mo['occ_num']
    return

  def get_extrapolation(self,r1,r2,mo_coeff,deg=1,grid1d=None):
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

  def order_using_analytical_overlap(self,fid_list=None,itype=None,deg=0,numproc=1,
                                     **kwargs):
    '''Performs an ordering routine using analytical overlap integrals between 
    molecular orbitals. Set fid_list to None to omit the reading of input files.
    
    If :literal:`deg` is set to a value larger than zero, the molecular orbital 
    coefficients are extrapolated with a polynomial of degree :literal:`deg`,
    before computing the molecular orbital overlap matrix.
    
    **Paramerters:**
    
    fid_list : list of str or None
      If not None, it contains the list of input file names.
    itype : str, choices={'auto', 'tar', 'molden', 'gamess', 'gaussian.log', 'gaussian.fchk'}
      Specifies the type of the input files.
    deg : None|int, optional
      - If deg is None, atomic orbitals of two successive geometries will be assumed 
        to be on the same positions.
      - If greater than zero, specifies the degree of the extrapolation polynomial
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
    '''
    
    if fid_list is not None:
      self.read(fid_list,itype=itype,**kwargs)
      
    display('\nStarting the ordering routine using the molecular orbital overlap...')

    iterate= list(range(1,len(self.geo_spec_all)))
    
    if deg is not None and deg > 0:
      display('\tThe molecular orbital coefficients will be extrapolated')
      display('\tusing a least squares polynomial fit of degree %d.' % deg)
      std = numpy.array([numpy.std(i-self.geo_spec_all[0]) for i in self.geo_spec_all])
   
    sym_sorted_keys = sorted(self.sym.keys())
    mo_overlap = [[] for i in sym_sorted_keys]
    index_list = [[] for i in sym_sorted_keys]
    for ik in sym_sorted_keys:
      s = self.sym[ik]
      shape = numpy.shape(self.mo_coeff_all[s])
      index_list[s] = numpy.ones((shape[0],shape[1]),dtype=int)
      index_list[s] *= numpy.arange(shape[1],dtype=int)
    
    c = 0
    t = time()
    for rr in iterate:
      r1 = rr-1
      r2 = rr

      if (deg is None) or (deg > 0 and r1 >= deg):
        ao_overlap = get_ao_overlap(self.geo_spec_all[r2],self.geo_spec_all[r2],self.ao_spec)
      else:
        ao_overlap = get_ao_overlap(self.geo_spec_all[r1],self.geo_spec_all[r2],self.ao_spec)

      cs = 0
      for ik in sym_sorted_keys:
        s = self.sym[ik]
        mo_coeff = self.mo_coeff_all[s]
        shape = numpy.shape(mo_coeff)
        if deg is not None and deg > 0 and r1 >= deg:
          mo_r1 = self.get_extrapolation(r1,r2,mo_coeff,grid1d=std,deg=deg)
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
        
        self.mo_coeff_all[s] = mo_coeff
        
        index = numpy.abs(index_list[s])[r2,:]
        self.mo_energy_all[s][r2,:] = self.mo_energy_all[s][r2,index]
        self.mo_occ_all[s][r2,:] = self.mo_occ_all[s][r2,index]
      
      c += 1
      #if not c % int(numpy.ceil(len(iterate)/10.)):
      display('\tFinished %d of %d geometries (%.1f s)' % (c, len(iterate), time()-t))
      t = time()
    
    tmp = []
    for i in mo_overlap:
      tmp.append(numpy.array(i))
    
    mo_overlap = tmp
    
    return index_list, mo_overlap

  def order_using_extrapolation(self,fid_list=None,itype=None,deg=1,
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
    '''
    
    mu = 5e-2

    # Read all input files
    if fid_list is not None:
      self.read(fid_list,itype=itype,**kwargs)

    radius = range(len(self.geo_spec_all)) #: We assume an equally spaced grid
    
    if deg < 2:
      function = self.order_mo
    else:
      function = self.order_mo_higher_deg
    
    if matrix is not None:
      display('\tOdering backward')
      matrix, index_list = function(matrix,index_list=self.index_list,backward=True,mu=mu,deg=deg)
      display('\tOdering forward')
      matrix, index_list = function(matrix,index_list=self.index_list,backward=False,mu=mu,deg=deg)
      return matrix, index_list
    
    index_list = [None for i in self.sym.keys()]
    
    for s,ii_s in self.sym.items():
      display('Starting ordering of MOs of self.symmetry %s' % s)
      
      shape = numpy.shape(self.mo_coeff_all[ii_s])
      
      matrix = self.mo_coeff_all[ii_s]
      if use_mo_values:
        display('\tComputing molecular orbitals at the nuclear positions')
        matrix = self.compute_mo_list(self.geo_spec_all,self.ao_spec,matrix,
                                 iter_drv=[None, 'x', 'y', 'z'])
      
      
      display('\tOdering backward')
      matrix, index_list[ii_s] = function(matrix,index_list=index_list[ii_s],backward=True,mu=mu,deg=deg)
      display('\tOdering forward')
      matrix, index_list[ii_s] = function(matrix,index_list=index_list[ii_s],backward=False,mu=mu,deg=deg)
      
      for rr in range(shape[0]):
        index = numpy.abs(index_list[ii_s])[rr,:]
        sign = (-1)**(index_list[ii_s][rr] < 0)
        self.mo_energy_all[ii_s][rr,:] = self.mo_energy_all[ii_s][rr,index]
        self.mo_occ_all[ii_s][rr,:] = self.mo_occ_all[ii_s][rr,index]
        self.mo_coeff_all[ii_s][rr,:,:] = sign[:,numpy.newaxis]*self.mo_coeff_all[ii_s][rr,index,:] # numpy.array(matrix,copy=True)
      
    return index_list

  def order_manually(self,matrix,i_0,i_1,r_range,using_sign=True):
    '''Performs the ordering manually.
    '''
    
    def sign(x): return -1 if x < 0 and using_sign else 1
    
    for rr in r_range:
      matrix[rr,[abs(i_0),abs(i_1)]] = sign(i_1)*matrix[rr,[abs(i_1),abs(i_0)]]
    
    return matrix

  def order_mo(self,mo,index_list=None,backward=True,mu=1e-1,use_factor=False,**kwargs):
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

  def order_mo_higher_deg(self,mo,index_list=None,backward=True,mu=1e-1,deg=2,**kwargs):
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

  def order_pm(self,x,y,backward=True,mu=1e-1,use_factor=False):
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

  def mo_list_parsing(self,indata=None):
    # Parses lists of mos to and from the native Orbkit format
    parameters = {'energy': self.mo_energy_all, 
                  'coeff': self.mo_coeff_all, 
                  'occ': self.mo_occ_all}
    if indata is None:
      outdata = {}
      for param in parameters:
        for i in range(len(parameters[param])):
          outdata[param+'#'+str(i)] = parameters[param][i]
          print(param,parameters[param][i].shape)
      return outdata
    else:
      order = {'energy': [], 
               'coeff': [], 
               'occ': []}
      param_tmp = {'energy': [], 
               'coeff': [], 
               'occ': []}
      for name in indata:
        param = name.split('#')[0]
        param_tmp[param].append(indata[name])
        order[param].append(int(name.split('#')[-1]))
      for param in parameters:
        sort = numpy.argsort(numpy.array(order[param],dtype=numpy.intc))
        for s in sort:
          parameters[param].append(param_tmp[param][s])
      return

  def todict(self):
    data = {}
    data['ao_spec'] = self.ao_spec.todict()
    data['geo_info'] = self.geo_info
    data['geo_spec_all'] = self.geo_spec_all
    data['mo_data'] = self.mo_list_parsing()
    data['sym'] = self.sym
    data['index_list'] = self.index_list
    data['parent_class_name'] = self.__module__ + '.' + self.__class__.__name__
    return data

  def construct_qc(self, all_mo=True):
    '''Converts all global variables to a list of `QCinfo` classes.
    '''
    self.QC = []
    ilumo = None
    for rr in range(len(self.geo_spec_all)):
      qc = QCinfo()
      qc.geo_spec = self.geo_spec_all[rr]
      qc.geo_info = self.geo_info
      qc.ao_spec = self.ao_spec
      qc.mo_spec = []
      for s,ii_s in self.sym.items():
        for i,coeffs in enumerate(self.mo_coeff_all[ii_s][rr]):
          qc.mo_spec.append({'coeffs': coeffs,
                             'energy' : self.mo_energy_all[ii_s][rr,i],
                             'occ_num' : self.mo_occ_all[ii_s][rr,i],
                             'sym': '%d.%s' % (i+1,s)})
      qc.ao_spec.update()
      qc.mo_spec = MOClass(qc.mo_spec)
      qc.mo_spec.update()
      if not all_mo:
        ilumo = max(ilumo or 0, qc.mo_spec.get_lumo())
      
      self.QC.append(qc)
    
    if not all_mo:
      for i in range(len(self.QC)):
        self.QC[i].mo_spec = self.QC[i].mo_spec[slice(None,ilumo)]
    
    return self.QC

  def compute_mo_list(self,ao_spec,mo_matrix,
                      iter_drv=[None, 'x', 'y', 'z']):
    '''Computes the values of the molecular orbitals and, if requested, their 
    derivatives at the nuclear positions for a complete 
    mo_matrix (shape=(Nfiles,NMO,NAO)).'''
    from orbkit.core import ao_creator
    
    shape = numpy.shape(mo_matrix)
    mo_list = numpy.zeros((shape[0],shape[1],4*numpy.shape(self.geo_spec_all)[1]))

    for rr in range(shape[0]):
      geo_spec = self.geo_spec_all[rr]
      x = geo_spec[:,0]
      y = geo_spec[:,1]
      z = geo_spec[:,2]
      N = len(x)
      for i,drv in enumerate(iter_drv):
        ao_list = ao_creator(geo_spec,self.ao_spec,
                             exp_list=False,
                             is_vector=True,
                             drv=drv,
                             x=x,y=y,z=z)
        for i_mo in range(shape[1]):
          for i_ao in range(shape[2]):
            mo_list[rr,i_mo,N*i+numpy.arange(N)] += mo_matrix[rr,i_mo,i_ao] * ao_list[i_ao,:]
      
      return mo_list

  def data_interp(self,x,y,xnew,k=3,der=0,s=0,**kwargs):
    '''Interpolates a dataset y(x) to y(xnew) using B-Splines of order k.'''
    from scipy import interpolate
    tck = interpolate.splrep(x,y,s=s,k=k)
    ynew = interpolate.splev(xnew,tck,der=der)
    
    return ynew

  def splrep_all(self,x,k=3,**kwargs):
    from scipy import interpolate
    
    geo_spec_tck  = []
    mo_coeff_tck  = []
    mo_energy_tck = []
    mo_occ_tck    = []
    
    shape = self.geo_spec_all.shape
    for i in range(shape[1]):
      geo_spec_tck.append([])
      for j in range(shape[2]):
        geo_spec_tck[-1].append(interpolate.splrep(x,self.geo_spec_all[:,i,j],
                                k=k,**kwargs))
    
    for i_mo in range(len(self.mo_coeff_all)):
      mo_coeff_tck.append([])
      mo_energy_tck.append([])
      mo_occ_tck.append([])
      shape = self.mo_coeff_all[i_mo].shape    
      for i in range(shape[1]):
        mo_coeff_tck[-1].append([])
        mo_energy_tck[-1].append(interpolate.splrep(x,self.mo_energy_all[i_mo][:,i],
                                 k=k,**kwargs))
        mo_occ_tck[-1].append(interpolate.splrep(x,self.mo_occ_all[i_mo][:,i],
                              k=k,**kwargs))
        for j in range(shape[2]):
          mo_coeff_tck[-1][-1].append(interpolate.splrep(x,
                                      self.mo_coeff_all[i_mo][:,i,j],
                                      k=k,**kwargs))

  def interpolate_all(self,x,xnew,k=3,**kwargs):
    '''Interpolates a dataset y(x) to y(xnew) using B-Splines of order k.'''
    from scipy import interpolate
    
    shape = list(self.geo_spec_all.shape)
    shape[0] = len(xnew)
    tmp = numpy.zeros(shape)  
    for i in range(shape[1]):
      for j in range(shape[2]):
        tmp[:,i,j] = self.data_interp(x,self.geo_spec_all[:,i,j],xnew,k=k,**kwargs)
    self.geo_spec_all = numpy.copy(tmp)
    
    for i_mo in range(len(self.mo_coeff_all)):
      shape = list(self.mo_coeff_all[i_mo].shape)
      shape[0] = len(xnew)
      tmp = numpy.zeros(shape)  
      for i in range(shape[1]):
        for j in range(shape[2]):
          tmp[:,i,j] = self.data_interp(x,self.mo_coeff_all[i_mo][:,i,j],xnew,k=k,**kwargs)
      self.mo_coeff_all[i_mo] = numpy.copy(tmp)
      
      shape = list(self.mo_energy_all[i_mo].shape)
      shape[0] = len(xnew)
      tmp = numpy.zeros(shape)  
      for i in range(shape[1]):
        tmp[:,i] = self.data_interp(x,self.mo_energy_all[i_mo][:,i],xnew,k=k,**kwargs)
      self.mo_energy_all[i_mo] = numpy.copy(tmp)
      
      shape = list(self.mo_occ_all[i_mo].shape)
      shape[0] = len(xnew)
      tmp = numpy.zeros(shape)  
      for i in range(shape[1]):
        tmp[:,i] = self.data_interp(x,self.mo_occ_all[i_mo][:,i],xnew,k=k,**kwargs)
      self.mo_occ_all[i_mo] = numpy.copy(tmp)

  def plot(self,mo_matrix,symmetry='1',title='All',x_label='index',
           y_label='MO coefficients',output_format='png',
           plt_dir='Plots',ylim=None,thresh=0.1,x0=0,grid=True,x_grid=None,**kwargs):
    '''Plots all molecular orbital coefficients of one self.symmetry.'''
    import pylab as plt
    from matplotlib.ticker import MultipleLocator
    import os
    
    display('Plotting data of self.symmetry %s to %s/' % (symmetry,plt_dir))
    if not os.path.exists(plt_dir):
      os.makedirs(plt_dir)
    
    if numpy.ndim(mo_matrix) == 2:
      mo_matrix = mo_matrix[:,numpy.newaxis,:]
    
    shape = numpy.shape(mo_matrix)
    
    
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
          curves.append(ax.plot(X,Y, '.-' ,linewidth=1.5))
      
      
      plt.xlabel(x_label, fontsize=16);
      plt.ylabel(y_label, fontsize=16);
      plt.title('%s: %d.%s'%  (title,i+1,symmetry))
      plt.ylim(ylim)
      
      plt.tight_layout()
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

  def show_selected_mos(self,selected_mos,r0=0,steps=1,select_slice='xz',where=0.0,
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
        ao_list = ao_creator(self.geo_spec_all[rr],self.ao_spec)
        mo.append(mo_creator(ao_list,self.mo_coeff_all[self.sym[j]][rr,int(i)-1,numpy.newaxis])[0].reshape(tuple(npts)))
      
      f, pics = self.contour_mult_mo(xyz[k[0]],xyz[k[1]],mo,
                      xlabel=select_slice[0],ylabel=select_slice[1],
                      title='MO:%s' % mo_sel,r0=r0)
      for i,pic in enumerate(pics):
        pic.plot(self.geo_spec_all[rr,:,k[1]],self.geo_spec_all[rr,:,k[0]],nuclear_pos,
                 markersize=10,markeredgewidth=2)

  def contour_mult_mo(self,x,y,mo,xlabel='x',ylabel='y',title='',r0=0):
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
