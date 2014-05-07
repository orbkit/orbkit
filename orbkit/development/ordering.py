#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
'''Module for reading multiple output files of quantum chemical software with
a subsequent ordering of the molecular orbital coefficients 
(e.g. for the preparation of an interpolation). 

Example for the Execution
-------------------------

# Create a list of input filenames
fid_list = []
for i in range(101):
  fid_list.append('~/molden_files/%03d.molden' % (path,i))

init_display(name = 'MO_ordering')  # Specify a filename for the oklog file 
                                    # (We want a oklog file but we have no 
                                    # options.outputname here.)

# Start the read and the subsequent ordering routine.
main_order(fid_list,itype='molden')
'''

from copy import deepcopy
import numpy
import os

from orbkit.read import main_read
from orbkit.display import display,init_display

def read(fid_list,itype='molden'):
  '''Read a list of input files and construct arrays for the 
  ordering/interpolation routine.
  
  **Parameters:**
  
    fid_list : list of str
      List of input file names.
    itype : str, choices={'molden', 'gamess', 'gaussian.log', 'gaussian.fchk'}
        Specifies the type of the input files.
  
  **Returns:**
  
    - Geo_Spec, geo_info, ao_spec, MO_coeff, MO_energy, MO_occ, sym
    
    Geo_Spec : numpy.ndarray, shape=((len(fid_list),)+shape(geo_spec)) 
      Contains all molecular geometries. 
    geo_info, ao_spec : 
      See `Central Variables`_ for details.
    MO_coeff, MO_energy, MO_occ : list of numpy.ndarray
      Contains all molecular orbital coefficients, energies, and occupations,
      respectively.
    sym : dict
      Contains the molecular orbital symmetries and the corresponding position
      in MO_coeff, MO_energy, and MO_occ respectively.
  '''

  # Construct arrays for the ordering/interpolation routine
  Geo_Spec = []
  MO_Spec = []
  MO_coeff = []
  MO_energy = []
  MO_occ = []
  # geo_info and ao_info need to stay unchanged
  geo_old = []
  ao_old = []

  sym_list = {}
  n_ao = {}
  n_r = len(fid_list)

  for i,filename in enumerate(fid_list):
    geo_spec, geo_info, ao_spec, mo_spec = main_read(filename,
                                                itype=itype,
                                                all_mo=True)
    # Geo Section
    if i > 0 and (geo_old != geo_info):
      raise IOError('geo_info has changed!')
    else:
      geo_old = deepcopy(geo_info)
    Geo_Spec.append(geo_spec)
    # AO Section
    if (i > 0 and not
        numpy.alltrue([numpy.allclose(ao_old[j]['coeffs'],ao_spec[j]['coeffs']) 
                      for j in range(len(ao_old))]
                      )):
      raise IOError('ao_spec has changed!')
    else:
      ao_old = deepcopy(ao_spec)
    # MO Section
    sym = {}
    MO_Spec.append(mo_spec)
    for mo in mo_spec:
      key = mo['sym'].split('.')
      if key[1] not in sym.keys():
        sym[key[1]] = 0
        n_ao[key[1]] = len(mo_spec[0]['coeffs'])
      sym[key[1]] += 1
    
    for k,it in sym.iteritems():
      if k in sym_list:
        sym_list[k] = max(sym_list[k],it)
      else: 
        sym_list[k] = it
      
        
    
  Geo_Spec = numpy.array(Geo_Spec)

  # Presorting of the MOs according to their symmetry

  sym = []
  for k,it in sym_list.iteritems():
    sym.append((k,len(sym)))
    MO_coeff.append(numpy.zeros((n_r,it,n_ao[k])))
    MO_energy.append(numpy.zeros((n_r,it)))
    MO_occ.append(numpy.zeros((n_r,it)))

  sym = dict(sym)

  for i,spec in enumerate(MO_Spec):
    for j,mo in enumerate(spec):
      index,k = mo['sym'].split('.')
      
      index = int(index)-1
      
      MO_coeff[sym[k]][i,index,:] = mo['coeffs']
      MO_energy[sym[k]][i,index] = mo['energy']
      MO_occ[sym[k]][i,index] = mo['occ_num']
  
  return Geo_Spec, geo_info, ao_spec, MO_coeff, MO_energy, MO_occ, sym

def plot(mo_matrix,symmetry='1',plt_dir='Plots',title=''):
  '''Plots the all MO coefficients of one symmetry.'''
  import pylab as plt
  if not os.path.exists(plt_dir):
    os.makedirs(plt_dir)
  
  shape = numpy.shape(mo_matrix)
  
  ## Plot MO_coeffs
  for ii in range(shape[1]):
    fig=plt.figure()
    plt.rc('xtick', labelsize=16) 
    plt.rc('ytick', labelsize=16)
    curves=[]
    for ij in range(shape[2]):
      Y = mo_matrix[:,ii,ij]
      curves.append(plt.plot(Y, '-' ,linewidth=1.5))
    
    x_label='${\\sf index}$';y_label='MO coefficients';
    plt.xlabel(x_label, fontsize=16);
    plt.ylabel(y_label, fontsize=16);
    plt.title('%s %d.%s'%  (title,ii+1,symmetry))
    #plt.xlim([0, 10])
    #lgd = plt.legend(prop={'size':8})
    output_fid = '%d.%s.png' % (ii+1,symmetry)# sym[ii]
    fig.savefig('%s/%s' % (plt_dir, output_fid) ,format='png')


def variance_check(mo_matrix,shape,start=None,stop=1,step=-1,mu=1e-1):
  Dev = []
  var_check = []
  for ii in range(shape[1]):
    for ij in range(shape[2]):
      dev = 0.0
      for rr in range(shape[0])[start:stop:step]:
        dev += (mo_matrix[rr+step,ii,ij]-mo_matrix[rr,ii,ij])**2
      if dev >= mu:
        if ii not in var_check:
          var_check.append(ii)
        Dev.append(dev)
  display('\tElements with variance mu >= %.1e: %d' % (mu,len(var_check)))
  return var_check

def data_interp(x,y,xnew):
  
  tck = interpolate.splrep(x,y,s=0)
  ynew = interpolate.splev(xnew,tck,der=0)
  
  return ynew

# Ordering routine 
def order_mo(mo_matrix,index_list=None,backward=True,mu=1e-1):
  # MO coeffs can only interchange between MOs (in dimension 1) -
  # coeff[radius,mo,ao] --
  #--
  #[start:stop:step]
  shape = numpy.shape(mo_matrix)
  
  if index_list == None:
    index_list = numpy.ones((shape[0],shape[1]),dtype=int)
    index_list *= numpy.arange(shape[1],dtype=int)
  
  if backward:
    st = [None,1,-1]
  else:
    st = [None,-2,1]
  
  # First check using variance 
  var_check = variance_check(mo_matrix,shape,start=st[0],stop=st[1],step=st[2],mu=mu)
  
  # Sortiere nur Elemente aus var_check (dev >= 0.1)
  x=2 # Wir wollen den naechsten Wert bei der linearen Extrapolation    
   
  bFirst = [True for i in range(len(var_check))]
  for rr in range(shape[0])[st[0]:st[1]:st[2]]:
    for ii_index,ii in enumerate(var_check): 
      diff = numpy.zeros((len(var_check),shape[2]))
      if bFirst[ii_index]:
        dev = 0.0
        for ij in range(shape[2]):
          dev += (mo_matrix[rr+st[2],ii,ij]-mo_matrix[rr,ii,ij])**2
        if dev < mu:
          bFirst[ii_index] = False
      if not bFirst[ii_index]:
        Min = [0,0]
        AMin = [0,0]
        for ii_s,sign in [(0,-1),(1,+1)]:
          for ij in range(shape[2]):
            m = (sign*mo_matrix[rr+st[2],ii,ij]-mo_matrix[rr,ii,ij])
            epol = abs(m * x + mo_matrix[rr,ii,ij])
            ik_index = 0
            for ik in var_check:
              diff[ik_index,ij] = (mo_matrix[rr+2*st[2],ik,ij]-epol)**2
              ik_index += 1    
          Min[ii_s] = numpy.min(numpy.sum(diff,1))
          AMin[ii_s] = numpy.argmin(numpy.sum(diff,1))
        argmin = AMin[numpy.argmin(Min)]
        i_0 = var_check[ii_index]
        i_1 = var_check[argmin]     
        if i_0 != i_1:
          mo_matrix[rr+2*st[2],[i_0,i_1],:] = mo_matrix[rr+2*st[2],[i_1,i_0],:]
          index_list[rr+2*st[2],[i_0,i_1]] = index_list[rr+2*st[2],[i_1,i_0]]
  
  return mo_matrix, index_list
  

def order_pm(x,y,backward=True,mu=1e-1):
  if backward:
    st = [None,1,-1]
  else:
    st = [None,-2,1]
  
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
    bFirst = [True for i in range(shape[1])]
    for rr in range(shape[0])[st[0]:st[1]:st[2]]:
      diff = numpy.zeros(2)
      for ii in range(shape[1]):
        if bFirst:
          dev = (y[rr+st[2],ii]-y[rr,ii])**2
          if dev < mu:
            bFirst[ii] = False
        if not bFirst[ii]:
          m = (y[rr+st[2],ii]-y[rr,ii])/(x[rr+st[2]]-x[rr])
          epol = m * (x[rr+2*st[2]]-x[rr]) + y[rr,ii]
          for ii_d in range(2):
            diff[ii_d] += (((-1)**ii_d * y[rr+2*st[2],ii])-epol)**2
      if numpy.argmin(numpy.abs(diff)) == 1:
        y[rr+2*st[2],:] = -y[rr+2*st[2],:]    
  else:
    display('Function order_pm only works for vectors and 2D matrices')
  return y

def main_order(fid_list,itype='molden',
         plt_dir=['/tmp/Plots/original','/tmp/Plots/ordered']):
  '''Perform all read, ordering(, and interpolation) routines.
  
  **Paramerters:**
  
  fid_list : list of str
    List of input file names.
  itype : str, choices={'molden', 'gamess', 'gaussian.log', 'gaussian.fchk'}
    Specifies the type of the input files.
  plt_dir : {None or list of str}, optional
    If not None, plots the molecular orbital coefficients before and after the
    ordering routine the folder specified.
  '''
  global Geo_Spec, geo_info, ao_spec, MO_coeff, MO_energy, MO_occ, sym
  # Read all input files
  Geo_Spec, geo_info, ao_spec, MO_coeff, MO_energy, MO_occ, sym = read(fid_list)

  radius = range(len(Geo_Spec)) #: We assume an equally spaced grid
  for s,ii_s in sym.iteritems():
    display('Starting ordering of MOs of symmetry %s' % s)
    
    if plt_dir is not None:
      display('\tPlotting molecular orbitals before ordering')
      plot(MO_coeff[ii_s],symmetry=s,plt_dir=plt_dir[0])
    
    shape = numpy.shape(MO_coeff[ii_s])
    
    # Save sign of MOs 
    display('\tSaving the sign of the molecular orbitals')
    coeff_abs = numpy.abs(MO_coeff[ii_s])
    coeff_phase = numpy.zeros(shape, numpy.bool)
    
    for rr in range(shape[0]):
      for ii in range(shape[1]):
        for ij in range(shape[2]):
          coeff_phase[rr,ii,ij] = MO_coeff[ii_s][rr,ii,ij] < 0
    
    mo_matrix = coeff_abs
    display('\tOdering backward')
    mo_matrix, index_list = order_mo(mo_matrix,index_list=None,backward=True,mu=1e-1)
    display('\tOdering forward')
    mo_matrix, index_list = order_mo(mo_matrix,index_list=index_list,backward=False,mu=1e-1)
    
    for rr in range(shape[0]):
      index = index_list[rr,:]
      MO_coeff[ii_s][rr,:,:] = MO_coeff[ii_s][rr,index,:]
      MO_energy[ii_s][rr,:] = MO_energy[ii_s][rr,index]
      #MO_index[ii_s][rr,:] = MO_index[ii_s][rr,index]
      MO_occ[ii_s][rr,:] = MO_occ[ii_s][rr,index]
    
    display('\tOdering the signs of the molecular orbitals')
    for ii in range(shape[1]):
      MO_coeff[ii_s][:,ii,:] = order_pm(radius,MO_coeff[ii_s][:,ii,:],backward=True)
      MO_coeff[ii_s][:,ii,:] = order_pm(radius,MO_coeff[ii_s][:,ii,:],backward=False)
    
    if plt_dir is not None:
      display('\tPlotting ordered molecular orbitals')
      plot(MO_coeff[ii_s],symmetry=s,plt_dir=plt_dir[1],title='Ordered:')


