# -*- coding: iso-8859-1 -*-

lgpl = '''
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

lgpl_short = '''This is orbkit.
Copyright (C) 2014 Gunter Hermann, Vincent Pohl, and Axel Schild. 
This program comes with ABSOLUTELY NO WARRANTY. 
This is free software, and you are welcome to redistribute it 
under certain conditions. Type '-l' for details.
'''

# Import general modules
import copy
import os
import optparse
import re
import string
import sys
import time

import numpy
from scipy import weave
from scipy import integrate

import multiprocessing as mp
from multiprocessing import Process
from multiprocessing import Pool

# Import orbkit modules
import orbkit_grid as grid
import orbkit_output 

def l_creator_atom(at_pos=None,exp_list=None,coeff_list=None,
	      pnum_list=None,x=None,y=None,z=None,N=None):
  ''' FUNCTION l_creator_atom
  Calculate all contracted atomic orbitals of ONE selected atom 
  
  We use scipy.weave.inline to run C++ Code within the python environment, see
  http://docs.scipy.org/doc/numpy/user/c-info.python-as-glue.html#inline-c-code
  '''
  
  if x == None: x = grid.x
  if y == None: y = grid.y
  if z == None: z = grid.z
  if N == None: N = numpy.array(grid.N_)
  
  
  #--- All AOs from one atom ---
  if exp_list == None or coeff_list == None or numpy.shape(at_pos) != (3,):
    orbkit_output.display(
      "If all AOs for one atom are requested, the variables " +
      "exp_list and coeff_list and pnum_list \nhave to be specify explicitly.\n" +
      "Moreover, geo_spec has to contain only one set of [x,y,z] coordinates:\n\t" +
      "l_creator(at_pos=[x,y,z],exp_list=exp,coeff_list=coeff,pnum_list=pnum)")
    return 255
  else:
    at_pos = numpy.array(at_pos)
    x = x - at_pos[0]
    y = y - at_pos[1]
    z = z - at_pos[2]
    exp_list = numpy.array(exp_list)
    coeff_list = numpy.array(coeff_list)
    pnum = numpy.shape(coeff_list)[0]
    pnum_list = numpy.array(pnum_list)
    which_ao = numpy.zeros((pnum,))

    ii_p = 0
    for ii in range(len(pnum_list)):
      for ij in range(pnum_list[ii]):
	which_ao[ii_p] = ii
	ii_p += 1
  
  ao_list = numpy.zeros((len(exp_list),len(x),len(y),len(z)))  
  
  ao_code = '''
  double Norm[pnum][Nao_list[0]];
  double X, Y, Z;
  int il; 
  int lx[Nao_list[0]], ly[Nao_list[0]], lz[Nao_list[0]];
  double rr, ao_l0, ao_xyz;
  
  for (int ii=0; ii<pnum; ii++)
  {
    il = which_ao[ii];
    lx[il] = EXP_LIST2(il,0);
    ly[il] = EXP_LIST2(il,1);
    lz[il] = EXP_LIST2(il,2);
    Norm[ii][il] = ao_norm(lx[il],ly[il],lz[il],&COEFF_LIST2(ii,0));
  }
  
  for (int i=0; i<Nao_list[1]; i++)
  {
    for (int j=0; j<Nao_list[2]; j++)
    {
      for (int k=0; k<Nao_list[3]; k++)
      {
	X = x[i];
	Y = y[j];
	Z = z[k];
	
	rr = pow(X,2)+pow(Y,2)+pow(Z,2);
	
	
	for (int ii=0; ii<pnum; ii++)
	{
	  ao_l0 = COEFF_LIST2(ii,1) * exp(-COEFF_LIST2(ii,0) * rr);
	  
	  il = which_ao[ii];
	  ao_xyz = xyz(X, Y, Z, lx[il], ly[il], lz[il]);
	  AO_LIST4(il,i,j,k) += Norm[ii][il] * ao_xyz * ao_l0;

	}
	
      }
    }
  }
  '''

  weave.inline(ao_code, ['x','y','z','N','exp_list','coeff_list',
			 'pnum','which_ao','ao_list'], 
	      support_code = cSupportCode.norm + cSupportCode.xyz,verbose = 1)
  
  return ao_list
  #--- l_creator_atom ---


def delta_l_creator_atom(at_pos=None,exp_list=None,coeff_list=None,
		    pnum_list=None,drv='x',x=None,y=None,z=None,N=None):
  ''' FUNCTION delta_l_creator_atom
  Calculate the derivative of ALL contracted atomic orbitals
  of ONE selected atom 
  with respect to a specific variable (e.g. drv = 'x' or drv = 0)
  
  We use scipy.weave.inline to run C++ Code within the python environment.
  http://docs.scipy.org/doc/numpy/user/c-info.python-as-glue.html#inline-c-code
  '''
  
  if x == None: x = grid.x
  if y == None: y = grid.y
  if z == None: z = grid.z
  if N == None: N = numpy.array(grid.N_)
  
  #--- With respect to which variable the derivative shall be computed? --- 
  if not isinstance(drv, (int, long)):
    drv = 'xyz'.find(drv)
  if drv == -1:		# Was the selection valid? If not drv='x'
    drv = 0
    orbkit_output.display("The selection of the derivative variable was not valid!" +
		    " (drv = 'x' or 'y' or 'z')")
    orbkit_output.display("Calculating the derivative with respect to x...")	


  #--- All AOs from one atom
  if exp_list == None or coeff_list == None or numpy.shape(at_pos) != (3,):
    orbkit_output.display(
      "If all AOs for one atom are requested, the variables " +
      "exp_list and coeff_list and pnum_list \nhave to be specify explicitly.\n" +
      "Moreover, geo_spec has to contain only one set of [x,y,z] coordinates:\n\t" +
      "l_creator(at_pos=[x,y,z],exp_list=exp,coeff_list=coeff,pnum_list=pnum)")
    return 255
  else:
    at_pos = numpy.array(at_pos)
    x = x - at_pos[0]
    y = y - at_pos[1]
    z = z - at_pos[2]
    exp_list = numpy.array(exp_list)
    coeff_list = numpy.array(coeff_list)
    pnum = numpy.shape(coeff_list)[0]
    pnum_list = numpy.array(pnum_list)
    which_ao = numpy.zeros((pnum,))

    ii_p = 0
    for ii in range(len(pnum_list)):
      for ij in range(pnum_list[ii]):
	which_ao[ii_p] = ii
	ii_p += 1

  ao_list = numpy.zeros((len(exp_list),len(x),len(y),len(z)))  


  ao_code = '''
  double Norm[pnum][Nao_list[0]];
  double X, Y, Z;
  int il; 
  int lx[Nao_list[0]], ly[Nao_list[0]], lz[Nao_list[0]];
  double rr, ao_l0, ao_xyz;
  
  for (int ii=0; ii<pnum; ii++)
  {
    il = which_ao[ii];
    lx[il] = EXP_LIST2(il,0);
    ly[il] = EXP_LIST2(il,1);
    lz[il] = EXP_LIST2(il,2);
    Norm[ii][il] = ao_norm(lx[il],ly[il],lz[il],&COEFF_LIST2(ii,0));
  }
  
  for (int i=0; i<Nao_list[1]; i++)
  {
    for (int j=0; j<Nao_list[2]; j++)
    {
      for (int k=0; k<Nao_list[3]; k++)
      {
	X = x[i];
	Y = y[j];
	Z = z[k];
	
	rr = pow(X,2)+pow(Y,2)+pow(Z,2);
	
	for (int ii=0; ii<pnum; ii++)
	{
	  ao_l0 = COEFF_LIST2(ii,1) * exp(-COEFF_LIST2(ii,0) * rr);
	  
	  il = which_ao[ii];
	  switch(drv)
	  {
	    case 0:
	    {
	      if (lx[il] == 0)
	      {
		ao_xyz = - 2 * COEFF_LIST2(ii,0) * xyz(X, Y, Z, lx[il]+1, ly[il], lz[il]);
	      }
	      else
	      {
		ao_xyz = lx[il] * xyz(X, Y, Z, lx[il]-1, ly[il], lz[il])  - 2 * COEFF_LIST2(ii,0) * xyz(X, Y, Z, lx[il]+1, ly[il], lz[il]);	    
	      }
	    } break;
	    
	    case 1:
	    {
	      if (ly[il] == 0)
	      {
		ao_xyz = - 2 * COEFF_LIST2(ii,0) * xyz(X, Y, Z, lx[il], ly[il]+1, lz[il]);
	      }
	      else
	      {
		ao_xyz = ly[il] * xyz(X, Y, Z, lx[il], ly[il]-1, lz[il])  - 2 * COEFF_LIST2(ii,0) * xyz(X, Y, Z, lx[il], ly[il]+1, lz[il]);	    
	      }
	    } break;
	    
	    case 2:
	    {
	      if (lz[il] == 0)
	      {
		ao_xyz = - 2 * COEFF_LIST2(ii,0) * xyz(X, Y, Z, lx[il], ly[il], lz[il]+1);
	      }
	      else
	      {
		ao_xyz = lz[il] * xyz(X, Y, Z, lx[il], ly[il], lz[il]-1)  - 2 * COEFF_LIST2(ii,0) * xyz(X, Y, Z, lx[il], ly[il], lz[il]+1);	    
	      }
	    } break;
	    
	    default:
	    {
	      std::cout << "False statement for derivative variable!" << std::endl;
	    }
	  }
	  AO_LIST4(il,i,j,k) += Norm[ii][il] * ao_xyz * ao_l0;
	
	}
	
      }
    }
  }
  
  '''

  weave.inline(ao_code, ['x','y','z','N','at_pos','exp_list','coeff_list',
			 'pnum','which_ao','ao_list','drv'], 
	      support_code = cSupportCode.norm + cSupportCode.xyz,verbose = 1)

  return ao_list
  #--- delta_l_creator_atom ---


def mo_creator_matrop(ao_list,mo_coeff,AO_Nr,x=None,y=None,z=None,vector=False):
  '''FUNCTION mo_creator
  Calculate the molecular orbitals for MO-coefficient matrix

  Parameters
  ----------
  ao_list : numpy array
	    containing all atomic orbitals on the specified grid
  mo_coeff: python list of numpy arrays (N_mo x N_ao)
	    containing all molecular orbital coefficients
  AO_Nr: python list of python tuples
	    containing information which element of mo_coeff corresponds to
	    which atomic orbital
  '''
  
  if x == None: x = grid.x
  if y == None: y = grid.y
  if z == None: z = grid.z
  
  if not vector:
    N = tuple(grid.N_)
  else:
    if len(x) != len(y) or len(x) != len(z):
      orbkit_output.display("Dimensions of x-, y-, and z- coordinate differ!")
      return 0
    else:
      N = (len(x),)
  
  mo_list = []
  i = 0
  for ii_s in range(len(mo_coeff)):
    shape = numpy.shape(mo_coeff[ii_s])
    for ii in range(shape[0]):
      mo_list.append(numpy.zeros(N))
      for jj in range(shape[1]):
        for jj_ao in range(len(ao_list)):
	  if AO_Nr[jj_ao] == (jj+1,ii_s+1):
	    mo_list[i] += mo_coeff[ii_s][ii,jj] * ao_list[jj_ao]
      i += 1
  return mo_list
  #--- mo_creator_matrop ---

###################################################
def ao_creator2(geo_spec,ao_spec,fid='tmp',x=None,y=None,z=None):
  '''FUNCTION ao_creator2 -> Will be removed/changed in a future release!
  The atomic orbital creator for save_ao
  
  save_ao: Save all AOs to disk and reload for every MO calculation
  '''
  if x == None: x = grid.x
  if y == None: y = grid.y
  if z == None: z = grid.z
  
  HDF5_f = h5py.File('%(f)s_ao.h5' % {'f': fid},'w')
  shape = (len(x),len(y),len(z))
  ii_ao = 0
  # Generalized AO creator
  for ii in range(len(ao_spec)):
    ao_list = l_creator(geo_spec,ao_spec,ii,x=x,y=y,z=z)    
    for kk in range(len(ao_list)):
      #--- Save the atomic orbitals as hdf5 file ---
      dset_id = 'ao:%(i)03d' % {'i': ii_ao}
      dset = HDF5_f.create_dataset(dset_id, shape, data=ao_list[kk])
      #with h5py.File(filename) as f: 
	#f['/ao'] = ao_list[kk]  
      ii_ao += 1
  HDF5_f.close()
  #--- Store and return the number of atomic orbitals ---
  ao_num = ii_ao
  
  return ao_num
  #--- ao_creator2 ---

def mo_creator2(ao_num,mo_spec,fid='tmp',x=None,y=None,z=None):
  '''FUNCTION mo_creator2 -> Will be removed/changed in a future release!
  The molecular orbital creator for save_ao
  
  save_ao: Save all AOs to disk and reload for every MO calculation
  '''
  if x == None: x = grid.x
  if y == None: y = grid.y
  if z == None: z = grid.z
  
  HDF5_f = h5py.File('%(f)s_ao.h5' % {'f': fid},'r')
  mo_list = []
  for ii in range(len(mo_spec)):
    mo_list.append(numpy.zeros((len(x),len(y),len(z))));
  
  mo_count = 0
  for ii_ao in range(ao_num): 
    dset_id = 'ao:%(i)03d' % {'i': ii_ao}
    ao_list = HDF5_f[dset_id]
    #filename = '%(f)s_ao_%(i)03d.h5' % {'f': fid, 'i': ii_ao}
    #with h5py.File(filename) as f: 
      #ao_list = f['/ao'].value 
    for ii in range(len(mo_spec)):
      mo_list[ii] += mo_spec[ii]['coeffs'][mo_count] * ao_list
    mo_count += 1
  HDF5_f.close()
  return mo_list
  #--- mo_creator2 ---

###################################################  
def mo_creator3(geo_spec,ao_spec,mo_spec,x=None,y=None,z=None):
  '''FUNCTION mo_creator3
  The molecular orbital creator for discard_ao
  
  discard_ao: Calculate MOs without keeping the AOs in the RAM
	      Dedicated for huge systems to save RAM
  '''
  if x == None: x = grid.x
  if y == None: y = grid.y
  if z == None: z = grid.z
  
  mo_list = []
  for ii in range(len(mo_spec)):
    mo_list.append(numpy.zeros((len(x),len(y),len(z))))

  mo_count = 0 
  #--- Generalized AO creator ---
  for ii in range(len(ao_spec)):
    ao_list = l_creator(geo_spec,ao_spec,ii,x=x,y=y,z=z)
    
    for kk in range(len(ao_list)):
      for jj in range(len(mo_spec)):
	mo_list[jj] += mo_spec[jj]['coeffs'][mo_count] * ao_list[kk]
      mo_count += 1
              
  return mo_list
  #--- mo_creator3 ---


'''

    if options.save_ao:
      fid = '%(f)s_x_%(i)0.4f' % {'f': options.outputname, 'i': x} 
      ao_count = ao_creator2(geo_spec,ao_spec,fid=fid,x=x,y=y,z=z,vector=vector)
      mo_list = mo_creator2(ao_count,mo_spec,fid=fid,x=x,y=y,z=z,vector=vector)
    elif options.discard_ao:
      mo_list = mo_creator3(geo_spec,ao_spec,mo_spec,x=x,y=y,z=z,vector=vector)
    else:
'''




################

####
def slice_rho_vector(zz):
  #--- FUNCTION slice_rho ------------------------------------------------------
  #--- Calculate the density or the derivative of the density with respect -----
  #--- to Spec['Derivative'] for one slice in x-direction ----------------------
  #-----------------------------------------------------------------------------
  try:
    #--- All desired information is stored in the Global variable Spec ---
    geo_spec = Spec['geo_spec']
    ao_spec = Spec['ao_spec']
    mo_spec = Spec['mo_spec']
    drv = Spec['Derivative']
    
    # Set up Grid
    x = grid.x[zz[0]:zz[1]]
    y = grid.y[zz[0]:zz[1]]
    z = grid.z[zz[0]:zz[1]]
    
    
    #--- Calculate the MOs and AOs for this slice ---
    if options.save_ao:
      fid = '%(f)s_sl_%(i)03d' % {'f': options.outputname, 'i': zz[0]} 
      ao_count = ao_creator2(geo_spec,ao_spec,fid=fid,x=x,y=y,z=z,vector=True)
      mo_list = mo_creator2(ao_count,mo_spec,fid=fid,x=x,y=y,z=z,vector=True)
    elif options.discard_ao:
      mo_list = mo_creator3(geo_spec,ao_spec,mo_spec,x=x,y=y,z=z,vector=True)
    else:
      ao_list = ao_creator(geo_spec,ao_spec,x=x,y=y,z=z,vector=True)
      mo_list = mo_creator(ao_list,mo_spec,x=x,y=y,z=z,vector=True)

    #--- Initialize a numpy array for the density ---
    rho = numpy.zeros((len(x),))
    
    #--- Initialize a numpy array for the norm of the MOs ---
    mo_norm = numpy.zeros((len(mo_list)))
    
    #--- Calculate the density and the norm ---
    for ii_mo in range(len(mo_list)): 
      mo_norm[ii_mo] = numpy.sum(numpy.square(mo_list[ii_mo]))
      rho += mo_spec[ii_mo]['occ_num'] * numpy.square(numpy.abs(mo_list[ii_mo]))
      
    if drv == None:
      #--- Return the density and the norm ---
      return rho, mo_norm
    else:
      #--- Initialize a numpy array for the derivative of the density ---
      delta_rho = numpy.zeros((len(drv),len(x)))
      for ii_d, d in enumerate(drv):
	#--- Calculate the derivatives of the AOs and MOs for this slice ---
	delta_ao_list = delta_ao_creator(geo_spec,ao_spec,drv=d,x=x,y=y,z=z,vector=True)
	delta_mo_list = mo_creator(delta_ao_list,mo_spec,x=x,y=y,z=z,vector=True)
	
	#--- Calculate the derivative of the density
	for ii_mo in range(len(delta_mo_list)): 
	  delta_rho[ii_d,:] += (mo_spec[ii_mo]['occ_num'] * 
			    2 * delta_mo_list[ii_mo]*mo_list[ii_mo])
      #--- Return the derivative of the density ---
      return rho, mo_norm, delta_rho
  except KeyboardInterrupt:
    #--- Catch keybord interrupt signal to prevent a hangup of the worker processes ---
    return 0
  #--- slice_rho ---

def rho_compute_vector(geo_spec,geo_info,ao_spec,mo_spec,delta_slice=1e4):
  '''Function for computing the density for a vector grid
  
  The 3D grid is divided into Slices, and the computational tasks
  are distributed to the worker processes
  '''
  
  #--- Specify the global variable containing all desired information needed ---
  #--- by the function slice_rho ---
  global Spec
  
  Spec = {'geo_spec': geo_spec, 'ao_spec': ao_spec, 
	      'mo_spec': mo_spec, 'Derivative': None}
  
  if len(grid.x) != len(grid.y) or len(grid.x) != len(grid.z):
      orbkit_output.display("Dimensions of x-, y-, and z- coordinate differ!")
      sys.exit(99)
  
  N = len(grid.x)
  sNum = int(numpy.floor(N/delta_slice)+1)
  #--- Make slices ---
  
  rho = numpy.zeros((N,))
  
  #--- Print information regarding the density calculation ---
  orbkit_output.display("\nStarting the density calculation...")
  
  #--- Initialize some additional user information ---
  status_old = 0
  s_old = 0
  t = [time.time()]
  
  #--- Initialize an array to store the norm of the MOs ---
  mo_norm = numpy.zeros((len(mo_spec)))
  
  #--- Start the worker processes --
  pool = Pool(processes=options.numproc)
  
  #--- Write the slices in x to an array zz ---
  zz=[]
  i = 0 
  
  for s in range(sNum):
    if i+delta_slice >= N:
      zz.append((numpy.array([i,N+1],dtype=int)))      
    else:
      zz.append((numpy.array([i,i+delta_slice],dtype=int)))
    i += delta_slice 
  
  #--- Compute the density slice by slice ---
  it = pool.imap(slice_rho_vector, zz)
  for s in range(sNum):
    rho[zz[s][0]:zz[s][1]], norm = it.next()
    mo_norm += norm
    
    #--- Print out the progress of the computation ---
    status = round(s*100/float(sNum))
    if not status % 20 and status != status_old: 
      t.append(time.time())
      orbkit_output.display("\tFinished %(f)d%% (%(s)d slices in %(t).3fs)" 
		      % {'f': status,
			 's': s + 1 - s_old,
			 't': t[-1]-t[-2]})
      status_old = status
      s_old = s + 1
  
  #--- Close the worker processes --
  pool.close()
  
  #--- Print the norm of the MOs ---
  orbkit_output.display("\nNorm of the MOs:")
  for ii_mo in range(len(mo_norm)):
    orbkit_output.display("\t%(m).6f\tMO %(n)s" 
		    % {'m':mo_norm[ii_mo]*grid.d3r, 'n':mo_spec[ii_mo]['sym']})
  
  #--- Print the number of electrons ---
  orbkit_output.display("We have " + str(numpy.sum(rho)*grid.d3r) + " electrons.")
  
  return rho
  #--- rho_compute ---



def delta_rho_compute(geo_spec, geo_info, ao_spec, mo_spec, drv=[0,1,2]):
  '''FUNCTION delta_rho_compute
  Function for computing the derivative of the density
  with respect to a specific variable (e.g. drv = ['x','y'] or drv = 0)
  
  The 3D grid is divided into Slices, and the computational tasks
  are distributed to the worker processes
  '''
  try:
    drv = list(drv)
  except TypeError: 
    drv = [drv]
  
  #--- Specify the global variable containing all desired information needed ---
  #--- by the function slice_rho ---
  global Spec
  Spec = {'geo_spec': geo_spec, 'ao_spec': ao_spec, 
	      'mo_spec': mo_spec, 'Derivative': drv}
  
  #--- Save the complete initial grid ---
  x = grid.x
  y = grid.y
  z = grid.z
  N = copy.deepcopy(grid.N_)
  Min = copy.deepcopy(grid.min_)
  Max = copy.deepcopy(grid.max_)
  
  #--- Make slices ---
  rho = numpy.zeros((len(grid.x),len(grid.y),len(grid.z)))
  delta_rho = numpy.zeros((len(drv),len(grid.x),len(grid.y),len(grid.z)))
  sDim = 0
  sNum = grid.N_[sDim]
  grid.N_[sDim] = 1
  
  #--- The number of worker processes is capped to the number of ---
  #--- grid points in x-direction. --- 
  if options.numproc > sNum: options.numproc = sNum
  
  #--- Print information regarding the density calculation ---
  orbkit_output.display("\nStarting the density calculation...")
  orbkit_output.display("The grid has been seperated in 2d slices and the")
  orbkit_output.display("calculation will be carried out with %d subprocesses." 
		  % options.numproc)
  orbkit_output.display("\nThere are %d contracted AOs and %d MOs to be calculated."
		  % (len(mo_spec[0]['coeffs']), len(mo_spec)))
  
  #--- Initialize some additional user information ---
  status_old = 0
  s_old = 0
  t = [time.time()]
  
  #--- Initialize an array to store the norm of the MOs ---
  mo_norm = numpy.zeros((len(mo_spec)))
    
  #--- Start the worker processes --
  pool = Pool(processes=options.numproc)
  
  #--- Write the slices in x to an array zz ---
  zz=[]
  for s in range(sNum):
    zz.append((numpy.array([x[s]])))
  
  #--- Compute the derivative of the density slice by slice ---

  it = pool.imap(slice_rho, zz)
  for s in range(sNum):
    rho[s,:,:], norm, d_rho = it.next()
    for ii_d in range(len(drv)):
      delta_rho[ii_d,s,:,:] = d_rho[ii_d,:,:,:]
    mo_norm += norm
    
    #--- Print out the progress of the computation ---
    status = round(s*100/float(sNum))
    if not status % 20 and status != status_old: 
      t.append(time.time())
      orbkit_output.display('\tFinished %(f)d%% (%(s)d slices in %(t).3fs)'
		      % {'f': status,
			 's': s + 1 - s_old,
			 't': t[-1]-t[-2]})
      status_old = status
      s_old = s + 1
  
  #--- Close the worker processes --
  pool.close()
  
  #--- Restore the complete initial grid ---
  grid.N_ = N
  grid.min_ = Min
  grid.max_ = Max
  grid.grid_init()

  #--- Print the norm of the MOs ---
  orbkit_output.display('\nNorm of the MOs:')
  for ii_mo in range(len(mo_norm)):
    orbkit_output.display('\t%(m).6f\tMO %(n)s' 
		    % {'m':mo_norm[ii_mo]*grid.d3r, 'n':mo_spec[ii_mo]['sym']})
  
  #--- Print the number of electrons ---
  orbkit_output.display('We have %f electrons.' % integration(rho))
  
  return rho, delta_rho
  #--- delta_rho_compute ---

def delta_rho_compute_no_slice(geo_spec, geo_info, ao_spec, mo_spec, drv=[0,1,2]):
  '''FUNCTION delta_rho_compute_no_slice
  Function for computing the derivative of the density
  with respect to a specific variable (e.g. drv = ['x','y'] or drv = 0)
  without slicing the grid
  '''
   try:
    drv = list(drv)
  except TypeError: 
    drv = [drv]
  
  #--- Initialize a numpy array for the derivative of the density ---
  ao_list = ao_creator(geo_spec,ao_spec)
  mo_list = mo_creator(ao_list,mo_spec)
  
  #--- Initialize a numpy array for the density
  rho = numpy.zeros((len(grid.x),len(grid.y),len(grid.z)))
  
  #--- Calculate the density and the norm
  for ii_mo in range(len(mo_list)): 
    rho += numpy.square(numpy.abs(mo_list[ii_mo])) * mo_spec[ii_mo]['occ_num']
  
  #--- Print the norm of the MOs ---
  orbkit_output.display('Norm of the MOs:')
  for ii_mo in range(len(mo_list)): 
    orbkit_output.display('\t%(m).6f\tMO %(n)s' % 
		    {'m':numpy.sum(mo_list[ii_mo]**2)*grid.d3r, 
		     'n':mo_spec[ii_mo]['sym']})
  
  #--- Print the number of electrons ---
  orbkit_output.display('We have %f electrons.' % integration(rho))
  
  #--- Print information ---
  orbkit_output.display('\nCalculating the derivative of the density...')
  delta_rho = numpy.zeros((len(drv),len(grid.x),len(grid.y),len(grid.z)))
  
  #--- Loop over spatial directions ---
  for d_ii in range(len(drv)):
    orbkit_output.display('\t...with respect to %s' % drv[d_ii])
    #--- Calculate the derivatives of the AOs and MOs ---
    delta_ao_list = delta_ao_creator(geo_spec,ao_spec,drv=drv[d_ii])
    delta_mo_list = mo_creator(delta_ao_list,mo_spec)
    
    
    #--- Calculate the derivative of the density
    for ii_mo in range(len(delta_mo_list)): 
      delta_rho[d_ii,:,:,:] += (mo_spec[ii_mo]['occ_num'] * 
				2 * delta_mo_list[ii_mo]*mo_list[ii_mo])
  
  return rho, delta_rho
  #--- delta_rho_compute_no_slice ---

