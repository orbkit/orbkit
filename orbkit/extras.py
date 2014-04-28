# -*- coding: iso-8859-1 -*-
'''Module for all additional features of orbkit.'''
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
import os
import sys

import numpy
from scipy import integrate

# Import orbkit modules
from orbkit import core,grid,output,options
from orbkit.display import display

def mo_select(mo_spec, fid_mo_list):
  '''Selects  molecular orbitals from an external file.

  **Parameters:**
   
    mo_spec :        
      See `Central Variables`_ for details.
    fid_mo_list : str
      Specifies the filename of the molecular orbitals list. 
      If fid_mo_list is 'all_mo', creates a list containing all molecular orbitals.

  **Supported Formats:**
  
    Integer List:
    
      .. literalinclude:: ../examples/MO_List_int.tab
            :language: bash
    
    List with Symmetry Labels:
    
      .. literalinclude:: ../examples/MO_List.tab
            :language: bash

  **Returns:**
  
    Dictionary with following Members:
      :mo: - List of molecular orbital labels.
      :mo_ii: - List of molecular orbital indices.
      :mo_spec: - Selected elements of mo_spec. See `Central Variables`_ for details.
      :mo_in_file: - List of molecular orbital labels within the fid_mo_list file.
      :sym_select: - If True, symmetry labels have been used. 
      
  '''
  mo_in_file = []
  all_mo = []
  selected_mo_spec = []
  selected_mo_ii = [] 
  sym_select = False
  
  if fid_mo_list.lower() == 'all_mo':
    selected_mo = numpy.array(numpy.arange(len(mo_spec))+1, dtype=numpy.str)
    mo_in_file = [selected_mo]
    selected_mo_spec = mo_spec
  else:  
    try:
      fid=open(fid_mo_list,'r')
      flines = fid.readlines()
      fid.close()
      for line in flines:
        integer = line.split()
        mo_in_file.append(integer)
        all_mo = all_mo + integer
    except:
      raise IOError('The selected mo-list (%(m)s) is not valid!' % 
                    {'m': fid_mo_list} + '\ne.g.\n\t1\t3\n\t2\t7\t9\n')
    
    #--- Check if the mos are specified ---
    #--- by symmetry (e.g. 1.1 in MOLPRO nomenclature) or ---
    #--- by the number in the molden file (e.g. 1) ---
    selected_mo=list(set(all_mo))
    try: # Try to convert selections into integer
      for i in selected_mo: 
        int(i)
    except ValueError:
      sym_select = True
      # Add '.' after molecular orbital number if missing
      for i in range(len(selected_mo)):
        if not '.' in selected_mo[i]:
          from re import search
          a = search(r'\d+', selected_mo[i]).group()
          if a == selected_mo[i]:
            selected_mo[i] = '%s.1' % a
          else:
            selected_mo[i] = selected_mo[i].replace(a, '%s.' % a)
    
    if sym_select:
      for k in range(len(mo_spec)):
        if mo_spec[k]['sym'] in selected_mo:
          selected_mo_spec.append(mo_spec[k])
          selected_mo_ii.append(mo_spec[k]['sym'])
      selected_mo_ii = numpy.array(selected_mo_ii)
    else:
      selected_mo = map(int, selected_mo)            
      selected_mo.sort()
      selected_mo = map(str, selected_mo)
      for k in range(len(mo_spec)):
        if str(k+1) in selected_mo: selected_mo_spec.append(mo_spec[k])
  
  return {'mo': selected_mo, 'mo_ii': selected_mo_ii,
          'mo_spec': selected_mo_spec, 
          'mo_in_file': mo_in_file, 'sym_select': sym_select}

def calc_mo(geo_spec, geo_info, ao_spec, mo_spec, fid_mo_list, 
            drv=None, vector=None, otype=None):
  '''Calculates and saves the selected molecular orbitals or the derivatives thereof.

  **Parameters:**
   
    geo_spec, geo_info, ao_spec, mo_spec :        
      See `Central Variables`_ for details.
    fid_mo_list : str
      Specifies the filename of the molecular orbitals list. 
      If fid_mo_list is 'all_mo', creates a list containing all molecular orbitals.
    otype : str or list of str, optional
      Specifies output file type. See :data:`otypes` for details.
    drv : int or string, {None, 'x', 'y', 'z', 0, 1, 2}, optional
      If not None, a derivative calculation of the atomic orbitals 
      is requested.
    vector : None or int, optional
      If not None, performs the computations on a vectorized grid, i.e., 
      with x, y, and z as vectors.
    
    
  **Returns:**
    mo_list :        
      See `Central Variables`_ for details.
    mo : dict 
      Contains information of the selected molecular orbitals and has following Members:
        :mo: - List of molecular orbital labels.
        :mo_ii: - List of molecular orbital indices.
        :mo_spec: - Selected elements of mo_spec. See `Central Variables`_ for details.
        :mo_in_file: - List of molecular orbital labels within the fid_mo_list file.
        :sym_select: - If True, symmetry labels have been used. 
      
  '''

  mo = mo_select(mo_spec, fid_mo_list)
  
  #--- Calculate the AOs and MOs ---
  mo_list = core.rho_compute(geo_spec,ao_spec,mo['mo_spec'],
                             calc_mo=True,drv=drv,vector=vector,
                             numproc=options.numproc)
  
  fid = '%s_MO' % (options.outputname)
  
  if not options.no_output:
    if 'h5' in otype:    
      output.main_output(mo_list,geo_info,geo_spec,data_id='MO',
                    outputname=fid,
                    mo_spec=mo['mo_spec'],drv=drv,is_mo_output=True)
    #--- Create Output ---
    for i,j in enumerate(mo['mo_spec']):   
      index = numpy.index_exp[:,i] if drv is not None else i
      output.main_output(mo_list[index],geo_info,geo_spec,
                      outputname='%s_%s' % (fid,j['sym']),
                      otype=otype,no_hdf5=True,drv=drv,
                      is_vector=(options.vector is not None))
  return mo_list, mo
  
def mo_set(geo_spec, geo_info, ao_spec, mo_spec, fid_mo_list, 
            drv=None, vector=None, otype=None):
  '''Calculates and saves the density or the derivative thereof 
  using selected molecular orbitals.
  
  **Parameters:**
   
    geo_spec, geo_info, ao_spec, mo_spec :        
      See `Central Variables`_ for details.
    fid_mo_list : str
      Specifies the filename of the molecular orbitals list. 
      If fid_mo_list is 'all_mo', creates a list containing all molecular orbitals.
    otype : str or list of str, optional
      Specifies output file type. See :data:`otypes` for details.
    drv : int or string, {None, 'x', 'y', 'z', 0, 1, 2}, optional
      If not None, a derivative calculation of the atomic orbitals 
      is requested.
    vector : None or int, optional
      If not None, performs the computations on a vectorized grid, i.e., 
      with x, y, and z as vectors.
    
  '''

  mo = mo_select(mo_spec, fid_mo_list)
  if 'h5' in otype:
    try:
      os.remove('%s.h5' % options.outputname)
    except OSError:
      pass
  
  return_data = []
  
  for i_file,j_file in enumerate(mo['mo_in_file']):
    display('\nStarting with the %d. element of the mo-List (%s)...\n\t' % 
                (i_file+1,fid_mo_list) + str(j_file) + 
                '\n\t(Only regarding existing and occupied mos.)\n')
    
    Spec = []
    for i_mo,j_mo in enumerate(mo['mo']):
      if j_mo in j_file: 
        if mo['sym_select']: 
          ii_mo = numpy.argwhere(mo['mo_ii'] == j_mo)
        else: 
          ii_mo = i_mo
        Spec.append(mo['mo_spec'][int(ii_mo)])
    
    data = core.rho_compute(geo_spec, ao_spec, Spec,
                            drv=drv,vector=vector,
                            numproc=options.numproc)
    if drv is None:
      rho = data
    else:
      rho, delta_rho = data
    
    if options.z_reduced_density:
      if vector is not None:
	display(
	'\nSo far, reducing the density is not supported for vectorized grids.\n')
      elif drv is not None:
	display(
	'\nSo far, reducing the density is not supported for the derivative of the density.\n')
      else:
	rho = integrate.simps(rho, grid.x, dx=grid.delta_[0], axis=0, even='avg')
	rho = integrate.simps(rho, grid.y, dx=grid.delta_[1], axis=0, even='avg')
    
    if not options.no_output:
      if 'h5' in otype:
	display('Saving to Hierarchical Data Format file (HDF5)...')
	fid = options.outputname
	group = '/mo_set:%03d' % (i_file+1)	
	display('\n\t%s.h5 in the group "%s"' % (fid,group))	
	output.HDF5_creator(rho,fid,geo_info,geo_spec,data_id='rho',
				   append=group,mo_spec=Spec)
	if options.drv is not None:
	  for i,j in enumerate(options.drv):
	    data_id = 'rho_d%s' % j
	    output.HDF5_creator(delta_rho[i],fid,geo_info,geo_spec,
				       data_id=data_id,data_only=True,
				       append=group,mo_spec=Spec)
      
      fid = '%s_%03d' % (options.outputname, i_file+1) 
      output.main_output(rho,geo_info,geo_spec,outputname=fid,
				otype=otype,no_hdf5=True,
				is_vector=(vector is not None))
      if options.drv is not None:
	for i,j in enumerate(options.drv):
	  fid = '%s_%03d_d%s' % (options.outputname, i_file+1, j) 
	  output.main_output(rho,geo_info,geo_spec,outputname=fid,
				    otype=otype,no_hdf5=True,
				    is_vector=(vector is not None))
	    
  return None
  #--- mo_select ---

def save_mo_hdf5(filename,geo_info,geo_spec,ao_spec,mo_spec,
                 x=None,y=None,z=None,N=None):
  '''Calculate and save selected MOs to an HDF5 File requiering only a
  small amount of the RAM.
  '''
  import h5py
  if x is None: x = grid.x
  if y is None: y = grid.y
  if z is None: z = grid.z
  if N is None: N = grid.N_
  
  #--- Initialize HDF5_File ---
  fid = filename if filename.endswith('.h5') else '%s.h5' % filename
  f = h5py.File(fid, 'w')
  
  f.create_dataset('x',(1,len(x)),data=x)
  f.create_dataset('y',(1,len(y)),data=y)
  f.create_dataset('z',(1,len(z)),data=z)
  
  f.create_dataset('geo_info',(numpy.shape(geo_info)),data=numpy.array(geo_info))
  f.create_dataset('geo_spec',(numpy.shape(geo_spec)),data=geo_spec)
  
  mo_name = []
  occ_num=[]
  energy=[]
  sym=[]
  
  for mo in mo_spec:
    mo_name.append(mo['sym'])
    occ_num.append(mo['occ_num'])
    energy.append(mo['energy'])
    sym.append(mo['sym'])
    
    dID = 'mo:%s' % mo['sym']
    dset = f.create_dataset(dID,N)
  
  
  f.create_dataset('mo:Content',data=numpy.array(mo_name))
  
  mo_info = f.create_group('mo_info')
  mo_info.create_dataset('occ_num',((1,len(mo_spec))),data=occ_num)
  mo_info.create_dataset('energy',((1,len(mo_spec))),data=energy)
  mo_info.create_dataset('sym',((1,len(mo_spec))),data=sym)
  
  f.close()
    
  
  #--- Slice the grid ---
  sDim = 0
  sNum = N[sDim]
  N[sDim] = 1
  xyz = [x,y,z]
  zz = [x,y,z]
  zz[sDim] = numpy.zeros(1)
  for ii_z in range(len(xyz[sDim])):
    zz[sDim][:] = xyz[sDim][ii_z]
    
    ao_list = core.ao_creator(geo_spec,ao_spec,x=zz[0],y=zz[1],z=zz[2],N=N)
    core.mo_creator(ao_list,mo_spec,x=zz[0],y=zz[1],z=zz[2],N=N,
                    vector=False,HDF5_save=fid,h5py=h5py,s=ii_z,
                    numproc=options.numproc)
  
  return 0

def atom2index(atom,geo_info=None):
  '''Converts a list of atom numbers to indices of :data:`geo_info`.'''
  if not (isinstance(atom,list) or isinstance(atom,numpy.ndarray)):
    atom = [atom]
    
  index = []
  
  if geo_info != None:
    geo_info = numpy.array(geo_info)[:,1]
    for a in atom:
      i = numpy.argwhere(geo_info.astype(int) == a)
      if len(i) != 1:
        raise ValueError('No or multiple occurence ' + 
                         'of the atom number %d in geo_info!' % a)
      index.append(int(i))
    
    index = numpy.array(index, dtype=int)
  
  else:
    try:
      index = numpy.array(atom,dtype=int)
    except ValueError:
      raise ValueError('Cannot convert atom to integer array!')
  
  return atom, index

def atom_projected_density(atom,geo_spec,ao_spec,mo_spec,geo_info=None,
                    bReturnmo=False,ao_list=None,mo_list=None,
                    x=None,y=None,z=None,N=None,is_vector=False):
  '''Computes the projected electron density with respect to the selected atoms.
  
  .. math::
  
        \\rho^a(x,y,z) = \sum_i occ_i * mo_i^a * mo_i
  
  
  **Parameters:**
  
    atom : int or list of int
      Specifies the atoms to which the projected electron density will be computed.  
    geo_spec, geo_info, ao_spec, mo_spec :  
      See `Central Variables`_ for details.
    bReturnmo : bool, optional
      If True, the atom projected molecular orbitals are additionally returned.

  **Returns:**
  
    rho_atom : list of numpy.ndarrays, shape=(len(atoms,) + N)
      Contains the atom projected electron density on a grid.
    mo_atom : list of numpy.ndarrays, shape=(len(atoms,NMO) + N)
      Contains the NMO=len(mo_spec) atom projected molecular orbitals on a grid. 
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
  
  atom, index = atom2index(atom,geo_info=geo_info)
  
  display('Computing the atom-projected density with respect to '+ 
          'the atom(s) (internal numbering)')
  outp = '\t['
  for i,a in enumerate(atom):
    if not i % 10 and i != 0:
      outp += '\n\t'
    outp += '%d,\t' % a
  display('%s]\n' % outp[:-2])
  
  display('\tCalculating ao_list & mo_list')
  if ao_list is None:
    ao_list = core.ao_creator(geo_spec,ao_spec,x=x,y=y,z=z,N=N,
                    is_vector=is_vector)
  if mo_list is None:
    mo_list = core.mo_creator(ao_list,mo_spec,x=x,y=y,z=z,N=N,
                              is_vector=is_vector)
    
  display('\tCalculating the atom-projected density')
  
  if bReturnmo: mo_atom = [[] for a in index]
  rho_atom = [numpy.zeros(N) for a in index]
  for i,a in enumerate(index):
    display('\t\tFinished %d of %d' % (i+1,len(index)))
    ao_index = []
    ao = []
    ll = 0
    for ii in ao_spec:
      for l in range(core.l_deg(l=ii['type'])):
        if ii['atom'] == a:
          ao_index.append(ll)
          ao.append(ao_list[ll])
        ll += 1
    
    for ii_mo,spec in enumerate(mo_spec):
      mo = numpy.zeros(N)
      for jj in range(len(ao_index)):
        mo += spec['coeffs'][ao_index[jj]] * ao[jj]
      rho_atom[i] += spec['occ_num'] * mo_list[ii_mo] * mo
      if bReturnmo: mo_atom[i].append(mo)
  
  if bReturnmo:
    display('Returning the atom-projected density and' +
    '\n\tthe atom-projected molecular orbitals')
    return rho_atom, mo_atom
  else:
    display('Returning the atom-projected density')
    return rho_atom
  #--- atom_projected_density ---

def compute_mulliken_charges(atom,geo_info,geo_spec,ao_spec,mo_spec,
            ao_list=None,mo_list=None,rho_atom=None,
            x=None,y=None,z=None,N=None,is_vector=False):
  '''Compute the Mulliken charges of the selected atoms using
  the respective projected electron densities.
  
  .. math::
  
    P^a = \int \\rho^a(x,y,z) d^3r
  
  **Parameters:**
  
    atom : int or list of int
      Compute the numerical Mulliken charges of which atom.
  
  **Returns:**
  
    geo_spec, geo_info, ao_spec, mo_spec :
      See the `Central Variables`_ in the manual for details.
        
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
  
  atom, index = atom2index(atom,geo_info=geo_info)
  
  rho_atom = atom_projected_density(index,geo_spec,ao_spec,mo_spec,
                    ao_list=ao_list,mo_list=mo_list,
                    x=x,y=y,z=z,N=N,is_vector=is_vector)
  
  display('\nNumerical Mulliken Charges:')
  if is_vector:
    display('Warning: You have applied a vectorized grid.\n' + 
    'The integration is not implemented for such a grid!\n' +
    'You will get wrong Mulliken Charges...')
  mulliken_charge = []  
  for i,a in enumerate(index):
    mulliken_charge.append(core.integration(rho_atom[i]))
    #--- Print the Mulliken Charges ---
    a = int(a)
    display('\tAtom %s (%s):\t%+0.4f' % 
            (geo_info[a][1],geo_info[a][0],mulliken_charge[-1]))
  
  return rho_atom, numpy.array(mulliken_charge)


def mo_transition_flux_density(i,j,geo_spec,ao_spec,mo_spec,drv='x',
                    ao_list=None,mo_list=None,
                    delta_ao_list=None,delta_mo_list=None,
                    x=None,y=None,z=None,N=None,is_vector=False):
  '''Calculate one component (e.g. drv='x') of the 
  transition electoronic flux density between the 
  molecular orbitals i and j.
  
  .. math::
  
        moTEFD_{i,j}(r) = <mo_i|\delta(r-r')\\nabla_x|mo_j>
    
  **Parameters:**
  
     i: int
        index of the first mo (Bra)
     j: int
        index of the second mo (Ket)
     drv: {0,1,2,'x','y','z'}
        The desired component of the vector field of the 
        transition electoronic flux density

  **Returns:**
  
     mo_tefd : numpy.ndarray
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
  
  if mo_list is None:
    if ao_list is None:
      display('\tComputing ao_list and ' +
                    'mo #%d, since it is not given.' % i)
      ao_list = core.ao_creator(geo_spec,ao_spec,x=x,y=y,z=z,N=N,
                    is_vector=is_vector)  
    else:
      display('\tComputing mo #%d, since it is not given.' % i)
    mo = core.mo_creator(ao_list,[mo_spec[i]],x=x,y=y,z=z,N=N,
                    is_vector=is_vector)[0]
  else:
    mo = mo_list[i]
  if delta_mo_list is None:
    if delta_ao_list is None:
      display('\tComputing delta_ao_list and the derivative of ' +
                    'mo #%d, since it is not given.' % j)
      delta_ao_list = core.ao_creator(geo_spec,ao_spec,drv=drv,
                        x=x,y=y,z=z,N=N,is_vector=is_vector)
    else:
      display('\tComputing mo #%d, since it is not given.' % j)
    delta_mo = core.mo_creator(delta_ao_list,[mo_spec[j]],x=x,y=y,z=z,N=N,
                    is_vector=is_vector)[0]
  else:
    delta_mo = delta_mo_list[i]
  
  return mo*delta_mo
  #--- mo_transition_flux_density ---
