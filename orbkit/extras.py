# -*- coding: iso-8859-1 -*-
'''Module for all additional features of orbkit.'''
'''
orbkit
Gunter Hermann, Vincent Pohl, and Axel Schild

Institut fuer Chemie und Biochemie, Freie Universitaet Berlin, 14195 Berlin, 
Germany

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

import numpy
from scipy import integrate

# Import orbkit modules
from orbkit import core,grid,output,options
from orbkit.display import display
from orbkit.read import mo_select

def calc_mo(qc, fid_mo_list, drv=None, vector=None, otype=None, ofid=None):
  '''Calculates and saves the selected molecular orbitals or the derivatives thereof.

  **Parameters:**
   
    qc.geo_spec, qc.geo_info, qc.ao_spec, qc.mo_spec :
      See :ref:`Central Variables` for details.
    fid_mo_list : str
      Specifies the filename of the molecular orbitals list or list of molecular
      orbital labels (cf. :mod:`orbkit.read.mo_select` for details). 
      If fid_mo_list is 'all_mo', creates a list containing all molecular orbitals.
    otype : str or list of str, optional
      Specifies output file type. See :data:`otypes` for details.
    ofid : str, optional
      Specifies output file name. If None, the filename will be based on
      :mod:`orbkit.options.outputname`.
    drv : int or string, {None, 'x', 'y', 'z', 0, 1, 2}, optional
      If not None, a derivative calculation of the atomic orbitals 
      is requested.
    vector : None or int, optional
      If not None, performs the computations on a vector grid, i.e., 
      with x, y, and z as vectors.
    
  **Returns:**
    mo_list : numpy.ndarray, shape=((NMO,) + N)
      Contains the NMO=len(qc.mo_spec) molecular orbitals on a grid.
    mo_info : dict 
      Contains information of the selected molecular orbitals and has following Members:
        :mo: - List of molecular orbital labels.
        :mo_ii: - List of molecular orbital indices.
        :mo_spec: - Selected elements of mo_spec. See :ref:`Central Variables` for details.
        :mo_in_file: - List of molecular orbital labels within the fid_mo_list file.
        :sym_select: - If True, symmetry labels have been used. 
      
  '''
  mo_info = mo_select(qc.mo_spec, fid_mo_list, strict=True)
  qc_select = qc.todict()
  qc_select['mo_spec'] = mo_info['mo_spec']
  
  # Calculate the AOs and MOs 
  mo_list = core.rho_compute(qc_select,calc_mo=True,drv=drv,vector=vector,
                             numproc=options.numproc)
  
  if otype is None:
    return mo_list, mo_info
  
  if ofid is None:
    ofid = '%s_MO' % (options.outputname)
  
  if not options.no_output:
    if 'h5' in otype:    
      output.main_output(mo_list,qc.geo_info,qc.geo_spec,data_id='MO',
                    outputname=ofid,
                    mo_spec=qc_select['mo_spec'],drv=drv,is_mo_output=True)
    # Create Output     
    cube_files = []
    for i,j in enumerate(qc_select['mo_spec']):
      outputname = '%s_%s' % (ofid,mo_info['mo'][i])
      comments = ('%s,Occ=%.1f,E=%+.4f' % (mo_info['mo'][i],
                                           j['occ_num'],
                                           j['energy']))
      index = numpy.index_exp[:,i] if drv is not None else i
      output_written = output.main_output(mo_list[index],
                                      qc.geo_info,qc.geo_spec,
                                      outputname=outputname,
                                      comments=comments,
                                      otype=otype,omit=['h5','vmd','mayavi'],
                                      drv=drv,
                                      is_vector=(options.vector is not None))
      
      for i in output_written:
        if i.endswith('.cb'):
          cube_files.append(i)
    
    if 'vmd' in otype:
      display('\nCreating VMD network file...' +
                      '\n\t%(o)s.vmd' % {'o': ofid})
      output.vmd_network_creator(ofid,cube_files=cube_files)
  if 'mayavi' in otype:
    datalabels = ['MO %(sym)s, Occ=%(occ_num).2f, E=%(energy)+.4f E_h' % 
                  i for i in qc_select['mo_spec']]
    if drv is not None:
      tmp = []
      for i in drv:
        for j in datalabels:
          tmp.append('d/d%s of %s' % (i,j))
      datalabels = tmp
    data = mo_list.reshape((-1,) + grid.get_shape())
    
    output.main_output(data,qc.geo_info,qc.geo_spec,
                       otype='mayavi',datalabels=datalabels)
  
  return mo_list, mo_info
  
def mo_set(qc, fid_mo_list, drv=None, laplacian=False, vector=None, 
           otype=None, ofid=None, return_all=True):
  '''Calculates and saves the density or the derivative thereof 
  using selected molecular orbitals.
  
  **Parameters:**
   
    qc.geo_spec, qc.geo_info, qc.ao_spec, qc.mo_spec :
      See :ref:`Central Variables` for details.
    fid_mo_list : str
      Specifies the filename of the molecular orbitals list or list of molecular
      orbital labels (cf. :mod:`orbkit.read.mo_select` for details). 
    otype : str or list of str, optional
      Specifies output file type. See :data:`otypes` for details.
    ofid : str, optional
      Specifies output file name. If None, the filename will be based on
      :mod:`orbkit.options.outputname`.
    drv : int or string, {None, 'x', 'y', 'z', 0, 1, 2}, optional
      If not None, a derivative calculation of the atomic orbitals 
      is requested.
    vector : None or int, optional
      If not None, performs the computations on a vector grid, i.e., 
      with x, y, and z as vectors.
    return_all : bool
      If False, no data will be returned.
  
  **Returns:**
    datasets : numpy.ndarray, shape=((NSET,) + N)
      Contains the NSET molecular orbital sets on a grid.
    delta_datasets : numpy.ndarray, shape=((NSET,NDRV) + N)
      Contains the NDRV NSET molecular orbital set on a grid. This is only 
      present if derivatives are requested.
    mo_info : dict 
      Contains information of the selected molecular orbitals and has following Members:
        :mo: - List of molecular orbital labels.
        :mo_ii: - List of molecular orbital indices.
        :mo_spec: - Selected elements of mo_spec. See :ref:`Central Variables` for details.
        :mo_in_file: - List of molecular orbital labels within the fid_mo_list file.
        :sym_select: - If True, symmetry labels have been used. 
  '''

  mo_info = mo_select(qc.mo_spec, fid_mo_list, strict=False)
  qc_select = qc.todict()
  
  if ofid is None:
    ofid = options.outputname
  if 'h5' in otype and os.path.exists('%s.h5' % ofid):
    raise IOError('%s.h5 already exists!' % ofid)
  
  datasets = []
  delta_datasets = []
  cube_files = []
  for i_file,j_file in enumerate(mo_info['mo_in_file']):
    display('Starting with the %d. element of the molecular orbital list (%s)...\n\t' % 
                (i_file+1,fid_mo_list) + str(j_file) + 
                '\n\t(Only regarding existing and occupied mos.)\n')
    
    qc_select['mo_spec'] = []
    for i_mo,j_mo in enumerate(mo_info['mo']):
      if j_mo in j_file:
        if mo_info['sym_select']: 
          ii_mo = numpy.argwhere(mo_info['mo_ii'] == j_mo)
        else: 
          ii_mo = i_mo
        qc_select['mo_spec'].append(mo_info['mo_spec'][int(ii_mo)])
    
    data = core.rho_compute(qc_select,
                            drv=drv,
                            laplacian=laplacian,
                            vector=vector,
                            numproc=options.numproc)
    datasets.append(data)
    if drv is None:
      rho = data
    elif options.laplacian:
      rho, delta_rho, laplacian_rho = data 
      delta_datasets.extend(delta_rho)
      delta_datasets.append(laplacian_rho)
    else:
      rho, delta_rho = data
      delta_datasets.append(delta_rho)
    
    if options.z_reduced_density:
      if vector is not None:
        display('\nSo far, reducing the density is not supported for vector grids.\n')
      elif drv is not None:
        display('\nSo far, reducing the density is not supported for the derivative of the density.\n')
      else:
        rho = integrate.simps(rho, grid.x, dx=grid.delta_[0], axis=0, even='avg')
        rho = integrate.simps(rho, grid.y, dx=grid.delta_[1], axis=0, even='avg')
    
    if not options.no_output:
      if 'h5' in otype:
        display('Saving to Hierarchical Data Format file (HDF5)...')
        group = '/mo_set:%03d' % (i_file+1)	
        display('\n\t%s.h5 in the group "%s"' % (ofid,group))	
        output.HDF5_creator(rho,ofid,qc.geo_info,qc.geo_spec,data_id='rho',
                            mode='w',group=group,mo_spec=qc_select['mo_spec'])
        if options.drv is not None:
          for i,j in enumerate(options.drv):
            data_id = 'rho_d%s' % j
            output.HDF5_creator(delta_rho[i],ofid,qc.geo_info,qc.geo_spec,
                                data_id=data_id,data_only=True,mode='a',
                                group=group,mo_spec=qc_select['mo_spec'])
          if options.laplacian:
            data_id = 'rho_laplacian' 
            output.HDF5_creator(laplacian_rho,ofid,qc.geo_info,qc.geo_spec,
                                data_id=data_id,data_only=True,mode='a',
                                group=group,mo_spec=qc_select['mo_spec'])
            
      fid = '%s_%03d' % (ofid, i_file+1) 
      cube_files.append('%s.cb' % fid)
      comments = ('mo_set:'+','.join(j_file))
      output.main_output(rho,qc.geo_info,qc.geo_spec,outputname=fid,
                         otype=otype,omit=['h5','vmd','mayavi'],
                         comments=comments,is_vector=(vector is not None))
      if options.drv is not None:
        for i,j in enumerate(options.drv):
          fid = '%s_%03d_d%s' % (ofid, i_file+1, j) 
          cube_files.append('%s.cb' % fid)
          comments = ('d%s_of_mo_set:' % j + ','.join(j_file))
          output.main_output(delta_rho[i],qc.geo_info,qc.geo_spec,outputname=fid,
                             otype=otype,omit=['h5','vmd','mayavi'],
                             comments=comments,is_vector=(vector is not None))  
        if options.laplacian:
          fid = '%s_%03d_laplacian' % (ofid, i_file+1) 
          cube_files.append('%s.cb' % fid)
          comments = ('laplacian_of_mo_set:' + ','.join(j_file))
          output.main_output(laplacian_rho,qc.geo_info,qc.geo_spec,outputname=fid,
                             otype=otype,omit=['h5','vmd','mayavi'],
                             comments=comments,is_vector=(vector is not None))  
          
  if 'vmd' in otype and options.vector is None:
      display('\nCreating VMD network file...' +
                      '\n\t%(o)s.vmd' % {'o': ofid})
      output.vmd_network_creator(ofid,cube_files=cube_files)
   
  datasets = numpy.array(datasets)
  if drv is None:
    if 'mayavi' in otype:
      output.main_output(datasets,qc.geo_info,qc.geo_spec,
                       otype='mayavi',datalabels=mo_info['mo_in_file'])
    return datasets, mo_info
  else:
    delta_datasets = numpy.array(delta_datasets)
    if 'mayavi' in otype:
      datalabels = []
      for i in mo_info['mo_in_file']:
        datalabels.extend(['d/d%s of %s' % (j,i) for j in drv])
        if options.laplacian:
          datalabels.append('laplacian of %s' % i)
      #output.view_with_mayavi(grid.x,grid.y,grid.z,
                               #delta_datasets.reshape((-1,) + tuple(grid.N_)),
                               #geo_spec=qc.geo_spec,
      output.main_output(delta_datasets.reshape((-1,) + grid.get_shape()),
                         qc.geo_info,qc.geo_spec,otype='mayavi',datalabels=datalabels)
    return datasets, delta_datasets, mo_info
  # mo_set 

def calc_ao(qc, drv=None, is_vector=False, otype=None, ofid=None):  
  '''Computes and saves all atomic orbital or a derivative thereof.
  
  **Parameters:**
   
    qc.geo_spec, qc.geo_info, qc.ao_spec :
      See :ref:`Central Variables` for details.
    otype : str or list of str, optional
      Specifies output file type. See :data:`otypes` for details.
    ofid : str, optional
      Specifies output file name. If None, the filename will be based on
      :mod:`orbkit.options.outputname`.
    drv : int or string, {None, 'x', 'y', 'z', 0, 1, 2}, optional
      If not None, a derivative calculation of the atomic orbitals 
      is requested.
    is_vector : bool, optional
      If True, performs the computations on a vector grid, i.e., 
      with x, y, and z as vectors.
    
  **Returns:**    
  
    ao_list : numpy.ndarray, shape=((NAO,) + N)
      Contains the computed NAO atomic orbitals on a grid. Is only returned, if
      otype is None.
  '''
  from omp_functions import run
  global global_val
  
  if ofid is None:
    ofid = options.outputname
  dstr = '' if drv is None else '_d%s' % drv
  
  ao_spec = []
  lxlylz = []
  datalabels = []
  for sel_ao in range(len(qc.ao_spec)):
    if 'exp_list' in qc.ao_spec[sel_ao].keys():
      l = qc.ao_spec[sel_ao]['exp_list']
    else:
      l = core.exp[core.lquant[qc.ao_spec[sel_ao]['type']]]
    lxlylz.extend(l)
    for i in l:
      ao_spec.append(qc.ao_spec[sel_ao].copy())
      ao_spec[-1]['exp_list'] = [i]
      datalabels.append('lxlylz=%s,atom=%d' %(i,ao_spec[-1]['atom']))

  global_val = {'qc':qc,'drv':drv,'is_vector':is_vector,'otype':otype,
                'ofid':ofid,'ao_spec':ao_spec}
  
  ao = run(get_ao,x=numpy.arange(len(ao_spec)).reshape((-1,1)),numproc=options.numproc,display=display)
  
  del global_val
  
  if otype is None or 'h5' in otype or 'mayavi' in otype:   
    ao_list = []
    for i in ao:
      ao_list.extend(i)      
    ao_list = numpy.array(ao_list)
    if otype is None:
      return ao_list
    
    if 'h5' in otype:
      import h5py
      fid = '%s_AO%s.h5' % (ofid,dstr)
      display('Saving to Hierarchical Data Format file (HDF5)...\n\t%s' % fid)
      output.hdf5_write(fid,mode='w',gname='general_info',
                        x=grid.x,y=grid.y,z=grid.z,
                        geo_info=qc.geo_info,geo_spec=qc.geo_spec,
                        lxlylz=numpy.array(lxlylz,dtype=numpy.int64),
                        aolabels=numpy.array(datalabels),
                        grid_info=numpy.array(is_vector,dtype=int))      
      for f in output.hdf5_open(fid,mode='a'):
        output.hdf5_append(ao_spec,f,name='ao_spec')
        output.hdf5_append(ao_list,f,name='ao_list')
      
    if 'mayavi' in otype and options.vector is None:
      output.view_with_mayavi(grid.x,grid.y,grid.z,
                               ao_list,
                               geo_spec=qc.geo_spec,
                               datalabels=datalabels)
  
  if 'vmd' in otype and options.vector is None:
    fid = '%s_AO%s' % (ofid,dstr)
    display('\nCreating VMD network file...' +
                    '\n\t%(o)s.vmd' % {'o': fid})
    output.vmd_network_creator(fid,
                            cube_files=['%s_AO%s_%03d.cb' % (ofid,
                                        dstr,x) for x in range(len(ao_spec))])
  return None
  # calc_ao

def save_mo_hdf5(filename,geo_info,geo_spec,ao_spec,mo_spec,is_vector=False,
                 x=None,y=None,z=None):
  '''Calculate and save selected MOs to an HDF5 File requiering only a
  small amount of the RAM.
  '''
  import h5py
  if x is None: x = grid.x
  if y is None: y = grid.y
  if z is None: z = grid.z  
  if not is_vector:
    N = (len(x), len(y), len(z))
  else:
    if len(x) != len(y) or len(x) != len(z):
      raise ValueError("Dimensions of x-, y-, and z- coordinate differ!")      
    else:
      N = (len(x),)
  
  # Initialize HDF5_File 
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
  
  for mo_info in mo_spec:
    mo_name.append(mo_info['sym'])
    occ_num.append(mo_info['occ_num'])
    energy.append(mo_info['energy'])
    sym.append(mo_info['sym'])
    
    dID = 'mo_info:%s' % mo_info['sym']
    dset = f.create_dataset(dID,N)
  
  
  f.create_dataset('mo_info:Content',data=numpy.array(mo_name))
  
  mo_info = f.create_group('mo_info')
  mo_info.create_dataset('occ_num',((1,len(mo_spec))),data=occ_num)
  mo_info.create_dataset('energy',((1,len(mo_spec))),data=energy)
  mo_info.create_dataset('sym',((1,len(mo_spec))),data=sym)
  
  f.close()
    
  
  # Slice the grid 
  sDim = 0
  N[sDim] = 1
  xyz = [x,y,z]
  zz = [x,y,z]
  zz[sDim] = numpy.zeros(1)
  for ii_z in range(len(xyz[sDim])):
    zz[sDim][:] = xyz[sDim][ii_z]
    
    ao_list = core.ao_creator(geo_spec,ao_spec,x=zz[0],y=zz[1],z=zz[2])
    core.mo_creator(ao_list,mo_spec,x=zz[0],y=zz[1],z=zz[2],
                    vector=False,HDF5_save=fid,h5py=h5py,s=ii_z,
                    numproc=options.numproc)
  
  return 0

def get_ao(x):
  ao_list = core.l_creator(global_val['qc'].geo_spec,global_val['ao_spec'],
                           int(x),exp_list=global_val['ao_spec'][x]['exp_list'],
                           is_vector=global_val['is_vector'],
                           drv=global_val['drv'])
  
  if global_val['otype'] is not None:
    drv = global_val['drv']
    comments = '%03d.lxlylz=%s,at=%d' %(x,
                                         global_val['ao_spec'][x]['exp_list'][0],
                                         global_val['ao_spec'][x]['atom'])
    output.main_output(ao_list[0] if drv is None else ao_list,
                       global_val['qc'].geo_info,
                       global_val['qc'].geo_spec,
                       comments=comments,
                       outputname='%s_AO_%03d' % (global_val['ofid'],x),
                       otype=global_val['otype'],
                       omit=['h5','vmd','mayavi'],
                       drv = None if drv is None else [drv],
                       is_vector=global_val['is_vector'])
  return ao_list

def atom2index(atom,geo_info=None):
  '''Converts a list of atom numbers to indices of :data:`geo_info`.
  
  **Parameters:**
  
    atom : int or list of int
      Specifies the indices of the atoms (counting from one).  
    geo_info : optional
      See :ref:`Central Variables` for details.
  
  **Returns:**
  
    atom : list of int
      Contains the indices of the atoms (counting from one). 
    index : numpy.ndarray, dtype=int
      Contains the indices of the atoms as occuring in qc.geo_info 
      (counting from zero). 
    
  '''
   
  if not (isinstance(atom,list) or isinstance(atom,numpy.ndarray)):
    atom = [atom]
    
  index = []
  
  if geo_info is not None:
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

def gross_atomic_density(atom,qc,
                    bReturnmo=False,ao_list=None,mo_list=None,
                    x=None,y=None,z=None,drv=None,is_vector=False):
  r'''Computes the gross atomic density with respect to the selected atoms.
  
  .. math::
  
    \rho^a = \sum_i^{N_{\rm MO}} {\rm occ}_i \cdot \varphi_i^a \cdot \varphi_i
  
  
  **Parameters:**
  
    atom : 'all' or int or list of int
      Specifies the atoms (counting from one) for which the gross atomic density
      will be computed.  
      If (atom == 'all') or (atom == -1), computes the gross atomic density  
      for all atoms.
    qc.geo_spec, qc.geo_info, qc.ao_spec, qc.mo_spec :
      See :ref:`Central Variables` for details.
    bReturnmo : bool, optional
      If True, the gross atomic molecular orbitals are additionally returned.

  **Returns:**
  
    rho_atom : list of numpy.ndarrays, shape=(len(atoms,) + N)
      Contains the atom gross atomic density on a grid.
    mo_atom : list of numpy.ndarrays, shape=(len(atoms,NMO) + N)
      Contains the NMO=len(mo_spec) gross atomic molecular orbitals on a grid. 
  '''  
  
  if x is None: x = grid.x
  if y is None: y = grid.y
  if z is None: z = grid.z
  
  if (atom == 'all' or atom == -1):
    atom = range(1, len(qc.geo_info) + 1)  
  
  if not is_vector:
    N = (len(x), len(y), len(z))
  else:
    if len(x) != len(y) or len(x) != len(z):
      raise ValueError("Dimensions of x-, y-, and z- coordinate differ!")      
    else:
      N = (len(x),)
    
  atom, index = atom2index(atom,geo_info=qc.geo_info)
  
  display('Computing the gross atomic density with respect to '+ 
          'the atom(s) (internal numbering)')
  outp = '\t['
  for i,a in enumerate(atom):
    if not i % 10 and i != 0:
      outp += '\n\t'
    outp += '%d,\t' % a
  display('%s]\n' % outp[:-2])
  
  display('\tCalculating ao_list & mo_list')
  if ao_list is None:
    ao_list = core.ao_creator(qc.geo_spec,qc.ao_spec,x=x,y=y,z=z,drv=drv,
                    is_vector=is_vector)
  if mo_list is None:
    mo_list = core.mo_creator(ao_list,qc.mo_spec,x=x,y=y,z=z,
                              is_vector=is_vector)
    
  display('\tCalculating the gross atomic density')
  
  if bReturnmo: mo_atom = [[] for a in index]
  rho_atom = [numpy.zeros(N) for a in index]
  for i,a in enumerate(index):
    display('\t\tFinished %d of %d' % (i+1,len(index)))
    ao_index = []
    ao = []
    ll = 0
    for ii in qc.ao_spec:
      for l in range(core.l_deg(l=ii['type'])):
        if ii['atom'] == a:
          ao_index.append(ll)
          ao.append(ao_list[ll])
        ll += 1
    
    for ii_mo,spec in enumerate(qc.mo_spec):
      mo_info = numpy.zeros(N)
      for jj in range(len(ao_index)):
        mo_info += spec['coeffs'][ao_index[jj]] * ao[jj]
      rho_atom[i] += spec['occ_num'] * mo_list[ii_mo] * mo_info
      if bReturnmo: mo_atom[i].append(mo_info)
  
  string = 'Returning the gross atomic density'
  if bReturnmo:
    display(string + ' and\n\tthe gross atomic molecular orbitals')
    return rho_atom, mo_atom
  else:
    display(string)
    return rho_atom
  # gross_atomic_density 

def numerical_mulliken_charges(atom,qc,
            ao_list=None,mo_list=None,rho_atom=None,
            x=None,y=None,z=None,is_vector=False):
  r'''Compute the Mulliken charges and gross atomic populations of the selected 
  atoms *numerically* using the respective gross atomic densities.
  
  .. math::
  
    P^a = \int \rho^a(x,y,z) {\rm d}^3r
  
  **Parameters:**
  
    atom : 'all' or int or list of int
      Specifies the atom (numbering starting from one) to which the numerical 
      Mulliken charges and gross populations will be computed. 
      If (atom == 'all') or (atom == -1), computes the charges and populations 
      for all atoms.
  
  **Returns:**
  
    rho_atom : list of numpy.ndarrays, shape=(len(atoms,) + N)
      Contains the atom gross atomic density on a grid.
    mulliken_num : dict, (present if not is_vector)
      Contains information of Mulliken charge analysis and has following members:
        :population: - Mulliken population for each atom.
        :charges: - Mulliken charges for each atom.
  '''
  
  if x is None: x = grid.x
  if y is None: y = grid.y
  if z is None: z = grid.z
  
  if (atom == 'all' or atom == -1 or atom == [-1]):
    atom = range(1, len(qc.geo_info) + 1) 
  
  atom, index = atom2index(atom,geo_info=qc.geo_info)
  
  rho_atom = gross_atomic_density(atom,qc,
                    ao_list=ao_list,mo_list=mo_list,
                    x=x,y=y,z=z,is_vector=is_vector)
  
  if is_vector:
    display('Warning: You have applied a vector grid.\n' + 
    'The integration is not implemented for such a grid!\n' +
    '\nReturning only the gross atomic densities...')
    return rho_atom
  
  GP_A = numpy.array([core.integration(i) for i in rho_atom])
  mulliken_num = {'population': GP_A,
                  'charge': numpy.array(qc.geo_info[:,2],dtype=float)-GP_A}  
  
  display('\nNumerical Mulliken Charges Q and Gross Atomic Populations GP:')
  for i,a in enumerate(index):
    a = int(a)
    display('\tAtom %s (%s):\tQ = %+0.4f ( GP = %0.4f )' % 
            (qc.geo_info[a][1],qc.geo_info[a][0],
            mulliken_num['charge'][i],mulliken_num['population'][i]))     
  return rho_atom, mulliken_num


def mo_transition_flux_density(i,j,qc,drv='x',
                    ao_list=None,mo_list=None,
                    delta_ao_list=None,delta_mo_list=None,
                    x=None,y=None,z=None,is_vector=False):
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
      ao_list = core.ao_creator(qc.geo_spec,qc.ao_spec,x=x,y=y,z=z,
                    is_vector=is_vector)  
    else:
      display('\tComputing mo #%d, since it is not given.' % i)
    mo = core.mo_creator(ao_list,[qc.mo_spec[i]],x=x,y=y,z=z,
                    is_vector=is_vector)[0]
  else:
    mo = mo_list[i]
  if delta_mo_list is None:
    if delta_ao_list is None:
      display('\tComputing delta_ao_list and the derivative of ' +
                    'mo #%d, since it is not given.' % j)
      delta_ao_list = core.ao_creator(qc.geo_spec,qc.ao_spec,drv=drv,
                        x=x,y=y,z=z,is_vector=is_vector)
    else:
      display('\tComputing mo #%d, since it is not given.' % j)
    delta_mo = core.mo_creator(delta_ao_list,[qc.mo_spec[j]],x=x,y=y,z=z,
                    is_vector=is_vector)[0]
  else:
    delta_mo = delta_mo_list[i]
  
  return mo*delta_mo
  # mo_transition_flux_density 
