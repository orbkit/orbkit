# -*- coding: iso-8859-1 -*-
'''Module for all additional features of ORBKIT.'''
'''
ORBKIT
Gunter Hermann, Vincent Pohl, Lukas Eugen Marsoner Steinkasserer, Axel Schild, and Jean Christophe Tremblay

Institut fuer Chemie und Biochemie, Freie Universitaet Berlin, 14195 Berlin, 
Germany

This file is part of ORBKIT.

ORBKIT is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or any later version.

ORBKIT is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with ORBKIT.  If not, see <http://www.gnu.org/licenses/>.
'''

# Import general modules
import os

import numpy
from scipy import integrate

# Import orbkit modules
from orbkit import core,grid, options
from orbkit.write import vmd_network_creator,main_output,hdf5_creator,hdf5_write,hdf5_append
from orbkit.display import display
from orbkit.orbitals import MOClass


def calc_mo(qc, fid_mo_list, drv=None, otype=None, ofid=None,
            numproc=None, slice_length=None):
  '''Calculates and saves the selected molecular orbitals or the derivatives thereof.
  **Parameters:**
   
    qc.geo_spec, qc.geo_info, qc.ao_spec, qc.mo_spec :
      See :ref:`Central Variables` for details.
    fid_mo_list : str
      Specifies the filename of the molecular orbitals list or list of molecular
      orbital labels (cf. :mod:`orbkit.read.mo_select` for details). 
      If fid_mo_list is 'all_mo', creates a list containing all molecular orbitals.
    drv : int or string, {None, 'x', 'y', 'z', 0, 1, 2}, optional
      If not None, a derivative calculation of the atomic orbitals 
      is requested.
    otype : str or list of str, optional
      Specifies output file type. See :data:`otypes` for details.
    ofid : str, optional
      Specifies output file name. If None, the filename will be based on
      :mod:`orbkit.options.outputname`.
    numproc : int, optional
      Specifies number of subprocesses for multiprocessing. 
      If None, uses the value from :mod:`options.numproc`.
    slice_length : int, optional
      Specifies the number of points per subprocess.
      If None, uses the value from :mod:`options.slice_length`.
    
  **Returns:**
    mo_list : numpy.ndarray, shape=((NMO,) + N)
      Contains the NMO=len(qc.mo_spec) molecular orbitals on a grid.
  '''
  mo_spec = qc.mo_spec.select(fid_mo_list, flatten_input=True)
  qc_select = qc.copy()
  qc_select.mo_spec = mo_spec

  slice_length = options.slice_length if slice_length is None else slice_length
  numproc = options.numproc if numproc is None else numproc

  # Calculate the AOs and MOs 
  mo_list = core.rho_compute(qc_select,
                             calc_mo=True,
                             drv=drv,
                             slice_length=slice_length,
                             numproc=numproc)
  
  if otype is None:
    return mo_list
  
  if ofid is None:
    ofid = '%s_MO' % (options.outputname)
  
  if not options.no_output:
    if 'h5' in otype:    
      main_output(mo_list,qc.geo_info,qc.geo_spec,data_id='MO',
                  outputname=ofid,
                  mo_spec=qc_select.mo_spec,drv=drv,is_mo_output=True)
    # Create Output     
    cube_files = []
    for i in range(len(qc_select.mo_spec)):
      outputname = '%s_%s' % (ofid,qc_select.mo_spec.selected_mo[i])
      comments = ('%s,Occ=%.1f,E=%+.4f' % (qc_select.mo_spec.selected_mo[i],
                                           qc_select.mo_spec.get_occ()[i],
                                           qc_select.mo_spec.get_eig()[i]))
      index = numpy.index_exp[:,i] if drv is not None else i
      output_written = main_output(mo_list[index],
                                   qc.geo_info,qc.geo_spec,
                                   outputname=outputname,
                                   comments=comments,
                                   otype=otype,omit=['h5','vmd','mayavi'],
                                   drv=drv)
      
      for i in output_written:
        if i.endswith('.cb'):
          cube_files.append(i)
    
    if 'vmd' in otype and cube_files != []:
      display('\nCreating VMD network file...' +
                      '\n\t%(o)s.vmd' % {'o': ofid})
      vmd_network_creator(ofid,cube_files=cube_files)
  if 'mayavi' in otype:
    datalabels = ['MO %(sym)s, Occ=%(occ_num).2f, E=%(energy)+.4f E_h' % 
                  i for i in qc_select.mo_spec]
    if drv is not None:
      tmp = []
      for i in drv:
        for j in datalabels:
          tmp.append('d/d%s of %s' % (i,j))
      datalabels = tmp
    data = mo_list.reshape((-1,) + grid.get_shape())
    
    main_output(data,qc.geo_info,qc.geo_spec,
                       otype='mayavi',datalabels=datalabels)
  
  return mo_list
  
def mo_set(qc, fid_mo_list, drv=None, laplacian=None,
           otype=None, ofid=None, return_all=True,
           numproc=None, slice_length=None):
  '''Calculates and saves the density or the derivative thereof 
  using selected molecular orbitals.
  
  **Parameters:**
   
    qc.geo_spec, qc.geo_info, qc.ao_spec, qc.mo_spec :
      See :ref:`Central Variables` for details.
    fid_mo_list : str
      Specifies the filename of the molecular orbitals list or list of molecular
      orbital labels (cf. :mod:`orbkit.orbitals.MOClass.select` for details). 
    otype : str or list of str, optional
      Specifies output file type. See :data:`otypes` for details.
    ofid : str, optional
      Specifies output file name. If None, the filename will be based on
      :mod:`orbkit.options.outputname`.
    drv : int or string, {None, 'x', 'y', 'z', 0, 1, 2}, optional
      If not None, a derivative calculation of the atomic orbitals 
      is requested.
    return_all : bool
      If False, no data will be returned.
    numproc : int, optional
      Specifies number of subprocesses for multiprocessing. 
      If None, uses the value from :mod:`options.numproc`.
    slice_length : int, optional
      Specifies the number of points per subprocess.
      If None, uses the value from :mod:`options.slice_length`.
  
  **Returns:**
    datasets : numpy.ndarray, shape=((NSET,) + N)
      Contains the NSET molecular orbital sets on a grid.
    delta_datasets : numpy.ndarray, shape=((NSET,NDRV) + N)
      Contains the NDRV NSET molecular orbital set on a grid. This is only 
      present if derivatives are requested.
  '''

  #Can be an mo_spec or a list of mo_spec
  # For later itteration we'll make it into a list here if it is not
  mo_info_list = qc.mo_spec.select(fid_mo_list)
  if isinstance(mo_info_list, MOClass):
    mo_info_list = [mo_info_list]
    
  drv = options.drv if drv is None else drv
  laplacian = options.laplacian if laplacian is None else laplacian
  slice_length = options.slice_length if slice_length is None else slice_length
  numproc = options.numproc if numproc is None else numproc
  
  if ofid is None:
    ofid = options.outputname
  if 'h5' in otype and os.path.exists('%s.h5' % ofid):
    raise IOError('%s.h5 already exists!' % ofid)
  
  datasets = []
  delta_datasets = []
  cube_files = []

  for i_file, mo_info in enumerate(mo_info_list):
    qc_select = qc.copy()
    qc_select.mo_spec = mo_info
    display('Starting with the molecular orbital list'
            + str(mo_info.selection_string) +
            '\n\t(Only regarding existing and occupied mos.)\n')
      
    data = core.rho_compute(qc_select,
                            drv=drv,
                            laplacian=laplacian,
                            slice_length=slice_length,
                            numproc=numproc)

    datasets.append(data)

    if drv is None:
      rho = data
    elif laplacian:
      rho, delta_rho, laplacian_rho = data 
      delta_datasets.extend(delta_rho)
      delta_datasets.append(laplacian_rho)
    else:
      rho, delta_rho = data
      delta_datasets.append(delta_rho)
    
    if options.z_reduced_density:
      if grid.is_vector:
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
        HDF5_creator(rho,ofid,qc.geo_info,qc.geo_spec,data_id='rho',
                            mode='w',group=group,mo_spec=qc_select['mo_spec'])
        if drv is not None:
          for i,j in enumerate(drv):
            data_id = 'rho_d%s' % j
            HDF5_creator(delta_rho[i],ofid,qc.geo_info,qc.geo_spec,
                                data_id=data_id,data_only=True,mode='a',
                                group=group,mo_spec=qc_select['mo_spec'])
          if laplacian:
            data_id = 'rho_laplacian' 
            HDF5_creator(laplacian_rho,ofid,qc.geo_info,qc.geo_spec,
                                data_id=data_id,data_only=True,mode='a',
                                group=group,mo_spec=qc_select['mo_spec'])
            
      fid = '%s_%03d' % (ofid, i_file+1) 
      cube_files.append('%s.cb' % fid)
      comments = ('mo_set:'+','.join(mo_info.selection_string))

      main_output(rho,qc.geo_info,qc.geo_spec,outputname=fid,
                         otype=otype,omit=['h5','vmd','mayavi'],
                         comments=comments)
      if drv is not None:
        for i,j in enumerate(drv):
          fid = '%s_%03d_d%s' % (ofid, i_file+1, j) 
          cube_files.append('%s.cb' % fid)
          comments = ('d%s_of_mo_set:' % j + ','.join(mo_info.selection_string))
          main_output(delta_rho[i],qc.geo_info,qc.geo_spec,outputname=fid,
                             otype=otype,omit=['h5','vmd','mayavi'],
                             comments=comments)  
        if laplacian:
          fid = '%s_%03d_laplacian' % (ofid, i_file+1) 
          cube_files.append('%s.cb' % fid)
          comments = ('laplacian_of_mo_set:' + ','.join(mo_info.selection_string))
          main_output(laplacian_rho,qc.geo_info,qc.geo_spec,outputname=fid,
                             otype=otype,omit=['h5','vmd','mayavi'],
                             comments=comments)  
          
  if 'vmd' in otype and cube_files != []:
      display('\nCreating VMD network file...' +
                      '\n\t%(o)s.vmd' % {'o': ofid})
      vmd_network_creator(ofid,cube_files=cube_files)
   
  datasets = numpy.array(datasets)
  if drv is None:
    if 'mayavi' in otype:
      main_output(datasets,qc.geo_info,qc.geo_spec,
                       otype='mayavi',datalabels=mo_info.selected_mo)
    print(datasets.shape)
    return datasets#, mo_info
  else:
    delta_datasets = numpy.array(delta_datasets)
    if 'mayavi' in otype:
      datalabels = []
      for i in mo_info.selected_mo:
        datalabels.extend(['d/d%s of %s' % (j,i) for j in drv])
        if laplacian:
          datalabels.append('laplacian of %s' % i)
      main_output(delta_datasets.reshape((-1,) + grid.get_shape()),
                         qc.geo_info,qc.geo_spec,otype='mayavi',datalabels=datalabels)
    return datasets, delta_datasets
  # mo_set 

def calc_ao(qc, drv=None, otype=None, ofid=None):  
  '''Computes and saves all atomic orbital or a derivative thereof.
  
  **Parameters:**
   
    qc.geo_spec, qc.geo_info, qc.ao_spec :
      See :ref:`Central Variables` for details.
    otype : str or list of str, optional
      Specifies output file type. See :data:`otypes` for details.
    ofid : str, optional
      Specifies output file name. If None, the filename will be based on
      :mod:`orbkit.options.outputname`.
    drv : int or string, {None, 'x', 'y', 'z', 'xx', 'xy', ...}, optional
      If not None, a derivative calculation of the atomic orbitals 
      is requested.
    
  **Returns:**    
  
    ao_list : numpy.ndarray, shape=((NAO,) + N)
      Contains the computed NAO atomic orbitals on a grid. Is only returned, if
      otype is None.
  '''
  from omp_functions import run
  
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
  
  global_args = {'geo_spec':qc.geo_spec,
                 'geo_info':qc.geo_info,
                 'ao_spec':ao_spec,
                 'drv':drv,
                 'x':grid.x,
                 'y':grid.y,
                 'z':grid.z,
                 'is_vector':grid.is_vector,
                 'otype':otype,
                 'ofid':ofid}
  display('Starting the computation of the %d atomic orbitals'%len(ao_spec)+
         (' using %d subprocesses.' % options.numproc if options.numproc > 1 
          else '.' )
          )
  ao = run(get_ao,x=numpy.arange(len(ao_spec)).reshape((-1,1)),
           numproc=options.numproc,display=display,
           initializer=initializer,global_args=global_args)
    
  if otype is None or 'h5' in otype or 'mayavi' in otype:   
    ao_list = []
    for i in ao:
      ao_list.extend(i)      
    ao_list = numpy.array(ao_list)
    if otype is None or options.no_output:
      return ao_list
    
    if 'h5' in otype:
      import h5py
      fid = '%s_AO%s.h5' % (ofid,dstr)
      display('Saving to Hierarchical Data Format file (HDF5)...\n\t%s' % fid)
      hdf5_write(fid,mode='w',gname='general_info',
                        x=grid.x,y=grid.y,z=grid.z,
                        geo_info=qc.geo_info,geo_spec=qc.geo_spec,
                        lxlylz=numpy.array(lxlylz,dtype=numpy.int64),
                        aolabels=numpy.array(datalabels),
                        grid_info=numpy.array(grid.is_vector,dtype=int))      
      for f in hdf5_open(fid,mode='a'):
        hdf5_append(ao_spec,f,name='ao_spec')
        hdf5_append(ao_list,f,name='ao_list')
      
    if 'mayavi' in otype:
      main_output(ao_list,qc.geo_info,qc.geo_spec,
                  otype='mayavi',datalabels=datalabels)

  if 'vmd' in otype and options.vector is None:
    fid = '%s_AO%s' % (ofid,dstr)
    display('\nCreating VMD network file...' +
                    '\n\t%(o)s.vmd' % {'o': fid})
    vmd_network_creator(fid,
                        cube_files=['%s_AO%s_%03d.cb' % (ofid,
                        dstr,x) for x in range(len(ao_spec))])
  # calc_ao

def initializer(gargs):
  global global_args
  global_args = gargs
  grid.set_grid(global_args['x'],
                global_args['y'],
                global_args['z'],
                is_vector=global_args['is_vector'])

def get_ao(x):
  ao_list = core.ao_creator(global_args['geo_spec'],
                            [global_args['ao_spec'][int(x)]],
                            drv=global_args['drv'])
  
  if global_args['otype'] is not None:
    drv = global_args['drv']
    comments = '%03d.lxlylz=%s,at=%d' %(x,
                                         global_args['ao_spec'][x]['exp_list'][0],
                                         global_args['ao_spec'][x]['atom'])
    
    main_output(ao_list[0] if drv is None else ao_list,
                global_args['geo_info'],
                global_args['geo_spec'],
                comments=comments,
                outputname='%s_AO_%03d' % (global_args['ofid'],x),
                otype=global_args['otype'],
                omit=['h5','vmd','mayavi'],
                drv = None if drv is None else [drv])
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
                    bReturnmo=False,ao_list=None,mo_list=None,drv=None):
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
    
  if (atom == 'all' or atom == -1):
    atom = range(1, len(qc.geo_info) + 1)  
  
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
    ao_list = core.ao_creator(qc.geo_spec,qc.ao_spec,drv=drv)
  if mo_list is None:
    mo_list = core.mo_creator(ao_list,qc.mo_spec)
  
  display('\tCalculating the gross atomic density')
  N = mo_list.shape[1:]
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
            ao_list=None,mo_list=None,rho_atom=None):
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
  
  if (atom == 'all' or atom == -1 or atom == [-1]):
    atom = range(1, len(qc.geo_info) + 1) 
  
  atom, index = atom2index(atom,geo_info=qc.geo_info)
  
  rho_atom = gross_atomic_density(atom,qc,
                    ao_list=ao_list,mo_list=mo_list)
  
  if grid.is_vector:
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
                    delta_ao_list=None,delta_mo_list=None):
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
  if mo_list is None:
    if ao_list is None:
      display('\tComputing ao_list and ' +
                    'mo #%d, since it is not given.' % i)
      ao_list = core.ao_creator(qc.geo_spec,qc.ao_spec)  
    else:
      display('\tComputing mo #%d, since it is not given.' % i)
    mo = core.mo_creator(ao_list,[qc.mo_spec[i]])[0]
  else:
    mo = mo_list[i]
  if delta_mo_list is None:
    if delta_ao_list is None:
      display('\tComputing delta_ao_list and the derivative of ' +
                    'mo #%d, since it is not given.' % j)
      delta_ao_list = core.ao_creator(qc.geo_spec,qc.ao_spec,drv=drv)
    else:
      display('\tComputing mo #%d, since it is not given.' % j)
    delta_mo = core.mo_creator(delta_ao_list,[qc.mo_spec[j]])[0]
  else:
    delta_mo = delta_mo_list[i]
  
  return mo*delta_mo
  # mo_transition_flux_density 
