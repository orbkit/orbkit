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

from .tools import *

# Import orbkit modules
from orbkit import core,grid,options
from orbkit.output import main_output
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
    drv : string or list of strings {None,'x','y', 'z', 'xx', 'xy', ...}, optional
      If not None, a derivative calculation of the molecular orbitals 
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
  mo_spec = qc.mo_spec[fid_mo_list]
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
    if '@' in options.outputname:
      outputname,group = options.outputname.split('@')
    else:
      outputname,group = options.outputname, ''
    outputname,autootype = os.path.splitext(outputname)
    ofid = '%s_MO%s@%s' % (outputname,autootype,group)
  if not options.no_output:
    output_written = main_output(mo_list,
                                 qc_select,
                                 outputname=ofid,
                                 datalabels=qc_select.mo_spec.get_labels(),
                                 dataindices=qc_select.mo_spec.get_indices(),
                                 otype=otype,
                                 drv=drv)
  
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
    drv : string or list of strings {None,'x','y', 'z', 'xx', 'xy', ...}, optional
      If not None, a derivative calculation is requested.
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
  # For later iteration we'll make it into a list here if it is not
  mo_info_list = qc.mo_spec.select(fid_mo_list, flatten_input=False)
    
  drv = options.drv if drv is None else drv
  laplacian = options.laplacian if laplacian is None else laplacian
  slice_length = options.slice_length if slice_length is None else slice_length
  numproc = options.numproc if numproc is None else numproc
  
  if ofid is None:
    ofid = options.outputname
  
  datasets = []
  datalabels = []
  delta_datasets = []
  delta_datalabels = []
  cube_files = []
  
  for i_file, mo_info in enumerate(mo_info_list):
    qc_select = qc.copy()
    qc_select.mo_spec = mo_info
    label = 'mo_set:'+mo_info.selection_string
    display('\nStarting with the molecular orbital list \n\t'
            + label +
            '\n\t(Only regarding existing and occupied mos.)\n')
    
    data = core.rho_compute(qc_select,
                            drv=drv,
                            laplacian=laplacian,
                            slice_length=slice_length,
                            numproc=numproc)

    if drv is None:
      rho = data
      delta_datasets = numpy.zeros((0,)+rho.shape)
    elif laplacian:
      rho, delta_rho, laplacian_rho = data 
      delta_datasets.extend(delta_rho)
      delta_datasets.append(laplacian_rho)
      delta_datalabels.extend(['d^2/d%s^2 %s' % (i,label) for i in 'xyz'])
      delta_datalabels.append('laplacian_of_' + label)
    else:
      rho, delta_rho = data
      delta_datasets.extend(delta_rho)
      delta_datalabels.extend(['d/d%s %s' % (i,label) for i in drv])
      
    
    datasets.append(rho)
    datalabels.append(label)
  
  datasets = numpy.array(datasets)
  delta_datasets = numpy.array(delta_datasets)
  delta_datalabels.append('mo_set')
  data = numpy.append(datasets,delta_datasets,axis=0)
  
  if not options.no_output:
    output_written = main_output(data,
                                 qc,
                                 outputname=ofid,
                                 datalabels=datalabels+delta_datalabels,
                                 otype=otype,
                                 drv=None)
  return data
  # mo_set 


def calc_ao(qc, drv=None, otype=None, ofid=None,
            numproc=None, slice_length=None):
  '''Computes and saves all atomic orbital or a derivative thereof.
  
  **Parameters:**
   
    qc.geo_spec, qc.geo_info, qc.ao_spec, qc.mo_spec :
      See :ref:`Central Variables` for details.
  drv : string or list of strings {None,'x','y', 'z', 'xx', 'xy', ...}, optional
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
    ao_list : numpy.ndarray, shape=((NAO,) + N)
      Contains the computed NAO atomic orbitals on a grid. 
  '''
  
  slice_length = options.slice_length if slice_length is None else slice_length
  numproc = options.numproc if numproc is None else numproc
  datalabels = qc.ao_spec.get_labels()
  
  ao_list = core.rho_compute(qc,
                             calc_ao=True,
                             drv=drv,
                             slice_length=options.slice_length,
                             numproc=options.numproc)
  
  if otype is None:
    return ao_list
  
  if ofid is None:
    ofid = '%s_AO' % (options.outputname)
  
  if not options.no_output:
    output_written = main_output(ao_list,
                                 qc,
                                 outputname=ofid,
                                 datalabels=qc.ao_spec.get_labels(),
                                 otype=otype,
                                 drv=drv)
  
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


def calc_jmo(qc, ij, drv=['x','y','z'], numproc=1, otype=None, ofid='',**kwargs):
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
  ij = numpy.asarray(ij)
  if ij.ndim == 1 and len(ij) == 2:
    ij.shape = (1,2)
  assert ij.ndim == 2
  assert ij.shape[1] == 2
  
  u, indices = numpy.unique(ij,return_inverse=True)
  indices.shape = (-1,2)

  mo_spec = qc.mo_spec[u]
  qc_select = qc.copy()
  qc_select.mo_spec = mo_spec
  labels = mo_spec.get_labels(format='short')
  
  mo_matrix = core.calc_mo_matrix(qc_select,drv=drv,
                                  numproc=numproc,**kwargs)
  jmo = numpy.zeros((len(drv),len(indices)) + grid.get_shape())
  datalabels = []
  for n,(i,j) in enumerate(indices):
    jmo[:,n] = - 0.5 * (mo_matrix[:,i,j] - mo_matrix[:,j,i])
    datalabels.append('j( %s , %s )'%(labels[i],labels[j]))
  
  if not options.no_output:
    output_written = main_output(jmo,
                                 qc,
                                 outputname=ofid,
                                 datalabels=datalabels,
                                 otype=otype,
                                 drv=drv)
  return jmo