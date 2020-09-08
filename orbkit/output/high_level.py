# -*- coding: iso-8859-1 -*-
'''Module for creating the requested output files.
'''
'''
orbkit
Gunter Hermann, Vincent Pohl, Lukas Eugen Marsoner Steinkasserer, Axel Schild, and Jean Christophe Tremblay

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
import numpy
import h5py
from os import path

# Import orbkit modules
from orbkit import grid, options
from orbkit.display import display
from orbkit.units import a0_to_aa

from .amira import amira_creator, hx_network_creator
from .cube import cube_creator
from .hdf5 import hdf5_creator, hdf5_append, hdf5_write
from .mayavi_interface import view_with_mayavi
from .pdb import pdb_creator
from .vmd import vmd_network_creator
from .xyz import xyz_creator
from .native import write_native

synonyms = {'auto':'auto',
            'h5':'h5', 'hdf5':'h5',
            'npz':'npz','numpy':'npz',
            'cube':'cube', 'cb':'cube',
            'cube.gz':'cube', 'cb.gz':'cube',
            'am':'am', 
            'hx':'hx',
            'vmd':'vmd',
            'mayavi':'mayavi',
            '':None
            }

def main_output(data,qc=None,outputname='data',otype='auto',gname='',
                drv=None,omit=[],datalabels='',dataindices=None,mode='w',**kwargs):
  '''Creates the requested output.
  
  **Parameters:**
  
  data : numpy.ndarray, shape=N, shape=((NDRV,) + N), shape=(n, (NDRV,) + N) or list of numpy.ndarrays
    Contains the output data. The shape (N) depends on the grid and the data, i.e.,
    3d for regular grid, 1d for vector grid. 
  qc : class or dict
    QCinfo class or dictionary containing the following attributes/keys.
    See :ref:`Central Variables` for details.
  outputname : str or list of str
    Contains the base name of the output file. If outputname contains @, string will be split and first 
    part interpreted as outputname and second as gname (cf. Parameters:gname). 
  otype : str or list of str, optional
    Contains the output file type. Possible options:
    'auto', 'h5', 'cb', 'am', 'hx', 'vmd', 'mayavi'
    
    If otype='native', a native input file will be written, the type of which may be specifie by
    ftype='numpy'.
  gname : str, optional
    For native, HDF5, or npz output, specifies the group, where the data will be stored.
  drv : None, list of str or list of list of str, optional
    If not None, a 4d(regular)/2d(vector) input data array will be expected
    with NDRV = len(drv). Specifies the file labels, i.e. e.g., data_d{drv}.cube for 4d array.
    For 5d arrays i.e., data_0_d{drv}.cube
  datalabels : list of str, optional
    If not empty, the output file types specified here are omitted.
  dataindices : list of int, optional
    If not empty, the output file types specified here are omitted.
  omit : list of str, optional
    If not empty, the output file types specified here are omitted.
  mode : str={'r', 'w', 'a'}, optional
    Specifies the mode used to open the file (native, HDF5, or npz). 
  
  **Note:**
  
    All additional keyword arguments are forwarded to the output functions.
  '''
  if otype is None or otype == []:
    return []
  
  if isinstance(outputname, str):
    if '@' in outputname:
      outputname,gname = outputname.split('@')
  if isinstance(otype, str):
    if otype == 'auto':
      outputname, otype = path.splitext(outputname)
      otype = otype[1:]
    otype = [otype]
  elif isinstance(otype, list) and len(otype) == 1:
    if otype[0] == 'auto':
      outputname, otype = path.splitext(outputname)
      otype = [otype[1:]]
  else:
    for iot in range(len(otype)):
      if otype[iot] == 'auto':
        outputname, tmp = path.splitext(outputname)
        if tmp != '':
          otype[iot] = tmp[1:]
  
  # Catch our native format before all else
  # We can't figure this out by the file ending alone
  # as we support hdf5 for output of both grid-based data
  # as well as our internal format

  output_written = []
  internals = [i for i in range(len(otype)) if otype[i] == 'native']
  
  if len(internals) > 0:
    if not isinstance(data,list):
      data = [data]

    if isinstance(outputname,str):
      outputname = [outputname for _ in data]

    if 'ftype' in kwargs.keys():
      if isinstance(kwargs['ftype'],str):
        ftype = [kwargs['ftype'] for _ in data]
    else:
      ftype = ['numpy' for _ in data]

    if 'group' in kwargs.keys():
      if isinstance(kwargs['group'],str):
        group = [kwargs['group'] for _ in data]
    else:
      group = [i.__class__.__name__.lower() for i in data]
    
    display('Writing native input file...' )
    for i, oname in enumerate(outputname):
      output_written.append(write_native(data[i], oname, ftype[i], mode=mode,
                                         gname=path.join(gname,group[i])))
    display('\n'.join(['\t' + i for i in output_written]))

  else:
    print_waring = False
    output_not_possible = (grid.is_vector and not grid.is_regular)
    
    # Shape shall be (Ndrv,Ndata,Nx,Ny,Nz) or (Ndrv,Ndata,Nxyz)
    data = numpy.array(data)
    dims = 1 if grid.is_vector else 3
    shape = data.shape
    
    if drv is not None and isinstance(drv,str):
      drv = [drv]
    
    if data.ndim < dims:
      output_not_possible = True
      display('data.ndim < ndim of grid')
    elif data.ndim == dims: # 3d data set
      data = data[numpy.newaxis,numpy.newaxis]
    elif data.ndim == dims + 1: # 4d data set
      if drv is not None:
        data = data[:,numpy.newaxis]
      else:
        data = data[numpy.newaxis]
    elif data.ndim == dims + 2: # 5d data set check if drv matches Ndrv
      if drv is None or len(drv) != data.shape[0]:
        drv = list(range(data.shape[0]))
    elif data.ndim > dims + 2:
      output_not_possible = True
      display('data.ndim > (ndim of grid) +2')
    
    if 'vmd' in otype and not ('cb' in otype or 'cube' in otype):
      otype.append('cube')
    if 'hx' in otype and not 'am' in otype:
      otype.append('am')
    
    otype = [i for i in otype if i not in omit]
    otype_synonyms = [synonyms[i] for i in otype]
    otype_ext = dict(zip(otype_synonyms,otype))
    
    # Convert the data to a regular grid, if possible
    is_regular_vector = (grid.is_vector and grid.is_regular)

    if is_regular_vector:  
      display('\nConverting the regular 1d vector grid to a 3d regular grid.')
      grid.vector2grid(*grid.N_)
      data = numpy.array(grid.mv2g(data))
    
    isstr = isinstance(outputname, str)
    if isinstance(datalabels, str):
      if data.shape[1] > 1:
        datalabels = numpy.array([str(idata) + ',' + datalabels 
                                  for idata in range(data.shape[1])])
      else:
        datalabels = numpy.array([datalabels])
    elif isinstance(datalabels, list):
      datalabels = numpy.array(datalabels)
    
    if drv is not None:
      fid = '%(f)s_d%(d)s.'
      datalabel_id = 'd/d%(d)s %(f)s'
      contents = {
        'axis:0': numpy.array(['d/d%s' % i if i is not None else str(i) 
                               for i in drv]),
        'axis:1': datalabels}
      it = enumerate(drv)
    elif data.shape[0] > 1:
      fid = '%(f)s_%(d)s.'
      datalabel_id = '%(d)s %(f)s'
      it = enumerate(data.shape[0])
      contents = {
        'axis:0': numpy.arange(data.shape[0]).astype(str),
        'axis:1': datalabels}
    else:
      fid = '%(f)s.'
      datalabel_id = '%(f)s'
      it = [(0,None)]
      if data.shape[1] > 1:
        contents = {'axis:0': datalabels}
      else:
        contents = datalabels
    
    cube_files = []
    all_datalabels = []
    for idrv,jdrv in it:
      datasetlabels = []
      for idata in range(data.shape[1]):
        if isstr:
          index = str(idata) if dataindices is None else str(dataindices[idata])
          f = {'f': outputname + '_' + index if data.shape[1] > 1 else outputname,
               'd':jdrv}
        else:
          f = {'f': outputname[idata], 'd':jdrv}
        c = {'f': datalabels[idata],'d':jdrv}
        datalabel = datalabel_id%c
        datasetlabels.append(datalabel)
        
        if 'am' in otype_synonyms and not print_waring:
          if output_not_possible: print_waring = True
          else: 
            filename = fid % f + otype_ext['am']
            display('\nSaving to ZIBAmiraMesh file...\n\t' + filename)
            amira_creator(data[idrv,idata],filename)
            output_written.append(filename)
        if 'hx' in otype_synonyms and not print_waring:
          if output_not_possible: print_waring = True
          else: 
            filename = fid % f + otype_ext['hx']
            display('\nCreating ZIBAmira network file...\n\t' + filename)
            hx_network_creator(data[idrv,idata],filename)
            output_written.append(filename)
        if 'cube' in otype_synonyms and not print_waring:
          if output_not_possible: print_waring = True
          elif qc is None: 
            display('\nFor cube file output `qc` is a required keyword parameter in `main_output`.')
          else: 
            filename = fid % f + otype_ext['cube']
            display('\nSaving to cube file...\n\t' + filename)
            cube_creator(data[idrv,idata],filename,qc.geo_info,qc.geo_spec,
                         comments=datalabel,
                         **kwargs)
            output_written.append(filename)
            cube_files.append(filename)
        
      all_datalabels.extend(datasetlabels)
    
    if 'vmd' in otype_synonyms and not print_waring:
      if output_not_possible: print_waring = True
      else: 
        filename = (outputname if isstr else outputname[-1]) +'.'+ otype_ext['vmd']
        display('\nCreating VMD network file...\n\t' + filename)
        vmd_network_creator(filename,cube_files=cube_files,**kwargs)
        output_written.append(filename)

    
    if 'h5' in otype_synonyms: 
      filename = (outputname if isstr else outputname[-1]) +'.'+ otype_ext['h5']
      display('\nSaving to Hierarchical Data Format file (HDF5)...\n\t' + filename)
      
      hdf5_creator(data.reshape(shape),filename,qcinfo=qc,gname=gname,
                   ftype='hdf5',contents=contents,mode=mode,**kwargs)
      output_written.append(filename)
    
    if 'npz' in otype_synonyms: 
      filename = (outputname if isstr else outputname[-1]) 
      display('\nSaving to a compressed .npz archive...\n\t' + filename+'.npz')
      hdf5_creator(data.reshape(shape),filename,qcinfo=qc,gname=gname,
                   ftype='numpy',contents=contents,mode=mode,**kwargs)
      output_written.append(filename)
    
    if 'mayavi' in otype_synonyms:
      if output_not_possible: print_waring = True
      else: 
        display('\nDepicting the results with MayaVi...\n\t')
        if drv == ['x','y','z'] or drv == [0,1,2]:
          is_vectorfield = True
          data = numpy.swapaxes(data,0,1)
          datalabels = datalabels
        else:
          is_vectorfield = False
          data = data.reshape((-1,)+grid.get_shape())
          datalabels = all_datalabels
        
        view_with_mayavi(grid.x,grid.y,grid.z,
                         data,
                         is_vectorfield=is_vectorfield,
                         geo_spec=qc.geo_spec,
                         datalabels=datalabels,**kwargs)
        
    if print_waring:
      display('For a non-regular vector grid (`if grid.is_vector and not grid.is_regular`)')
      display('only HDF5 is available as output format...')
      display('Skipping all other formats...')
    
    if is_regular_vector:
      # Convert the data back to a regular vector grid
      grid.grid2vector()
  
  return output_written

