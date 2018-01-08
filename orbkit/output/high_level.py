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

from .amira import amira_creator
from .cube import cube_creator
from .hdf5 import hdf5_creator, hdf5_append, hdf5_write, hx_network_creator
from .mayavi_interface import view_with_mayavi
from .pdb import pdb_creator
from .vmd import vmd_network_creator
from .xyz import xyz_creator
from .native import write_native

synonyms = {'h5':'h5', 'hdf5':'h5',
            'cube':'cb', 'cb':'cb',
            'am':'am', 
            'hx':'hx',
            'vmd':'vmd',
            'mayavi':'mayavi'
            }

def main_output(data,geo_info=None,geo_spec=None,outputname='new',otype='auto',
                drv=None,omit=[],datalabels='',**kwargs):
  '''Creates the requested output.
  
  **Parameters:**
  
  data : numpy.ndarray, shape=N, shape=((NDRV,) + N), shape=(n, (NDRV,) + N) or list of numpy.ndarrays
    Contains the output data. The shape (N) depends on the grid and the data, i.e.,
    3d for regular grid, 1d for vector grid. 
  geo_info, geo_spec : 
    See :ref:`Central Variables` for details.
  outputname : str or list of str
    Contains the base name of the output file.
  otype : str or list of str
    Contains the output file type. Possible options:
    'auto', 'h5', 'cb', 'am', 'hx', 'vmd', 'mayavi'
  drv : None, list of str or list of list of str optional
    If not None, a 4d(regular)/2d(vector) input data array will be expected
    with NDRV = len(drv).
  omit : list of str, optional
    If not empty, the output file types specified here are omitted.
  
  **Note:**
  
    All additional keyword arguments are forwarded to the output functions.
  '''

  if isinstance(otype, str):
    if otype == 'auto':
      filename, otype = path.splitext(outputname)
      otype = otype[1:]
    otype = [otype]
  else:
    for iot in range(len(otype)):
      if otype[iot] == 'auto':
        otype[iot] = otype[iot].split('.')[-1]

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
      group = ['' for _ in data]

    for i, oname in enumerate(outputname):
      write_native(data[i], oname, ftype[i])
      if group[i] != '':
        output_written.append('{0}_{1}.{2}'.format(oname, group[i], ftype[i]))
      else:
        output_written.append('{0}.{1}'.format(oname, ftype[i]))

  else:
    print_waring = False
    output_not_possible = (grid.is_vector and not grid.is_regular)
    
    # Shape shall be (Ndrv,Ndata,Nx,Ny,Nz) or (Ndrv,Ndata,Nxyz)
    data = numpy.array(data)
    dims = 1 if grid.is_vector else 3
    
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
    elif data.ndim > dims + 2:
      output_not_possible = True
      display('data.ndim > (ndim of grid) +2')
    
    if 'vmd' in otype and not 'cb' in otype:
      otype.append('cb')
    
    otype = [i for i in otype if i not in omit]
    otype_synonyms = [synonyms[i] for i in otype]
    
    if otype is None or otype == []:
      return output_written 

    # Convert the data to a regular grid, if possible
    is_regular_vector = (grid.is_vector and grid.is_regular)

    if is_regular_vector:  
      display('\nConverting the regular 1d vector grid to a 3d regular grid.')
      grid.vector2grid(*grid.N_)
      #new_data = grid.mv2g(data) #[]
      #for idata in range(data.shape[1]):
        #new_data.append(grid.mv2g(data=data[:,idata]))
      data = numpy.array(grid.mv2g(data))

    isstr = isinstance(outputname, str)
    isstr_datalabels = isinstance(datalabels, str)
    
    if drv is not None:
      fid = '%(f)s_d%(d)s'
      datalabel_id = 'd/d%(d)s %(f)s'
      it = enumerate(drv)
    else:
      fid = '%(f)s'
      datalabel_id = '%(f)s'
      it = [(0,None)]
    
    
    if 'h5' in otype: 
      for idrv,jdrv in it:
        f = fid % {'f':outputname if isstr else outputname[-1],'d':jdrv}
        display('\nSaving to Hierarchical Data Format file (HDF5)...' +
                '\n\t{0}.h5'.format(f))
        hdf5_creator(data[idrv],f,geo_info,geo_spec,**kwargs)
        output_written.append('{0}.h5'.format(f))
    
    cube_files = []
    all_datalabels = []
    for idrv,jdrv in it:
      for idata in range(data.shape[1]):
        if isstr:
          f = {'f': outputname + '_' + str(idata) if data.shape[1] > 1 else outputname,
               'd':jdrv}
        else:
          f = {'f': outputname[idata], 'd':jdrv}
        if isstr_datalabels:
          c = {'f': str(idata) + ',' + datalabels if data.shape[1] > 1 else datalabels,
               'd':jdrv}
        else: 
          c = {'f': datalabels[idata],'d':jdrv}
        datalabel = datalabel_id%c
        all_datalabels.append(datalabel)
        
        if 'am' in otype or 'hx' in otype and not print_waring:
          if output_not_possible: print_waring = True
          else: 
            display('\nSaving to ZIBAmiraMesh file...' +
                           '\n\t%(o)s.am' % {'o': fid % f})
            amira_creator(data[idrv,idata],(fid % f))
            output_written.append('%s.am' % (fid % f))
        if 'hx' in otype and not print_waring:
          if output_not_possible: print_waring = True
          else: 
            # Create Amira network incl. Alphamap
            display('\nCreating ZIBAmira network file...')
            hx_network_creator(data[idrv,idata],(fid % f))
            output_written.append('%s.hx' % (fid % f))
        if 'cb' in otype or 'vmd' in otype and not print_waring:
          if output_not_possible: print_waring = True
          else: 
            display('\nSaving to .cb file...' +
                            '\n\t%(o)s.cb' % {'o': fid % f})
            cube_creator(data[idrv,idata],(fid % f),geo_info,geo_spec,
                         comments=datalabel,
                         **kwargs)
            output_written.append('%s.cb' % (fid % f))
            cube_files.append('%s.cb' % (fid % f))

    if 'vmd' in otype and not print_waring:
      if output_not_possible: print_waring = True
      else: 
        # Create VMD network 
        display('\nCreating VMD network file...' +
                        '\n\t%(o)s.vmd' % {'o': fid % f})        
        vmd_network_creator(outputname if isstr else outputname[-1],
                            cube_files=cube_files,**kwargs)
        output_written.append('%s.vmd' % (fid % f))

    if 'mayavi' in otype:
      data = data.reshape((-1,)+grid.get_shape())
      if output_not_possible: print_waring = True
      else: view_with_mayavi(grid.x,grid.y,grid.z,data,geo_spec=geo_spec,datalabels=all_datalabels,**kwargs)
        
    if print_waring:
      display('For a non-regular vector grid (`if grid.is_vector and not grid.is_regular`)')
      display('only HDF5 is available as output format...')
      display('Skipping all other formats...')
    
    if is_regular_vector:
      # Convert the data back to a regular vector grid
      grid.grid2vector()
  
  return output_written

