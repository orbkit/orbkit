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

# Import orbkit modules
from orbkit import grid, options
from orbkit.display import display
from orbkit.units import a02aa

from .amira import amira_creator
from .cube import cube_creator
from .hdf5 import hdf5_creator, hdf5_append, hdf5_write, hx_network_creator
from .mayavi import view_with_mayavi
from .pdb import pdb_creator
from .vmd import vmd_network_creator
from .xyz import xyz_creator

def main_output(data,geo_info,geo_spec,outputname='new',otype='auto',
                drv=None,omit=[],**kwargs):
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

  print_waring = False
  output_written = []
  data = numpy.array(data)
  if not drv:
    ex_dim = 0
  else:
    ex_dim = 1
  if grid.is_vector:
    if data.ndim == ex_dim + 1:
      data = data[numpy.newaxis]
  else:
    if data.ndim == ex_dim + 3:
      data = data[numpy.newaxis]

  if isinstance(otype, str):
    if otype == 'auto':
      otype = otype.split('.')[-1]
    otype = [otype]
  else:
    for iot in range(len(otype)):
      if otype[iot] == 'auto':
        otype[iot] = otype[iot].split('.')[-1]
  
  if 'vmd' in otype and not 'cb' in otype:
    otype.append('cb')
  
  otype = [i for i in otype if i not in omit]
  
  if otype is None or otype == []:
    return output_written 

  # Convert the data to a regular grid, if possible
  output_not_possible = (grid.is_vector and not grid.is_regular)
  is_regular_vector = (grid.is_vector and grid.is_regular)

  if is_regular_vector:  
    display('\nConverting the regular 1d vector grid to a 3d regular grid.')
    grid.vector2grid(*grid.N_)
    new_data = []
    for idata in range(data.shape[0]):
      new_data.append(grid.mv2g(data=data[idata]))
    data = numpy.array(new_data)

  if isinstance(outputname, str):
    isstr = True

  if 'h5' in otype:
    if not isstr:
      display('Recived list of output names - writing each dateset' + 
              '\n\t to separate HDF5 file')
      for idata in range(data.shape[0]):
        fid = '%(f)s'
        it = [(0,None)]

        display('\nSaving to Hierarchical Data Format file (HDF5)...' +
                '\n\t{0}.h5'.format(fid))
        hdf5_creator(data,fid,geo_info,geo_spec,**kwargs)
        output_written.append('{0}.h5'.format(fid))

    else:
      display('\nSaving to Hierarchical Data Format file (HDF5)...' +
              '\n\t{0}.h5'.format(outputname))
      hdf5_creator(data,outputname,geo_info,geo_spec,**kwargs)
      output_written.append('{0}.h5'.format(outputname))

  for idata in range(data.shape[0]):
    if drv is not None:
      fid = '%(f)s_d%(d)s'
      it = enumerate(drv)
    else:
      fid = '%(f)s'
      it = [(0,None)]

    if isstr:
      f = {'f': outputname + '_' + str(idata)}
    else:
      f = {'f': outputname[idata]}

    for i,j in it:
      f['d'] = j
      if 'am' in otype or 'hx' in otype and not print_waring:
        if output_not_possible: print_waring = True
        else: 
          display('\nSaving to ZIBAmiraMesh file...' +
                         '\n\t%(o)s.am' % {'o': fid % f})
          amira_creator(data[idata],(fid % f))
          output_written.append('%s.am' % (fid % f))
      if 'hx' in otype and not print_waring:
        if output_not_possible: print_waring = True
        else: 
          # Create Amira network incl. Alphamap
          display('\nCreating ZIBAmira network file...')
          hx_network_creator(data[idata],(fid % f))
          output_written.append('%s.hx' % (fid % f))
      if 'cb' in otype or 'vmd' in otype and not print_waring:
        if output_not_possible: print_waring = True
        else: 
          display('\nSaving to .cb file...' +
                          '\n\t%(o)s.cb' % {'o': fid % f})
          cube_creator(data[idata],(fid % f),geo_info,geo_spec,**kwargs)
          output_written.append('%s.cb' % (fid % f))
         #else: output_creator(data[idata],(fid % f),geo_info,geo_spec)  # Axel's cube files

  if 'vmd' in otype and not print_waring:
    if output_not_possible: print_waring = True
    else: 
      # Create VMD network 
      display('\nCreating VMD network file...' +
                      '\n\t%(o)s.vmd' % {'o': fid % f})        
      vmd_network_creator((fid % f),cube_files=['%s.cb' % (fid % f)],**kwargs)
      output_written.append('%s.vmd' % (fid % f))

  if 'mayavi' in otype:
    if output_not_possible: print_waring = True
    else: view_with_mayavi(grid.x,grid.y,grid.z,data,geo_spec=geo_spec,**kwargs)
      
  if print_waring:
    display('For a non-regular vector grid (`if grid.is_vector and not grid.is_regular`)')
    display('only HDF5 is available as output format...')
    display('Skipping all other formats...')
  
  if is_regular_vector:
    # Convert the data back to a regular vector grid
    grid.grid2vector()
  
  return output_written

