# -*- coding: iso-8859-1 -*-
'''Module for creating the requested output files.
'''
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
import numpy
import h5py

# Import orbkit modules
from orbkit import grid, options
from orbkit.display import display

def main_output(data,geo_info,geo_spec,outputname='new',otype='h5',
                drv=None,omit=[],**kwargs):
  '''Creates the requested output.
  
  **Parameters:**
  
  data : numpy.ndarray, shape=N or shape=((NDRV,) + N)
    Contains the output data. The shape (N) depends on the grid and the data, i.e.,
    3d for regular grid, 1d for vector grid. 
  geo_info, geo_spec : 
    See :ref:`Central Variables` for details.
  outputname : str
    Contains the base name of the output file.
  otype : str or list of str
    Contains the output file type. Possible options:
    'h5', 'cb', 'am', 'hx', 'vmd', 'mayavi'
  drv : None or list of str, optional
    If not None, a 4d(regular)/2d(vector) input data array will be expected
    with NDRV = len(drv).
  omit : list of str, optional
    If not empty, the input file types specified here are omitted.
  
  **Note:**
  
    All additional keyword arguments are forwarded to the output functions.
  '''
  print_waring = False
  output_written = []
  if isinstance(otype,str):
    otype = [otype]
  
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
    data = grid.mv2g(data=data)
  
  if 'mayavi' in otype:
    if output_not_possible: print_waring = True
    else: view_with_mayavi(grid.x,grid.y,grid.z,data,geo_spec=geo_spec,**kwargs)
  
  if drv is not None:
    fid = '%(f)s_d%(d)s'
    it = enumerate(drv)
  else:
    fid = '%(f)s'
    it = [(0,None)]
    data = [data]
  f = {'f': outputname}
  
  for i,j in it:
    f['d'] = j
    d = data[i]
    if 'h5' in otype:
      display('\nSaving to Hierarchical Data Format file (HDF5)...' +
              '\n\t%(o)s.h5' % {'o': fid % f})
      HDF5_creator(d,(fid % f),geo_info,geo_spec,**kwargs)
      output_written.append('%s.h5' % (fid % f))
    if 'am' in otype or 'hx' in otype and not print_waring:
      if output_not_possible: print_waring = True
      else: 
        display('\nSaving to ZIBAmiraMesh file...' +
                     '\n\t%(o)s.am' % {'o': fid % f})
        amira_creator(d,(fid % f))
        output_written.append('%s.am' % (fid % f))
    if 'hx' in otype and not print_waring:
      if output_not_possible: print_waring = True
      else: 
        # Create Amira network incl. Alphamap
        display('\nCreating ZIBAmira network file...')
        hx_network_creator(data,(fid % f))
        output_written.append('%s.hx' % (fid % f))
    if 'cb' in otype or 'vmd' in otype and not print_waring:
      if output_not_possible: print_waring = True
      else: 
        display('\nSaving to .cb file...' +
                      '\n\t%(o)s.cb' % {'o': fid % f})
        cube_creator(d,(fid % f),geo_info,geo_spec,**kwargs)
        output_written.append('%s.cb' % (fid % f))
       #else: output_creator(d,(fid % f),geo_info,geo_spec)  # Axel's cube files
    if 'vmd' in otype and not print_waring:
      if output_not_possible: print_waring = True
      else: 
        # Create VMD network 
        display('\nCreating VMD network file...' +
                      '\n\t%(o)s.vmd' % {'o': fid % f})        
        vmd_network_creator((fid % f),cube_files=['%s.cb' % (fid % f)],**kwargs)
        output_written.append('%s.vmd' % (fid % f))
      
  if print_waring:
    display('For a non-regular vector grid (`if grid.is_vector and not grid.is_regular`)')
    display('only HDF5 is available as output format...')
    display('Skipping all other formats...')
  
  if is_regular_vector:
    # Convert the data back to a regular vector grid
    grid.grid2vector()
  
  return output_written

def cube_creator(rho,filename,geo_info,geo_spec,comments='',**kwargs):
  '''Creates a plain text Gaussian cube file. 
  
  **Parameters:**
  
  rho : numpy.ndarray, shape=N
    Contains the output data.
  filename : str
    Contains the base name of the output file.
  geo_info, geo_spec : 
    See :ref:`Central Variables` for details.
  comments : str, optional
    Specifies the second (comment) line of the cube file.  
  '''
  
  # Open an empty file 
  fid = open('%(f)s.cb' % {'f': filename}, 'w')
  
  # Write the type and the position of the atoms in the header 
  string = 'orbkit calculation\n'
  string += ' %(f)s\n'  % {'f': comments}
  # How many atoms 
  string += ('%(at)d' % {'at': len(geo_info)}).rjust(5)
  # Minima
  for ii in range(3):
    string += ('%(min)0.6f' % {'min': grid.min_[ii]}).rjust(12)

  for ii in range(3):
    string += '\n'
    string += ('%(N)d' % {'N': grid.N_[ii]}).rjust(5)
    for jj in range(3):
      if jj == ii: 
        string += ('%(dr)0.6f' % {'dr': grid.delta_[ii]}).rjust(12)
      else:
        string += ('%(dr)0.6f' % {'dr': 0}).rjust(12)
  
  for ii in range(len(geo_info)):
    string += '\n'
    string += ('%(N)d' % {'N': round(float(geo_info[ii][2]))}).rjust(5)
    string += ('%(ch)0.6f' % {'ch': float(geo_info[ii][1])}).rjust(12)
    for jj in range(3):
      string += ('%(r)0.6f' % {'r': geo_spec[ii][jj]}).rjust(12)
  string += '\n'
  for rr in range(len(grid.x)):
    for ss in range(len(grid.y)):
      for tt in range(len(grid.z)):
        string += ('%(rho).6E' % {'rho': rho[rr,ss,tt]}).rjust(13)
        if (tt % 6 == 5): 
          string += '\n'
      string += '\n'
  
  
  fid.write(string)
  
  # Close the file 
  fid.close()

def vmd_network_creator(filename,cube_files=None,render=False,iso=(-0.01,0.01),
                        abspath=False,**kwargs):
  '''Creates a VMD script file from a list of cube files provided.
  
  **Parameters:**
  
  filename : str
    Contains the base name of the output file.
  cube_files : None or list of str
    Specifies the cube files which serve as input for the VMD script.
    If None, searches the directory for '.cb' and '.cube' files.
  render : bool
    If True, the VMD script will automatically create '.tga' files for each 
    cube file.
  iso : tuple
    Specifies the isovalue for the blue and the red isosurface, respectively.
  abspath : bool
    If True, the paths of the cube files will be expanded to absolute file paths.
  '''
  from os import path,listdir
  import linecache
  from orbkit import vmd_network_draft
  if cube_files is None:
    display('No list of cube (.cb or .cube) filenames provided. Checking the directory' + 
            ' of the outputfile...')
    cube_files = []
    for fid in listdir(path.dirname(filename)):
      if fid.endswith('.cb') or fid.endswith('.cube'):
        cube_files.append(fid)
    if cube_files == []:
      raise IOError('Could not find valid cube files in %s' % path.dirname(filename))
  elif isinstance(cube_files,str):
    cube_files = [cube_files]
  elif not isinstance(cube_files,list):
    raise IOError('`cube_files` has to be a list of strings.')
  
  title = []
  mo = ''
  for i,f in enumerate(cube_files):
    title = linecache.getline(f,2)
    if title.split() == []:
      title = path.splitext(path.basename(f))[0]
    else:
      title = title.replace('\n','').replace(' ','')
    linecache.clearcache()
    pid = path.abspath(f) if abspath else path.relpath(f,path.dirname(filename))
    mo += vmd_network_draft.mo_string % {
                                  'c': i, 
                                  'n1': pid,
                                  'n2': title, 
                                  'isored': iso[0], 
                                  'isoblue': iso[1], 
                                  'render': '' if render else '#'
                                  }
  
  f = open('%(f)s.vmd' % {'f': filename},'w')
  f.write(vmd_network_draft.vmd_string % {'mo':mo})
  f.close()

def hdf5_open(fid,mode='w'):
  '''Open an HDF5 file.
  
  **Parameters:**
  
  fid : str
    Contains the filename of the HDF5 file.
  mode : str, optional
    Specifies the mode used to open the file. ('r', 'w', 'a')
  
  **Usage:**
  
  for HDF5_file in hdf5_open(fid,mode='w'):
    # Do something with HDF5_file...
  '''
  import h5py
  HDF5_file = h5py.File(fid, mode)
  try:
    yield HDF5_file
  finally:
    HDF5_file.close()

def hdf5_append(x,group,name='data'):
  '''Automatically append data to an open HDF5 file/group.
  
  **Parameters:**
  
  s : numpy.ndarray, list, dict, int, float, str
    Input data. Not supported: None
  group : h5py.File or h5py.Group
    The HDF5 file/group where the data will be appended.
  name : string, optional
    Specifies the group/dataset name in the HDF5 file/group. If empty,
    root directory of HDF5 file/group is chosen.  
  '''
  if isinstance(x,numpy.ndarray):
    if x.dtype.type is numpy.unicode_:
      x = numpy.asarray(x,dtype=numpy.string_)
    h5_dset = group.create_dataset(name,numpy.shape(x),data=x)
  elif isinstance(x,list):
    if name != '':
      subgroup = group.create_group(name)
    else:
      subgroup = group
    num = len(x)
    subgroup.attrs['num'] = num
    keys = range(num)
    for ii in keys:
      hdf5_append(x[ii],subgroup,name='%d' % ii)
  elif isinstance(x,dict):
    if name != '':
      subgroup = group.create_group(name)
    else:
      subgroup = group
    keys = x.keys()
    for ii in keys:
      hdf5_append(x[ii],subgroup,name=ii)
  else:
    group.attrs[name] = x
  

def hdf52dict(group,HDF5_file):
  '''Automatically convert the data stored in an HDF5 file/group 
  to a python dictionary.
    
  **Parameters:**
  
  group : str
    Specifies the group/dataset name in the HDF5 file/group.
  HDF5_file : h5py.File or h5py.Group
    The source HDF5 file/group.
  '''
  try:
    # The selected group is a dataset 
    x = HDF5_file['%s' % group][()]
  except:
    # Read all members 
    members = list(HDF5_file['%(g)s' % {'g':group}])
    try:
      members = numpy.array(members, dtype = numpy.int64)
      x = []
      for mm in numpy.sort(members):
        x.append(hdf52dict('%(g)s/%(m)s' % {'g':group, 'm': mm},HDF5_file))
    except (ValueError,OverflowError):
      x = {}
      for mm in members:
        x[mm] = hdf52dict('%(g)s/%(m)s' % {'g':group, 'm': mm},HDF5_file)
      attrs = HDF5_file['%(g)s' % {'g':group}].attrs
      for ii_a in attrs:
        x[ii_a] = attrs[ii_a]
  return x

def hdf5_write(fid,mode='w',gname='',**kwargs):
  '''Writes all keyword arguments to an HDF5 file.
  
  **Parameters:**
  
  fid : str
    Contains the filename of the HDF5 file.
  mode : str, optional
    Specifies the mode used to open the file. ('r', 'w', 'a')
  gname : str
    Specifies the group name in the HDF5 file, where the kwargs will be appended
    to.
  
  **Note:**
  
    All additional keyword arguments specify the input data.
  '''
  for f in hdf5_open(fid,mode=mode):
    if gname is not '':
      group = f.create_group(gname)
      for i,j in kwargs.items(): 
        group.create_dataset(i,numpy.shape(j),data=j)
    else:
      for i,j in kwargs.items(): 
        f.create_dataset(i,numpy.shape(j),data=j)

def HDF5_creator(data,filename,geo_info,geo_spec,mode='w',data_id='rho',group=None,
            data_only=False,ao_spec=None,mo_spec=None,is_mo_output=False,
            x=None,y=None,z=None,**kwargs):
  '''Creates an HDF5 file (Hierarchical Data Format) output.
  
  **Parameters:**
  
  data : numpy.ndarray
    Contains the output data.
  filename : str
    Contains the base name of the output file.
  geo_info, geo_spec : 
    See :ref:`Central Variables` for details.
  data_id : str, optional
    Specifies name of the dataset in the HDF5 file.
  group : None of str
    If not None, output data will be stored in this group.
  data_only : bool
    Specifies if only the dataset `data` should be saved.
  ao_spec, mo_spec : optional
    If not None, these data sets will be saved additionally.
    (cf. :ref:`Central Variables` for details)
    If mo_spec is not None, some information about the molecular orbitals 
    will be saved additionally. 
  is_mo_output : bool
    If True, information about the dataset will be saved to 'MO:Content'.
  x,y,z : numpy.ndarray, optional
    If not None, these variables will replace grid.x, grid.y, and/or grid.z,
    respectively.
  '''
  import h5py
  HDF5_file = h5py.File('%s.h5' % filename, mode)
  if group is None:
    f = HDF5_file
  else:
    HDF5_file.require_group(group)
    f = HDF5_file[group]
    
  if is_mo_output:    
    mo_name = []
    for i,j in enumerate(data):
      try:
        mo_name.append(mo_spec[i]['sym'])
      except TypeError:
        mo_name.append(str(i))
      
    dset = f.create_dataset('MO:Content',data=numpy.array(mo_name,dtype='S'))
  
  if options.z_reduced_density:
    dset = f.create_dataset(data_id,numpy.shape(data),data=data)
  else:
    dset = f.create_dataset(data_id,numpy.shape(data),data=data)
  
  if data_only:
    HDF5_file.close()
    return
  
  if x is None: x = grid.x
  if y is None: y = grid.y
  if z is None: z = grid.z
  
  if options.z_reduced_density:
    dset = f.create_dataset('z',(1,len(z)),data=z)
  else:
    dset = f.create_dataset('x',(1,len(x)),data=x)
    dset = f.create_dataset('y',(1,len(y)),data=y)
    dset = f.create_dataset('z',(1,len(z)),data=z)
  
  if ao_spec is not None:
    hdf5_append(ao_spec,f,name='ao_spec')
    if mo_spec is not None:
      hdf5_append(mo_spec,f,name='mo_spec')
  if mo_spec is not None:
    MO_info = f.create_group('MO_info')
    occ_num=[]
    energy=[]
    sym=[]
    for ii in range(len(mo_spec)):
      occ_num.append(mo_spec[ii]['occ_num'])
      energy.append(mo_spec[ii]['energy'])
      sym.append(numpy.string_(mo_spec[ii]['sym']))
    dset = MO_info.create_dataset('occ_num',((1,len(mo_spec))),data=occ_num)
    dset = MO_info.create_dataset('energy',((1,len(mo_spec))),data=energy)
    dset = MO_info.create_dataset('sym',((1,len(mo_spec))),data=sym)
  
  dset = f.create_dataset('geo_info',(numpy.shape(geo_info)),data=numpy.array(geo_info,dtype='S'))
  dset = f.create_dataset('geo_spec',(numpy.shape(geo_spec)),data=geo_spec)
  
  HDF5_file.close()

def hx_network_creator(rho,filename):
  '''Creates a ZIBAmira hx-network file including a colormap file (.cmap)
  adjusted to the density for the easy depiction of the density.
  '''
  from orbkit.hx_network_draft import hx_network
  # Create a .cmap colormap file using the default values 
  display('\tCreating ZIBAmira colormap file...\n\t\t%(f)s.cmap' % 
                {'f': filename})
  colormap_creator(rho,filename)
  
  # Create a .hx network file based on the file orbkit.hx_network_draft.py 
  display('\tCreating ZIBAmira network file...\n\t\t%(f)s.hx' % 
                {'f': filename})
  # Open an empty file
  fid = open('%(f)s.hx' % {'f': filename},'w')
  
  filename = filename.split('/')[-1]
  # Copy the content of the draft file and replace the keywords 
  fid.write(hx_network.replace("FILENAME",filename)) 
  
  # Close the file 
  fid.close()  

def amira_creator(rho,filename):
  '''Creates a ZIBAmira mesh file. (plain text)
  
  **Parameters:**
  
  rho : numpy.ndarray, shape=N
    Contains the output data.
  filename : str
    Contains the base name of the output file.
  '''
  # Open an empty file 
  fid = open('%(f)s.am' % {'f': filename},'w')

  # usage:
  #     - open Amira
  #     - left-click File -> Open Data
  #     - choose the sys.argv[2].am
  #     - Press OK
  #     - right-click on the new Data Object -> Compute -> Arithmetic
  #     - choose the Arithmetic object
  #     - select  Expr -> A and Result type -> regular and Apply
  #     - use the new data object to display your density as usual,
  
  # Write Header 
  fid.write('# AmiraMesh 3D ASCII 2.0\n\n\n')
  fid.write('define Lattice %(Nx)d %(Ny)d %(Nz)d\n' % 
                    {'Nx': grid.N_[0],'Ny': grid.N_[1],'Nz': grid.N_[2]})
  fid.write('define Coordinates %(N)d\n\n' % {'N': numpy.sum(grid.N_)})
  fid.write('Parameters {\n')
  fid.write('    Content "%(Nx)dx%(Ny)dx%(Nz)d float, uniform coordinates",\n' %
                    {'Nx': grid.N_[0],'Ny': grid.N_[1],'Nz': grid.N_[2]})
  fid.write('    BoundingBox %(xmin)f %(xmax)f %(ymin)f %(ymax)f %(zmin)f %(zmax)f,\n' %
            {'xmin': grid.min_[0],'xmax': grid.max_[0],
             'ymin': grid.min_[1],'ymax': grid.max_[1],
             'zmin': grid.min_[2],'zmax': grid.max_[2]})
  fid.write('    CoordType "uniform"\n}\n\n')
  fid.write('Lattice { float Data } @1\n')
  fid.write('# Data section follows\n@1\n')
  
  # Write density information to .am file
  string = ''
  for tt in range(len(grid.z)):
    for ss in range(len(grid.y)):
      for rr in range(len(grid.x)): 
        string += '%g\n' % rho[rr,ss,tt]
  
  fid.write(string)
  
  # Close the file 
  fid.close()

def colormap_creator(rho,filename,n_peaks=5,start=0.01,stop=0.999,peak_width=0.1):
  '''Creates a .cmap colormap for ZIBAmira adjusted to the density.

  **Default:** Isosurface values between 1% and 99.9% of the total density.
  '''
  
  # Where do we have the start and stop percentage? 
  rho_min, rho_max = determine_rho_range(rho,start=start,stop=stop)
  
  # Compute the distance between two isosurface values 
  delta_peak =(rho_max-rho_min)/(n_peaks-1)
  
  # Open a cmap file 
  fid = open('%(f)s.cmap' % {'f': filename}, 'w')
  
  # Write the header 
  fid.write('<!DOCTYPE Colormap>\n')
  fid.write('<ColormapVisage2.0 Name="%(f)s">\n' % {'f': filename})
  fid.write('  <Graph Active="1" Type="0" Name="">\n')
  
  # Initialize a counter for the contour values 
  counter = 0
  
  # Initialize a string for the contours 
  c_str = ('    <Control Opacity="%(o)f" Number="%(c)d" Blue="%(v)f"' + 
            ' Red="%(v)f" Green="%(v)f" Value="%(p)f"/>\n' )
  
  # Write the initial value at zero with zero opacity 
  fid.write(c_str % {'o': 0, 'c': counter, 'v': 0, 'p': 0})
  counter += 1
  
  # Loop over the contour values 
  for ii in range(n_peaks):
    # Calculate the value for the isosurface peak and two values 
    # next to the peak with zero opacity 
    peak = rho_min+delta_peak*(ii+1)
    peak_minus = peak * (1 - peak_width/2.)
    peak_plus = peak * (1 + peak_width/2.)
    
    # Calculate a value for the opacity and the color
    value = 1-(float(ii+1)/(n_peaks+1)*0.9)
    
    # Write the peak 
    # Left edge of the peak 
    fid.write(c_str % {'o': 0, 'c': counter, 'v': value, 'p': peak_minus})
    counter += 1
    # The peak 
    fid.write(c_str % {'o': 1-value, 'c': counter, 'v': value, 'p': peak})
    counter += 1
    # Right edge of the peak 
    fid.write(c_str % {'o': 0, 'c': counter, 'v': value, 'p': peak_plus})
    counter += 1
  
  # Write a final value at 1.5 * (final peak value) with zero opacity 
  fid.write(c_str % {'o': 0, 'c': counter, 'v': 0, 'p': peak_plus*1.5})
  
  # Finalize the file 
  fid.write('  </Graph>\n')
  fid.write('</ColormapVisage2.0>')
  
  # Close the file
  fid.close()

def determine_rho_range(rho,start=0.01,stop=0.999):
  '''Get a range for the isosurface values for the colormap_creator.'''
  # Sort the density values 
  a=numpy.reshape(rho,-1)
  a=a[numpy.argsort(a)]
  
  # Where do we have the start and stop percentage? 
  n1=int(start*len(a))
  n2=int(stop*len(a))
  rho_min=a[n1]
  rho_max=a[n2]  
  return rho_min, rho_max

def colormap_creator_peaks(filename,peaks,peak_width=0.02,peak_minus=None,
  peak_plus=None,alpha=0.2,rgb=0.2):  
  '''Creates a ZIBAmira colomap for selected data values.
  
  **Parameters:**
  
  filename : str
    Specifies the filename of the colormap.
  peaks : list
    Determines the values for peaks in the colormap.
  peak_width : float
    Specifies the width of of the peaks in the colormap.
  peak_min : None or float, optional
    Specifies the lower boundary of the colomap. (Peak with no hight.)
    If None, is set to the smallest data value minus 2*peak_width.
  peak_min : None or float, optional
    Specifies the upper boundary of the colomap. (Peak with no hight.)
    If None, is set to the larges data value plus 2*peak_width.
  alpha : float, data_range={0..1}
    Determines the opacity of the peak. (alpha = 1 - opacity)
  rgb : float or list or numpy.ndarray, data_range={0..1}
    If float or shape=(len(peaks),), specifies a the grey tone.
    Else, specifies the color in [red,green,blue].
  
  '''
  peaks = numpy.sort(peaks)
  if peak_minus is None:
    peak_minus = peaks[0]-(2*peak_width)
  if peak_plus is None:
    peak_plus = peaks[-1]+(2*peak_width)
  
  if isinstance(rgb,(float,int)):
    rgb = numpy.zeros((len(peaks),3)) + rgb
  else:
    rgb = numpy.array(rgb, dtype=float)
    if rgb.ndim == 1:
      if len(rgb) == 3:
        rgb = numpy.zeros((len(peaks),3)) + rgb[numpy.newaxis,:]
      elif len(rgb) == len(peaks):
        rgb = numpy.zeros((len(peaks),3)) + rgb[:,numpy.newaxis]
      else:
        raise ValueError("Wrong shape of 'rgb'")
    elif not (rgb.ndim == 2 and rgb.shape == (len(peaks),3)):
        raise ValueError("Wrong shape of 'rgb'")
  
  # Open a cmap file 
  fid = open('%(f)s.cmap' % {'f': filename}, 'w')
  
  # Write the header 
  fid.write('<!DOCTYPE Colormap>\n')
  fid.write('<ColormapVisage2.0 Name="%(f)s">\n' % {'f': filename})
  fid.write('  <Graph Active="1" Type="0" Name="">\n')
  
  # Initialize a counter for the contour values 
  counter = 0
  
  # Initialize a string for the contours 
  c_str = ('    <Control Opacity="%(o)f" Number="%(c)d" Blue="%(b)f"' + 
            ' Red="%(r)f" Green="%(g)f" Value="%(p)f"/>\n' )
  
  # Write the initial value at zero with zero opacity 
  fid.write(c_str % {'o': 0, 'c': counter, 'r': 0, 'g': 0, 'b': 0, 
                     'p': peak_minus})
  counter += 1
  
  for i,p in enumerate(peaks):
    # Write the peak 
    # Left edge of the peak
    fid.write(c_str % {'o': 0, 'c': counter, 'p': p-peak_width/2., 
                       'r': rgb[i,0], 'g': rgb[i,1], 'b': rgb[i,2]})
    counter += 1
    # The peak 
    fid.write(c_str % {'o': 1-alpha, 'c': counter, 'p': p, 
                       'r': rgb[i,0], 'g': rgb[i,1], 'b': rgb[i,2]})
    counter += 1
    # Right edge of the peak
    fid.write(c_str % {'o': 0, 'c': counter, 'p': p+peak_width/2., 
                       'r': rgb[i,0], 'g': rgb[i,1], 'b': rgb[i,2]})
    counter += 1
  
  # Write a final value at 1.5 * (final peak value) with zero opacity 
  fid.write(c_str % {'o': 0, 'c': counter, 'r': 0, 'g': 0, 'b': 0, 
                     'p': peak_plus})
  
  # Finalize the file 
  fid.write('  </Graph>\n')
  fid.write('</ColormapVisage2.0>')
  
  # Close the file
  fid.close()
def meshgrid2(*arrs):
  '''adapted from:
  http://stackoverflow.com/a/1830192
  '''
  arrs = tuple(reversed(arrs)) 
  lens = map(len, arrs)
  dim = len(arrs)
  
  sz = 1
  for s in lens:
      sz*=s
  
  ans = []    
  for i, arr in enumerate(arrs):
    slc = [1]*dim
    slc[i] = lens[i]
    arr2 = numpy.asarray(arr).reshape(slc)
    for j, sz in enumerate(lens):
      if j!=i:
        arr2 = arr2.repeat(sz, axis=j) 
    ans.append(arr2)
  
  return tuple(ans[::-1])

def view_with_mayavi(x,y,z,data,geo_spec=None,datalabels=None,
                     iso_min=1e-4,iso_val=0.01,iso_max=10.0):
  ''' Creates an interactive mayavi dialog showing isosurface plots of the input
  data. 
  
  Components adapted from:
  http://stackoverflow.com/a/1830192
  http://docs.enthought.com/mayavi/mayavi/auto/example_mlab_interactive_dialog.html
  
  **Parameters:**
  
  x,y,z : numpy.ndarray, 1-dim
    Contains the grid.
  data : numpy.ndarray, shape=(len(x), len(y),len(z)) or shape=(N, len(x), len(y),len(z))
    Contains the output data.
  geo_spec : 
    See :ref:`Central Variables` for details.
    If not None, the atom positions will be drawn additionally.
  datalabels : None or list of str
    Contains information about the plotted data with len(datalabels) == len(data).
  '''
  try:
    from enthought.traits.api import HasTraits, Range, Instance, on_trait_change, Bool, Str,List,Button
    from enthought.traitsui.api import View, Item, Group, ListStrEditor,HSplit
  except ImportError:
    from traits.api import HasTraits, Range, Instance, on_trait_change, Bool, Str,List,Button
    from traitsui.api import View, Item, Group, ListStrEditor,HSplit
  
  try:
    from enthought.mayavi import mlab
    from enthought.mayavi.core.api import PipelineBase
    from enthought.mayavi.core.ui.api import MayaviScene, SceneEditor, MlabSceneModel
  except ImportError:
    from mayavi import mlab
    from mayavi.core.api import PipelineBase
    from mayavi.core.ui.api import MayaviScene, SceneEditor, MlabSceneModel
  
  from copy import deepcopy
  data = numpy.array(data)
  
  if data.ndim == 3:
    data = data[numpy.newaxis]
  #else:
    #try:
      #data = data.reshape((-1,) + tuple(grid.N_))
    #except ValueError:
  #if data.ndim != 4:
    #raise ValueError('`data` has to be a ``numpy.array`` with four dimensions')
  if datalabels is not None and len(datalabels) != len(data):
    raise ValueError('`datalabels` has to be a list of strings with the same' +
                     'length as `data`.')
  if datalabels is not None:
    datalabels = ['%03d: %s' % (i,j) for i,j in enumerate(datalabels)]
  
  Z,Y,X = meshgrid2(z,y,x)
  
  class MyModel(HasTraits):  
      select  = Range(0, len(data)-1, 0)
      last_select = deepcopy(select)
      iso_value  = Range(iso_min, iso_max, iso_val,mode='logslider')
      opacity    = Range(0, 1.0, 0.85)
      show_atoms = Bool(True)
      label = Str()
      available = List(Str)
      available = datalabels
      
      prev_button = Button('Previous')
      next_button = Button('Next')
      
      scene = Instance(MlabSceneModel, ())
      
      plot_atoms = Instance(PipelineBase)
      plot0 = Instance(PipelineBase)
      
      # When the scene is activated, or when the parameters are changed, we
      # update the plot.
      @on_trait_change('select,iso_value,show_atoms,opacity,label,scene.activated')
      def update_plot(self):
        #if self.select < len(data)-1:
          #self.select += 1
        
        if self.plot0 is None:          
          src = mlab.pipeline.scalar_field(X,Y,Z,data[self.select])
          self.plot0 = self.scene.mlab.pipeline.iso_surface(\
                      src, contours= [-self.iso_value,self.iso_value], opacity=self.opacity,colormap='blue-red',vmin=-1e-8,vmax=1e-8)
          lut = self.plot0.module_manager.scalar_lut_manager.lut.table.to_array()
          self.plot0.module_manager.scalar_lut_manager.lut.table = lut[::-1]
          self.plot0.contour.scene.background = (1,1,1)
        elif self.select != self.last_select:
          self.plot0.mlab_source.set(scalars=data[self.select])
        self.plot0.contour.contours = [-self.iso_value,self.iso_value]
        self.plot0.actor.property.opacity = self.opacity
        self.last_select = deepcopy(self.select)
        if datalabels is not None:
          self.label = datalabels[self.select]
        if geo_spec is not None: 
          if self.plot_atoms is None:            
            self.plot_atoms = self.scene.mlab.points3d(geo_spec[:,0],geo_spec[:,1],geo_spec[:,2])
          self.plot_atoms.visible = self.show_atoms
      
      def _prev_button_fired(self):
        if self.select > 0:
          self.select -= 1
      def _next_button_fired(self):
        if self.select < len(data)-1:
          self.select += 1
      
      
      # The layout of the dialog created
      items = (Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                      height=400, width=600, 
                      show_label=False),)
      items0 = ()
      if len(data) > 1:
        items0 += (Group('select',
                        HSplit(Item('prev_button', show_label=False),
                               Item('next_button', show_label=False)
                              )),)               
      items0 += (Group('iso_value', 'opacity','show_atoms'
                      ),)
              
      if datalabels is not None:
        if len(datalabels) > 1:        
            items1 = (Item('available', 
                       editor=ListStrEditor(title='Available Data',editable=False),
                       show_label=False,style='readonly',width=300
                       ),)
            items0 = HSplit(items0,items1)
        items += (Group('_',Item('label',label='Selected Data',style='readonly', show_label=True),'_'),
               items0,)
      else:
        items += items0
      view = View(*items,
                  resizable=True
                  )
  
  my_model = MyModel()
  my_model.configure_traits()

def pdb_creator(geo_info,geo_spec,filename='new',charges=None,comments='',
    angstrom=True):
  '''Creates a plain text file with information concerning the molecular 
  structure in Protein Data Bank (PDB) format. 
  
  **Parameters:**
  
  geo_info, geo_spec : 
    See :ref:`Central Variables` for details.
  filename : str
    Contains the base name of the output file.
  charges : numpy.ndarray, shape=(Natoms,), optional
    Contains a partial charge for each atom.
  comments : str, optional
    Specifies the second (comment) line of the pdb file.
  angstrom : bool, optional
    If True, conversion of molecular coordinates from Bohr radii to Angstrom.
  '''
  
  # Conversion factor form Angstrom in Bohr radii
  aa_to_au = 1/0.52917720859 if angstrom else 1.
  
  # Check for charges
  if charges is None:
    charges = numpy.zeros(len(geo_spec))
    
  # Open an empty file 
  fid = open('%(f)s.pdb' % {'f': filename},'w')
  
  # Write HEADER, TITLE and AUTHOR
  fid.write('HEADER\n')
  fid.write('TITLE    %s\n' % comments)
  fid.write('AUTHOR    orbkit\n') 
  
  # Write ATOM records
  string = ''
  for il in range(len(geo_spec)):
    string = 'ATOM    %s' % (il+1)
    while len(string) < 13:
      string += ' '
    string += '%s' % (geo_info[il][0])
    while len(string) < 17:
      string += ' '
    xyz = geo_spec[il]/aa_to_au
    string += '             %s %s %s        ' % (('%.9f' % xyz[0])[:7],
                                                 ('%.9f' % xyz[1])[:7],
                                                 ('%.9f' % xyz[2])[:7])
    string += '%s        ' % (('%.9f' % charges[il])[:6])
    fid.write(str(string) + '\n')
      
  # Write MASTER and END line
  fid.write('MASTER        0    0    0    0    0    0    0    0 ' + 
            '%s    0    0    0\nEND' % ('%s'.rjust(4) % len(geo_spec)))
  
  # Close the file 
  fid.close()

def xyz_creator(geo_info,geo_spec,filename='new',charges=None,comments='',
    angstrom=True):
  '''Creates a xyz file containing the molecular coordinates. 
  
  **Parameters:**
  
  geo_info, geo_spec : 
    See :ref:`Central Variables` for details.
  filename : str
    Contains the base name of the output file.
  charges : numpy.ndarray, shape=(Natoms,), optional
    Contains a partial charge for each atom.
  comments : str, optional
    Specifies the second (comment) line of the xyz file.
  angstrom : bool, optional
    If True, conversion of molecular coordinates from Bohr radii to Angstrom.
  '''
  
  # Conversion factor form Angstrom in Bohr radii
  aa_to_au = 1/0.52917720859 if angstrom else 1.
  
  # Open an empty file 
  fid = open('%(f)s.xyz' % {'f': filename},'w')
  
  # Write number of atoms and a comment line
  fid.write('%d\n%s\n' % (len(geo_spec),comments))
  
  # Write Cartesian coordinates of molecular structure
  string = ''
  for il in range(len(geo_spec)):
    string += '%-2s' % geo_info[il][0]
    for i in range(3):
      string += ' %22.15f'  % (geo_spec[il][i]/aa_to_au)
    if charges is not None:
      string += ' %22.15f'  % charges[il]
    string += '\n'
  fid.write(string)
  
  # Close the file
  fid.close()
