import numpy

from orbkit import options, grid
from orbkit.display import display
from .tools import colormap_creator

def hdf5_creator(data,filename,geo_info,geo_spec,mode='w',data_id='rho',group=None,
            data_only=False,ao_spec=None,mo_spec=None,is_mo_output=False,
            x=None,y=None,z=None,**kwargs):
  '''Creates an hdf5 file (Hierarchical Data Format) output.
  
  **Parameters:**
  
  data : numpy.ndarray
    Contains the output data.
  filename : str
    Contains the base name of the output file.
  geo_info, geo_spec : 
    See :ref:`Central Variables` for details.
  data_id : str, optional
    Specifies name of the dataset in the hdf5 file.
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
  hdf5_file = h5py.File('%s.h5' % filename, mode)
  if group is None:
    f = hdf5_file
  else:
    hdf5_file.require_group(group)
    f = hdf5_file[group]
    
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
    hdf5_file.close()
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
  
  hdf5_file.close()

def hx_network_creator(rho,filename):
  '''Creates a ZIBAmira hx-network file including a colormap file (.cmap)
  adjusted to the density for the easy depiction of the density.
  '''
  from orbkit.hx_network_draft import hx_network
  # Create a .cmap colormap file using the default values 
  display('\tCreating ZIBAmira colormap file...\n\t\t%(f)s.cmap' % 
                {'f': filename})
  
  AssertionError (rho.shape != tuple(grid.N_)), 'The grid does not fit the data.'
  
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

def hdf5_open(fid,mode='w'):
  '''Open an hdf5 file.
  
  **Parameters:**
  
  fid : str
    Contains the filename of the hdf5 file.
  mode : str, optional
    Specifies the mode used to open the file. ('r', 'w', 'a')
  
  **Usage:**
  
  for hdf5_file in hdf5_open(fid,mode='w'):
    # Do something with hdf5_file...
  '''
  import h5py
  hdf5_file = h5py.File(fid, mode)
  try:
    yield hdf5_file
  finally:
    hdf5_file.close()

def hdf5_append(x,group,name='data'):
  '''Automatically append data to an open hdf5 file/group.
  
  **Parameters:**
  
  group : h5py.File or h5py.Group
    The hdf5 file/group where the data will be appended.
  name : string, optional
    Specifies the group/dataset name in the hdf5 file/group. If empty,
    root directory of hdf5 file/group is chosen.  
  '''
  from orbkit.orbitals import MOClass, AOClass

  if isinstance(x, MOClass) or isinstance(x, AOClass):
    x = x.todict()
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
    if x is None:
      x = str('None')
    group.attrs[name] = x
  

def hdf52dict(group,hdf5_file):
  '''Automatically convert the data stored in an hdf5 file/group 
  to a python dictionary.
    
  **Parameters:**
  
  group : str
    Specifies the group/dataset name in the hdf5 file/group.
  hdf5_file : h5py.File or h5py.Group
    The source hdf5 file/group.
  '''
  try:
    # The selected group is a dataset 
    x = hdf5_file['%s' % group][()]
  except:
    # Read all members 
    members = list(hdf5_file['%(g)s' % {'g':group}])
    try:
      members = numpy.array(members, dtype = numpy.int64)
      x = []
      for mm in numpy.sort(members):
        x.append(hdf52dict('%(g)s/%(m)s' % {'g':group, 'm': mm},hdf5_file))
    except (ValueError,OverflowError):
      x = {}
      for mm in members:
        x[mm] = hdf52dict('%(g)s/%(m)s' % {'g':group, 'm': mm},hdf5_file)
      attrs = hdf5_file['%(g)s' % {'g':group}].attrs
      for ii_a in attrs:
        x[ii_a] = attrs[ii_a]
  return x

def hdf5_write(fid,mode='w',gname='',**kwargs):
  '''Writes all keyword arguments to an hdf5 file.
  
  **Parameters:**
  
  fid : str
    Contains the filename of the hdf5 file.
  mode : str, optional
    Specifies the mode used to open the file. ('r', 'w', 'a')
  gname : str
    Specifies the group name in the hdf5 file, where the kwargs will be appended
    to.
  
  **Note:**
  
    All additional keyword arguments specify the input data.
  '''
  for f in hdf5_open(fid,mode=mode):
    if gname is not '':
      group = f.create_group(gname)
    else:
      group = f
    for key,data in kwargs.items():
      if isinstance(data,(list,numpy.ndarray)):
        data = numpy.array(data)
        if data.dtype in ['<U1', '<U4', 'O']:
          data = numpy.array(data, dtype='S10')
        group.create_dataset(key,numpy.shape(data),data=data)
      else:
        if data is None or type(data) == bool:
          data = str(data)
        group.attrs[key] = data







