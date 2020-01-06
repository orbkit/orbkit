import os
import numpy

from orbkit import options, grid
from orbkit.display import display
from orbkit.orbitals import AOClass, MOClass

from .native import write_native

def hdf5_creator(data, filename, qcinfo=None, gname='', ftype='hdf5', mode='w',attrs={},**kwargs):
  '''Creates an hdf5 file (Hierarchical Data Format) output.
  
  **Parameters:**
  
  data : numpy.ndarray
    Contains the output data.
  filename : str
    Contains the base name of the output file.
  qcinfo : 
    Stores all information from qcinfo class (see :ref:`Central Variables` for
    details).
  gname : str, optional
    Specifies the path, where the data is stored.
  mode : str, optional
    Specifies the mode used to open the file. ('r', 'w', 'a')
  attrs : dict or dict of dicts, optional
    Supply datasets with attributes. If dict, entries will be added to 'data'.
  grid : module or class, global
    Contains the grid, i.e., grid.x, grid.y, and grid.z.
  '''
  
  if ftype.lower() in ['hdf5', 'h5']:
    write = hdf5_write
    attributes = hdf5_attributes
  elif ftype.lower() in ['numpy', 'npz']:
    write = npz_write
    attributes = npz_attributes
  else:
    raise NotImplementedError('File format {0} not implemented for writing.'.format(ftype.lower()))
  
  # Save data and all other kwargs
  write(filename,mode=mode,gname=gname,data=data,**kwargs)
  
  # Save grid
  write(filename,mode='a',gname=os.path.join(gname,'grid'),
             x=grid.x,y=grid.y,z=grid.z,
             is_vector=grid.is_vector,is_regular=grid.is_regular)
  
  # Set attributes
  attributes(filename,gname=gname,**attrs)
  
  if qcinfo is not None:
    write_native(qcinfo, outputname=filename, ftype=ftype, mode='a', 
                 gname=os.path.join(gname,'qcinfo'))  

def npz_write(filename,gname='',mode='w',compress=True,**namedict):
  '''Function borrowed from numpy.lib.npyio._savez to allow for groups and different modes.
  '''
  import zipfile

  if not filename.endswith('npz'):
    filename += '.npz'
  
  if compress:
    compression = zipfile.ZIP_DEFLATED
  else:
    compression = zipfile.ZIP_STORED
  
  zipf = numpy.lib.npyio.zipfile_factory(filename, mode=mode, compression=compression)

  # Stage arrays in a temporary file on disk, before writing to zip.

  # Import deferred for startup time improvement
  import tempfile
  # Since target file might be big enough to exceed capacity of a global
  # temporary directory, create temp file side-by-side with the target file.
  file_dir, file_prefix = os.path.split(filename) if  numpy.lib.npyio._is_string_like(filename) else (None, 'tmp')
  fd, tmpfile = tempfile.mkstemp(prefix=file_prefix, dir=file_dir, suffix='-numpy.npy')
  os.close(fd)
  try:
    for key, val in namedict.items():
      fname = key + '.npy'
      fid = open(tmpfile, 'wb')
      try:
        if val is None:
            continue
        numpy.lib.format.write_array(fid, numpy.asanyarray(val), allow_pickle=False)
        fid.close()
        fid = None
        zipf.write(tmpfile, arcname=os.path.join(gname,fname))
      except IOError as exc:
        raise IOError("Failed to write to %s: %s" % (tmpfile, exc))
      finally:
        if fid:
          fid.close()
  finally:
    os.remove(tmpfile)

  zipf.close()

def npz_attributes(filename, gname='', **attrs):
  npz_write(filename,gname=gname+'_attrs',mode='a',**attrs)

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
  if isinstance(fid,h5py.File):
    yield fid
  
  if not (fid.endswith('h5') or fid.endswith('hdf5')):
    fid += '.h5'
  
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
    if gname != '':
      group = f.create_group(gname) if gname not in f.keys() else f[gname]
    else:
      group = f
    for key,data in kwargs.items():
      if isinstance(data,(list,numpy.ndarray)):
        data = numpy.array(data)
        if data.dtype.type is numpy.unicode_:
          data = numpy.asarray(data,dtype=numpy.string_)
        elif data.dtype in ['O']:
          data = numpy.array(data, dtype='S10')
        group[key] = data
      elif isinstance(data,dict):
        hdf5_append(data,group,name=key)
      else:
        if data is None or type(data) == bool:
          data = str(data)
        group.attrs[key] = data

def hdf5_attributes(filename, gname='', **attrs):
  warning = []
  for f in hdf5_open(filename,mode='a'):
    for key,value in attrs.items():
      if isinstance(value,dict):
        for k,v in value.items():
          if k in f[gname].keys():
            f[os.path.join(gname,key)].attrs[k] = v
          else:
            warning.append(key)
      else:
        if key in f[gname].keys():
          f[os.path.join(gname,'data')].attrs[key] = value
        else:
          warning.append('data')
  if warning:
    display('Warning: The following attributes are not datasets in '+filename)
    display(', '.join([os.path.join(gname,i) for i in warning]))
