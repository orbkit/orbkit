'''Input module for Orbkit native format
'''
import numpy

def import_module(name):
  lclass = len(name.split('.')[-1]) + 1
  module = __import__(name[:-lclass])
  return getattr(module, name[-lclass+1:])

def str_to_bool(entry):
  if isinstance(entry, (str,numpy.bytes_)):
    for revtype in [True, False, None]:
      if str(entry) == str(revtype):
        entry = revtype
  return entry

def parse_group(group):
  group_data = {}
  for key in group.keys():
    entry = group[key][...]
    if entry.dtype in (numpy.bytes_, '|S10'):
      entry = numpy.array(entry, dtype=str)
      if len(entry.shape) == 1:
        if all([e in [True, False, None] for e in entry]):
          entry = numpy.array(entry, dtype=bool)
    group_data[key] = entry
  for key in group.attrs.keys():
    entry = group.attrs[key]
    if isinstance(entry, (str,numpy.bytes_)):
      entry = str_to_bool(entry)
    group_data[key] = entry
  return group_data

# ``all_mo`` and ``spin`` make no sense here but defining them anyways makes ``main_read`` cleaner.
def read_native(fname, all_mo=None, spin=None, **kwargs):

  ftype = fname.split('.')[-1]

  if ftype.lower() in ['numpy', 'npz']:
    from numpy import load
    data = load(fname)
  elif ftype.lower() in ['hdf5', 'h5']:
    import h5py
    fd = h5py.File(fname, 'r')
    data = {}
    for name in fd.keys():
      if isinstance(fd[name], h5py._hl.group.Group):
        data[name] = parse_group(fd[name])
      else:
        entry = fd[name][...]
        if entry.dtype in (numpy.bytes_, '|S10'):
          entry = numpy.array(entry, dtype=str)
          if len(entry.shape) == 1:
            if all([e in [True, False, None] for e in entry]):
              entry = numpy.array(entry, dtype=bool)
        data[name] = entry
    for name in fd.attrs.keys():
      entry = fd.attrs[name]
      if isinstance(entry, (str,numpy.bytes_)):
        entry = str_to_bool(entry)
      data[name] = entry
  else:
    raise NotImplementedError('File format {0} not implemented for reading.'.format(ftype.lower()))

  parent_module_instance = import_module(str(data['parent_class_name']))
  return parent_module_instance(data)






