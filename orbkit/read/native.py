'''Input module for Orbkit native format
'''

def import_module(name):
  lclass = len(name.split('.')[-1]) + 1
  module = __import__(name[:-lclass])
  return getattr(module, name[-lclass+1:])

def str_to_bool(entry):
  if isinstance(entry, str):
    for revtype in [True, False, None]:
      if entry == str(revtype):
        entry = revtype
  return entry

# ``all_mo`` and ``spin`` make no sense here but defining them anyways makes ``main_read`` cleaner.
def read_native(fname, all_mo=None, spin=None):

  ftype = fname.split('.')[-1]

  if ftype.lower() in ['numpy', 'npz']:
    from numpy import load
    data = load(fname)
  elif ftype.lower() in ['hdf5', 'h5']:
    import h5py
    fd = h5py.File(fname, 'r')
    data = {}
    for key in fd.keys():
      entry = fd[key][...]
      if isinstance(entry, str):
        entry = str_to_bool(entry)
      data[key] = entry
    for key in fd.attrs.keys():
      entry = fd.attrs[key]
      if isinstance(entry, str):
        entry = str_to_bool(entry)
      data[key] = entry
  else:
    raise NotImplementedError('File format {0} not implemented for reading.'.format(ftype.lower()))

  parent_module_instance = import_module(str(data['parent_class_name']))
  return parent_module_instance(data)
