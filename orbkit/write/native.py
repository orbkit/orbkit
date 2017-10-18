'''Output module for Orbkit native format
'''
import time

def write_native(data, outputname='new', ftype='numpy', group=''):
  '''Creates the requested output.
  
  **Parameters:**

    data : 
      Instance of the class to be written to fine. This must support a **todict** function.
      outputname : str or list of str
      Contains the base name of the output file.
    ftype, str : 
      Data can be written to ``.npz`` or ``hdf5`` files.
    group : str, optional
      Name for ``hdf5`` group
  '''

  if not callable(getattr(data, 'todict')):
    raise NotImplementedError('Supplied class instance does not support todict()')

  if ftype.lower() not in ['numpy', 'npz', 'hdf5', 'h5']:
    raise NotImplementedError('Ony npz and hdf5 are currently supportet.')

  odata = data.todict()
  odata['date'] = time.strftime("%Y-%m-%d") 
  odata['time'] = time.strftime("%H:%M:%S")

  if ftype.lower() in ['numpy', 'npz']:
    from numpy import savez_compressed as save
    save(outputname + '.npz', **odata)
  elif ftype.lower() in ['hdf5', 'h5']:
    from orbkit.write import hdf5_write
    hdf5_write(outputname + '.' + ftype.lower(), 'w', group, **odata)
  else:
    raise NotImplementedError('File format {0} not implemented for writing.'.format(ftype.lower()))
