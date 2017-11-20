'''Output module for Orbkit native format
'''
import time

def recursive_dict_resolution(indict):
  resolved = False
  dictchain = [indict]
  keychain = ['']
  i = 0
  while not resolved or i > 0:
    resolved = True
    keypop = []
    for key in dictchain[i].keys():
      if isinstance(dictchain[i][key], dict):
        dictchain.append(dictchain[i][key])
        keychain.append(key)
        keypop.append(key)
        resolved = False
    for key in keypop:
      dictchain[i].pop(key)
    if not resolved:
      i += 1
    else:
      i -= 1
  outdict = {}
  for i,j in zip(keychain, dictchain):
    outdict[i] = j
  return outdict
  
def write_native(data, outputname='new', ftype='numpy'):
  '''Creates the requested output.
  
  **Parameters:**

    data : 
      Instance of the class to be written to fine. This must support a **todict** function.
      outputname : str or list of str
      Contains the base name of the output file.
    ftype, str : 
      Data can be written to ``.npz`` or ``hdf5`` files.
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
    from orbkit.output import hdf5_write
    odata = recursive_dict_resolution(odata)
    for i,key in enumerate(odata.keys()):
      if i == 0:
        mode = 'w'
      else:
        mode = 'a'
      hdf5_write(outputname + '.' + ftype.lower(), mode, gname=key, **odata[key])
  else:
    raise NotImplementedError('File format {0} not implemented for writing.'.format(ftype.lower()))






