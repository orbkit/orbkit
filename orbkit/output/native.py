'''Output module for Orbkit native format
'''
import time
from os.path import join,splitext

def unravel_dicts(indict):
  '''Unravels encapsulated dictionaries stemming from class-subclass structures'''
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
  
def write_native(outdata, outputname, ftype='numpy', gname='qcinfo', mode='w',
                 **add_data):
  '''Creates the requested output.
  
  **Parameters:**

    outdata : 
      Instance of the class to be written to fine. This must support a **todict** function.
      outputname : str or list of str
      Contains the base name of the output file.
    ftype, str : 
      Data can be written to ``.npz`` or ``hdf5`` files.
  '''

  if not callable(getattr(outdata, 'todict')):
    raise NotImplementedError('Supplied class instance does not support todict()')

  if ftype.lower() in ['hdf5', 'h5']:
    from . import hdf5_write
    write = hdf5_write
    if not (outputname.endswith('hdf5') or outputname.endswith('h5')):
      outputname += '.hdf5'
  elif ftype.lower() in ['numpy', 'npz']:
    from . import npz_write
    write = npz_write
    if not (outputname.endswith('numpy') or outputname.endswith('npz')):
      outputname += '.npz'
  else:
    raise NotImplementedError('File format {0} not implemented for writing.'.format(ftype.lower()))
    #raise NotImplementedError('Only .npz and .hdf5 are currently supported.')
  
  odata = outdata.todict()
  odata.update(add_data)
  odata['date'] = time.strftime("%Y-%m-%d") 
  odata['time'] = time.strftime("%H:%M:%S")
  odata = unravel_dicts(odata)
  for i,key in enumerate(odata.keys()):
    if i == 0:
      current_mode = mode
    else:
      current_mode = 'a'
    write(outputname, mode=current_mode, gname=join(gname,key), **odata[key])
  
  if gname != '':
    return '{0}@{1}'.format(outputname,gname)
  else:
    return outputname





