import numpy

from orbkit.units import a0_to_aa

def xyz_creator(geo_info,geo_spec,filename='new',mode='w',charges=None,comments='',
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
  
  # Open an empty file 
  fid = open('%s.xyz' % filename,mode)
  
  # Write number of atoms and a comment line
  fid.write('%d\n%s\n' % (len(geo_spec),comments))
  
  # Write Cartesian coordinates of molecular structure
  string = ''
  for il in range(len(geo_spec)):
    string += '%-2s' % geo_info[il][0]
    for i in range(3):
      xyz = geo_spec[il]
      if angstrom:
        xyz *= a0_to_aa
      string += ' %22.15f'  % (xyz[i])
    if charges is not None:
      string += ' %22.15f'  % charges[il]
    string += '\n'
  fid.write(string)
  
  # Close the file
  fid.close()
