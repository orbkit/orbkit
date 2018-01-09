import numpy

from orbkit.units import a0_to_aa

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
    xyz = geo_spec[il]
    if angstrom:
      xyz *= a0_to_aa
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
