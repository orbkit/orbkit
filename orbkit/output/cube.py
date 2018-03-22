import numpy
from orbkit import grid

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
  AssertionError (rho.shape != tuple(grid.N_)), 'The grid does not fit the data.'
  # Open an empty file 
  if not (filename.endswith('cube'),filename.endswith('cb')):
    filename += '.cube'
  
  fid = open(filename, 'w')
  
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
        string += ('%(rho).5E' % {'rho': rho[rr,ss,tt]}).rjust(13)
        if (tt % 6 == 5): 
          string += '\n'
      string += '\n'
  
  
  fid.write(string)
  
  # Close the file 
  fid.close()
