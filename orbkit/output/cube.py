import gzip
import numpy
from orbkit import grid

def cube_creator(data,filename,geo_info,geo_spec,comments='',labels=None,**kwargs):
  '''Creates a plain text Gaussian cube file. 
  
  **Parameters:**
  
  data : numpy.ndarray, shape=N
    Contains the output data.
  filename : str
    Contains the base name of the output file.
  geo_info, geo_spec : 
    See :ref:`Central Variables` for details.
  comments : str, optional
    Specifies the second (comment) line of the cube file.  
  '''
  # Shape shall be (Ndrv,Ndata,Nx,Ny,Nz) or (Ndrv,Ndata,Nxyz)
  data = numpy.array(data)
  dims = 3
  ndim = data.ndim
  shape = data.shape
    
  if data.ndim < dims:
    raise AssertionError('data.ndim < ndim of grid')
  elif data.ndim == dims: # 3d data set
    data = data[numpy.newaxis]
  elif data.ndim > dims + 1:
    raise AssertionError('data.ndim > (ndim of grid) +2')
  
  if labels is not None:
    if labels == True or labels == 'auto':
      labels = list(range(len(data)))
    assert len(labels) == len(data)
    try: 
      labels = map(int,labels)
    except ValueError:
      raise AssertionError('labels has to be list of integers.')

  assert data.shape[1:] == tuple(grid.N_), 'The grid does not fit the data.'
  
  if not any([filename.endswith(ext)
              for ext in ['cube', 'cb', 'cube.gz', 'cb.gz']]):
    filename += '.cube'
  
  # Write the type and the position of the atoms in the header
  string = 'orbkit calculation\n'
  string += ' %(f)s\n'  % {'f': comments}
  # How many atoms 
  string += ('%(at)d' % {'at': (-1)**(labels is not None)*len(geo_info)}).rjust(5)
  # Minima
  for ii in range(3):
    string += ('%(min)0.6f' % {'min': grid.min_[ii]}).rjust(12)
  # Number of data points per grid point
  if len(data) > 1:
    string += ('%d' % len(data)).rjust(12)
  for ii in range(3):
    string += '\n'
    string += ('%(N)d' % {'N': grid.N_[ii]}).rjust(5)
    for jj in range(3):
      if jj == ii: 
        string += ('%(dr)0.6f' % {'dr': grid.delta_[ii]}).rjust(12)
      else:
        string += ('%(dr)0.6f' % {'dr': 0}).rjust(12)
  
  string += '\n'
  for ii in range(len(geo_info)):
    string += ('%(N)d' % {'N': round(float(geo_info[ii][2]))}).rjust(5)
    string += ('%(ch)0.6f' % {'ch': float(geo_info[ii][1])}).rjust(12)
    for jj in range(3):
      string += ('%(r)0.6f' % {'r': geo_spec[ii][jj]}).rjust(12)
    string += '\n'
  # If exist, write labels
  if labels is not None:
    string += ('%(N)d' % {'N': len(data)}).rjust(5)
    c = 0
    for i,j in enumerate(labels):
      c += 1
      string += str(j).rjust(5)
      if (c  % 9 == 8):
        string += '\n'
    string += '\n'
  # Write data
  for rr in range(len(grid.x)):
    for ss in range(len(grid.y)):
      c = 0
      for tt in range(len(grid.z)):
        for dd in data[:,rr,ss,tt]:
          string += ('%(data).5E' % {'data': dd}).rjust(13)
          if (c % 6 == 5): 
            string += '\n'
          c += 1
      string += '\n'
  
  if filename.endswith('gz'):
    with gzip.open(filename, 'wb') as f:
      f.write(string.encode('utf-8'))
  else:
    with open(filename, 'w') as f:
      f.write(string)
