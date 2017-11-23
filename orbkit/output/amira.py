import numpy
from orbkit import grid

def amira_creator(data,filename):
  '''Creates a ZIBAmira mesh file. (plain text)
  
  **Parameters:**
  
  data : numpy.ndarray, shape=N or (3,)+N
    Contains the output data.
  filename : str
    Contains the base name of the output file.
  '''
  if data.ndim == 3:
    typename = 'double'
    formatter = '{0:.15e}\n'
    data = numpy.swapaxes(data,0,2).reshape((-1,1)) #ravel(order='F')[:,numpy.newaxis]
  elif data.ndim == 4:
    typename = 'double[3]'
    data = numpy.array([j.ravel(order='F') for j in data])
    formatter = '{0:.15e} {1:.15e} {2:.15e}\n'
  else:
    raise IOError("amira_creator only supports 3D or 4D data.")
  
  N = tuple(grid.N_)
  
  # Open an empty file 
  fid = open('%s.am' % filename,'w')
  # Write Header 
  fid.write('# AmiraMesh 3D ASCII 2.0\n\n\n')
  fid.write('define Lattice %d %d %d\n' % N)
  fid.write('Parameters {\n')
  fid.write('    Content "%dx%dx%d %s, uniform coordinates",\n' % 
            (N + (typename,)))
  fid.write('    BoundingBox %(xmin)f %(xmax)f %(ymin)f %(ymax)f %(zmin)f %(zmax)f,\n' %
            {'xmin': grid.min_[0],'xmax': grid.max_[0],
             'ymin': grid.min_[1],'ymax': grid.max_[1],
             'zmin': grid.min_[2],'zmax': grid.max_[2]})
  fid.write('    CoordType "uniform"\n}\n\n')
  fid.write('Lattice { %s Data } @1\n'%typename)
  fid.write('# Data section follows\n@1\n')
  for i in data:
    fid.write(formatter.format(*i))
  #for tt in range(len(grid.z)):
    #for ss in range(len(grid.y)):
      #for rr in range(len(grid.x)): 
        #string += '%g\n' % rho[rr,ss,tt]
  fid.close()

def amira_creator_old(rho,filename):
  '''Creates a ZIBAmira mesh file. (plain text)
  
  **Parameters:**
  
  rho : numpy.ndarray, shape=N
    Contains the output data.
  filename : str
    Contains the base name of the output file.
  '''
  # Open an empty file 
  fid = open('%(f)s.am' % {'f': filename},'w')

  # usage:
  #     - open Amira
  #     - left-click File -> Open Data
  #     - choose the sfilename.am
  #     - Press OK
  #     - right-click on the new Data Object -> Compute -> Arithmetic
  #     - choose the Arithmetic object
  #     - select  Expr -> A and Result type -> regular and Apply
  #     - use the new data object to display your density as usual,
  
  # Write Header 
  fid.write('# AmiraMesh 3D ASCII 2.0\n\n\n')
  fid.write('define Lattice %(Nx)d %(Ny)d %(Nz)d\n' % 
                    {'Nx': grid.N_[0],'Ny': grid.N_[1],'Nz': grid.N_[2]})
  fid.write('define Coordinates %(N)d\n\n' % {'N': numpy.sum(grid.N_)})
  fid.write('Parameters {\n')
  fid.write('    Content "%(Nx)dx%(Ny)dx%(Nz)d float, uniform coordinates",\n' %
                    {'Nx': grid.N_[0],'Ny': grid.N_[1],'Nz': grid.N_[2]})
  fid.write('    BoundingBox %(xmin)f %(xmax)f %(ymin)f %(ymax)f %(zmin)f %(zmax)f,\n' %
            {'xmin': grid.min_[0],'xmax': grid.max_[0],
             'ymin': grid.min_[1],'ymax': grid.max_[1],
             'zmin': grid.min_[2],'zmax': grid.max_[2]})
  fid.write('    CoordType "uniform"\n}\n\n')
  fid.write('Lattice { float Data } @1\n')
  fid.write('# Data section follows\n@1\n')
  
  # Write density information to .am file
  string = ''
  for tt in range(len(grid.z)):
    for ss in range(len(grid.y)):
      for rr in range(len(grid.x)): 
        string += '%g\n' % rho[rr,ss,tt]
  
  fid.write(string)
  
  # Close the file 
  fid.close()
