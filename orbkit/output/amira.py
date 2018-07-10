import os
import numpy

from orbkit import grid
from orbkit.display import display

from .tools import colormap_creator

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
  filename += '.am' if not filename.endswith('.am') else ''
  fid = open(filename,'w')
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

def hx_network_creator(rho,filename):
  '''Creates a ZIBAmira hx-network file including a colormap file (.cmap)
  adjusted to the density for the easy depiction of the density.
  '''
  from .hx_network_draft import hx_network
  filename += '.hx' if not filename.endswith('.hx') else ''
    
  # Create a .cmap colormap file using the default values 
  display('\tCreating ZIBAmira colormap file...\n\t\t' + filename.replace('.hx','.cmap'))
  
  assert (rho.shape != tuple(grid.N_)), 'The grid does not fit the data.'
  
  colormap_creator(rho,filename.replace('.hx','.cmap'))
  
  # Create a .hx network file based on the file orbkit.hx_network_draft.py 
  # Open an empty file
  
  fid = open(filename,'w')
  
  # Copy the content of the draft file and replace the keywords 
  fid.write(hx_network.replace("FILENAME",os.path.splitext(os.path.basename(filename))[0])) 
  
  # Close the file 
  fid.close()  
