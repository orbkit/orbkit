import gzip
import numpy
from orbkit import grid
from skimage import measure

def obj_creator(data,filename,geo_info,geo_spec,iso=(-0.01,0.01),**kwargs):
  '''Creates a plain text 3D surface file (obj). 
  
  **Parameters:**
  
  data : numpy.ndarray, shape=N
    Contains the output data.
  filename : str
    Contains the base name of the output file.
  geo_info, geo_spec : 
    See :ref:`Central Variables` for details.
  iso : tuple
    Specifies the isovalue for negative and posititve surfaces, respectively.
  '''
  # Shape shall be (Ndrv,Ndata,Nx,Ny,Nz) or (Ndrv,Ndata,Nxyz)
  data = numpy.array(data)
  dims = 3
  ndim = data.ndim
  shape = data.shape
    
  if data.ndim < dims:
    raise AssertionError('data.ndim < ndim of grid')
  if data.ndim > dims + 1:
    raise AssertionError('data.ndim > (ndim of grid) +2')

  if shape != tuple(grid.N_):
    raise AssertionError('The grid does not fit the data.')

  # Conversion of 3d array and storage as surface file
  basisname = filename.replace('.obj.gz', '').replace('.obj', '')
  spacing = tuple(numpy.concatenate(grid.delta_))
  for i in iso:
        
    # Data conversion    
    verts, faces, normals, values = measure.marching_cubes(data,
                                                           level=i,
                                                           spacing=spacing,
                                                           step_size=1)
    faces = faces + 1
    
    # Define filename for plus and minus surfaces
    pm = 'plus' if i > 0 else 'minus'
    fid = f'{basisname}_{pm}.obj'

    # Write vertices
    string = '# Number of vertices: %d\n' % len(verts)
    for item in verts:
      string += 'v %0.6g %0.6g %0.6g\n' % (item[0],item[1],item[2])
          
    # Write normals
    string += '# Number of normals: %d\n' % len(normals)
    for item in normals:
      string += ('vn %.3g' % item[0]).replace('0.', '.')
      string += (' %0.3g' % item[1]).replace('0.', '.')
      string += (' %0.3g\n' % item[2]).replace('0.', '.')
      
    # Write faces
    string += '# Number of faces: %d\n' % len(faces)
    for item in faces:
      string += 'f {0}//{0} {1}//{1} {2}//{2}\n'.format(item[0],item[1],item[2])

    if filename.endswith('gz'):
      with gzip.open(fid + '.gz', 'wb') as f:
        f.write(string.encode('utf-8'))
    else:
      with open(fid, 'w') as f:
        f.write(string)
