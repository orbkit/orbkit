from orbkit.display import display

def vmd_network_creator(filename,cube_files=None,render=False,exit=False,
                        iso=(-0.01,0.01),abspath=False,mo_options='',**kwargs):
  '''Creates a VMD script file from a list of cube files provided.

  **Parameters:**

  filename : str
    Contains the base name of the output file.
  cube_files : None or list of str
    Specifies the cube files which serve as input for the VMD script.
    If None, searches the directory for '.cb' and '.cube' files.
  render : bool
    If True, the VMD script will automatically create '.tga' files for each
    cube file.
  exit : bool
    If True, exit vmd after rendering.
  iso : tuple
    Specifies the isovalue for the blue and the red isosurface, respectively.
  abspath : bool
    If True, the paths of the cube files will be expanded to absolute file paths.
  vmd_options : strings
    Further commands to include in vmd script.
  '''
  from os import path,listdir
  import linecache
  from orbkit import vmd_network_draft
  
  filename += '.vmd' if not filename.endswith('vmd') else ''
  
  if cube_files is None:
    display('No list of cube (.cb or .cube) filenames provided. Checking the directory' +
            ' of the outputfile...')
    cube_files = []
    for fid in listdir(path.dirname(filename)):
      if fid.endswith('.cb') or fid.endswith('.cube'):
        cube_files.append(fid)
    if cube_files == []:
      raise IOError('Could not find valid cube files in %s' % path.dirname(filename))
  elif isinstance(cube_files,str):
    cube_files = [cube_files]
  elif not isinstance(cube_files,list):
    raise IOError('`cube_files` has to be a list of strings.')

  title = []
  mo = ''
  for i,f in enumerate(cube_files):
    title = linecache.getline(f,2)
    if title.split() == []:
      title = path.splitext(path.basename(f))[0]
    else:
      title = title.replace('\n','').replace(' ','')
    linecache.clearcache()
    pid = path.abspath(f) if abspath else path.relpath(f,path.dirname(filename))
    mo += vmd_network_draft.mo_string % {
                                  'c': i,
                                  'n1': pid,
                                  'n2': title,
                                  'isored': iso[0],
                                  'isoblue': iso[1],
                                  'render': '' if render else '#',
                                  'mo_options' : mo_options,
                                  }
    
  with open(filename) as f:
    f.write(vmd_network_draft.vmd_string % {'mo':mo})
    if exit:
      f.write('exit\n')
  f.close()
