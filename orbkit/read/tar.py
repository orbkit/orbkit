'''
Tar interface for ORBKIT readers
'''

import tarfile
import numpy

from .tools import find_itype

def is_tar_file(infile):
  itype = ''
  if '.' in infile:
      if 'tar' not in infile.split('.')[-2:]:
        return None
      else:
        return True
  else:
    return None

def get_all_files_from_tar(infile, sort=False):
  files = []
  itypes = []

  fd = tarfile.open(infile, 'r')
  for tarinfo in fd:
    if tarinfo.isreg():
      one_file = fd.extractfile(tarinfo)
      itype = find_itype(one_file)
      one_file = fd.extractfile(tarinfo)
      itypes.append(itype)
      files.append(one_file)
  files = numpy.array(files)
  itypes = numpy.array(itypes)
  if sort:
    filenames = numpy.array([f.name for f in files], dtype=str)
    si = numpy.argsort(filenames)
    files = files[si]
    itypes = itypes[si]
    del filenames, si
  return files, itypes

#The ci_descriptor attribute is a hack which becomes
#outdated as soon as find_itype recognizes CI Files
#If this is done the main Ci reader should probably also
#automatically recognize filetypes in the same way the
#"normal" high_level reader does now.
def get_file_from_tar(infile, index=0, ci_descriptor=False):
  one_file = None
  i = 0
  fd = tarfile.open(infile, 'r')
  for tarinfo in fd:
    if tarinfo.isreg():
      if i == index:
        one_file = fd.extractfile(tarinfo)
        itype = None
        if not ci_descriptor:
          itype = find_itype(one_file)
        one_file = fd.extractfile(tarinfo)
        break
      i += 1
  if not one_file:
    raise ValueError('File {} not found'.format(i))
  return one_file, itype


