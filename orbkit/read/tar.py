'''
Tar interface for ORBKIT readers
'''

import tarfile

from orbkit.read.high_level import find_itype

def is_tar_file(infile):
  itype = ''
  if '.' in infile:
      if 'tar' not in infile.split('.')[-2:]:
        return None
      else:
        return True
  else:
    return None

def get_all_files_from_tar(infile):
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
  return files, itypes

def get_file_from_tar(infile, index=0):
  i = 0
  fd = tarfile.open(infile, 'r')
  for tarinfo in fd:
    if tarinfo.isreg() and i == index:
      one_file = fd.extractfile(tarinfo)
      itype = find_itype(one_file)
      one_file = fd.extractfile(tarinfo)
      break
  return one_file, itype
