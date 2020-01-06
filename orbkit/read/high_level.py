'''
High level interface to Orbkit reading functions
'''

import re

from .native import read_native
from orbkit.read.molden import read_molden
from orbkit.read.gamess import read_gamess
from orbkit.read.gaussian_fchk import read_gaussian_fchk
from orbkit.read.gaussian_log import read_gaussian_log
from orbkit.read.aomix import read_aomix
from orbkit.read.wfx import read_wfx
from orbkit.read.wfn import read_wfn
from orbkit.read.native import read_native
from orbkit.read.cclib_parser import read_with_cclib
from orbkit.display import display
from orbkit.analytical_integrals import check_mo_norm

from .tools import find_itype, descriptor_from_file

readers = {'molden': read_molden,
           'aomix': read_aomix,
           'gamess': read_gamess, 
           'gaussian.log': read_gaussian_log, 'gaussian_log': read_gaussian_log, 
           'gaussian.fchk': read_gaussian_fchk, 'fchk': read_gaussian_fchk,
           'wfn': read_wfn,
           'wfx': read_wfx,
           'cclib': read_with_cclib,
           'native': read_native
          }                        #: Specifies possible input types.

def main_read(fname, all_mo=False, spin=None, itype='auto', check_norm=False, **kwargs):
  '''
  This is the high-lever interface for the
  orbkit reading routines.
  
  **Parameters:**
  
  fname: str, file descriptor
    Specifies the filename for the input file.
    fname can also be used with a file descriptor instad of a filename.
  all_mo : bool, optional
    If True, all molecular orbitals are returned.
  spin : {None, 'alpha', or 'beta'}, optional
    If not None, returns exclusively 'alpha' or 'beta' molecular orbitals.
  itype : str, optional
    Can be used to manually specify the input filetype.
  check_norm : bool, optional
    If True, ORBKIT verifies that molecular orbitals are orthonormal.
    
  **Note:**
  
    All additional keyword arguments are forwarded to the reading functions.
    
  **Returns:**
  
    qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, mo_spec, etot :
    See :ref:`Central Variables` for details.
  '''

  if isinstance(fname, str):
    filename = fname
  else:
    filename = fname.name

  if itype == 'auto':
    itype = find_itype(fname)
 
  display('Loading data from {0} type file {1}\n'.format(itype, filename))

  qc = readers[itype](fname, all_mo=all_mo, spin=spin, **kwargs)

  if check_norm:
    deviation = check_mo_norm(qc)
    if deviation >= 1e-5:
      raise ValueError('Bad molecular orbital norm: {0:%4e}'.format(deviation))

  return qc
