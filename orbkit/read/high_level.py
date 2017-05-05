'''
High level interface to Orbkit reading functions
'''

import re
from orbkit.read.molden import read_molden
from orbkit.read.gamess import read_gamess
from orbkit.read.gaussian_fchk import read_gaussian_fchk
from orbkit.read.gaussian_log import read_gaussian_log
from orbkit.read.aomix import read_aomix
from orbkit.read.wfx import read_wfx
from orbkit.read.wfn import read_wfn
from orbkit.read.cclib import read_with_cclib
from orbkit.display import display

readers = {'fchk': read_gaussian_fchk, 'wfx': read_wfx, 'wfn': read_wfn, 
           'cclib': read_with_cclib, 'molden': read_molden, 'gamess': gaussian_log, 
           'aomix': read_aomix}


def read(fname, all_mo=False,spin=None, ignore_molden=False, cclib_parser= None, **kwargs):
  '''
  This is the high-lever interface for the
  orbkit reading routines.
  
  **Parameters:**
  
  fname: str
    Specifies the filename for the input file.
  all_mo : bool, optional
    If True, all molecular orbitals are returned.
  spin : {None, 'alpha', or 'beta'}, optional
    If not None, returns exclusively 'alpha' or 'beta' molecular orbitals.
  cclib_parser : str
    If a cclib file is read, specifies the cclib.parser.
  ignore_molden: bool, optional
    Molden files are problematic for automatic detection 
    as they might be appended to other types of files.
    Default bahavior is to first check if the file is a molden file
    which will lead to other output being disregarded.
    To ignore the possible presence of molden data use ignore_molden = True.
    
  **Note:**
  
    All additional keyword arguments are forwarded to the reading functions.
    
  **Returns:**
  
    qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, mo_spec, etot :
    See :ref:`Central Variables` for details.
  '''
  
  filetype = find_filetype(fname, ignore_molden)
  
  if filetype == 'cclib':
    assert cclib_parser not None, 'For cclib files a parser must be specified'
  
  display('Loading data from {0} type file {1}\n'.format(filetype, fname))
  
  if filetype =! 'cclib':
    return readers[filetype](fname, all_mo=all_mo, spin=spin, **kwargs)
  else:
    return readers[filetype](fname, cclib_parser=cclib_parser, all_mo=all_mo, spin=spin, **kwargs)
  

def find_filetype(fname, ignore_molden=False):
  '''
  This function is used by the high-level read
  to determine what reader to use.
  Filetypes are determined either by extension or
  by looking for magic strings within the files.
  
  **Parameters:**
  
  fname: str
    Specifies the filename for the input file.
  ignore_molden: bool, optional
    Molden files are problematic for automatic detection 
    as they might be appended to other types of files.
    Default bahavior is to first check if the file is a molden file
    which will lead to other output being disregarded.
    To ignore the possible presence of molden data use ignore_molden = True.
    
  **Returns:**
  
  filetype, str
    Filetype for the file specified by fname
  
  filetypes determied from file extension:
    - .fchk
    - .wfx
    - .wfn
    - .cclib
    
  filetypes determined from magic strings:
    - Molden
    - Gamess US
    - Gaussian log
    - AOMix
  '''
  
  extensions = ['fchk', 'wfx', 'wfn', 'cclib']
  if fname.split('.')[-1] in extensions:
    return fname.split('.')[-1]

  molden_regex = re.compile(r"\[[ ]{,}Molden[ ]+Format[ ]{,}\]")
  gamess_regex = re.compile(r"GAMESS") #This might be too weak - Can someone who knows Gamess please check?
  gaussian_regex = re.compile(r"Copyright[,\s\(\)c0-9]+Gaussian\s{,},\s+Inc.")
  aomix_regex = re.compile(r"\[[ ]{,}AOMix[ ]+Format[ ]{,}\]")
  
  regexes = {'molden': molden_regex, 'gamess': gamess_regex, 'gaussian_log': gaussian_regex, 'aomix': aomix_regex}
  
  if ignore_molden:
    filetypes = ['gamess', 'gaussian_log', 'aomix']
  else:
    filetypes = ['molden', 'gamess', 'gaussian_log', 'aomix']  
    
  with open(fname, 'r') as fd:
    text = fd.read()
    for regname in filetypes:
      if regname != 'molden':
          text = text[:1e3] #Only Molden files can be appended so the entire file has to be scanned.
      if regexes['regname'].search(text):
        return regname
  raise NotImplementedError('File format not reccognized or reader not implemented!')

