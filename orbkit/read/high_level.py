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

readers = {'fchk': read_gaussian_fchk, 'wfx': read_wfx, 'wfn': read_wfn,
           'molden': read_molden, 'gamess': gaussian_log, 'aomix': read_aomix}

def find_filetype(fname, ignore_molden=False):
  '''
  This function is used by the high-level read
  to determine what reader to use.
  Filetypes are determined either by extension or
  by looking for magic strings within the files.
  cclib files require the specification of a parser 
  and are not supported
  
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
    
  filetypes determined from magic strings:
    - Molden
    - Gamess US
    - Gaussian log
    - AOMix
  '''
  if 'cclib' in fname.split('.')[-1]:
    print('Higl-level interface does not currently support cclib files')
    print('Use the dedicated read_cclib() function')
    exit()
  
  extensions = ['fchk', 'wfx', 'wfn']
  if fname.split('.')[-1] in extensions:
    return fname.split('.')[-1]

  molden_regex = re.compile(r"\[[ ]{,}Molden[ ]+Format[ ]{,}\]")
  gamess_regex = re.compile(r"GAMESS") #This might be too weak - Can someone who knows Gamess please check?
  gaussian_regex = re.compile(r"r"Copyright[,\s\(\)c0-9]+Gaussian\s{,},\s+Inc.")
  aomix_regex = re.compile(r"\[[ ]{,}AOMix[ ]+Format[ ]{,}\]")
  
  regexes = {'molden': molden_regex, 'gamess': gamess_regex, 'gaussian_log': gaussian_regex, 'aomix': aomix_regex}
  
  if ignore_molden:
    filetypes = ['gamess', 'gaussian_log', 'aomix']
  else:
    filetypes = ['molden', 'gamess', 'gaussian_log', 'aomix']  
    
  with open(fname, 'r') as fd:
    text = fd.read()
    for regname in filetypes:
      if regexes['regname'].search(text):
        return regname
