'''
Some tools needed by Orbkit functions
'''

from os import path
import re
import numpy

from orbkit.display import display
from orbkit.tools import *
from orbkit.units import u_to_me
from orbkit import options
from orbkit import analytical_integrals

nist_mass = None
# Standard atomic masses as "Linearized ASCII Output", see http://physics.nist.gov
nist_file = path.join(path.dirname(path.realpath(__file__)),
                      '../supporting_data/Atomic_Weights_NIST.html')
# see http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii2&isotype=some

def descriptor_from_file(filename, index=0, ci_descriptor=False):
  from .tar import is_tar_file, get_file_from_tar

  if is_tar_file(filename):
    fname, _ = get_file_from_tar(filename, index=index, ci_descriptor=ci_descriptor)
  else:
    fname = open(filename, 'r')
  return fname

def check_mo_norm(qc):
  geo_spec = qc.geo_spec
  ao_spec = qc.ao_spec

  aoom = analytical_integrals.get_ao_overlap(geo_spec, geo_spec, ao_spec)
  moom = analytical_integrals.get_mo_overlap_matrix(mo_spec, mo_spec, aoom, numproc=options.numproc)

  deviation = numpy.linalg.norm(moom - numpy.eye(len(moom)))
  return deviation

def basissize(types, cart_spher):
  assert cart_spher in ['cartesian', 'pure'], 'Must specify eigther cartesian or pure spherical harmonics'
  size = 0
  for i in range(len(types)):
    l = spdfg_to_l[types[i]]
    if cart_spher == 'cartesian':
      size += (l+1)*(l+2)/2
    else:
      size += 2*l+1
  return size

def read_nist():
  '''Reads and converts the atomic masses from the "Linearized ASCII Output", 
  see http://physics.nist.gov.
  '''
  global nist_mass

  f = open(nist_file,'r')
  flines = f.readlines()
  f.close()
  
  nist_mass = []
  index = None
  new = True
  
  def rm_brackets(text,rm=['(',')','[',']']):
    for i in rm:
      text = text.replace(i,'')
    return text
  
  for line in flines:
    thisline = line.split()
    if 'Atomic Number =' in line:
      i = int(thisline[-1]) - 1
      new = (i != index)
      if new:
        nist_mass.append(['',0])
      index = i
    elif 'Atomic Symbol =' in line and new:
      nist_mass[index][0] = thisline[-1]
    elif 'Standard Atomic Weight =' in line and new:
      nist_mass[index][1] = float(rm_brackets(thisline[-1]))

def spin_check(spin,restricted,has_alpha,has_beta):
  '''Check if `spin` keyword is valid.
  '''
  if spin is not None:
    if restricted:
      raise IOError('The keyword `spin` is only supported for unrestricted calculations.')    
    if spin != 'alpha' and spin != 'beta':
      raise IOError('`spin=%s` is not a valid option' % spin)
    elif spin == 'alpha' and has_alpha:
      display('Reading only molecular orbitals of spin alpha.')
    elif spin == 'beta' and has_beta:
      display('Reading only molecular orbitals of spin beta.')
    elif (not has_alpha) and (not has_beta):
      raise IOError(
          'Molecular orbitals in the input file do not contain `Spin=` keyword')
    elif ((spin == 'alpha' and not has_alpha) or 
          (spin == 'beta' and not has_beta)):
      raise IOError('You requested `%s` orbitals, but None of them are present.'
                    % spin)

def standard_mass(atom):
  '''Returns the standard atomic mass of a given atom.
    
  **Parameters:**
  
  atom : int or str
    Contains the name or atomic number of the atom.
  
  **Returns:**
  
  mass : float
    Contains the atomic mass in atomic units.
  '''
  if nist_mass is None:
    read_nist()  
  try:
    atom = int(atom) - 1
    return nist_mass[atom][1] * u_to_me
  except ValueError:
    return dict(nist_mass)[atom.title()] * u_to_me
    
def get_atom_symbol(atom):
  '''Returns the atomic symbol of a given atom.
    
  **Parameters:**
  
  atom : int or str
    Contains the atomic number of the atom.
  
  **Returns:**
  
  symbol : str
    Contains the atomic symbol.
  '''
  if nist_mass is None:
    read_nist()  
  try:
    atom = int(atom) - 1
    return nist_mass[atom][0]
  except ValueError:    
    return atom.upper()

def set_ao_spherical(ao_spec,p=[1,0]):
  ao_spec.spherical = True
  for i,ao in enumerate(ao_spec):
    ii = ao['type']
    l = lquant[ii]
    for m in (range(0,l+1) if l != 1 else p):
      ao_spec[i]['ao_spherical'].append((l,m))
      if m != 0:
        ao_spec[i]['ao_spherical'].append((l,-m))
    for m in (range(1,l+1) if l != 1 else p):
      if m != 0:
        ao_spherical.append([i,(l,-m)])
  return

def find_itype(fname, extension=None):
  '''
  This function is used by the high-level read
  to determine what reader to use.
  Filetypes are determined either by extension or
  by looking for magic strings within the files.
  
  **Parameters:**
  
  fname: str, file descriptor
    Specifies the filename for the input file.
    fname can also be used with a file descriptor instad of a filename.
  extension: str, optional
    If extension is not None it will be used to attempt determining
    filetypes by extension in place of fname.
    
  **Returns:**
  
  filetype, str
    Filetype for the file specified by fname
  
  filetypes determied from file extension:
    - .fchk
    - .wfx
    - .wfn
    - .npz
    - .hdf5 / .h5
    
  filetypes determined from magic strings:
    - Molden
    - Gamess US
    - Gaussian log
    - AOMix
  '''

  #Don't like this was_str stuff... but I don't know what to do
  if isinstance(fname, str):
    filename = fname
    fname = descriptor_from_file(filename, index=0)
    was_str = True
  else:
    filename = fname.name
    was_str = False

  if not extension:
    extension = filename.split('.')[-1]
  if extension.lower() in ['fchk', 'wfx', 'wfn']:
    return extension
  elif extension.lower() in ['numpy', 'npz', 'hdf5', 'h5']:
    return 'native'
  
  molden_regex = re.compile(r"\[[ ]{,}[Mm]olden[ ]+[Ff]ormat[ ]{,}\]")
  gamess_regex = re.compile(r"[Gg][Aa][Mm][Ee][Ss][Ss]") #This might be too weak - Can someone who knows Gamess please check?
  gaussian_regex = re.compile(r"[Cc]opyright[,\s\(\)c0-9]+[Gg]aussian\s{,},\s+Inc.")
  aomix_regex = re.compile(r"\[[ ]{,}[Aa][Oo][Mm]ix[ ]+[Ff]ormat[ ]{,}\]")
  
  regexes = {'molden': molden_regex, 'gamess': gamess_regex, 'gaussian_log': gaussian_regex, 'aomix': aomix_regex}
  
  itypes = ['molden', 'gamess', 'gaussian_log', 'aomix']  
  
  from io import TextIOWrapper
  if isinstance(fname, TextIOWrapper):
    text = fname.read()
  else:
    text = fname.read().decode("iso-8859-1")

  for regname in itypes:
    if regexes[regname].search(text):
      return regname

  if was_str:
    fname.close()
  
  raise NotImplementedError('File format not reccognized or reader not implemented!')
