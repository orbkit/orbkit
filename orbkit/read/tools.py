'''
Some tools needed by Orbkit functions
'''

from orbkit.display import display
from os import path
import numpy
from orbkit.tools import lquant
from orbkit.units import u2me

nist_mass = None
# Standard atomic masses as "Linearized ASCII Output", see http://physics.nist.gov
nist_file = path.join(path.dirname(path.realpath(__file__)),
                      '../supporting_data/Atomic_Weights_NIST.html')
# see http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii2&isotype=some

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
    return nist_mass[atom][1] * u2me
  except ValueError:
    return dict(nist_mass)[atom.title()] * u2me
    
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

def get_ao_spherical(ao_spec,p=[1,0]):
  ao_spherical = []
  for i,ao in enumerate(ao_spec):
    ii = ao['type']
    l = lquant[ii]
    for m in (range(0,l+1) if l != 1 else p):
      ao_spherical.append([i,(l,m)])
      if m != 0:
        ao_spherical.append([i,(l,-m)])
    #for m in (range(1,l+1) if l != 1 else p):
      #if m != 0:
        #ao_spherical.append([i,(l,-m)])
  return ao_spherical

def mo_select(mo_spec, fid_mo_list, strict=False):
  '''Selects molecular orbitals from an external file or a list of molecular 
  orbital labels.

  **Parameters:**
   
    mo_spec :        
      See :ref:`Central Variables` for details.
    strict : bool, optional
      If True, orbkit will follow strictly the fid_mo_list, i.e., the order of 
      the molecular orbitals will be kept and multiple occurrences of items 
      will evoke multiple calculations of the respective molecular orbitals. 
    fid_mo_list : str, `'all_mo'`, or list
      | If fid_mo_list is a str, specifies the filename of the molecular orbitals list.
      | If fid_mo_list is 'all_mo', creates a list containing all molecular orbitals.
      | If fid_mo_list is a list, provides a list (or a list of lists) of molecular 
        orbital labels.

  **Supported Formats:**
  
    Integer List (Counting from **ONE**!)::
    
      1       2       3
      5       4
      homo    lumo+2:lumo+4
    
    List with Symmetry Labels::
    
      1.1     2.1     1.3
      1.1     4.1
      4.1     2.3     2.1
  
  **Returns:**
  
    Dictionary with following Members:
      :mo: - List of molecular orbital labels.
      :mo_ii: - List of molecular orbital indices.
      :mo_spec: - Selected elements of mo_spec. See :ref:`Central Variables` for details.
      :mo_in_file: - List of molecular orbital labels within the fid_mo_list file.
      :sym_select: - If True, symmetry labels have been used. 
  
  ..attention:
    
    For **unrestricted** calculations, orbkit adds `_a` (alpha) or `_b` (beta) to
    the symmetry labels, e.g., `1.1_a`. 
    If you have specified the option `spin=alpha` or `spin=beta`, only the 
    alpha or the beta orbitals are taken into account for the counting 
    within the Integer List.
  '''
  import re
  display('\nProcessing molecular orbital list...')
  
  mo_in_file = []
  selected_mo = []
  sym_select = False
  
  def assign_selected_mo(selected_mo,mo_spec,strict=False,  
                          what=lambda x,y: y[x]['sym']):
    selected_mo_spec = []
    selected_mo_ii = [] 
    for i in selected_mo:
      is_present = False
      for k in range(len(mo_spec)):
        if (what(k,mo_spec) == i):
          is_present = True
          if strict or (i not in selected_mo_ii):
            selected_mo_spec.append(mo_spec[k])
            selected_mo_ii.append(what(k,mo_spec))
      if not is_present:
        raise IOError('Cannot find %s in mo_spec' % i)
    selected_mo_ii = numpy.array(selected_mo_ii)
    return selected_mo_spec,selected_mo_ii
  
  def expr2int(expr):
    if isinstance(expr,int):
      return expr
    x = 0
    for i in re.findall(r'\d+|[+-]\d+',expr):
      x += int(i)
    return x
  
  def get_selection(selected_mo):
    mo_occup = numpy.array([i['occ_num'] for i in mo_spec])
    homo = (mo_occup>0.).nonzero()[0][-1]   + 1 # molden numbering
    lumo = (mo_occup>0.).nonzero()[0][-1]+1 + 1 # molden numbering
    mo_energy = numpy.array([i['energy'] for i in mo_spec])
    last_bound = sum(mo_energy<=0.0)            # molden numbering
    sel = []
    for i in selected_mo:
      i = i.lower().replace('homo',str(homo)).replace('lumo',str(lumo))
      i = i.replace('last_bound',str(last_bound))
      if ':' in i:
        k = [1,len(mo_spec)+1,1]
        i = i.split(':')
        for ik,j in enumerate(i):
          if j != '': k[ik] = j
        i = list(range(*[expr2int(j) for j in k]))
        sel.extend(i)
      else:
        sel.append(int(i))
    return sel
  
  if isinstance(fid_mo_list,str) and fid_mo_list.lower() == 'all_mo':
    selected_mo = numpy.array(numpy.arange(len(mo_spec))+1, dtype=numpy.str)
    mo_in_file = [selected_mo]
    selected_mo_spec = mo_spec
    selected_mo_ii = numpy.array([i['sym'] for i in selected_mo_spec])
  else:
    if isinstance(fid_mo_list,str) and not path.exists(fid_mo_list):
      if ',' in fid_mo_list:
        fid_mo_list = fid_mo_list.split(',')
      else:
        fid_mo_list = [fid_mo_list]
    if isinstance(fid_mo_list, list):
      for i in fid_mo_list:
        if not isinstance(i, list):
          i = i.split(',') if isinstance(i,str) else [i]
        selected_mo.extend(list(map(str,i)))
        mo_in_file.append(list(map(str,i)))
    else:
      try:
        fid=open(fid_mo_list,'r')
        flines = fid.readlines()
        fid.close()
        for line in flines:
          integer = line.replace(',',' ').split()
          mo_in_file.append(integer)
          selected_mo.extend(integer)
      except:
        raise IOError('The selected mo-list (%(m)s) is not valid!' % 
                      {'m': fid_mo_list} + '\ne.g.\n\t1\t3\n\t2\t7\t9\n')
    
    # Print some information
    for i,j in enumerate(mo_in_file):
      display('\tLine %d: %s' % (i+1,', '.join(j)))
    
    # Check if the molecular orbitals are specified by symmetry 
    # (e.g. 1.1 in MOLPRO nomenclature) or 
    # by the number in the input file (e.g. 1)
    
    try: # Try to convert selections into integer
      for i in selected_mo:
        if isinstance(i,int):
          continue
        i = i.replace('homo','1').replace('lumo','2').replace('last_bound','3')
        for r in ['-','+',':']:
          i = i.replace(r,'')
        int(i)
    except ValueError:
      sym_select = True
      errors = []
      for i in range(len(selected_mo)):
        if not '.' in selected_mo[i]:
          errors.append(i)      
      if errors:
        err = [selected_mo[i] for i in errors]
        raise IOError('`%s` are no valid labels according '% ', '.join(err) +
                      'to the MOLPRO nomenclature, e.g., `5.1` or `5.A1`.' +
                      '\n\tHint: You cannot mix integer numbering and MOLPRO\'s ' +
                      'symmetry labels')
    
    if sym_select:
      what = lambda x,y: y[x]['sym']
      selected_mo_spec,selected_mo_ii = assign_selected_mo(selected_mo,
                                                           mo_spec,
                                                           strict=strict,
                                                           what=what)
    else:
      selected_mo = get_selection(selected_mo)
      
      if not strict:
        selected_mo = list(map(int, selected_mo))            
        selected_mo.sort()
      selected_mo = list(map(str, selected_mo))
      what = lambda x,y: str(x+1)
      selected_mo_spec,selected_mo_ii = assign_selected_mo(selected_mo,
                                                           mo_spec,
                                                           strict=strict,
                                                           what=what)
      selected_mo = selected_mo_ii
      for i in range(len(mo_in_file)):
        mo_in_file[i] = list(map(str, get_selection(mo_in_file[i])))
    
    # Print some information
    display('\nThe following orbitals will be considered...')
    for i,j in enumerate(mo_in_file):
      display('\tLine %d: %s' % (i+1,', '.join(j)))
  
  display('')
  return {'mo': selected_mo, 'mo_ii': selected_mo_ii,
          'mo_spec': selected_mo_spec, 
          'mo_in_file': mo_in_file, 'sym_select': sym_select}
