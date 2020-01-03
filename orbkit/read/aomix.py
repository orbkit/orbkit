import re
import numpy

from orbkit.qcinfo import QCinfo
from orbkit.orbitals import AOClass, MOClass
from orbkit.tools import l_deg, lquant
from orbkit.display import display

from .tools import descriptor_from_file, spin_check

def read_aomix(fname, all_mo=False, spin=None, i_md=-1, interactive=True,
               created_by_tmol=True, **kwargs):
  '''Reads all information desired from a aomix file.
  
  **Parameters:**
  
  fname : str, file descriptor
    Specifies the filename for the input file.
    fname can also be used with a file descriptor instad of a filename.
  all_mo : bool, optional
    If True, all molecular orbitals are returned.
  spin : {None, 'alpha', or 'beta'}, optional
    If not None, returns exclusively 'alpha' or 'beta' molecular orbitals.
  i_md : int, default=-1
    Selects the `[AOMix Format]` section of the output file.
  interactive : bool
    If True, the user is asked to select the different sets.
  created_by_tmol : bool
    If True and if Cartesian basis set is found, the molecular orbital 
    coefficients will be converted.
  
  **Returns:**
  
  qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, mo_spec, etot :
    See :ref:`Central Variables` for details.
  '''
  
  aomix_regex = re.compile(r"\[[ ]{,}[Aa][Oo][Mm]ix[ ]+[Ff]ormat[ ]{,}\]")

  if isinstance(fname, str):
    filename = fname
    fname = descriptor_from_file(filename, index=0)
  else:
    filename = fname.name

  from io import TextIOWrapper
  if isinstance(fname, TextIOWrapper):
    flines = fname.readlines()       # Read the WHOLE file into RAM
  else:
    magic = 'This is an Orbkit magic string'
    text = fname.read().decode("iso-8859-1").replace('\n','\n{}'.format(magic))
    flines = text.split(magic)
    flines.pop()
  
  # Is this really a aomix file? 
  if not '[AOMix Format]\n' in flines:
    raise IOError('The input file %s is no valid aomix file!\n\nIt does'  % filename+
          ' not contain the keyword: [AOMix Format]\n')
  
  def check_sel(count,i,interactive=False):
    if count == 0:
      raise IndexError
    elif count == 1:
      return 0
    message = '\tPlease give an integer from 0 to %d: ' % (count-1)
    
    try:
      if interactive:
        i = int(input(message))
      i = range(count)[i]
    except (IndexError,ValueError):
      raise IOError(message.replace(':','!'))
    else:
      display('\tSelecting the %s' %
             ('last element.' if (i == count-1) else 'element %d.' % i))
    return i
  
  has_alpha = []
  has_beta = []
  restricted = []
  count = 0
  # Go through the file line by line 
  for il in range(len(flines)):
    line = flines[il]            # The current line as string
    
    # Check the file for keywords 
    if aomix_regex.search(line):
      count += 1
      has_alpha.append(False)
      has_beta.append(False)
      restricted.append(False)
    if 'Spin' in line and 'alpha' in line.lower():
      has_alpha[-1] = True
    if 'Spin' in line and 'beta' in line.lower():
      has_beta[-1] = True
    if 'Occup' in line:
      restricted[-1] = restricted[-1] or (float(line.split('=')[1]) > 1.+1e-4)
  
  if count == 0:
    raise IOError('The input file %s is no valid aomix file!\n\nIt does' % filename +
            ' not contain the keyword: [AOMix Format]\n')
  else:
    if count > 1:
      display('\nContent of the aomix file:')
      display('\tFound %d [AOMix Format] keywords, i.e., ' % count + 
              'this file contains %d aomix files.' % count)
    i_md = check_sel(count,i_md,interactive=interactive)
    
  spin_check(spin,restricted[i_md],has_alpha[i_md],has_beta[i_md])
  
  # Set a counter for the AOs 
  basis_count = 0

  # Declare synonyms for molden keywords 
  synonyms = {'Sym': 'sym',
              'Ene': 'energy',
              'Occup': 'occ_num',
              'Spin': 'spin'
             }
  MO_keys = synonyms.keys()
  
  lxlylz = []
  count = 0
  start_reading = False
  # Go through the file line by line 
  for il in range(len(flines)):
    line = flines[il]              # The current line as string
    thisline = line.split()        # The current line split into segments
    
    # Check the file for keywords 
    if '[aomix format]' in line.lower():
      # A new file begins 
      # Initialize the variables 
      if i_md == count:
        qc = QCinfo()
        qc.ao_spec = AOClass([])
        qc.mo_spec = MOClass([])
        sec_flag = False           # A Flag specifying the current section 
        start_reading = True       # Found the selected section
      else:
        start_reading = False
      count += 1
      continue
    if start_reading:
      if '[SCF Energy / Hartree]' in line:
        try:
          qc.etot = float(flines[il+1].split()[0])
        except IndexError:
          pass
      elif '[atoms]' in line.lower():
        # The section containing information about 
        # the molecular geometry begins 
        sec_flag = 'geo_info'
        angstrom = 'Angs' in line
      elif '[gto]' in line.lower():
        # The section containing information about 
        # the atomic orbitals begins 
        sec_flag = 'ao_info'
        bNew = True                  # Indication for start of new AO section
      elif '[mo]' in line.lower():
        # The section containing information about 
        # the molecular orbitals begins 
        sec_flag = 'mo_info'
        bNew = True                  # Indication for start of new MO section
      elif '[sto]' in line.lower():
        # The orbkit does not support Slater type orbitals 
        raise IOError('orbkit does not work for STOs!\nEXIT\n')
      else:
        # Check if we are in a specific section 
        if sec_flag == 'geo_info':
          # Geometry section 
          qc.geo_info.append(thisline[0:3])
          qc.geo_spec.append([float(ii) for ii in thisline[3:
              ]])
        if sec_flag == 'ao_info':
          # Atomic orbital section 
          def check_int(i):
            try:
              int(i)
              return True
            except ValueError:
              return False
          
          if thisline == []:
            # There is a blank line after every AO 
            bNew = True
          elif bNew:
            # The following AOs are for which atom? 
            bNew = False
            at_num = int(thisline[0]) - 1
            ao_num = 0
          elif len(thisline) == 3 and check_int(thisline[1]):
            # AO information section 
            # Initialize a new dict for this AO 
            ao_num = 0               # Initialize number of atomic orbiatls 
            ao_type = thisline[0]    # Which type of atomic orbital do we have
            pnum = int(thisline[1])  # Number of primatives
            # Calculate the degeneracy of this AO and increase basis_count 
            for i_ao in ao_type:
              # Calculate the degeneracy of this AO and increase basis_count 
              basis_count += l_deg(lquant[i_ao])
              qc.ao_spec.append({'atom': at_num,
                              'type': i_ao,
                              'pnum': pnum,
                              #'ao_spherical': None,
                              'coeffs': numpy.zeros((pnum, 2))
                              })
          else:
            # Append the AO coefficients 
            coeffs = numpy.array(line.replace('D','e').split(), dtype=numpy.float64)
            for i_ao in range(len(ao_type)):
              qc.ao_spec[-len(ao_type)+i_ao]['coeffs'][ao_num,:] = [coeffs[0],
                                                                coeffs[1+i_ao]]
            ao_num += 1
        if sec_flag == 'mo_info':
          # Molecular orbital section 
          if '=' in line:
            # MO information section 
            if bNew:
              # Create a numpy array for the MO coefficients and 
              # for backward compability create a simple counter for 'sym'
              qc.mo_spec.append({'coeffs': numpy.zeros(basis_count),
                                'sym': '%d.1' % (len(qc.mo_spec)+1)})
              bNew = False
            # Append information to dict of this MO 
            info = line.replace('\n','').replace(' ','')
            info = info.split('=')
            if info[0] in MO_keys:               
              if info[0] == 'Spin':
                info[1] = info[1].lower()
              elif info[0] != 'Sym':
                info[1] = float(info[1])
              elif not '.' in info[1]:
                from re import search
                a = search(r'\d+', info[1]).group()
                if a == info[1]:
                  info[1] = '%s.1' % a
                else:
                  info[1] = info[1].replace(a, '%s.' % a, 1)
              qc.mo_spec[-1][synonyms[info[0]]] = info[1]
          else:
            if ('[' or ']') in line:
              # start of another section that is not (yet) read
              sec_flag = None
            else:
              # Append the MO coefficients 
              bNew = True            # Reset bNew
              index = int(thisline[0])-1
              try: 
                # Try to convert coefficient to float 
                qc.mo_spec[-1]['coeffs'][index] = float(thisline[-1])
                if len(qc.mo_spec) == 1:
                  lxlylz.append(thisline[-2])
              except ValueError:
                # If it cannot be converted print error message 
                raise ValueError('Error in coefficient %d of MO %s!' % (index, 
                qc.mo_spec[-1]['sym']) + '\nSetting this coefficient to zero...')

  # Check usage of same atomic basis sets
  for ii in range(len(lxlylz)):
    s = lxlylz[ii]
    exp = [0,0,0]
    c_last = None
    for jj in s[1:]:
      try:
        c = int(jj)
        exp[c_last] += (c-1)
      except ValueError:
        for kk,ll in enumerate('xyz'):
          if jj == ll:
            exp[kk] += 1
            c_last = kk
    lxlylz[ii] = exp

  count = 0
  for i,j in enumerate(qc.ao_spec):    
    l = l_deg(lquant[j['type']])
    j['lxlylz'] = []
    for i in range(l):
      j['lxlylz'].append((lxlylz[count][0],
                          lxlylz[count][1],
                          lxlylz[count][2]))
      count += 1
    j['lxlylz'] = numpy.array(j['lxlylz'],dtype=numpy.int64)
  
  # For Cartesian basis sets in Turbomole, the molecular orbital coefficients 
  # have to be converted.
  is_tmol_cart = not (len(qc.mo_spec) % len(qc.mo_spec[0]['coeffs'])) 
  
  # Are all MOs requested for the calculation? 
  if not all_mo:
    for i in range(len(qc.mo_spec))[::-1]:
      if qc.mo_spec[i]['occ_num'] < 0.0000001:
        del qc.mo_spec[i] 
  
  # Modify qc.mo_spec to support spin
  qc.select_spin(restricted[i_md],spin=spin)
  
  # Convert geo_info and geo_spec to numpy.ndarrays
  qc.format_geo(is_angstrom=angstrom)
  
  if is_tmol_cart and created_by_tmol:
    display('\nFound a Cartesian basis set in the AOMix file.')
    display('We assume that this file has been created by Turbomole.')
    display('Applying a conversion to the molecular orbital coefficients, ')
    display('in order to get normalized orbitals.')
    
    # Convert MO coefficients
    def dfact(n):
      if n <= 0:
        return 1
      else:
        return n * dfact(n-2)

    mo = qc.mo_spec.get_coeffs()
    for i,j in enumerate(qc.ao_spec.get_lxlylz()):
      norm = (dfact(2*j[0] - 1) * dfact(2*j[1] - 1) * dfact(2*j[2] - 1))
      j = sum(j)
      if j >1: 
        mo[:,i] *= numpy.sqrt(norm)   
    for ii in range(len(qc.mo_spec)):
      qc.mo_spec[ii]['coeffs'] = mo[ii]

  qc.mo_spec.update()
  qc.ao_spec.update()
  return qc
