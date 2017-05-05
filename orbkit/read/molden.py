'''
New Molden interface
'''

import re

def read_molden(filename, all_mo=False, spin=None, i_md=-1, interactive=True,
                **kwargs):
  '''Reads all information desired from a molden file.
  
  **Parameters:**
  
    filename : str
      Specifies the filename for the input file.
    all_mo : bool, optional
      If True, all molecular orbitals are returned.
    spin : {None, 'alpha', or 'beta'}, optional
      If not None, returns exclusively 'alpha' or 'beta' molecular orbitals.
    i_md : int, default=-1
      Selects the `[Molden Format]` section of the output file.
    interactive : bool
      If True, the user is asked to select the different sets.
  
  **Returns:**
  
    qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, mo_spec, etot :
        See :ref:`Central Variables` for details.
  '''

  fid    = open(filename,'r')      # Open the file
  flines = fid.readlines()         # Read the WHOLE file into RAM
  fid.close()                      # Close the file
  
  def check_sel(count,i,interactive=False):
    if count == 0:
      raise IndexError
    elif count == 1:
      return 0
    message = '\tPlease give an integer from 0 to {0}: '.format(count-1)
    
    try:
      if interactive:
        i = int(raw_input(message))
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
  cartesian_basis = []
  mixed_warning = []
  by_orca = []
  count = 0
  # Go through the file line by line 
  for il in range(len(flines)):
    line = flines[il]            # The current line as string
    
    # Check the file for keywords 
    if '[molden format]' in line.lower():
      count += 1
      has_alpha.append(False)
      has_beta.append(False)
      restricted.append(False)
      cartesian_basis.append(True)
      mixed_warning.append(False)
      by_orca.append(False)
    if 'orca' in line.lower():
      by_orca[-1] = True
    if '[5d]' in line.lower() or '[5d7f]' in line.lower():
      cartesian_basis[-1] = False
    if '[5d10f]'  in line.lower():
      mixed_warning[-1] = '5D, 10F'
      cartesian_basis[-1] = False
    if '[7f]'  in line.lower():
      mixed_warning[-1] = '6D, 7F'
      cartesian_basis[-1] = True
    if 'Spin' in line and 'alpha' in line.lower():
      has_alpha[-1] = True
    if 'Spin' in line and 'beta' in line.lower():
      has_beta[-1] = True
    if 'Occup' in line:
      restricted[-1] = restricted[-1] or (float(line.split('=')[1]) > 1.+1e-4)
  
  if count == 0:
    display('The input file %s is no valid molden file!\n\nIt does' % filename +
          ' not contain the keyword: [Molden Format]\n')
    raise IOError('Not a valid input file')
  else:
    if count > 1:
      display('\nContent of the molden file:')
      display('\tFound %d [Molden Format] keywords, i.e., ' % count + 
              'this file contains %d molden files.' % count)
    i_md = check_sel(count,i_md,interactive=interactive)
  
  if spin is not None:
    if restricted[i_md]:
      raise IOError('The keyword `spin` is only supported for unrestricted calculations.')    
    if spin != 'alpha' and spin != 'beta':
      raise IOError('`spin=%s` is not a valid option' % spin)
    elif spin == 'alpha' and has_alpha[i_md]:
      display('Reading only molecular orbitals of spin alpha.')
    elif spin == 'beta' and has_beta[i_md]:
      display('Reading only molecular orbitals of spin beta.')
    elif (not has_alpha[i_md]) and (not has_beta[i_md]):
      raise IOError(
           'Molecular orbitals in `molden` file do not contain `Spin=` keyword')
    elif ((spin == 'alpha' and not has_alpha[i_md]) or 
          (spin == 'beta' and not has_beta[i_md])):
      raise IOError('You requested `%s` orbitals, but None of them are present.'
                    % spin)
  
  # Set a counter for the AOs 
  basis_count = 0
  sym = {}

  # Declare synonyms for molden keywords 
  synonyms = {'Sym': 'sym',
              'Ene': 'energy',
              'Occup': 'occ_num',
              'Spin': 'spin'
             }
  MO_keys = synonyms.keys()
  
  count = 0
  max_l = 0
  start_reading = False
  # Go through the file line by line 
  for il in range(len(flines)):
    line = flines[il]              # The current line as string
    thisline = line.split()        # The current line split into segments
    
    # Check the file for keywords 
    if '[molden format]' in line.lower():
      # A new file begins 
      # Initialize the variables 
      if i_md == count:
        qc = QCinfo()
        sec_flag = False           # A Flag specifying the current section 
        start_reading = True       # Found the selected section
      else:
        start_reading = False
      count += 1
      continue
    if start_reading:
      if '_ENERGY=' in line:
        try:
          qc.etot = float(thisline[1])
        except IndexError:
          pass
      elif '[atoms]' in line.lower():
        # The section containing information about 
        # the molecular geometry begins 
        sec_flag = 'geo_info'
        if 'Angs' in line:
          # The length are given in Angstroem 
          # and have to be converted to Bohr radii --
          aa_to_au = 1/0.52917720859
        else:
          # The length are given in Bohr radii 
          aa_to_au = 1.0
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
        display('orbkit does not work for STOs!\nEXIT\n');
        raise IOError('Not a valid input file')
      elif '[' in line:
        sec_flag = None
      else:
        # Check if we are in a specific section 
        if sec_flag == 'geo_info' and thisline != []:
          # Geometry section 
          qc.geo_info.append(thisline[0:3])
          qc.geo_spec.append([float(ii)*aa_to_au for ii in thisline[3:]])
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
              basis_count += l_deg(lquant[i_ao],cartesian_basis=cartesian_basis[i_md])
              max_l = max(max_l,lquant[i_ao])
              qc.ao_spec.append({'atom': at_num,
                              'type': i_ao,
                              'pnum': -pnum if by_orca[i_md] else pnum,
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
                try:
                  a = search(r'\d+', info[1]).group()
                  if a == info[1]:
                    info[1] = '%s.1' % a
                  elif info[1].startswith(a):
                    info[1] = info[1].replace(a, '%s.' % a,1)
                  else:
                    raise AttributeError
                except AttributeError:
                  if info[1] not in sym.keys(): sym[info[1]] = 1
                  else: sym[info[1]] += 1
                  info[1] = '%d.%s' % (sym[info[1]],info[1]) 
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
                qc.mo_spec[-1]['coeffs'][index] = float(thisline[1])
              except ValueError:
                # If it cannot be converted print error message 
                display('Error in coefficient %d of MO %s!' % (index, 
                  qc.mo_spec[-1]['sym']) + '\nSetting this coefficient to zero...')
  
  # Spherical basis?
  if not cartesian_basis[i_md]:    
    qc.ao_spherical = get_ao_spherical(qc.ao_spec,p=[1,0])
  if max_l > 2 and mixed_warning[i_md]:
    display('='*80)
    display('The input file %s contains ' % filename +
            'mixed spherical and Cartesian function (%s).' %  mixed_warning[i_md] + 
            'ORBKIT does not support these basis functions yet. '+
            'Pleas contact us, if you need this feature!')    
    display('='*80)
  # Are all MOs requested for the calculation? 
  if not all_mo:
    for i in range(len(qc.mo_spec))[::-1]:
      if qc.mo_spec[i]['occ_num'] < 0.0000001:
        del qc.mo_spec[i]
  
  # Only molecular orbitals of one spin requested?
  if spin is not None:
    for i in range(len(qc.mo_spec))[::-1]:
      if qc.mo_spec[i]['spin'] != spin:
        del qc.mo_spec[i]

  if restricted[i_md]:
    # Closed shell calculation
    for mo in qc.mo_spec:
      del mo['spin']
  else:
    # Rename MOs according to spin
    for mo in qc.mo_spec:
      mo['sym'] += '_%s' % mo['spin'][0]
  
  # Orca uses for all molecular orbitals the same name 
  sym = [i['sym'] for i in qc.mo_spec]
  if sym[1:] == sym[:-1]:
    sym = sym[0].split('.')[-1]
    for i in range(len(qc.mo_spec)):
      qc.mo_spec[i]['sym'] = '%d.%s' % (i+1,sym)
  
  # Convert geo_info and geo_spec to numpy.ndarrays
  qc.format_geo()
  
  # Check the normalization
  from orbkit.analytical_integrals import get_ao_overlap,get_lxlylz
  norm = numpy.diagonal(get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec))
  
  if sum(numpy.abs(norm-1.)) > 1e-8:
    display('The atomic orbitals are not normalized correctly, renormalizing...\n')
    if not by_orca[i_md]: 
      j = 0
      for i in range(len(qc.ao_spec)):
        qc.ao_spec[i]['coeffs'][:,1] /= numpy.sqrt(norm[j])
        for n in range(l_deg(lquant[qc.ao_spec[i]['type']],cartesian_basis=True)):
          j += 1
    else:
      qc.ao_spec[0]['N'] = 1/numpy.sqrt(norm[:,numpy.newaxis])
  
    if cartesian_basis[i_md]:
      from orbkit.cy_overlap import ommited_cca_norm
      cca = ommited_cca_norm(get_lxlylz(qc.ao_spec))
      for mo in qc.mo_spec:
        mo['coeffs'] *= cca
  
  return qc
  # read_molden 


def spin_count(text):
  spin_regex = re.compile(r"\[[ ]{,}Molden[ ]+Format[ ]{,}\].*Spin[=\s\(]+\b{0}\b\s+\)".format(spin))
  spin_regex.findall(text)
  
