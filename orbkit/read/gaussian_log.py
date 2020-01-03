import numpy

from orbkit.qcinfo import QCinfo
from orbkit.orbitals import AOClass, MOClass
from orbkit.tools import l_deg, lquant
from orbkit.display import display

from .tools import descriptor_from_file

def read_gaussian_log(fname,all_mo=False,spin=None,orientation='standard',
                      i_link=-1,i_geo=-1,i_ao=-1,i_mo=-1,interactive=True,
                      **kwargs):
  '''Reads all information desired from a Gaussian .log file.

  **Parameters:**
  
    fname: str, file descriptor
      Specifies the filename for the input file.
      fname can also be used with a file descriptor instad of a filename.
    all_mo :  bool, optional
      If True, all molecular orbitals are returned.
    spin : {None, 'alpha', or 'beta'}, optional
      If not None, returns exclusively 'alpha' or 'beta' molecular orbitals.
    orientation : string, choices={'input', 'standard'}, optional
      Specifies orientation of the molecule in Gaussian nomenclature. [#first]_ 
    i_link : int, default=-1
      Selects the file for linked Gaussian jobs.
    i_geo : int, default=-1
      Selects the geometry section of the output file.
    i_ao : int, default=-1
      Selects the atomic orbital section of the output file.
    i_mo : int, default=-1
      Selects the molecular orbital section of the output file.
    interactive : bool
      If True, the user is asked to select the different sets.
  
  **Returns:**
  
    qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, mo_spec, etot :
        See :ref:`Central Variables` for details.

.. [#first] Attention: The MOs in the output are only valid for the standard orientation!

  '''
  
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
  
  # Search the file the specific sections
  count = {'link': 0, 'geometry': 0, 'geometry_input': 0, 'atomic orbitals': 0, 
           'molecular orbitals': [], 'state': []}

  def check_sel(count,i,interactive=False,default=-1):
    if count == 0:
      raise IndexError
    elif count == 1:
      return 0
    message = '\tPlease give an integer from 0 to {0} (default: {0}): '.format(count-1)
    
    try:
      if interactive:
        i = input(message)
        i = default if i == '' else int(i)
      i = range(count)[i]
    except (IndexError,ValueError):
      raise IOError(message.replace(':','!'))
    else:
      display('\tSelecting the %s' %
             ('last element.' if (i == count-1) else 'element %d.' % i))
    return i
  
  # Go through the file line by line 
  for il in range(len(flines)):
      line = flines[il]            # The current line as string  
      # Check the file for keywords 
      if ' Entering Link 1' in line:
        count['link'] += 1
  
  try:
    display('\tFound %d linked GAUSSIAN files.' % count['link'])
    i_link = check_sel(count['link'],i_link,interactive=interactive)
  except IndexError:
    raise IOError('Found no `Entering Link 1` keyword!')
  
  cartesian_basis = True
  c_link = 0
  # Go through the file line by line 
  for il in range(len(flines)):
      line = flines[il]            # The current line as string
      thisline = line.split()      # The current line split into segments
      
      
      # Check the file for keywords 
      if ' Entering Link 1' in line:
        c_link += 1
      if i_link == (c_link-1):
        if ' orientation:' in line:
          if '%s orientation:' % orientation in line.lower():
            count['geometry'] += 1
          if 'input orientation:' in line.lower():
            count['geometry_input'] += 1
        elif 'Standard basis:' in line or 'General basis read from cards:' in line:
          # Check if a cartesian basis has been applied
          if '(5D, 7F)' in line:
            cartesian_basis = False
          elif '(6D, 10F)' not in line:
            raise IOError('Please apply a Spherical Harmonics (5D, 7F) or '+
                          'a Cartesian Gaussian Basis Set (6D, 10F)!')
        elif 'AO basis set in the form of general basis input' in line:
          count['atomic orbitals'] += 1
        elif 'The electronic state is ' in line:
          count['state'].append(thisline[-1][:-1])
        elif 'Orbital Coefficients:' in line:
          mo_type = thisline[0]
          if mo_type != 'Beta':
            count['molecular orbitals'].append(mo_type)
          else:
            count['molecular orbitals'][-1] = 'Alpha&Beta'
 
  display('\nContent of the GAUSSIAN .log file:')
  display('\tFound %d geometry section(s). (%s orientation)' % 
         (count['geometry'], orientation))
  try:
    i_geo = check_sel(count['geometry'],i_geo,interactive=interactive)
  except IndexError:
    count['geometry'] = count['geometry_input']
    orientation = 'input'
    display('\Looking for "Input orientation": \n' + 
            '\tFound %d geometry section(s). (%s orientation)' % 
            (count['geometry'], orientation))
    try:
      i_geo = check_sel(count['geometry'],i_geo,interactive=interactive)
    except IndexError:
      raise IOError('Found no geometry section!'+
                    ' Are you sure this is a GAUSSIAN .log file?')
    
  
  try:
    display('\tFound %d atomic orbitals section(s) %s.' % 
            (count['atomic orbitals'],
            '(6D, 10F)' if cartesian_basis else '(5D, 7F)'))
    i_ao = check_sel(count['atomic orbitals'],i_ao,interactive=interactive)
  except IndexError:
    raise IOError('Write GFINPUT in your GAUSSIAN route section to print' + 
                  ' the basis set information!')
  
  try:
    display('\tFound the following %d molecular orbitals section(s):' % 
            len(count['molecular orbitals']))
  except IndexError:
    raise IOError('Write IOP(6/7=3) in your GAUSSIAN route section to print\n' + 
                  ' all molecular orbitals!')
  for i,j in enumerate(count['molecular orbitals']):
    string = '\t\tSection %d: %s Orbitals'% (i,j)
    try:
      string += ' (electronic state: %s)' % count['state'][i]
    except IndexError:
      pass
    display(string)
  i_mo = check_sel(len(count['molecular orbitals']),i_mo,interactive=interactive)
  
  if spin is not None:
    if spin != 'alpha' and spin != 'beta':
      raise IOError('`spin=%s` is not a valid option' % spin)
    else:
      display('Reading only molecular orbitals of spin %s.' % spin)

  # Set a counter for the AOs 
  basis_count = 0
  
  # Initialize some variables 
  sec_flag = None
  skip = 0
  c_link = 0
  c_geo = 0
  c_ao = 0
  c_mo = 0
  c_sao = 0
  old_ao = -1
  orb_sym = []
  qc = QCinfo()
  qc.ao_spec = AOClass([])
  qc.mo_spec = MOClass([])
  index = []
  
  # Go through the file line by line 
  for il in range(len(flines)):
    line = flines[il]              # The current line as string
    thisline = line.split()        # The current line split into segments
    
    # Check the file for keywords 
    if ' Entering Link 1' in line:
      c_link += 1
    if i_link == (c_link-1):
      if '%s orientation:' % orientation in line.lower():
        # The section containing information about 
        # the molecular geometry begins 
        if i_geo == c_geo:
          qc.geo_info = []
          qc.geo_spec = []
          sec_flag = 'geo_info'
        c_geo += 1
        skip = 4
      elif 'Standard basis:' in line or 'General basis read from cards:' in line:
        # Check if a cartesian basis has been applied
        if '(5D, 7F)' in line:
          cartesian_basis = False
        elif '(6D, 10F)' not in line:
          raise IOError('Please apply a Spherical Harmonics (5D, 7F) or '+
                        'a Cartesian Gaussian Basis Sets (6D, 10F)!')
      elif 'AO basis set in the form of general basis input' in line:
        # The section containing information about 
        # the atomic orbitals begins
        if i_ao == c_ao:
          qc.ao_spec = AOClass([])
          if not cartesian_basis:
            qc.ao_spec.spherical = True
          sec_flag = 'ao_info'
          basis_count = 0
        c_ao += 1
        bNew = True                  # Indication for start of new AO section
      elif 'Orbital symmetries:' in line:
          sec_flag = 'mo_sym'
          add = ''
          orb_sym = []
      elif 'Orbital Coefficients:' in line:
        # The section containing information about 
        # the molecular orbitals begins 
        if (i_mo == c_mo):
          sec_flag = 'mo_info'
          mo_type = count['molecular orbitals'][i_mo]
          qc.mo_spec = MOClass([])
          offset = 0
          add = ''
          orb_spin = []
          if orb_sym == []:
            if 'Alpha' in mo_type:
              add = '_a'
              orb_spin = ['alpha'] * basis_count
            orb_sym = ['A1'+add] * basis_count
            if 'Beta' in mo_type:
              add = '_b'
              orb_spin += ['beta'] * basis_count
              orb_sym += ['A1'+add] * basis_count
          for i in range(len(orb_sym)):
            # for numpy version < 1.6 
            c = ((numpy.array(orb_sym[:i+1]) == orb_sym[i]) != 0).sum()
            # for numpy version >= 1.6 this could be used:
            #c = numpy.count_nonzero(numpy.array(orb_sym[:i+1]) == orb_sym[i])
            qc.mo_spec.append({'coeffs': numpy.zeros(basis_count),
                            'energy': 0.,
                            'sym': '%d.%s' % (c,orb_sym[i])})
            if orb_spin != []:
              qc.mo_spec[-1]['spin'] = orb_spin[i]
        if mo_type != 'Beta':
          c_mo += 1
        bNew = True                  # Indication for start of new MO section
      elif 'E(' in line:
        qc.etot = float(line.split('=')[1].split()[0])
      else:
        # Check if we are in a specific section 
        if sec_flag == 'geo_info':
          if not skip:
            qc.geo_info.append([thisline[1],thisline[0],thisline[1]])
            qc.geo_spec.append([float(ij) for ij in thisline[3:]])
            if '-----------' in flines[il+1]:
              sec_flag = None
          else:
            skip -= 1
        if sec_flag == 'ao_info':
          # Atomic orbital section 
          if ' ****' in line: 
            # There is a line with stars after every AO 
            bNew = True
            # If there is an additional blank line, the AO section is complete
            if flines[il+1].split() == []:
              sec_flag = None
          elif bNew:
            # The following AOs are for which atom? 
            bNew = False
            at_num = int(thisline[0]) - 1
            ao_num = 0
          elif len(thisline) == 4:
            # AO information section 
            # Initialize a new dict for this AO 
            ao_num = 0               # Initialize number of atomic orbiatls 
            ao_type = thisline[0].lower()   # Type of atomic orbital            
            pnum = int(thisline[1])  # Number of primatives
            for i_ao in ao_type:
              # Calculate the degeneracy of this AO and increase basis_count 
              basis_count += l_deg(lquant[i_ao],cartesian_basis=cartesian_basis)
              qc.ao_spec.append({'atom': at_num,
                              'type': i_ao,
                              'pnum': pnum,
                              'coeffs': numpy.zeros((pnum, 2))
                              })
              if not cartesian_basis:
                qc.ao_spec[-1]['lm'] = []
          else:
            # Append the AO coefficients 
            coeffs = numpy.array(line.replace('D','e').split(), dtype=numpy.float64)
            for i_ao in range(len(ao_type)):
              qc.ao_spec[-len(ao_type)+i_ao]['coeffs'][ao_num,:] = [coeffs[0],
                                                                coeffs[1+i_ao]]
            ao_num += 1
        if sec_flag == 'mo_sym':
          if 'electronic state' in line:
            sec_flag = None
          else:
            info = line[18:].replace('(','').replace(')','').split()
            if 'Alpha' in line:
              add = '_a'
            elif 'Beta' in line:
              add = '_b'
            for i in info:
              orb_sym.append(i + add)   
        if sec_flag == 'mo_info':
          # Molecular orbital section 
          info = line[:21].split()          
          if info == []:
            coeffs = line[21:].split()
            if bNew:
              index = [offset+i for i in range(len(coeffs))]
              bNew = False
            else:
              for i,j in enumerate(index):
                qc.mo_spec[j]['occ_num'] = int('O' in coeffs[i])
                if mo_type not in 'Alpha&Beta':
                  qc.mo_spec[j]['occ_num'] *= 2
          elif 'Eigenvalues' in info:
            coeffs = line[21:].replace('-',' -').split()
            if mo_type == 'Natural':
              key = 'occ_num'
            else:
              key = 'energy'
            for i,j in enumerate(index):
              qc.mo_spec[j][key] = float(coeffs[i])
          else:
            try:
              int(info[0])
            except ValueError:
              for j in list(range(index[-1]+1,len(qc.mo_spec)))[::-1]:
                del qc.mo_spec[j]
              sec_flag = None
              orb_sym = []
              bNew = True
              continue
            coeffs = line[21:].replace('-',' -').split()
            if not cartesian_basis and offset == 0:
              if old_ao != line[:14].split()[-1] or len(line[:14].split()) == 4:
                old_ao = line[:14].split()[-1]
                c_sao += 1
              i = c_sao-1
              l = lquant[line[13].lower()] 
              m = line[14:21].replace(' ', '').lower()
              p = 'yzx'.find(m) if len(m) == 1 else -1
              if p != -1:
                m = p - 1
              elif m == '':
                m = 0
              else:
                m = int(m)
                
              qc.ao_spec[i]['lm'].append((l,m))
            for i,j in enumerate(index):
              qc.mo_spec[j]['coeffs'][int(info[0])-1] = float(coeffs[i])
            if int(info[0]) == basis_count:
              bNew = True
              offset = index[-1]+1
              if index[-1]+1 == len(orb_sym):
                sec_flag = None
                orb_sym = []
  
  # Are all MOs requested for the calculation? 
  if not all_mo:
    for i in range(len(qc.mo_spec))[::-1]:
      if qc.mo_spec[i]['occ_num'] < 0.0000001:
        del qc.mo_spec[i]
    
  if spin is not None:
    if orb_spin == []:
      raise IOError('You requested `%s` orbitals, but None of them are present.'
                    % spin)
    else:
      for i in range(len(qc.mo_spec))[::-1]:
        if qc.mo_spec[i]['spin'] != spin:
          del qc.mo_spec[i]
  
  # Convert geo_info and geo_spec to numpy.ndarrays
  qc.format_geo(is_angstrom=True)
  qc.mo_spec.update()
  qc.ao_spec.update()
  return qc
