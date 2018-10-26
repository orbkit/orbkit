import numpy
import re

from orbkit.qcinfo import QCinfo
from orbkit.orbitals import AOClass, MOClass
from orbkit.display import display
from orbkit.tools import l_deg, lquant, orbit

from .tools import descriptor_from_file

'''
New Molden interface
'''

FLAGS_SPH  = ['5d',  '7f',  '9g']
FLAGS_CART = ['6d', '10f', '15g']
# regex to match all lines like "[5d7f]"
regex_flagline = re.compile(r'\[((' + '|'.join(FLAGS_SPH + FLAGS_CART) + r')+)\]')
regex_flag = re.compile('(\d+[dfg])')
regex_molden = re.compile(r'\[[ ]{,}[Mm]olden[ ]+[Ff]ormat[ ]{,}\]')

def read_molden(fname, all_mo=False, spin=None, i_md=-1, interactive=True,
                **kwargs):
  '''Reads all information desired from a molden file.
  
  **Parameters:**
  
    fname: str, file descriptor
      Specifies the filename for the input file.
      fname can also be used with a file descriptor instad of a filename.
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

  if 'index' not in kwargs.keys():
    kwargs['index'] = 0

  if isinstance(fname, str):
    filename = fname
    fname = descriptor_from_file(filename, index=kwargs['index'])
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
  spherical_basis = []
  cartesian_basis = []
  by_orca = []
  count = 0
  # Go through the file line by line 
  for line in flines:
    # Check the file for keywords
    if regex_molden.search(line):
      count += 1
      has_alpha.append(False)
      has_beta.append(False)
      restricted.append(False)
      spherical_basis.append([])
      cartesian_basis.append([])
      by_orca.append(False)
    if 'orca' in line.lower():
      by_orca[-1] = True
    # check whether line contains flags for spherical/cartesian basis functions
    m = regex_flagline.match(line.lower())
    if m:
      # get list of all flags in line
      flags = regex_flag.findall(m.group(1))
      # check whether cartesian or spherical
      for flag in flags:
        if flag in FLAGS_SPH:
          spherical_basis[-1].append(flag)
        if flag in FLAGS_CART:
          cartesian_basis[-1].append(flag)
      # check for ambiguous flags
      lsph = [l[-1] for l in spherical_basis[-1]]
      lcart = [l[-1] for l in cartesian_basis[-1]]
      if set(lsph) & set(lcart):
          raise IOError('The input file {} contains ambiguous flags for spherical and cartesian basis functions: {}'.format(filename, ', '.join(spherical_basis[-1]+cartesian_basis[-1])))
    if 'Spin' in line and 'alpha' in line.lower():
      has_alpha[-1] = True
    if 'Spin' in line and 'beta' in line.lower():
      has_beta[-1] = True
    if 'Occup' in line:
      restricted[-1] = restricted[-1] or (float(line.split('=')[1]) > 1.+1e-4)
  
  if count == 0:
    raise IOError('The input file %s is no valid molden file!\n\nIt does' % filename +
            ' not contain the keyword: [Molden Format]\n')
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
  for line in flines:
    thisline = line.split()        # The current line split into segments
    
    # Check the file for keywords 
    if regex_molden.search(line):
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
          angstrom = True
        else:
          angstrom = False
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
      elif '[' in line:
        sec_flag = None
      else:
        # Check if we are in a specific section 
        if sec_flag == 'geo_info' and thisline != []:
          # Geometry section 
          qc.geo_info.append(thisline[0:3])
          qc.geo_spec.append([float(ii) for ii in thisline[3:]])
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
            ao_type = thisline[0].lower()    # Which type of atomic orbital do we have
            pnum = int(thisline[1])  # Number of primatives
            # Calculate the degeneracy of this AO and increase basis_count 
            for i_ao in ao_type:
              # Calculate the degeneracy of this AO and increase basis_count 
              # TODO: check for mixed sph/cart basis 
              cart = True
              if cartesian_basis[i_md]:
                cart = i_ao in [f[-1] for f in cartesian_basis[i_md]]
              if spherical_basis[i_md]:
                cart = not i_ao in [f[-1] for f in spherical_basis[i_md]]
              basis_count += l_deg(lquant[i_ao], cartesian_basis=cart)
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
                try:
                  a = re.search(r'\d+', info[1]).group()
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
                raise ValueError('Error in coefficient %d of MO %s!' % (index, 
                      qc.mo_spec[-1]['sym']) + '\nSetting this coefficient to zero...')

  # check for mixed spherical/cartesian basis functions
  if max_l >= 2:
    # remove flags for unused angular momentum
    l = orbit[2:max_l+1]
    sph = [f for f in spherical_basis[i_md] if f[-1] in l]
    cart = [f for f in cartesian_basis[i_md] if f[-1] in l]
    if sph and cart:
      raise IOError('''The input file {} contains mixed spherical and Cartesian function ({}).
                  ORBKIT does not support these basis functions yet.
                  Pleas contact us, if you need this feature!'''.format(filename, ', '.join(sph+cart))) 

  # Spherical basis?
  if spherical_basis[i_md]:
    qc.ao_spec.set_lm_dict(p=[1,0])

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
  qc.format_geo(is_angstrom=angstrom)
  
  # Check the normalization
  from orbkit.analytical_integrals import get_ao_overlap
  spher_tmp = qc.ao_spec.spherical
  qc.ao_spec.spherical = False 
  norm = numpy.diagonal(get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec))
  qc.ao_spec.spherical = spher_tmp 
  if max(numpy.abs(norm-1.)) > 1e-5:
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
      cca = ommited_cca_norm(qc.ao_spec.get_lxlylz())
      for mo in qc.mo_spec:
        mo['coeffs'] *= cca

  qc.mo_spec.update()
  qc.ao_spec.update()
  return qc
