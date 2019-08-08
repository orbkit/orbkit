import numpy
import re
from .tools import *
from orbkit.qcinfo import QCinfo
from orbkit.display import display
from orbkit.orbitals import AOClass, MOClass

def read_gamess(fname, all_mo=False, spin=None, read_properties=False,
                **kwargs):
  '''Reads all information desired from a Gamess-US output file.

  **Parameters:**

  fname : str, file descriptor
    Specifies the filename for the input file.
    fname can also be used with a file descriptor instad of a filename.
  all_mo : bool, optional
      If True, all molecular orbitals are returned.

  **Returns:**

    qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, mo_spec, etot :
        See :ref:`Central Variables` for details.
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

  # Initialize the variables
  qc = QCinfo()
  qc.ao_spec = AOClass([])
  qc.mo_spec = MOClass([])
  has_alpha = False                  # Flag for alpha electron set
  has_beta = False                   # Flag for beta electron set
  restricted = True                  # Flag for restricted calculation
  sec_flag = None                    # A Flag specifying the current section
  is_pop_ana = True                  # Flag for population analysis for ground state
  keyword = [' ATOM      ATOMIC                      COORDINATES','']
                                     # Keywords for single point calculation and
                                     # geometry optimization
  mokey = 'EIGENVECTORS'             # Keyword for MOs
  unrestopt = False                  # Flag for unrestricted optimization
  bopt = False                       # Flag for geometry optimization
  sym={}                             # Symmetry of MOs
  geo_skip = 1                       # Number of lines to skip in geometry section

  for il in range(len(flines)):
    line = flines[il]               # The current line as string
    thisline = line.split()         # The current line split into segments
    # Check the file for keywords
    if 'RUNTYP=OPTIMIZE' in line:
      keyword = [' COORDINATES OF ALL ATOMS ARE',
                 '***** EQUILIBRIUM GEOMETRY LOCATED *****']
      geo_skip = 2
      bopt = True
      if 'SCFTYP=UHF' in line:
        mokey = ' SET ****'
        restricted = False
      else:
        mokey = 'EIGENVECTORS'

    elif keyword[0] in line and keyword[1] in flines[il-1]:
      # The section containing information about
      # the molecular geometry begins
      sec_flag = 'geo_info'
      atom_count = 0  # Counter for Atoms
      angstrom = not '(BOHR)' in line

    elif 'ATOMIC BASIS SET' in line:
      # The section containing information about
      # the atomic orbitals begins
      sec_flag = 'ao_info'
      ao_skip = 6                     # Number of lines to skip
      AO = []                         # Atomic orbitals

    elif '----- ALPHA SET ' in line:
      # The section for alpha electrons
      has_alpha = True
      has_beta = False
      restricted = False

    elif '----- BETA SET ' in line:
      # The section for alpha electrons
      restricted = False
      has_alpha = False
      has_beta = True

    elif mokey in line and len(thisline) < 3:
      # The section containing information about
      # the molecular orbitals begins
      sec_flag = 'mo_info'
      mo_skip = 1
      len_mo = 0                      # Number of MOs
      init_mo = False                 # Initialize new MO section
      info_key = None                 # A Flag specifying the energy and symmetry section
      lxlylz = []
      if 'ALPHA' in line:
        has_alpha = True
        mo_skip = 0
      elif 'BETA' in line:
        has_beta = True
        has_alpha = False
        mo_skip = 0

    elif 'NATURAL ORBITALS' in line and len(thisline) <= 3:
      display('The natural orbitals are not extracted.')

    elif ' NUMBER OF OCCUPIED ORBITALS (ALPHA)          =' in line:
      occ = []                        # occupation number of molecular orbitals
      occ.append(int(thisline[-1]))
    elif ' NUMBER OF OCCUPIED ORBITALS (BETA )          =' in line:
      occ.append(int(thisline[-1]))
#      elif 'ECP POTENTIALS' in line:
#        sec_flag = 'ecp_info'
#        ecp = ''
    elif ' NUMBER OF OCCUPIED ORBITALS (ALPHA) KEPT IS    =' in line:
      occ = []                        # occupation number of molecular orbitals
      occ.append(int(thisline[-1]))
    elif ' NUMBER OF OCCUPIED ORBITALS (BETA ) KEPT IS    =' in line:
      occ.append(int(thisline[-1]))
    elif 'NUMBER OF STATES REQUESTED' in line and read_properties:
      # get the number of excited states and initialize variables for
      # transition dipole moment and energies
      exc_states = int(line.split('=')[1])  # Number of excited states
      # Dipole moments matrix: Diagonal elements -> permanent dipole moments
      # Off-diagonal elements -> transition dipole moments
      qc.dipole_moments = numpy.zeros(((exc_states+1),(exc_states+1),3))
      # Multiplicity of ground and excited states
      qc.states['multiplicity'] = numpy.zeros(exc_states+1)
      # Energies of ground and excited states
      qc.states['energy'] = numpy.zeros(exc_states+1)
      qc.states['energy'][0]           = qc.etot
      qc.states['multiplicity'][0]     = gs_multi
      dm_flag = None                  # Flag specifying the dipole moments section
    elif 'TRANSITION DIPOLE MOMENTS' in line and read_properties:
      # Section containing energies of excited states
      sec_flag = 'dm_info'
        # Energy and Multiplicity for ground state
    elif 'SPIN MULTIPLICITY' in line and read_properties:
      # odd way to get gound state multiplicity
      gs_multi = int(line.split()[3])
    elif 'FINAL' in line and read_properties:
      # get (last) energy
      qc.etot = float(line.split()[4])
    elif 'TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS' in line and is_pop_ana == True and read_properties:
      # Read Mulliken and Lowdin Atomic Populations
      sec_flag = 'pop_info'
      pop_skip = 1
      is_pop_ana == False
      qc.pop_ana['Lowdin'] = []
      qc.pop_ana['Mulliken'] = []
    else:
      # Check if we are in a specific section
      if sec_flag == 'geo_info':
        if not geo_skip:
          if len(line) < 2:
            sec_flag = None
          else:
            qc.geo_info.append([thisline[0],atom_count+1,thisline[1]])
            qc.geo_spec.append([float(ii) for ii in thisline[2:]])
            atom_count += 1
        elif geo_skip:
          geo_skip -= 1

      elif sec_flag == 'ao_info':
        if not ao_skip:
          if ' TOTAL NUMBER OF BASIS SET SHELLS' in line:
            sec_flag = None
          else:
            if len(thisline) == 1:
              # Read atom type
              at_type = thisline[0]
              AO.append([])
              new_ao = False
            elif len(thisline) == 0 and new_ao == False:
              new_ao = True
            else:
              coeffs = [float(ii) for ii in thisline[3:]]
              if new_ao:
                ao_type = thisline[1].lower().replace('l','sp')
                for i_ao,t_ao in enumerate(ao_type):
                  AO[-1].append({'atom_type': at_type,
                                'type': t_ao,
                                'pnum': 1,
                                'coeffs': [[coeffs[0],coeffs[1+i_ao]]]})
                new_ao = False
              else:
                for i_ao in range(len(ao_type)):
                  AO[-1][-len(ao_type)+i_ao]['coeffs'].append([coeffs[0],
                                                              coeffs[1+i_ao]])
                  AO[-1][-len(ao_type)+i_ao]['pnum'] += 1
        elif ao_skip:
          ao_skip -= 1
      elif sec_flag == 'mo_info':
        if not mo_skip:
          if 'END OF' in line and 'CALCULATION' in line or '-----------' in line:
            sec_flag = None
            has_alpha = False
            has_beta = False
          else:
            if thisline == []:
              info_key = None
              init_mo = True
              try:
                int(flines[il+1].split()[0])
              except ValueError:
                sec_flag = None
                init_mo = False
            elif init_mo:
              init_len = len(thisline)
              lxlylz = []
              for ii in range(len(thisline)):
                if has_alpha == True or has_beta == True:
                  qc.mo_spec.append({'coeffs': [],
                                  'energy': 0.0,
                                  'occ_num': 0.0,
                                  'sym': '',
                                  'spin': ''
                                  })
                else:
                  qc.mo_spec.append({'coeffs': [],
                  'energy': 0.0,
                  'occ_num': 0.0,
                  'sym': ''
                  })
              init_mo = False
              info_key = 'energy'
            elif len(thisline) == init_len and info_key == 'energy':
              for ii in range(init_len,0,-1):
                qc.mo_spec[-ii]['energy'] = float(thisline[init_len-ii])
              info_key = 'symmetry'
            elif len(thisline) == init_len and info_key == 'symmetry':
              for ii in range(init_len,0,-1):
                len_mo += 1
                a = thisline[init_len-ii]
                if a not in sym.keys(): sym[a] = 1
                else: sym[a] = len_mo
                if has_alpha:
                  qc.mo_spec[-ii]['sym'] = '%d.%s_a' % (sym[a], thisline[init_len-ii])
                  qc.mo_spec[-ii]['spin'] = 'alpha'
                elif has_beta:
                  qc.mo_spec[-ii]['sym'] = '%d.%s_b' % (sym[a], thisline[init_len-ii])
                  qc.mo_spec[-ii]['spin'] = 'beta'
                else:
                  qc.mo_spec[-ii]['sym'] = '%d.%s' % (sym[a], thisline[init_len-ii])
              info_key = 'coeffs'
            elif thisline != [] and info_key == 'coeffs':
              lxlylz.append((line[11:17]))
              for ii, m in enumerate(re.finditer('-?\d+\.\d+', line[16:])):
                  qc.mo_spec[-init_len+ii]['coeffs'].append(float(m.group()))
        elif mo_skip:
          mo_skip -= 1
      elif sec_flag == 'ecp_info':
        if 'THE ECP RUN REMOVES' in line:
            sec_flag = None
        elif 'PARAMETERS FOR' in line:
          if line[17:25].split()[0] != ecp:
            ecp = line[17:25].split()[0]
            zcore = float(line[51:55].split()[0])
            ii_geo = int(line[35:41].split()[0])-1
            qc.geo_info[ii_geo][2] = str(float(qc.geo_info[ii_geo][2]) - zcore)
          else:
            ii_geo = int(line[35:41].split()[0])-1
            qc.geo_info[ii_geo][2] = str(float(qc.geo_info[ii_geo][2]) - zcore)

      elif sec_flag == 'dm_info':
        # instead of giving the output in a useful human and machine readable
        # way, gamess output syntax differs for transitions involving the
        # ground state compared to transitions between excited states...
        if 'GROUND STATE (SCF) DIPOLE=' in line:
          # ground state dipole is in debye...convert to atomic units
          for ii in range(3):
            qc.dipole_moments[0][0][ii] = float(thisline[ii+4])*0.393430307
        if 'EXPECTATION VALUE DIPOLE MOMENT FOR EXCITED STATE' in line:
          state = (int(line.replace('STATE', 'STATE ').split()[7]))
          dm_flag = 'state_info'
        if 'TRANSITION FROM THE GROUND STATE TO EXCITED STATE' in line:
          state = [0,
                    int(line.replace('STATE', 'STATE ').split()[8])]
          dm_flag = 'transition_info'
        if 'TRANSITION BETWEEN EXCITED STATES' in line:
          state = [int(thisline[4]),
                    int(line.replace('AND', 'AND ').split()[6])]
          dm_flag = 'transition_info'
        if 'NATURAL ORBITAL OCCUPATION NUMBERS FOR EXCITED STATE' in line:
          sec_flag = None
          dm_flag = None
        if dm_flag == 'state_info':
          if 'STATE MULTIPLICITY' in line:
            qc.states['multiplicity'][state] = int(line.split('=')[1])
          if 'STATE ENERGY' in line:
            qc.states['energy'][state] = float(line.split('=')[1])
          if 'STATE DIPOLE' and 'E*BOHR' in line:
            for ii in range(3):
              qc.dipole_moments[state][state][ii] = float(thisline[ii+3])
        elif dm_flag == 'transition_info':
          if 'TRANSITION DIPOLE' and 'E*BOHR' in line:
            for ii in range(3):
              qc.dipole_moments[state[0]][state[1]][ii] = float(thisline[ii+3])
              qc.dipole_moments[state[1]][state[0]][ii] = float(thisline[ii+3])
      elif sec_flag == 'pop_info':
        if not pop_skip:
          if  line == '\n':
            sec_flag = None
          else:
            qc.pop_ana = {}
            qc.pop_ana['Lowdin'].append(float(thisline[5]))
            qc.pop_ana['Mulliken'].append(float(thisline[3]))
        elif pop_skip:
          pop_skip -= 1

  # Check usage of same atomic basis sets
  basis_set = {}
  for ii in range(len(AO)):
    if not AO[ii][0]['atom_type'] in basis_set.keys():
      basis_set[AO[ii][0]['atom_type']] = AO[ii]
    else:
      for jj in range(len(AO[ii])):
        if AO[ii][jj]['coeffs'] !=  basis_set[AO[ii][0]['atom_type']][jj]['coeffs']:
          raise IOError('Different basis sets for the same atom.')
  # Numpy array
  for ii in basis_set.keys():
    for jj in range(len(basis_set[ii])):
      basis_set[ii][jj]['coeffs'] = numpy.array(basis_set[ii][jj]['coeffs'])

  for kk in range(len(qc.mo_spec)):
    qc.mo_spec[kk]['coeffs'] = numpy.array(qc.mo_spec[kk]['coeffs'])

  # Complement atomic basis sets
  for kk in range(len(qc.geo_info)):
    for ll in range(len(basis_set[qc.geo_info[kk][0]])):
      qc.ao_spec.append({'atom': qc.geo_info[kk][1]-1,
                         'type': basis_set[qc.geo_info[kk][0]][ll]['type'],
                         'pnum': basis_set[qc.geo_info[kk][0]][ll]['pnum'],
                         'coeffs': basis_set[qc.geo_info[kk][0]][ll]['coeffs'],
                         'lxlylz': None
                         })
  # Reconstruct exponents list for ao_spec
  count = 0
  for i,j in enumerate(qc.ao_spec):
    l = l_deg(lquant[j['type']])
    j['lxlylz'] = []
    for i in range(l):
      j['lxlylz'].append((lxlylz[count].lower().count('x'),
                          lxlylz[count].lower().count('y'),
                          lxlylz[count].lower().count('z')))
      count += 1
    j['lxlylz'] = numpy.array(j['lxlylz'],dtype=numpy.int64)

  if restricted:
    for ii in range(len(qc.mo_spec)):
      if occ[0] and occ[1]:
        qc.mo_spec[ii]['occ_num'] += 2.0
        occ[0] -= 1
        occ[1] -= 1
      if not occ[0] and occ[1]:
        qc.mo_spec[ii]['occ_num'] += 1.0
        occ[1] -= 1
      if not occ[1] and occ[0]:
        qc.mo_spec[ii]['occ_num'] += 1.0
        occ[0] -= 1

  if restricted == False:
    for ii in range(len(qc.mo_spec)):
      if qc.mo_spec[ii]['spin'] == 'alpha' and occ[0] > 0:
        qc.mo_spec[ii]['occ_num'] += 1.0
        occ[0] -= 1
        has_alpha = True
      elif qc.mo_spec[ii]['spin'] == 'beta' and occ[1] > 0:
        qc.mo_spec[ii]['occ_num'] += 1.0
        occ[1] -= 1
        has_beta = True

  if spin is not None:
    if restricted:
      raise IOError('The keyword `spin` is only supported for unrestricted calculations.')
    if spin != 'alpha' and spin != 'beta':
      raise IOError('`spin=%s` is not a valid option' % spin)
    elif spin == 'alpha' and has_alpha == True:
      display('Reading only molecular orbitals of spin alpha.')
    elif spin == 'beta' and has_beta == True:
      display('Reading only molecular orbitals of spin beta.')
    elif (not has_alpha) and (not has_beta):
      raise IOError(
            'No spin molecular orbitals available')
    elif ((spin == 'alpha' and not has_alpha) or
          (spin == 'beta' and not has_beta)):
      raise IOError('You requested `%s` orbitals, but None of them are present.'
                    % spin)

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

  # Convert geo_info and geo_spec to numpy.ndarrays
  qc.format_geo(is_angstrom=angstrom)

  qc.mo_spec.update()
  qc.ao_spec.update()
  return qc
