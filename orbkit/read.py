# -*- coding: iso-8859-1 -*-
'''Module for reading the output files of quantum chemical software. 
'''
'''
orbkit
Gunter Hermann, Vincent Pohl, and Axel Schild

Institut fuer Chemie und Biochemie, Freie Universitaet Berlin, 14195 Berlin, Germany

This file is part of orbkit.

orbkit is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or any later version.

orbkit is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with orbkit.  If not, see <http://www.gnu.org/licenses/>.
'''
import os
import copy
import numpy

from orbkit.core import l_deg, lquant, orbit, exp, exp_wfn
from orbkit.display import display
from orbkit.qcinfo import QCinfo, get_atom_symbol

def main_read(filename,itype='molden',all_mo=False,spin=None,cclib_parser=None,
              **kwargs):
  '''Calls the requested read function.
  
  **Parameters:**
  
    filename : str
      Specifies the filename for the input file.
    itype : str, choices={'molden', 'gamess', 'gaussian.log', 'gaussian.fchk', 'aomix', 'cclib'}
      Specifies the type of the input file.
    all_mo : bool, optional
      If True, all molecular orbitals are returned.
    spin : {None, 'alpha', or 'beta'}, optional
      If not None, returns exclusively 'alpha' or 'beta' molecular orbitals.
    cclib_parser : str
      If itype is 'cclib', specifies the cclib.parser.
  
  **Returns:**
  
    qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, mo_spec, etot :
          See :ref:`Central Variables` for details.
  
  **Note:**
  
    All additional keyword arguments are forwarded to the reading functions.
  '''
  display('Opened \n\t%s\n' % filename)
  
  # What kind of input file has to be read?
  reader = {'molden': read_molden, 
            'gamess': read_gamess, 
            'gaussian.log': read_gaussian_log, 
            'gaussian.fchk': read_gaussian_fchk,
            'aomix' : read_aomix,
            'cclib' : read_with_cclib,
            'wfn': read_wfn,
            'wfx': read_wfx}
  
  display('Loading %s file...' % itype)
  
  # Return required data
  return reader[itype](filename, all_mo=all_mo, spin=spin, 
                       cclib_parser=cclib_parser,**kwargs)
  # main_read 
  
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
  
  # Is this really a molden file? 
  if not '[Molden Format]\n' in flines:
    display('The input file %s is no valid molden file!\n\nIt does'  % filename+
          ' not contain the keyword: [Molden Format]\n')
    raise IOError('Not a valid input file')
  
  def check_sel(count,i,interactive=False):
    if count == 0:
      raise IndexError
    elif count == 1:
      return 0
    message = '\tPlease give an integer from 0 to %d: ' % (count-1)
    
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
  count = 0
  # Go through the file line by line 
  for il in range(len(flines)):
    line = flines[il]            # The current line as string
    
    # Check the file for keywords 
    if '[Molden Format]' in line:
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

  # Declare synonyms for molden keywords 
  synonyms = {'Sym': 'sym',
              'Ene': 'energy',
              'Occup': 'occ_num',
              'Spin': 'spin'
             }
  MO_keys = synonyms.keys()
  
  count = 0
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
              basis_count += l_deg(lquant[i_ao])
              qc.ao_spec.append({'atom': at_num,
                              'type': i_ao,
                              'pnum': pnum,
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
                  info[1] = info[1].replace(a, '%s.' % a)
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
    
  # Convert geo_info and geo_spec to numpy.ndarrays
  qc.format_geo()
    
  return qc
  # read_molden 

def read_gamess(filename, all_mo=False, spin=None, read_properties=False,
                **kwargs):
  '''Reads all information desired from a Gamess-US output file.
  
  **Parameters:**
  
    filename : str
      Specifies the filename for the input file.
    all_mo :   bool, optional
      If True, all molecular orbitals are returned.
  
  **Returns:**
  
    qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, mo_spec, etot :
        See :ref:`Central Variables` for details.
  '''
  
  fid    = open(filename,'r')      # Open the file
  flines = fid.readlines()         # Read the WHOLE file into RAM
  fid.close()                      # Close the filename
  
  # Initialize the variables 
  qc = QCinfo()
  has_alpha = False                  # Flag for alpha electron set
  has_beta = False                   # Flag for beta electron set
  restricted = True                  # Flag for restricted calculation
  sec_flag = None                    # A Flag specifying the current section
  is_pop_ana = True                  # Flag for population analysis for ground state
  keyword = [' ATOM      ATOMIC                      COORDINATES','\n'] 
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
        mokey = 'MOLECULAR ORBITALS'
      
    elif keyword[0] in line and keyword[1] in flines[il-1]:
      # The section containing information about 
      # the molecular geometry begins 
      sec_flag = 'geo_info'
      atom_count = 0  # Counter for Atoms
      if '(BOHR)' in line:
        # The length are given in Bohr radii
        aa_to_au = 1.0
      else:
        # The length are given in Angstroem 
        # and have to be converted to Bohr radii
        aa_to_au = 1/0.52917720859
        
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
      
    elif mokey in line:
      # The section containing information about 
      # the molecular orbitals begins 
      sec_flag = 'mo_info'
      mo_skip = 1
      len_mo = 0                      # Number of MOs 
      init_mo = False                 # Initialize new MO section
      info_key = None                 # A Flag specifying the energy and symmetry section
      exp_list = []
      if 'ALPHA' in line:
        has_alpha = True
        mo_skip = 0
      elif 'BETA' in line:
        has_beta = True
        has_alpha = False
        mo_skip = 0
    
    elif 'NATURAL ORBITALS' in line:
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
          if line == '\n':
            sec_flag = None
          else:
            qc.geo_info.append([thisline[0],atom_count+1,thisline[1]])
            qc.geo_spec.append([float(ii)*aa_to_au for ii in thisline[2:]])
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
              exp_list = []
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
              exp_list.append((line[11:17]))
              for ii in range(init_len,0,-1):
                qc.mo_spec[-ii]['coeffs'].append(float(line[16:].split()[init_len-ii]))
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
  for ii in basis_set.iterkeys():
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
                         'exp_list': None
                         })
  # Reconstruct exponents list for ao_spec
  count = 0
  for i,j in enumerate(qc.ao_spec):    
    l = l_deg(lquant[j['type']])
    j['exp_list'] = []
    for i in range(l):
      j['exp_list'].append((exp_list[count].lower().count('x'),
                            exp_list[count].lower().count('y'),
                            exp_list[count].lower().count('z')))
      count += 1
    j['exp_list'] = numpy.array(j['exp_list'],dtype=numpy.int64)
  
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
  qc.format_geo()

  return qc
  # read_gamess 

def read_gaussian_fchk(filename, all_mo=False, spin=None, **kwargs):
  '''Reads all information desired from a Gaussian FChk file. 

  **Parameters:**
  
    filename : str
      Specifies the filename for the input file.
    all_mo : bool, optional
      If True, all molecular orbitals are returned.
  
  **Returns:**
  
    qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, mo_spec, etot :
        See :ref:`Central Variables` for details.
  '''
  
  fid    = open(filename,'r')   # Open the file
  flines = fid.readlines()      # Read the WHOLE file into RAM
  fid.close()                   # Close the file
  
  # Is this an unrestricted calculation?
  has_beta = False
  is_6D = False
  is_10F = False
  for line in flines:
    if 'beta mo coefficients' in line.lower():
      has_beta = True
    if 'Pure/Cartesian d shells' in line:
      is_6D = int(line.split()[-1]) == 1
    if 'Pure/Cartesian f shells' in line:
      is_10F = int(line.split()[-1]) == 1
  
  cartesian_basis = (is_6D and is_10F)
  if ((not is_6D) and is_10F) or (is_6D and (not is_10F)):
    raise IOError('Please apply a Spherical Harmonics (5D, 7F) or '+
                          'a Cartesian Gaussian Basis Set (6D, 10F)!')
  
  if spin is not None:
    if spin != 'alpha' and spin != 'beta':
      raise IOError('`spin=%s` is not a valid option' % spin)
    elif has_beta:
      display('Reading only molecular orbitals of spin %s.' % spin)
    else:
      raise IOError('The keyword `spin` is only supported for unrestricted calculations.')
  restricted = (not has_beta)
  
  sec_flag = None
  
  el_num = [0,0]
  mo_i0 = {'alpha': 0, 'beta': 0}
  what = 'alpha'
  index = 0
  at_num = 0
  
  ao_num = 0 
  switch = 0
  qc = QCinfo()
  qc.geo_info = [[],[],[]]
  if not cartesian_basis:
    qc.ao_spherical = []
  # Set a counter for the AOs 
  basis_count = 0
  
  # Go through the file line by line 
  for il in range(len(flines)):
    line = flines[il]         # The current line as string
    thisline = line.split()   # The current line split into segments
    
    # Check the file for keywords 
    if 'Number of alpha electrons' in line:
      el_num[0] = int(thisline[5]) 
    elif 'Number of beta electrons' in line:
      el_num[1] = int(thisline[5])
    elif 'Number of basis functions' in line:
      basis_number = int(thisline[5])
    elif 'Atomic numbers'  in line:
      sec_flag = 'geo_info'
      index = 0
      at_num = int(thisline[-1])
      count = 0
      qc.geo_info[1] = list(range(1,at_num+1))
    elif 'Nuclear charges' in line:
      sec_flag = 'geo_info'
      index = 2
      at_num = int(thisline[-1])
      count = 0
    elif 'Total Energy' in line:
      qc.etot = float(thisline[3])
    elif 'Current cartesian coordinates' in line:
      at_num = int(thisline[-1])/3
      sec_flag = 'geo_pos'
      qc.geo_spec = []
      count = 0
      xyz = []
    elif 'Shell types' in line:
      sec_flag = 'ao_info'
      index = 'type'
      ao_num = int(thisline[-1])
      count = 0
      if qc.ao_spec == []:
        for ii in range(ao_num):
          qc.ao_spec.append({})
    elif 'Number of primitives per shell' in line:
      sec_flag = 'ao_info'
      index = 'pnum'
      ao_num = int(thisline[-1])
      count = 0
      if qc.ao_spec == []:
        for ii in range(ao_num):
          qc.ao_spec.append({})
    elif 'Shell to atom map' in line:
      sec_flag = 'ao_info'
      index = 'atom'
      ao_num = int(thisline[-1])
      count = 0
      if qc.ao_spec == []:
        for ii in range(ao_num):
          qc.ao_spec.append({})
    elif 'Primitive exponents' in line:
      sec_flag = 'ao_coeffs'
      ao_num = int(thisline[-1])
      count = 0
      switch = 0
      index = 0
      if qc.ao_spec == []:
        raise IOError('Shell types need to be defined before the AO exponents!')
      if not 'coeffs' in qc.ao_spec[0].keys():
        for ii in range(len(qc.ao_spec)):
          pnum = qc.ao_spec[ii]['pnum']
          qc.ao_spec[ii]['coeffs'] = numpy.zeros((pnum, 2))
    elif 'Contraction coefficients' in line:
      sec_flag = 'ao_coeffs'
      ao_num = int(thisline[-1])
      count = 0
      switch = 1
      index = 0
      if qc.ao_spec == []:
        raise IOError('Shell types need to be defined before the AO exponents!')
      if not 'coeffs' in qc.ao_spec[0].keys():
        for ii in range(len(qc.ao_spec)):
          pnum = qc.ao_spec[ii]['pnum']
          qc.ao_spec[ii]['coeffs'] = numpy.zeros((pnum, 2))
    elif 'Orbital Energies' in line:
      sec_flag = 'mo_eorb'
      mo_num = int(thisline[-1])      
      mo_i0[thisline[0].lower()] = len(qc.mo_spec)#index
      if restricted:
        if el_num[0] == el_num[1]:
          i = el_num[0]
          occ = 2
        else:
          i = el_num[0 if 'Alpha' in line else 1]
          occ = 1
      else:
        i = el_num[0 if 'Alpha' in line else 1]
        occ = 1      
      for ii in range(mo_num):
        qc.mo_spec.append({'coeffs': numpy.zeros(basis_count),
                        'energy': 0.0,
                        'occ_num': float(occ if ii < i else 0),
                        'sym': '%i.1' % (ii+1),
                        'spin':thisline[0].lower()
                        })
    elif 'MO coefficients' in line:
      sec_flag = 'mo_coeffs'
      count = 0
      index = 0
      mo_num = int(thisline[-1])
      what = thisline[0].lower()
    else:
      # Check if we are in a specific section 
      if sec_flag == 'geo_info':
        for ii in thisline:
          qc.geo_info[index].append(ii)
          count += 1
          if count == at_num:
            sec_flag = None
      elif sec_flag == 'geo_pos':
        for ii in thisline:
          xyz.append(float(ii))
          if len(xyz) == 3:
            qc.geo_spec.append(xyz)
            xyz = []
            count += 1
            if count == at_num:
              sec_flag = None
      elif sec_flag == 'ao_info':
        for ii in thisline:
          ii = int(ii)
          if index is 'type':
            ii = orbit[abs(ii)]
            l = lquant[ii]
            basis_count += l_deg(l,cartesian_basis=cartesian_basis)
            if not cartesian_basis:
              for m in (range(0,l+1) if l != 1 else [1,0]):
                qc.ao_spherical.append([count,(l,m)])
                if m != 0:
                  qc.ao_spherical.append([count,(l,-m)])
          elif index is 'atom':
            ii -= 1
          qc.ao_spec[count][index] = ii
          count += 1
          if count == ao_num:
            sec_flag = None
      elif sec_flag == 'ao_coeffs':
        for ii in thisline:
          qc.ao_spec[index]['coeffs'][count,switch] = float(ii)
          count += 1
          ao_num -= 1
          if count == qc.ao_spec[index]['pnum']:
            index += 1
            count = 0
        if not ao_num:
          sec_flag = None
      elif sec_flag == 'mo_eorb':
        for ii in thisline:
          qc.mo_spec[count]['energy'] = float(ii)
          count += 1
          if index != 0 and not count % basis_count:
            sec_flag = None
      elif sec_flag == 'mo_coeffs':
        for ii in thisline:    
          qc.mo_spec[mo_i0[what]+index]['coeffs'][count] = float(ii)
          count += 1
          if count == basis_count:
            count = 0
            index += 1
          if index != 0 and not index % basis_count:
            sec_flag = None
  
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
  
  if restricted:
    # Closed shell calculation
    for mo in qc.mo_spec:
      del mo['spin']
  else:
    # Rename MOs according to spin
    for mo in qc.mo_spec:
      mo['sym'] += '_%s' % mo['spin'][0]
  
  # Check for natural orbital occupations
  energy_sum = sum([abs(i['energy']) for i in qc.mo_spec])
  if energy_sum < 0.0000001:
    display('Attention!\n\tThis FChk file contains natural orbitals. '+
            '(There are no energy eigenvalues.)\n\t' + 
            'In this case, Gaussian does not print the respective natural' +
            'occupation numbers!' )
  
  qc.geo_info = numpy.array(qc.geo_info).T
  # Convert geo_info and geo_spec to numpy.ndarrays
  qc.format_geo()
  
  return qc
  # read_gaussian_fchk 

def read_gaussian_log(filename,all_mo=False,spin=None,orientation='standard',
                      i_link=-1,i_geo=-1,i_ao=-1,i_mo=-1,interactive=True,
                      **kwargs):
  '''Reads all information desired from a Gaussian .log file.

  **Parameters:**
  
    filename : str
      Specifies the filename for the input file.
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
  
    qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, ao_spherical, mo_spec, etot :
        See :ref:`Central Variables` for details.

.. [#first] Attention: The MOs in the output are only valid for the standard orientation!

  '''
  
  fid    = open(filename,'r')      # Open the file
  flines = fid.readlines()         # Read the WHOLE file into RAM
  fid.close()                      # Close the file
  
  
  # Search the file the specific sections
  count = {'link': 0, 'geometry': 0, 'geometry_input': 0, 'atomic orbitals': 0, 
           'molecular orbitals': [], 'state': []}

  def check_sel(count,i,interactive=False):
    if count == 0:
      raise IndexError
    elif count == 1:
      return 0
    message = '\tPlease give an integer from 0 to %d: ' % (count-1)
    
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
        elif 'AO basis set' in line:
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
  
  aa_to_au = 1/0.52917720859  # conversion factor for Angstroem to bohr radii

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
      elif 'AO basis set' in line:
        # The section containing information about 
        # the atomic orbitals begins
        if i_ao == c_ao:
          qc.ao_spec = []
          if not cartesian_basis:
            qc.ao_spherical = []
          sec_flag = 'ao_info'
        c_ao += 1
        basis_count = 0
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
          qc.mo_spec = []
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
            qc.geo_spec.append([aa_to_au*float(ij) for ij in thisline[3:]])
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
          else:
            # Append the AO coefficients 
            coeffs = numpy.array(line.replace('D','e').split(), dtype=numpy.float64)
            for i_ao in range(len(ao_type)):
              qc.ao_spec[-len(ao_type)+i_ao]['coeffs'][ao_num,:] = [coeffs[0],
                                                                coeffs[1+i_ao]]
            ao_num += 1
        if sec_flag == 'mo_sym':
          if 'The electronic state is' in line:
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
            coeffs = line[21:].replace('-',' -').split()
            if not cartesian_basis and offset == 0:
              if old_ao != line[:14].split()[-1]:
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
              qc.ao_spherical.append([i,(l,m)])
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
  qc.format_geo()
  return qc
  # read_gaussian_log 


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

def read_aomix(filename, all_mo=False, spin=None, i_md=-1, interactive=True,
               created_by_tmol=True, **kwargs):
  '''Reads all information desired from a aomix file.
  
  **Parameters:**
  
    filename : str
      Specifies the filename for the input file.
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
  
  fid    = open(filename,'r')      # Open the file
  flines = fid.readlines()         # Read the WHOLE file into RAM
  fid.close()                      # Close the file
  
  # Is this really a aomix file? 
  if not '[AOMix Format]\n' in flines:
    display('The input file %s is no valid aomix file!\n\nIt does'  % filename+
          ' not contain the keyword: [AOMix Format]\n')
    raise IOError('Not a valid input file')
  
  def check_sel(count,i,interactive=False):
    if count == 0:
      raise IndexError
    elif count == 1:
      return 0
    message = '\tPlease give an integer from 0 to %d: ' % (count-1)
    
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
  count = 0
  # Go through the file line by line 
  for il in range(len(flines)):
    line = flines[il]            # The current line as string
    
    # Check the file for keywords 
    if '[AOMix Format]' in line:
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
    display('The input file %s is no valid aomix file!\n\nIt does' % filename +
          ' not contain the keyword: [AOMix Format]\n')
    raise IOError('Not a valid input file')
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
  
  exp_list = []
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
      else:
        # Check if we are in a specific section 
        if sec_flag == 'geo_info':
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
              basis_count += l_deg(lquant[i_ao])
              qc.ao_spec.append({'atom': at_num,
                              'type': i_ao,
                              'pnum': pnum,
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
                  info[1] = info[1].replace(a, '%s.' % a)
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
                  exp_list.append(thisline[-2])
              except ValueError:
                # If it cannot be converted print error message 
                display('Error in coefficient %d of MO %s!' % (index, 
                  qc.mo_spec[-1]['sym']) + '\nSetting this coefficient to zero...')

  # Check usage of same atomic basis sets
  for ii in range(len(exp_list)):
    s = exp_list[ii]
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
    exp_list[ii] = exp

  count = 0
  for i,j in enumerate(qc.ao_spec):    
    l = l_deg(lquant[j['type']])
    j['exp_list'] = []
    for i in range(l):
      j['exp_list'].append((exp_list[count][0],
                      exp_list[count][1],
                      exp_list[count][2]))
      count += 1
    j['exp_list'] = numpy.array(j['exp_list'],dtype=numpy.int64)
  
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
  qc.format_geo()
  
  if is_tmol_cart and created_by_tmol:
    display('\nFound a Cartesian basis set in the AOMix file.')
    display('We assume that this file has been created by Turbomole.')
    display('Applying a conversion to the molecular orbital coefficients, ')
    display('in order to get normalized orbitals.')
    
    # Convert MO coefficients
    from orbkit import cSupportCode
    from orbkit.analytical_integrals import create_mo_coeff, get_lxlylz
    def dfact(n):
      if n <= 0:
        return 1;
      else:
        return n * dfact(n-2)

    mo = create_mo_coeff(qc.mo_spec)
    for i,j in enumerate(get_lxlylz(qc.ao_spec)):
      norm = (dfact(2*j[0] - 1) * dfact(2*j[1] - 1) * dfact(2*j[2] - 1))
      j = sum(j)
      if j >1: 
        mo[:,i] *= numpy.sqrt(norm)   
    for ii in range(len(qc.mo_spec)):
      qc.mo_spec[ii]['coeffs'] = mo[ii]
  
  return qc
  # read_aomix
  

def read_wfx(filename, all_mo=False, spin=None, **kwargs):
  '''Reads all information desired from a wfn file.
  
  **Parameters:**
  
    filename : str
      Specifies the filename for the input file.
    all_mo : bool, optional
      If True, all molecular orbitals are returned.
    spin : {None, 'alpha', or 'beta'}, optional
      If not None, returns exclusively 'alpha' or 'beta' molecular orbitals.
  
  **Returns:**
  
    qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, mo_spec, etot :
        See :ref:`Central Variables` for details.
  '''

  # Initialize the variables 
  qc = QCinfo()
  exp_list = []
  for j in exp_wfn:         
    exp_list.extend(j)
  exp_list = numpy.array(exp_list,dtype=numpy.int64)
  
  fid    = open(filename,'r')      # Open the file
  flines = fid.readlines()         # Read the WHOLE file into RAM
  fid.close()                      # Close the file
  
  is_valid = False
  for il in range(len(flines)):
    if '<Keywords>' in flines[il] and 'GTO' in flines[il+1]:
      is_valid = True
  
  if not is_valid: raise IOError('No valid .wfx file!\nMissing:\n' + 
                                 '<Keywords>\n  GTO\n</Keywords>')
  
  sec_flag = None                 # A Flag specifying the current section
  at_num = None
  mo_num = None
  ao_num = None
  restricted = True
  count = 0
  
  # Go through the file line by line 
  for il in range(len(flines)):
    line = flines[il]             # The current line as string
    
    if '<Number of Nuclei>' in line:
      at_num = int(flines[il+1])
      qc.geo_info = [[None, i+1, None] for i in range(at_num)]
      qc.geo_spec = []
    elif '<Nuclear Names>' in line:
      if not at_num: raise IOError('`<Number of Nuclei>` has to be found ' +
                                   'before `<Nuclear Names>`.')
      for i in range(at_num):
        qc.geo_info[i][0] = flines[il+i+1].replace(' ','').replace('\n','')
    elif '<Atomic Numbers>' in line:
      if not at_num: raise IOError('`<Number of Nuclei>` has to be found ' +
                                   'before `<Atomic Numbers>`.')
      for i in range(at_num):
        qc.geo_info[i][2] = flines[il+i+1].replace(' ','').replace('\n','')
    elif '<Nuclear Cartesian Coordinates>' in line:
      if not at_num: raise IOError('`<Number of Nuclei>` has to be found ' +
                                   'before `<Nuclear Cartesian Coordinates>`.')
      for i in range(at_num):
        qc.geo_spec.append(flines[il+i+1].split())
    elif '<Number of Primitives>' in line:
      ao_num = int(flines[il+1])
      qc.ao_spec = [{'atom': None,
                    'pnum': -1,
                    'coeffs': None,
                    'exp_list': None,
                    } for i in range(ao_num)]
    elif '<Primitive Centers>' in line:
      sec_flag = 'ao_center'
      count = 0
    elif '<Primitive Types>' in line:
      sec_flag = 'ao_type'
      count = 0
    elif '<Primitive Exponents>' in line:
      sec_flag = 'ao_exp'
      count = 0
    elif '<Number of Occupied Molecular Orbitals>' in line:    
      mo_num = int(flines[il+1])
      qc.mo_spec = [{'coeffs': numpy.zeros(ao_num),
                     'energy': None,
                     'occ_num': None,
                     'spin': None,
                     'sym': '%s.1' % (i+1)
                    } for i in range(mo_num)]
    elif '<Molecular Orbital Occupation Numbers>' in line:
      for i in range(mo_num):
        qc.mo_spec[i]['occ_num'] = float(flines[il+1+i])
    elif '<Molecular Orbital Energies>' in line:
      for i in range(mo_num):
        qc.mo_spec[i]['energy'] = float(flines[il+1+i])
    elif '<Molecular Orbital Spin Types>' in line:
      for i in range(mo_num):
        qc.mo_spec[i]['spin'] = (flines[il+1+i].replace(' ','').replace('\n','')
                                 ).replace('and','_').lower()
        restricted = restricted and ('_' in qc.mo_spec[i]['spin'])
    elif '<MO Number>' in line:
      index = int(flines[il+1])-1
      for i in range(ao_num):
        qc.mo_spec[index]['coeffs'][i] = float(flines[il+3+i])
    elif '</' in line:
      sec_flag = None
    elif sec_flag is not None:
      if sec_flag == 'ao_center':
        for i in line.split():
          qc.ao_spec[count]['atom'] = int(i)-1
          count += 1
      if sec_flag == 'ao_type':
        for i in line.split():
          qc.ao_spec[count]['exp_list'] = exp_list[int(i)-1][numpy.newaxis]
          count += 1
      if sec_flag == 'ao_exp':
        for i in line.split():
          qc.ao_spec[count]['coeffs'] = numpy.array([[float(i),1.0]])
          count += 1
      
  
  has_alpha = any([i['spin'] == 'alpha' for i in qc.mo_spec])
  has_beta  = any([i['spin'] == 'beta'  for i in qc.mo_spec])
  
  spin_check(spin,restricted,has_alpha,has_beta)
  qc.select_spin(restricted,spin=spin)
  
  # Remove numbers from atom names
  for i in qc.geo_info:
    i[0] = ''.join([k for k in i[0] if not k.isdigit()])
  # Convert geo_info and geo_spec to numpy.ndarrays
  qc.format_geo()
  return qc
      
  
def read_wfn(filename, all_mo=False, spin=None, **kwargs):
  '''Reads all information desired from a wfn file.
  
  **Parameters:**
  
    filename : str
      Specifies the filename for the input file.
    all_mo : bool, optional
      If True, all molecular orbitals are returned.
  
  **Returns:**
  
    qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, mo_spec, etot :
        See :ref:`Central Variables` for details.
  '''
  if spin is not None:
    raise IOError('The option `spin` is not supported for the `.wfn` reader.')
  
  # Initialize the variables 
  qc = QCinfo()
  sec_flag = None                 # A Flag specifying the current section
  is_wfn = False                  # Check type of file
  ao_num = 0                      # Number of AO
  mo_num = 0                      # Number of MO
  at_num = 0                      # Number of atoms
  c_type = 0                      # Counting variable for AO type
  c_exp = 0                       # Counting variable for AO exponents
  exp_list = []
  for j in exp_wfn:         
    exp_list.extend(j)
  exp_list = numpy.array(exp_list,dtype=numpy.int64)
  
  with open(filename) as fileobject:
    for line in fileobject:
      thisline = line.split()      # The current line split into segments
      # Check the file for keywords 
      if 'GAUSSIAN' in line or 'GTO' in line:
        if len(thisline) == 8:
          mo_num = int(thisline[1])
          ao_num = int(thisline[4])
          at_num = int(thisline[6])
          sec_flag = 'geo_info'
      elif 'CENTRE ASSIGNMENTS' in line:
        thisline = line[20:].split()
        for i in range(len(thisline)):
          qc.ao_spec.append({'atom': int(thisline[i])-1,
                'pnum': -1,
                'coeffs': None,
                'exp_list': None,
                })
      elif 'TYPE ASSIGNMENTS' in line:
        thisline = line[18:].split()
        for i in range(len(thisline)):
          qc.ao_spec[c_type]['exp_list'] = exp_list[int(thisline[i])-1][numpy.newaxis]
          c_type += 1
      elif 'EXPONENTS' in line:
        thisline = line.replace('EXPONENTS','').replace('D','E').split()
        for i in thisline:
          qc.ao_spec[c_exp]['coeffs'] = numpy.array([[float(i),1.0]])
          c_exp += 1
      elif 'MO' in line and 'OCC NO =' in line and 'ORB. ENERGY =' in line:
        qc.mo_spec.append({'coeffs': numpy.zeros(ao_num),
                'energy': float(line[25:].split()[7]),
                'occ_num': float(line[25:].split()[3]),
                'sym': '%s.1' % thisline[1]
                })
        sec_flag = 'mo_info'
        c_mo = 0                        # Counting variable for MOs
      else:
        if sec_flag == 'geo_info':
          if not at_num:
            sec_flag = None
          elif at_num:
            qc.geo_info.append([thisline[0],thisline[-7][:-1],thisline[-1]]) 
            qc.geo_spec.append([float(ii) for ii in thisline[-6:-3]])
            at_num -= 1
        elif sec_flag == 'mo_info':
          for i in thisline:
            if (c_mo) < ao_num:
              qc.mo_spec[-1]['coeffs'][c_mo] = numpy.array(
                                                      float(i.replace('D','E')))
              c_mo += 1
            if (c_mo) == ao_num:
              sec_flag = None
  
  # Remove numbers from atom names
  for i in qc.geo_info:
    i[0] = ''.join([k for k in i[0] if not k.isdigit()])
  # Convert geo_info and geo_spec to numpy.ndarrays
  qc.format_geo()
  
  return qc

def read_with_cclib(filename, cclib_parser=None, all_mo=False, spin=None, 
                    **kwargs):
  '''Reads all information desired using cclib.
  
  **Parameters:**
  
    filename : str
      Specifies the filename for the input file.
    cclib_parser : str
      If itype is 'cclib', specifies the cclib.parser.
    all_mo : bool, optional
      If True, all molecular orbitals are returned.
    spin : {None, 'alpha', or 'beta'}, optional
      If not None, returns exclusively 'alpha' or 'beta' molecular orbitals.
  
  **Returns:**
  
    qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, mo_spec, etot :
        See :ref:`Central Variables` for details.
  '''
  
  if not isinstance(cclib_parser,str):
    raise IOError('cclib requires the specification of parser, e.g., ' + 
                  'cclib_parser="Gaussian".')
  
  if cclib_parser == 'Molpro':
    display('\nThe Molpro basis set is not properly read by the cclib parser.')
    display('Please create a molden file with Molpro, i.e., ' + 
            '\n\tput,molden,output.molden,NEW;\n')
  
  exec('from cclib.parser import %s as parser' % cclib_parser)
  ccData = parser(filename).parse()
  return convert_cclib(ccData, all_mo=all_mo, spin=spin)

def convert_cclib(ccData, all_mo=False, spin=None):
  '''Converts a ccData class created by cclib to an instance of
  orbkit's QCinfo class.

  **Parameters:**
  
    ccData : class
      Contains the input data created by cclib.
    all_mo : bool, optional
      If True, all molecular orbitals are returned.
    spin : {None, 'alpha', or 'beta'}, optional
      If not None, returns exclusively 'alpha' or 'beta' molecular orbitals.


  **Returns:**

    qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, mo_spec, etot :
          See :ref:`Central Variables` for details.
  '''
  aa_to_au = 1/0.52917720859
  # Initialize the variables 
  qc = QCinfo()
  
  # Converting all information concerning atoms and geometry
  qc.geo_spec = ccData.atomcoords[0] * aa_to_au
  for ii in range(ccData.natom):
    symbol = get_atom_symbol(atom=ccData.atomnos[ii])
    qc.geo_info.append([symbol,str(ii+1),str(ccData.atomnos[ii])])
  
  # Convert geo_info and geo_spec to numpy.ndarrays
  qc.format_geo()
  
  # Converting all information about atomic basis set

  for ii in range(ccData.natom):
    for jj in range(len(ccData.gbasis[ii])):
      pnum = len(ccData.gbasis[ii][jj][1])
      qc.ao_spec.append({'atom': ii,
                  'type': str(ccData.gbasis[ii][jj][0]).lower(),
                  'pnum':  pnum,
                  'coeffs': numpy.zeros((pnum, 2))
                  })
      for kk in range(pnum):
        qc.ao_spec[-1]['coeffs'][kk][0] = ccData.gbasis[ii][jj][1][kk][0]
        qc.ao_spec[-1]['coeffs'][kk][1] = ccData.gbasis[ii][jj][1][kk][1]
        
  # Reconstruct exponents list for ao_spec
  cartesian_basis = True
  for i in ccData.aonames:
    if '+' in i or '-' in i:
      cartesian_basis = False

  if not cartesian_basis:
      qc.ao_spherical = []
  
  count = 0
  for i,j in enumerate(qc.ao_spec):
    l = l_deg(lquant[j['type']],cartesian_basis=cartesian_basis)
    if cartesian_basis:
      j['exp_list'] = []
      
    for ll in range(l):
      if cartesian_basis:
        j['exp_list'].append((ccData.aonames[count].lower().count('x'),
                              ccData.aonames[count].lower().count('y'),
                              ccData.aonames[count].lower().count('z')))
      else:
        m = ccData.aonames[count].lower().split('_')[-1]
        m = m.replace('+',' +').replace('-',' -').replace('s','s 0').split(' ') 
        p = 'yzx'.find(m[0][-1])
        if p != -1:
          m = p - 1
        else:
          m = int(m[-1])
        qc.ao_spherical.append([i,(lquant[j['type']],m)])
      count += 1
  
  # Converting all information about molecular orbitals
  ele_num = numpy.sum(ccData.atomnos) - numpy.sum(ccData.coreelectrons) - ccData.charge
  ue = (ccData.mult-1)
  
  # Check for natural orbitals and occupation numbers
  is_natorb = False
  if hasattr(ccData,'nocoeffs'):
    if not hasattr(ccData,'nooccnos'):
      raise IOError('There are natural orbital coefficients (`nocoeffs`) in the cclib' + 
                    ' ccData, but no natural occupation numbers (`nooccnos`)!')
    is_natorb = True 
  
  restricted = (len(ccData.mosyms) == 1)
  if spin is not None:
    if spin != 'alpha' and spin != 'beta':
      raise IOError('`spin=%s` is not a valid option' % spin)
    elif restricted:
      raise IOError('The keyword `spin` is only supported for unrestricted calculations.')
    else:
      display('Converting only molecular orbitals of spin %s.' % spin)
  
  sym = {}
  if len(ccData.mosyms) == 1:
    add = ['']
    orb_sym = [None]
  else:
    add = ['_a','_b']      
    orb_sym = ['alpha','beta']
  
  for ii in range(ccData.nmo):    
    for i,j in enumerate(add):
      a = '%s%s' % (ccData.mosyms[i][ii],j)
      if a not in sym.keys(): sym[a] = 1
      else: sym[a] += 1
      if is_natorb:
        occ_num = ccData.nooccnos[ii]
      elif not restricted:
        occ_num = 1.0 if ii <= ccData.homos[i] else 0.0
      elif ele_num > ue:
        occ_num = 2.0
        ele_num -= 2.0
      elif ele_num > 0.0 and ele_num <= ue: 
        occ_num = 1.0
        ele_num -= 1.0
        ue -= 1.0
      else:
        occ_num = 0.0
      qc.mo_spec.append({'coeffs': (ccData.nocoeffs if is_natorb else ccData.mocoeffs[i])[ii],
              'energy': 0.0 if is_natorb else ccData.moenergies[i][ii],
              'occ_num': occ_num,
              'sym': '%d.%s' %(sym[a],a)
              })
      if orb_sym[i] is not None:
        qc.mo_spec[-1]['spin'] = orb_sym[i]
        if spin is not None and spin != orb_sym[i]:
          del qc.mo_spec[-1]
    
  # Are all MOs requested for the calculation? 
  if not all_mo:
    for i in range(len(qc.mo_spec))[::-1]:
      if qc.mo_spec[i]['occ_num'] < 0.0000001:
        del qc.mo_spec[i]

  return qc


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
  
  def get_selection(selected_mo):
    mo_occup = numpy.array([i['occ_num'] for i in mo_spec])
    homo = (mo_occup>0.).nonzero()[0][-1]   + 1 # molden numbering
    lumo = (mo_occup>0.).nonzero()[0][-1]+1 + 1 # molden numbering
    sel = []
    for i in selected_mo:
      i = i.lower().split(':')
      if len(i) == 1:
        sel.append(eval(i[0]))
      else:
        i = range(*[eval(j) for j in i])
        sel.extend(i)
    return sel
  
  if isinstance(fid_mo_list,str) and fid_mo_list.lower() == 'all_mo':
    selected_mo = numpy.array(numpy.arange(len(mo_spec))+1, dtype=numpy.str)
    mo_in_file = [selected_mo]
    selected_mo_spec = mo_spec
    selected_mo_ii = numpy.array([i['sym'] for i in selected_mo_spec])
  else:
    if isinstance(fid_mo_list,str) and not os.path.exists(fid_mo_list):
      if ',' in fid_mo_list:
        fid_mo_list = fid_mo_list.split(',')
      else:
        fid_mo_list = [fid_mo_list]
    if isinstance(fid_mo_list, list):
      for i in fid_mo_list:
        if not isinstance(i, list):
          i = i.split(',') if isinstance(i,str) else [i]
        selected_mo.extend(map(str,i))
        mo_in_file.append(map(str,i))
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
      display('\tLine %d: %s' % (i+1,str(j)))
    
    # Check if the molecular orbitals are specified by symmetry 
    # (e.g. 1.1 in MOLPRO nomenclature) or 
    # by the number in the input file (e.g. 1)
    
    try: # Try to convert selections into integer
      for i in selected_mo:
        if isinstance(i,int):
          continue
        i = i.replace('homo','1').replace('lumo','2')
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
        selected_mo = map(int, selected_mo)            
        selected_mo.sort()
      selected_mo = map(str, selected_mo)
      what = lambda x,y: str(x+1)
      selected_mo_spec,selected_mo_ii = assign_selected_mo(selected_mo,
                                                           mo_spec,
                                                           strict=strict,
                                                           what=what)
      selected_mo = selected_mo_ii
      for i in range(len(mo_in_file)):
        mo_in_file[i] = map(str, get_selection(mo_in_file[i]))
    
    # Print some information
    display('\nThe following orbitals will be considered...')
    for i,j in enumerate(mo_in_file):
      display('\tLine %d: %s' % (i+1,str(j)))
  
  display('')
  return {'mo': selected_mo, 'mo_ii': selected_mo_ii,
          'mo_spec': selected_mo_spec, 
          'mo_in_file': mo_in_file, 'sym_select': sym_select}
