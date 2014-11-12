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

import copy
import numpy

from orbkit.core import l_deg, lquant, orbit, exp
from orbkit.display import display
from orbkit.qcinfo import QCinfo

def main_read(filename,itype='molden',all_mo=False,**kwargs):
  '''Calls the requested read function.
  
  **Parameters:**
  
    filename : str
      Specifies the filename for the input file.
    itype : str, choices={'molden', 'gamess', 'gaussian.log', 'gaussian.fchk'}
      Specifies the type of the input file.
    all_mo : bool, optional
      If True, all molecular orbitals are returned.
  
  **Returns:**
  
    qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, mo_spec, etot :
          See `Central Variables`_ for details.
  
  **Note:**
  
  
  '''
  display('Opened \n\t%s\n' % filename)
  
  # What kind of input file has to be read?
  reader = {'molden': read_molden, 
            'gamess': read_gamess, 
            'gaussian.log': read_gaussian_log, 
            'gaussian.fchk': read_gaussian_fchk,
            'wfn': read_wfn}
  display('Loading %s file...' % itype)
  
  # Return required data
  return reader[itype](filename, all_mo=all_mo,**kwargs)
  # main_read 

def read_molden(filename, all_mo=False, i_md=-1, interactive=True):
  '''Reads all information desired from a molden file.
  
  **Parameters:**
  
    filename : str
      Specifies the filename for the input file.
    all_mo : bool, optional
      If True, all molecular orbitals are returned.
    i_mo : int, default=-1
      Selects the `[Molden Format]` section of the output file.
    interactive : bool
      If True, the user is asked to select the different sets.
  
  **Returns:**
  
    qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, mo_spec, etot :
        See `Central Variables`_ for details.
  '''

  fid    = open(filename,'r')      # Open the file
  flines = fid.readlines()         # Read the WHOLE file into RAM
  fid.close()                      # Close the file
  
  # Is this really a molden file? 
  if not '[Molden Format]\n' in flines:
    display('The input file %s is no valid molden file!\nIt does' +
          ' not contain the keyword: [Molden Format]\nEXIT\n' % filename)
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
  
  count = 0
  # Go through the file line by line 
  for il in range(len(flines)):
      line = flines[il]            # The current line as string
      
      # Check the file for keywords 
      if '[Molden Format]' in line:
        count += 1
  
  if count == 0:
    display('The input file %s is no valid molden file!\nIt does' % filename +
          ' not contain the keyword: [Molden Format]\nEXIT\n' )
    raise IOError('Not a valid input file')
  else:
    if count > 1:
      display('\nContent of the molden file:')
      display('\tFound %d [Molden Format] keywords, i.e., ' % count + 
              'this file contains %d molden files.' % count)
    i_md = check_sel(count,i_md,interactive=interactive)
  
  
  # Set a counter for the AOs 
  basis_count = 0

  # Declare synonyms for molden keywords 
  synonyms = {'Sym': 'sym',
              'Ene': 'energy',
              'Occup': 'occ_num',
             }
  MO_keys = synonyms.keys()
  
  # Go through the file line by line 
  for il in range(len(flines)):
    line = flines[il]              # The current line as string
    thisline = line.split()        # The current line split into segments
    
    # Check the file for keywords 
    if '[Molden Format]' in line:
      # A new file begins 
      # Initialize the variables 
      qc = QCinfo()
      sec_flag = False             # A Flag specifying the current section 
    elif '_ENERGY=' in line:
      try:
        qc.etot = float(thisline[1])
      except IndexError:
        pass
    elif '[Atoms]' in line:
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
    elif '[GTO]' in line:
      # The section containing information about 
      # the atomic orbitals begins 
      sec_flag = 'ao_info'
      bNew = True                  # Indication for start of new AO section
    elif '[MO]' in line:
      # The section containing information about 
      # the molecular orbitals begins 
      sec_flag = 'mo_info'
      bNew = True                  # Indication for start of new MO section
    elif '[STO]' in line:
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
            if info[0] != 'Sym':
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
    mo_range   = copy.deepcopy(qc.mo_spec)
    qc.mo_spec = []
    for ii_mo in mo_range:
      if ii_mo['occ_num'] > 0.0000001:
        qc.mo_spec.append(ii_mo)
  
  # Convert geo_info and geo_spec to numpy.ndarrays
  qc.geo_info = numpy.array(qc.geo_info)
  qc.geo_spec = numpy.array(qc.geo_spec)
  
  return qc
  # read_molden 

def read_gamess(filename, all_mo=False,read_properties=False):
  '''Reads all information desired from a Gamess-US output file.
  
  **Parameters:**
  
    filename : str
      Specifies the filename for the input file.
    all_mo :   bool, optional
      If True, all molecular orbitals are returned.
  
  **Returns:**
  
    qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, mo_spec, etot :
        See `Central Variables`_ for details.
  '''
  # Initialize the variables 
  qc = QCinfo()
  sec_flag = None                 # A Flag specifying the current section
  occ = []                        # occupation number of molecular orbitals
  is_pop_ana = True               # Flag for population analysis for ground state

  # Go through the file line by line 
  with open(filename) as fileobject:
    for line in fileobject:
      thisline = line.split()      # The current line split into segments
      # Check the file for keywords 
      if ' ATOM      ATOMIC                      COORDINATES' in line:
        # The section containing information about 
        # the molecular geometry begins 
        sec_flag = 'geo_info'
        geo_skip = 1     # Number of lines to skip
        atom_count = 0  # Counter for Atoms
        if '(BOHR)' in line:
          # The length are given in Angstroem 
          # and have to be converted to Bohr radii 
          aa_to_au = 1.0
        else:
          # The length are given in Bohr radii 
          aa_to_au = 1/0.52917720859
      
      elif 'ATOMIC BASIS SET' in line:
        # The section containing information about 
        # the atomic orbitals begins 
        sec_flag = 'ao_info'
        ao_skip = 6                     # Number of lines to skip
        AO = []                         # Atomic orbitals
      elif '          EIGENVECTORS' in line:
        # The section containing information about 
        # the molecular orbitals begins 
        sec_flag = 'mo_info'
        mo_skip = 1
        init_mo = False             # Initialize new MO section
        info_key = None             # A Flag specifying the energy and symmetry section
        len_mo = 0                  # Number of MOs 
        sym={}                      # Symmetry of MOs
      elif ' NUMBER OF OCCUPIED ORBITALS (ALPHA)          =' in line:
        occ.append(int(thisline[-1]))
      elif ' NUMBER OF OCCUPIED ORBITALS (BETA )          =' in line:
        occ.append(int(thisline[-1]))
      elif '          ECP POTENTIALS' in line:
        occ = []
      elif ' NUMBER OF OCCUPIED ORBITALS (ALPHA) KEPT IS    =' in line:
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
            if 'END OF' in line and 'CALCULATION' in line:
              sec_flag = None
            else:
              if thisline == []:
                init_mo = True
                info_key = None
              elif init_mo:
                init_len = len(thisline)
                for ii in range(len(thisline)):
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
                  a= thisline[init_len-ii]
                  if a not in sym.keys(): sym[a] = 1
                  else: sym[a] += 1
                  qc.mo_spec[-ii]['sym'] = '%d.%s' % (sym[a], thisline[init_len-ii])
                info_key = 'coeffs'
              elif thisline != [] and info_key == 'coeffs':
                for ii in range(init_len,0,-1):
                  qc.mo_spec[-ii]['coeffs'].append(float(line[16:].split()[init_len-ii]))
          elif mo_skip:
            mo_skip -= 1
            
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
                         'coeffs': basis_set[qc.geo_info[kk][0]][ll]['coeffs']
                         })

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
  
  # Convert geo_info and geo_spec to numpy.ndarrays
  qc.geo_info = numpy.array(qc.geo_info)
  qc.geo_spec = numpy.array(qc.geo_spec)
  
  return qc
  # read_gamess 

def read_gaussian_fchk(filename, all_mo=False):
  '''Reads all information desired from a Gaussian FChk file. 

  **Parameters:**
  
    filename : str
      Specifies the filename for the input file.
    all_mo : bool, optional
      If True, all molecular orbitals are returned.
  
  **Returns:**
  
    qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, mo_spec, etot :
        See `Central Variables`_ for details.
  '''

  fid    = open(filename,'r')   # Open the file
  flines = fid.readlines()      # Read the WHOLE file into RAM
  fid.close()                   # Close the file
  
  sec_flag = None
  
  el_num = [0,0]
  index = 0
  at_num = 0
  ao_num = 0 
  switch = 0
  qc = QCinfo()
  
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
    #elif 'Number of basis functions' in line:
      #basis_count = int(thisline[5])
    elif 'Atomic numbers'  in line:
      sec_flag = 'geo_info'
      index = 0
      at_num = int(thisline[-1])
      count = 0
      if qc.geo_info == []:
        for ii in range(at_num):
          qc.geo_info.append(['',ii,''])
    elif 'Integer atomic weights' in line:
      sec_flag = 'geo_info'
      index = 2
      at_num = int(thisline[-1])
      count = 0
      if qc.geo_info == []:
        for ii in range(at_num):
          qc.geo_info.append(['',ii,''])
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
      index = 0
      count = len(qc.mo_spec)
      if el_num[0] == el_num[1]:
        i = el_num[0]
        occ = 2
      else:
        i = el_num[0 if 'Alpha' in line else 1]
        occ = 1
      for ii in range(mo_num):
        qc.mo_spec.append({'coeffs': numpy.zeros(basis_count),
                        'energy': 0.0,
                        'occ_num': float(occ if ii < i else 0)
                        })
        qc.mo_spec[-1]['sym'] = '%i.1' % len(qc.mo_spec)
    elif 'MO coefficients' in line:
      sec_flag = 'mo_coeffs'
      count = 0
      mo_num = int(thisline[-1])
    else:
      # Check if we are in a specific section 
      if sec_flag == 'geo_info':
        for ii in thisline:
          qc.geo_info[count][index] = ii
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
            basis_count += l_deg(lquant[ii])
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
          qc.mo_spec[index]['coeffs'][count] = float(ii)
          count += 1
          if count == basis_count:
            count = 0
            index += 1
          if index != 0 and not index % basis_count:
            sec_flag = None
  
  # Are all MOs requested for the calculation? 
  if not all_mo:
    mo_range = copy.deepcopy(qc.mo_spec)
    qc.mo_spec = []
    for ii_mo in mo_range:
      if ii_mo['occ_num'] > 0.0000001:
        qc.mo_spec.append(ii_mo)
  
  # Convert geo_info and geo_spec to numpy.ndarrays
  qc.geo_info = numpy.array(qc.geo_info)
  qc.geo_spec = numpy.array(qc.geo_spec)
  
  return qc
  # read_gaussian_fchk 

def read_gaussian_log(filename,all_mo=False,orientation='standard',
                      i_geo=-1,i_ao=-1,i_mo=-1,interactive=True):
  '''Reads all information desired from a Gaussian .log file.

  **Parameters:**
  
    filename : str
      Specifies the filename for the input file.
    all_mo :  bool, optional
      If True, all molecular orbitals are returned.
    orientation : string, choices={'input', 'standard'}, optional
      Specifies orientation of the molecule in Gaussian nomenclature. [#first]_ 
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
        See `Central Variables`_ for details.

.. [#first] Attention: The MOs in the output are only valid for the standard orientation!

  '''
  
  fid    = open(filename,'r')      # Open the file
  flines = fid.readlines()         # Read the WHOLE file into RAM
  fid.close()                      # Close the file
  
  # Search the file the specific sections
  count = {'geometry': 0, 'geometry_input': 0, 'atomic orbitals': 0, 
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
      thisline = line.split()      # The current line split into segments
      
      # Check the file for keywords 
      if ' orientation:' in line:
        if '%s orientation:' % orientation in line.lower():
          count['geometry'] += 1
        if 'input orientation:' in line.lower():
          count['geometry_input'] += 1
      elif 'Standard basis:' in line or 'General basis read from cards:' in line:
        # Check if a cartesian basis has been applied
        if '(6D, 10F)' not in line:
          raise IOError('Please apply a Cartesian Gaussian Basis Sets (6D, 10F)!')
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
    display('\tFound %d atomic orbitals section(s).' % count['atomic orbitals'])
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
  
  aa_to_au = 1/0.52917720859  # conversion factor for Angstroem to bohr radii

  # Set a counter for the AOs 
  basis_count = 0
  
  # Initialize some variables 
  sec_flag = None
  skip = 0
  c_geo = 0
  c_ao = 0
  c_mo = 0
  orb_sym = []
  qc = QCinfo()
  index = []
  
  # Go through the file line by line 
  for il in range(len(flines)):
    line = flines[il]              # The current line as string
    thisline = line.split()        # The current line split into segments
    
    # Check the file for keywords 
    if '%s orientation:' % orientation in line.lower():
      # The section containing information about 
      # the molecular geometry begins 
      if i_geo == c_geo:
        qc.geo_info = []
        qc.geo_spec = []
        sec_flag = 'geo_info'
      c_geo += 1
      skip = 4
    elif 'Standard basis:' in line:
      # Check if a cartesian basis has been applied
      if '(6D, 10F)' not in line:
        raise IOError('Please apply a Cartesian Gaussian Basis Sets (6D, 10F)!')
    elif 'AO basis set' in line:
      # The section containing information about 
      # the atomic orbitals begins
      if i_ao == c_ao:
        qc.ao_spec = []
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
        if orb_sym == []:
          if 'Alpha' in mo_type:
            add = '(a)'
          orb_sym = ['A1'+add] * basis_count
          if 'Beta' in mo_type:
            add = '(b)'
            orb_sym += ['A1'+add] * basis_count
        for i in range(len(orb_sym)):
          # for numpy version < 1.6 
          c = ((numpy.array(orb_sym[:i+1]) == orb_sym[i]) != 0).sum()
          # for numpy version >= 1.6 this could be used:
          #c = numpy.count_nonzero(numpy.array(orb_sym[:i+1]) == orb_sym[i])
          qc.mo_spec.append({'coeffs': numpy.zeros(basis_count),
                          'energy': 0.,
                          'sym': '%d.%s' % (c,orb_sym[i])})
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
      if sec_flag == 'mo_sym':
        if 'The electronic state is' in line:
          sec_flag = None
        else:
          info = line[18:].replace('(','').replace(')','').split()
          if 'Alpha' in line:
            add = '(a)'
          elif 'Beta' in line:
            add = '(b)'
          for i in info:
            orb_sym.append(i + add)   
      if sec_flag == 'mo_info':
        # Molecular orbital section 
        info = line[:21].split()
        coeffs = line[21:].split()
        if info == []:
          if bNew:
            index = [offset+i for i in range(len(coeffs))]
            bNew = False
          else:
            for i,j in enumerate(index):
              qc.mo_spec[j]['occ_num'] = int('O' in coeffs[i])
              if mo_type not in 'Alpha&Beta':
                qc.mo_spec[j]['occ_num'] *= 2
        elif 'Eigenvalues' in info:
          if mo_type == 'Natural':
            key = 'occ_num'
          else:
            key = 'energy'
          for i,j in enumerate(index):
            qc.mo_spec[j][key] = float(coeffs[i])
        else:
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
    mo_range = copy.deepcopy(qc.mo_spec)
    qc.mo_spec = []
    for ii_mo in mo_range:
      if ii_mo['occ_num'] > 0.0000001:
        qc.mo_spec.append(ii_mo)
  
  # Convert geo_info and geo_spec to numpy.ndarrays
  qc.geo_info = numpy.array(qc.geo_info)
  qc.geo_spec = numpy.array(qc.geo_spec)
  
  return qc
  # read_gaussian_log 
  

def read_wfn(filename, all_mo=False):
  '''Reads all information desired from a wfn file.
  
  **Parameters:**
  
    filename : str
      Specifies the filename for the input file.
    all_mo : bool, optional
      If True, all molecular orbitals are returned.
  
  **Returns:**
  
    qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, mo_spec, etot :
        See `Central Variables`_ for details.
  '''
  aa_to_au = 1/0.52917720859
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
  for j in exp:         
    exp_list.extend(j)

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
                'pnum': 1,
                'coeffs': [],
                'exp_list': [] 
                })
      elif 'TYPE ASSIGNMENTS' in line:
        thisline = line[18:].split()
        for i in range(len(thisline)):
          qc.ao_spec[c_type]['exp_list'] = exp_list[int(thisline[i])-1]
          c_type += 1
      elif 'EXPONENTS' in line:
        thisline = line[11:].split()
        for i in range(len(thisline)):
          qc.ao_spec[c_exp]['coeffs'] = numpy.array([[float(thisline[i]),1.0]])
          c_exp += 1
      elif 'MO' in line and 'OCC NO =' in line and 'ORB. ENERGY =' in line:
        qc.mo_spec.append({'coeffs': numpy.zeros(ao_num),
                'energy': line[25:].split()[7],
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
            qc.geo_info.append([thisline[0],thisline[1],thisline[9]])   ###FIXME Bohr or Angstroem
            qc.geo_spec.append([float(ii) for ii in thisline[4:7]])
            at_num -= 1
        elif sec_flag == 'mo_info':
          for i in thisline:
            if (c_mo) < ao_num:
              qc.mo_spec[-1]['coeffs'][c_mo] = numpy.array(float(i))
              c_mo += 1
            if (c_mo) == ao_num:
              sec_flag = None
  
  # Convert geo_info and geo_spec to numpy.ndarrays
  qc.geo_info = numpy.array(qc.geo_info)
  qc.geo_spec = numpy.array(qc.geo_spec)
  
  return qc


def mo_select(mo_spec, fid_mo_list):
  '''Selects molecular orbitals from an external file.

  **Parameters:**
   
    mo_spec :        
      See `Central Variables`_ for details.
    fid_mo_list : str, `'all_mo'`, or list
      | If fid_mo_list is a str, specifies the filename of the molecular orbitals list.
      | If fid_mo_list is 'all_mo', creates a list containing all molecular orbitals.
      | If fid_mo_list is a list, provides a list (or a list of lists) of molecular 
        orbital labels.

  **Supported Formats:**
  
    Integer List:
    
      .. literalinclude:: ../examples/MO_List_int.tab
            :language: bash
    
    List with Symmetry Labels:
    
      .. literalinclude:: ../examples/MO_List.tab
            :language: bash

  **Returns:**
  
    Dictionary with following Members:
      :mo: - List of molecular orbital labels.
      :mo_ii: - List of molecular orbital indices.
      :mo_spec: - Selected elements of mo_spec. See `Central Variables`_ for details.
      :mo_in_file: - List of molecular orbital labels within the fid_mo_list file.
      :sym_select: - If True, symmetry labels have been used. 
      
  '''
  mo_in_file = []
  all_mo = []
  selected_mo = []
  selected_mo_spec = []
  selected_mo_ii = [] 
  sym_select = False
  
  if isinstance(fid_mo_list,str) and fid_mo_list.lower() == 'all_mo':
    selected_mo = numpy.array(numpy.arange(len(mo_spec))+1, dtype=numpy.str)
    mo_in_file = [selected_mo]
    selected_mo_spec = mo_spec
  else:  
    if isinstance(fid_mo_list, list):
      if any(isinstance(i, list) for i in fid_mo_list):
        for i in fid_mo_list:
          if isinstance(i, list):
            selected_mo.extend(i)
          else:
            selected_mo.append(i)
            i = [i]
          mo_in_file.append(map(str,i))
      else:
        selected_mo.extend(fid_mo_list)
        mo_in_file.append(map(str, fid_mo_list))
    else:
      try:
        fid=open(fid_mo_list,'r')
        flines = fid.readlines()
        fid.close()
        for line in flines:
          integer = line.replace(',',' ').split()
          mo_in_file.append(integer)
          all_mo = all_mo + integer
        selected_mo=list(set(all_mo))
      except:
        raise IOError('The selected mo-list (%(m)s) is not valid!' % 
                      {'m': fid_mo_list} + '\ne.g.\n\t1\t3\n\t2\t7\t9\n')
      
    # Check if the molecular orbitals are specified by symmetry 
    # (e.g. 1.1 in MOLPRO nomenclature) or 
    # by the number in the molden file (e.g. 1)
    
    try: # Try to convert selections into integer
      for i in selected_mo: 
        int(i)
    except ValueError:
      sym_select = True
      # Add '.' after molecular orbital number if missing
      for i in range(len(selected_mo)):
        if not '.' in selected_mo[i]:
          from re import search
          a = search(r'\d+', selected_mo[i]).group()
          if a == selected_mo[i]:
            selected_mo[i] = '%s.1' % a
          else:
            selected_mo[i] = selected_mo[i].replace(a, '%s.' % a)
    
    if sym_select:
      for k in range(len(mo_spec)):
        if mo_spec[k]['sym'] in selected_mo:
          selected_mo_spec.append(mo_spec[k])
          selected_mo_ii.append(mo_spec[k]['sym'])
      selected_mo_ii = numpy.array(selected_mo_ii)
    else:
      selected_mo = map(int, selected_mo)            
      selected_mo.sort()
      selected_mo = map(str, selected_mo)
      for k in range(len(mo_spec)):
        if str(k+1) in selected_mo: selected_mo_spec.append(mo_spec[k])
  
  return {'mo': selected_mo, 'mo_ii': selected_mo_ii,
          'mo_spec': selected_mo_spec, 
          'mo_in_file': mo_in_file, 'sym_select': sym_select}