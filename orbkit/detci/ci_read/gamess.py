import numpy
from copy import copy

from orbkit.units import ev_to_ha
from orbkit.display import display
from orbkit.qcinfo import CIinfo
from orbkit.read.tools import descriptor_from_file

from .tools import multiplicity

def gamess_tddft(fname,select_state=None,threshold=0.0,**kwargs):
  '''Reads GAMESS-US TDDFT output. 
  
  **Parameters:**
  
    fname: str, file descriptor
      Specifies the filename for the input file.
      fname can also be used with a file descriptor instad of a filename.
    select_state : None or list of int, optional
      If not None, specifies the states to be read (0 corresponds to the ground 
      state), else read all electronic states.
    threshold : float, optional
      Specifies a read threshold for the CI coefficients.
  
  **Returns:**
  
    ci : list of CIinfo class instances
      See :ref:`Central Variables` for details.
  '''
  display('\nReading data of TDDFT calculation from GAMESS-US...')
  # Initialize variables
  ci = []
  ci_flag = False
  prttol = False
  init_state = False
  rhfspin = 0
  spin = 'Unknown'
  
  if isinstance(select_state,int): select_state = [select_state]

  if isinstance(fname, str):
    filename = fname
    fname = descriptor_from_file(filename, index=0, ci_descriptor=True)
  else:
    filename = fname.name

  for line in fname:
    thisline = line.split()             # The current line split into segments
    #--- Check the file for keywords ---
    # Initialize Hartree-Fock ground state
    if 'NUMBER OF ELECTRONS' in line:
      nel = int(thisline[-1])
    elif 'SPIN MULTIPLICITY' in line:
      rhfspin = int(thisline[-1])
    elif 'SINGLET EXCITATIONS' in line:
      spin = 'Singlet'
    #elif ' FINAL RHF ENERGY IS' in line and (select_state is None or 0 in select_state):
    elif ' FINAL' in line and ' ENERGY IS' in line and (select_state is None or 0 in select_state):
        ci.append(CIinfo(method='tddft'))
        ci[-1].info   = []
        ci[-1].coeffs = []
        ci[-1].occ    = []
        ci[-1].occ.append([0,0])
        ci[-1].coeffs.append(1.0)
        ci[-1].info = {'state': '0',
                       'energy': float(thisline[4]),
                       'fileinfo': filename,
                       'read_threshold': threshold,
                       'spin': spin,
                       'nel': nel}
    # Initialize new excited state
    elif 'STATE #' in line and 'ENERGY =' in line:
      if select_state is None or int(thisline[2]) in select_state:
        init_state = True
        tddft_skip = 8
        ci.append(CIinfo(method='cis'))
        ci[-1].info   = []
        ci[-1].coeffs = []
        ci[-1].occ    = []
        ci[-1].info = {'state': thisline[2],
                       'energy': float(thisline[-2])*ev_to_ha + ci[0].info['energy'],
                       'fileinfo': filename,
                       'read_threshold': threshold,
                       'spin': 'Unknown',
                       'nel': nel}
    if init_state == True and line != '\n' and 'WARNING:' not in line:
      if not tddft_skip:
        if 'NON-ABELIAN' in line or 'SUMMARY' in line or 'SYMMETRY' in line or 'STATE #' in line:
          init_state = False
        else:
          if abs(float(thisline[2])) > threshold:
            ci[-1].occ.append(thisline[:2])
            ci[-1].coeffs.append(thisline[2])
      elif tddft_skip:
        tddft_skip -= 1
          
  fname.close()
          
  #--- Calculating norm of CI states
  display('\nIn total, %d states have been read.' % len(ci)) 
  display('Norm of the states:')
  for i in range(len(ci)):
    j = numpy.array(ci[i].coeffs,dtype=float)
    norm = numpy.sum(j**2)
    ci[i].coeffs = j
    # Write Norm to log-file
    display('\tState %s:\tNorm = %0.8f (%d Coefficients)' % (ci[i].info['state'],norm, len(ci[i].coeffs)))
    # Transform to numpy arrays
    ci[i].occ = numpy.array([s for s in ci[i].occ],dtype=numpy.intc)-1
  
  return ci

def gamess_cis(filename,select_state=None,threshold=0.0,**kwargs):
  '''Reads GAMESS-US CIS output. 
  
  **Parameters:**
  
    filename : str
      Specifies the filename for the input file.
    select_state : None or list of int, optional
      If not None, specifies the states to be read (0 corresponds to the ground 
      state), else read all electronic states.
    threshold : float, optional
      Specifies a read threshold for the CI coefficients.
  
  **Returns:**
  
    ci : list of CIinfo class instances
      See :ref:`Central Variables` for details.
  '''
  display('\nReading data of CIS calculation from GAMESS-US...')
  # Initialize variables
  ci = []
  ci_flag = False
  prttol = False
  init_state = False
  rhfspin = 0
  min_c = -1
  
  if isinstance(select_state,int): select_state = [select_state]
  with open(filename) as fileobject:
    for line in fileobject:
      thisline = line.split()             # The current line split into segments
      #--- Check the file for keywords ---
      # Initialize Hartree-Fock ground state
      if 'NUMBER OF ELECTRONS' in line and "=" in line:
        nel = int(thisline[-1])
      elif 'SPIN MULTIPLICITY' in line:
        rhfspin = int(thisline[-1])
      elif ' FINAL RHF ENERGY IS' in line and (select_state is None or 0 in select_state):
          ci.append(CIinfo(method='cis'))
          ci[-1].info   = []
          ci[-1].coeffs = []
          ci[-1].occ    = []
          ci[-1].occ.append([0,0])
          ci[-1].coeffs.append(1.0)
          ci[-1].info = {'state': '0',
                         'energy': float(thisline[4]),
                         'fileinfo': filename,
                         'read_threshold': threshold,
                         'spin': multiplicity()[rhfspin],
                         'nel': nel}
      # Printing parameter
      elif ' PRINTING CIS COEFFICIENTS LARGER THAN' in line:
        min_c = float(thisline[-1])
      # Initialize new excited state
      elif ' EXCITED STATE ' in line and 'ENERGY=' and 'SPACE SYM' in line:
        if select_state is None or int(thisline[2]) in select_state:
          init_state = True
          cis_skip = 6
          ci.append(CIinfo(method='cis'))
          ci[-1].info   = []
          ci[-1].coeffs = []
          ci[-1].occ    = []
          ci[-1].info = {'state': thisline[2],
                         'energy': float(thisline[4]),
                         'fileinfo': filename,
                         'read_threshold': threshold,
                         'spin': multiplicity()[int(2*float(thisline[7])+1)],
                         'nel': nel}
      if init_state == True:
        if not cis_skip:
          if '----------------------------------------------' in line:
            init_state = False
          else:
            if abs(float(thisline[2])) > threshold:
              ci[-1].occ.append(thisline[:2])
              ci[-1].coeffs.append(thisline[2])
        elif cis_skip:
          cis_skip -= 1
          
  #--- Calculating norm of CI states
  display('\nIn total, %d states have been read.' % len(ci)) 
  display('Norm of the states:')
  for i in range(len(ci)):
    j = numpy.array(ci[i].coeffs,dtype=float)
    norm = numpy.sum(j**2)
    ci[i].coeffs = j
    # Write Norm to log-file
    display('\tState %s:\tNorm = %0.8f (%d Coefficients)' % (ci[i].info['state'],norm, len(ci[i].coeffs)))
    # Transform to numpy arrays
    ci[i].occ = numpy.array([s for s in ci[i].occ],dtype=numpy.intc)-1
  display('')
  if min_c > threshold:
    display('\nInfo:'+
       '\n\tSmallest coefficient (|c|=%f) is larger than the read threshold (%f).' 
       %(min_c,threshold) + 
       '\n\tUse `PRTTOL=0.0` in the `$CIS` input card to print ' +
       'all CI coefficients.\n')
  
  return ci
