import numpy

from orbkit.display import display
from orbkit.qcinfo import CIinfo
from orbkit.read.tools import descriptor_from_file
from orbkit.units import ev_to_ha
from .tools import multiplicity

def gaussian_tddft(fname,select_state=None,threshold=0.0,**kwargs):
  '''Reads Gaussian16 TDDFT output. 
  
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
  display('\nReading data of TDDFT calculation from Gaussian...')
  # Initialize variables
  ci = []
  ci_flag = False
  prttol = False
  init_state = False
  rhfspin = 0
  spin = 'Unknown'
  nel = 0
  deex = []
  
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
    if ' SCF Done:  E(' in line:
        ci.append(CIinfo(method='tddft'))
        ci[-1].info   = []
        ci[-1].coeffs = []
        ci[-1].occ    = []
        ci[-1].occ.append([0,0])
        ci[-1].coeffs.append(1.0)
        ci[-1].info = {'state': '0',
                       'energy': float(thisline[4]),
		       'energy_nm': 0.0,
                       'fileinfo': filename,
                       'read_threshold': threshold,
                       'spin': spin,
		       'f_0i': 0.0}
    # Initialize new excited state
    elif ' Excited State' in line and 'eV' in line and 'nm' in line:
      if select_state is None or int(thisline[2].replace(':',' ')) in select_state:
        init_state = True
        tddft_skip = 1
        ci.append(CIinfo(method='tddft'))
        ci[-1].info   = []
        ci[-1].coeffs = []
        ci[-1].occ    = []
        ci[-1].info = {'state': thisline[2][:-1],
                       'energy': float(thisline[-6])*ev_to_ha + ci[0].info['energy'],
		       'energy_nm': float(thisline[-4]),
                       'fileinfo': filename,
                       'read_threshold': threshold,
                       'spin': thisline[3].split('-')[0],
		       'f_0i': float(thisline[8].replace('=',' ').split()[-1])}
        deex.append([])
    if init_state == True:
      if not tddft_skip:
        if thisline == [] or '->' not in line and '<-' not in line:
          init_state = False
        else:
          if '->' in line:
            thisline = line.replace('->','-> ').split()
            if abs(float(thisline[-1])) > threshold:
              tmp_occ = [thisline[0],thisline[2]]
              ci[-1].occ.append(tmp_occ)
              ci[-1].coeffs.append(float(thisline[-1])*numpy.sqrt(2))
          elif '<-' in line:
            deex[-1].append(float(thisline[-1])*numpy.sqrt(2))
      elif tddft_skip:
        tddft_skip -= 1

  fname.close()
  deex = numpy.array(deex)

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
