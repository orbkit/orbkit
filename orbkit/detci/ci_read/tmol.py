import numpy
from copy import copy

from orbkit.display import display
from orbkit.qcinfo import CIinfo
from orbkit.read.tools import descriptor_from_file

from .tools import orthonorm

def tmol_tddft(fname,nmoocc=None,nforbs=0,select_state=None,threshold=0.0,
               bortho=False,**kwargs):
  '''Reads Turbomole TD-DFT output (``sing_a`` file). 
  
  Hint: Only implemented for singlet TD-DFT computations in C_1 symmetry, the
  output data file of which is usually called ``sing_a``.
  
  **Parameters:**
  
    fname: str, file descriptor
      Specifies the filename for the input file.
      fname can also be used with a file descriptor instad of a filename.
    nmoocc : int
      Specifies the number of non-frozen occupied orbitals.
    nforbs : int,optional
      Specifies the number of frozen orbitals.
    select_state : None or list of int, optional
      If not None, specifies the states to be read (0 corresponds to the ground 
      state), else read all electronic states. 
    threshold : float, optional
      Specifies a read threshold for the TD-DFT coefficients.
    bortho : bool
      If True, an orthonormalization of  the TD-DFT coefficients is perfomed
      with Gram-Schmidt.
  
  **Returns:**
  
    ci : list of CIinfo class instances
      See :ref:`Central Variables` for details.
  '''
  display('\nReading data of TD-DFT calculation from Turbomole...')
  
  nel = (nmoocc + nforbs) * 2
  lumo = (nmoocc + nforbs)
  states = []
  if isinstance(select_state,int): select_state = [select_state]

  if isinstance(fname, str):
    filename = fname
    fname = descriptor_from_file(filename, index=0, ci_descriptor=True)
  else:
    filename = fname.name

  start = False
  for i,line in enumerate(fname):
    if '$tensor space dimension' in line:
      tspace = int(line.split()[-1])
    if 'eigenvalue' in line:
      start = True
      l = line.split()
      states.append({'eigenvalue': float(l[-1].replace('D','E')),
                     'coeffs': []})
    elif start:
      for i in range((len(line)-1)/20):
        states[-1]['coeffs'].append(float(line[i*20:(i+1)*20].replace('D','E')))
  
  for i in range(len(states)):
    if len(states[i]['coeffs']) == tspace*2:
      states[i]['coeffs'] = numpy.array(states[i]['coeffs']).reshape((2,nmoocc,-1))
      states[i]['xia'] = 0.5 * (
                          (states[i]['coeffs'][0] + states[i]['coeffs'][1])
                          )
      #states[i]['yia'] = 0.5 * (
                          #(states[i]['coeffs'][0] - states[i]['coeffs'][1])
                          #)
    elif len(states[i]['coeffs']) == tspace:
      states[i]['xia'] = numpy.array(states[i]['coeffs']).reshape((nmoocc,-1))
    else:
      raise IOError('Shape of coefficient matrix does not match with tensor space.')
    del states[i]['coeffs']

  coeffs = numpy.zeros((len(states)+1,numpy.prod(states[1]['xia'].shape)+1))
  for i in range(len(states)):
    coeffs[i+1][1:] = states[i]['xia'].reshape((-1,))
  coeffs[0,0] = 1.0
  
  # Set up configurations
  ci = []
  for i in range(len(states)+1):
    if select_state is None or i in select_state:
      ci.append(CIinfo(method='tddft'))
      ci[-1].coeffs = coeffs[i]
      ci[-1].occ    = []
      ci[-1].info = {'state': str(i),
                     'energy': 0.0 if i==0 else numpy.abs(states[i-1]['eigenvalue']),
                     'fileinfo': filename,
                     'read_threshold': threshold,
                     'spin': 'Unknown',
                     'nel': nel}
      if not i:
        ci[-1].occ.append([-1,-1])
      else:
        ci[-1].occ.append([0,0])
        for jj in range(states[i-1]['xia'].shape[0]):
          for kk in range(states[i-1]['xia'].shape[1]):
            ci[-1].occ.append([jj+nforbs,kk+lumo])
      ci[-1].occ = numpy.array(ci[-1].occ,dtype=numpy.intc)

  # Gram-Schmidt
  display('Orthonormalizing the TD-DFT coefficients with Gram-Schmidt...\n')
  ci = orthonorm(ci)

  if bortho:
    for c in ci:
      c.apply_threshold(threshold,keep_length=True)
    ci = orthonorm(ci)
  else:
    for c in ci:
      c.apply_threshold(threshold,keep_length=False)
  
  for c in ci:
      c.apply_threshold(0.0,keep_length=False)
    
  #--- Calculating norm of CI states
  display('\nIn total, %d states have been read.' % len(ci)) 
  display('Norm of the states:')
  for i in ci:
    display(str(i))
  
  return ci

def tmol_escf(filename,ene_unit='nm'):
  '''Reads Turbomole TD-DFT output (``escf`` file). 
  
  **Parameters:**
  
    filename : str
      Specifies the filename for the input file.
    ene_unit : str
      Specifies in which units the excitation energies are read out.
      [nm,E_h,eV,cm-1]
  
  **Returns:**
  
    ex_ene,osc_str : Excitation energies and oscillator strength for excited states
    from Turbomole TD-DFT output. 
  '''

  fid    = open(filename,'r')      # Open the file
  flines = fid.readlines()         # Read the WHOLE file into RAM
  fid.close()                      # Close the file
  
  # Is this really a escf file? 
  if not '                                e s c f\n' in flines:
    raise IOError('The input file %s is no valid escf file!\n\nIt does'  % filename+
          ' not contain the keyword: e s c f\n')
  
  # Initialize variables
  osc_str = []
  ex_ene = []
  ncount = 0
  sec_flag = None
  
  # Select the excitation energy
  if ene_unit == 'nm':
    ene_str = 'Excitation energy / nm:'
  elif ene_unit == 'E_h':
    ene_str = 'Excitation energy:'
  elif ene_unit == 'eV':
    ene_str = 'Excitation energy / eV:'
  elif ene_unit == 'cm-1':
    ene_str = 'Excitation energy / cm^(-1):'
  else: 
    display('The units selection for the excitation energies is not valid.' +
            'The excitation energies are read out in units of Hartree.\n')
    ene_str = 'Excitation energy:'

  for il in range(len(flines)):
    line = flines[il]            # The current line as string
    thisline = line.split()      # The current line split into segments
    # Check the file for keywords 
    if ene_str in line:
      ex_ene.append(float(thisline[-1]))
    elif 'Oscillator strength:' in line:
      sec_flag = 'osc'
      count = 1
    else:
      if sec_flag == 'osc':
        if count == 0:
          osc_str.append(float(thisline[-1]))
          ncount += 1
          sec_flag = None
        else:
          count -= 1

  return numpy.array(ex_ene),numpy.array(osc_str)
