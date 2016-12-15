from ..display import display
import numpy
from copy import copy
from ..qcinfo import QCinfo,CIinfo

def main_ci_read(qc,filename,itype='psi4_detci',threshold=0.0,
                 select=None,nforbs=0,bortho=False,
                 **kwargs):
  '''Reads determinant CI calculation. 
  
  Supported input files (``itype``):
  
  itype           QC-Program  Level of Theory
  ==============  ==========  =====================
  'psi4_detci'    PSI4        All CI calculations
  'gamess_cis'    GAMESS-US   CIS
  'tmol_tddft'    TURBOMOLE   TD-DFT
  'molpro_mcscf'  MOLPRO      MCSCF
  ==============  ==========  =====================
  
  **Parameters:**
  
    qc : class QCinfo
      See :ref:`Central Variables` for details.
    filename : str
      Specifies the filename for the input file.
    itype : str, choices={'psi4_detci', 'gamess_cis', 'tmol_tddft', 'molpro_mcscf'}
      Specifies the type of the input file.
    threshold : float, optional
      Specifies a read threshold for the CI coefficients.
    select : None or (list of) int (counting from zero), optional
      if None, reads all states, else ...
      | TURBOMOLE & GAMESS-US: Specifies the STATE to be read (0 -> ground state)
      | PSI4 & MOLPRO: Specifies the CALCULATION to be read
    nforbs : int, optional, TURBOMOLE-specific
      Specifies the number of frozen orbitals.
    bortho : bool, optional, TURBOMOLE-specific
      If True, an orthonormalization of  the TD-DFT coefficients is perfomed
      with Gram-Schmidt.
  
  **Returns:**
  
    qc : class QCinfo
      New instance of QCinfo. See :ref:`Central Variables` for details.
      qc.mo_spec is possibly reordered.
    ci : list of CIinfo class instances 
      See :ref:`Central Variables` for details.
  '''
  display('Opened \n\t%s\n' % filename)  
  assert isinstance(qc,QCinfo), '`qc` has to be an instance of the QCinfo class'
  
  # What kind of input file has to be read?
  reader = {'psi4_detci': psi4_detci, 
            'gamess_cis': gamess_cis, 
            'tmol_tddft': tmol_tddft, 
            'molpro_mcscf': molpro_mcscf}
  
  display('Loading %s file...' % itype)
  if itype not in reader.keys():
    display('Available reader (`itype`) are:\n  ' + ', '.join(reader.keys()))
    raise NotImplementedError("itype='%s' not implemented!"%itype)
  
  kwargs['nmoocc'] = qc.get_nmoocc()
  
  ci = reader[itype](filename,select_state=select,threshold=threshold, 
                     select_run=select,                 # PSI4/MOLPRO specific
                     nforbs=nforbs,bortho=bortho,       # TURBOMOLE specific
                     **kwargs)
  
  # Get a copy of qc
  qc = qc.copy()
  if itype in ['tmol_tddft','gamess_cis']: # CIS-like
    moocc = qc.get_mo_occ()
  elif itype in ['psi4_detci','molpro_mcscf']: # detCI-like
    # Reorder qc.mo_spec
    irreps=None if 'irreps' not in ci[0].info.keys() else ci[0].info['irreps']                                                
    closed,active,external = molpro_mo_order_ci(ci[0].info['occ_info'],
                                                qc.mo_spec,
                                                irreps=irreps
                                                )
    qc.mo_spec = closed+active+external
    moocc = numpy.zeros(len(closed),dtype=numpy.intc) + 2
  
  # Add moocc to CI class instances
  for i in ci:
    i.set_moocc(moocc)
  
  return ci

def psi4_detci(filename,select_run=None,threshold=0.0,**kwargs):
  '''Reads PSI4 DETCI output. 
  
  **Parameters:**
  
    filename : str
      Specifies the filename for the input file.
    select_run : (list of) int
      Specifies the DETCI calculation (D E T C I) to be read. 
      For the selected DETCI calculation, all electronic states will be read.
    threshold : float, optional
      Specifies a read threshold for the CI coefficients.
  
  **Returns:**
  
    ci : list of CIinfo class instances 
      See :ref:`Central Variables` for details.
  
  ..#ATTENTION: Changed return value to list of CI classes
  '''
  display('\nReading data of DETCI calculation from PSI4...')
  count = 0
  numci = []
  # Go through the file line by line 
  with open(filename) as fileobject:
    for line in fileobject:
      if 'D E T C I' in line:
        count += 1
        numci.append(0)
      if '* ROOT' in line:
        numci[-1] += 1
  
  if count == 0:
    display('The input file %s does not contain any DETCI calculations!\n' % (filename) + 
            'It does not contain the keyword:\n\tD E T C I')
    raise IOError('Not a valid input file')
  else:
    display('The input file %s contains %d DETCI calculation(s)' % (filename,count))
    display('with %s root(s)%s.'%(', '.join(map(str,numci)),
                               ', respectively' if len(numci)>1 else ''))
    
    if select_run is None:
      select_run = numpy.arange(count)   
    if isinstance(select_run,int) and 0 <= select_run < count:
      select_run = [select_run]
      ci = [[]]
    elif isinstance(select_run,(list,numpy.ndarray)):
      ci = []
      for i in select_run:
        ci.append([])
        if not isinstance(i,int):
          raise IOError(str(i) + ' is not a valid selection for select_run')
    else: 
      raise IOError(str(select_run) + ' is a not valid selection for select_run')    
  select_run = numpy.array(select_run)
  
  display('\n\tYour selection (counting from zero): %s' % 
          ', '.join(map(str,select_run)))
  
  general_information = {'fileinfo': filename,
                         'read_threshold': threshold}
  
  ci_skip = 0
  count = 0
  count_runs = 0
  min_c = 1
  orbs_in_ci = 0
  start_reading = False
  
  occ_types = ['core','closed','active','external']
  synonyms = {'frozen docc': 'core',
              'restricted docc': 'closed',
              'ras ': 'active',
              'active': 'active',
              'restricted uocc': 'external',
              'frozen uocc': 'external'}
  irreps = []
  fist_active_mo = []
  index_active_mo = []
  occ_symbols = {'A': numpy.array([1,0]),
                 'B': numpy.array([0,1]),
                 'X': numpy.array([1,1])}
  with open(filename) as fileobject:
    for line in fileobject:
      thisline = line.split()             # The current line split into segments
      
      if 'Running in ' in line and 'symmetry.' in line: 
        nIRREP = point_groups[thisline[2].lower()]
        rhf_occ = numpy.zeros(nIRREP,dtype=numpy.intc)
      #--- RHF occupation
      elif 'Final Occupation by Irrep:' in line:
        irreps = fileobject.next().split()
        line = fileobject.next()
        c_occ = line.replace(',','').split()[2:-1]
        for ii in range(len(c_occ)):
          rhf_occ[ii] = int(c_occ[ii])
      #--- A DETCI Calculation starts ---
      elif 'D E T C I' in line:
        occ_info = {}
        for i in occ_types:
          occ_info[i] = numpy.zeros(nIRREP,dtype=numpy.intc)
        info = {}
        num_roots = 0
        method = 'detci'
        i = numpy.argwhere(select_run == count_runs)
        if len(i) > 0:
          index_run = int(i)
          start_reading = True
        count_runs += 1
      elif start_reading:
        if 'NUM ROOTS     =' in line:
          num_roots = int(thisline[3])
        elif 'FCI          =' in line and line.endswith('yes'):
          method = 'fci'
        elif 'REF SYM       =' in line:
          info['sym'] = thisline[3] 
          if info['sym'] != 'auto':
            info['sym'] = irreps[int(info['sym'])]
        elif 'S             =' in  line:
          info['spin'] = multiplicity[int(2*float(thisline[2])+1)]
        elif 'NUM ALP      =' in line:
          info['nel'] = int(thisline[3]) + int(thisline[-1])
        elif 'ORBS IN CI   =' in line: 
          orbs_in_ci = int(thisline[-1])
        elif '=' in line and any([i in line.lower() for i in synonyms.keys()]):
          for i in synonyms.keys():
            if i in line.lower():
              occ_info[synonyms[i]] += numpy.array(thisline[-nIRREP:],
                                                   dtype=numpy.intc)
        elif '* ROOT' in line and 'CI total energy' in line:
          ci[index_run].append(CIinfo(method=method))
          ci[index_run][-1].info = copy(general_information) 
          ci[index_run][-1].info['fileinfo'] += '@%d' % index_run
          ci[index_run][-1].info['irreps'] = irreps
          ci[index_run][-1].info['state'] = '%s.%s'%(thisline[2],info['sym'])
          ci[index_run][-1].info['energy'] = float(thisline[7])
          ci[index_run][-1].info['spin'] = info['spin']
          ci[index_run][-1].info['nel'] = info['nel']
          ci[index_run][-1].info['occ_info'] = occ_info
          closed = []
          active = {}
          c = 0
          for i in range(nIRREP):
            a = occ_info['core'][i] + occ_info['closed'][i]
            for b in range(occ_info['active'][i]):
              active['%d%s' % (a+b+1,irreps[i])] = c
              c += 1
        elif 'most important determinants' in line:
          num_det = int(thisline[1])
          ci[index_run][-1].coeffs = [] # numpy.zeros(num_det)
          ci[index_run][-1].occ = [] # numpy.zeros((numpy.sum(occ_info['active']),2),
                                          #dtype=numpy.intc)
          fileobject.next()
          def rm(s,w='*(,)'):
            for i in w:
              s = s.replace(i,' ')
            return s
          for i in range(num_det):            
            thisline = rm(fileobject.next()).split()
            if thisline == []:
              break
            min_c = min(min_c,abs(float(thisline[1])))
            if abs(float(thisline[1])) > threshold:
              #INFO: beta alpha has to have the opposite sign of alpha beta
              f = -1. if int(thisline[2]) > int(thisline[3]) else 1.
              ci[index_run][-1].coeffs.append(f * float(thisline[1]))
              ci[index_run][-1].occ.append(numpy.zeros((numpy.sum(occ_info['active']),2),
                                            dtype=numpy.intc))
              for j in thisline[4:]:
                ci[index_run][-1].occ[-1][active[j[:-1]]] = occ_symbols[j[-1]]
          ci[index_run][-1].coeffs = numpy.array(ci[index_run][-1].coeffs)
          ci[index_run][-1].occ = numpy.array(ci[index_run][-1].occ,dtype=numpy.intc)
        elif 'A good bug is a dead bug' in line:
          start_reading = False
  #--- Calculating norm of CI states
  display('\nIn total, %d states have been read.' % sum([len(i) for i in ci])) 
  display('Norm of the states:')
  for i in range(len(ci)):
    for j in range(len(ci[i])):
      norm = sum(ci[i][j].coeffs**2)
      display('\tState %s (%s):\tNorm = %0.8f (%d Coefficients)' % 
              (ci[i][j].info['state'],ci[i][j].info['spin'],
               norm,len(ci[i][j].coeffs)))
  display('')
  if threshold and min_c > threshold:
    display('\nInfo:'+
       '\n\tSmallest coefficient (|c|=%f) larger than the read threshold (%f).' 
       %(min_c,threshold) + 
       '\n\tUse `set num_dets_print -1` in the PSI4 input file to print ' +
       'all CI coefficients.\n')
  
  ci_new = []
  for i in ci:
    ci_new.extend(i)
  
  return ci_new

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
      if 'NUMBER OF ELECTRONS' in line:
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
                         'spin': multiplicity[rhfspin],
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
                         'spin': multiplicity[int(2*float(thisline[7])+1)],
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

def tmol_tddft(filename,nmoocc=None,nforbs=0,select_state=None,threshold=0.0,
               bortho=False,**kwargs):
  '''Reads Turbomole TD-DFT output (``sing_a`` file). 
  
  Hint: Only implemented for singlet TD-DFT computations in C_1 symmetry, the
  output data file of which is usually called ``sing_a``.
  
  **Parameters:**
  
    filename : str
      Specifies the filename for the input file.
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
  with open(filename,'r') as f:
    start = False
    for i,line in enumerate(f):
      if 'eigenvalue' in line:
        start = True
        l = line.split()
        states.append({'eigenvalue': float(l[-1].replace('D','E')),
                       'coeffs': []})
      elif start:
        for i in range((len(line)-1)/20):
          states[-1]['coeffs'].append(float(line[i*20:(i+1)*20].replace('D','E')))
  
  for i in range(len(states)):
    states[i]['coeffs'] = numpy.array(states[i]['coeffs']).reshape((2,nmoocc,-1))
    states[i]['xia'] = 0.5 * (
                        (states[i]['coeffs'][0] + states[i]['coeffs'][1])
                        )
    states[i]['yia'] = 0.5 * (
                        (states[i]['coeffs'][0] - states[i]['coeffs'][1])
                        )
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
                     'energy': 0.0 if i==0 else numpy.abs(states[i+1]['eigenvalue']),
                     'fileinfo': filename,
                     'read_threshold': threshold,
                     'spin': 'Unknown',
                     'nel': nel}
      if not i:
        ci[-1].occ.append([-1,-1])
        for jj in range(states[i]['xia'].shape[0]):
          for kk in range(states[i]['xia'].shape[1]):
            ci[-1].occ.append([0,0])
      else:
        ci[-1].occ.append([0,0])
        for jj in range(states[i]['xia'].shape[0]):
          for kk in range(states[i]['xia'].shape[1]):
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
turbomole_tddft = tmol_tddft

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
    display('The input file %s is no valid escf file!\n\nIt does'  % filename+
          ' not contain the keyword: e s c f\n')
    raise IOError('Not a valid input file')
  
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

def molpro_mcscf(filename,select_run=0,threshold=0.0,**kwargs):
  '''Reads MOLPRO MCSCF output. 
  
  **Parameters:**
  
    filename : str
      Specifies the filename for the input file.
    select_run : (list of) int
      Specifies the MCSCF calculation (1PROGRAM * MULTI) to be read. 
      For the selected MCSCF calculation, all electronic states will be read.
    threshold : float, optional
      Specifies a read threshold for the CI coefficients.
  
  **Returns:**
  
    ci : list of CIinfo class instances 
      See :ref:`Central Variables` for details.
  
  ..#ATTENTION: Changed return value to list of CI classes
  '''
  display('\nReading data of MCSCF calculation from MOLPRO...')
  method = 'mcscf'
  available = {'mcscf': 'MULTI'}
  single_run_selected = isinstance(select_run,int)
  count = 0
  # Go through the file line by line 
  with open(filename) as fileobject:
    for line in fileobject:
      if '1PROGRAM * %s' % available[method] in line:
        count += 1
  
  if count == 0:
    display('The input file %s does not contain any %s calculations!\n' % (filename,method) + 
            'It does not contain the keyword:\n\t1PROGRAM * %s' % available[method])
    raise IOError('Not a valid input file')
  else:
    display('The input file %s contains %d %s calculation(s).' % (filename,count,method) )
    if select_run is None:
      select_run = numpy.arange(count)   
    elif isinstance(select_run,int) and 0 <= select_run < count:
      display('\tYour selection (counting from zero): %d' % select_run )
      select_run = [select_run]
      ci = [[]]
    elif isinstance(select_run,list):
      ci = []
      for i in select_run:
        ci.append([])
        if not isinstance(i,int):
          raise IOError(str(i) + ' is not a valid selection for select_run')
    else: 
      raise IOError(str(select_run) + ' is a not valid selection for select_run')    
  select_run = numpy.array(select_run)
  
  display('\tYour selection (counting from zero): %s' % 
          ', '.join(map(str,select_run)))
  
  general_information = {'fileinfo': filename,
                         'read_threshold': threshold}
  
  ci_skip = 0
  count = 0
  count_runs = 0
  min_c = 1
  start_reading = False
  sec_flag = False

  info_flag = False
  info_sep = False
  info_split = False
    
  occ_types = ['core','closed','active','external']
  
  with open(filename) as fileobject:
    for line in fileobject:
      thisline = line.split()             # The current line split into segments
      
      if '1PROGRAM *' in line:
        start_reading = False
      #--- Number of IRREPs
      if 'Point group' in line or '_PGROUP' in line: 
        nIRREP = point_groups[thisline[-1].lower()]
        rhf_occ = numpy.zeros(nIRREP,dtype=numpy.intc)
      #--- RHF occupation
      elif 'Final occupancy:' in line:
        c_occ = line.split()[2:]
        for ii in range(len(c_occ)):
          rhf_occ[ii] = int(c_occ[ii])
      #--- A MCSCF Calculation starts ---
      elif '1PROGRAM * %s' % available[method] in line:
        occ_info = {}
        for i in occ_types:
          occ_info[i] = numpy.zeros(nIRREP,dtype=numpy.intc)
        state_info = []
        i = numpy.argwhere(select_run == count_runs)
        if len(i) > 0:
          index_run = int(i)
          start_reading = True
          count = 0
          old = 0
        count_runs += 1
      elif start_reading:
        #--- Active space ---
        if 'Number of ' in line and 'orbitals:' in line:          
          line = line.replace('(','').replace(')','').replace('-shell','')   
          c_occ = numpy.array(line.split()[-nIRREP:],dtype=numpy.intc)
          occ_info[line.split()[2]] += c_occ
        elif 'State symmetry' in line:
          fileobject.next()
          thisline = fileobject.next()
          if 'State symmetry' in thisline:
            fileobject.next()
            thisline = fileobject.next()
          thisline = thisline.replace('=',' ').split()
          data = {'nel': thisline[3],
                  'spin': thisline[6],
                  'sym': thisline[9]}
          thisline = fileobject.next().split()
          state_info.extend([data for i in range(int(thisline[-1]))])
        elif '!%s state' % method in line.lower() and 'energy' in line.lower():
          ci[index_run].append(CIinfo(method=method))
          info = state_info[count]
          thisline = line.lower().replace('state','state ').split()
          ci[index_run][count].info = copy(general_information) 
          ci[index_run][count].info['fileinfo'] += '@%d' % index_run
          ci[index_run][count].info['state'] = thisline[2]
          ci[index_run][count].info['energy'] = float(thisline[4])
          ci[index_run][count].info['spin'] = info['spin']
          ci[index_run][count].info['nel'] = info['nel']
          ci[index_run][count].info['occ_info'] = occ_info
          count += 1
        elif 'CI vector' in line:
          sec_flag = 'mcscf'
          info_split = '     '
          ci_skip = 3
          info = thisline[-1]
          count = old
          first = True
        if not ci_skip:
          if line == '\n' or '/EOF' in line:
            sec_flag = False
          elif sec_flag != False:
            split_line = filter(None, line.split(info_split))
            if len(split_line) > 1:
              occupation = numpy.zeros((numpy.sum(occ_info['active']),2),dtype=numpy.intc)
              for i,j in enumerate(split_line[0].replace(' ','')):
                if j == '2': 
                  occupation[i,:] = 1
                elif j == 'a':
                  occupation[i,0] = 1
                elif j == 'b':
                  occupation[i,1] = 1
              c0 = 0
            coeffs = split_line[-1].split()
            for i,j in enumerate(coeffs):
              if first:
                old += 1
              min_c = min(min_c,abs(float(j)))
              if abs(float(j)) > threshold:
                ci[index_run][count + c0 + i].coeffs.append(float(j))
                ci[index_run][count + c0 + i].occ.append(occupation)
            c0 += len(coeffs)
            first = False
        else:
          ci_skip -= 1 
  
  #--- Calculating norm of CI states
  display('\nIn total, %d states have been read.' % sum([len(i) for i in ci])) 
  display('Norm of the states:')
  for i in range(len(ci)):
    for j in range(len(ci[i])):
      ci[i][j].coeffs = numpy.array(ci[i][j].coeffs)
      ci[i][j].occ = numpy.array(ci[i][j].occ,dtype=numpy.intc)
      norm = sum(ci[i][j].coeffs**2)
        # Write Norm to log-file
      display('\tState %s (%s):\tNorm = %0.8f (%d Coefficients)' % 
              (ci[i][j].info['state'],ci[i][j].info['spin'],
               norm,len(ci[i][j].coeffs)))
  display('')
  if min_c > threshold:
    display('\nInfo:'+
       '\n\tSmallest coefficient (|c|=%f) larger than the read threshold (%f).' 
       %(min_c,threshold) + 
       '\n\tUse `gthresh, printci=0.0` in the MOLPRO input file to print ' +
       'all CI coefficients.\n')
  
  #if single_run_selected and len(ci) == 1:
    #ci = ci[0]
  ci_new = []
  for i in ci:
    ci_new.extend(i)
  
  return ci_new

def molpro_mo_order_ci(occ_info,mo_spec,irreps=None,nIRREP=None,order_sym=False):
  '''Orders the molecular orbitals according to the occupation information 
  of a CI class. (For processing MCSCF MOLPRO output)  
  
  **Parameters:**
  
    occ_info : ci.info["occ_info"]
      Contains a dict with the number of orbitals in the different 
      IRREPS for "core", "closed", "active", and "external" orbitals.
    mo_spec : list of dictionaries 
      See :ref:`Central Variables` for details.
  
  **Returns:**
  
    closed,active,external : list of mo_spec
      Contains mo_spec assigned to the "closed", "active", and "external" 
      orbitals.
  '''
  if nIRREP is None: nIRREP = len(occ_info['closed']) 
  if irreps is None: 
    irreps = dict(zip(map(str,range(1,nIRREP+1)),range(nIRREP)))
  else:
    nIRREP = len(irreps)
    irreps = dict(zip(irreps,range(nIRREP)))
  
  mo_index = [[] for i in range(nIRREP)]
  for i,i_mo in enumerate(mo_spec):
    info = i_mo['sym'].split('.')
    mo_index[irreps[info[1]]].append((i,int(info[0])))
  
  mo_sym = [[] for i in range(nIRREP)]
  for i in range(nIRREP):
    if mo_index[i] == []:
      continue
    mo_index[i] = numpy.array(mo_index[i])
    j = numpy.argsort(mo_index[i],axis=0)[:,1]
    mo_index[i] = mo_index[i][j]
    mo_sym[i] = [mo_spec[j] for j,k in mo_index[i]]
  if order_sym:
    return mo_sym
  
  closed = []
  active = []
  external = []
  for i in range(nIRREP):
    a = occ_info['core'][i] + occ_info['closed'][i]
    closed.extend(mo_sym[i][:a])
    b = a+occ_info['active'][i]
    active.extend(mo_sym[i][a:b])
    external.extend(mo_sym[i][b:])
  return closed,active,external

def orthonorm(ci):
  '''Orthonomalize CI coefficients after Grahm-Schmidt 
  
  **Parameters:**
  
    ci: 
  
  **Returns:**
  
    ci: 
      Contains the orthonormalized CI coefficients for each Slater-determinant
  '''
  # Grahm-Schmidt orthonormalization
  display('Orthonormalization of CI-coefficients with Gram-Schmidt...\n')
  for ii in range(len(ci)):
    scal = numpy.dot(ci[ii].coeffs,ci[ii].coeffs)
    ci[ii].coeffs /= numpy.sqrt(scal)
    for jj in range(ii+1,len(ci)):
      scal = numpy.dot(ci[ii].coeffs,ci[jj].coeffs)
      ci[jj].coeffs -= ci[ii].coeffs*scal
      
  return ci

point_groups = {'c1': 1, 'cs': 2, 'c2': 2, 'ci': 2, 
                'c2v': 4, 'c2h': 4, 'd2': 4, 'd2h': 8}

multiplicity = ['Unknown','Singlet','Doublet','Triplet','Quartet','Quintet']