import numpy
from copy import copy

from orbkit.display import display
from orbkit.qcinfo import CIinfo
from orbkit.read.tools import descriptor_from_file

from .tools import point_groups

def molpro_mcscf(fname,select_run=0,threshold=0.0,**kwargs):
  '''Reads MOLPRO MCSCF output. 
  
  **Parameters:**
  
    fname: str, file descriptor
      Specifies the filename for the input file.
      fname can also be used with a file descriptor instad of a filename.
    select_run : (list of) int
      Specifies the MCSCF calculation (1PROGRAM * MULTI) to be read. 
      For the selected MCSCF calculation, all electronic states will be read.
    threshold : float, optional
      Specifies a read threshold for the CI coefficients.
  
  **Returns:**
  
    ci : list of CIinfo class instances 
      See :ref:`Central Variables` for details.
  

    ..ATTENTION: Changed return value to list of CI classes
  '''
  display('\nReading data of MCSCF calculation from MOLPRO...')
  method = 'mcscf'
  available = {'mcscf': 'MULTI'}
  single_run_selected = isinstance(select_run,int)
  count = 0
  numci = []

  if isinstance(fname, str):
    filename = fname
    fname = descriptor_from_file(filename, index=0, ci_descriptor=True)
  else:
    filename = fname.name

  from io import TextIOWrapper
  if isinstance(fname, TextIOWrapper):
    flines = fname.readlines()      # Read the WHOLE file into RAM
  else:
    magic = 'This is an Orbkit magic string'
    text = fname.read().decode("iso-8859-1").replace('\n','\n{}'.format(magic))
    flines = text.split(magic)
    flines.pop()

  # Go through the file line by line 
  for line in flines:
    if '1PROGRAM * %s' % available[method] in line:
      count += 1
      numci.append(0)
    if '!%s state' % method in line.lower() and 'energy' in line.lower():
      numci[-1] += 1
          
  if count == 0:
    display('The input file %s does not contain any DETCI calculations!\n' % (filename) + 
            'It does not contain the keyword:\n\t1PROGRAM * MULTI')
    raise IOError('Not a valid input file')
  else:
    string = ', '.join(map(str,numci)).rsplit(', ',1)
    string = ' and '.join(string) if len(string[0].split(',')) < 2 else ', and '.join(string)
    display('The input file %s contains' % (filename))
    display('%d MCSCF calculation(s) with %s root(s)%s.'%(count,string,
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
  
  flines = iter(flines)
  for line in flines:
    thisline = line.split()             # The current line split into segments
    
    if '1PROGRAM *' in line:
      start_reading = False
    #--- Number of IRREPs
    if 'Point group' in line or '_PGROUP' in line: 
      nIRREP = point_groups()[thisline[-1].lower()]
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
        next(flines)
        thisline = next(flines)
        if 'State symmetry' in thisline:
          next(flines)
          thisline = next(flines)
        thisline = thisline.replace('=',' ').split()
        data = {'nel': thisline[3],
                'spin': thisline[6],
                'sym': thisline[9]}
        thisline = next(flines).split()
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
          split_line = list(filter(None, line.split(info_split)))
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
