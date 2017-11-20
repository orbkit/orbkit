from copy import copy
import numpy
import sys

from orbkit.display import display
from orbkit.qcinfo import CIinfo
from orbkit.read.tools import descriptor_from_file

from .tools import point_groups, multiplicity

def psi4_detci(fname,select_run=None,threshold=0.0,**kwargs):
  '''Reads PSI4 DETCI output. 
  
  **Parameters:**
  
    fname: str, file descriptor
      Specifies the filename for the input file.
      fname can also be used with a file descriptor instad of a filename.
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

  for line in flines:
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
    string = ', '.join(map(str,numci)).rsplit(', ',1)
    string = ' and '.join(string) if len(string[0].split(',')) < 2 else ', and '.join(string)
    display('The input file %s contains' % (filename))
    display('%d DETCI calculation(s) with %s root(s)%s.'%(count, string,
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
        if not isinstance(i,(int, numpy.int, numpy.int64, numpy.intc)):
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
  synonyms = {'dropped docc': 'core',
              'frozen docc': 'core',
              'restricted docc': 'closed',
              'ras ': 'active',
              'active': 'active',
              'restricted uocc': 'external',
              'frozen uocc': 'external',
              'dropped uocc': 'external'  }
  irreps = []
  fist_active_mo = []
  index_active_mo = []
  occ_symbols = {'A': numpy.array([1,0]),
                 'B': numpy.array([0,1]),
                 'X': numpy.array([1,1])}

  flines = iter(flines)
  for line in flines:
    thisline = line.split()             # The current line split into segments
    
    if 'Running in ' in line and 'symmetry.' in line: 
      nIRREP = point_groups()[thisline[2].lower()]
      rhf_occ = numpy.zeros(nIRREP,dtype=numpy.intc)
    #--- RHF occupation
    elif 'Final Occupation by Irrep:' in line:
      if sys.version_info.major == 2:
        irreps = next(flines).split()
        line = next(flines)
      else:
        irreps = next(flines).split()
        line = next(flines)
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
        info['sym'] = thisline[-1]
        if info['sym'] != 'auto':
          info['sym'] = irreps[int(info['sym'])]
      elif ' S     ' in line and ' =' in  line:
        info['spin'] = multiplicity()[int(2*float(thisline[2])+1)]
      elif 'Electrons    =' in line:
        info['nel'] = int(thisline[-1])
      elif 'ORBS IN CI   =' in line: 
        orbs_in_ci = int(thisline[-1])
      elif  any([i in line.lower() for i in synonyms.keys()]):
        for i in synonyms.keys():
          if i in line.lower():
            occ_info[synonyms[i]] += numpy.array(thisline[-nIRREP:],
                                                 dtype=numpy.intc)
      elif ('* ROOT' in line and 'CI total energy' in line) or 'CI Root' in line and 'energy' in line:
        if '*' in line:
          state = '%s.%s'%(thisline[2],info['sym'])
          energy = float(thisline[7])
        else:
          state = '%s.%s'%(int(thisline[2])+1,info['sym'])
          energy = float(thisline[-1])
        ci[index_run].append(CIinfo(method=method))
        ci[index_run][-1].info = copy(general_information) 
        ci[index_run][-1].info['fileinfo'] += '@%d' % index_run
        ci[index_run][-1].info['irreps'] = irreps
        ci[index_run][-1].info['state'] = state
        ci[index_run][-1].info['energy'] = energy
        ci[index_run][-1].info['spin'] = info['spin']
        ci[index_run][-1].info['nel'] = sum(rhf_occ)
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
        if sys.version_info.major == 2:
          next(flines)
        else:
          next(flines)
        def rm(s,w='*(,)'):
          for i in w:
            s = s.replace(i,' ')
          return s
        for i in range(num_det):  
          if sys.version_info.major == 2:
            thisline = rm(next(flines)).split()
          else:
            thisline = rm(next(flines)).split()
          if thisline == []:
            break
          min_c = min(min_c,abs(float(thisline[1])))
          if abs(float(thisline[1])) > threshold:
            ci[index_run][-1].coeffs.append(float(thisline[1]))
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
