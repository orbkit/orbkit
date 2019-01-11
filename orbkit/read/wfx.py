import numpy

from orbkit.qcinfo import QCinfo
from orbkit.orbitals import AOClass, MOClass
from orbkit.tools import orbit, exp_wfn

from .tools import descriptor_from_file, spin_check

def read_wfx(fname, all_mo=False, spin=None, **kwargs):
  '''Reads all information desired from a wfn file.
  
  **Parameters:**
  
    fname: str, file descriptor
      Specifies the filename for the input file.
      fname can also be used with a file descriptor instad of a filename.
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
  qc.ao_spec = AOClass([])
  qc.mo_spec = MOClass([])
  lxlylz = []
  for j in exp_wfn:         
    lxlylz.extend(j)
  lxlylz = numpy.array(lxlylz,dtype=numpy.int64)
  
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
      qc.ao_spec = AOClass([{'atom': None,
                    'pnum': -1,
                    'coeffs': None,
                    'lxlylz': None,
                    #'lm': None
                    } for i in range(ao_num)])
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
      qc.mo_spec = MOClass([{'coeffs': numpy.zeros(ao_num),
                     'energy': None,
                     'occ_num': None,
                     'spin': None,
                     'sym': '%s.1' % (i+1)
                    } for i in range(mo_num)])
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
          qc.ao_spec[count]['lxlylz'] = lxlylz[int(i)-1][numpy.newaxis]
          qc.ao_spec[count]['type'] = orbit[sum(lxlylz[int(i)-1])]
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

  qc.mo_spec.update()
  qc.ao_spec.update()
  return qc
