from orbkit.qcinfo import QCinfo
from orbkit.orbitals import AOClass, MOClass
from orbkit.tools import *
import numpy

from .tools import descriptor_from_file

def read_wfn(fname, all_mo=False, spin=None, **kwargs):
  '''Reads all information desired from a wfn file.
  
  **Parameters:**
  
    fname: str, file descriptor
      Specifies the filename for the input file.
      fname can also be used with a file descriptor instad of a filename.
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
  qc.ao_spec = AOClass([])
  qc.mo_spec = MOClass([])
  sec_flag = None                 # A Flag specifying the current section
  is_wfn = False                  # Check type of file
  ao_num = 0                      # Number of AO
  mo_num = 0                      # Number of MO
  at_num = 0                      # Number of atoms
  c_type = 0                      # Counting variable for AO type
  c_exp = 0                       # Counting variable for AO exponents
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

  for line in flines:
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
              'lxlylz': None,
              #'lm': None
              })
    elif 'TYPE ASSIGNMENTS' in line:
      thisline = line[18:].split()
      for i in range(len(thisline)):
        qc.ao_spec[c_type]['lxlylz'] = lxlylz[int(thisline[i])-1][numpy.newaxis]
        qc.ao_spec[c_type]['type'] = orbit[sum(lxlylz[int(i)-1])]
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
 
  if isinstance(fname, str):
    fname.close()                    # Leave existing file descriptors alive
  
  # Remove numbers from atom names
  for i in qc.geo_info:
    i[0] = ''.join([k for k in i[0] if not k.isdigit()])
  # Convert geo_info and geo_spec to numpy.ndarrays
  qc.format_geo(is_angstrom=False)
  
  qc.mo_spec.update()
  qc.ao_spec.update()
  return qc
