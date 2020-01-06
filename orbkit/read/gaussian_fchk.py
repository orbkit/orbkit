import numpy
from copy import deepcopy

from orbkit.qcinfo import QCinfo
from orbkit.orbitals import AOClass, MOClass
from orbkit.display import display
from orbkit.tools import orbit, lquant

from .tools import descriptor_from_file

def read_gaussian_fchk(fname, all_mo=False, spin=None, **kwargs):
  '''Reads all information desired from a Gaussian FChk file. 

  **Parameters:**
  
  fname : str, file descriptor
    Specifies the filename for the input file.
    fname can also be used with a file descriptor instad of a filename.
  all_mo : bool, optional
      If True, all molecular orbitals are returned.
  
  **Returns:**
  
    qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, mo_spec, etot :
        See :ref:`Central Variables` for details.
  '''
  
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
  ao_sp_coeffs = {}
  switch = 0
  qc = QCinfo()
  qc.geo_info = [[],[],[]]

  if not cartesian_basis:
    qc.ao_spec.spherical = True
  
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
          if not cartesian_basis:
            qc.ao_spec[-1]['lm'] = []
    elif 'Number of primitives per shell' in line:
      sec_flag = 'ao_info'
      index = 'pnum'
      ao_num = int(thisline[-1])
      count = 0
      if qc.ao_spec == []:
        for ii in range(ao_num):
          qc.ao_spec.append({})
          if not cartesian_basis:
            qc.ao_spec[-1]['lm'] = []
    elif 'Shell to atom map' in line:
      sec_flag = 'ao_info'
      index = 'atom'
      ao_num = int(thisline[-1])
      count = 0
      if qc.ao_spec == []:
        for ii in range(ao_num):
          qc.ao_spec.append({})
          if not cartesian_basis:
            qc.ao_spec[-1]['lm'] = []
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
      if 'P(S=P)' not in line:
        sec_flag = 'ao_coeffs'  
      else:
        sec_flag = 'ao_sp_coeffs'  
        ao_sp_coeffs = {0: []}
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
      mo_i0[thisline[0].lower()] = len(qc.mo_spec)
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
        qc.mo_spec.append({'coeffs': numpy.zeros(basis_number),
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
          if index == 'type':         
            ii = orbit[abs(ii)]
            l = lquant[ii]
            if not cartesian_basis:
              for m in (range(0,l+1) if l != 1 else [1,0]):
                #Used to be append so last one should be fine
                qc.ao_spec[count]['lm'].append((l,m))
                if m != 0:
                  qc.ao_spec[count]['lm'].append((l,-m))
          elif index == 'atom':
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
      elif sec_flag == 'ao_sp_coeffs':
        for ii in thisline:
          ao_sp_coeffs[index].append(float(ii))
          count += 1
          ao_num -= 1
          if count == qc.ao_spec[index]['pnum']:
            index += 1
            ao_sp_coeffs[index] = []
            count = 0
        if not ao_num:
          sec_flag = None
      elif sec_flag == 'mo_eorb':
        for ii in thisline:
          qc.mo_spec[count]['energy'] = float(ii)
          count += 1
          if index != 0 and not count % basis_number:
            sec_flag = None
      elif sec_flag == 'mo_coeffs':
        for ii in thisline:    
          qc.mo_spec[mo_i0[what]+index]['coeffs'][count] = float(ii)
          count += 1
          if count == basis_number:
            count = 0
            index += 1
          if index != 0 and not index % basis_number:
            sec_flag = None
  
  # Look for SP atomic orbitals
  if ao_sp_coeffs:
    ao_new = []
    for i,ao in enumerate(qc.ao_spec):
      if ao['type'] == 'p' and sum(numpy.abs(ao_sp_coeffs[i])) > 0:
        ao_new.append(deepcopy(ao))
        ao_new[-1]['type'] = 's'
        ao_new.append(ao)
        ao_new[-1]['type'] = 'p'
        ao_new[-1]['coeffs'][:,1] = numpy.array(ao_sp_coeffs[i])        
      else:
        ao_new.append(ao)
    qc.ao_spec = ao_new   
    
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
  qc.format_geo(is_angstrom=False)
  
  qc.mo_spec = MOClass(qc.mo_spec)
  qc.mo_spec.update()
  qc.ao_spec.update()
  return qc
