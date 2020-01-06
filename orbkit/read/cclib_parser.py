import numpy

from orbkit.qcinfo import QCinfo
from orbkit.orbitals import AOClass, MOClass
from orbkit.units import aa_to_a0, ev_to_ha
from orbkit.display import display
from orbkit.tools import l_deg, lquant, get_atom_symbol
from importlib import import_module

def read_with_cclib(filename, cclib_parser=None, all_mo=False, spin=None, 
                    **kwargs):
  '''Reads all information desired using cclib.
  
  **Parameters:**
  
    filename : str
      Specifies the filename for the input file.
    cclib_parser : str
      If itype is 'cclib', specifies the cclib.parser.
    all_mo : bool, optional
      If True, all molecular orbitals are returned.
    spin : {None, 'alpha', or 'beta'}, optional
      If not None, returns exclusively 'alpha' or 'beta' molecular orbitals.
  
  **Returns:**
  
    qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, mo_spec, etot :
        See :ref:`Central Variables` for details.
  '''
  
  #Maybe we actually don't need this
  #Can someone check if cclib can handle
  #file descriptors?
  assert isinstance(filename, str)

  if not isinstance(cclib_parser,str):
    raise IOError('cclib requires the specification of parser, e.g., ' + 
                  'cclib_parser="Gaussian".')
  
  if cclib_parser == 'Molpro':
    display('\nThe Molpro basis set is not properly read by the cclib parser.')
    display('Please create a molden file with Molpro, i.e., ' + 
            '\n\tput,molden,output.molden,NEW;\n')
  
  parsedic = {'Gaussian': 'gaussianparser', 'Gamess': 'gamessparser',
              'Orca': 'orcaparser'}
  module = import_module('cclib.parser.{}'.format(parsedic[cclib_parser]))
  if cclib_parser != 'Gaussian':
    cclib_parser = cclib_parser.upper()
  parser = getattr(module,cclib_parser)(filename)
  ccData = parser.parse()
  return convert_cclib(ccData, all_mo=all_mo, spin=spin)

def convert_cclib(ccData, all_mo=False, spin=None):
  '''Converts a ccData class created by cclib to an instance of
  orbkit's QCinfo class.

  **Parameters:**
  
    ccData : class
      Contains the input data created by cclib.
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
  
  # Converting all information concerning atoms and geometry
  qc.geo_spec = ccData.atomcoords[0] * aa_to_a0
  for ii in range(ccData.natom):
    symbol = get_atom_symbol(atom=ccData.atomnos[ii])
    qc.geo_info.append([symbol,str(ii+1),str(ccData.atomnos[ii])])
  
  # Convert geo_info and geo_spec to numpy.ndarrays
  qc.format_geo()
  
  # Converting all information about atomic basis set

  for ii in range(ccData.natom):
    for jj in range(len(ccData.gbasis[ii])):
      pnum = len(ccData.gbasis[ii][jj][1])
      qc.ao_spec.append({'atom': ii,
                  'type': str(ccData.gbasis[ii][jj][0]).lower(),
                  'pnum':  pnum,
                  'coeffs': numpy.zeros((pnum, 2))
                  })
      for kk in range(pnum):
        qc.ao_spec[-1]['coeffs'][kk][0] = ccData.gbasis[ii][jj][1][kk][0]
        qc.ao_spec[-1]['coeffs'][kk][1] = ccData.gbasis[ii][jj][1][kk][1]
  
  if hasattr(ccData,'aonames'):
    # Reconstruct exponents list for ao_spec
    cartesian_basis = True
    for i in ccData.aonames:
      if '+' in i or '-' in i:
        cartesian_basis = False

    if not cartesian_basis:
        qc.ao_spec.spherical = True
    
    count = 0
    for i,ao in enumerate(qc.ao_spec):
      l = l_deg(lquant[ao['type']],cartesian_basis=cartesian_basis)
      if cartesian_basis:
        ao['lxlylz'] = []
      else:
        ao['lm'] = []
      for ll in range(l):
        if cartesian_basis:
          ao['lxlylz'].append((ccData.aonames[count].lower().count('x'),
                               ccData.aonames[count].lower().count('y'),
                               ccData.aonames[count].lower().count('z')))
        else:
          m = ccData.aonames[count].lower().split('_')[-1]
          m = m.replace('+',' +').replace('-',' -').replace('s','s 0').split(' ') 
          p = 'yzx'.find(m[0][-1])
          if p != -1:
            m = p - 1
          else:
            m = int(m[-1])
          ao['lm'].append((lquant[ao['type']],m))
        count += 1
  
  # Converting all information about molecular orbitals
  ele_num = numpy.sum(ccData.atomnos) - numpy.sum(ccData.coreelectrons) - ccData.charge
  ue = (ccData.mult-1)
  
  # Check for natural orbitals and occupation numbers
  is_natorb = False
  if hasattr(ccData,'nocoeffs'):
    if not hasattr(ccData,'nooccnos'):
      raise IOError('There are natural orbital coefficients (`nocoeffs`) in the cclib' + 
                    ' ccData, but no natural occupation numbers (`nooccnos`)!')
    is_natorb = True
  
  restricted = (len(ccData.mosyms) == 1)
  if spin is not None:
    if spin != 'alpha' and spin != 'beta':
      raise IOError('`spin=%s` is not a valid option' % spin)
    elif restricted:
      raise IOError('The keyword `spin` is only supported for unrestricted calculations.')
    else:
      qc.mo_spec.spinpola
      display('Converting only molecular orbitals of spin %s.' % spin)
  
  sym = {}
  if len(ccData.mosyms) == 1:
    add = ['']
    orb_sym = [None]
  else:
    add = ['_a','_b']      
    orb_sym = ['alpha','beta']
  
  nmo = ccData.nmo if hasattr(ccData,'nmo') else len(ccData.mocoeffs[0])  
  for ii in range(nmo):    
    for i,j in enumerate(add):
      a = '%s%s' % (ccData.mosyms[i][ii],j)
      if a not in sym.keys(): sym[a] = 1
      else: sym[a] += 1
      if is_natorb:
        occ_num = ccData.nooccnos[ii]
      elif not restricted:
        occ_num = 1.0 if ii <= ccData.homos[i] else 0.0
      elif ele_num > ue:
        occ_num = 2.0
        ele_num -= 2.0
      elif ele_num > 0.0 and ele_num <= ue: 
        occ_num = 1.0
        ele_num -= 1.0
        ue -= 1.0
      else:
        occ_num = 0.0
        
      qc.mo_spec.append({'coeffs': (ccData.nocoeffs if is_natorb else ccData.mocoeffs[i])[ii],
              'energy': 0.0 if is_natorb else ccData.moenergies[i][ii]*ev_to_ha,
              'occ_num': occ_num,
              'sym': '%d.%s' %(sym[a],a)
              })
      if orb_sym[i] is not None:
        qc.mo_spec[-1]['spin'] = orb_sym[i]
        if spin is not None and spin != orb_sym[i]:
          del qc.mo_spec[-1]
  
  # Use default order for atomic basis functions if aonames is not present
  if not hasattr(ccData,'aonames'):
    display('The attribute `aonames` is not present in the parsed data.')
    display('Using the default order of basis functions.')
    
    # Check which basis functions have been used
    c_cart = sum([l_deg(l=ao['type'], cartesian_basis=True) for ao in qc.ao_spec])
    c_sph = sum([l_deg(l=ao['type'], cartesian_basis=False) for ao in qc.ao_spec])
    
    c = qc.mo_spec.get_coeffs().shape[-1]
    if c != c_cart and c == c_sph: # Spherical basis
      qc.ao_spec.set_lm_dict(p=[0,1])
    elif c != c_cart:
      display('Warning: The basis set type does not match with pure spherical ' +
              'or pure Cartesian basis!') 
      display('Please specify qc.ao_spec["lxlylz"] and/or qc.ao_spec["lm"] by your self.')
  
  # Are all MOs requested for the calculation? 
  if not all_mo:
    for i in range(len(qc.mo_spec))[::-1]:
      if qc.mo_spec[i]['occ_num'] < 0.0000001:
        del qc.mo_spec[i]

  qc.mo_spec.update()
  qc.ao_spec.update()
  return qc
