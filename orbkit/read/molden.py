from collections import defaultdict
import numpy
import re

from orbkit.qcinfo import QCinfo
from orbkit.orbitals import AOClass, MOClass
from orbkit.display import display
from orbkit.tools import l_deg, lquant, orbit

from .tools import descriptor_from_file

'''
New Molden interface
'''

### define regular expressions
# re.I = re.IGNORECASE
regex_molden = re.compile(r'\[\s*molden\s+format\s*\]', re.I)
re_float = r'([-+]?\d+\.?\d*[ed]?[+-]?\d+|NaN)'

# [Atoms] section (geo_info)
regex_atoms = re.compile(r'\[atoms\]\s*\(?(angs|au)\)?', re.I)
regex_atom = re.compile(r'\s*([a-z]+)\s+(\d+)\s+(\d+)\s+' + r'\s+'.join((re_float,)*3), re.I)

# [GTO] section (ao_info)
regex_basis = re.compile(r'\s*(\d+)\s+(\d+)$', re.I)
regex_contraction = re.compile(r'\s*([a-z]+)\s+(\d+)\s+(\d+(\.\d+)?)($)', re.I)
regex_primitive = re.compile(r'\s*'+r'\s+'.join((re_float,)*2), re.I)

# flags for use of spherical/cartesian basis functions
FLAGS_SPH  = ['5d',  '7f',  '9g']
FLAGS_CART = ['6d', '10f', '15g']

# regex to match all lines like "[5d7f]"
regex_flagline = re.compile(r'\[((' + '|'.join(FLAGS_SPH + FLAGS_CART) + r')+)\]', re.I)
regex_flag = re.compile('(\d+[dfg])', re.I)

# [MO] section (mo_info)
regex_sym = re.compile(r'\s*sym\s*=\s*(\S+)', re.I)
regex_energy = re.compile(r'\s*ene(?:rgy)?\s*=\s*' + re_float, re.I)
regex_spin = re.compile(r'\s*spin\s*=\s*(alpha|beta)', re.I)
regex_occu = re.compile(r'\s*occup\s*=\s*' + re_float, re.I)
regex_coeff = re.compile(r'\s*(\d+)\s+' + re_float, re.I)


def read_molden(fname, all_mo=False, spin=None, i_md=-1, interactive=True,
                **kwargs):
  '''Reads all information desired from a molden file.

  **Parameters:**

    fname : str, file descriptor
      Specifies the filename for the input file.
      fname can also be used with a file descriptor instad of a filename.
    all_mo : bool, optional
      If True, all molecular orbitals are returned.
    spin : {None, 'alpha', or 'beta'}, optional
      If not None, returns exclusively 'alpha' or 'beta' molecular orbitals.
    i_md : int, default=-1
      Selects the `[Molden Format]` section of the output file.
    interactive : bool
      If True, the user is asked to select the different sets.

  **Returns:**

    qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, mo_spec, etot :
        See :ref:`Central Variables` for details.
  '''

  if 'index' not in kwargs.keys():
    kwargs['index'] = 0

  if isinstance(fname, str):
    fd = descriptor_from_file(fname, index=kwargs['index'])
  else:
    fd = fname
    fname = fd.name

  ### read the whole file into RAM
  # TODO: optimize for large files

  molden = fd.read()
  if isinstance(molden, bytes):
    molden = molden.decode()

  ### find number of [Molden Format] entries and figure our which one to use
  entries = [m.start() for m in regex_molden.finditer(molden)]
  count = len(entries)

  if count == 0:
    raise IOError('The input file {:s} is no valid molden file!\n\nIt does'.format(fname) +
            ' not contain the keyword: [Molden Format]\n')

  if count > 1:
    display('\nContent of the molden file:')
    display('\tFound {:d} [Molden Format] keywords, i.e., '.format(count) +
            'this file contains {:d} molden files.'.format(count))

    if interactive:
      message = '\tPlease give an integer from 0 to {0}: '.format(count-1)
      from builtins import input  # Python2 compatibility

      while 1:
        try:
          i_md = int(input(message))
        except ValueError:
          print('An Integer is required!')
        else:
          if i_md >= count or i_md < -count:
              # invalid index
              continue
          break

    i_md = list(range(count))[i_md]

    # log selected index
    display('\tSelecting the element with index {:d}.'.format(i_md))

    # select molden entry
    start = entries[i_md]
    end = (entries + [None])[i_md+1]
    molden = molden[start:end]

  molden = molden.splitlines()

  ### parse [Atoms] and [GTO] section
  qc = QCinfo()
  has_alpha = False
  has_beta = False
  restricted = False
  spherical_basis = []  # found flags for spherical basis
  cartesian_basis = []  # found flags for cartesian basis
  angular = []          # angular momentum actually used
  by_orca = False

  for iline, line in enumerate(molden):

    if 'orca' in line.lower():
      by_orca = True
      continue

    if '_ENERGY=' in line:
      try:
        qc.etot = float(line.split()[1])
      except IndexError:
        pass
      continue

    # [Atoms] section (geo_info)
    m = regex_atoms.match(line)
    if m:
      angstrom = 'angs' == m.group(1).lower()
      continue

    m = regex_atom.match(line)
    if m:
      qc.geo_info.append(list(m.groups()[:3]))
      qc.geo_spec.append([float(f) for f in m.groups()[3:]])
      continue

    # [GTO] section (ao_info)
    if '[sto]' in line.lower():
      # orbkit does not support Slater type orbitals
      raise IOError('orbkit does not work for STOs!\nEXIT\n')

    m = regex_basis.match(line)
    if m:
      at_num = int(m.group(1)) - 1
      #ao_num = 0
      continue

    # check spherical/cartesian flags
    m = regex_flagline.match(line.lower())
    if m:
      # get list of all flags in line
      flags = regex_flag.findall(m.group(1))
      # check whether cartesian or spherical
      for flag in flags:
        if flag in FLAGS_SPH:
          spherical_basis.append(flag)
        if flag in FLAGS_CART:
          cartesian_basis.append(flag)

    m = regex_contraction.match(line)
    if m:
      ao_num = 0               # Initialize number of atomic orbitals
      ao_type = m.group(1).lower()    # angular momentum
      pnum = int(m.group(2))          # Number of primatives

      for l in ao_type:
        qc.ao_spec.append({'atom': at_num,
                          'type': l,
                          'pnum': -pnum if by_orca else pnum,
                          'coeffs': numpy.zeros((pnum, 2))
                          })
        if not l in angular:
            angular.append(l)
      continue

    m = regex_primitive.match(line)
    if m:
      # split line as regex only captures the first two floats, and there may be more
      coeffs = numpy.array(line.lower().replace('d','e').split(), dtype=numpy.float64)
      for i_ao in range(len(ao_type)):
        qc.ao_spec[-len(ao_type)+i_ao]['coeffs'][ao_num,:] = [coeffs[0],
                                                              coeffs[1+i_ao]]
      ao_num += 1
      continue

    if '[mo]' in line.lower():
      break

  ### checks for cartesion/spherical basis

  # check for mixed spherical/cartesian basis functions
  max_l = max(lquant[l] for l in angular)
  if max_l >= 2:
    # remove flags for unused angular momentum
    l = orbit[2:max_l+1]
    sph = [f for f in spherical_basis if f[-1] in l]
    cart = [f for f in cartesian_basis if f[-1] in l]
    if sph and cart:
      raise IOError('''The input file {} contains mixed spherical and Cartesian function ({}).
                  ORBKIT does not support these basis functions yet.
                  Pleas contact us, if you need this feature!'''.format(
                      fname, ', '.join(sph+cart)))

    # check for ambiguous spherical/cartesian flags
    sph = [l[-1] for l in sph]
    cart = [l[-1] for l in cart]
    if set(sph) & set(cart):
        raise IOError('The input file {} contains ambiguous flags for spherical and cartesian basis functions: {}'.format(fname, ', '.join(spherical_basis+cartesian_basis)))

    cartesian = not bool(sph)

  else:
    cartesian = True    # does not matter for s and p orbitals

  # count number of basis functions
  basis_count = 0
  for AO in qc.ao_spec:
      l = AO['type']
      # TODO: check for mixed sph/cart basis
      basis_count += l_deg(lquant[l], cartesian_basis=cartesian)

  ### parse [MO] section (mo_info)
  newMO = False
  MO_sym = None
  MO_spin = None
  MO_energy = None
  MO_occ = None
  sym = defaultdict(int)     # counter for MOs per IRREP

  for line in molden[iline:]:

    m = regex_coeff.match(line)
    if m:

      if newMO:

        # infer incomplete data
        MO_spin = MO_spin or 'alpha'
        m2 = re.search(r'\d+', MO_sym)
        if m2:
            a = m2.group()
            if MO_sym == a:
                MO_sym = '{:s}.1'.format(a)
            elif MO_sym.startswith(a):
                MO_sym.replace(a, '{:s}.'.format(a), 1)
            else:
                sym[a] += 1
                MO_sym = '{:d}.{:s}'.format(sym[a], MO_sym)
        MO_sym = MO_sym or '%d.1' % (len(qc.mo_spec)+1)

        # create a new MO entry
        qc.mo_spec.append({'coeffs': numpy.zeros(basis_count),
                           'sym' : MO_sym,
                           'energy' : MO_energy,
                           'occ_num' : MO_occ,
                           'spin' : MO_spin,
                        })

        # reset variables
        newMO = False
        MO_sym = None
        MO_spin = None
        MO_energy = None
        MO_occ = None

      # parse and store current coefficient
      iMO = int(m.group(1)) - 1
      coeff = float(m.group(2))
      if numpy.isnan(coeff):
        display('Warning: coefficient {:d} of MO {:s} is NaN! Using zero instead'.format(
            iMO, qc.mo_spec[-1]['sym']))
      else:
        qc.mo_spec[-1]['coeffs'][iMO] = coeff
      continue

    newMO = True
    m = regex_sym.match(line)
    if m:
      MO_sym = m.group(1)
      continue

    m = regex_energy.match(line)
    if m:
      MO_energy = m.group(1)
      continue

    m = regex_spin.match(line)
    if m:
      MO_spin = m.group(1).lower()
      has_alpha = has_alpha or MO_spin == 'alpha'
      has_beta = has_beta or MO_spin == 'beta'
      continue

    m = regex_occu.match(line)
    if m:
      MO_occ = float(m.group(1))
      restricted = restricted or (MO_occ > 1.0001)
      continue

  ### post checks and clean up

  if spin is not None:
    if restricted:
      raise IOError('The keyword `spin` is only supported for unrestricted calculations.')
    if spin != 'alpha' and spin != 'beta':
      raise IOError('`spin=%s` is not a valid option' % spin)
    elif spin == 'alpha' and has_alpha:
      display('Reading only molecular orbitals of spin alpha.')
    elif spin == 'beta' and has_beta:
      display('Reading only molecular orbitals of spin beta.')
    elif (not has_alpha) and (not has_beta):
      raise IOError(
           'Molecular orbitals in `molden` file do not contain `Spin=` keyword')
    elif ((spin == 'alpha' and not has_alpha) or
          (spin == 'beta' and not has_beta)):
      raise IOError('You requested `%s` orbitals, but None of them are present.' % spin)

  # Spherical basis?
  if spherical_basis:
    qc.ao_spec.set_lm_dict(p=[1,0])

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

  # Orca uses for all molecular orbitals the same name
  sym = [i['sym'] for i in qc.mo_spec]
  if sym[1:] == sym[:-1]:
    sym = sym[0].split('.')[-1]
    for i in range(len(qc.mo_spec)):
      qc.mo_spec[i]['sym'] = '%d.%s' % (i+1,sym)

  # Convert geo_info and geo_spec to numpy.ndarrays
  qc.format_geo(is_angstrom=angstrom)

  # Check the normalization
  from orbkit.analytical_integrals import get_ao_overlap
  spher_tmp = qc.ao_spec.spherical
  qc.ao_spec.spherical = False
  norm = numpy.diagonal(get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec))
  qc.ao_spec.spherical = spher_tmp
  if max(numpy.abs(norm-1.)) > 1e-5:
    display('The atomic orbitals are not normalized correctly, renormalizing...\n')
    if not by_orca:
      j = 0
      for i in range(len(qc.ao_spec)):
        qc.ao_spec[i]['coeffs'][:,1] /= numpy.sqrt(norm[j])
        for n in range(l_deg(lquant[qc.ao_spec[i]['type']],cartesian_basis=True)):
          j += 1
    else:
      qc.ao_spec[0]['N'] = 1/numpy.sqrt(norm[:,numpy.newaxis])

    if cartesian_basis:
      from orbkit.cy_overlap import ommited_cca_norm
      cca = ommited_cca_norm(qc.ao_spec.get_lxlylz())
      for mo in qc.mo_spec:
        mo['coeffs'] *= cca

  qc.mo_spec.update()
  qc.ao_spec.update()
  return qc
