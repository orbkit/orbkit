import numpy
import re
from itertools import chain
from .symmetry import irrep_mul
from . import cy_mo_integrals

# chemists' notation:
# [ij|kl] = \int\int dr1 dr2 \phi_i^*(r1) \phi_j(r1) 1/r12 \phi_k^*(r2) \phi_l(r2)

class FCIDUMP(object):
  '''Stores 1- and 2-electron integrals.
  '''
  def __init__(self,
        Norb=None, Nelec=None, spin=1,
        OrbSym=None, ISym=1, restricted=True,
    ):
    '''
    **Parameters:**

    Norb : int
      number of orbitals
    Nelec : int
      number of electrons (can probably be ignored?)
    spin : int
      spin multiplicity 2S+1
    OrbSym : list of int
      IRREP assigned to each orbital; len(OrbSym)==Norb; if None all orbitals are assumed to be of IRREP 1
    ISym : int
      symmetry of the state (can probably be ignored?)
    '''

    # general parameters
    self.restricted = restricted
    self.initialized = False

    # system parameters
    self.Norb = Norb
    self.Nelec = Nelec
    self.spin = spin
    self.OrbSym = numpy.array(OrbSym)
    self.ISym = ISym

    if self.Norb:
      self.init_arrays()
    self.nuclear_repulsion = 0

  def set_OrbSym(self, nmopi):
    '''Set IRREPs for orbitals.

    **Parameters:**

    nmopi : list of int
      number of orbitals per IRREP
    '''
    self.OrbSym = list(chain(*[[irrep+1]*nmo for irrep, nmo in enumerate(nmopi)]))
    self.OrbSym = numpy.array(self.OrbSym)

  def init_arrays(self):

    assert self.Norb, 'require number of orbitals \'Norb\''

    if not self.Nelec:
      self.Nelec = self.Norb
    if not numpy.any(self.OrbSym):
      self.OrbSym = numpy.ones(self.Norb, dtype=int)

    self.H = numpy.zeros((self.Norb,self.Norb))
    self.G = numpy.zeros((self.Norb,self.Norb,self.Norb,self.Norb))

    if not self.restricted:
      self.Hbeta = numpy.zeros(self.H.shape)
      self.Gbeta = numpy.zeros(self.G.shape)
      self.Galphabeta = numpy.zeros(self.G.shape)

    self.initialized = True

  def apply_numerical_zero(self, numerical_zero=1e-15):

    self.G[numpy.absolute(self.G) < numerical_zero] = 0.0
    self.H[numpy.absolute(self.H) < numerical_zero] = 0.0

    if not self.restricted:
      self.Gbeta[numpy.absolute(self.Gbeta) < numerical_zero] = 0.0
      self.Galphabeta[numpy.absolute(self.Galphabeta) < numerical_zero] = 0.0
      self.Hbeta[numpy.absolute(self.Hbeta) < numerical_zero] = 0.0

    if numpy.absolute(self.nuclear_repulsion) < numerical_zero:
      self.nuclear_repulsion = 0.0

  def store(self, filename='FCIDUMP', **kwargs):
    '''Write FCIDUMP file.'''
    store(self, filename, **kwargs)

  ########################
  ###  access to data  ###
  ########################

  def getitem(self, indices, spin='alpha'):
    assert self.initialized, 'Matrices not yet initalized'

    if len(indices) == 4:
      i, j, k, l = indices
    elif len(indices) == 2:
      i, j = indices
      k = l = -1

    if k == l == -1:
      if i == j == -1:
        return self.nuclear_repulsion
      return self.get_H(i, j, spin)

    return self.get_G(i, j, k, l, spin)

  def setitem(self, indices, value, spin='alpha'):
    assert self.initialized, 'Matrices not yet initalized'

    # set indices
    if len(indices) == 4:
      i, j, k, l = cy_mo_integrals.normalize_indices(*indices)
    elif len(indices) == 2:
      i, j = cy_mo_integrals.normalize_indices(*indices)
      k = l = -1

    # store integral
    if k == l == -1:
      if i == j == -1:
        self.nuclear_repulsion = value
      else:
        if self.restricted or spin == 'alpha':
          self.H[i,j] = value
        elif spin == 'beta':
          self.Hbeta[i,j] = value
        else:
          raise ValueError('Invalid spin '+spin)
    else:
      if self.restricted or spin == 'alpha':
        self.G[i,j,k,l] = value
      elif spin == 'beta':
        self.Gbeta[i,j,k,l] = value
      elif spin == 'alphabeta':
        self.Galphabeta[i,j,k,l] = value
      else:
        raise ValueError('Invalid spin '+spin)

  def set_H(self, H, irrep, spin='alpha'):
    '''Set Hcore block for given irrep.'''
    # TODO: normalize indices
    assert self.initialized, 'Matrices not yet initalized'
    assert irrep in self.OrbSym, 'no orbital(s) specified for IRREP %s' %irrep

    orbitals = self.OrbSym == irrep
    orbitals = orbitals[:,numpy.newaxis] & orbitals[numpy.newaxis,:]

    if self.restricted or spin == 'alpha':
      self.H[orbitals] = H.flatten()
    else:
      self.Hbeta[orbitals] = H.flatten()

  def set_G(self, V, irrepi, irrepj, irrepk, irrepl, spin='alpha'):
    '''Set ERI block for given irreps.

    **Parameters:**

    V : numpy.ndarray
      array of integrals (chemists notation)
    irrepi, irrepj, irrepk, irrepl : int
      IRREP for indices i, j, k and l respectively
    '''
    # TODO: normalize indices

    assert self.initialized, 'Matrices not yet initalized'
    for irrep in (irrepi, irrepj, irrepk, irrepl):
      assert irrep in self.OrbSym, 'no orbital(s) specified for IRREP %s' %irrep

    orbitalsi = self.OrbSym == irrepi
    orbitalsj = self.OrbSym == irrepj
    orbitalsk = self.OrbSym == irrepk
    orbitalsl = self.OrbSym == irrepl

    orbitals1 = orbitalsi[:,numpy.newaxis] & orbitalsj[numpy.newaxis,:]
    orbitals2 = orbitalsk[:,numpy.newaxis] & orbitalsl[numpy.newaxis,:]
    orbitals = orbitals1[:,:,numpy.newaxis,numpy.newaxis] & orbitals2[numpy.newaxis,numpy.newaxis,:,:]

    if self.restricted or spin == 'alpha':
      self.G[orbitals] = V.flatten()
    elif spin == 'beta':
      self.Gbeta[orbitals] = V.flatten()
    else:
      self.Galphabeta[orbitals] = V.flatten()

  def get_H(self, i=-1, j=-1, spin='alpha', full=False):
    '''One-electron integrals in spatial basis.
    Returns only one matrix element if indices i and j are given.
    '''
    assert spin in ('alpha', 'beta')
    i, j = int(i), int(j)

    if not (i == j == -1):
      # return only one matrix element
      i, j = cy_mo_integrals.normalize_indices(i, j)
      if self.restricted or spin == 'alpha':
        return self.H[i,j]
      return self.Hbeta[i,j]

    # return (full) matrix
    if self.restricted or spin == 'alpha':
      H = self.H.copy()
    else:
      H = self.Hbeta.copy()
    if full:
      tril = numpy.tril_indices(self.Norb, k=-1)
      Hl = H.transpose()
      H[tril] = Hl[tril]
    return H

  def get_h(self, P=-1, Q=-1, full=False):
    '''One-electron integrals in spin orbital basis.
    Returns only one matrix element if indices P and Q are given.
    '''

    if not (P == Q == -1):
      ## return only one matrix element
      # check spin: <alpha beta | alpha beta> = 0
      if P%2 != Q%2:
        return 0
      # map spin orbitals to spatial orbitals using integer division
      p, q = P//2, Q//2
      # get spin
      spin = ('alpha', 'beta')[int(P)%2]
      # return integral
      return self.get_H(p,q, spin)

    # generate (full) matrix
    F = self.get_F()
    if self.restricted:
      Fbeta = F
    else:
      Fbeta = self.get_F(spin='beta')
    f = cy_mo_integrals.transform_F(F, Fbeta)
    if full:
      tril = numpy.tril_indices(2*self.Norb, k=-1)
      fl = f.transpose()
      f[tril] = fl[tril]
    return f

  # create alias
  get_F = get_H
  get_f = get_h

  def get_G(self, i=-1, j=-1, k=-1, l=-1, spin='alpha', full=False):
    '''Two-electron integrals in spatial basis.
    Returns only one matrix element if indices i, j, k and l are given.
    '''
    # (ij|kl) = \int\int dr1 dr2 \phi_i^*(r1) \phi_j(r1) 1/r12 \phi_k^*(r2) \phi_l(r2)
    assert spin in ('alpha', 'beta', 'alphabeta')
    i, j, k, l = int(i), int(j), int(k), int(l)

    if not (i == j == k == l == -1):
      # return only one matrix element
      i, j, k, l = cy_mo_integrals.normalize_indices(i, j, k, l)
      if self.restricted or spin == 'alpha':
        return self.G[i,j,k,l]
      elif spin == 'beta':
        return self.Gbeta[i,j,k,l]
      else:
        return self.Galphabeta[i,j,k,l]

    # return (full) matrix
    if self.restricted or spin == 'alpha':
      g = self.G.copy()
    elif spin == 'beta':
      g = self.Gbeta.copy()
    else:
      g = self.Galphabeta.copy()
    if full:
      cy_mo_integrals.expand_G(g)
    return g

  def get_g(self, P=-1, Q=-1, R=-1, S=-1, full=False):
    '''Two-electron integrals in spin orbital basis.
    Returns only one matrix element if indices P, Q, R and S are given.
    '''

    if not (P == Q == R == S == -1):
      ## return only one matrix element
      # check spin: <alpha beta | alpha beta> = 0
      if (P%2 != Q%2) or (R%2 != S%2):
        return 0
      # map spin orbitals to spatial orbitals using integer division
      p, q = P//2, Q//2
      r, s = R//2, S//2
      # get spin
      spin1 = ('alpha', 'beta')[P%2]
      spin2 = ('alpha', 'beta')[R%2]
      spin = ('alphabeta', spin1)[spin1==spin2]
      # return integral
      return self.get_G(p,q,r,s, spin)

    # generate (full) matrix
    G = self.get_G()
    if self.restricted:
      Gbeta = G
      Galphabeta = G
    else:
      Gbeta = self.get_G(spin='beta')
      Galphabeta = self.get_G(spin='alphabeta')
    g = cy_mo_integrals.transform_G(G, Gbeta, Galphabeta)
    if full:
      cy_mo_integrals.expand_G(g)
    return g

  #############################
  ###  change active space  ###
  #############################

  def reduce_active_space(self, core=0, occ=None):
    '''Reduce to active orbitals.

    **Parameters:**

    occ : int or list of ints
      A single integer selects the first n MOs. Supply a list of integers to select for each IRREP indivdually.
    core : int or list of ints
      Defines core orbitals.
    '''
    assert self.restricted, 'unrestricted not yet implemented'

    if isinstance(core, int):
      core = range(core)
    else:
      tmp = []
      for i, c in enumerate(core):
        OrbsIrrep = numpy.where(self.OrbSym==(i+1))[0]
        if c > OrbsIrrep.size:
          raise IndexError('core=%i for IRREP %i requested, but only %i orbitals are available' %(c, i+1, OrbsIrrep.size))
        tmp.append(OrbsIrrep[:c])
      core = sorted(chain(*tmp))

    if occ is None:
      occ = self.Norb
    if isinstance(occ, int):
      occ = range(occ)
    else:
      tmp = []
      for i, c in enumerate(occ):
        OrbsIrrep = numpy.where(self.OrbSym==(i+1))[0]
        if c > OrbsIrrep.size:
          raise IndexError('occ=%i for IRREP %i requested, but only %i orbitals are available' %(c, i+1, OrbsIrrep.size))
        tmp.append(OrbsIrrep[:c])
      occ = sorted(chain(*tmp))

    not_core = sorted(set(range(self.Norb)).difference(core))
    occ = sorted(set(occ).difference(core).intersection(range(self.Norb)))

    if not occ:
      raise ValueError('No occupied orbitals left')

    # expand to full matrices
    self.G = self.get_G(full=True)

    # index convention:
    # i, j : not_core
    # r, s : core

    # update nuclear repulsion
    # NR += \sum_r H_rr + 2*\sum_rs G_rrss - \sum_rs G_rssr
    self.nuclear_repulsion += 2*numpy.sum(self.H[core,core])
    self.nuclear_repulsion += 2*numpy.sum(self.G[core,core,:,:][:,core,core])
    self.nuclear_repulsion -= numpy.sum(self.G[core,:,:,core][:,core,core])

    # update H
    # H_ij += 2*\sum_r G_ijrr - \sum_r G_irrj
    Hind = numpy.ix_(not_core, not_core)
    self.H[Hind] += 2*numpy.sum(self.G[:,:,core,core][not_core,:,:][:,not_core,:], axis=2)
    self.H[Hind] -= numpy.sum(self.G[:,core,core,:][not_core,:,:][:,:,not_core], axis=1)

    # cut Matrices
    self.H = self.H[numpy.ix_(occ,occ)]
    self.G = self.G[numpy.ix_(occ,occ,occ,occ)]

    if not self.restricted:
      self.Hbeta = self.Hbeta[numpy.ix_(occ,occ)]
      self.Gbeta = self.Gbeta[numpy.ix_(occ,occ,occ,occ)]
      self.Galphabeta = self.Galphabeta[numpy.ix_(occ,occ,occ,occ)]

    # update system parameters
    self.Nelec -= 2*len(core)
    self.Norb = len(occ)
    self.OrbSym = self.OrbSym[occ]

  #################################
  ###  change orbital ordering  ###
  #################################

  def change_order(self, order):
    '''Reorder orbitals by given order.'''
    assert len(set(order)) == len(order), 'duplicate entries in order detected'
    assert max(order) <= self.Norb, 'order contains index exceeding number of orbitals'

    if len(order) < self.Norb:
      # add missing indices in original order
      order.extend(set(range(self.Norb)) - set(order))

    # expand to full matrices
    self.H = self.get_H(full=True)
    self.G = self.get_G(full=True)

    if not self.restricted:
      self.Hbeta = self.get_H(spin='beta', full=True)
      self.Gbeta = self.get_G(spin='beta', full=True)
      self.Galphabeta = self.get_G(spin='alphabeta', full=True)

    # reorder
    self.H = self.H[numpy.ix_(order, order)]
    self.G = self.G[numpy.ix_(order, order, order, order)]

    if not self.restricted:
      self.Hbeta = self.Hbeta[numpy.ix_(order, order)]
      self.Gbeta = self.Gbeta[numpy.ix_(order, order, order, order)]
      self.Galphabeta = self.Galphabeta[numpy.ix_(order, order, order, order)]

    self.OrbSym = self.OrbSym[order]

  def order_irreps(self):
    '''Reorder orbitals in groups of IRREPs'''
    order = []
    for i in set(self.OrbSym):
      order.extend(numpy.where(self.OrbSym == i)[0])
    self.change_order(order)

###############################
###  load/store ASCII file  ###
###############################

def load(filename='FCIDUMP'):
  '''Read Integrals from an ASCII file. Return a FCIDUMP instance.'''

  ## load file
  with open(filename) as fd:
    filedata = fd.readlines()

  ## parse head
  Norb = int(re.search('NORB\s*=\s*([0-9]*)', filedata[0], re.IGNORECASE).group(1))
  Nelec = int(re.search('NELEC\s*=\s*([0-9]*)', filedata[0], re.IGNORECASE).group(1))
  spin = int(re.search('MS2\s*=\s*([0-9]*)', filedata[0], re.IGNORECASE).group(1)) + 1
  restricted = True

  for iline, line in enumerate(filedata):
    if line == ' /\n':
      start = iline + 1
      break
    elif 'ISYM' in line:
      ISym = int(line.strip(' ISYM=,\n'))
    elif 'IUHF' in line:
      restricted = False

  OrbSym = re.search('ORBSYM=([\d,\s]*)ISYM', ''.join(filedata[:start]), re.IGNORECASE).group(1)
  OrbSym = [int(m.group()) for m in re.finditer('\d', OrbSym)]

  ## create Integrals object
  integrals = FCIDUMP(
    Norb=Norb, Nelec=Nelec, spin=spin,
    OrbSym=OrbSym, ISym=ISym,
    restricted=restricted,
  )

  ## parse lines
  spin = 'alpha'
  for line in filedata[start:]:

    fields = line.split()
    integral = float(fields[0])
    i, j, k, l = map(int, fields[1:5])

    if k == l == 0:
      if i == j == 0:
        # nuclear repulsion
        integrals.nuclear_repulsion = integral
        spin = {'alpha':'beta', 'beta':'alphabeta', 'alphabeta':'alpha'}[spin]
      else:
        # Hcore
        integrals.setitem((i-1,j-1), integral, spin)
    else:
      # Gijkl
      integrals.setitem((i-1,j-1,k-1,l-1), integral, spin)

  ## return integrals object
  return integrals

def store(integrals, filename='FCIDUMP', numerical_zero=0, molpro_format=True):
  '''
  Store all integrals in an ASCII file.

  **Parameters**:

  filename : string
  numerical_zero : float
    all values below this threshold will be omitted in file
  molpro_format : boolean
    if True, shift decimal point from scientific notation one to the left, so that all numbers start with 0
  '''

  assert integrals.initialized, 'Matrices not initalized'

  def format_line(I_ijkl, i=-1, j=-1, k=-1, l=-1):

    I_ijkl = ' % 19.16E' %(I_ijkl)
    if molpro_format:
      # shift dot and increment exponent
      mantisse = '%18.16f' %(float(I_ijkl[2:20])/10)
      exp = '%0+3d' %(int(I_ijkl[21:24])+1)
      I_ijkl = I_ijkl[:2] + mantisse + 'E' + exp + I_ijkl[24:]

    return I_ijkl + '%4d%4d%4d%4d' %(i+1, j+1, k+1, l+1)

  with open(filename, 'w') as fd:

    ## head
    fd.write(' &FCI NORB=%3d,NELEC=%3d,MS2=%2d,\n' %(integrals.Norb, integrals.Nelec, integrals.spin-1))
    fd.write('  ORBSYM=%s,\n' %(','.join(map(str, integrals.OrbSym))))
    fd.write('  ISYM=%s,\n' %integrals.ISym)
    if not integrals.restricted:
      fd.write('  IUHF=1,\n')
    fd.write(' /')

    ## write two-electron integrals
    for spin in ('alpha', 'beta', 'alphabeta'):

      for i in range(integrals.Norb):
        for j in range(integrals.Norb):

          # for real basis functions
          if i < j:
            continue

          for k in range(integrals.Norb):
            for l in range(integrals.Norb):

              # for real basis functions
              if k < l:
                continue

              # exchange of electronic coordinates
              if spin != 'alphabeta':
                if (i < k) or (i == k and j < l):
                  continue

              # check symmetry
              if irrep_mul(*integrals.OrbSym[[i,j,k,l]]) != 1:
                 continue

              Vijkl = integrals.getitem((i,j,k,l), spin)
              if not numerical_zero or abs(Vijkl) > numerical_zero:
                fd.write('\n' + format_line(Vijkl, i, j, k, l))

      if integrals.restricted:
        break
      else:
        fd.write('\n' + format_line(0, -1, -1, -1, -1))

    ## write one-electron integrals
    for spin in ('alpha', 'beta'):

      for i in range(integrals.Norb):
        for j in numpy.where(integrals.OrbSym==integrals.OrbSym[i])[0]:

          # for real basis functions
          if i < j:
            continue

          tij = integrals.getitem((i,j), spin)
          if not numerical_zero or abs(tij) > numerical_zero:
            fd.write('\n' + format_line(tij, i, j))

      if integrals.restricted:
        break
      else:
        fd.write('\n' + format_line(0, -1, -1, -1, -1))

    ## write nuclear repulsion
    fd.write('\n' + format_line(integrals.nuclear_repulsion))
