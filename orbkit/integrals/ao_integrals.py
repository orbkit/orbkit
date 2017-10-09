import numpy
import ctypes

from ..tools import lquant
from itertools import chain

try:
  from numpy import moveaxis
except:
  from orbkit.tools import moveaxis

# search and load libcint from $PATH
import os
libcint = None
for dirname in os.environ['PATH'].split(':'):
  if os.path.isfile(os.path.join(dirname, 'libcint.so')):
    libcint = ctypes.cdll.LoadLibrary(os.path.abspath(os.path.join(dirname, 'libcint.so')))
    break
if libcint is None:
  raise ImportError("libcint not found: please add to PATH environment variable")

# set return type
libcint.CINTgto_norm.restype = ctypes.c_double

class CINTOpt(ctypes.Structure):
  _fields_ = [
    ('index_xyz_array', ctypes.POINTER(ctypes.POINTER(ctypes.c_int))),
    ('prim_offset', ctypes.POINTER(ctypes.c_int)),
    ('non0ctr', ctypes.POINTER(ctypes.c_int)),
    ('non0idx', ctypes.POINTER(ctypes.POINTER(ctypes.c_int))),
    ('non0coeff', ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),
    ('expij', ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),
    ('rij', ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),
    ('cceij', ctypes.POINTER(ctypes.POINTER(ctypes.c_int))),
    ('tot_prim', ctypes.c_int),
  ]

############################
###  external interface  ###
############################

# TODO:
# - reordering for spherical
# - non-hermitian operators (?)
# - operators including gradients (return tensors)

class AOIntegrals():
  '''Interface to calculate AO Integrals with libcint.
  https://github.com/sunqm/libcint
  '''
  def __init__(self, qc, cartesian=True):

    # general parameters
    self.qc = qc
    self.clear_all_blocks()

    # libcint parameters
    ienv = 0
    self.env = [0]*ienv
    self.atm = []
    self.bas = []
    self.cache_norm = {}  # cache renormalization factors for cartesian integrals

    # build atm
    self.env.extend(qc.geo_spec.flatten())
    for center in qc.geo_info:
      # tuple for each center: (charge, coords offset, 0, 0, 0, 0)
      self.atm.append((float(center[2]), ienv, 0, 0, 0, 0))
      ienv += 3

    # build bas
    for ii, ao in enumerate(qc.ao_spec):

      # TODO: column-wise coefficients for contractions of same exponents
      l = lquant[ao['type']]
      if l > 6:
        raise ValueError('maximum angular moment supported by libcint is l=6')
      exps, coeffs = ao['coeffs'].T

      # scale coefficients
      coeffs = rescale_coeffs_libcint(exps, coeffs, l)

      # add exponents and coefficients to env
      self.env.extend(exps.tolist())
      self.env.extend(coeffs)

      # tuple for each contraction:
      #  (center, l, #prim, #contr, kappa, prim offset, contr offset, 0)
      self.bas.append((ao['atom'], l, ao['pnum'], 1, 1, ienv, ienv+ao['pnum'], 0))
      ienv += 2*ao['pnum']

    # store as arrays
    self.env = numpy.array(self.env)
    self.atm = numpy.array(self.atm, dtype=numpy.int32)
    self.bas = numpy.array(self.bas, dtype=numpy.int32)

    # store ctypes
    self.c_env = self.env.ctypes.data_as(ctypes.c_void_p)
    self.c_atm = self.atm.ctypes.data_as(ctypes.c_void_p)
    self.c_bas = self.bas.ctypes.data_as(ctypes.c_void_p)

    # some further parameters
    self._cartesian = cartesian
    self.Nshells = self.bas.shape[0]
    self.Norb = self.count_contractions()
    self.natm = ctypes.c_int(self.atm.shape[0])
    self.nbas = ctypes.c_int(self.bas.shape[0])
    self.shell_dims = numpy.array(self._get_dims(*range(self.Nshells)))
    self.shell_offsets = numpy.zeros((self.Nshells,), dtype=numpy.int)
    self.shell_offsets[1:] = numpy.cumsum(self.shell_dims[:-1])

  @property
  def cartesian(self):
    return self._cartesian

  @cartesian.setter
  def cartesian(self, value):
    self._cartesian = value
    self.Norb = self.count_contractions()
    self.shell_dims = numpy.array(self._get_dims(*range(self.Nshells)))
    self.shell_offsets = numpy.zeros((self.Nshells,), dtype=numpy.int)
    self.shell_offsets[1:] = numpy.cumsum(self.shell_dims[:-1])
    self.clear_all_blocks()

  def count_contractions(self):
    '''Counts the number of contracted gaussians.'''
    if self.cartesian:
      fun = libcint.CINTcgto_cart
    else:
      fun = libcint.CINTcgto_spheric
    ncntr = 0
    for i in range(self.Nshells):
      ncntr += fun(i, self.c_bas)
    return ncntr

  def count_primitives(self):
    '''Counts the number of primitive gaussians.'''
    if self.cartesian:
      fun = libcint.CINTcgto_cart
    else:
      fun = libcint.CINTcgto_spheric
    nprim = 0
    for i in range(self.Nshells):
      nprim += self.bas[i,2]*self.bas[i,3]*fun(i, self.c_bas)
    return nprim

  ########################
  ###  some Shortcuts  ###
  ########################

  def overlap(self, **kwargs):
    '''Shortcut to calculate overlap integrals <i|j>.'''
    return self.int1e('ovlp', **kwargs)

  def kinetic(self, **kwargs):
    '''Shortcut to calculate kinetic energy integrals <i|-0.5*\nabla^2|j>.'''
    return self.int1e('kin', **kwargs)

  def Vne(self, **kwargs):
    '''Shortcut to calculate electron-nuclear repulsion integrals \sum_a <i|-Z_a/r|j>.'''
    return self.int1e('nuc', **kwargs)

  def Hcore(self, **kwargs):
    '''Shortcut to calculate one-electron Hamiltonian H_core = T_e + V_ne.'''
    T = self.int1e('kin', **kwargs)
    V = self.int1e('nuc', **kwargs)
    if isinstance(T, list):
      return [t+v for t,v in zip(T, V)]
    else:
      return T + V

  def Vee(self, **kwargs):
    '''Shortcut to calculate electron-electron repulsion integrals.'''
    return self.int2e('', **kwargs)

  ########################################
  ###  set blocks of AO and MO ranges  ###
  ########################################

  def add_AO_block_1e(self, AOrange=None, AOrangei=None, AOrangej=None, shell=False):
    '''Specify block of 1-electron AO integrals to be calculated.

      **Parameters:**

        AOrangei, AOrangej : lists or range objects of integers
          Indices i and j specifing the desired block of AO integrals. If omitted, the whole range is take for the corrresponding index.
        AOrange : lists or range objects of integers
          Set same range for i and j.
        shell : boolean
          If True, indices i and j specify AO shells, instead of AO basis functions.
    '''
    if AOrange is not None:
      AOrangei = AOrange
      AOrangej = AOrange
    if AOrangei is None and AOrangej is None:
      return
    if shell:
      AOrangei = [self._shell2ao(s) for s in AOrangei] if AOrangei is not None else range(self.Norb)
      AOrangej = [self._shell2ao(s) for s in AOrangej] if AOrangej is not None else range(self.Norb)
    else:
      AOrangei = AOrangei if AOrangei is not None else range(self.Norb)
      AOrangej = AOrangej if AOrangej is not None else range(self.Norb)
    self.AO_blocks_1e.append((AOrangei, AOrangej))

  def add_AO_block_2e(self, AOrange=None, AOrangei=None, AOrangej=None, AOrangek=None, AOrangel=None, shell=False):
    '''Specify block of 2-electron AO integrals to be calculated.

      **Parameters:**

        AOrangei, AOrangej, AOrangek, AOrangel : lists or range objects of integers
          Indices i, j, k and l specifing the desired block of AO integrals. If omitted, the whole range is take for the corrresponding index.
        AOrange : lists or range objects of integers
          Set same range for i, j, k and l.
        shell : boolean
          If True, indices i, j, k and l specify AO shells, instead of AO basis functions.
    '''
    if AOrange is not None:
      AOrangei = AOrange
      AOrangej = AOrange
      AOrangek = AOrange
      AOrangel = AOrange

    if AOrangei is None and AOrangej is None and AOrangek is None and AOrangel is None:
      return
    if shell:
      AOrangei = [self._shell2ao(s) for s in AOrangei] if AOrangei is not None else range(self.Norb)
      AOrangej = [self._shell2ao(s) for s in AOrangej] if AOrangej is not None else range(self.Norb)
      AOrangek = [self._shell2ao(s) for s in AOrangek] if AOrangek is not None else range(self.Norb)
      AOrangel = [self._shell2ao(s) for s in AOrangel] if AOrangel is not None else range(self.Norb)
    else:
      AOrangei = AOrangei if AOrangei is not None else range(self.Norb)
      AOrangej = AOrangej if AOrangej is not None else range(self.Norb)
      AOrangek = AOrangek if AOrangek is not None else range(self.Norb)
      AOrangel = AOrangel if AOrangel is not None else range(self.Norb)
    self.AO_blocks_2e.append((AOrangei, AOrangej, AOrangek, AOrangel))

  def add_MO_block_1e(self, MOrangei=None, MOrangej=None):
    '''Specify block of 2-electron MO integrals to be calculated.

      **Parameters:**

        MOrangei, MOrangei : lists or range objects of integers
          Indices i and j specifing the desired block of MO integrals. If omitted, the whole range is take for the corrresponding index.
    '''
    if MOrangei is None and MOrangej is None:
      return
    MOrangei = MOrangei if MOrangei is not None else range(self.Norb)
    MOrangej = MOrangej if MOrangej is not None else range(self.Norb)
    self.MO_blocks_1e.append((MOrangei, MOrangej))

  def add_MO_block_2e(self, MOrangei=None, MOrangej=None, MOrangek=None, MOrangel=None):
    '''Specify block of 2-electron MO integrals to be calculated.

      **Parameters:**

        MOrangei, MOrangej : lists or range objects of integers
          Indices i, j, k and l specifing the desired block of MO integrals. If omitted, the whole range is take for the corrresponding index.
    '''
    if MOrangei is None and MOrangej is None and MOrangek is None and MOrangel is None:
      return
    MOrangei = MOrangei if MOrangei is not None else range(self.Norb)
    MOrangej = MOrangej if MOrangej is not None else range(self.Norb)
    MOrangek = MOrangek if MOrangek is not None else range(self.Norb)
    MOrangel = MOrangel if MOrangel is not None else range(self.Norb)
    self.MO_blocks_2e.append((MOrangei, MOrangej, MOrangek, MOrangel))

  def clear_all_blocks(self):
    self.clear_AO_blocks_1e()
    self.clear_AO_blocks_2e()
    self.clear_MO_blocks_1e()
    self.clear_MO_blocks_2e()

  def clear_AO_blocks_1e(self):
    self.AO_blocks_1e = []

  def clear_AO_blocks_2e(self):
    self.AO_blocks_2e = []

  def clear_MO_blocks_1e(self):
    self.MO_blocks_1e = []

  def clear_MO_blocks_2e(self):
    self.MO_blocks_2e = []

  #######################################
  ###  calculate requested integrals  ###
  #######################################

  def int1e(self, operator, asMO=False, max_dims=0):
    '''Calculates all one-electron integrals <i|operator|j>.

      Use add_AO_block_1e() or add_MO_block_1e() if only a subset of AO or MO integrals is needed.

      **Parameters:**

        operator : string
          Base name of function/integral in libcint.
        asMO : boolean
          if True, transform from AO to MO basis.
        max_dims : int
          If > 0, calculate AO Integrals in slices containing no more than max_dims AOs (but at least one shell). Slower, but requires less memory.

      **Returns:**

        (Hermitian) 2D array of integrals if all or only one block of AOs/MOs is requested, otherwise a list of 2D arrays in the order of blocks specified.
    '''

    # get list of all required AO shells
    shells = None
    if asMO:
      if not self.MO_blocks_1e:
        AOs = None
        shells = [range(self.Nshells), range(self.Nshells)]
        shell_offsets = self.shell_offsets
      else:
        MOs = [sorted(set(chain(*blocks))) for blocks in zip(*self.MO_blocks_1e)]
        MOs = sorted(set(chain(*MOs)))
        AOs = numpy.where(~numpy.isclose(self.qc.mo_spec.coeffs[MOs,:], 1e-15).all(axis=0))[0]
        AOs = [AOs, AOs]  # same AO shells for i and j
    else:
      if not self.AO_blocks_1e:
        AOs = None
        shells = [range(self.Nshells), range(self.Nshells)]
        shell_offsets = self.shell_offsets
      else:
        AOs = [numpy.array(sorted(set(chain(*blocks)))) for blocks in zip(*self.AO_blocks_1e)]
    if not shells:
      shells = [sorted(set([self._ao2shell(ao) for ao in aos])) for aos in AOs]

    # set up libcint function
    if self.cartesian:
      ext = 'cart'
    else:
      ext = 'sph'
    operator = 'cint1e_%s_%s' %(operator, ext)
    operator = getattr(libcint, operator)
    operator.restype = ctypes.c_void_p

    # loop over shells
    if not asMO or not max_dims:
      # single slice
      if AOs is None:
        mat = numpy.zeros((self.Norb, self.Norb))
      else:
        mat = numpy.zeros((len(AOs[0]), len(AOs[1])))
      for i in shells[0]:
        for j in shells[1]:

          # hermitian
          if j < i:
            continue

          res = self._libcint1e(operator, i, j)
          ii, jj = self.shell_offsets[[i, j]]
          di, dj = res.shape

          if AOs:
            # only keep required AOs
            AOi = numpy.array(sorted(set(range(ii, ii+di)).intersection(AOs[0]))) - ii
            AOj = numpy.array(sorted(set(range(jj, jj+dj)).intersection(AOs[1]))) - jj
            res = res[numpy.ix_(AOi, AOj)]
            di, dj = res.shape
            ii = numpy.where(AOs[0]>=ii)[0][0]
            jj = numpy.where(AOs[1]>=jj)[0][0]

          mat[ii:ii+di,jj:jj+dj] = res
          mat[jj:jj+dj,ii:ii+di] = res.T

      if not asMO:
        if not self.AO_blocks_1e:
          return mat
        else:
          results = []
          for block in self.AO_blocks_1e:
            # indices need to be updated to account for discarded AOs
            maski, maskj = [numpy.array([False]*self.Norb)]*2
            maski[block[0]] = True
            maskj[block[1]] = True
            results.append(mat[numpy.ix_(maski[AOs[0]], maskj[AOs[1]])])
      else:
        results = []
        if AOs is None:
          AOrange = range(self.Norb)
        else:
          AOrange = AOs[0]
        if not self.MO_blocks_1e:
          return ao2mo(mat, self.qc.mo_spec.coeffs[:,AOrange])
        for block in self.MO_blocks_1e:
          results.append(ao2mo(
            mat, self.qc.mo_spec.coeffs[:,AOrange],
            MOrangei=block[0], MOrangej=block[1],
          ))

    else:
      # slice into smaller pieces to save some memory (asMO=True)
      if not self.MO_blocks_1e:
        results = [numpy.zeros((self.Norb, self.Norb))]
      else:
        results = [numpy.zeros((len(block[0]), len(block[1]))) for block in self.MO_blocks_1e]
      slices = self._get_slices(max_dims, shells[0])
      for s in slices:
        if AOs is None:
          Norb = numpy.sum(self.shell_dims[s[0]:s[-1]+1])
          mat_s = numpy.zeros((Norb, self.Norb))
        else:
          AOs_slice = set(chain(*[self._shell2ao(i) for i in s]))
          AOs_slice = sorted(AOs_slice.intersection(AOs[0]))
          mat_s = numpy.zeros((len(AOs_slice), len(AOs[1])))
        ss = self.shell_offsets[s[0]]
        if AOs is not None:
          ss = numpy.where(AOs[0]>=ss)[0][0]
        for i in s:
          for j in shells[1]:

            # hermitian
            if j < i and j in s:
              continue

            res = self._libcint1e(operator, i, j)
            ii, jj = self.shell_offsets[[i, j]]
            di, dj = res.shape

            if AOs:
              # only keep required AOs
              AOi = numpy.array(sorted(set(range(ii, ii+di)).intersection(AOs[0]))) - ii
              AOj = numpy.array(sorted(set(range(jj, jj+dj)).intersection(AOs[1]))) - jj
              res = res[numpy.ix_(AOi, AOj)]
              di, dj = res.shape
              ii = numpy.where(AOs[0]>=ii)[0][0]
              jj = numpy.where(AOs[1]>=jj)[0][0]

            ii_s = ii - ss
            mat_s[ii_s:ii_s+di,jj:jj+dj] = res
            if j in s and j != i:
              jj_s = jj - ss
              mat_s[jj_s:jj_s+dj,ii:ii+di] = res.T

        if AOs is None:
          AOrange = range(self.Norb)
          AOslice = range(self.shell_offsets[s[0]], self.shell_offsets[s[-1]]+self.shell_dims[s[-1]])
        else:
          AOrange = AOs[0]
          AOslice = [numpy.where(AOs[0]>=ao)[0][0] for ao in AOs_slice]
        if not self.MO_blocks_1e:
          results[0] += _ao2mo_slice(mat_s, self.qc.mo_spec.coeffs[:,AOrange], AOslice=AOslice)
        for i, block in enumerate(self.MO_blocks_1e):
          results[i] += _ao2mo_slice(
            mat_s, self.qc.mo_spec.coeffs[:,AOrange], AOslice=AOslice,
            MOrangei=block[0], MOrangej=block[1],
          )

    # return results
    if len(results) == 1:
      return results[0]
    return results

  def int2e(self, operator, asMO=False, max_dims=0):
    '''Calculates all two-electron integrals <ij|operator|kl>.

      Use add_AO_block_2e() or add_MO_block_2e() if only a subset of AO or MO integrals is needed.

      **Parameters:**

        operator : string
          Base name of function/integral in libcint.
        asMO : boolean
          if True, transform from AO to MO basis.
        MOrangei, MOrangej, MOrangek, MOrangel : list|range object|None
          Only transform selected MOs for indices i, j, k and l respectively.
        MOrange : list|range object|None
          sets MOrangei to MOrangel at the same range.
        max_dims : int
          If > 0, calculate AO Integrals in slices containing no more than max_dims AOs (but at least one shell). Slower, but requires less memory.

      **Returns:**

        4D array of integrals.
    '''

    # get list of all required AO shells
    shells = None
    if asMO:
      if not self.MO_blocks_2e:
        AOs = None
        shells = [range(self.Nshells), range(self.Nshells), range(self.Nshells), range(self.Nshells)]
      else:
        MOs = [list(set(chain(*blocks))) for blocks in zip(*self.MO_blocks_2e)]
        MOs = list(set(chain(*MOs)))
        AOs = numpy.where(~numpy.isclose(self.qc.mo_spec.coeffs[MOs,:], 1e-15).all(axis=0))[0]
        AOs = [AOs, AOs, AOs, AOs]  # same AO shells for i, j, k and l
    else:
      if not self.AO_blocks_2e:
        AOs = None
        shells = [range(self.Nshells), range(self.Nshells), range(self.Nshells), range(self.Nshells)]
      else:
        AOs = [numpy.array(list(set(chain(*blocks)))) for blocks in zip(*self.AO_blocks_2e)]
    if not shells:
      shells = [list(set([self._ao2shell(ao) for ao in aos])) for aos in AOs]

    # set up libcint function and optimizer
    if self.cartesian:
      ext = 'cart'
    else:
      ext = 'sph'
    if operator:
      operator = 'cint2e_%s_%s' %(operator, ext)
    else:
      operator = 'cint2e_%s' %ext
    optimizer = '_'.join((operator, 'optimizer'))

    operator = getattr(libcint, operator)
    operator.restype = ctypes.c_void_p

    # initialize optimizer
    opt = CINTOpt()
    fopt = getattr(libcint, optimizer)
    fopt(ctypes.byref(opt), self.c_atm, self.natm, self.c_bas, self.nbas, self.c_env)
    #opt = ctypes.POINTER(ctypes.c_void_p)()  # disables optimizer

    # loop over shells
    if not asMO or not max_dims:
      if AOs is None:
        mat = numpy.zeros((self.Norb, self.Norb, self.Norb, self.Norb))
      else:
        mat = numpy.zeros((len(AOs[0]), len(AOs[1]), len(AOs[2]), len(AOs[3])))
      for i in shells[0]:
        for j in shells[1]:
          for k in shells[2]:
            for l in shells[3]:

              # hermitian
              if j < i or l < k:
                continue

              # exchange of electronic coordinates
              if (i < k) or (i == k and j < l):
                continue

              res = self._libcint2e(operator, i, j, k, l, opt)
              di, dj, dk, dl = res.shape
              ii, jj, kk, ll = self.shell_offsets[[i, j, k, l]]

              if AOs:
                # only keep required AOs
                AOi = numpy.array(sorted(set(range(ii, ii+di)).intersection(AOs[0]))) - ii
                AOj = numpy.array(sorted(set(range(jj, jj+dj)).intersection(AOs[1]))) - jj
                AOk = numpy.array(sorted(set(range(kk, kk+dk)).intersection(AOs[2]))) - kk
                AOl = numpy.array(sorted(set(range(ll, ll+dl)).intersection(AOs[3]))) - ll
                res = res[numpy.ix_(AOi, AOj, AOk, AOl)]
                di, dj, dk, dl = res.shape
                ii = numpy.where(AOs[0]>=ii)[0][0]
                jj = numpy.where(AOs[1]>=jj)[0][0]
                kk = numpy.where(AOs[2]>=kk)[0][0]
                ll = numpy.where(AOs[3]>=ll)[0][0]

              # store/replicate results
              mat[ii:ii+di,jj:jj+dj,kk:kk+dk,ll:ll+dl] = res

              mat[jj:jj+dj,ii:ii+di,kk:kk+dk,ll:ll+dl] = numpy.swapaxes(res, 0, 1)
              mat[ii:ii+di,jj:jj+dj,ll:ll+dl,kk:kk+dk] = numpy.swapaxes(res, 2, 3)
              mat[jj:jj+dj,ii:ii+di,ll:ll+dl,kk:kk+dk] = moveaxis(res, (0,1,2,3), (1,0,3,2))

              mat[kk:kk+dk,ll:ll+dl,ii:ii+di,jj:jj+dj] = moveaxis(res, (0,1,2,3), (2,3,0,1))
              mat[ll:ll+dl,kk:kk+dk,ii:ii+di,jj:jj+dj] = moveaxis(res, (0,1,2,3), (2,3,1,0))
              mat[kk:kk+dk,ll:ll+dl,jj:jj+dj,ii:ii+di] = moveaxis(res, (0,1,2,3), (3,2,0,1))
              mat[ll:ll+dl,kk:kk+dk,jj:jj+dj,ii:ii+di] = moveaxis(res, (0,1,2,3), (3,2,1,0))

      if not asMO:
        if not self.AO_blocks_2e:
          results = [mat]
        else:
          results = []
          for block in self.AO_blocks_2e:
            # indices need to be updated to account for discarded AOs
            maski, maskj, maskk, maskl = [numpy.array([False]*self.Norb)]*4
            maski[block[0]] = True
            maskj[block[1]] = True
            maskk[block[2]] = True
            maskl[block[3]] = True
            results.append(
              mat[numpy.ix_(maski[AOs[0]], maskj[AOs[1]], maskk[AOs[2]], maskl[AOs[3]])]
            )
      else:
        results = []
        if AOs is None:
          AOrange = range(self.Norb)
        else:
          AOrange = AOs[0]
        if not self.MO_blocks_2e:
          results = [ao2mo(mat, self.qc.mo_spec.coeffs[:,AOrange])]
        for block in self.MO_blocks_2e:
          results.append(ao2mo(
            mat, self.qc.mo_spec.coeffs[:,AOrange],
            MOrangei=block[0], MOrangej=block[1], MOrangek=block[2], MOrangel=block[3],
          ))

    else:
      # slice into smaller pieces to save some memory (asMO=True)
      if not self.MO_blocks_2e:
        results = [numpy.zeros((self.Norb, self.Norb, self.Norb, self.Norb))]
      else:
        results = [numpy.zeros((len(block[0]), len(block[1]), len(block[2]), len(block[3]))) for block in self.MO_blocks_2e]
      slices = self._get_slices(max_dims, shells[0])
      for s in slices:
        if AOs is None:
          Norb = numpy.sum(self.shell_dims[s[0]:s[-1]+1])
          mat_s = numpy.zeros((Norb, self.Norb, self.Norb, self.Norb))
        else:
          AOs_slice = set(chain(*[self._shell2ao(i) for i in s]))
          AOs_slice = sorted(AOs_slice.intersection(AOs[0]))
          mat_s = numpy.zeros((len(AOs_slice), len(AOs[1]), len(AOs[2]), len(AOs[3])))
        ss = self.shell_offsets[s[0]]
        if AOs is not None:
          ss = numpy.where(AOs[0]>=ss)[0][0]
        for i in s:
          for j in shells[1]:
            for k in shells[2]:
              for l in shells[3]:

                ## hermitian
                if j in s and j < i:
                  continue
                if l < k:
                  continue

                # exchange of electronic coordinates
                if k in s and k < i:
                  continue
                if l in s and i == k and l < j:
                  continue

                res = self._libcint2e(operator, i, j, k, l, opt)
                ii, jj, kk, ll = self.shell_offsets[[i, j, k, l]]
                di, dj, dk, dl = res.shape

                if AOs:
                  # only keep required AOs
                  AOi = numpy.array(sorted(set(range(ii, ii+di)).intersection(AOs[0]))) - ii
                  AOj = numpy.array(sorted(set(range(jj, jj+dj)).intersection(AOs[1]))) - jj
                  AOk = numpy.array(sorted(set(range(kk, kk+dk)).intersection(AOs[2]))) - kk
                  AOl = numpy.array(sorted(set(range(ll, ll+dl)).intersection(AOs[3]))) - ll
                  res = res[numpy.ix_(AOi, AOj, AOk, AOl)]
                  di, dj, dk, dl = res.shape
                  ii = numpy.where(AOs[0]>=ii)[0][0]
                  jj = numpy.where(AOs[1]>=jj)[0][0]
                  kk = numpy.where(AOs[2]>=kk)[0][0]
                  ll = numpy.where(AOs[3]>=ll)[0][0]

                # store/replicate results
                ii_s = ii - ss
                mat_s[ii_s:ii_s+di,jj:jj+dj,kk:kk+dk,ll:ll+dl] = res
                mat_s[ii_s:ii_s+di,jj:jj+dj,ll:ll+dl,kk:kk+dk] = numpy.swapaxes(res, 2, 3)

                if j in s:
                  jj_s = jj - ss
                  mat_s[jj_s:jj_s+dj,ii:ii+di,kk:kk+dk,ll:ll+dl] = numpy.swapaxes(res, 0, 1)
                  mat_s[jj_s:jj_s+dj,ii:ii+di,ll:ll+dl,kk:kk+dk] = moveaxis(res, (0,1,2,3), (1,0,3,2))

                if k in s:
                  kk_s = kk - ss
                  mat_s[kk_s:kk_s+dk,ll:ll+dl,ii:ii+di,jj:jj+dj] = moveaxis(res, (0,1,2,3), (2,3,0,1))
                  mat_s[kk_s:kk_s+dk,ll:ll+dl,jj:jj+dj,ii:ii+di] = moveaxis(res, (0,1,2,3), (3,2,0,1))

                if l in s:
                  ll_s = ll - ss
                  mat_s[ll_s:ll_s+dl,kk:kk+dk,ii:ii+di,jj:jj+dj] = moveaxis(res, (0,1,2,3), (2,3,1,0))
                  mat_s[ll_s:ll_s+dl,kk:kk+dk,jj:jj+dj,ii:ii+di] = moveaxis(res, (0,1,2,3), (3,2,1,0))

        if AOs is None:
          AOrange = range(self.Norb)
          AOslice = range(self.shell_offsets[s[0]], self.shell_offsets[s[-1]]+self.shell_dims[s[-1]])
        else:
          AOrange = AOs[0]
          AOslice = [numpy.where(AOs[0]>=ao)[0][0] for ao in AOs_slice]
        if not self.MO_blocks_2e:
          results[0] += _ao2mo_slice(mat_s, self.qc.mo_spec.coeffs[:,AOrange], AOslice=AOslice)
        for i, block in enumerate(self.MO_blocks_2e):
          results[i] += _ao2mo_slice(
            mat_s, self.qc.mo_spec.coeffs[:,AOrange], AOslice=AOslice,
            MOrangei=block[0], MOrangej=block[1], MOrangek=block[2], MOrangel=block[3],
          )

    # release optimizer
    libcint.CINTdel_optimizer(ctypes.byref(opt))

    # return results
    if len(results) == 1:
      return results[0]
    return results

  ############################
  ###  internal functions  ###
  ############################

  def _shell2ao(self, i):
    '''Returns all AO indices for given shell.'''
    if i == self.Nshells - 1:
      return range(self.shell_offsets[i], self.Norb)
    return range(self.shell_offsets[i], self.shell_offsets[i+1])

  def _ao2shell(self, i):
    '''Returns shell index for given AO.'''
    return numpy.where(self.shell_offsets <= i)[0][-1]

  def _get_slices(self, max_dims, shells=None):
    '''Slices shells in separate ranges.'''
    if not shells:
      shells = range(self.Nshells)
    shells = list(shells)
    slices = []
    _slice = shells[:1]
    for s in shells[1:]:
      if sum(self.shell_dims[_slice+[s]]) > max_dims:
        slices.append(_slice)
        _slice = [s]
      else:
        _slice.append(s)
    slices.append(_slice)
    return slices

  def _get_dims(self, *indices):
    '''Returns number of basis functions for shell indices'''
    if self.cartesian:
      fun = libcint.CINTcgto_cart
    else:
      fun = libcint.CINTcgto_spheric
    dims = []
    for ind in indices:
      dims.append( fun(ctypes.c_int(ind), self.c_bas) )
    if len(indices) == 1:
      return dims[0]
    return dims

  def _norm_cart_shell(self, i):
    '''Calculates normalization factors for shell i to rescale cartesian integrals.'''

    # check if cached
    S = self.cache_norm.get(i, None)

    if S is None:
      # calculate square root of self overlap
      di = self._get_dims(i)
      mat = (ctypes.c_double * di*di)()
      shls = (ctypes.c_int * 2)(i, i)
      libcint.cint1e_ovlp_cart.restype = ctypes.c_void_p
      libcint.cint1e_ovlp_cart(mat, shls, self.c_atm, self.natm, self.c_bas, self.nbas, self.c_env)
      S = numpy.reshape(mat, (di, di)).T
      S = numpy.sqrt(numpy.diag(S))

      # add to cache
      self.cache_norm[i] = S

    return S

  def _libcint1e(self, operator, i, j):
    '''Calls libcint to evaluate 1-electron integrals over shells i and j.

    **Parameters:**

      operator : libcint function to call.
      i, j : Shells for bra and ket vectors.

    **Returns:**

      2D array of integrals.

    '''

    di, dj = self.shell_dims[[i, j]]
    mat = (ctypes.c_double * di*dj)()
    shls = (ctypes.c_int * 2)(i, j)

    operator(mat, shls, self.c_atm, self.natm, self.c_bas, self.nbas, self.c_env)
    mat = numpy.reshape(mat, (dj, di)).T

    # cartesian integrals need to be rescaled according to overlap matrix
    if self.cartesian:
      mat /= self._norm_cart_shell(i).reshape(di,1)
      mat /= self._norm_cart_shell(j).reshape(1,dj)

    # switch order of basis functions (angl>1)
    angl = self.bas[i][1]
    if angl > 1:
      order = get_order(angl, self.cartesian)
      mat = mat[order,:]
    angl = self.bas[j][1]
    if angl > 1:
      order = get_order(angl, self.cartesian)
      mat = mat[:,order]

    return mat

  def _libcint2e(self, operator, i, j, k, l, opt):
    '''Calls libcint to evaluate 2-electron integrals over shells i, j, k and l.

    **Parameters:**

      operator : libcint function to call.
      i, j, k, l : Shells for bra and ket vectors.
      opt : libcint optimizer to be passed to operator.

    **Returns:**

      4D array of integrals.

    '''

    di, dj, dk, dl = self.shell_dims[[i, j, k, l]]
    mat = (ctypes.c_double * di*dj*dk*dl)()
    shls = (ctypes.c_int * 4)(i, j, k, l)

    operator(mat, shls, self.c_atm, self.natm, self.c_bas, self.nbas, self.c_env, opt)
    mat = numpy.reshape(mat, (dl, dk, dj, di))
    mat = moveaxis(mat, (0, 1, 2, 3), (3, 2, 1, 0))

    # cartesian integrals need to be rescaled according to overlap matrix
    if self.cartesian:
      mat /= self._norm_cart_shell(i).reshape(di,1,1,1)
      mat /= self._norm_cart_shell(j).reshape(1,dj,1,1)
      mat /= self._norm_cart_shell(k).reshape(1,1,dk,1)
      mat /= self._norm_cart_shell(l).reshape(1,1,1,dl)

    # switch order of basis functions (angl>1)
    angl = self.bas[i][1]
    if angl > 1:
      order = get_order(angl, self.cartesian)
      mat = mat[order,:,:,:]
    angl = self.bas[j][1]
    if angl > 1:
      order = get_order(angl, self.cartesian)
      mat = mat[:,order,:,:]
    angl = self.bas[k][1]
    if angl > 1:
      order = get_order(angl, self.cartesian)
      mat = mat[:,:,order,:]
    angl = self.bas[l][1]
    if angl > 1:
      order = get_order(angl, self.cartesian)
      mat = mat[:,:,:,order]

    return mat

#######################################################
###  rescaling and transformations of coefficients  ###
#######################################################

def ao2mo(mat, coeffs, MOrange=None, MOrangei=None, MOrangej=None, MOrangek=None, MOrangel=None):
  '''Transforms array of one- or two-electron integrals from AO to MO basis.

  **Parameters:**

    mat : 2 or 4 dimensional numpy.ndarray
      integrals to be transformed
    coeffs : 2 dimensional numpy.ndarry
      MO coefficients
    MOrangei, MOrangej, MOrangek, MOrangel : list|range object|None
      Only transform selected MOs for indices i, j, k and l respectively.
    MOrange: list|range object|None
      Set same range for all indices.

  **Returns:**

    numpy.ndarray

  '''
  assert len(mat.shape) in (2, 4), "'mat' musst be of size 2 or 4."

  if MOrange is not None:
    MOrangei = MOrange
    MOrangej = MOrange
    MOrangek = MOrange
    MOrangel = MOrange
  if MOrangei is None:
    MOrangei = range(coeffs.shape[0])
  if MOrangej is None:
    MOrangej = range(coeffs.shape[0])
  if MOrangek is None:
    MOrangek = range(coeffs.shape[0])
  if MOrangel is None:
    MOrangel = range(coeffs.shape[0])

  # discard zero columns in MO coeffs
  if len(mat.shape) == 2:
    MOs = list(set(chain(MOrangei, MOrangej)))
  elif len(mat.shape) == 4:
    MOs = list(set(chain(MOrangei, MOrangej, MOrangek, MOrangel)))
  AOs = numpy.where(~numpy.isclose(coeffs[MOs,:], 1e-15).all(axis=0))[0]
  coeffs = coeffs[:,AOs]

  if len(mat.shape) == 2:
    # 1-electron integrals
    mat = mat[numpy.ix_(AOs, AOs)]
    return numpy.dot(coeffs[MOrangei,:], numpy.dot(mat, coeffs[MOrangej,:].T))
  elif len(mat.shape) == 4:
    # 2-electron integrals
    mat = mat[numpy.ix_(AOs, AOs, AOs, AOs)]
    mat = numpy.tensordot(mat, coeffs[MOrangei,:], axes=(0, 1))
    mat = numpy.tensordot(mat, coeffs[MOrangej,:], axes=(0, 1))
    mat = numpy.tensordot(mat, coeffs[MOrangek,:], axes=(0, 1))
    mat = numpy.tensordot(mat, coeffs[MOrangel,:], axes=(0, 1))
    return mat

def _ao2mo_slice(mat, coeffs, AOslice, MOrange=None, MOrangei=None, MOrangej=None, MOrangek=None, MOrangel=None):
  '''Transforms array of one- or two-electron integrals from AO to MO basis. AO sliced version.'''
  assert len(mat.shape) in (2, 4), "'mat' musst be of size 2 or 4."

  if MOrange is not None:
    MOrangei = MOrange
    MOrangej = MOrange
    MOrangek = MOrange
    MOrangel = MOrange
  if MOrangei is None:
    MOrangei = range(coeffs.shape[0])
  if MOrangej is None:
    MOrangej = range(coeffs.shape[0])
  if MOrangek is None:
    MOrangek = range(coeffs.shape[0])
  if MOrangel is None:
    MOrangel = range(coeffs.shape[0])

  # discard zero columns in MO coeffs
  if len(mat.shape) == 2:
    MOs = list(set(chain(MOrangei, MOrangej)))
  elif len(mat.shape) == 4:
    MOs = list(set(chain(MOrangei, MOrangej, MOrangek, MOrangel)))
  AOs = numpy.where(~numpy.isclose(coeffs[MOs,:], 1e-15).all(axis=0))[0]
  AOs_ = [ao for ao in AOs if ao in AOslice]
  AOs_mat = [ao-AOslice[0] for ao in AOs if ao in AOslice]

  if not AOs_mat:
    # no relevant AOs in this slice
    if len(mat.shape) == 2:
      return numpy.zeros((len(MOrangei), len(MOrangej)))
    elif len(mat.shape) == 4:
      return numpy.zeros((len(MOrangei), len(MOrangej), len(MOrangek), len(MOrangel)))

  if len(mat.shape) == 2:
    # 1-electron integrals
    mat = mat[numpy.ix_(AOs_mat,AOs)]
    return numpy.dot(coeffs[numpy.ix_(MOrangei, AOs_)], numpy.dot(mat, coeffs[numpy.ix_(MOrangej, AOs)].T))
  elif len(mat.shape) == 4:
    # 2-electron integrals
    mat = mat[numpy.ix_(AOs_mat, AOs, AOs, AOs)]
    mat = numpy.tensordot(mat, coeffs[numpy.ix_(MOrangei, AOs_)], axes=(0, 1))
    coeffs = coeffs[:,AOs]
    mat = numpy.tensordot(mat, coeffs[MOrangej,:], axes=(0, 1))
    mat = numpy.tensordot(mat, coeffs[MOrangek,:], axes=(0, 1))
    mat = numpy.tensordot(mat, coeffs[MOrangel,:], axes=(0, 1))
    return mat

def rescale_coeffs_libcint(exps, coeffs, l):
  return [ c*libcint.CINTgto_norm(l, ctypes.c_double(e)) for e, c in zip(exps, coeffs) ]

def get_order(l, cartesian):
  '''Returns indices to transform libcint order to orbkit order.'''
  if l == 0:
    return (0,)
  if l == 1:
    return (0,1,2)
  order = {
    # l : (cartesian, spherical)
    2 : ((0,3,5,1,2,4), range(2)),
    3 : ((0,6,9,3,1,2,5,8,7,4), range(3)),
    4 : ((0,10,14,1,2,6,11,9,13,3,5,12,4,7,8), range(4)),
  }
  if cartesian:
    return order[l][0]
  else:
    return order[l][1]
