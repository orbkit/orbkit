import numpy
import ctypes

from ..tools import lquant
from ..display import display
from itertools import chain, combinations_with_replacement
import time

try:
  from numpy import moveaxis
except:
  from orbkit.tools import moveaxis

# search and load libcint from $PATH
import os
libcint = None
for dirname in os.environ['PATH'].split(':'):
  if os.path.isfile(os.path.join(dirname, 'lib', 'libcint.so')):
    libcint = ctypes.cdll.LoadLibrary(os.path.abspath(os.path.join(dirname, 'lib', 'libcint.so')))
    break
if libcint is None:
  raise ImportError("libcint not found: please add to PATH environment variable")

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

# set return type
libcint.CINTgto_norm.restype = ctypes.c_double

class struct:
  pass

def match_order(a, ref):
  '''Determines order such that a[order,:] == ref.'''
  # inverse argsort of ref
  if ref.ndim == 1:
    fwd = numpy.argsort(ref)
  else:
    fwd = numpy.lexsort(ref.T, axis=0)
  inv = numpy.empty_like(fwd)
  inv[fwd] = numpy.arange(fwd.size)
  # argsort of a
  if a.ndim == 1:
    arg = numpy.argsort(a)
  else:
    arg = numpy.lexsort(a.T, axis=0)
  # now order is given by arg[inv]
  return arg[inv]

############################
###  external interface  ###
############################

# TODO:
# - non-hermitian operators (?)
# - operators including gradients (return tensors)

class AOIntegrals():
  '''Interface to calculate AO Integrals with libcint.
  https://github.com/sunqm/libcint
  '''
  def __init__(self, qc):

    # general parameters
    self.qc = qc
    self.clear_all_blocks()

    # libcint parameters
    ienv = 0
    self.basis = struct()
    self.env = [0]*ienv
    self.atm = []
    self.bas = []
    self.basis.cart_norm = []  # renormalization factors for cartesian integrals

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

    # convert to arrays
    self.env = numpy.array(self.env)
    self.atm = numpy.array(self.atm, dtype=numpy.int32)
    self.bas = numpy.array(self.bas, dtype=numpy.int32)
    self.basis.bas = self.bas

    # store ctypes
    self.basis.c_env = self.env.ctypes.data_as(ctypes.c_void_p)
    self.basis.c_atm = self.atm.ctypes.data_as(ctypes.c_void_p)
    self.basis.c_bas = self.bas.ctypes.data_as(ctypes.c_void_p)

    # some further parameters
    self.basis.cartesian = not qc.ao_spec.spherical
    self.basis.Nshells = self.basis.bas.shape[0]
    self.basis.Norb = self.count_contractions()
    self.basis.natm = ctypes.c_int(self.atm.shape[0])
    self.basis.nbas = ctypes.c_int(self.bas.shape[0])
    self.basis.shell_dims = numpy.array(self._get_dims(*range(self.basis.Nshells)))
    self.basis.shell_offsets = _get_shell_offsets(self.basis.shell_dims)
    self.basis.order = [self._get_order(i) for i in range(self.Nshells)]
    if self.basis.cartesian:
      self._norm_cart_shell()

  @property
  def cartesian(self):
    return self.basis.cartesian

  @cartesian.setter
  def cartesian(self, value):
    self.basis.cartesian = value
    self.basis.Norb = self.count_contractions()
    self.basis.shell_dims = numpy.array(self._get_dims(*range(self.Nshells)))
    self.basis.shell_offsets = _get_shell_offsets(self.basis.shell_dims)
    self.basis.order = [self._get_order(i) for i in range(self.Nshells)]
    if value:
      self._norm_cart_shell()
    self.clear_all_blocks()

  @property
  def Norb(self):
    return self.basis.Norb

  @property
  def Nshells(self):
    return self.basis.Nshells

  def count_contractions(self):
    '''Counts the number of contracted gaussians.'''
    if self.cartesian:
      fun = libcint.CINTcgto_cart
    else:
      fun = libcint.CINTcgto_spheric
    ncntr = 0
    for i in range(self.Nshells):
      ncntr += fun(i, self.basis.c_bas)
    return ncntr

  def count_primitives(self):
    '''Counts the number of primitive gaussians.'''
    if self.cartesian:
      fun = libcint.CINTcgto_cart
    else:
      fun = libcint.CINTcgto_spheric
    nprim = 0
    for i in range(self.Nshells):
      nprim += self.bas[i,2]*self.bas[i,3]*fun(i, self.basis.c_bas)
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
    if all((AOrangei is None, AOrangej is None)):
      return
    if shell:
      AOrangei = list(chain(*[self._shell2ao(s) for s in AOrangei])) if AOrangei is not None else range(self.Norb)
      AOrangej = list(chain(*[self._shell2ao(s) for s in AOrangej])) if AOrangej is not None else range(self.Norb)
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
    if all((AOrangei is None, AOrangej is None, AOrangek is None, AOrangel is None)):
      return
    if shell:
      AOrangei = list(chain(*[self._shell2ao(s) for s in AOrangei])) if AOrangei is not None else range(self.Norb)
      AOrangej = list(chain(*[self._shell2ao(s) for s in AOrangej])) if AOrangej is not None else range(self.Norb)
      AOrangek = list(chain(*[self._shell2ao(s) for s in AOrangek])) if AOrangek is not None else range(self.Norb)
      AOrangel = list(chain(*[self._shell2ao(s) for s in AOrangel])) if AOrangel is not None else range(self.Norb)
    else:
      AOrangei = AOrangei if AOrangei is not None else range(self.Norb)
      AOrangej = AOrangej if AOrangej is not None else range(self.Norb)
      AOrangek = AOrangek if AOrangek is not None else range(self.Norb)
      AOrangel = AOrangel if AOrangel is not None else range(self.Norb)
    self.AO_blocks_2e.append((AOrangei, AOrangej, AOrangek, AOrangel))

  def add_MO_block_1e(self, MOrange=None, MOrangei=None, MOrangej=None):
    '''Specify block of 2-electron MO integrals to be calculated.

      **Parameters:**

      MOrangei, MOrangej : lists or range objects of integers
        Indices i and j specifing the desired block of MO integrals. If omitted, the whole range is take for the corrresponding index.
      MOrange : lists or range objects of integers
        Set same range for i and j.
    '''
    if MOrange is not None:
      MOrangei = MOrange
      MOrangej = MOrange
    if all((MOrangei is None, MOrangej is None)):
      return
    MOrangei = MOrangei if MOrangei is not None else range(self.Norb)
    MOrangej = MOrangej if MOrangej is not None else range(self.Norb)
    self.MO_blocks_1e.append((MOrangei, MOrangej))

  def add_MO_block_2e(self, MOrange=None, MOrangei=None, MOrangej=None, MOrangek=None, MOrangel=None):
    '''Specify block of 2-electron MO integrals to be calculated.

      **Parameters:**

      MOrangei, MOrangej : lists or range objects of integers
        Indices i, j, k and l specifing the desired block of MO integrals. If omitted, the whole range is take for the corrresponding index.
      MOrange : lists or range objects of integers
        Set same range for i, j, k and l.
    '''
    if MOrange is not None:
      MOrangei = MOrange
      MOrangej = MOrange
      MOrangek = MOrange
      MOrangel = MOrange
    if all((MOrangei is None, MOrangej is None, MOrangek is None, MOrangel is None)):
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
        AOs = [range(self.Norb)]*2
        shells = [range(self.Nshells), range(self.Nshells)]
      else:
        MOs = [sorted(set(chain(*blocks))) for blocks in zip(*self.MO_blocks_1e)]
        MOs = sorted(set(chain(*MOs)))
        AOs = numpy.where(~numpy.isclose(self.qc.mo_spec.coeffs[MOs,:], 1e-15).all(axis=0))[0]
        AOs = [AOs, AOs]  # same AO indicies for i and j
    else:
      if not self.AO_blocks_1e:
        AOs = [range(self.Norb)]*2
        shells = [range(self.Nshells), range(self.Nshells)]
      else:
        AOs = [numpy.array(sorted(set(chain(*blocks)))) for blocks in zip(*self.AO_blocks_1e)]
    if not shells:
      shells = [sorted(set([self._ao2shell(ao) for ao in aos])) for aos in AOs]


    if not asMO or not max_dims:

      # calculate whole AO matrix
      mat = _libcint1e(self.basis, operator, *shells)

      # remove orbitals not presents in AOs
      masks = []
      for i in range(2):
        # mask of all AOs
        mask = numpy.array([False]*self.Norb)
        # only enable request AOs
        mask[AOs[i]] = True
        # cut skipped shells
        mask = mask[list(chain(*[self._shell2ao(s) for s in shells[i]]))]
        masks.append(mask)
      mat = mat[numpy.ix_(*masks)]

      if not asMO:
        # get AO blocks
        if not self.AO_blocks_1e:
          results = [mat]
        else:
          results = []
          for block in self.AO_blocks_1e:
            # indices need to be updated to account for discarded AOs
            maski, maskj = [numpy.array([False]*self.Norb)]*2
            maski[block[0]] = True
            maskj[block[1]] = True
            results.append(mat[numpy.ix_(maski[AOs[0]], maskj[AOs[1]])])
      else:
        # convert to MO (blocks)
        results = []
        if not self.MO_blocks_1e:
          results.append(ao2mo(mat, self.qc.mo_spec.coeffs[:,AOs[0]]))
        for block in self.MO_blocks_1e:
          results.append(ao2mo(
            mat, self.qc.mo_spec.coeffs[:,AOs[0]],
            MOrangei=block[0], MOrangej=block[1],
          ))

    else:

      # slice into smaller pieces to save some memory (asMO=True)
      if not self.MO_blocks_1e:
        results = [numpy.zeros((self.Norb, self.Norb))]
      else:
        results = [numpy.zeros((len(block[0]), len(block[1]))) for block in self.MO_blocks_1e]

      # get slices by shells (avoid splitting AOs belonging to the same shell)
      slices = self._get_slices(max_dims, shells[0])

      for s in slices:
        shells_sliced = [s, shells[1]]

        # calculate AO matrix slice
        mat = _libcint1e(self.basis, operator, *shells_sliced)

        # remove orbitals not presents in AOs
        masks = []
        for i in range(2):
          # mask of all AOs
          mask = numpy.array([False]*self.Norb)
          # only enable request AOs
          mask[AOs[i]] = True
          # cut skipped shells
          mask = mask[list(chain(*[self._shell2ao(s_) for s_ in shells_sliced[i]]))]
          masks.append(mask)
        mat = mat[numpy.ix_(*masks)]

        # convert to MO (blocks)
        # all AOs in current slice
        AOslice = sorted(set(chain(*[self._shell2ao(i) for i in s])).intersection(AOs[0]))
        # match indices with respect to AOs[0]
        AOslice = [numpy.where(numpy.array(AOs[0])>=ao)[0][0] for ao in AOslice]

        if not self.MO_blocks_1e:
          results[0] += _ao2mo_slice(mat, self.qc.mo_spec.coeffs[:,AOs[0]], AOslice=AOslice)
        for i, block in enumerate(self.MO_blocks_1e):
          results[i] += _ao2mo_slice(
            mat, self.qc.mo_spec.coeffs[:,AOs[0]], AOslice=AOslice,
            MOrangei=block[0], MOrangej=block[1],
          )

    # return results
    if len(results) == 1:
      return results[0]
    return results

  def int2e(self, operator, asMO=False, max_dims=0, max_mem=0):
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
      max_mem : float
        Rough memory limit in MB, to determine max_dims automatically. Note: Shells with high angular momentum quantum number may exceed the limit, if choosen to small.

      **Returns:**

      4D array of integrals.
    '''

    # get list of all required AO shells
    shells = None
    if asMO:
      if not self.MO_blocks_2e:
        AOs = [range(self.Norb)]*4
        shells = [range(self.Nshells), range(self.Nshells), range(self.Nshells), range(self.Nshells)]
      else:
        MOs = [sorted(set(chain(*blocks))) for blocks in zip(*self.MO_blocks_2e)]
        MOs = sorted(set(chain(*MOs)))
        AOs = numpy.where(~numpy.isclose(self.qc.mo_spec.coeffs[MOs,:], 1e-15).all(axis=0))[0]
        AOs = [AOs, AOs, AOs, AOs]  # same AO indices for i, j, k and l
    else:
      if not self.AO_blocks_2e:
        AOs = [range(self.Norb)]*4
        shells = [range(self.Nshells), range(self.Nshells), range(self.Nshells), range(self.Nshells)]
      else:
        AOs = [numpy.array(sorted(set(chain(*blocks)))) for blocks in zip(*self.AO_blocks_2e)]
    if not shells:
      shells = [sorted(set([self._ao2shell(ao) for ao in aos])) for aos in AOs]

    if max_mem:
      # get memory requirements for AO matrix (j,k,l)
      shape = []
      for shell in shells:
        Norb = numpy.sum(self.basis.shell_dims[shell[0]:shell[-1]+1])
        shape.append(Norb)
      mem_jkl = numpy.prod(shape[1:])*8/1000**2
      # get memory requirements for MO matrix (i,j,k,l)
      if asMO:
        block = self.MO_blocks_2e[0]
        mem_mo = [numpy.prod([len(i) for i in block])*8/1000.**2 for block in self.MO_blocks_2e]
        max_mem -= sum(mem_mo)
      # calculate max_dims (estimate 2*mem_jkl to account for ao2mo)
      max_dims = int(max_mem/(2*mem_jkl))
      if max_dims > shape[0]:
        # only 1 slice required
        max_dims = 0

    if not asMO or not max_dims:

      # calculate whole AO matrix
      display('calculating 2-electron AO matrix for {:4d} shells along first index'.format(len(shells[0])))
      T0 = time.time()
      mat = _libcint2e(self.basis, operator, *shells)
      T1 = time.time()
      display('finished in {:8.2f} sec'.format(T1-T0))

      # remove orbitals not presents in AOs
      masks = []
      for i in range(4):
        # mask of all AOs
        mask = numpy.array([False]*self.Norb)
        # only enable request AOs
        mask[AOs[i]] = True
        # cut skipped shells
        mask = mask[list(chain(*[self._shell2ao(s) for s in shells[i]]))]
        masks.append(mask)
      mat = mat[numpy.ix_(*masks)]

      if not asMO:
        # get AO blocks
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
        # convert to MO (blocks)
        results = []
        if not self.MO_blocks_2e:
          results = [ao2mo(mat, self.qc.mo_spec.coeffs[:,AOs[0]])]
        display('ao2mo for MO block(s) {}'.format(len(self.MO_blocks_2e)))
        T0 = time.time()
        for ib, block in enumerate(self.MO_blocks_2e):
          results.append(ao2mo(
            mat, self.qc.mo_spec.coeffs[:,AOs[0]],
            MOrangei=block[0], MOrangej=block[1], MOrangek=block[2], MOrangel=block[3],
          ))
        T1 = time.time()
        display('\tfinished in {:8.2f} sec'.format(T1-T0))

    else:

      # slice into smaller pieces to save some memory (asMO=True)
      if not self.MO_blocks_2e:
        results = [numpy.zeros((self.Norb, self.Norb, self.Norb, self.Norb))]
      else:
        results = [numpy.zeros((len(block[0]), len(block[1]), len(block[2]), len(block[3]))) for block in self.MO_blocks_2e]

      # get slices by shells (avoid splitting AOs belonging to the same shell)
      slices = self._get_slices(max_dims, shells[0])

      display('AO slices requested using max_dims={:d}'.format(max_dims))
      display('Number of AOs required/total:    {:4d}/{:4d}'.format(len(AOs[0]), self.Norb))
      display('Number of shells required/total: {:4d}/{:4d}'.format(len(shells[0]), self.Nshells))
      display('AO slices generated: {:d}'.format(len(slices)))

      T0 = time.time()
      for i_, slice_ in enumerate(slices):
        display('\tcalculating slice {:2d}: {:4d} AOs, {:4d} shells'.format(
          i_+1,                                     # slice number
          sum(self.basis.shell_dims[slice_]),       # number of AOs in slice
          len(slice_),                              # number of shells in slice
        ))

        t0 = time.time()
        shells_sliced = [slice_] + shells[1:]

        # calculate AO matrix slice
        mat = _libcint2e(self.basis, operator, *shells_sliced)
        t1 = time.time()
        mat_bytes = mat.nbytes / 1000**2

        # remove orbitals not presents in AOs
        masks = []
        for i in range(4):
          # mask of all AOs
          mask = numpy.array([False]*self.Norb)
          # only enable request AOs
          mask[AOs[i]] = True
          # cut skipped shells
          mask = mask[sorted(chain(*[self._shell2ao(s_) for s_ in shells_sliced[i]]))]
          masks.append(mask)
        mat = mat[numpy.ix_(*masks)]

        display('\t\tAO integrals {:8.2f} sec, {:8.2f} MB'.format(
          t1-t0,     # time
          mat_bytes, # MB mat
        ))

        # convert to MO (blocks)
        # all AOs in current slice
        AOslice = sorted(set(chain(*[self._shell2ao(i) for i in slice_])).intersection(AOs[0]))
        # match indices with respect to AOs[0]
        AOslice = [numpy.where(numpy.array(AOs[0])>=ao)[0][0] for ao in AOslice]

        if not self.MO_blocks_2e:
          res, mem = _ao2mo_slice(mat, self.qc.mo_spec.coeffs[:,AOs[0]], AOslice=AOslice)
          results[0] += res
        for i, block in enumerate(self.MO_blocks_2e):
          mem = 0
          res, mem_ = _ao2mo_slice(
            mat, self.qc.mo_spec.coeffs[:,AOs[0]], AOslice=AOslice,
            MOrangei=block[0], MOrangej=block[1], MOrangek=block[2], MOrangel=block[3],
          )
          results[i] += res
          mem = max(mem, mem_)

        t2 = time.time()
        display('\t\tao2mo        {:8.2f} sec, {:8.2f} MB'.format(
          t2-t1,         # time
          mem / 1000**2   # MB ao2mo
        ))

      T1 = time.time()
      display('\tMO integrals:        {:8.2f} sec, {:8.2f} MB'.format(
        T1-T0,
        sum([r.nbytes for r in results]) / 1000**2,   # MB results
      ))

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
      return range(self.basis.shell_offsets[i], self.Norb)
    return range(self.basis.shell_offsets[i], self.basis.shell_offsets[i+1])

  def _ao2shell(self, i):
    '''Returns shell index for given AO.'''
    return numpy.where(self.basis.shell_offsets <= i)[0][-1]

  def _get_slices(self, max_dims, shells=None):
    '''Slices shells in separate ranges.'''
    if not shells:
      shells = range(self.Nshells)
    shells = list(shells)
    slices = []
    _slice = shells[:1]
    for s in shells[1:]:
      if sum(self.basis.shell_dims[_slice+[s]]) > max_dims:
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
      dims.append( fun(ctypes.c_int(ind), self.basis.c_bas) )
    if len(indices) == 1:
      return dims[0]
    return dims

  def _get_order(self, i):
    '''Returns indices to transform libcint order to orbkit order for given shell.'''
    l = self.basis.bas[i][1]
    if l == 0:
      return (0,)
    if self.cartesian:
      order_orbkit = self.qc.ao_spec.get_lxlylz()[self.qc.ao_spec.get_assign_lxlylz_to_cont()==i,:]
      order_libcint = []
      for item in combinations_with_replacement('xyz', l):
        order_libcint.append([item.count('x'), item.count('y'), item.count('z')])
      order_libcint = numpy.array(order_libcint)
    else:
      order_orbkit = numpy.array(self.qc.ao_spec.get_lm())[self.qc.ao_spec.get_assign_lm_to_cont()==i,1]
      if l == 1:
        order_libcint = numpy.array([1,-1,0])
      else:
        order_libcint = numpy.array(range(-l,l+1))
    return match_order(order_libcint, order_orbkit)

  def _norm_cart_shell(self):
    '''Calculates normalization factors for each shell. Required to rescale cartesian integrals.'''
    self.basis.cart_norm = []
    for i in range(self.basis.Nshells):
      # calculate square root of self overlap
      di = self._get_dims(i)
      mat = (ctypes.c_double * di*di)()
      shls = (ctypes.c_int * 2)(i, i)
      libcint.cint1e_ovlp_cart.restype = ctypes.c_void_p
      libcint.cint1e_ovlp_cart(mat, shls, self.basis.c_atm, self.basis.natm, self.basis.c_bas, self.basis.nbas, self.basis.c_env)
      S = numpy.reshape(mat, (di, di)).T
      S = numpy.sqrt(numpy.diag(S))
      self.basis.cart_norm.append(S)

def _libcint1e(basis, operator, shells_i, shells_j):
  '''Calculate AO integral matrix.

    **Parameters:**

    basis : structure keeping all relevent system data
    operator : libcint function basename to call
    shells_i, shells_j : shells of requested AO integrals along indices i and j

    **Returns:**

    2D array of AO integrals.
  '''

  # shell offsets (account for missing/skipped shells)
  shell_dims = basis.shell_dims[shells_i]
  off = _get_shell_offsets(shell_dims)
  shell_offsets_i = {I:off[i] for i, I in enumerate(shells_i)}

  shell_dims = basis.shell_dims[shells_j]
  off = _get_shell_offsets(shell_dims)
  shell_offsets_j = {I:off[i] for i, I in enumerate(shells_j)}

  # set up libcint function
  if basis.cartesian:
    ext = 'cart'
  else:
    ext = 'sph'
  operator = 'cint1e_%s_%s' %(operator, ext)
  operator = getattr(libcint, operator)
  operator.restype = ctypes.c_void_p

  # init AO array
  shape = []
  for shell in (shells_i, shells_j):
    Norb = numpy.sum(basis.shell_dims[shell[0]:shell[-1]+1])
    shape.append(Norb)
  mat = numpy.zeros(shape)

  # loop over shells
  for i in shells_i:
    for j in shells_j:

      # hermitian
      if i < j and i in shells_j and j in shells_i:
        continue

      # calc shell integrals
      di, dj = basis.shell_dims[[i, j]]

      buf = (ctypes.c_double * di*dj)()
      shls = (ctypes.c_int * 2)(i, j)

      operator(buf, shls, basis.c_atm, basis.natm, basis.c_bas, basis.nbas, basis.c_env)
      buf = numpy.reshape(buf, (dj, di)).T

      # cartesian integrals need to be rescaled according to overlap matrix
      if basis.cartesian:
        buf /= basis.cart_norm[i].reshape(di,1)
        buf /= basis.cart_norm[j].reshape(1,dj)

      # switch order of basis functions (angl>1)
      if di > 1:
        buf = buf[basis.order[i],:]
      if dj > 1:
        buf = buf[:,basis.order[j]]

      # store
      ii = shell_offsets_i[i]
      jj = shell_offsets_j[j]
      mat[ii:ii+di,jj:jj+dj] = buf
      if i > j and i in shells_j and j in shells_i:
        ii = shell_offsets_i[j]
        jj = shell_offsets_j[i]
        mat[ii:ii+dj,jj:jj+di] = buf.T

  return mat

def _libcint2e(basis, operator, shells_i, shells_j, shells_k, shells_l):
  '''Calculate AO integral matrix.

    **Parameters:**

    basis : structure keeping all relevent system data
    operator : libcint function basename to call
    shells_i, shells_j, shells_k, shells_l : shells of requested AO integrals along indices i, j, k and l

    **Returns:**

    4D array of AO integrals.
  '''

  # shell offsets (account for missing/skipped shells)
  shell_dims = basis.shell_dims[shells_i]
  off = _get_shell_offsets(shell_dims)
  shell_offsets_i = {I:off[i] for i, I in enumerate(shells_i)}

  shell_dims = basis.shell_dims[shells_j]
  off = _get_shell_offsets(shell_dims)
  shell_offsets_j = {I:off[i] for i, I in enumerate(shells_j)}

  shell_dims = basis.shell_dims[shells_k]
  off = _get_shell_offsets(shell_dims)
  shell_offsets_k = {I:off[i] for i, I in enumerate(shells_k)}

  shell_dims = basis.shell_dims[shells_l]
  off = _get_shell_offsets(shell_dims)
  shell_offsets_l = {I:off[i] for i, I in enumerate(shells_l)}

  # set up libcint function and optimizer
  if basis.cartesian:
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
  fopt(ctypes.byref(opt), basis.c_atm, basis.natm, basis.c_bas, basis.nbas, basis.c_env)
  #opt = ctypes.POINTER(ctypes.c_void_p)()  # disables optimizer

  # init AO array
  shape = []
  for shell in (shells_i, shells_j, shells_k, shells_l):
    Norb = numpy.sum(basis.shell_dims[shell[0]:shell[-1]+1])
    shape.append(Norb)
  mat = numpy.zeros(shape)

  # loop over shells
  for i in shells_i:
    for j in shells_j:
      for k in shells_k:
        for l in shells_l:

          # hermitian
          if i < j and i in shells_j and j in shells_i:
            continue
          if k < l and k in shells_l and l in shells_k:
            continue

          # exchange of electronic coordinates
          if i < k and i in shells_k and k in shells_i and j < l and j in shells_l and l in shells_j:
            continue
          if i < l and i in shells_l and l in shells_i and i == k and j > l and j in shells_l and l in shells_j:
            continue

          di, dj, dk, dl = basis.shell_dims[[i, j, k, l]]
          buf = (ctypes.c_double * di*dj*dk*dl)()
          shls = (ctypes.c_int * 4)(i, j, k, l)

          operator(buf, shls, basis.c_atm, basis.natm, basis.c_bas, basis.nbas, basis.c_env, opt)
          buf = numpy.reshape(buf, (dl, dk, dj, di))
          buf = moveaxis(buf, (0, 1, 2, 3), (3, 2, 1, 0))

          # cartesian integrals need to be rescaled according to overlap matrix
          if basis.cartesian:
            buf /= basis.cart_norm[i].reshape(di,1,1,1)
            buf /= basis.cart_norm[j].reshape(1,dj,1,1)
            buf /= basis.cart_norm[k].reshape(1,1,dk,1)
            buf /= basis.cart_norm[l].reshape(1,1,1,dl)

          # switch order of basis functions (angl>1)
          if di > 1:
            buf = buf[basis.order[i],:,:,:]
          if dj > 1:
            buf = buf[:,basis.order[j],:,:]
          if dk > 1:
            buf = buf[:,:,basis.order[k],:]
          if dl > 1:
            buf = buf[:,:,:,basis.order[l]]

          # store/replicate results

          # original block
          ii = shell_offsets_i[i]
          jj = shell_offsets_j[j]
          kk = shell_offsets_k[k]
          ll = shell_offsets_l[l]
          mat[ii:ii+di,jj:jj+dj,kk:kk+dk,ll:ll+dl] = buf

          # switch i <-> j (hermitian)
          if i > j and i in shells_j and j in shells_i:
            ii = shell_offsets_i[j]
            jj = shell_offsets_j[i]
            mat[ii:ii+dj,jj:jj+di,kk:kk+dk,ll:ll+dl] = numpy.swapaxes(buf, 0, 1)

          # switch k <-> l (hermitian)
          if k > l and k in shells_l and l in shells_k:
            ii = shell_offsets_i[i]
            jj = shell_offsets_j[j]
            kk = shell_offsets_k[l]
            ll = shell_offsets_l[k]
            mat[ii:ii+di,jj:jj+dj,kk:kk+dl,ll:ll+dk] = numpy.swapaxes(buf, 2, 3)

            # switch i <-> j and k <-> l (hermitian)
            if i > j and i in shells_j and j in shells_i:
              ii = shell_offsets_i[j]
              jj = shell_offsets_j[i]
              mat[ii:ii+dj,jj:jj+di,kk:kk+dl,ll:ll+dk] = moveaxis(buf, (0,1,2,3), (1,0,3,2))

          # switch i <-> k and j <-> l (exchange of electronic coordinates)
          if i > k and i in shells_k and k in shells_i and j > l and j in shells_l and l in shells_j:
            ii = shell_offsets_i[k]
            jj = shell_offsets_j[l]
            kk = shell_offsets_k[i]
            ll = shell_offsets_l[j]
            mat[ii:ii+dk,jj:jj+dl,kk:kk+di,ll:ll+dj] = moveaxis(buf, (0,1,2,3), (2,3,0,1))

            # switch  additionally i <-> j (hermitian)
            if i > j and i in shells_l and j in shells_k:
              kk = shell_offsets_k[j]
              ll = shell_offsets_l[i]
              mat[ii:ii+dk,jj:jj+dl,kk:kk+dj,ll:ll+di] = moveaxis(buf, (0,1,2,3), (3,2,0,1))

            # switch additionally k <-> l (hermitian)
            if k > l and k in shells_j and l in shells_i:
              ii = shell_offsets_i[l]
              jj = shell_offsets_j[k]
              kk = shell_offsets_k[i]
              ll = shell_offsets_l[j]
              mat[ii:ii+dl,jj:jj+dk,kk:kk+di,ll:ll+dj] = moveaxis(buf, (0,1,2,3), (2,3,1,0))

              # switch  additionally i <-> j and k <-> l (hermitian)
              if i > j and i in shells_l and j in shells_k:
                kk = shell_offsets_k[j]
                ll = shell_offsets_l[i]
                mat[ii:ii+dl,jj:jj+dk,kk:kk+dj,ll:ll+di] = moveaxis(buf, (0,1,2,3), (3,2,1,0))

  # release optimizer
  libcint.CINTdel_optimizer(ctypes.byref(opt))

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
  MOrange: list or range object or None
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
      return numpy.zeros((len(MOrangei), len(MOrangej), len(MOrangek), len(MOrangel))), 0

  if len(mat.shape) == 2:
    # 1-electron integrals
    mat = mat[numpy.ix_(AOs_mat,AOs)]
    return numpy.dot(coeffs[numpy.ix_(MOrangei, AOs_)], numpy.dot(mat, coeffs[numpy.ix_(MOrangej, AOs)].T))
  elif len(mat.shape) == 4:
    # 2-electron integrals
    mem = mat.nbytes
    mat = mat[numpy.ix_(AOs_mat, AOs, AOs, AOs)]
    mem = max(mem, mat.nbytes)
    mat = numpy.tensordot(mat, coeffs[numpy.ix_(MOrangej,AOs)], axes=(1, 1))
    mem = max(mem, mat.nbytes)
    mat = numpy.tensordot(mat, coeffs[numpy.ix_(MOrangek,AOs)], axes=(1, 1))
    mem = max(mem, mat.nbytes)
    mat = numpy.tensordot(mat, coeffs[numpy.ix_(MOrangel,AOs)], axes=(1, 1))
    mem = max(mem, mat.nbytes)
    mat = numpy.tensordot(mat, coeffs[numpy.ix_(MOrangei, AOs_)], axes=(0, 1))
    mat = numpy.rollaxis(mat, 3, 0)
    mem = max(mem, mat.nbytes)
    return mat, mem

def rescale_coeffs_libcint(exps, coeffs, l):
  return [ c*libcint.CINTgto_norm(l, ctypes.c_double(e)) for e, c in zip(exps, coeffs) ]

def _get_shell_offsets(shell_dims):
  offsets = numpy.zeros((len(shell_dims),), dtype=numpy.int)
  offsets[1:] = numpy.cumsum(shell_dims[:-1])
  return offsets
