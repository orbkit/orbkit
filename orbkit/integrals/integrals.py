import numpy
import ctypes

from ..tools import lquant

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
# - symmetry

class AOIntegrals():
  '''Interface to calculate AO Integrals with libcint.
  https://github.com/sunqm/libcint
  '''
  def __init__(self, qc, cartesian=True):
    self.qc = qc
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

    # get some parameters
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
    '''Shortcut to calculate one-electron Hamiltonian H_core = T_e + V_ne '''
    return self.int1e('kin', **kwargs) + self.int1e('nuc', **kwargs)

  def Vee(self, **kwargs):
    '''Shortcut to calculate electron-electron repulsion integrals.'''
    return self.int2e('',**kwargs)

  def int1e(self, operator, asMO=True, MOrange=None, max_dims=0):
    '''Calculates all one-electron integrals <i|operator|j>.

      **Parameters:**

        operator : string
          Base name of function/integral in libcint.
        asMO : boolean
          if True, transform from AO to MO basis.
        MOrange: list|range object|None
          only transform selected MOs
        max_dims : int
          if set, calculate AO in slices no larger than max_dims dimensions (but at least one shell). Slower, but needs less memory.

      **Returns:**

        Hermitian 2D array of integrals.
    '''
    #TODO: non-hermitian operators
    if self.cartesian:
      ext = 'cart'
    else:
      ext = 'sph'

    operator = 'cint1e_%s_%s' %(operator, ext)
    operator = getattr(libcint, operator)
    operator.restype = ctypes.c_void_p

    if not asMO or not max_dims:
      # single slice
      mat = numpy.zeros((self.Norb, self.Norb))
      for i in range(self.Nshells):
        for j in range(i,self.Nshells):
          ii, jj = self.shell_offsets[[i, j]]
          res = self._libcint1e(operator, i, j)
          di, dj = res.shape
          mat[ii:ii+di,jj:jj+dj] = res
          mat[jj:jj+dj,ii:ii+di] = res.transpose()
          jj += dj
        ii += di

      if asMO:
        mat = ao2mo(mat, self.qc.mo_spec.coeffs, MOrange=MOrange)

    else:
      # slice into smaller pieces to save some memory
      if MOrange:
        mat = numpy.zeros((len(MOrange), len(MOrange)))
      else:
        mat = numpy.zeros((self.Norb, self.Norb))
      slices = self._get_slices(max_dims)
      for s in slices:
        Norb = numpy.sum(self.shell_dims[s[0]:s[-1]+1])
        mat_s = numpy.zeros((Norb, self.Norb))
        for i in s:
          for j in range(self.Nshells):
            ii = self.shell_offsets[i] - self.shell_offsets[s[0]]
            jj = self.shell_offsets[j]
            res = self._libcint1e(operator, i, j)
            di, dj = res.shape
            mat_s[ii:ii+di,jj:jj+dj] = res

        if asMO:
          AOrange = range(self.shell_offsets[s[0]], self.shell_offsets[s[-1]]+self.shell_dims[s[-1]])
          mat_s = ao2mo(mat_s, self.qc.mo_spec.coeffs, AOrange=AOrange, MOrange=MOrange)
          mat += mat_s

    return mat

  def int2e(self, operator, asMO=True, MOrange=None, max_dims=0):
    '''Calculates all two-electron integrals <ij|operator|kl>.

      **Parameters:**

        operator : string
          Base name of function/integral in libcint.
        asMO : boolean
          if True, transform from AO to MO basis.
        MOrange: list|range object|None
          only transform selected MOs
        max_dims : int
          if set, calculate AO in slices no larger than max_dims dimensions (but at least one shell). Slower, but needs less memory.

      **Returns:**

        4D array of integrals.
    '''

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

    if not asMO or not max_dims:
      # single slice
      mat = numpy.zeros((self.Norb, self.Norb, self.Norb, self.Norb))
      for i in range(self.Nshells):
        for j in range(i,self.Nshells):      # hermitian
          for k in range(i,self.Nshells):    # exchange of electronic coordinates
            for l in range(k,self.Nshells):  # hermitian

              # get number of contractions and offsets for given shells
              di, dj, dk, dl = self.shell_dims[[i, j, k, l]]
              ii, jj, kk, ll = self.shell_offsets[[i, j, k, l]]

              # calculate integrals
              res = self._libcint2e(operator, i, j, k, l, opt)

              # store/replicate results
              mat[ii:ii+di,jj:jj+dj,kk:kk+dk,ll:ll+dl] = res

              mat[jj:jj+dj,ii:ii+di,kk:kk+dk,ll:ll+dl] = numpy.swapaxes(res, 0, 1)
              mat[ii:ii+di,jj:jj+dj,ll:ll+dl,kk:kk+dk] = numpy.swapaxes(res, 2, 3)
              mat[jj:jj+dj,ii:ii+di,ll:ll+dl,kk:kk+dk] = moveaxis(res, (0,1,2,3), (1,0,3,2))

              mat[kk:kk+dk,ll:ll+dl,ii:ii+di,jj:jj+dj] = moveaxis(res, (0,1,2,3), (2,3,0,1))
              mat[ll:ll+dl,kk:kk+dk,ii:ii+di,jj:jj+dj] = moveaxis(res, (0,1,2,3), (2,3,1,0))
              mat[kk:kk+dk,ll:ll+dl,jj:jj+dj,ii:ii+di] = moveaxis(res, (0,1,2,3), (3,2,0,1))
              mat[ll:ll+dl,kk:kk+dk,jj:jj+dj,ii:ii+di] = moveaxis(res, (0,1,2,3), (3,2,1,0))

      if asMO:
        mat = ao2mo(mat, self.qc.mo_spec.coeffs, MOrange=MOrange)

    else:
      # slice into smaller pieces to save some memory
      if MOrange:
        mat = numpy.zeros((len(MOrange), len(MOrange), len(MOrange), len(MOrange)))
      else:
        mat = numpy.zeros((self.Norb, self.Norb, self.Norb, self.Norb))
      slices = self._get_slices(max_dims)
      for s in slices:
        Norb = numpy.sum(self.shell_dims[s[0]:s[-1]+1])
        mat_s = numpy.zeros((Norb, self.Norb, self.Norb, self.Norb))
        for i in s:
          for j in range(self.Nshells):
            for k in range(self.Nshells):
              for l in range(k, self.Nshells):  # hermitian

                # get number of contractions and offsets for given shells
                di, dj, dk, dl = self.shell_dims[[i, j, k, l]]
                ii = self.shell_offsets[i] - self.shell_offsets[s[0]]
                jj = self.shell_offsets[j]
                kk = self.shell_offsets[k]
                ll = self.shell_offsets[l]

                # calculate integrals
                res = self._libcint2e(operator, i, j, k, l, opt)
                mat_s[ii:ii+di,jj:jj+dj,kk:kk+dk,ll:ll+dl] = res
                mat_s[ii:ii+di,jj:jj+dj,ll:ll+dl,kk:kk+dk] = numpy.swapaxes(res, 2, 3)

        if asMO:
          AOrange = range(self.shell_offsets[s[0]], self.shell_offsets[s[-1]]+self.shell_dims[s[-1]])
          mat_s = ao2mo(mat_s, self.qc.mo_spec.coeffs, AOrange=AOrange, MOrange=MOrange)
          mat += mat_s

    # release optimizer
    libcint.CINTdel_optimizer(ctypes.byref(opt))

    return mat

  ##############################
  ###  Interface to libcint  ###
  ##############################

  def _get_slices(self, max_dims):
    slices = []
    start, end = 0, 0
    while end < self.Nshells:
      end = numpy.where(numpy.cumsum(self.shell_dims[start:]) > max_dims)[0]
      if len(end) > 0:
        end = end[0] + start
      else:
        end = start + 1
      if end == start:
        end += 1
      slices.append(range(start,end))
      start = end
    return slices

  def _get_dims(self, *indices):
    '''Returns number of basis functions for shell indices'''
    if self.cartesian:
      fun = libcint.CINTcgto_cart
    else:
      fun = libcint.CINTcgto_spheric
    dims = []
    for ind in indices:
      dims.append( fun(ind, self.c_bas) )
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
      S = numpy.reshape(mat, (di, di)).transpose()
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
    mat = numpy.reshape(mat, (dj, di)).transpose()

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

def ao2mo(mat, coeffs, AOrange=None, MOrange=None):
  '''Transforms array of one- or two-electron integrals from AO to MO basis.

  **Parameters:**

    mat : 2 or 4 dimensional numpy.ndarray
      integrals to be transformed
    coeffs : 2 dimensional numpy.ndarry
      MO coefficients
    AOrange: list|range object|None
      only transform an AO slice
    MOrange: list|range object|None
      only transform selected MOs

  **Returns:**
    numpy.ndarray
  '''
  if MOrange is not None:
    coeffs = coeffs[MOrange,:]
  if len(mat.shape) == 2:
    # 1-electron integrals
    if AOrange is None:
      return numpy.dot(coeffs, numpy.dot(mat, coeffs.transpose()))
    else:
      return numpy.dot(coeffs[:,AOrange], numpy.dot(mat, coeffs.transpose()))
  elif len(mat.shape) == 4:
    # 2-electron integrals
    if AOrange is None:
      mat = numpy.tensordot(mat, coeffs, axes=(0, 1))  # i
      mat = numpy.tensordot(mat, coeffs, axes=(0, 1))  # j
      mat = numpy.tensordot(mat, coeffs, axes=(0, 1))  # k
      mat = numpy.tensordot(mat, coeffs, axes=(0, 1))  # l
    else:
      mat = numpy.tensordot(mat, coeffs[:,AOrange], axes=(0, 1))  # i
      mat = numpy.tensordot(mat, coeffs, axes=(0, 1))             # j
      mat = numpy.tensordot(mat, coeffs, axes=(0, 1))             # k
      mat = numpy.tensordot(mat, coeffs, axes=(0, 1))             # l
    return mat
  raise ValueError("'mat' musst be of size 2 or 4.")

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
