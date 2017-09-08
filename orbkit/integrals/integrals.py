from orbkit.tools import lquant
import numpy
import ctypes

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

############################
###  external interface  ###
############################

class AOIntegrals():
  '''Interface to calculate AO Integrals with libcint.
  https://github.com/sunqm/libcint
  '''
  def __init__(self, qc):
    self.qc = qc
    self.cartesian = True # TODO: spherical
    ienv = 0
    self.env = [0]*ienv
    self.atm = []
    self.bas = []

    # build atm
    self.env.extend(qc.geo_spec.flatten())
    for center in qc.geo_info:
      # tuple for each center: (charge, coords offset, 0, 0, 0, 0)
      self.atm.append((float(center[2]), ienv, 0, 0, 0, 0))
      ienv += 3

    # build bas
    for ao in qc.ao_spec:

        # TODO: column-wise coefficients for contractions of same exponents
        l = lquant[ao['type']]
        exps, coeffs = zip(*ao['coeffs'])
        coeffs = rescale_coeffs(exps, coeffs, l)

        # add exponents and coefficients to env
        self.env.extend(exps)
        self.env.extend(coeffs)

        # tuple for each contraction:
        #  (center, l, #prim, #contr, kappa, prim offset, contr offset, 0)
        self.bas.append((ao['atom'], l, ao['pnum'], 1, 1, ienv, ienv+ao['pnum'], 0))
        ienv += 2*ao['pnum']

    self.env = numpy.array(self.env)
    self.atm = numpy.array(self.atm, dtype=numpy.int32)
    self.bas = numpy.array(self.bas, dtype=numpy.int32)

    self.c_env = self.env.ctypes.data_as(ctypes.c_void_p)
    self.c_atm = self.atm.ctypes.data_as(ctypes.c_void_p)
    self.c_bas = self.bas.ctypes.data_as(ctypes.c_void_p)

    self.natm = ctypes.c_int(self.atm.shape[0])
    self.nbas = ctypes.c_int(self.bas.shape[0])
    self.Norb = self.count()

  def count(self):
    '''Counts the number of contracted gaussians.'''
    if self.cartesian:
      fun = libcint.CINTcgto_cart
    else:
      fun = libcint.CINTcgto_spheric
    ncntr = 0
    for i in range(self.nbas.value):
      ncntr += fun(i, self.c_bas)
    return ncntr

  def count_primitives(self):
    '''Counts the number of primitive gaussians.'''
    if self.cartesian:
      fun = libcint.CINTcgto_cart
    else:
      fun = libcint.CINTcgto_spheric
    nprim = 0
    for i in range(self.nbas.value):
      nprim += self.bas[i,2]*self.bas[i,3]*fun(i, self.c_bas)
    return nprim
    return numpy.sum(self.bas[:,2]*self.bas[:,3])

  def get_dims(self, *indices):
    '''Returns number of basis functions for shell indices'''
    if self.cartesian:
      fun = libcint.CINTcgto_cart
    else:
      fun = libcint.CINTcgto_spheric
    dims = []
    for ind in indices:
      dims.append( fun(ind, self.c_bas) )
    return dims

  def overlap(self, asMO=True):
    '''Shortcut to calculate overlap integrals <i|j>.'''
    return self.int1e('ovlp', asMO)

  def kinetic(self, asMO=True):
    '''Shortcut to calculate kinetic energy integrals <i|-0.5*\nabla^2|j>.'''
    return self.int1e('kin', asMO)

  def Vne(self, asMO=True):
    '''Shortcut to calculate electron-nuclear repulsion integrals \sum_a <i|-Z_a/r|j>.'''
    return self.int1e('nuc', asMO)

  def Vee(self, asMO=True):
    '''Shortcut to calculate electron-electron repulsion integrals.'''
    return self.int2e('', asMO)

  def int1e(self, operator, asMO=True):
    '''Calculates one-electron integrals <i|operator|j>.

      **Parameters:**

        operator : Base name of function/integral in libcint.
        asMO : if True, transform from AO to MO basis.

      **Returns:**

        Hermitian 2D array of integrals.
    '''

    mat = numpy.zeros((self.Norb, self.Norb))

    ii = 0
    for i in range(self.nbas.value):
      jj = ii
      for j in range(i,self.nbas.value):
        res = self.libcint1e(operator, i, j)
        di, dj = res.shape
        mat[ii:ii+di,jj:jj+dj] = res
        mat[jj:jj+dj,ii:ii+di] = res.transpose()
        jj += dj
      ii += di

    if asMO:
      mat = ao2mo(mat, self.qc.mo_spec.coeffs)

    return mat

  def int2e(self, operator, asMO=True):
    '''Calculates two-electron integrals <ij|operator|kl>.

      **Parameters:**

        operator : Base name of function/integral in libcint.
        asMO : if True, transform from AO to MO basis.

      **Returns:**

        4D array of integrals.
    '''

    mat = numpy.zeros((self.Norb, self.Norb, self.Norb, self.Norb))

    ii = 0
    for i in range(self.nbas.value):
      jj = ii
      for j in range(i,self.nbas.value):
        kk = 0
        for k in range(self.nbas.value):
          ll = kk
          for l in range(k,self.nbas.value):

            # get number of contractions for given shell
            di, dj, dk, dl = self.get_dims(i, j, k, l)

            # exchange of electronic coordinates
            if (i > k):
              ll += dl
              continue

            # calculate integrals
            res = self.libcint2e(operator, i, j, k, l)

            # store/replicate results
            mat[ii:ii+di,jj:jj+dj,kk:kk+dk,ll:ll+dl] = res

            mat[jj:jj+dj,ii:ii+di,kk:kk+dk,ll:ll+dl] = numpy.swapaxes(res, 0, 1)
            mat[ii:ii+di,jj:jj+dj,ll:ll+dl,kk:kk+dk] = numpy.swapaxes(res, 2, 3)
            mat[jj:jj+dj,ii:ii+di,ll:ll+dl,kk:kk+dk] = moveaxis(res, (0,1,2,3), (1,0,3,2))

            mat[kk:kk+dk,ll:ll+dl,ii:ii+di,jj:jj+dj] = moveaxis(res, (0,1,2,3), (2,3,0,1))
            mat[ll:ll+dl,kk:kk+dk,ii:ii+di,jj:jj+dj] = moveaxis(res, (0,1,2,3), (2,3,1,0))
            mat[kk:kk+dk,ll:ll+dl,jj:jj+dj,ii:ii+di] = moveaxis(res, (0,1,2,3), (3,2,0,1))
            mat[ll:ll+dl,kk:kk+dk,jj:jj+dj,ii:ii+di] = moveaxis(res, (0,1,2,3), (3,2,1,0))

            ll += dl
          kk += dk
        jj += dj
      ii += di


    if asMO:
      mat = ao2mo(mat, self.qc.mo_spec.coeffs)

    return mat

  def libcint1e(self, operator, i, j):
    '''Calls libcint to evaluate 1-electron integrals over shells i and j.

    **Parameters:**

      operator : Name of function/integral in libcint.
      i, j : Shells for bra and ket vectors.

    **Returns:**

      2D matrix of integrals.

    '''

    # get libcint function
    if self.cartesian:
      ext = 'cart'
    else:
      ext = 'sph'

    foperator = 'cint1e_%s_%s' %(operator, ext)
    fun = getattr(libcint, foperator)

    # call libcint
    di, dj = self.get_dims(i, j)
    buf = (ctypes.c_double * di*dj)()
    shls = (ctypes.c_int * 2)(i, j)

    fun.restype = ctypes.c_void_p
    fun(buf, shls, self.c_atm, self.natm, self.c_bas, self.nbas, self.c_env)

    return numpy.array(buf).reshape(di, dj)

  def libcint2e(self, operator, i, j, k, l):
    '''Calls libcint to evaluate 2-electron integrals over shells i, j, k and l.

    **Parameters:**

      operator : Name of function/integral in libcint.
      i, j, k, l : Shells for bra and ket vectors.

    **Returns:**

      4D matrix of integrals.

    '''

    # get libcint function
    if self.cartesian:
      ext = 'cart'
    else:
      ext = 'sph'

    if operator:
      foperator = 'cint2e_%s_%s' %(operator, ext)
    else:
      foperator = 'cint2e_%s' %ext

    fun = getattr(libcint, foperator)
    opt = ctypes.POINTER(ctypes.c_void_p)()  # optimizer disabled

    # call libcint
    di, dj, dk, dl = self.get_dims(i, j, k, l)
    buf = (ctypes.c_double * di*dj*dk*dl)()
    shls = (ctypes.c_int * 4)(i, j, k, l)

    fun.restype = ctypes.c_void_p
    fun(buf, shls, self.c_atm, self.natm, self.c_bas, self.nbas, self.c_env, opt)

    # return output
    buf = numpy.array(buf).reshape(dl, dk, dj, di)
    buf = numpy.moveaxis(buf, (0, 1, 2, 3), (3, 2, 1, 0))
    return buf

def ao2mo(mat, coeffs):
  '''Transforms matrix of one- or two-electron integrals from AO to MO basis.'''
  if len(mat.shape) == 2:
    # 1-electron integrals
    return numpy.dot(coeffs, numpy.dot(mat, coeffs.transpose()))
  elif len(mat.shape) == 4:
    # 2-electron integrals
    coeffs = coeffs.transpose()
    D = numpy.tensordot(coeffs, coeffs, 0)
    mat = numpy.tensordot(mat, D, ((1,3), (0,2)))
    mat = numpy.tensordot(D, mat, ((0,2), (0,1)))
    return numpy.swapaxes(mat, 1, 2)
  raise ValueError("'mat' musst be of size 2 or 4.")

###########################################
###  Code parts from pyscf              ###
###  BSD 2-clause "Simplified" License  ###
###  pyscf/gto/mole.py                  ###
###########################################

from scipy.special import gamma

def _gaussian_int(n, alpha):
  r'''int_0^inf x^n exp(-alpha x^2) dx'''
  n1 = (n + 1) * .5
  return gamma(n1) / (2. * alpha**n1)

def gto_norm(l, expnt):
  r'''Normalized factor for GTO radial part   :math:`g=r^l e^{-\alpha r^2}`

  .. math::

    \frac{1}{\sqrt{\int g^2 r^2 dr}}
    = \sqrt{\frac{2^{2l+3} (l+1)! (2a)^{l+1.5}}{(2l+2)!\sqrt{\pi}}}

  Ref: H. B. Schlegel and M. J. Frisch, Int. J. Quant.  Chem., 54(1995), 83-87.

  Args:
    l (int):
      angular momentum
    expnt :
      exponent :math:`\alpha`

  Returns:
    normalization factor

  Examples:

  >>> print gto_norm(0, 1)
  2.5264751109842591
  '''
  if l >= 0:
    #f = 2**(2*l+3) * math.factorial(l+1) * (2*expnt)**(l+1.5) \
    #        / (math.factorial(2*l+2) * math.sqrt(math.pi))
    #return math.sqrt(f)
    return 1/numpy.sqrt(_gaussian_int(l*2+2, 2*expnt))
  else:
    raise ValueError('l should be > 0')

def rescale_coeffs(es, cs, angl):
  """rescale coefficients"""

  if not isinstance(es, numpy.ndarray):
    es = numpy.array(es)
  if not isinstance(cs, numpy.ndarray):
    cs = numpy.array(cs)
  if len(cs.shape) < 2:
    cs = cs.reshape(-1,1)

  cs = numpy.einsum('pi,p->pi', cs, gto_norm(angl, es))
  ee = es.reshape(-1,1) + es.reshape(1,-1)
  ee = _gaussian_int(angl*2+2, ee)
  s1 = 1/numpy.sqrt(numpy.einsum('pi,pq,qi->i', cs, ee, cs))
  cs = numpy.einsum('pi,i->pi', cs, s1)

  return cs.flatten().tolist()
