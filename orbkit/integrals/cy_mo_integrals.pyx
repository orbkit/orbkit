import numpy
cimport numpy

def normalize_indices(*indices):
  '''Transforms indices for repeating integral.s'''
  cdef int i, j, k, l
  if len(indices) == 2:
    i, j = indices
    if i > j:
      i, j = j, i
    return i, j
  elif len(indices) == 4:
    i, j, k, l = indices
    # for real basis functions
    if i > j:
      i, j = j, i
    if k > l:
      k, l = l, k
    # exchange of electronic coordinates
    if (i > k) or (i == k and j > l):
      i, j, k, l = k, l, i, j
    return i, j, k, l
  raise ValueError("indices must be of lenght 2 or 4")

def indices_2e(int Norb, disable_exchange=False):
  '''Generator for normalized indices for 2-electron integrals.'''
  cdef int i, j, k, l
  for i in range(Norb):
    for j in range(i, Norb):
      for k in range(Norb):
        for l in range(k, Norb):
          # exchange of electronic coordinates
          if not disable_exchange and ((i > k) or (i == k and j > l)):
            continue
          yield i, j, k, l

def redundant_indices(int i, int j, int k, int l):
  '''Return redundant indices for normalized 2-electron integrals.'''
  # TODO: skip for equal indices, e.g. i=j=k=l=0
  yield j, i, k, l
  yield i, j, l, k
  yield j, i, l, k
  yield k, l, i, j
  yield l, k, i, j
  yield k, l, j, i
  yield l, k, j, i

def expand_G(numpy.ndarray[numpy.float_t, ndim=4] G):
  '''copy upper triangular part into lower triangular part'''
  cdef int i, j, k, l
  for i, j, k, l in indices_2e(G.shape[0]):
    for indices in redundant_indices(i, j, k, l):
      G[indices] = G[i,j,k,l]

def transform_F(numpy.ndarray[numpy.float_t, ndim=2] Falpha, numpy.ndarray[numpy.float_t, ndim=2] Fbeta):
  '''Transforms one-electron integrals from spatial to spin orbital basis.'''

  cdef int P, Q, p, q
  cdef int Norb = 2*Falpha.shape[0]
  cdef numpy.ndarray[numpy.float_t, ndim=2] f = numpy.zeros((Norb,Norb))

  for P in range(Norb):
    for Q in range(P, Norb):
      # check spin: <alpha beta | alpha beta> = 0
      if P%2 != Q%2:
        f[P,Q] = 0.0
      else:
        # map spin orbitals to spatial orbitals using integer division
        p, q = P//2, Q//2
        if P%2:
          f[P,Q] = Fbeta[p,q]
        else:
          f[P,Q] = Falpha[p,q]
  return f

def transform_G(numpy.ndarray[numpy.float_t, ndim=4] Galpha, numpy.ndarray[numpy.float_t, ndim=4] Gbeta, numpy.ndarray[numpy.float_t, ndim=4] Galphabeta):
  '''Transforms two-electron integrals from spatial to spin orbital basis.'''

  cdef int P, Q, R, S, p, q, r, s
  cdef int Norb = 2*Galpha.shape[0]
  cdef numpy.ndarray[numpy.float_t, ndim=4] g = numpy.zeros((Norb,Norb,Norb,Norb))
  cdef str spin1, spin2

  for P, Q, R, S in indices_2e(Norb):
    # check spin: <alpha beta | alpha beta> = 0
    if (P%2 != Q%2) or (R%2 != S%2):
      g[P,Q,R,S] = 0.0
    else:
      # map spin orbitals to spatial orbitals using integer division
      p, q, r, s = normalize_indices(P//2, Q//2, R//2, S//2)
      spin1 = ('alpha', 'beta')[P%2]
      spin2 = ('alpha', 'beta')[R%2]
      spin1 = ('alphabeta', spin1)[spin1==spin2]
      if spin1 == 'alpha':
        g[P,Q,R,S] = Galpha[p,q,r,s]
      elif spin1 == 'beta':
        g[P,Q,R,S] = Gbeta[p,q,r,s]
      else:
        g[P,Q,R,S] = Galphabeta[p,q,r,s]

  return g
