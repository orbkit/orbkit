import cython
from cython.parallel import prange
 
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

cdef extern from "math.h" nogil:
  double sin(double x)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef public sinwt(double t,
            const double omega,
            np.ndarray[double, ndim=1, mode="c"] x):
  '''Simple df/dt = sin(w*t + x) function intended for testing.'''
  cdef int i
  def np.ndarray[double, ndim=1, mode="c"] y
  for i in prange(x.shape[0], nogil=True):
    y[i] = sin(t * omega + x[i])
  return y
