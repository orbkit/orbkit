import cython
 
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np
from scipy.linalg.blas cimport caxpy

cdef extern from cy_functions.h:
  cpdef np.ndarray[double, ndim=1, mode="c"] sinwt(double t, const double omega, np.ndarray[double, ndim=1, mode="c"] x)

propagators = {'sinwt': sinwt}

@cython.boundscheck(False)
@cython.wraparound(False)
cdef RungeKutta(int nsteps, double dt, char* prop,
                np.ndarray[double, ndim=1, mode="c"] y,
                np.ndarray[double, ndim=1, mode="c"] x):
  integrator = propagators[prop]
  cdef int i, j
  cdef double dt2, dt6
  dt2 = dt * 0.5
  dt6 = dt / 6.
  cdef np.ndarray[double, ndim=1] yt = \
                 np.zeros(y.shape[0], dtype=np.float64)
  cdef np.ndarray[double, ndim=1] k2 = \
                  np.zeros(y.shape[0], dtype=np.float64)
  cdef np.ndarray[double, ndim=1] k3 = \
                  np.zeros(y.shape[0], dtype=np.float64)
  cdef np.ndarray[double, ndim=1] k4 = \
                  np.zeros(y.shape[0], dtype=np.float64)

  for i in range(nsteps):
    yt = caxpy(y + dt * integrator)
