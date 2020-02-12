#cython: language_level=3
import cython
 
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

cdef extern from "math.h":
    double sin(double x)
    double cos(double x)

@cython.boundscheck(False)
@cython.wraparound(False)
def grid2vector(np.ndarray[double, ndim=1, mode="c"] x not None,
                   np.ndarray[double, ndim=1, mode="c"] y not None,
                   np.ndarray[double, ndim=1, mode="c"] z not None):
  cdef int npts = x.shape[0] * y.shape[0] * z.shape[0]
  cdef np.ndarray[double, ndim=2, mode="c"] xyz = np.zeros([3,npts],
                                                            dtype=np.float64)
  cdef int i,j,k
  cdef int c = 0
  for i in range(x.shape[0]):
    for j in range(y.shape[0]):
      for k in range(z.shape[0]):
        xyz[0,c] = x[i]
        xyz[1,c] = y[j]
        xyz[2,c] = z[k]
        c += 1
  
  return xyz

@cython.boundscheck(False)
@cython.wraparound(False)
def vector2grid(np.ndarray[double, ndim=1, mode="c"] x not None,
                   np.ndarray[double, ndim=1, mode="c"] y not None,
                   np.ndarray[double, ndim=1, mode="c"] z not None,
                   int Nx, int Ny, int Nz):
  cdef int npts = x.shape[0]
  cdef np.ndarray[double, ndim=1, mode="c"] xnew = np.zeros([Nx],
                                                            dtype=np.float64)
  cdef np.ndarray[double, ndim=1, mode="c"] ynew = np.zeros([Ny],
                                                            dtype=np.float64)
  cdef np.ndarray[double, ndim=1, mode="c"] znew = np.zeros([Nz],
                                                            dtype=np.float64)
  cdef int i,j,k
  # x runs slowest
  for i in range(Nx):
    xnew[i] = x[i * Ny * Nz]
  for j in range(Ny):
    ynew[j] = y[j * Nz]
  # z runs fastest
  for k in range(Nz):
    znew[k] = z[k]
  
  return xnew,ynew,znew

@cython.boundscheck(False)
@cython.wraparound(False)
def sph2cart(np.ndarray[double, ndim=1, mode="c"] r not None,
                np.ndarray[double, ndim=1, mode="c"] theta not None,
                np.ndarray[double, ndim=1, mode="c"] phi not None):
  cdef int npts = r.shape[0] * theta.shape[0] * phi.shape[0]
  cdef np.ndarray[double, ndim=2, mode="c"] xyz = np.zeros([3,npts],
                                                            dtype=np.float64)
  cdef int i,j,k
  
  cdef int c = 0
  for i in range(r.shape[0]):
    for j in range(theta.shape[0]):
      for k in range(phi.shape[0]):
        xyz[0,c] = r[i] * sin(theta[j]) * cos(phi[k])
        xyz[1,c] = r[i] * sin(theta[j]) * sin(phi[k])
        xyz[2,c] = r[i] * cos(theta[j])
        c += 1
  
  return xyz

@cython.boundscheck(False)
@cython.wraparound(False)
def cyl2cart(np.ndarray[double, ndim=1, mode="c"] r not None,
                np.ndarray[double, ndim=1, mode="c"] phi not None,
                np.ndarray[double, ndim=1, mode="c"] zed not None):
  cdef int npts = r.shape[0] * phi.shape[0] * zed.shape[0]
  cdef np.ndarray[double, ndim=2, mode="c"] xyz = np.zeros([3,npts],
                                                            dtype=np.float64)
  cdef int i,j,k
  
  cdef int c = 0
  for i in range(r.shape[0]):
    for j in range(phi.shape[0]):
      for k in range(zed.shape[0]):
        xyz[0,c] = r[i] * cos(phi[j])
        xyz[1,c] = r[i] * sin(phi[j])
        xyz[2,c] = zed[k]
        c += 1
  
  return xyz
  