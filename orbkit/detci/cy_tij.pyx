import cython
from cython.parallel import prange
 
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

'''Implementation of cy_ci using the reduced density matrix'''

@cython.boundscheck(False)
@cython.wraparound(False)
def get_rho(int i,int j, 
           np.ndarray[dm, ndim=2, mode="c"] dm not None,
           np.ndarray[double, ndim=2, mode="c"] molist not None):

  # Initialize the variables
  cdef int slen = j-i
  cdef np.ndarray[double, ndim=1, mode="c"] rho = np.zeros([slen],dtype=np.float64)
  cdef int sta,stb
  cdef int iz,jz,ks,x
  cdef double citmp = 0.0
  cdef int nz = len(zero[0])
  cdef int ns = len(sing[0])

  for iz in range(nz):
    for jz in range(len(zero[0][iz])):
      sta = zero[1][iz][jz]
      citmp = zero[0][iz][jz]   
      for x in range(slen):
        rho[x] += citmp*molist[sta,(x+i)]*molist[sta,(x+i)]
  for ks in range(ns):
    sta = sing[1][ks][0]
    stb = sing[1][ks][1]
    citmp = sing[0][ks]
    for x in range(slen):
      rho[x] += citmp*molist[sta,(x+i)]*molist[stb,(x+i)]
  
  return rho 
