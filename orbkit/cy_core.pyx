#cython: language_level=3
import cython
 
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

cdef extern from "math.h":
    double sqrt(double x)
    double pow(double x,int y)

cdef extern from "c_grid-based.h":
  void c_lcreator(double* ao_list, int* lxlylz, double* coeff_list,
                    double* at_pos, double* x, double* y, double* z, 
                    int npts, int ao_num, int pnum, int drv, int is_normalized)

cdef extern from "c_support.h":
  double ao_norm(int l,int m,int n,double alpha, int is_normalized)
  double get_ao_xyz(double X, double Y, double Z, int lx, int ly, int lz, double alpha, int drv)

def aonorm(int lx,int ly,int lz,double alpha, int is_normalized):
  return ao_norm(lx,ly,lz,alpha,is_normalized)

def aoxyz(double x, double y, double z, int lx, int ly, int lz, double alpha, int drv):
  return get_ao_xyz(x, y, z, lx,ly,lz,alpha,drv)

@cython.boundscheck(False)
@cython.wraparound(False)
def lcreator(np.ndarray[double, ndim=2, mode="c"] ao_list       not None,
             np.ndarray[int,    ndim=2, mode="c"] lxlylz        not None,
             np.ndarray[double, ndim=2, mode="c"] coeff_list    not None, 
             np.ndarray[double, ndim=1, mode="c"] at_pos        not None,
             np.ndarray[double, ndim=1, mode="c"] x             not None,
             np.ndarray[double, ndim=1, mode="c"] y             not None,
             np.ndarray[double, ndim=1, mode="c"] z             not None,
             int ao_num,
             int pnum,
             int drv,
             int is_normalized):  
  """
  lcreator(ao_list,lxlylz,coeff_list,at_pos,x,y,z,ao_num,pnum,drv,is_normalized)
  """
  cdef int npts = x.shape[0]
  
  c_lcreator(&ao_list[0,0],&lxlylz[0,0],&coeff_list[0,0],
                         &at_pos[0],&x[0],&y[0],&z[0],npts,ao_num,pnum,
                         drv,is_normalized)

@cython.boundscheck(False)
@cython.wraparound(False)
def aocreator(np.ndarray[int,    ndim=2, mode="c"] lxlylz       not None,
              np.ndarray[int,    ndim=1, mode="c"] assign       not None,
              np.ndarray[double, ndim=2, mode="c"] ao_coeffs    not None, 
              np.ndarray[int,    ndim=1, mode="c"] pnum_list    not None,
              np.ndarray[double, ndim=2, mode="c"] geo_spec     not None, 
              np.ndarray[int,    ndim=1, mode="c"] atom_indices not None,
              np.ndarray[double, ndim=1, mode="c"] x            not None,
              np.ndarray[double, ndim=1, mode="c"] y            not None,
              np.ndarray[double, ndim=1, mode="c"] z            not None,
              int drv,
              int is_normalized):  
  """
  aocreator(lxlylz,assign,ao_coeffs,pnum_list,geo_spec,atom_indices,x,y,z,drv,is_normalized)
  """
  cdef int npts = x.shape[0]
  cdef int ao_num = lxlylz.shape[0]
  cdef np.ndarray[double, ndim=2, mode="c"] ao_list = np.zeros([ao_num,npts],
                                                               dtype=np.float64)
  cdef int i
  cdef int c_ao = 0 # counter for aos
  cdef int c_p = 0  # counter for primitves 
  for i in range(assign.shape[0]):    
    lcreator(ao_list[c_ao:,:],lxlylz[c_ao:,:],ao_coeffs[c_p:,:],
             geo_spec[atom_indices[i],:],x,y,z,assign[i],pnum_list[i],
             drv,is_normalized)
    c_ao += assign[i]
    c_p += pnum_list[i]
  return ao_list
  
@cython.boundscheck(False)
@cython.wraparound(False)
def mocreator(np.ndarray[double, ndim=2, mode="c"] ao_list      not None,
              np.ndarray[double, ndim=2, mode="c"] mo_coeffs    not None,):  
  """
  mocreator(ao_list,mo_coeffs)
  """
  cdef int ao_num = ao_list.shape[0]
  cdef int npts   = ao_list.shape[1]
  cdef int mo_num = mo_coeffs.shape[0]
  cdef double value = 0.0
  cdef np.ndarray[double, ndim=2, mode="c"] mo_list = np.zeros([mo_num,npts],
                                                               dtype=np.float64)
  cdef int i,j,k
  for i in range(mo_num):
    for j in range(npts):
      value = 0.0
      for k in range(ao_num):
        value += mo_coeffs[i,k] * ao_list[k,j]
      mo_list[i,j] = value
  
  return mo_list
