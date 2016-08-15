import cython

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

cdef extern from "math.h":
    double sqrt(double x)

cdef extern from "c_support.h":
  double ao_norm(int l,int m,int n,double alpha, int is_normalized)
  int doublefactorial(int n)

cdef extern from "c_non-grid-based.h":
  ctypedef struct S_Primitive:
          double alpha
          int l[3]
          double R[3]
  double get_overlap(S_Primitive *pA, S_Primitive *pB)

@cython.boundscheck(False)
@cython.wraparound(False)
def ommited_cca_norm(np.ndarray[int, ndim=2, mode="c"] lxlylz not None):
  """
  Get basis function normalization ommited in CCA standard. 
  
  Adaped from Psi4 Program code (FCHKWriter::write src/lib/libmints/writer.cc):

    // For Cartesian functions we need to add a basis function normalization
    // constant of
    //      _______________________________
    //     / (2lx-1)!! (2ly-1)!! (2lz-1)!!
    //    /  -----------------------------
    //  \/             (2l-1)!!
    //
    // which is omitted in the CCA standard, adopted by Psi4.
  """
  cdef int ao_num = lxlylz.shape[0]
  cdef double divisor = 0.0
  cdef np.ndarray[double, ndim=1, mode="c"] norm = np.zeros([ao_num],
                                                               dtype=np.float64)
  for i in range(ao_num):
    divisor = doublefactorial(2*(lxlylz[i,0] + lxlylz[i,1] + lxlylz[i,2]) - 1)
    norm[i] = sqrt( doublefactorial(2*lxlylz[i,0] - 1) 
                  * doublefactorial(2*lxlylz[i,1] - 1) 
                  * doublefactorial(2*lxlylz[i,2] - 1) 
                  / divisor
                  )
  return norm
  
@cython.boundscheck(False)
@cython.wraparound(False)
def aooverlap(np.ndarray[double, ndim=2, mode="c"] ra           not None,
              np.ndarray[double, ndim=2, mode="c"] rb           not None,
              np.ndarray[int,    ndim=2, mode="c"] lxlylz_a     not None,
              np.ndarray[int,    ndim=2, mode="c"] lxlylz_b     not None,
              np.ndarray[double, ndim=2, mode="c"] coeff_list   not None, 
              np.ndarray[int,    ndim=2, mode="c"] index        not None,
              int drv, int is_normalized):  
  """
  aooverlap(ra,rb,lxlylz_a,lxlylz_b,coeff_list,index,is_normalized)
  drv : type int
    0 - No derivative  
    1 - d/dx
    2 - d/dy
    3 - d/dz
  """
  cdef int nindex = index.shape[0]
  cdef int ao_num = lxlylz_a.shape[0]
  cdef np.ndarray[double, ndim=2, mode="c"] aoom = np.zeros([ao_num,ao_num],
                                                               dtype=np.float64)
  cdef np.ndarray[double, ndim=1, mode="c"] norm = np.zeros([nindex],
                                                               dtype=np.float64)
  cdef int i,j
  cdef int i_l, j_l, i_ao, j_ao,rr  
  cdef S_Primitive A,B
  
  for i in range(nindex):
    i_l = index[i,1]
    norm[i] = ao_norm(lxlylz_a[i_l,0],lxlylz_a[i_l,1],lxlylz_a[i_l,2],
                      coeff_list[i,0],is_normalized)
  
  for i in range(nindex):
    for j in range(nindex):
      i_ao = index[i,0]
      j_ao = index[j,0]
      i_l = index[i,1]
      j_l = index[j,1]
      
      for rr in range(3):
        A.R[rr] = ra[i_ao,rr]
        B.R[rr] = rb[j_ao,rr]
        A.l[rr] = lxlylz_a[i_l,rr]
        B.l[rr] = lxlylz_b[j_l,rr]
      
      A.alpha = coeff_list[i,0]
      B.alpha = coeff_list[j,0]
      if drv <= 0:
        aoom[i_l,j_l] += (coeff_list[i,1] * coeff_list[j,1] * norm[i] * norm[j]
                          * get_overlap(&A, &B))
      elif B.l[drv-1] == 0:
        B.l[drv-1] = lxlylz_b[j_l,drv-1] + 1
        aoom[i_l,j_l] += ((-2 * B.alpha) * coeff_list[i,1] * coeff_list[j,1] 
                          * norm[i] * norm[j] * get_overlap(&A, &B))
      else:
        B.l[drv-1] = lxlylz_b[j_l,drv-1] - 1
        aoom[i_l,j_l] += (lxlylz_b[j_l,drv-1] *coeff_list[i,1] * coeff_list[j,1] 
                          * norm[i] * norm[j] * get_overlap(&A, &B))
        B.l[drv-1] = lxlylz_b[j_l,drv-1] + 1
        aoom[i_l,j_l] += ((-2 * B.alpha) * coeff_list[i,1] * coeff_list[j,1] 
                          * norm[i] * norm[j] * get_overlap(&A, &B))
  return aoom
  
@cython.boundscheck(False)
@cython.wraparound(False)
def mooverlap(np.ndarray[double, ndim=1, mode="c"] mo_a           not None,
              np.ndarray[double, ndim=1, mode="c"] mo_b           not None,
              np.ndarray[double, ndim=2, mode="c"] aoom           not None
              ):
  cdef int nao = aoom.shape[0]
  cdef int mao = aoom.shape[1]
  cdef double overlap = 0.0
  cdef int k,l
  
  for k in range(nao):
    for l in range(nao):
      overlap += mo_a[k] * mo_b[l] * aoom[k,l]
  
  return overlap

@cython.boundscheck(False)
@cython.wraparound(False)
def mooverlapmatrix(np.ndarray[double, ndim=2, mode="c"] mo_a     not None,
                    np.ndarray[double, ndim=2, mode="c"] mo_b     not None,
                    np.ndarray[double, ndim=2, mode="c"] aoom     not None,
                    int i_start,
                    int i_end):
  cdef int nao = aoom.shape[0]
  cdef int mao = aoom.shape[1]
  cdef int nmo_a = mo_a.shape[0]
  cdef int nmo_b = mo_b.shape[0]
  cdef np.ndarray[double, ndim=2, mode="c"] moom = np.zeros([i_end-i_start,nmo_b],
                                                               dtype=np.float64)
  cdef double tmpsum
  cdef int i,j,k,l

  for i in range(i_start,i_end):
    for j in range(nmo_b):
      tmpsum = 0.0
      for k in range(nao):
        for l in range(mao):
          tmpsum += mo_a[i,k] * mo_b[j,l] * aoom[k,l]
      moom[i-i_start,j] = tmpsum
  
  return moom
