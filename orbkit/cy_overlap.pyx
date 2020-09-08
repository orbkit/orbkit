#cython: language_level=3
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
def tmol_aomix_norm(np.ndarray[int, ndim=2, mode="c"] lxlylz not None):
  """
  Get basis function normalization used by Turbomole and orca
    // For Cartesian functions we need to add a basis function normalization
    // constant of
    //      _______________________________
    //    \/ (2lx-1)!! (2ly-1)!! (2lz-1)!!
  """
  cdef int ao_num = lxlylz.shape[0]
  cdef double divisor = 0.0
  cdef np.ndarray[double, ndim=1, mode="c"] norm = np.zeros([ao_num],
                                                               dtype=np.float64)
  for i in range(ao_num):
    norm[i] = sqrt( doublefactorial(2*lxlylz[i,0] - 1) 
                  * doublefactorial(2*lxlylz[i,1] - 1) 
                  * doublefactorial(2*lxlylz[i,2] - 1) 
                  )
  return norm
  
  
  
@cython.boundscheck(False)
@cython.wraparound(False)
def aooverlap(np.ndarray[double, ndim=2, mode="c"] geo_spec_a   not None,
              np.ndarray[double, ndim=2, mode="c"] geo_spec_b   not None,
              np.ndarray[int,    ndim=2, mode="c"] lxlylz_a     not None,
              np.ndarray[int,    ndim=2, mode="c"] lxlylz_b     not None,
              np.ndarray[int,    ndim=1, mode="c"] assign       not None,
              np.ndarray[double, ndim=2, mode="c"] ao_coeffs    not None, 
              np.ndarray[int,    ndim=1, mode="c"] pnum_list    not None, 
              np.ndarray[int,    ndim=1, mode="c"] atom_indices not None,
              int drv,
              int is_normalized): 
  """
  aooverlap(geo_spec_a,geo_spec_b,lxlylz_a,lxlylz_b,assign,ao_coeffs,pnum_list,atom_indices,drv,is_normalized)
  drv : type int
    0 - No derivative  
    1 - d/dx
    2 - d/dy
    3 - d/dz
  """
  cdef int nindex = sum(pnum_list*assign)
  cdef int ao_num = lxlylz_a.shape[0]
  cdef np.ndarray[double, ndim=2, mode="c"] aoom = np.zeros([ao_num,ao_num],
                                                               dtype=np.float64)
  cdef np.ndarray[double, ndim=1, mode="c"] norm = np.zeros([nindex],
                                                               dtype=np.float64)
  cdef np.ndarray[int,    ndim=2, mode="c"] index = np.zeros([nindex,3],
                                                               dtype=np.intc)
  cdef int i, j, rr
  cdef int c_ao = 0 # counter for contracted aos
  cdef int c_p = 0  # counter for primitves of each contracted ao
  cdef int c = 0 # absolute counter for primitves
  
  cdef int i_ao, i_p # indices for primitves of each contracted ao
  cdef int i_l, j_l  # absolute indices for primitves
  
  cdef S_Primitive A, B # structs for primitves
  
  for i in range(assign.shape[0]):    
    for i_ao in range(assign[i]):   
      for i_p in range(pnum_list[i]):  
        norm[c] = ao_norm(lxlylz_a[c_ao+i_ao,0],lxlylz_a[c_ao+i_ao,1],lxlylz_a[c_ao+i_ao,2],
                          ao_coeffs[c_p+i_p,0],is_normalized)
        index[c,0] = c_p+i_p
        index[c,1] = c_ao+i_ao
        index[c,2] = atom_indices[i]
        c += 1
    c_ao += assign[i]
    c_p += pnum_list[i]
  
  # Rearange the data
  ao_coeffs  = np.copy(ao_coeffs[index[:,0]])
  lxlylz_a   = np.copy(lxlylz_a[index[:,1]])
  lxlylz_b   = np.copy(lxlylz_b[index[:,1]])
  geo_spec_a = np.copy(geo_spec_a[index[:,2]])
  geo_spec_b = np.copy(geo_spec_b[index[:,2]])
  
  
  for i in range(nindex):
    for j in range(nindex):
      i_l = index[i,1]
      j_l = index[j,1]
      
      for rr in range(3):
        A.R[rr] = geo_spec_a[i,rr]
        B.R[rr] = geo_spec_b[j,rr]
        A.l[rr] = lxlylz_a[i,rr]
        B.l[rr] = lxlylz_b[j,rr]
      
      A.alpha = ao_coeffs[i,0]
      B.alpha = ao_coeffs[j,0]
      
      if drv <= 0:
        aoom[i_l,j_l] += (ao_coeffs[i,1] * ao_coeffs[j,1] * norm[i] * norm[j]
                          * get_overlap(&A, &B))
      elif B.l[drv-1] == 0:
        B.l[drv-1] = lxlylz_b[j,drv-1] + 1
        aoom[i_l,j_l] += ((-2 * B.alpha) * ao_coeffs[i,1] * ao_coeffs[j,1] 
                          * norm[i] * norm[j] * get_overlap(&A, &B))
      else:
        B.l[drv-1] = lxlylz_b[j,drv-1] - 1
        aoom[i_l,j_l] += (lxlylz_b[j,drv-1] *ao_coeffs[i,1] * ao_coeffs[j,1] 
                          * norm[i] * norm[j] * get_overlap(&A, &B))
        B.l[drv-1] = lxlylz_b[j,drv-1] + 1
        aoom[i_l,j_l] += ((-2 * B.alpha) * ao_coeffs[i,1] * ao_coeffs[j,1] 
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
