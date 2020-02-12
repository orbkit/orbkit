#cython: language_level=3
import cython
 
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

cdef extern from "math.h":
    double sqrt(double x)

@cython.boundscheck(False)
@cython.wraparound(False)
def cis_ab(int i,int j, 
           np.ndarray[double, ndim=1, mode="c"] acoeffs not None,
           np.ndarray[int, ndim=2, mode="c"] aocc not None,
           np.ndarray[double, ndim=1, mode="c"] bcoeffs not None,
           np.ndarray[int, ndim=2, mode="c"] bocc not None,
           np.ndarray[int, ndim=1, mode="c"] moocc not None):
  
  # Initialize the variables
  cdef double citmp = 0.0
  cdef double f=0.0
  cdef int ia,ib,im
  cdef int nbcoeffs = bcoeffs.shape[0]
  cdef int nmoocc = moocc.shape[0]
  zero = []     #: Identical Slater-determinants
  zc = []       #: Prefactor (occupation * CI coefficient)
  zm = []       #: Indices of occupied orbitals
  singles = []  #: Effective single excitation
  sc = []       #: Product of CI coefficients
  sm = []       #: Indices of the two molecular orbitals 
  cdef double sqrt2 = sqrt(2.0) 
  for ia in range(i,j):
    if acoeffs[ia] != 0.0:
      for ib in range(nbcoeffs):
        if bcoeffs[ib] != 0.0:
          if (aocc[ia,0] == bocc[ib,0] and aocc[ia,1] == bocc[ib,1]):
            si = []
            zo = []
            citmp = acoeffs[ia]*bcoeffs[ib]
            for im in range(nmoocc):
              if (aocc[ia,0] == im): f = (moocc[im]-1.0)*citmp
              elif (aocc[ia,1] == im):f = (moocc[im]+1.0)*citmp
              else: f = (moocc[im])*citmp
              if f != 0:
                si.append(im)
                zo.append(f)
            zc.append(zo)
            zm.append(si)
          elif (aocc[ia,0] == bocc[ib,0]):
            si = []
            sc.append(1.0*acoeffs[ia]*bcoeffs[ib])
            si.append(aocc[ia,1])
            si.append(bocc[ib,1])
            sm.append(si)
          elif (aocc[ia,1] == bocc[ib,1]):
            si = []
            sc.append(-1.0*acoeffs[ia]*bcoeffs[ib])
            si.append(bocc[ib,0])
            si.append(aocc[ia,0])
            sm.append(si)
          elif (aocc[ia,0] == -1):
            si = []
            sc.append(-sqrt2*bcoeffs[ib])
            si.append(bocc[ib,0])
            si.append(bocc[ib,1])
            sm.append(si)
          elif (bocc[ib,0] == -1):
            si = []
            sc.append(-sqrt2*acoeffs[ia])
            si.append(aocc[ia,1])
            si.append(aocc[ia,0])
            sm.append(si)

  zero.append(zc)
  zero.append(zm)

  singles.append(sc)
  singles.append(sm)

  return zero,singles

@cython.boundscheck(False)
@cython.wraparound(False)
def mcscf_ab(int i,int j, 
             np.ndarray[double, ndim=1, mode="c"] acoeffs not None,
             np.ndarray[int, ndim=3, mode="c"] aocc not None,
             np.ndarray[double, ndim=1, mode="c"] bcoeffs not None,
             np.ndarray[int, ndim=3, mode="c"] bocc not None,
             np.ndarray[int, ndim=1, mode="c"] moocc not None,
             int ab_sorting):
  '''Sorting of alpha and beta strings:
  if ab_sorting: Alternating alpha and beta orbitals (MOLPRO style of orbital sorting)
  else:          First alpha than beta orbitals (PSI4 style of orbital sorting)
  '''
  # Initialize the variables
  cdef double citmp = 0.0
  cdef int ia,ib,im,occ,itmp,jtmp,nflips
  cdef int nbcoeffs = bcoeffs.shape[0]
  cdef int nmoocc = moocc.shape[0] #: Number of core orbitals
  cdef int nactive = aocc.shape[1] #: Number of active orbitals
  cdef int alpha,beta
  cdef int warning = 0
  cdef int am,ap #: creation, annihilation 
  zero = []     #: Identical Slater-determinants
  zc = []       #: Prefactor (occupation * CI coefficient)
  zm = []       #: Indices of occupied orbitals
  singles = []  #: Effective single excitation
  sc = []       #: Product of CI coefficients
  sm = []       #: Indices of the two molecular orbitals 

  for ia in range(i,j):
    if acoeffs[ia] != 0.0:
      for ib in range(nbcoeffs):
        if bcoeffs[ib] != 0.0:
          alpha = 0
          beta = 0
          for im in range(nactive):
            alpha += abs(aocc[ia,im,0] - bocc[ib,im,0])
            beta += abs(aocc[ia,im,1] - bocc[ib,im,1])
          
          if (alpha+beta == 0):
            si = []
            zo = []
            citmp = acoeffs[ia]*bcoeffs[ib]
            for im in range(nmoocc):
              si.append(im)
              zo.append(moocc[im]*citmp)
            for im in range(nactive):
              occ = (aocc[ia,im,0]+aocc[ia,im,1])
              if occ != 0:
                si.append(nmoocc+im)
                zo.append(occ*citmp)
            zc.append(zo)
            zm.append(si)
          
          elif (alpha == 2 and beta == 0):
            am = -1
            ap = -1
            for im in range(nactive):
              alpha = aocc[ia,im,0] - bocc[ib,im,0]
              if alpha == -1: am = im
              elif alpha == 1: ap = im
            itmp = min(am,ap)
            jtmp = max(am,ap)
            if ab_sorting:
              nflips = aocc[ia,itmp,1]
              for im in range(itmp+1,jtmp):
                nflips += aocc[ia,im,0] + aocc[ia,im,1]
            else:
              nflips = 0
              for im in range(itmp+1,jtmp):
                nflips += aocc[ia,im,0]
            sc.append(pow(-1.0,nflips)*acoeffs[ia]*bcoeffs[ib])
            sm.append([ap+nmoocc, am+nmoocc])
            warning = min(min(am,ap),warning)
          
          elif (beta == 2 and alpha == 0):
            am = -1
            ap = -1
            for im in range(nactive):
              beta = aocc[ia,im,1] - bocc[ib,im,1]
              if beta == -1: am = im
              elif beta == 1: ap = im
            itmp = min(am,ap)
            jtmp = max(am,ap)
            if ab_sorting:
              nflips = aocc[ia,jtmp,0]
              for im in range(itmp+1,jtmp):
                nflips += aocc[ia,im,0] + aocc[ia,im,1]
            else:
              nflips = 0
              for im in range(itmp+1,jtmp):
                nflips += aocc[ia,im,1]
            sc.append(pow(-1.0,nflips)*acoeffs[ia]*bcoeffs[ib])
            sm.append([ap+nmoocc, am+nmoocc])
            warning = min(min(am,ap),warning)
  
  zero.append(zc)
  zero.append(zm)

  singles.append(sc)
  singles.append(sm)
  
  if warning < 0:
    print('Electron number mismatch between selected states!')
  return zero,singles
