from __future__ import division
import numpy, scipy
from copy import deepcopy
from itertools import product

from orbkit.detci import occ_check

try:
  from UserList import UserList
except ImportError:
  from collections import UserList

class DM:
  def __init__(self, zero, sing, qc):
    zero = numpy.array(zero)
    self.qc = qc
    self.Tij = numpy.matrix(numpy.zeros((len(self.qc.mo_spec),len(self.qc.mo_spec))))
    #Ai and Ei contain eigenvalues and eigenvectors of Tij
    self.Ai = None
    self.Vi = None
    #build the trace of Tij
    zero0 = zero[0].flatten()
    zero1 = zero[1].flatten()
    for imo, it in enumerate(zero1):
      self.Tij[it,it] += zero0[imo]
    #build the extradiagonal elements of Tij
    for imo, it in enumerate(sing[1]):
      self.Tij[it[0],it[1]] += sing[0][imo]

  def diagonalize(self):
    if self.Ai is not None:
      self.Ai, self.Vi = scipy.linalg.eigh(self.Tij)
      self.Ai = require(self.Ai,dtype='f')
      self.Vi = require(self.Vi,dtype='f')

  def get_eigenvalues(self):
    return deepcopy(self.Ai)

  def get_entropy(self):
    '''Computes the von Neumann Entropy as S = -tr( Ai * log(Ai) ) where Ai is the vector
       of eigenvalues of Tij.'''
    self.diagonalize()
    Ai = deepcopy(self.Ai)
    if not self.qc.mo_spec.spinpolarized:
      Ai /= 2
    Ai = Ai[numpy.where(Ai > 0)]
    return -1 * sum(Ai * numpy.log(Ai))

class DMstates(UserList):
  def __init__(self, seq = [], fullci=None, qc=None):
    UserList.__init__(self, seq)
    self.qc = qc
    for i,j in product(range(len(fullci)), range(len(fullci))):
      zero, sing = occ_check.compare(fullci[i], fullci[j])
      self.append(DM(zero=zero, sing=sing, qc=qc))

class CIFock(UserList):
  '''This is a wrapper for CIinfo instances. Its main use
(for now) consists in the implementation of reduce density matrix based features'''
  def __init__(self, seq = [], qc=None):
    UserList.__init__(self, seq)
    self.qc = qc
    self.T = DMstates(fullci=self, qc=self.qc)











