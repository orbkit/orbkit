import numpy
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
    self.Tij = numpy.zeros((len(self.qc.mo_spec),len(self.qc.mo_spec)))
    #build the trace of Tij
    zero0 = zero[0].flatten()
    zero1 = zero[1].flatten()
    for imo, it in enumerate(zero1):
      self.Tij[it,it] += zero0[imo]
    #build the extradiagonal elements of Tij
    for imo, it in enumerate(sing[1]):
      self.Tij[it[0],it[1]] += sing[0][imo]

  def get_Diag(self):
    return numpy.diagonal(self.Tij)

class DMstates(UserList):
  def __init__(self, seq = [], fullci=None, qc=None):
    UserList.__init__(self, seq)
    self.qc = qc
    for i,j in product(range(len(fullci)), range(len(fullci))):
      zero, sing = occ_check.compare(fullci[i], fullci[j])
      self.append(DM(zero=zero, sing=sing, qc=qc))

class FullCI(UserList):
  '''This is a wrapper for CIinfo instances. Its main use
(for now) consists in the implementation of reduce density matrix based features'''
  def __init__(self, seq = [], qc=None):
    UserList.__init__(self, seq)
    self.qc = qc
    self.T = DMstates(fullci=self, qc=self.qc)











