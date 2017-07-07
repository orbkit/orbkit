'''
Helper functions for systematic testing of Orbkit
'''

import numpy
from orbkit import QCinfo
from orbkit.orbitals import AOClass, MOClass

def equal(a, b, tol=1e-5):
  if isinstance(a, list):
    a = numpy.array(a)
  if isinstance(b, list):
    b = numpy.array(b)
  if isinstance(a, (int, float)) and isinstance(b, (int, float)):
    assert abs(a-b) <= tol
  elif isinstance(a, numpy.ndarray) and isinstance(b, numpy.ndarray):
    assert numpy.allclose(a, b, rtol=tol*1e2, atol=tol)
  elif isinstance(a, QCinfo) and isinstance(b, QCinfo):
    assert a == b
  elif isinstance(a, AOClass) and isinstance(b, AOClass):
    assert a == b
  elif isinstance(a, MOClass) and isinstance(b, MOClass):
    assert a == b
  elif isinstance(a, dict) and isinstance(b, dict):
    for key in a.keys():
      assert numpy.allclose(a[key], b[key], rtol=tol*1e2, atol=tol)
  else:
    raise TypeError('Unsupported variable type - must be float, int, np.adarray, QCinfo, AOClass, and MOClass')
