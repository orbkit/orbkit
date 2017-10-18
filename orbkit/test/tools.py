'''
Helper functions for systematic testing of Orbkit
'''

import numpy
from orbkit import QCinfo
from orbkit.orbitals import AOClass, MOClass

non_array_types = (int, float, numpy.float64, numpy.int, numpy.intc, numpy.int64)

def equal(a, b, tol=1e-5):
  if isinstance(a, list):
    a = numpy.array(a)
  if isinstance(b, list):
    b = numpy.array(b)
  if isinstance(a, non_array_types) and isinstance(b, non_array_types):
    assert abs(a-b) <= tol
  elif isinstance(a, numpy.ndarray) and isinstance(b, numpy.ndarray):
    if a.dtype in non_array_types and b.dtype in non_array_types:
      assert numpy.allclose(a, b, rtol=tol*1e2, atol=tol)
    else:
      for i,j in zip(a.flatten(), b.flatten()):
        assert i == j
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
