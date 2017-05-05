'''
Helper functions for systematic testing of Orbkit
'''

import numpy as np

def equal(a, b, tol=1e-5):
  if isinstance(a, (int, float)) and isinstance(b, (int, float)):
    assert abs(a-b) <= tol
  elif isinstance(a, np.ndarray) and isinstance(b, np.ndarray):
    assert np.allclose(a, b, rtol=tol*1e2, atol=tol)
  else:
    raise TypeError('Unsupported variable type - must be float, int, or np.adarray')
