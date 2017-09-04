import numpy

class Integrator:
  '''Wrapper class for Cython integrators'''
  def __init__(self, method, dt, y, prop, nsteps):
    self.method = method
    self.dt = dt
    self.y = y
    self.prop = prop
    self.nsteps = nsteps

  def propagate(self):
