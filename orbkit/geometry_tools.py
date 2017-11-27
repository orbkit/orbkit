from __future__ import division
import numpy
from itertools import permutations, product
from .tools import require
from .display import display

def cellpar_to_cell(*args):
  '''Trubomole takes cell vector lengths and angles instead of a 3x3 cell

  **Parameters:**

    cell parameters and angles : floats
    You must first provide unit cell lengths BEFORE the angles between them.

  **Returns:**

    unit cell matrix : numpy.ndarray, shape (3,3)
  '''
  if len(args) < 1:
    raise ValueError('At least one value must be supplied to build a valid unit cell')

  cell = numpy.zeros((3,3))

  if len(args) not in [1,3,6]:
    raise ValueError('Wrong number of parameters supplied for cell construction')

  if len(args) == 1:
    display('Building 1D periodi cell')
    param = {'a': 0}
  elif len(args) == 3:
    display('Building 2D periodi cell')
    param = {'a': 0, 'b': 1, 'gamma': 2}
  else:
    display('Building 3D periodi cell')
    param = {'a': 0, 'b': 1, 'c': 2, 'alpha': 3, 'beta': 4, 'gamma': 5}

  cellpar = {}
  for p in param:
    if p in ['a','b','c']:
      cellpar[p] = float(args[param[p]])
    else:
      cellpar[p] = float(args[param[p]])

  # Code adopted from ase
  # Handle orthorhombic cells separately to avoid rounding errors
  cos_angle = {}

  eps = 2 * numpy.spacing(90.0, dtype=numpy.float64)  # around 1.4e-14
  for angle in cellpar.keys():
    if angle in ['alpha', 'beta']:
      if abs(abs(cellpar[angle]) - 90) < eps:
        cos_angle[angle] = 0.0
      else:
        cos_angle[angle] = numpy.cos(cellpar[angle] * numpy.pi / 180)
    elif angle == 'gamma':
      if abs(cellpar[angle] - 90) < eps:
          cos_angle[angle] = 0.0
          sin_gamma = 1.0
      elif abs(cellpar[angle] + 90) < eps:
          cos_angle[angle] = 0.0
          sin_gamma = -1.0
      else:
          cos_angle[angle] = numpy.cos(cellpar[angle] * numpy.pi / 180)
          sin_gamma = numpy.sin(cellpar[angle] * numpy.pi / 180)

  # Build the cell vectors
  cell = numpy.zeros((3,3))
  cell[0,0] = cellpar['a']
  if 'b' in cellpar.keys():
    cell[1,:] = cellpar['b'] * numpy.array([cos_angle['gamma'], sin_gamma, 0])
  if 'c' in cellpar.keys():
    cx = cos_angle['beta']
    cy = (cos_angle['alpha'] - cos_angle['beta'] * cos_angle['gamma']) / sin_gamma
    cz = (1. - cx * cx - cy * cy)**(1/2)
    cell[2,:] = cellpar['c'] * numpy.array([cx, cy, cz])
  return require(cell, dtype='f')

def norm(a, b, cell, disp=[0,0,0]):
  '''calculates the distance between an atom in the primitive cell
     and another atom in an arbitrary cell'''
  return numpy.linalg.norm(a - (b + numpy.dot(cell, disp)))

def build_computational_supercell(cell, coord, pbc_exp=[2,2,2], return_mapping=True):
  '''Builds up the supercell used in the AO overlap calculations
  and returns the mapping needed to project the supercell back
  to the primitive cell.

  **Parameters:**

  coord : numpy.ndarray
    Specifies the positions of the atoms.
  cell : numpy.ndarray
    Defines the unit cell. 
  pbc_exp : list, optioal
    Used to specify the maximum expansion of the unit cell along periodic directions.
    If not specified, only nearest-neighbor cells will be considered.
  return_mapping : bool, optioal
    Return mapping from supercell to primitive cell as well as the number of primitive cells
    within the supercell.

  **Returns:**

  comp_supercell : numpy.ndarray
    Contains the position of the atoms in the periodic supercell used to calculate
    AO overlap integrals
  supercell_mapping : list
    Contains the mapping between atoms in the supercell and atoms in the primitive cell.
    Only returned if ``return_mapping = True``
  nc : int
    Number of primitive cells within the supercell.
    Only returned if ``return_mapping = True``
  '''

  pbc = numpy.array([sum(v) for v in cell_a], dtype=bool)

  supercells = [numpy.zeros(3, dtype=numpy.intc)]
  for i, me in enumerate([pbc_exp[i] for i in range(3)]):
    disp_vec = numpy.zeros(3, dtype=numpy.intc)
    disp_vec[i] = 1
    supercells.extend([disp_vec*exp for exp in range(-me,me+1) if exp != 0])

  comp_supercell = []
  supercell_mapping = []
  for i in range(len(data)):
    map_to_i = []
    for j in range(len(data)):
      for c in supercells:
        if norm(i,j,c) <= cutoff:
          comp_supercell.append(self._data[j] + numpy.dot(self._cell, c))
          map_to_i.append(j)
    supercell_mapping.append(map_to_i)

  if return_mapping:
    return comp_supercell, supercell_mapping, nc
  else:
    return comp_supercell






