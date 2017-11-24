from __future__ import division
import numpy
from .tools import require

def norm(a, b, cell, disp=[0,0,0]):
  '''calculates the distance between an atom in the primitive cell
     and another atom in an arbitrary cell'''
  return numpy.linalg.norm(a - (b + numpy.dot(cell, disp)))

def build_computational_supercell(cell, coord, cell, pbc_exp=[1,1,1], return_mapping=True):
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






