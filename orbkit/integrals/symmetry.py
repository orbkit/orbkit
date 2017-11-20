# -*- coding: utf-8 -*-
import numpy
from functools import reduce

## Ordering of IRREPs:
# Molpro Ordering: http://www.molpro.net/info/current/doc/manual/node36.html
# Cotton Ordering: http://www.psicode.org/psi4manual/master/psithonmol.html#table-irrepordering

irrep_labels = {
  # cotton ordering
  'd2h' : ('', 'Ag', 'B1g', 'B2g', 'B3g', 'Au', 'B1u', 'B2u', 'B3u'),
  'c2v' : ('', 'A1', 'A2', 'B1', 'B2'),
  'c2h' : ('', 'Ag', 'Bg', 'Au', 'Bu'),
   'd2' : ('', 'A', 'B1', 'B2', 'B3'),
   'cs' : ('', 'A′', 'A″'),
   'c2' : ('', 'A', 'B'),
   'ci' : ('', 'Ag', 'Au'),
   'c1' : ('', 'A',),
}

cotton2molpro = {
  'd2h' : [0, 1, 4, 6, 7, 8, 5, 3, 2],
  'c2v' : [0, 1, 4, 2, 3],
  'c2h' : [0, 1, 4, 2, 3],
  'd2'  : [0, 1, 4, 3, 2],
  'cs'  : [0, 1, 2],
  'c2'  : [0, 1, 2],
  'ci'  : [0, 1, 2],
  'c1'  : [0, 1],
}

point_groups = irrep_labels.keys()
Nirreps = {sym:len(irrep_labels[sym])-1 for sym in point_groups}

def irrep_label(point_group, irrep, ordering='molpro'):
  point_group = point_group.lower()
  if ordering == 'molpro':
    irrep = cotton2molpro[point_group].index(irrep)
  return irrep_labels[point_group][irrep]

# works for all IRREPs and both orderings:
multiplication_table = numpy.array([
  (1, 2, 3, 4, 5, 6, 7, 8),
  (2, 1, 4, 3, 6, 5, 8, 7),
  (3, 4, 1, 2, 7, 8, 5, 6),
  (4, 3, 2, 1, 8, 7, 6, 5),
  (5, 6, 7, 8, 1, 2, 3, 4),
  (6, 5, 8, 7, 2, 1, 4, 3),
  (7, 8, 5, 6, 3, 4, 1, 2),
  (8, 7, 6, 5, 4, 3, 2, 1),
])

def irrep_mul(*irreps):
  '''Returns the product of a list of irreps.'''
  return reduce(lambda i,j: multiplication_table[i-1, j-1], irreps)

labels_psi4 = {
  'd2h' : ('', 'Ag', 'B1g', 'B2g', 'B3g', 'Au', 'B1u', 'B2u', 'B3u'),
  'c2v' : ('', 'A1', 'A2', 'B1', 'B2'),
  'c2h' : ('', 'Ag', 'Bg', 'Au', 'Bu'),
   'd2' : ('', 'A', 'B1', 'B2', 'B3'),
   'cs' : ('', 'A\'', 'A"'),
   'c2' : ('', 'A', 'B'),
   'ci' : ('', 'Ag', 'Au'),
   'c1' : ('', 'A',),
}

def parse_irrep(label, point_group):
  '''Returns IRREP index.'''

  # try Molpro
  if label.isdigit():
    return int(label)

  # try Psi4
  labels = labels_psi4[point_group]
  if label in labels:
    return labels.index(label)

  raise ValueError('unkown symmetry label {} for point group {}'.format(label, point_group.capitalize()))

