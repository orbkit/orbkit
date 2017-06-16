from orbkit.detci.ci_read import psi4_detci, molpro_mcscf
from orbkit.test.tools import equal
from orbkit import options
import numpy
import os, inspect

tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
folder = os.path.join(tests_home, 'outputs_for_testing')

options.quiet = True

#I don't know what elese to check here...
#Someone elese please add an occupation number check

filetypes = ['psi4', 'molpro', 'molpro.tar']

readers = {'psi4': psi4_detci,
           'molpro': molpro_mcscf}

files = {'psi4': 'h2o.out',
         'tar': 'lih.out.tar.gz',
         'molpro': 'lih.out'}

ref_coeff = {'psi4': numpy.array([0.976736, -0.093030, -0.076255, 0.067785]),
             'molpro': numpy.array([0.9460981, -0.1770487, 0.1770487, -0.0378771])}

for ftype in filetypes:
  if '.' in ftype:
    ftype, filename = ftype.split('.')
    ciclass = readers[ftype](os.path.join(folder, files[filename]))
  ciclass = readers[ftype](os.path.join(folder, files[ftype]))
  equal(ciclass[0].coeffs[:4], ref_coeff[ftype])
