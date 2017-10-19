from orbkit.read.high_level import main_read
from orbkit.test.tools import equal
from orbkit.qcinfo import QCinfo
import numpy
from orbkit import options
import os, inspect

options.quiet = True

tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
folder = os.path.join(tests_home, '../outputs_for_testing/gamess')
filepath = os.path.join(folder, 'formaldehyde.log')

qc = main_read(filepath, all_mo=True)

item = qc.mo_spec.get_coeffs()
for i in range(item.shape[0]):
  item[i] = numpy.zeros(item.shape[1]) + i
qc.mo_spec.set_coeffs(item)
equal(qc.mo_spec.get_coeffs(), item)

item = qc.mo_spec.get_occ()
qc.mo_spec.set_occ(item)
equal(qc.mo_spec.get_occ(), item)

item = qc.mo_spec.get_eig()
qc.mo_spec.set_eig(item)
equal(qc.mo_spec.get_eig(), item)

item = qc.mo_spec.get_sym()
qc.mo_spec.set_sym(item)
equal(qc.mo_spec.get_sym(), item)
