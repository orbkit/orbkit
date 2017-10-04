from orbkit.read.high_level import main_read
from orbkit.test.tools import equal
from orbkit.qcinfo import QCinfo
import numpy
from orbkit import options
import os, inspect

options.quiet = True

tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
folder = os.path.join(tests_home, '../outputs_for_testing/gaussian')
filepath = os.path.join(folder, 'h2o_uhf_cart.inp.log')

qc = main_read(filepath, all_mo=True)

item = qc.ao_spec.get_contspher()
qc.ao_spec.set_contspher(item)
equal(qc.ao_spec.get_contspher(), item)

item = qc.ao_spec.get_cont2atoms()
qc.ao_spec.set_cont2atoms(item)
equal(qc.ao_spec.get_cont2atoms(), item)

item = qc.ao_spec.get_prim2cont()
qc.ao_spec.set_prim2cont(item)
equal(qc.ao_spec.get_prim2cont(), item)

item = qc.ao_spec.get_lxlylz()
for i in range(item.shape[0]):
  item[i] = numpy.zeros(item.shape[1]) + i
qc.ao_spec.set_lxlylz(item)
equal(qc.ao_spec.get_lxlylz(), item)

item = qc.ao_spec.get_pao()
for i in range(item.shape[0]):
  item[i] = numpy.zeros(item.shape[1]) + i
qc.ao_spec.set_pao(item)
equal(qc.ao_spec.get_pao(), item)


