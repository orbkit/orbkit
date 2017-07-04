import os, inspect
import numpy

from orbkit import multiple_files as mult
from orbkit import options
from orbkit.test.tools import equal

options.quiet = True

tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
folder = os.path.join(tests_home, 'outputs_for_testing')
filepath = os.path.join(folder, 'NaCl_molden_files.tar.gz')

mult.read(filepath, all_mo=True, nosym=False)

index_list, mo_overlap = mult.order_using_analytical_overlap(None)

#numpy.savez('refdata_read_multiple_files.npz', index_list=index_list[0][0], mo_overlap=mo_overlap[0][0])

filepath = os.path.join(tests_home, 'refdata_read_multiple_files.npz')
with open(filepath, 'r') as fd:
  refdata = numpy.load(fd)
  equal(index_list[0][0], refdata['index_list'])
  equal(mo_overlap[0][0], refdata['mo_overlap'])
