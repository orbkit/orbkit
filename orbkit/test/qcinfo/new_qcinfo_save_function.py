from orbkit.read.high_level import main_read
from orbkit.test.tools import equal
from orbkit.qcinfo import QCinfo
from orbkit import options
import os, inspect

options.quiet = True

tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
folder = os.path.join(tests_home, '../read/outputs_for_testing')
filepath = os.path.join(folder, 'NaCl_molden_files.tar.gz')

qc = main_read(filepath)
filepath = os.path.join(tests_home, 'tmp_save.npz')
qc.save(filepath)

qc_new = QCinfo(filepath)

for i in range(len(qc.ao_spec)):
  equal(qc.ao_spec[i]['coeffs'], qc_new.ao_spec[i]['coeffs'])

for prop in ['coeffs', 'energy', 'occ_num']:
  for i in range(len(qc.mo_spec)):
    equal(qc.mo_spec[i][prop], qc_new.mo_spec[i][prop])

equal(qc.geo_spec, qc_new.geo_spec)

os.remove(filepath)
