from orbkit.read.high_level import main_read
from orbkit.test.tools import equal
from orbkit.qcinfo import QCinfo
from orbkit import options
import os, inspect

options.quiet = True

tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
folder = os.path.join(tests_home, '../read/outputs_for_testing')
filepath = os.path.join(folder, 'NaCl_molden_files.tar.gz')

qc_old = main_read(filepath)
filepath = os.path.join(tests_home, 'tmp_save.npz')
qc_old.save(filepath)

qc_new = QCinfo(filepath)
equal(qc_old, qc_new)

os.remove(filepath)
