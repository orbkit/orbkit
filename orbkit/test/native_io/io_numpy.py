from orbkit.read.high_level import main_read
from orbkit.output.high_level import main_output
from orbkit.grid import set_grid,get_shape
from orbkit.test.tools import equal
from orbkit.qcinfo import QCinfo
from orbkit import options
import os, inspect, tempfile, numpy

options.quiet = True

tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
folder = os.path.join(tests_home, '../outputs_for_testing')
filepath = os.path.join(folder, 'NaCl_molden_files.tar.gz')

qc_old = main_read(filepath)

tests_home = tempfile.gettempdir()
filepath = os.path.join(tests_home, 'tmp.npz')

main_output(qc_old, outputname=filepath, otype='native', ftype='numpy')

qc_new = main_read(filepath)
equal(qc_old, qc_new)

os.remove(filepath)

# Test restart from standard npz data output
set_grid(0,0,0,is_vector=False)
main_output(numpy.zeros(get_shape()),qc=qc_old,outputname=filepath)

qc_new = main_read(filepath)
equal(qc_old, qc_new)

os.remove(filepath)