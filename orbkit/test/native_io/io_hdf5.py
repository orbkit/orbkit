import os, shutil, inspect, tempfile, numpy

from orbkit.read.high_level import main_read
from orbkit.output.high_level import main_output
from orbkit.grid import set_grid,get_shape
from orbkit.test.tools import equal
from orbkit.qcinfo import QCinfo
from orbkit import options

options.quiet = True

tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
folder = os.path.join(tests_home, '../outputs_for_testing')
filepath = os.path.join(folder, 'NaCl_molden_files.tar.gz')

qc_old = main_read(filepath)

test_dir = tempfile.mkdtemp()
testfile = os.path.join(test_dir, 'tmp.hdf5')

main_output(qc_old, outputname=testfile, otype='native', ftype='hdf5')

qc_new = main_read(testfile)

equal(qc_old, qc_new)

# Test restart from standard HDF5 data output
set_grid(0,0,0,is_vector=False)
main_output(numpy.zeros(get_shape()),qc=qc_old,outputname=testfile)

qc_new = main_read(filepath)
equal(qc_old, qc_new)

shutil.rmtree(test_dir)
