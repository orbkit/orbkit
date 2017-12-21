import numpy
from orbkit import read
from orbkit import options
from orbkit import output
from orbkit.test.tools import equal

import os, inspect
tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))

qc = read.main_read(os.path.join(tests_home, 'H2O.molden'), all_mo=True)
output.molden_writer(qc, 'h2o_new.molden')
