import numpy
import os, inspect

from orbkit import read
from orbkit.analytical_integrals import get_dipole_moment
from orbkit.test.tools import equal

tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
folder = os.path.join(tests_home, '../read/outputs_for_testing')
filepath = os.path.join(folder, 'h2o_rhf_sph.molden')
qc = read.main_read(filepath, all_mo=True)
dip = get_dipole_moment(qc,component=['x','y','z'])

equal(dip, [0.00000000e+00,  -1.01130147e-16,   8.17259184e-01])

#Slightly move one atom and calculate dipoles again
qc.geo_spec[1] += numpy.array([1,1,0])
dip = get_dipole_moment(qc,component=['x','y','z'])

equal(dip, [-0.21150947, -0.66894653,  1.09887067])
