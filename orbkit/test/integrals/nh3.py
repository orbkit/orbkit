import numpy
from orbkit import read, integrals
from orbkit import options
from orbkit.test.tools import equal

import os, inspect
tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))

options.quiet = True

qc = read.main_read(os.path.join(tests_home, 'nh3.mold'), all_mo=True)
ao = integrals.AOIntegrals(qc)

S = ao.overlap(asMO=0)
S_mo = ao.overlap(asMO=1)
T = ao.kinetic(asMO=0)
V = ao.Vne(asMO=0)
Vee = ao.Vee(asMO=0)
Vee_mo = ao.Vee(asMO=1)

ref = numpy.load(os.path.join(tests_home, 'nh3.npz'))

equal(S, ref['S'])
equal(S_mo, ref['S_mo'])
equal(T, ref['T'])
equal(V, ref['V'])
equal(Vee, ref['Vee'])
equal(Vee_mo, ref['Vee_mo'])
