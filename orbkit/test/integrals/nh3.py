import numpy
from orbkit import read, integrals
from orbkit import options
from orbkit.test.tools import equal

import os, inspect
tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))

options.quiet = True

## NH3 STO-3G RHF
qc = read.main_read(os.path.join(tests_home, 'nh3.mold'), all_mo=True)
Nelec = int(qc.get_elec_charge())
ao = integrals.AOIntegrals(qc)

## test basic 1- and 2-electron integrals
S = ao.overlap(asMO=0)
S_mo = ao.overlap(asMO=1)
T = ao.kinetic(asMO=0)
V = ao.Vne(asMO=0)
Vee = ao.Vee(asMO=0)
Vee_mo = ao.Vee(asMO=1)

ref = numpy.load(os.path.join(tests_home, 'nh3.npz'))

equal(S, ref['S'])
equal(S_mo, numpy.eye(ao.Norb))
equal(T, ref['T'])
equal(V, ref['V'])
equal(Vee, ref['Vee'])
equal(Vee_mo, ref['Vee_mo'])

## test Hartree-Fock energy
P = 2*numpy.dot(qc.mo_spec.coeffs[:Nelec//2,:].T, qc.mo_spec.coeffs[:Nelec//2,:])
G = P[None,None,:,:]*Vee - numpy.swapaxes(0.5*P.T[None,:,:,None]*Vee, 1, 3)
G = numpy.sum(G, axis=(2, 3))
Hcore = T + V
F = Hcore + G
EHF = 0.5*numpy.tensordot(P, Hcore+F, axes=2)

Vnn = qc.nuclear_repulsion

# compare
equal(Vnn, 11.73717604)
equal(EHF+Vnn, -55.455419778557)
