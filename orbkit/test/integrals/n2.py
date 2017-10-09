import numpy
from orbkit import read, options
import orbkit.integrals
from orbkit.test.tools import equal

import os, inspect
tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))

options.quiet = True

## N2 cc-pVDZ RHF
qc = read.main_read(os.path.join(tests_home, 'n2.mold'), all_mo=True)
ao = orbkit.integrals.AOIntegrals(qc, cartesian=True)
Nelec = int(qc.get_elec_charge())

# MO overlap
S = ao.overlap(asMO=1)
equal(S, numpy.eye(ao.Norb))

# calculate Hartree-Fock energy from AOs
Hcore = ao.Hcore(asMO=0)
Vee = ao.Vee(asMO=0)
Vnn = qc.nuclear_repulsion

P = 2*numpy.dot(qc.mo_spec.coeffs[:Nelec//2,:].T, qc.mo_spec.coeffs[:Nelec//2,:])
G = P[None,None,:,:]*Vee - numpy.swapaxes(0.5*P.T[None,:,:,None]*Vee, 1, 3)
G = numpy.sum(G, axis=(2, 3))
F = Hcore + G
E = 0.5*numpy.tensordot(P, Hcore+F, axes=2)

equal(E+Vnn, -108.92022806)

# calculate Hartree-Fock energy from MOs
Hcore_mo = orbkit.integrals.ao2mo(Hcore, qc.mo_spec.coeffs, MOrange=range(Nelec//2))
Vee_mo = orbkit.integrals.ao2mo(Vee, qc.mo_spec.coeffs, MOrange=range(Nelec//2))

E = 2*numpy.trace(Hcore_mo[:Nelec//2,:Nelec//2])
J = numpy.trace(numpy.trace(Vee_mo[:Nelec//2,:Nelec//2,:Nelec//2,:Nelec//2], axis1=2, axis2=3))
K = numpy.trace(numpy.trace(Vee_mo[:Nelec//2,:Nelec//2,:Nelec//2,:Nelec//2], axis1=1, axis2=2))
E += 2*J - K

equal(E+Vnn, -108.92022806)

# test AO sclicing
ao.add_MO_block_1e(MOrangei=range(Nelec//2), MOrangej=range(Nelec//2))
Hcore = ao.Hcore(asMO=1)
Hs = ao.Hcore(asMO=1, max_dims=3)
equal(Hcore, Hs)

ao.add_MO_block_2e(MOrangei=range(Nelec//2), MOrangej=range(Nelec//2), MOrangek=range(Nelec//2), MOrangel=range(Nelec//2))
V = ao.Vee(asMO=1)
Vs = ao.Vee(asMO=1, max_dims=3)
equal(V, Vs)

# test MOranges
c = ao.Norb//2
ao.clear_all_blocks()
ao.add_MO_block_1e(MOrangei=range(c), MOrangej=range(c, ao.Norb))
ao.add_MO_block_1e(MOrangei=range(c, ao.Norb), MOrangej=range(c, ao.Norb))
S = ao.overlap(asMO=1)
equal(S[0], numpy.zeros(S[0].shape))
equal(S[1], numpy.eye(S[1].shape[0]))

# test FCIDUMP generator without symmetry
fcidump = orbkit.integrals.generate_fcidump(qc, '')
fcidump_ref = orbkit.integrals.load_fcidump(os.path.join(tests_home, 'FCIDUMP_N2'))

equal(fcidump.nuclear_repulsion, fcidump_ref.nuclear_repulsion)
equal(fcidump.get_H(full=True), fcidump_ref.get_H(full=True))
equal(fcidump.get_G(full=True), fcidump_ref.get_G(full=True))

# test FCIDUMP generator with symmetry
qc = read.main_read(os.path.join(tests_home, 'n2_c2v.mold'), all_mo=True)
fcidump = orbkit.integrals.generate_fcidump(qc, '', sym='c2v')
fcidump_ref = orbkit.integrals.load_fcidump(os.path.join(tests_home, 'FCIDUMP_N2_c2v'))

equal(fcidump.nuclear_repulsion, fcidump_ref.nuclear_repulsion)
equal(fcidump.get_H(full=True), fcidump_ref.get_H(full=True))
equal(fcidump.get_G(full=True), fcidump_ref.get_G(full=True))
