import numpy
from orbkit import read, options
from orbkit import libcint_interface as integrals
from orbkit.test.tools import equal

import os, inspect
tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))

options.quiet = True

## N2 cc-pVDZ RHF
folder = os.path.join(tests_home, '../outputs_for_testing/molpro')
qc = read.main_read(os.path.join(folder, 'n2.mold'), all_mo=True)
ao = integrals.AOIntegrals(qc)
Nelec = int(abs(qc.get_charge(nuclear=False)))

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
Hcore_mo = integrals.ao2mo(Hcore, qc.mo_spec.coeffs, MOrange=range(Nelec//2))
Vee_mo = integrals.ao2mo(Vee, qc.mo_spec.coeffs, MOrange=range(Nelec//2))

E = 2*numpy.trace(Hcore_mo[:Nelec//2,:Nelec//2])
J = numpy.trace(numpy.trace(Vee_mo[:Nelec//2,:Nelec//2,:Nelec//2,:Nelec//2], axis1=2, axis2=3))
K = numpy.trace(numpy.trace(Vee_mo[:Nelec//2,:Nelec//2,:Nelec//2,:Nelec//2], axis1=1, axis2=2))
E += 2*J - K

equal(E+Vnn, -108.92022806)

# test AO sclicing
ao.add_MO_block_1e(MOrange=range(Nelec//2))
Hcore = ao.Hcore(asMO=1)
Hs = ao.Hcore(asMO=1, max_dims=3)
equal(Hcore, Hs)

ao.add_MO_block_2e(MOrange=range(Nelec//2))
V = ao.Vee(asMO=1)
Vs = ao.Vee(asMO=1, max_dims=3)
equal(V, Vs)

# test AO ranges
ao.clear_all_blocks()
S = ao.overlap(asMO=0)
rangei = [5,6,7] + [12,13]
rangej = [4,5]   + [10,11,12,13,14,15]
ao.add_AO_block_1e(AOrangei=rangei, AOrangej=rangej)
Sr = ao.overlap(asMO=0)
equal(S[numpy.ix_(rangei,rangej)], Sr)

# test MOranges
c = ao.Norb//2
ao.clear_all_blocks()
ao.add_MO_block_1e(MOrangei=range(c), MOrangej=range(c, ao.Norb))
ao.add_MO_block_1e(MOrangei=range(c, ao.Norb), MOrangej=range(c, ao.Norb))
S = ao.overlap(asMO=1)
equal(S[0], numpy.zeros(S[0].shape))
equal(S[1], numpy.eye(S[1].shape[0]))

# test FCIDUMP generator without symmetry
fcidump = integrals.generate_fcidump(qc, '')
fcidump_ref = integrals.load_fcidump(os.path.join(tests_home, 'FCIDUMP_N2'))

equal(fcidump.nuclear_repulsion, fcidump_ref.nuclear_repulsion)
equal(fcidump.get_H(full=True), fcidump_ref.get_H(full=True))
equal(fcidump.get_G(full=True), fcidump_ref.get_G(full=True))

# test FCIDUMP generator with symmetry
qc = read.main_read(os.path.join(folder, 'n2_c2v.mold'), all_mo=True)
fcidump = integrals.generate_fcidump(qc, '', sym='c2v')
fcidump_ref = integrals.load_fcidump(os.path.join(tests_home, 'FCIDUMP_N2_c2v'))

equal(fcidump.nuclear_repulsion, fcidump_ref.nuclear_repulsion)
equal(fcidump.get_H(full=True), fcidump_ref.get_H(full=True))
equal(fcidump.get_G(full=True), fcidump_ref.get_G(full=True))
