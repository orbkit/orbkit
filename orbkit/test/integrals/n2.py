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
Nelec = 14

# MO overlap
S = ao.overlap(asMO=1)
equal(S, numpy.eye(ao.Norb))

# calculate Hartree-Fock energy from AOs
Hcore = ao.kinetic(asMO=0) + ao.Vne(asMO=0)
Vee = ao.Vee(asMO=0)
Vnn = qc.nuclear_repulsion

P = numpy.zeros(Hcore.shape)
for i in range(ao.Norb):
  for j in range(i,ao.Norb):
    for k in range(Nelec//2):
      P[i,j] += 2*qc.mo_spec.coeffs[k,i]*qc.mo_spec.coeffs[k,j]
    P[j,i] = P[i,j]

G = numpy.zeros(Hcore.shape)
for i in range(ao.Norb):
  for j in range(i,ao.Norb):
    for k in range(ao.Norb):
      for l in range(ao.Norb):
        G[i,j] += P[k,l]*(Vee[i,j,k,l]-0.5*Vee[i,l,k,j])
    G[j,i] = G[i,j]

F = Hcore + G

E = 0
for i in range(ao.Norb):
  E += 0.5*P[i,i]*(Hcore[i,i]+F[i,i])
  for j in range(i+1,ao.Norb):
      E += P[i,j]*(Hcore[i,j]+F[i,j])

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
V = ao.Vee(asMO=1)
Vs = ao.Vee(asMO=1, max_dims=3)
equal(V, Vs)

V = ao.Vee(asMO=1, MOrange=range(Nelec//2))
Vs = ao.Vee(asMO=1, MOrange=range(Nelec//2), max_dims=3)
equal(V, Vs)
