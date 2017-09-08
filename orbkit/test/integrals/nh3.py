import numpy
from scipy.linalg import norm
from orbkit import read, integrals
from orbkit import options
from orbkit.test.tools import equal

import os, inspect
tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))

options.quiet = True

## NH3 STO-3G RHF
qc = read.main_read(os.path.join(tests_home, 'nh3.mold'), all_mo=True)
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
equal(S_mo, ref['S_mo'])
equal(T, ref['T'])
equal(V, ref['V'])
equal(Vee, ref['Vee'])
equal(Vee_mo, ref['Vee_mo'])

## test Hartree-Fock energy
occupation = [2, 2, 2, 2, 2, 0, 0, 0]

# density matrix
P = numpy.zeros(S.shape)
for i in range(ao.Norb):
    for j in range(i,ao.Norb):
        for k in range(ao.Norb):
            P[i,j] += occupation[k]*qc.mo_spec.coeff[k,i]*qc.mo_spec.coeff[k,j]
        P[j,i] = P[i,j]

# build Fock operator
G = numpy.zeros(S.shape)
for i in range(ao.Norb):
    for j in range(i,ao.Norb):
        for k in range(ao.Norb):
            for l in range(ao.Norb):
                G[i,j] += P[k,l]*(Vee[i,j,k,l]-0.5*Vee[i,l,k,j])
        G[j,i] = G[i,j]

Hcore = T + V
F = Hcore + G

# calculate HF energy
EHF = 0
for i in range(ao.Norb):
  # diagonal elements
  EHF += 0.5*P[i,i]*(Hcore[i,i]+F[i,i])
  # offdiagonal elements
  for j in range(i+1,ao.Norb):
      EHF += P[i,j]*(Hcore[i,j]+F[i,j])

# calculate nuclear repulsion
Vnn = 0
Natoms = qc.geo_spec.shape[0]
for a in range(Natoms):
  Za = float(qc.geo_info[a,2])
  Ra = qc.geo_spec[a,:].astype(float)
  for b in range(a+1, Natoms):
        Zb = float(qc.geo_info[b,2])
        Rb = qc.geo_spec[b,:].astype(float)
        Vnn += Za*Zb / norm(Ra-Rb)

# compare
equal(Vnn, 11.73717604)
equal(EHF+Vnn, -55.455419778557)
