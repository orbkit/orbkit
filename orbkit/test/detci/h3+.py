import numpy
import os, inspect

from orbkit import grid, main_read, core
from orbkit.analytical_integrals import get_ao_overlap,get_mo_overlap_matrix
from orbkit.analytical_integrals import get_ao_dipole_matrix
from orbkit.test.tools import equal
from orbkit import options
from orbkit import detci

options.quiet = True

sin = numpy.sin
cos = numpy.cos
sqrt = numpy.sqrt

# Code for the complex superposition state
# Psi_pm = 1/sqrt(2) (\psi_0 + 1/sqrt(2) (psi_1 \pm i psi_2) exp(-i\DeltaE t/hbar))
# with \DeltaE = E_1 - E_0 = E_2 - E_0
pm = +1
get_left = lambda rho,DeltaE,tau: DeltaE/sqrt(2) * (-rho[(0,1)]*sin(tau) +pm* rho[(0,2)]*cos(tau))
get_j = lambda j,tau: -pm*1/2.*j[(1,2)] + 1./sqrt(2.)*j[(0,1)]*sin(tau) -pm* 1./sqrt(2.)*j[(0,2)]*cos(tau)
# where tau=\DeltaE t/hbar

# Set options
numproc = options.numproc
slice_length = 1e2
ci_readthreshold = 0.0

tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
folder = os.path.join(tests_home, 'read/outputs_for_testing')
file_ci = os.path.join(folder, 'h3+_fci_cc-pVTZ.out')
file_molden = file_ci + '.default.molden'

qc = main_read(file_molden, all_mo=True)
qc,ci = detci.ci_read.main_ci_read(qc,file_ci,itype='psi4_detci',
                                   threshold=ci_readthreshold)

grid.min_   = [-2.5,-2.5, 0.0]
grid.max_   = [ 2.5, 2.5, 0.0]
grid.delta_ = [ 0.1, 0.1, 0.1]
grid.grid_init()    

molist = core.rho_compute(qc,
                          calc_mo=True,
                          slice_length=slice_length,
                          drv=[None,'x','y','z','xx','yy','zz'],
                          numproc=numproc)
molistdrv = molist[1:4]                                # \nabla of MOs
molistdrv2 = molist[-3:]                               # \nabla^2 of MOs
molist = molist[0]                                     # MOs

aoom = get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec,
                      drv=[None,'x','y','z'])

dm_aoom = get_ao_dipole_matrix(qc,component=['x','y','z'])

moom = get_mo_overlap_matrix(qc.mo_spec,qc.mo_spec,aoom[0],
                             numproc=numproc)          # <m|n>
omr = numpy.zeros((3,) + moom.shape)                   # <m|r|n>
omv = numpy.zeros((3,) + moom.shape)                   # <m|\nabla|n>

for i in range(3):
  omr[i] = get_mo_overlap_matrix(qc.mo_spec,qc.mo_spec,dm_aoom[i],numproc=numproc)
  omv[i] = get_mo_overlap_matrix(qc.mo_spec,qc.mo_spec,aoom[i+1],numproc=numproc)



rho_01 = []
j_01 = []
nabla_j_01 = []
mu_01 = []
for a in range(len(ci)):
  for b in range(a+1,len(ci)):
    name = ci[a].info['state'] + ' -> ' + ci[b].info['state']

    zero,sing = detci.occ_check.compare(ci[a],ci[b],numproc=numproc)
    
    rho_01.append(detci.ci_core.rho(zero,sing,molist,slice_length=slice_length,numproc=numproc))

    j_01.append(detci.ci_core.jab(zero,sing,molist,molistdrv,
                                  slice_length=slice_length,numproc=numproc))

    nabla_j_01.append(-numpy.sum(detci.ci_core.jab(zero,sing,molist,molistdrv2,
                                  slice_length=slice_length,numproc=numproc),
                                  axis=0)) # Sum over the three components
    
    mu_01.append(detci.ci_core.mu(ci[a],ci[b],qc,zero,sing,omr,omv))


filepath = os.path.join(tests_home, 'refdata_h3+.npz')

#numpy.savez('refdata_h3+.npz', rho_01=rho_01, j_01=j_01, nabla_j_01=nabla_j_01, mu_01=mu_01)
refdata = numpy.load(filepath)
equal(rho_01, refdata['rho_01'])
equal(j_01, refdata['j_01'])
equal(nabla_j_01, refdata['nabla_j_01'])
equal(mu_01, refdata['mu_01'])






