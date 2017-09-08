'''
This example file shows for a CIS calculations of a LiH molecule 
how to use detCI@ORBKIT to compute and visualize several quantities 
related to the electronic continuity equation for 
multi-determinant CI wavefunctions.

This example is based on the second example of the detCI@ORBKIT publication.

It is shown...

  1a. How to read the quantum chemistry output of a PSI4 CI calculation.
      The input and output files for PSI4 may be found in examples/detci/h3plus
     
  1b. How to prepare all subsequent calculations by
      computing the molecular orbitals, the derivatives thereof and
      the analytical expectation values: <a|b>, <a|r|b>, <a|nabla|b> 

  2a. How to compare the determinants of different electronic states in order
      to identify identical determinants and formal single excitations, 
      which are required for evaluation of one-electron operators.
     
  2b. How to compute compute the transition electron densities, 
      transition electronic flux densities, and transition dipole moments
  
  3b. How to compute and visualize the time evolution the electronic flux density
      as streamline plot and the difference density for a superposition state.

'''
import numpy
import matplotlib.pyplot as plt

# Import the required orbkit modules
from orbkit import grid, read, core
from orbkit.analytical_integrals import get_ao_overlap,get_mo_overlap_matrix
from orbkit.analytical_integrals import get_ao_dipole_matrix

from orbkit import detci

sin = numpy.sin
cos = numpy.cos

# Code for the superposition state
# Psi = 1/sqrt(2) (psi_0 + psi_e  exp(-i(\DeltaE t/hbar + phase))
# with \DeltaE = E_e - E_0 
phase = numpy.pi
get_rho = lambda rho_0,rho_1,rho_01,tau: 1/2. * (rho_0+rho_1) + rho_01*cos(tau+phase)
get_yield = lambda rho_01,tau: rho_01*(cos(tau+phase) - cos(phase))
get_left = lambda rho_01,DeltaE,tau: -DeltaE * rho_01*sin(tau+phase) 
get_j = lambda j_01,tau: j_01*sin(tau+phase)
# where tau=\DeltaE t/hbar

# Define Some Constants
eV = 27.211384
Debye = 2.5415800000001716
Ampere = 0.006623619

# Set some options
numproc = 4
slice_length = 1e2
ci_readthreshold = 0.0

print('''
==========================================
Reading the Quantum Chemistry Output Files
==========================================
''')
# 1a.
fid_psi4 = 'psi4_data/lih_cis_aug-cc-pVTZ.out'
fid_molden = fid_psi4 + '.default.molden'
qc = read.main_read(fid_molden,all_mo=True,interactive=False)
qc,ci = detci.ci_read.main_ci_read(qc,fid_psi4,itype='psi4_detci',
                                threshold=ci_readthreshold)

print('='*80)
print('Preparing All Subsequent Calculations')
print('='*80+'\n')
# 1b.
print('Setting up the grid...')
grid.min_   = [-3.9, 0.0,-5.9]
grid.max_   = [ 3.9, 0.0, 9.9]
grid.delta_ = [ 0.2, 0.0, 0.2]
grid.grid_init()    
print(grid.get_grid())

print('Computing the Molecular Orbitals and the Derivatives Thereof...\n')
molist = core.rho_compute(qc,
                          calc_mo=True,
                          slice_length=slice_length,
                          drv=[None,'x','y','z','xx','yy','zz'],
                          numproc=numproc)
molistdrv = molist[1:4]                                # \nabla of MOs
molistdrv2 = molist[-3:]                               # \nabla^2 of MOs
molist = molist[0]                                     # MOs

print('\nComputing the Analytical Overlaps of the Molecular Orbitals...\n')
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

print('''
==========================================
  Starting the detCI@ORBKIT Computations
==========================================
''')
a = 0 # Ground State 
b = 1 # First Excited State

'''
Uncomment the following lines to compute the electron density, electron number, 
and permanent dipole moments of the ground and first excited state
'''
# zero,sing = detci.occ_check.compare(ci[a],ci[a],numproc=numproc)
# en_0 = detci.ci_core.enum(zero,sing,moom)
# mu_0 = detci.ci_core.mu(ci[a],ci[a],qc,zero,sing,omr,omv)
# rho_0 = detci.ci_core.rho(zero,sing,molist,slice_length=slice_length,numproc=numproc)
# print('\tNorm of state %s: %f' % (ci[a].info['state'],en_0))
# zero,sing = detci.occ_check.compare(ci[b],ci[b],numproc=numproc)
# en_1 = detci.ci_core.enum(zero,sing,moom)
# mu_0 = detci.ci_core.mu(ci[b],ci[b],qc,zero,sing,omr,omv)
# rho_1 = detci.ci_core.rho(zero,sing,molist,slice_length=slice_length,numproc=numproc)
# print('\tNorm of state %s: %f' % (ci[b].info['state'],en_1))

name = ci[a].info['state'] + ' -> ' + ci[b].info['state']
print('\n'+'='*len(name)+'\n'+name+'\n'+'='*len(name))

# Compare the occupation patterns
# 2a.
zero,sing = detci.occ_check.compare(ci[a],ci[b],numproc=numproc)

# 2b.
print('\nComputing the transition density...')
rho_01 = detci.ci_core.rho(zero,sing,molist,slice_length=slice_length,numproc=numproc)
print('\nComputing the transition flux density...')
# detCI@ORBKIT returns the imaginary part of the flux density,
# since the real part vanishes for the real-valued electronic eigenstates.
j_01 = detci.ci_core.jab(zero,sing,molist,molistdrv,
                         slice_length=slice_length,numproc=numproc)
print('and its divergence...')
nabla_j_01 = -numpy.sum(detci.ci_core.jab(zero,sing,molist,molistdrv2,
                        slice_length=slice_length,numproc=numproc),
                        axis=0) # Sum over the three components

mu_01 = detci.ci_core.mu(ci[a],ci[b],qc,zero,sing,omr,omv)
print('\nAnalytical transition dipole moments:')
print('\t%+.04f D\t%+.04f D\t%+.04f D (in length gauge)' % tuple(mu_01[0]*Debye))
print('\t%+.04f D\t%+.04f D\t%+.04f D (from velocity gauge)' % tuple(mu_01[-1]*Debye))  
print('\t%+.04f i mA\t%+.04f i mA\t%+.04f i mA (in velocity gauge)' % tuple(mu_01[1]*1e3*Ampere))  

print('''
==========================================
          Plotting the Results
==========================================
''')
zv,xv = numpy.meshgrid(grid.z,grid.x)
levels = numpy.linspace(-.08,.08,num=60)
  
print('Showing the time evolution for the super position state')
# 3a.
names = [c.info['state'] for c in ci]
string = '\Psi = 1/\sqrt{2} (\psi_{%s} + \psi_{%s}  e^{-i(\Delta E t/\hbar + phase)})' % tuple(names[:2])
print('\t'+string)

fig, axs = plt.subplots(nrows=1, ncols=3,sharey=True,sharex=True)
fig.subplots_adjust(hspace=0.3, wspace=0)
axs[1].set_title(r'$%s$' % string,fontsize=10)
axs[0].set_ylabel(r'$y\ [a_0]$')
for ax in axs:
  ax.axis('equal')
  ax.set_xlabel(r'$x\ [a_0]$')

for i,t in enumerate([2*numpy.pi/4.,2*numpy.pi/2.,2*3*numpy.pi/4]):
  axs[i].text(0.05,0.95,r'$t=%.2f\tau$' % (t/(2*numpy.pi)),
                 va='top',ha='left',transform=axs[i].transAxes,fontsize=8)
  
  # Plot the left and the right side of the continuity equation
  contour = axs[i].contourf(xv,zv,
                            get_yield(rho_01,t).squeeze(),
                            levels=levels, cmap='RdGy', extend="both")
  # Plot the corresponding electronic flux density
  # Taking only each 6th data point
  jx,jy,jz = get_j(j_01,t).squeeze()
  speed = numpy.sqrt(jx**2 + jy**2 + jz**2)
  if speed.sum() != 0. and speed.max() > 1e-2:
    lw = 2e2*speed
    speed[speed>speed.max()/20.] = speed.max()/20.
    Qs = axs[i].streamplot(xv.T,zv.T,jx.T,jz.T,color=speed.T, density=1, cmap='Blues', 
                          arrowstyle='fancy')

plt.show()
