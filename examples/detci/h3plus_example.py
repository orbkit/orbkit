'''
This example file shows for a Full CI calculations of a H3+ molecule 
how to use detCI@ORBKIT to compute and visualize several quantities 
related to the electronic continuity equation for 
multi-determinant CI wavefunctions.

This example is based on the first example of the detCI@ORBKIT publication.

It is shown...

  1a. How to read the quantum chemistry output of a PSI4 CI calculation.
      The input and output files for PSI4 may be found in examples/detci/h3plus
     
  1b. How to prepare all subsequent calculations by
      computing the molecular orbitals, the derivatives thereof and
      the analytical expectation values: <a|b>, <a|r|b>, <a|nabla|b> 

  2a. How to compare the determinants of different electronic states in order
      to identify identical determinants and formal single excitations, 
      which are required for evaluation of one-electron operators.
     
  2b. How to compute the transition electron densities, 
      transition electronic flux densities, and transition dipole moments
  
  3a. How to visualize the time independent quantities: 
        - electronic flow (time-derivative of the density)
        - electronic flux density
        - divergence of the electronic flux density
  
  3b. How to compute and visualize the time evolution of complex superposition state:
        - electronic flow (time-derivative of the density)
        - electronic flux density
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
sqrt = numpy.sqrt

# Code for the complex superposition state
# Psi_pm = 1/sqrt(2) (\psi_0 + 1/sqrt(2) (psi_1 \pm i psi_2) exp(-i\DeltaE t/hbar))
# with \DeltaE = E_1 - E_0 = E_2 - E_0
pm = +1
get_left = lambda rho,DeltaE,tau: DeltaE/sqrt(2) * (-rho[(0,1)]*sin(tau) +pm* rho[(0,2)]*cos(tau))
get_j = lambda j,tau: -pm*1/2.*j[(1,2)] + 1./sqrt(2.)*j[(0,1)]*sin(tau) -pm* 1./sqrt(2.)*j[(0,2)]*cos(tau)
# where tau=\DeltaE t/hbar

# Define Constants
eV = 27.211384
Debye = 2.5415800000001716
Ampere = 0.006623619

# Set options
numproc = 4
slice_length = 1e2
ci_readthreshold = 0.0

print('''
==========================================
Reading the Quantum Chemistry Output Files
==========================================
''')
# 1a.
fid_psi4 = 'psi4_data/h3+_fci_cc-pVTZ.out'
fid_molden = fid_psi4 + '.default.molden'
qc = read.main_read(fid_molden,all_mo=True,interactive=False)
qc,ci = detci.ci_read.main_ci_read(qc,fid_psi4,itype='psi4_detci',
                                   threshold=ci_readthreshold)

print('='*80)
print('Preparing All Subsequent Calculations')
print('='*80+'\n')
# 1b.
print('Setting up the grid...')
grid.min_   = [-2.5,-2.5, 0.0]
grid.max_   = [ 2.5, 2.5, 0.0]
grid.delta_ = [ 0.1, 0.1, 0.1]
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
# Create a dictionary for the results
rho_01 = {}
j_01 = {}
nabla_j_01 = {}
mu_01 = {}
for a in range(len(ci)):
  for b in range(a+1,len(ci)):
    name = ci[a].info['state'] + ' -> ' + ci[b].info['state']
    print('\n'+'='*len(name)+'\n'+name+'\n'+'='*len(name))
    
    # Compare the occupation patterns
    # 2a.
    zero,sing = detci.occ_check.compare(ci[a],ci[b],numproc=numproc)
    
    # 2b.
    print('\nComputing the transition density...')
    rho_01[(a,b)] = detci.ci_core.rho(zero,sing,molist,slice_length=slice_length,numproc=numproc)
    print('\nComputing the transition flux density...')
    # The transition density between electronic eigenstates is purely imaginary.
    j_01[(a,b)] = detci.ci_core.jab(zero,sing,molist,molistdrv,
                                  slice_length=slice_length,numproc=numproc)
    print('and its divergence...')
    nabla_j_01[(a,b)] = -numpy.sum(detci.ci_core.jab(zero,sing,molist,molistdrv2,
                                  slice_length=slice_length,numproc=numproc),
                                  axis=0) # Sum over the three components
    
    mu_01[(a,b)] = detci.ci_core.mu(ci[a],ci[b],qc,zero,sing,omr,omv)
    print('\nAnalytical transition dipole moments:')
    print('\t%+.04f D\t%+.04f D\t%+.04f D (in length gauge)' % tuple(mu_01[(a,b)][0]*Debye))
    print('\t%+.04f D\t%+.04f D\t%+.04f D (from velocity gauge)' % tuple(mu_01[(a,b)][-1]*Debye))  
    print('\t%+.04f i mA\t%+.04f i mA\t%+.04f i mA (in velocity gauge)' % tuple(mu_01[(a,b)][1]*1e3*Ampere))  

print('''
==========================================
          Plotting the Results
==========================================
''')
yv,xv  = numpy.meshgrid(grid.y,grid.x)
levels = numpy.linspace(-0.2,0.2,num=50)

print('Showing the time-independent convergence...')
# 3a.
fig, axs = plt.subplots(nrows=len(ci), ncols=3,sharey=True,sharex=True)
fig.subplots_adjust(hspace=0.3, wspace=0)
axs[0,0].set_aspect('equal',share=True)

# Plot axis labels
for i in range(len(axs)):
  axs[i,0].set_ylabel(r'$y\ [a_0]$')
  for ax in axs[i]:
    # Plot the molecule
    ax.plot(*qc.geo_spec[[0,1,2,0],:2].T,color='grey',lw=0.5)
for ax in axs[-1]:
  ax.set_xlabel(r'$x\ [a_0]$')

a = 0 #: Considering all transitions from the ground state
for i,b in enumerate(range(1,len(ci))):
  name = ci[a].info['state'] + ' -> ' + ci[b].info['state']
  axs[i,1].set_title('Convergence of ' + name,fontsize=10)
  deltaE_01 = ci[b].info['energy'] - ci[a].info['energy']
  
  for j,label in enumerate([r'$\dot{\rho} = -\frac{\Delta E_{ge}}{\hbar} \rho_{ge}$',
                            r'$\vec{j}_{ge}$',
                            r'$-\nabla \cdot \vec{j}_{ge}$']):
    axs[i,j].text(0.05,0.95,label,va='top',ha='left',transform=ax.transAxes,fontsize=8)
  
  # Plot the left and the right side of the continuity equation
  contour = axs[i,0].contourf(xv,yv,-deltaE_01*rho_01[(a,b)].squeeze(),
                              levels=levels, cmap='RdBu', extend="both")
  contour = axs[i,2].contourf(xv,yv,nabla_j_01[(a,b)].squeeze(),
                              levels=levels, cmap='RdBu', extend="both")
  
  # Plot the corresponding electronic flux density
  Q = axs[i,1].quiver(xv[::6,::6],yv[::6,::6],
                      j_01[(a,b)][0,::6,::6].squeeze(),
                      j_01[(a,b)][1,::6,::6].squeeze(),
                      scale=0.7,units='width')
  

print('Showing the time evolution for the super position state')
# 3a.
names = [c.info['state'] for c in ci]
string = '\Psi = 1/\sqrt{2} [\psi_{%s} + 1/\sqrt{2} (\psi_{%s} + i \psi_{%s}) e^{-i(\Delta E t/hbar + phase)}]' % tuple(names)
print('\t'+string)
print('with Delta E = %.2f eV' % (deltaE_01*eV))
axs[-1,1].set_title(r'$%s$' % string,fontsize=10)

for i,t in enumerate([0,numpy.pi/2.,numpy.pi]):
  axs[-1,i].text(0.05,0.95,r'$t=%.1f\tau$' % (t/(2*numpy.pi)),
                 va='top',ha='left',transform=axs[-1,i].transAxes,fontsize=8)
  
  # Plot the left and the right side of the continuity equation
  contour = axs[-1,i].contourf(xv,yv,
                               get_left(rho_01,deltaE_01,t).squeeze(),
                               levels=levels, cmap='RdBu', extend="both")
  # Plot the corresponding electronic flux density
  # Taking only each 6th data point
  jx,jy,jz = get_j(j_01,t).squeeze()
  Q = axs[-1,i].quiver(xv[::6,::6],yv[::6,::6],jx[::6,::6],jy[::6,::6],scale=0.7,units='width')

plt.show()
