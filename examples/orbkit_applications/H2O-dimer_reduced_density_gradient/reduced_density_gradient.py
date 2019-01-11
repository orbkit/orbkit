'''
This script shows how to implement a new grid based quantity.

Here, we implement the reduced density gradient as proposed by 
Johnson et al., J. Am. Chem. Soc. 132, 6498 (2010).
'''

from orbkit import read, grid, display, core, output

'''
Identify the non-covalent interaction regions 
using the reduced density gradient
'''

# open molden file and read parameters
qc = read.main_read('h2o_dimer.molden',itype='molden')

# Adjust the grid parameters to the molecular structure
grid.adjust_to_geo(qc,extend=1.0,step=0.1)
# initialize grid
grid.grid_init()
# print grid information
display.display(grid.get_grid())

# Run orbkit
rho,delta_rho = core.rho_compute(qc,drv=['x','y','z'],numproc=4)

# Compute the reduced density gradient 
from numpy import sqrt,pi
s = (1./(2*(3*pi**2)**(1/3.)) * sqrt((delta_rho**2).sum(axis=0))/(rho**(4/3.)))

# Apply a density cutoff, i.e., consider only values for rho < 0.05 a.u.
s[rho>0.05] = 1e3

# Write a cube file and vmd script for the reduced density gradient 
output.main_output(s,
                   qc,                          # atomic information
                   outputname='reduced_drho',
                   otype=['cb','vmd'],
                   iso=(-0.5,0.5)               # isocontour value of s = 0.5 a.u.
                   )

'''
Classify the interaction types
using the second derivatives of the density 
("the three eigenvalues of the electronic Hessian matrix")
'''
# Compute the second derivatives of the density
rho,delta2_rho = core.rho_compute(qc,drv=['xx','yy','zz'],numproc=4)

# Sort the three components of the derivatives
# The sign of the second value determines, if the interaction is bonding (negative)
# non-bonding (positive)
delta2_rho.sort(axis=0)

# Plot of the reduced density gradient as a function of the electron density 
# multiplied by the sign of the second largest component
import pylab as plt
from numpy import sign
plt.plot((sign(delta2_rho[1])*rho).reshape((-1,)),s.reshape((-1),),'.')
plt.ylabel('Reduced Density Gradient')
plt.xlabel(r'Electron Density $\times$ Sign of Second Hessian Eigenvalue')
plt.ylim(0,2)
plt.xlim(-0.05,0.05)
plt.show()
