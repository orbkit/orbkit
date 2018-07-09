import numpy
import os, inspect

from orbkit import read
from orbkit.core import rho_compute
from orbkit.test.tools import equal
from orbkit import grid
from orbkit import options
from orbkit import units

options.quiet = True
dx = 4.970736
grid.x = numpy.arange(3)*dx-4.970736
grid.y = numpy.arange(3)*dx-4.970736
grid.z = numpy.arange(3)*dx-4.732975
grid.is_initialized = True
grid.is_regular = True
grid.is_vector = False

tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
folder = os.path.join(tests_home, '../outputs_for_testing/gaussian')
filepath = os.path.join(folder, 'h2o_rhf_sph.fchk')
qc = read.main_read(filepath, all_mo=False)

rho,delta_rho = rho_compute(qc,drv='xyz')
rho2,delta_rho2,laplacian = rho_compute(qc,laplacian=True)
equal(rho,rho2)
#equal(delta_rho,delta_rho2) #BUG

refdata = numpy.genfromtxt(os.path.join(folder, 'h2o_rhf_sph_grad.cube'), 
                           skip_header=9).reshape((-1,))
refrho = numpy.zeros_like(rho)
refdrho = numpy.zeros((3,) + rho.shape)
c = 0
for i in range(len(grid.x)):
  for j in range(len(grid.y)):
    for  k in range(len(grid.z)):
      refrho[i,j,k] = refdata[c]
      c += 1
      for l in range(3):
        refdrho[l,i,j,k] = refdata[c]
        c += 1

reflaplace = numpy.genfromtxt(os.path.join(folder, 'h2o_rhf_sph_laplace.cube'), 
                           skip_header=9).reshape(laplacian.shape)

equal(rho,refrho)
equal(delta_rho,refdrho, tol=1e-03)
equal(laplacian,reflaplace)

qc.mo_spec = qc.mo_spec[1:] # Gaussian ignores core density for computing Total SCF Density

data = rho_compute(qc)

refdata = numpy.genfromtxt(os.path.join(folder, 'h2o_rhf_sph.cube'), 
                           skip_header=9).reshape(data.shape)
equal(data,refdata)

qc.mo_spec = qc.mo_spec['homo']
data = rho_compute(qc,calc_mo=True)

refdata = numpy.genfromtxt(os.path.join(folder, 'h2o_rhf_sph_mo.cube'), 
                           skip_header=10).reshape(data.shape)
equal(data,refdata)