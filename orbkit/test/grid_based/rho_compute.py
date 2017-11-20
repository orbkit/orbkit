import numpy
import os, inspect

from orbkit import read
from orbkit.core import rho_compute
from orbkit.test.tools import equal
from orbkit import grid
from orbkit import options

options.quiet = True

tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
folder = os.path.join(tests_home, '../outputs_for_testing/molpro')
filepath = os.path.join(folder, 'h2o_rhf_sph.molden')
qc = read.main_read(filepath, all_mo=True)

grid.adjust_to_geo(qc, extend=2.0, step=1)
grid.grid_init(is_vector=False,force=True)

drv = [None,'x','y','z','xx','xy','xz','yy','yz','zz']
data = []
for i in range(2):
  if i: grid.grid2vector()
  data.append([
      rho_compute(qc,slice_length=0),
      rho_compute(qc,numproc=options.numproc),
      rho_compute(qc,laplacian=True,slice_length=0)[-1],
      rho_compute(qc,laplacian=True,numproc=options.numproc)[-1],
      rho_compute(qc,calc_mo=True,drv=drv,slice_length=0),
      rho_compute(qc,calc_mo=True,drv=drv,numproc=options.numproc)
      ])

data[1] = [grid.mv2g(d=i) for i in data[1]]

for i in range(len(data[0])):
  equal(data[0][i],data[1][i])

filepath = os.path.join(tests_home, 'refdata_rho_compute.npz')
#numpy.savez(filepath, data=data[0])
refdata = numpy.load(filepath)
in_dic = {0: 'zero', 1: 'one', 2: 'two', 3: 'three', 4: 'four'}
for i in range(5):
  equal(data[0][i], refdata[in_dic[i]])
