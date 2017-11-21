import numpy
import os, inspect

from orbkit import read
from orbkit.analytical_integrals import get_dipole_moment
from orbkit.test.tools import equal
from orbkit import options

'''Reference dipole moments for 'h2o_rhf_cart' (2015-09-30):

\tMolpro 2012 (.out):        0.00000000     0.00000000     0.81739430
\tGaussian 09 (.fchk):     -4.58602321E-17 -3.46944695E-18 -8.17393478E-01
\tTurbomole 6.5 (dscf.log): 0.000000     0.000000     0.817391
'''

ref_dip = {'molpro': [ 0.00000000,        0.00000000,      0.81739430],
           'gaussian': [-4.58602321E-17, -3.46944695E-18, -8.17393478E-01],
           'turbomole': [0.000000,        0.000000,        0.817391],
           }

options.quiet = True
tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
output_folder = os.path.join(tests_home, '../outputs_for_testing')

tests = ['h2o_rhf_cart','h2o_rhf_sph','h2o_uhf_cart','h2o_uhf_sph']

ok_opt = ['molden',
          'gaussian.log',
          'cclib',
          'gaussian.fchk',
          'aomix']

folder = ['molpro',
          'gaussian',
          'gaussian',
          'gaussian',
          'turbomole']

fileext = ['.molden', 
           '.inp.log', 
           '.inp.log', 
           '.fchk', 
           '/aomix.in']


for i in range(len(tests)):
  for j in range(len(folder)):

    skip = False
    if ok_opt[j] == 'cclib':
      try:
        __import__(ok_opt[j])
      except ImportError:
        skip = True

    if not skip:
      fid = os.path.join(output_folder,'%s/%s%s'%(folder[j],tests[i],fileext[j]))
      
      if 'uhf' in tests[i] and folder[j] == 'molpro':
        # Read the alpha input file
        qc = read.main_read(fid,itype=ok_opt[j],
                            all_mo=True,spin=None,i_md=0,interactive=False)
        # Extend the beta input file
        qc_b = read.main_read(fid,itype=ok_opt[j],
                            all_mo=True,spin=None,i_md=1,interactive=False)
        qc.mo_spec.extend(qc_b.mo_spec)
        qc.mo_spec.update()
      else:
        qc = read.main_read(fid ,itype=ok_opt[j],interactive=False,
                            all_mo=True,cclib_parser='Gaussian')
      
      dip = get_dipole_moment(qc,component=['x','y','z'])
      
      equal(dip, ref_dip[folder[j]])
    
'''
Old tests
'''

tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
folder = os.path.join(tests_home, '../outputs_for_testing/molpro')
filepath = os.path.join(folder, 'h2o_rhf_sph.molden')
qc = read.main_read(filepath, all_mo=True)
dip = get_dipole_moment(qc,component=['x','y','z'])

ref_dip = [0.00000000e+00,  -1.01130147e-16,   8.17259184e-01]
equal(dip, ref_dip)

qc.geo_spec += numpy.array([1,1,0])
dip = get_dipole_moment(qc,component=['x','y','z'])
equal(dip, ref_dip)

#Slightly move one atom and calculate dipoles again
qc.geo_spec[1] += numpy.array([1,1,0])
dip = get_dipole_moment(qc,component=['x','y','z'])
equal(dip, [0.25214432, -0.20529275,  1.09887067])
