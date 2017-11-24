import numpy
import os, inspect

from orbkit import read
from orbkit import analytical_integrals as ai
from orbkit import options

from orbkit.test.tools import equal

options.quiet = True
tol = 1e-4
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
      ao_overlap_matrix = ai.get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec)
      # Loop over spin if unrestricted
      for spin in ['alpha','beta']:
        mos = qc.mo_spec['all_mo ' + spin]
        if mos != []:
          moom = ai.get_mo_overlap_matrix(mos,mos,ao_overlap_matrix,numproc=options.numproc)
          equal(moom, numpy.eye(len(moom)),tol=tol)
