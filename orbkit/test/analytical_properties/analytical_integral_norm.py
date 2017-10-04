import numpy
import os, inspect

from orbkit import read
from orbkit import analytical_integrals as ai
from orbkit import options

from orbkit.test.tools import equal

options.quiet = True

tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
folder = os.path.join(tests_home, '../outputs_for_testing/molpro')
filepath = os.path.join(folder, 'h2o_rhf_sph.molden')
qc = read.main_read(filepath, all_mo=True)

ao_overlap_matrix = ai.get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec)

moom = ai.get_mo_overlap_matrix(qc.mo_spec,qc.mo_spec,ao_overlap_matrix,numproc=options.numproc)

equal(moom, numpy.eye(len(moom)))
