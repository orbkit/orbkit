from orbkit.read.high_level import main_read
from orbkit.test.tools import equal
from orbkit.qcinfo import QCinfo
from orbkit.orbitals import MOClass
from orbkit import options
import os, inspect

options.quiet = True

tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
folder = os.path.join(tests_home, '../read/outputs_for_testing')
filepath = os.path.join(folder, 'h2o_uhf_cart.inp.log')

qc = main_read(filepath, all_mo=True)
equal(qc.mo_spec.get_homo(), 4)
equal(qc.mo_spec.get_lumo(), 5)

usecases = {'homo-1:lumo+2': [3,4,5,6],
            'homo,lumo': [4,5],
            'homo': [4],
            'homo,lumo:lastbound+4': [4,5,6,7], #Last bound is HOMO
            'all_mo': range(len(qc.mo_spec)),
            '1.A1_a 2.A1_a 1.B2_a 3.A1_a': [0,1,2,3],
            '0 1 2 3': [0,1,2,3],
            '0       1            2 3': [0,1,2,3], #Test spaces and tabs
            '0,1 2:4': [0,1,2,3],
            '::2': range(0,len(qc.mo_spec),2),
            ':2': [0,1],
            '2:': range(2,len(qc.mo_spec)),
            '0:5:2': [0,2,4],
            '1.A1_a 2.A1_a 1.B2_a 3.A1_a alpha': [0,1,2,3],
            'homo,lumo alpha': [4,5],
            'all_mo alpha': list(range(len(qc.mo_spec)//2))
           }

for case in usecases:
  refmo = MOClass([qc.mo_spec[i] for i in usecases[case]])
  refmo.update()
  equal(qc.mo_spec.select(case, flatten_input=True), refmo)

# Test lists of strings
refmo = MOClass([qc.mo_spec[i] for i in [0,1,2,3]])
refmo.update()
equal(qc.mo_spec.select(['1.A1_a', '2.A1_a', '1.B2_a', '3.A1_a', 'alpha'], flatten_input=True), refmo)

refmo = MOClass([qc.mo_spec[i] for i in [0,1,2,3]])
refmo.update()
equal(qc.mo_spec.select([['1.A1_a', '2.A1_a'], ['1.B2_a', '3.A1_a', 'alpha']], flatten_input=True), refmo)

# Test lists of integers
refmo = MOClass([qc.mo_spec[i] for i in range(12)])
refmo.update()
equal(qc.mo_spec.select(list(range(12)), flatten_input=True), refmo)

refmo = MOClass([qc.mo_spec[i] for i in range(12)])
refmo.update()
equal(qc.mo_spec.select([list(range(5)), list(range(5,12))], flatten_input=True), refmo)






