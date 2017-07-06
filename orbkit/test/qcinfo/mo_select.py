from orbkit.read.high_level import main_read
from orbkit.test.tools import equal
from orbkit.qcinfo import QCinfo
from orbkit.orbitals import MOClass
from orbkit import options
import os, inspect

options.quiet = True

tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
folder = os.path.join(tests_home, '../read/outputs_for_testing')
filepath = os.path.join(folder, 'h2o_rhf_sph.molden')

qc = main_read(filepath, all_mo=True)

equal(qc.mo_spec.get_homo(), 4)
equal(qc.mo_spec.get_lumo(), 5)

usecases = {'homo-1:lumo+2': [3,4,5,6],
            'homo,lumo': [4,5],
            'homo': [4],
            'homo,lumo:lastbound+4': [4,5,6,7], #Last bound is HOMO
            'all_mo': range(len(qc.mo_spec)),
            '1.1,2.1,1.3,3.1': [0,1,2,3],
            '1.1 2.1 1.3 3.1': [0,1,2,3],
            '0 1 2 3': [0,1,2,3],
            '0       1            2 3': [0,1,2,3], #Test spaces and tabs
            '0,1 2:4': [0,1,2,3],
            '::2': range(0,len(qc.mo_spec),2),
            ':2': [0,1],
            '0:5:2': [0,2,4],
           }

for case in usecases:
  refmo = MOClass([qc.mo_spec[i] for i in usecases[case]])
  refmo.update()
  print(qc.mo_spec.select(case).sel_mo)
  equal(qc.mo_spec.select(case), refmo)

