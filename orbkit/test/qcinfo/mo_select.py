from orbkit.read.high_level import main_read
from orbkit.test.tools import equal
import numpy
from orbkit.qcinfo import QCinfo
from orbkit.orbitals import MOClass
from orbkit import options
import os, inspect
import warnings
warnings.filterwarnings("ignore")

options.quiet = True

tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
folder = os.path.join(tests_home, '../outputs_for_testing/gaussian')
filepath = os.path.join(folder, 'h2o_uhf_cart.inp.log')

qc = main_read(filepath, all_mo=True)
equal(qc.mo_spec.get_homo(), 9)
equal(qc.mo_spec.get_lumo(), 10)
equal(qc.mo_spec.get_lastbound(), qc.mo_spec.get_homo())

homo = qc.mo_spec.get_homo()
lumo = qc.mo_spec.get_lumo()
lastbound = qc.mo_spec.get_lastbound()

alpha = numpy.array(qc.mo_spec.alpha_index)
beta = numpy.array(qc.mo_spec.beta_index)

usecases = {'homo-1:lumo+2': range(homo-1,lumo+2),
            'homo,lumo': [homo,lumo],
            'homo': [homo],
            'homo,lumo:lastbound+4': [homo-1] + list(range(lumo,lastbound+4)), #Last bound is HOMO
            'all_mo': range(len(qc.mo_spec)),
            '1.A1 alpha': [0],
            '1.A1 2.A1 1.B2 3.A1': range(8),
            '0 1 2 3': range(4),
            '0       1            2 3': range(4), #Test spaces and tabs
            '0,1 2:4': range(4),
            '::2': range(0,len(qc.mo_spec),2),
            ':2': range(2),
            '2:': range(2,len(qc.mo_spec)),
            '0:5:2': list(range(0,5,2)),
            '1.A1 2.A1 1.B2 3.A1 alpha': alpha[:4],
            '1.A1 2.A1 1.B2 3.A1 beta': beta[:4],
            '*.A1': [i for i,s in enumerate(qc.mo_spec.get_sym()) if 'A1' in s],
            '*.A1 1.B2 2.B2': [i for i,s in enumerate(qc.mo_spec.get_sym()) if 'A1' in s or s in ['1.B2', '2.B2']],
            '*.A1 beta': [i for i,s in enumerate(qc.mo_spec.get_sym()) if 'A1' in s and i in beta],
            'homo,lumo alpha': alpha[[homo//2,lumo//2]],
            'all_mo alpha': alpha
           }

for case in usecases:
  refmo = MOClass([qc.mo_spec[i] for i in usecases[case]])
  refmo.update()
  equal(qc.mo_spec[case], refmo)

# Test lists of strings
refmo = MOClass([qc.mo_spec[i] for i in alpha[:4]])
refmo.update()
equal(qc.mo_spec[['1.A1', '2.A1', '1.B2', '3.A1', 'all_alpha']], refmo)

refmo = MOClass([qc.mo_spec[i] for i in [0, 1, 2, 3, 4, 6]])
refmo.update()
equal(qc.mo_spec[[['1.A1', '2.A1'], ['1.B2', '3.A1', 'alpha']]], refmo)

refmo = MOClass([qc.mo_spec[i] for i in range(12)])
refmo.update()
equal(qc.mo_spec[list(map(str,range(12)))], refmo)

# Test lists of integers
refmo = MOClass([qc.mo_spec[i] for i in range(12)])
refmo.update()
equal(qc.mo_spec[range(12)], refmo)

refmo = MOClass([qc.mo_spec[i] for i in range(12)])
refmo.update()
equal(qc.mo_spec[[list(range(5)), list(range(5,12))]], refmo)

# Test mixed lists
refmo = MOClass([qc.mo_spec[i] for i in list(alpha[:8])])
refmo.update()
equal(qc.mo_spec[[['1.A1','2.A1','1.B2','3.A1', 'all_alpha'] + list(alpha[4:8])]], refmo)

refmo = MOClass([qc.mo_spec[i] for i in list(alpha[:12])])
refmo.update()
equal(qc.mo_spec[[['1.A1','2.A1','1.B2','3.A1'] + ['all_alpha'] + list(alpha[4:8])+ list(map(str,list(alpha[8:12])))]], refmo)







