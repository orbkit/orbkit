from orbkit.read.tools import find_itype
import os, inspect

tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
folder = os.path.join(tests_home, 'outputs_for_testing')

files = {'fchk': 'h2o_rhf_cart.fchk',
         'gaussian_log': 'h2o_rhf_cart.inp.log',
         'molden': 'h2o_rhf_sph.molden',
         'aomix': 'aomix.in',
         'wfx': '1.wfx',
         'wfn': 'water_gamess-us.wfn'}

for fname in files:
  assert fname == find_itype(os.path.join(folder, files[fname]))
