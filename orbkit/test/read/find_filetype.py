from orbkit.read.tools import find_itype
import os, inspect

tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
folder = os.path.join(tests_home, '../outputs_for_testing')

files = {'fchk': 'gaussian/h2o_rhf_cart.fchk',
         'gaussian_log': 'gaussian/h2o_rhf_cart.inp.log',
         'gamess': 'gamess/formaldehyde.log',
         'molden': 'molpro/h2o_rhf_sph.molden',
         'aomix': 'turbomole/h2o_rhf_sph/aomix.in',
         'wfx': 'orca/1.wfx',
         'wfn': 'gamess/water_gamess-us.wfn'}

for fname in files:
  assert fname == find_itype(os.path.join(folder, files[fname]))
