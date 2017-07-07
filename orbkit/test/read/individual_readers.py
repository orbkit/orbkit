from orbkit.read.molden import read_molden
#from orbkit.read.gamess import read_gamess
from orbkit.read.gaussian_fchk import read_gaussian_fchk
from orbkit.read.gaussian_log import read_gaussian_log
from orbkit.read.aomix import read_aomix
from orbkit.read.wfx import read_wfx
from orbkit.read.wfn import read_wfn
from orbkit.read.cclib import read_with_cclib
from orbkit.test.tools import equal
from orbkit import options
import numpy
import os, inspect

tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
folder = os.path.join(tests_home, 'outputs_for_testing')

options.quiet = True

#Gamess Reader is not tested at this point
#We need to test cclib as well

files = {'fchk': 'h2o_rhf_cart.fchk',
         'gaussian_log': 'h2o_rhf_cart.inp.out',
         'molden': 'h2o_rhf_sph.molden',
         'aomix': 'aomix.in',
         'wfx': '1.wfx',
         'tar.gz': 'NaCl_molden_files.tar.gz',
         'wfn': 'water_gamess-us.wfn'}

readers = {'fchk': read_gaussian_fchk,
           'gaussian_log': read_gaussian_log,
           'molden': read_molden,
           'aomix': read_aomix,
           'wfx': read_wfx,
           'tar.gz': read_molden,
           'wfn': read_wfn}

refeigen = {'fchk': numpy.array([-2.05704730E+01, -1.27755361E+00, -6.41980010E-01, -5.33424332E-01]),
            'gaussian_log': numpy.array([-20.57047,  -1.27755,  -0.64198,  -0.53342]),
            'molden': numpy.array([-20.5692, -1.2772, -0.6416, -0.5331]),
            'aomix': numpy.array([-20.57047,  -1.27755,  -0.64198,  -0.53342]),
            'wfx': numpy.array([-2.06529440471407E+01, -1.17607391494576E+00, -5.09563645324521E-01, -4.41339273759071E-01]),
            'tar.gz': numpy.array([-104.9469, -40.5349, -10.6711, -8.1357]),
            'wfn': numpy.array([-20.56804075, -1.27715750, -0.64186203, -0.53341030])}

refexp = {'fchk': numpy.array([1.30100000E+01, 1.22000000E-01, 7.27000000E-01, 1.17200000E+04]),
          'gaussian_log': numpy.array([0.1301000000E+02, 0.1220000000E+00, 0.7270000000E+00, 0.1172000000E+05]),
          'molden': numpy.array([0.1301000000E+02, 0.1301000000E+02, 0.7270000000E+00, 0.1301000000E+02]),
          'aomix': numpy.array([0.13010000000000E+02,  0.12200000000000E+0,  0.72700000000000E+00,  0.11720000000000E+05]),
          'wfx': numpy.array([2.26617677850000E+03, 3.40870101910000E+02, 7.73631351670000E+01, 2.14796449400000E+01]),
          'tar.gz': numpy.array([0.4230000000E+06, 0.4230000000E+06, 0.4230000000E+06, 0.4230000000E+06]),
          'wfn': numpy.array([1.1720000E+04, 1.7590000E+03, 4.0080000E+02, 1.1370000E+02])}

refcontrac = {'fchk': numpy.array([3.34987264E-02, 1.00000000E+00, 1.00000000E+00, 7.11864434E-04]),
              'gaussian_log': numpy.array([0.3349872639E-01, 0.1000000000E+01, 0.1000000000E+01, 0.7118644339E-03]),
              'molden': numpy.array([0.1968498999E-01, 0.0000000000E+00, 0.1000000000E+01, 0.1968498999E-01]),
              'aomix': numpy.array([0.33498726389998E-01,  0.10000000000000E+01,  0.10000000000000E+01,  0.70964594651845E-03]),
              'tar.gz': numpy.array([0.1806182306E-04, -0.4406529065E-05, 0.6630189916E-06, 0.0000000000E+00]),
              'wfx': numpy.array([1, 1, 1, 1]),
              'wfn': numpy.array([1, 1, 1, 1])}

#I'm not shure that .fchk files are read correctly. Seems good but I
#don't really know... Can somemone please check?

for fname in files:
  qcinfo = readers[fname](os.path.join(folder, files[fname]))
  e_list = numpy.zeros(4, dtype=float)
  coeffs_list = numpy.zeros(4, dtype=float)
  contrac_list = numpy.zeros(4, dtype=float)
  for i in range(4):
    e_list[i] = qcinfo.mo_spec[i]['energy']
    coeffs_list[i] = qcinfo.ao_spec[i]['coeffs'][0,0]
    contrac_list[i] = qcinfo.ao_spec[i]['coeffs'][0,1]
  for ref in [[e_list, refeigen], [coeffs_list, refexp], [contrac_list, refcontrac]]:
    equal(ref[0], ref[1][fname])





