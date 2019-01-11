from orbkit.read.molden import read_molden
from orbkit.read.gamess import read_gamess
from orbkit.read.gaussian_fchk import read_gaussian_fchk
from orbkit.read.gaussian_log import read_gaussian_log
from orbkit.read.aomix import read_aomix
from orbkit.read.wfx import read_wfx
from orbkit.read.wfn import read_wfn
from orbkit.read.cclib_parser import read_with_cclib
from orbkit.test.tools import equal
from orbkit import options
import numpy
import os, inspect

tests_home = os.path.dirname(inspect.getfile(inspect.currentframe()))
folder = os.path.join(tests_home, '../outputs_for_testing')

options.quiet = True

files = {
         'fchk': 'gaussian/h2o_rhf_cart.fchk',
         'cclib': 'gaussian/h2o_rhf_cart.inp.log',
         'gamess': 'gamess/formaldehyde.log',
         'gaussian_log': 'gaussian/h2o_rhf_cart.inp.log',
         'molpro_molden': 'molpro/h2o_rhf_sph.molden',
         'psi4_molden': 'psi4/lih_cis_aug-cc-pVTZ.out.default.molden',
         'aomix': 'turbomole/h2o_rhf_cart/aomix.in',
         'wfx': 'orca/1.wfx',
         'tar.gz': 'NaCl_molden_files.tar.gz',
         'wfn': 'gamess/water_gamess-us.wfn',
         }

readers = {'fchk': read_gaussian_fchk,
           'gaussian_log': read_gaussian_log,
           'cclib': read_with_cclib,
           'gamess': read_gamess,
           'molpro_molden': read_molden,
           'psi4_molden': read_molden,
           'aomix': read_aomix,
           'wfx': read_wfx,
           'tar.gz': read_molden,
           'wfn': read_wfn}

refeigen = {'fchk': [-2.05704730E+01, -1.27755361E+00, -6.41980010E-01, -5.33424332E-01],
            'gaussian_log': [-20.57047,  -1.27755,  -0.64198,  -0.53342],
            'cclib': [-20.57047,  -1.27755,  -0.64198,  -0.53342],
            'gamess': [-20.5775,  -11.3375,  -1.4227,  -0.8676],
            'molpro_molden': [-20.5692, -1.2772, -0.6416, -0.5331],
            'psi4_molden': [-2.4492131459, -0.2996316473, -0.0077681166, 0.0160817783],
            'aomix': [-20.57047,  -1.27755,  -0.64198,  -0.53342],
            'wfx': [-2.06529440471407E+01, -1.17607391494576E+00, -5.09563645324521E-01, -4.41339273759071E-01],
            'tar.gz': [-104.9469, -40.5349, -10.6711, -8.1357],
            'wfn': [-20.56804075, -1.27715750, -0.64186203, -0.53341030]}

refexp = {'fchk': [1.30100000E+01, 1.22000000E-01, 7.27000000E-01, 1.17200000E+04],
          'gaussian_log': [0.1301000000E+02, 0.1220000000E+00, 0.7270000000E+00, 0.1172000000E+05],
          'cclib': [0.1301000000E+02, 0.1220000000E+00, 0.7270000000E+00, 0.1172000000E+05],
          'gamess': [8236.0000000,  8236.0000000,  0.9059000,  0.1285000],
          'molpro_molden': [0.1301000000E+02, 0.1301000000E+02, 0.7270000000E+00, 0.1301000000E+02],
          'psi4_molden': [5988.0000000000, 5988.0000000000, 0.0750900000, 0.0283200000],
          'aomix': [0.13010000000000E+02,  0.12200000000000E+0,  0.72700000000000E+00,  0.11720000000000E+05],
          'wfx': [2.26617677850000E+03, 3.40870101910000E+02, 7.73631351670000E+01, 2.14796449400000E+01],
          'tar.gz': [0.4230000000E+06, 0.4230000000E+06, 0.4230000000E+06, 0.4230000000E+06],
          'wfn': [1.1720000E+04, 1.7590000E+03, 4.0080000E+02, 1.1370000E+02]}

refcontrac = {'fchk': [3.34987264E-02, 1.00000000E+00, 1.00000000E+00, 7.11864434E-04],
              'gaussian_log': [0.3349872639E-01, 0.1000000000E+01, 0.1000000000E+01, 0.7118644339E-03],
              'cclib': [0.3349872639E-01, 0.1000000000E+01, 0.1000000000E+01, 0.7118644339E-03],
              'gamess': [0.000542430189, -0.000196392234, 1.000000000000, 1.000000000000],
              'molpro_molden': [0.1968498999E-01, 0.0000000000E+00, 0.1000000000E+01, 0.1968498999E-01],
              'aomix': [0.33498726389998E-01,  0.10000000000000E+01,  0.10000000000000E+01,  0.70964594651845E-03],
              'tar.gz': [0.1806182306E-04, -0.4406529065E-05, 0.6630189916E-06, 0.0000000000E+00],
              'wfx': [1, 1, 1, 1],
              'wfn': [1, 1, 1, 1]}

refgeo = {'tar.gz': [ 0.        ,  0.        , -2.54176518],
          'wfx': [ 0.        ,  0.        , -0.22204202],
          'wfn': [ 0.        ,  0.        , -0.13302525],
          'gamess': [ 0.        ,  0.        , -0.99504027],
          'gaussian_log': [ 0.        ,  1.68211892, -0.95104436],
          'aomix': [ 0.        , -1.68211949,  1.05577926],
          'cclib': [ 0.        ,  1.68211892, -0.95104436],
          'fchk': [ -6.16297582e-32,   1.68211949e+00,  -9.51043615e-01],
          'molpro_molden': [ 0.        , -1.68211948,  1.05577926],
          'psi4_molden': [ 0. , 0.,  -0.386892149961],
          }

for fname in files:

  skip = False
  if fname == 'cclib':
    try:
      __import__(fname)
    except ImportError:
      skip = True

  if not skip:
    qcinfo = readers[fname](os.path.join(folder, files[fname]),cclib_parser='Gaussian',all_mo=True)
    assert qcinfo.ao_spec.get_ao_num() == qcinfo.mo_spec.get_coeffs().shape[1]
    e_list = numpy.zeros(4, dtype=float)
    coeffs_list = numpy.zeros(4, dtype=float)
    contrac_list = numpy.zeros(4, dtype=float)
    geo_list = qcinfo.geo_spec[0]
    for i in range(4):
      e_list[i] = qcinfo.mo_spec[i]['energy']
      coeffs_list[i] = qcinfo.ao_spec[i]['coeffs'][0,0]
      contrac_list[i] = qcinfo.ao_spec[i]['coeffs'][0,1]
    for ref in [[e_list, refeigen], [coeffs_list, refexp], [contrac_list, refcontrac], [geo_list,refgeo]]:
      if fname in ref[1].keys():
        equal(ref[0], ref[1][fname])



