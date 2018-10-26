'''
Module for systematic testing of Orbkit functions
'''

import unittest
from orbkit import __file__
from orbkit.display import display
import os
import sys
from orbkit import options

tests = ['units',
         'read/find_filetype',
         'read/read_multiple_input_files',
         'read/individual_readers',
         'detci/read/detci_readers',
         'detci/h3+',
         'native_io/io_numpy',
         'native_io/io_hdf5',
         'qcinfo/mo_select',
         'libcint_interface/libcint_basic',
         'libcint_interface/libcint_n2',
         'qcinfo/mo_spec',
         'analytical_properties/analytical_integral_norm',
         'analytical_properties/dipole',
         'grid_based/rho_compute',
         'grid_based/cube_files']

def check_import(filename):
  # This defines which modules should be checked for which folders
  # We try both to import the module and to find it within the PATH variable
  impdict = {'hdf5': 'h5py', 'libcint': 'lib/libcint'}
  refmod = None
  for module in impdict.keys():
    if module in filename:
      refmod = module
      break
  if refmod:
    run_test = False
    try:
      __import__(impdict[refmod])
      run_test = True
    except ImportError:
      pass
    for dirname in os.environ['PATH'].split(':'):
      if os.path.isfile(os.path.join(dirname, '{}.so'.format(impdict[refmod]))):
        run_test = True
        break
    return run_test
  else:
    return True

# Add skip conditions to tests here
def skip_tests():
  def decorator(func):
    def wrapper(self):
      can_import = False
      if check_import(self.filename):
        can_import = True
      if not can_import:
        self.skipTest('libcint not found')
      else:
        func(self)
    return wrapper
  return decorator

#This class is a modification of the ASE class of the same name
class ScriptTestCase(unittest.TestCase):
  def __init__(self, filename):
    unittest.TestCase.__init__(self, 'testfile')
    self.filename = filename + '.py'

  @skip_tests()
  def testfile(self):
    with open(self.filename) as fd:
      exec(compile(fd.read(), self.filename, 'exec'), {})

  def id(self):
    return self.filename

  def __str__(self):
    return self.filename.split('/')[-1][:-3]


def test():

  if options.numproc == 1:
    display('----------------------------------------------------------------------')
    display('              Testing serial ORBKIT functionality                     ')
    display('----------------------------------------------------------------------\n')
  else:
    display('----------------------------------------------------------------------')
    display('             Testing parallel ORBKIT functionality                    ')
    display('----------------------------------------------------------------------\n')


  ts = unittest.TestSuite()
  tsrun = unittest.TextTestRunner(verbosity=2)

  tests_home = __file__.split('__')[0] + 'test/'

  for test in tests:
    ts.addTest(ScriptTestCase(filename=os.path.abspath(tests_home + test)))

  results = tsrun.run(ts)

  sys.exit(len(results.errors + results.failures))
