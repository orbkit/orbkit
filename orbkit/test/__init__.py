'''
Module for systematic testing of Orbkit functions
'''

import unittest
import subprocess
from orbkit import __file__
from orbkit.display import display
import os
import sys
import shutil
from orbkit import options

tests = ['units',
         'read/find_filetype',
         'read/read_multiple_input_files',
         'read/individual_readers',
         'detci/read/detci_readers',
         'detci/h3+',
         'native_io/numpy',
         'native_io/hdf5',
         'qcinfo/mo_select',
         'qcinfo/mo_spec',
         'qcinfo/ao_spec',
         'analytical_properties/analytical_integral_norm',
         'analytical_properties/dipole',
         'grid_based/rho_compute']

def check_import(filename):
  # This defines which modules should be checked for which folders
  impdict = {'hdf5': 'h5py'}
  refmod = None
  for module in impdict.keys():
    if module in filename:
      refmod = module
      break
  if refmod:
    if sys.version_info.major == 2:
      import imp
      imp_check = imp.find_module
    else:
      import importlib
      if sys.version_info.minor <= 3:
        imp_check = importlib.find_loader
      else:
        imp_check = importlib.util.find_spec
    try:
      imp_check(impdict[refmod])
      return True
    except ImportError:
      return False
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
        self.skipTest('skipped')
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
