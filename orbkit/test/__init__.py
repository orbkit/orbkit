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

tests = ['units',
         'analytical_properties/analytical_integral_norm',
         'read/find_filetype',
         'read/read_multiple_input_files',
         'read/specific_readers',
         'detci/read/detci_readers']

def clean(testdir):
  if os.path.isdir(testdir):
    shutil.rmtree(testdir)
  return

#This class is a modification of the ASE class of the same name
class ScriptTestCase(unittest.TestCase):
    def __init__(self, filename):
        unittest.TestCase.__init__(self, 'testfile')
        self.filename = filename + '.py'

    def testfile(self):
        try:
            with open(self.filename) as fd:
                exec(compile(fd.read(), self.filename, 'exec'), {})
        except KeyboardInterrupt:
            raise RuntimeError('Keyboard interrupt')

    def id(self):
        return self.filename

    def __str__(self):
        return self.filename.split('/')[-1]

    def id(self):
        return self.filename


def test():

  ts = unittest.TestSuite()

  display('----------------------------------------------------------------------')
  display('                     Testing ORBKIT functionality                     ')
  display('----------------------------------------------------------------------\n')

  tsrun = unittest.TextTestRunner(verbosity=2)

  tests_home = __file__.split('__')[0] + 'test/'

  for test in tests:
    ts.addTest(ScriptTestCase(filename=os.path.abspath(tests_home + test)))

  testdir = tests_home + 'test_tmp'
  clean(testdir)
  os.mkdir(testdir)
  os.chdir(testdir)

  results = tsrun.run(ts)
  clean(testdir)

  sys.exit(len(results.errors + results.failures))





