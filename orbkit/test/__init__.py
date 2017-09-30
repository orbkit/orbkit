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
         'qcinfo/new_qcinfo_save_function',
         'qcinfo/mo_select',
         'integrals/nh3',
         'integrals/n2',
         'qcinfo/mo_spec',
         'qcinfo/ao_spec',
         'analytical_properties/analytical_integral_norm',
         'analytical_properties/dipole',
         'grid_based/rho_compute']

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
        return self.filename.split('/')[-1][:-3]

    def id(self):
        return self.filename


def test():

  ts = unittest.TestSuite()

  if options.numproc == 1:
    display('----------------------------------------------------------------------')
    display('              Testing serial ORBKIT functionality                     ')
    display('----------------------------------------------------------------------\n')
  else:
    display('----------------------------------------------------------------------')
    display('             Testing parallel ORBKIT functionality                    ')
    display('----------------------------------------------------------------------\n')

  tsrun = unittest.TextTestRunner(verbosity=2)

  tests_home = __file__.split('__')[0] + 'test/'

  for test in tests:
    ts.addTest(ScriptTestCase(filename=os.path.abspath(tests_home + test)))

  results = tsrun.run(ts)

  sys.exit(len(results.errors + results.failures))
