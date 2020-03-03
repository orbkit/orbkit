#!/usr/bin/env python
import sys
from setuptools import find_packages, setup, Extension
from Cython.Distutils import build_ext
from os.path import join

import numpy

with open('README.rst', 'r') as fh:
  # Remove header
  for line in fh:
    if 'ORBKIT' in line:
      next(fh)
      next(fh)
      break
  
  long_description = fh.read()


if len(sys.argv) and 'bdist_wheel' in sys.argv[1:]:
  print("OpenMP for detCI disabled")
  ompcompileflags = []
  omplinkflags = []
else:
  ompcompileflags = ['-fopenmp']
  omplinkflags = ['-fopenmp']


setup(
  name='orbkit',
  packages=find_packages(),
  include_package_data=True,
  version='1.1.0.dev1',
  license='lgpl-3.0',
  description='A Toolbox for Post-Processing Quantum Chemical Wavefunction Data',
  long_description=long_description,
  long_description_content_type='text/x-rst',
  url='https://github.com/orbkit/orbkit',
  classifiers=[
      'Development Status :: 5 - Production/Stable',
      'Operating System :: OS Independent',
      'License :: OSI Approved :: '
      'GNU Lesser General Public License v3 or later (LGPLv3+)',
      'Programming Language :: Python :: 3',
      'Topic :: Scientific/Engineering :: Physics',
      'Topic :: Scientific/Engineering :: Chemistry',
      'Topic :: Scientific/Engineering :: Visualization'
  ],
  cmdclass={'build_ext': build_ext},
  entry_points={'console_scripts': ['orbkit = orbkit.main:run_standalone']},
  build_requires=['numpy'],
  install_requires=['numpy', 'scipy', 'matplotlib', 'h5py', 'setuptools'],
  ext_modules=[
    Extension(
      'orbkit.cy_grid',
      sources=['orbkit/cy_grid.pyx'],
      include_dirs=[numpy.get_include()],
    ),
    Extension(
      'orbkit.cy_core',
      sources=[
        'orbkit/cy_core.pyx', 'orbkit/c_grid-based.c',
        'orbkit/c_support.c'
      ],
      include_dirs=[numpy.get_include()],
      depends=[join('orbkit', '*.h')],
    ),
    Extension(
      'orbkit.cy_overlap',
      sources=[
        'orbkit/cy_overlap.pyx', 'orbkit/c_non-grid-based.c',
        'orbkit/c_support.c'
      ],
      include_dirs=[numpy.get_include()],
      depends=[join('orbkit', '*.h')],
    ),
    # detCI@ORBKIT
    Extension(
      'orbkit.detci.cy_occ_check',
      sources=['orbkit/detci/cy_occ_check.pyx'],
      include_dirs=[numpy.get_include()],
    ),
    Extension(
      'orbkit.detci.cy_ci',
      sources=['orbkit/detci/cy_ci.pyx'],
      extra_compile_args=ompcompileflags,
      extra_link_args=omplinkflags,
      include_dirs=[numpy.get_include()],
    ),
    # Libcint
    Extension(
      'orbkit.libcint_interface.cy_mo_integrals',
      sources=['orbkit/libcint_interface/cy_mo_integrals.pyx'],
      include_dirs=[numpy.get_include()],
    ),
  ],
)
