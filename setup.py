#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("orbkit.cy_core",
                    sources=["orbkit/cy_core.pyx", "orbkit/c_support.c"],
                             include_dirs=[numpy.get_include()]),
                   Extension("orbkit.cy_grid",
                    sources=["orbkit/cy_grid.pyx"],
                             include_dirs=[numpy.get_include()])],
)
