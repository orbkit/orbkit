.. image:: doc/orbkit_small.png
    :align: center

ORBKIT
======

ORBKIT is a parallel Python program package for post-processing 
wave function data from output files of quantum chemical programs.

The computational capabilities of ORBKIT range from grid-based quantities, e.g., molecular orbitals or 
electron density, to non grid-based quantities, for instance, Mulliken population charges or
analytical overlap integrals between molecular orbitals. 
There are several options and features to control the respective calculations, like grid types and parameters. 
The required data can be extracted from MOLPRO (Molden File Format), 
TURBOMOLE (AOMix file format), GAMESS-US, PROAIMS/AIMPAC (wfn/wfx file format), and Gaussian (.log File and Formatted Checkpoint File)
output files. Futhermore, an interface to cclib, a parser for quantum chemical logfiles, is provided.

**NEW:** `detCI\@ORBKIT`__ extends ORBKIT's functionality to multi-determinantal wave functions.

__ http://orbkit.github.io/dev/detci/index.html

ORBKIT's documentation may be found at `http://orbkit.github.io <http://orbkit.github.io/dev>`_


Support
-------

If you need help for the usage of ORBKIT, please do not hesitate to contact the 
ORBKIT support team via 

https://github.com/orbkit/orbkit/issues

Citation
--------

If you use ORBKIT in your work, please cite:

  Gunter Hermann, Vincent Pohl, Jean Christophe Tremblay, Beate Paulus, Hans-Christian Hege, and Axel Schild,
  `"ORBKIT: A Modular Python Toolbox for Cross-Platform Postprocessing of Quantum Chemical Wavefunction Data" <http://dx.doi.org/10.1002/jcc.24358>`_, 
  *J. Comput. Chem.* **2016**, *37*, 1511-1520.

If you use detCI\@ORBKIT in your work, please additionally cite:

  Vincent Pohl, Gunter Hermann, and Jean Christophe Tremblay,
  `"An Open-Source Framework for Analyzing N-Electron Dynamics. I. Multideterminantal Wave Functions" <http://dx.doi.org/10.1002/jcc.24792>`_, 
  *J. Comput. Chem.* **2017**, *38*, 1515-1527.

  Vincent Pohl, Gunter Hermann, and Jean Christophe Tremblay,
  `"An Open-Source Framework for Analyzing N-Electron Dynamics. II. Hybrid Density Functional Theory/Configuration Interaction Methodology" <http://dx.doi.org/10.1002/jcc.24896>`_, 
  *J. Comput. Chem.* **2017**, `DOI:10.1002/jcc.24896 <http://dx.doi.org/10.1002/jcc.24896>`_.

The papers are also freely available on arXiv (`ORBKIT <https://arxiv.org/abs/1601.03069>`_, `detCI\@ORBKIT_I <https://arxiv.org/abs/1701.06885>`_, and `detCI\@ORBKIT_II <https://arxiv.org/abs/1704.08137>`_) and a BibTex file may be
found in `doc/orbkit.bib <doc/orbkit.bib>`_.



Installation Requirements
-------------------------

For a proper execution of ORBKIT, the following Python modules are required:

1) Python 2.6 - 2.7, Python 3.x (http://www.python.org)
2) Cython (http://cython.org/)
3) NumPy Library of high-level mathematical functions (http://www.numpy.org/)
4) SciPy Library of algorithms and mathematical tools (http://www.scipy.org/)
5) h5py Interface to the HDF5 binary data format (http://www.h5py.org/)
6) Mayavi Tool for 3D scientific data visualization (optional, http://code.enthought.com/projects/mayavi/)

The package h5py is not mandatory but strongly recommended.

Installation
------------

ORBKIT needs to be installed manually, i.e.,
the Cython modules need to be pre-compiled and some 
environment variables need to be set. 
In the following, we describe this procedure exemplary 
for the different platforms.

Linux and Mac OS X
..................

The manual installation of ORBKIT is simple and can 
be carried out using ``bash`` as follows:

Choose the directory, where you want to install ORBKIT. Open a terminal window, 
e.g. ``gnome-terminal``, and navigate to this directory. In this example we 
will use the home directory. If you use a different directory simply replace 
``$HOME`` by your preferred folder throughout the whole section::

    $ cd $HOME

Get a copy of ORBKIT, either with git or using a zip archive. It is strongly
recommended to use git, since this version always contains the newest 
bug fixes and features. If git is not available on your system, the newest 
version can additionally be cloned from https://github.com/orbkit/orbkit.

  * Using git:

    Clone the repository::

        $ git clone https://github.com/orbkit/orbkit.git
  * **OR:** Using a zip archive:

    Download the latest ORBKIT release and extract the file::

        $ wget https://github.com/orbkit/orbkit/archive/cython.zip
        $ unzip orbkit-cython.zip
        $ mv orbkit-cython orbkit

Set an environment variable to this directory::

    $ export ORBKITPATH=$HOME/orbkit

Now, you have to build to ORBKIT::

    $ cd $ORBKITPATH
    $ python setup.py build_ext --inplace clean

In order to use ORBKIT, you have to add the ORBKIT directory to your ``$PYTHONPATH``
environment variable either *temporarily* by typing::

    $ export PYTHONPATH=$PYHONPATH:$ORBKITPATH

or permanently by adding these two lines to your ~/.bashrc file::

    $ export ORBKITPATH=$HOME/orbkit
    $ export PYTHONPATH=$PYHONPATH:$ORBKITPATH

To use ORBKIT as a standalone program, you have to modify your 
$PATH variable in the same way::

    $ export PATH=$PATH:$ORBKITPATH/tools

Windows
.......

We have tested ORBKIT on Windows using the free Visual Studio 2015 Community Edition 
(https://www.visualstudio.com/en-us/downloads/download-visual-studio-vs.aspx)
and the free version of the Python environment Entought Canopy 
(https://www.enthought.com/products/canopy/). 

Download and unzip the newest version of ORBKIT (or use git and clone the newest version):
from 
  
  https://github.com/orbkit/orbkit/archive/cython.zip

In the following, we assume that ORBKIT can be found at ``C:\orbkit``

Install Visual Studio 2015 including the Python-Tools for Visual Studio.
After installing Canopy (and using it as your default Python environment), 
install the required Python packages using the graphical package manager. 

If you are using the 64-bit version of Canopy (Python), please start the
``VS2013 x64 Native Tools Command Prompt``. For 32-bit, start the 
``VS2013 x86 Native Tools Command Prompt``.

Navigate to the ORBKIT folder::

  > cd C:\orbkit

Set some environment variables and build ORBKIT::

  > SET DISTUTILS_USE_SDK=1
  > SET MSSdk=1
  > python setup.py build_ext --inplace --compiler=msvc clean

Finally, you have to set the PYTHONPATH and the PATH variables to use ORBKIT.

Licence Note
------------

ORBKIT is free software: you can redistribute it and/or modify it under the 
terms of the GNU Lesser General Public License as published by the Free Software 
Foundation, either version 3 of the License, or any later version.

ORBKIT is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along 
with ORBKIT. If not, see <http://www.gnu.org/licenses/>.

Copyright (C) 2017, Gunter Hermann, Vincent Pohl, Lukas Eugen Marsoner Steinkasserer, Axel Schild, and Jean Christophe Tremblay.