.. _installation-instructions:

Installation Instructions
=========================

In this section, we present the manual installation of the development version of ORBKIT. For the stable release, refer to `https://orbkit.github.io/install.html <https://orbkit.github.io/install.html>`__.

.. contents:: Table of Contents:
  :local:
  :depth: 2

Prerequisites
-------------

For a proper execution of ORBKIT, the following Python modules are
required:

1) Python_ 2.6 - 2.7, 3.X
2) Cython_
3) NumPy_ Library of high-level mathematical functions
4) SciPy_ Library of algorithms and mathematical tools
5) h5py_ Interface to the HDF5 binary data format
6) mayavi_ Tool for 3D scientific data visualization (optional)

The package h5py is not mandatory but strongly recommended.

.. _Python: http://www.python.org
.. _Cython: http://cython.org/
.. _SciPy: http://www.scipy.org/
.. _NumPy: http://www.numpy.org/
.. _h5py: http://www.h5py.org/
.. _mayavi: http://docs.enthought.com/mayavi/mayavi/index.html

Instructions
------------

ORBKIT needs to be installed manually, i.e.,
the Cython modules need to be pre-compiled and some 
environment variables need to be set. 
In the following, we describe this procedure exemplary 
for the different platforms.

Irrespective of the system you're using, you should test
the code after completing the installation process.

* First, to test the serial version of the code run::

    $ orbkit test

* If tests run successfully also test the parallel version::

    $ orbkit test -p 2

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
bug fixes and features. 

  * Using git:

    Clone the repository::

        $ git clone https://github.com/orbkit/orbkit.git
    
    Go to the orbkit folder and checkout the development branch::
    
        $ cd orbkit
        $ git checkout development
        $ cd ..

  * **OR:** Using a zip archive:

    Download the latest ORBKIT (development) release and extract the file::

        $ wget https://github.com/orbkit/orbkit/archive/development.zip
        $ unzip orbkit-development.zip
        $ mv orbkit-development orbkit

Set an environment variable to this directory::

    $ export ORBKITPATH=$HOME/orbkit

Now, you have to build to ORBKIT::

    $ cd $ORBKITPATH
    $ python setup.py build_ext --inplace clean

In order to use ORBKIT, you have to add the ORBKIT directory to your ``$PYTHONPATH``
environment variable either temporarily by typing::

    $ export PYTHONPATH=$PYHONPATH:$ORBKITPATH

or permanently by adding these lines to your ~/.bashrc file.

To use ORBKIT as a standalone program, you have to modify your 
$PATH variable in the same way::

    $ export PATH=$PATH:$ORBKITPATH/tools

or permanently by adding these lines to your ~/.bashrc file.

Windows
.......

We have tested ORBKIT on Windows using the free Visual Studio 2015 Community Edition 
(https://www.visualstudio.com/en-us/downloads/download-visual-studio-vs.aspx)
and the free version of the Python environment Entought Canopy 
(https://www.enthought.com/products/canopy/). 

Download and unzip the newest (development) version of ORBKIT (or use git and clone the newest version):
from 
  
  https://github.com/orbkit/orbkit/archive/development.zip

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
