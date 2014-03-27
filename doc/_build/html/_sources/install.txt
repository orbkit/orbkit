Installation Instructions for orbkit
====================================

For a proper execution of orbkit, the following Python modules are
required:

1) Python_ 2.6 - 2.7
2) SciPy_ Library of algorithms and mathematical tools
3) NumPy_ Library of high-level mathematical functions
4) h5py_ Interface to the HDF5 binary data format

The package h5py is not mandatory but strongly recommended.

.. _Python: http://www.python.org
.. _SciPy: http://www.scipy.org/
.. _NumPy: http://www.numpy.org/
.. _h5py: http://www.h5py.org/

Installation
------------

Currently, no setup.py file is available. So, if you want to use orbkit as a 
standalone program or a Python module, a manual installation is prerequisite. 
An instruction for linux using bash follows in the subsequent section.

Two ways will be described: An installation via git and an installation 
using the tarball provided.

First, choose the directory, where you want to install orbkit. 
In this example we will use the home directory.
Second, follow the subsequent instructions either for git or for tarball. 

- Using git:
  Clone the repository and set a path variable for orbkit to that directory:

.. code-block:: bash

    $ cd $HOME
    $ git clone git://git.code.sf.net/p/orbkit/code orbkit
    $ export ORBKITPATH=$HOME/orbkit

- Using a tarball:
  Download the latest orbkit release, extract the file (`v0.2.0` has to be replaced
  by your version number), and set a path variable for orbkit to that directory:

.. code-block:: bash

    $ cd $HOME
    $ wget http://sourceforge.net/projects/orbkit/files/latest/download 
    $ tar xzvf orbkit.v0.2.0.tar.gz
    $ export ORBKITPATH=$HOME/orbkit

**orbkit Modules:**

In order to use orbkit, you must add the orbkit directory to your PYTHONPATH 
environment variable either temporarily by typing

.. code-block:: bash

    $ export PYTHONPATH=$PYHONPATH:$ORBKITPATH

or permanently by adding this line to your ~/.bashrc file.

To use the orbkit standalone program, you have to modify additionally your PATH 
variable in the same way:

.. code-block:: bash

    $ export PATH=$PATH:$ORBKITPATH/tools
