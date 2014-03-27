orbkit
======
Copyright (C) 2014, Gunter Hermann, Vincent Pohl, and Axel Schild.

orbkit is a parallel Python program package for post-processing 
wavefunction data extracted from output files of MOLPRO (Molden File Format), 
TURBOMOLE (Molden file format), GAMESS US, and Gaussian (Formatted Checkpoint File). 

Computational capabilities using arbitrary grids:

- Atomic orbitals, molecular orbitals, electron density, 
  and the spatial derivatives of these quantities
- Reduced electron density
- Electron density for selected molecular orbitals
- Atom-projected electron density
- Molecular orbital transition electronic flux density

If you use orbkit in your work, please cite it as follows:

Gunter Hermann, Vincent Pohl, Axel Schild: orbkit: A Toolbox for Post-Processing 
Quantum Chemical Wavefunction Data, available via http://sourceforge.net/p/orbkit (2014).

Website: http://sourceforge.net/p/orbkit

Licence Note
============

orbkit is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or any later version.

orbkit is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with orbkit.  If not, see <http://www.gnu.org/licenses/>.

Installation Requirements
=========================

For a proper execution of orbkit, the following Python modules are
required:

1) Python 2.6 - 2.7 (http://www.python.org)
2) SciPy Library of algorithms and mathematical tools (http://www.scipy.org/)
3) NumPy Library of high-level mathematical functions (http://www.numpy.org/)
4) h5py Interface to the HDF5 binary data format (http://www.h5py.org/)

The package h5py is not mandatory but strongly recommanded.

Installation
============

Currently, no setup.py file is available. So, if you want to use orbkit as a 
standalone program or a python module, a manual installation is prerequisite. 
An instruction for linux using bash follows in the subsequent section.

Two methods will be described. The installation via git and the installation 
using a tarball provided.

First choose and navigate the directory, where you want to install orbkit. 
In this example we will use the home directory. 

- Using git:
  Clone the repository and set a path variable for orbkit to that directory:

    $ cd $HOME
    $ git clone git://git.code.sf.net/p/orbkit/code orbkit
    $ export ORBKITPATH=$HOME/orbkit

- Using a tarball:
  Download the latest orbkit release, extract the file (`v0.2.0` has to be replaced
  by the version number), and set a path variable for orbkit to that directory:

    $ cd $HOME
    $ wget http://sourceforge.net/projects/orbkit/files/latest/download 
    $ tar xzvf orbkit.v0.2.0.tar.gz
    $ export ORBKITPATH=$HOME/orbkit

orbkit Modules
--------------

In order to use orbkit, you add the orbkit dirrectory to your PYTHONPATH 
enviroment variable either temporarily by typing

    $ export PYTHONPATH=$PYHONPATH:$ORBKITPATH

or permanently by adding this line to your ~/.bashrc file.

To use the orbkit standalone program you have to modify your PATH variable as
well:

    $ export PATH=$PATH:$ORBKITPATH/tools

Support
=======

If you need help for the usage of orbkit, do not hesitate to contact the 
orbkit support team via http://sourceforge.net/projects/orbkit/support