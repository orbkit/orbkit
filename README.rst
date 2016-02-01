orbkit
======

orbkit is a parallel Python program package for post-processing wavefunction 
data extracted from output files of MOLPRO (Molden File Format), TURBOMOLE 
(AOMix file format), GAMESS US, and Gaussian (Formatted Checkpoint File). 

Computational capabilities:

- Atomic orbitals, molecular orbitals, electron density, and the spatial derivatives of these quantities
- Reduced electron density
- Electron density for selected molecular orbitals
- Mulliken and Löwdin population analysis
- Analytical overlap integrals between atomic and molecular orbitals
- Atom-projected electron density
- Molecular orbital transition electronic flux density

If you use orbkit in your work, please cite it as follows:

Gunter Hermann, Vincent Pohl, Axel Schild: orbkit: A Toolbox for Post-Processing
Quantum Chemical Wavefunction Data, available via 
http://sourceforge.net/p/orbkit (2016).

Website: http://orbkit.github.io

Copyright (C) 2016, Gunter Hermann, Vincent Pohl, and Axel Schild.

Licence Note
============

orbkit is free software: you can redistribute it and/or modify it under the 
terms of the GNU Lesser General Public License as published by the Free Software 
Foundation, either version 3 of the License, or any later version.

orbkit is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along 
with orbkit. If not, see <http://www.gnu.org/licenses/>.

Installation Requirements
=========================

For a proper execution of orbkit, the following Python modules are required:

1) Python 2.6 - 2.7 (http://www.python.org) 
   Python 3.x probably works, too, but is officially not (yet) supported
2) SciPy Library of algorithms and mathematical tools (http://www.scipy.org/)
3) Weave (scipy.weave) for parallel computing (orbkit will not work without this)
4) NumPy Library of high-level mathematical functions (http://www.numpy.org/)
5) h5py Interface to the HDF5 binary data format (http://www.h5py.org/)
6) Mayavi Tool for 3D scientific data visualization (optional)

The package h5py is not mandatory but strongly recommended.

Installation
============

orbkit has to be manually installed. This is a simple procedure and can 
be carried out in Linux using bash as follows:

Choose the directory, where you want to install orbkit. Open a terminal window, 
e.g. gnome-terminal, and navigate to this directory. In this example we 
will use the home directory. If you use a different directory simply replace 
$HOME by your preferred folder throughout the whole section.

    $ cd $HOME

Get a copy of orbkit, either with git or using a tarball. It is strongly
recommended to use git, since this version always contains the newest 
bug fixes and features. If git is not available on your system, the newest 
version can additionally be cloned from https://github.com/orbkit/orbkit.

  * Using git:

    Clone the repository:

        $ git clone https://github.com/orbkit/orbkit.git

  * Using a tarball:

    Download the latest orbkit release and extract the file ("0.4.0" has to be 
    replaced by your version number):

        $ wget http://sourceforge.net/projects/orbkit/files/latest/download 
        $ tar xzvf orbkit.v0.4.0.tar.gz

In order to use orbkit, you have to add the orbkit directory to your $PYTHONPATH
environment variable either temporarily by typing

    $ export ORBKITPATH=$HOME/orbkit
    $ export PYTHONPATH=$PYHONPATH:$ORBKITPATH

or permanently by adding these lines to your ~/.bashrc file.

To use orbkit as a standalone program, you have to modify your 
$PATH variable in the same way:

    $ export ORBKITPATH=$HOME/orbkit
    $ export PATH=$PATH:$ORBKITPATH/tools

Documentation
=============

orbkit's documentation may be found at http://orbkit.github.io

Support
=======

If you need help for the usage of orbkit, contact the orbkit support team via 
https://github.com/orbkit/orbkit/issues
