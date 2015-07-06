.. _installation-instructions:

Installation Instructions for orbkit
====================================

For a proper execution of orbkit, the following Python modules are
required:

1) Python_ 2.6 - 2.7
2) SciPy_ Library of algorithms and mathematical tools
3) NumPy_ Library of high-level mathematical functions
4) h5py_ Interface to the HDF5 binary data format
5) mayavi_ Tool for 3D scientific data visualization (optional)

The package h5py is not mandatory but strongly recommended.

.. _Python: http://www.python.org
.. _SciPy: http://www.scipy.org/
.. _NumPy: http://www.numpy.org/
.. _h5py: http://www.h5py.org/
.. _mayavi: http://docs.enthought.com/mayavi/mayavi/index.html

Installation
------------

orbkit has to be manually installed. This is a simple procedure and can 
be carried out in Linux using ``bash`` as follows:

Choose the directory, where you want to install orbkit. Open a terminal window, 
e.g. ``gnome-terminal``, and navigate to this directory. In this example we 
will use the home directory. If you use a different directory simply replace 
``$HOME`` by your preferred folder throughout the whole section.

    .. code-block:: bash

        $ cd $HOME

Get a copy of orbkit, either with `git`_ or using a `tarball`_. It is strongly
recommended to use `git`_, since this version always contains the newest 
bug fixes and features. If git is not available on your system, the newest 
version can additionally be cloned from http://sourceforge.net/p/orbkit/code.

  .. _git:

  * Using **git**:

    Clone the repository:

    .. code-block:: bash

        $ git clone http://git.code.sf.net/p/orbkit/code orbkit

  .. _tarball:

  * Using a **tarball**:

    Download the latest orbkit release and extract the file (``v0.2.0`` has to be 
    replaced by your version number):

    .. code-block:: bash

        $ wget http://sourceforge.net/projects/orbkit/files/latest/download 
        $ tar xzvf orbkit.v0.2.0.tar.gz

In order to use orbkit, you have to add the orbkit directory to your ``PYTHONPATH``
environment variable either temporarily by typing

.. code-block:: bash

    $ export ORBKITPATH=$HOME/orbkit
    $ export PYTHONPATH=$PYHONPATH:$ORBKITPATH

or permanently by adding these lines to your ~/.bashrc file.

To use orbkit as a standalone program, you have to modify your 
``PATH`` variable in the same way:

.. code-block:: bash

    $ export ORBKITPATH=$HOME/orbkit
    $ export PATH=$PATH:$ORBKITPATH/tools
