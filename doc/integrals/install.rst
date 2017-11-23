.. _`Install libcint`:

Installation of libcint
=======================

Libcint needs to be installed separately, it is distributed under the BSD-2-Clause and available at GitHub.

First, download the source,

.. code-block:: bash

    $ git clone http://github.com/sunqm/libcint.git

then build build it

.. code-block:: bash

    $ cd libcint
    $ mkdir build; cd build
    $ cmake -DCMAKE_INSTALL_PREFIX:PATH=$HOME/libcint/ ..
    $ make install

and finally, add the library to the PATH environment variable by adding:

.. code-block:: bash

    $ export PATH=$PATH:$HOME/libcint

to your `~/.bashrc`.

.. hint::

    There is also an X86 optimized version of libcint, which can be installed in a similar way: https://github.com/sunqm/qcint
