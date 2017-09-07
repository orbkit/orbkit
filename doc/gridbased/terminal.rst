.. _`terminal_interface`:

Usage via the Terminal
======================

For an overview of all available features, a general help text can be
called in the terminal by using the command:

.. code-block:: bash

    $ orbkit -h

.. note::

  ORBKIT creates a logfile (:literal:`.oklog`) containing all information printed
  during the calculation.


.. contents:: Table of Contents:
  :local:
  :depth: 1

.. _io:

Input/Output
------------

For the standard density calculation including all occupied molecular 
orbitals on the default grid, the following command must be entered in the
console: 

.. code-block:: bash

    $ orbkit -i INPUT 

The only option required by ORBKIT is the definition of the input file with 
:literal:`-i INPUT`. Unless otherwise stated, ORBKIT will attempt to automatically
determine the filetype of the input file and use HDF5 file as output type (:literal:`--otype=h5`). 
As output name (:literal:`-o OUTPUTNAME`) it assumes 
the base name of :literal:`INPUT`. The file extension of the output is 
automatically added by the program. In order to produce different output types 
within one run, multiple calls of the option :literal:`--otype=OTYPE` are possible.

The available input file types (``--itype=ITYPE``) are: **auto** (default), **molden** (Molden file), 
**aomix** (AOMix file), **gamess** (GAMESS-US output file), **gaussian.log** 
(GAUSSIAN output file), **gaussian.fchk** (GAUSSIAN formatted checkpoint file),
**wfn** and **wfx** files. 

You can choose between standard Gaussian_ cube files (:literal:`--otype=cb`), 
HDF5_ files (:literal:`--otype=h5`), or ZIBAmira_ Mesh files (:literal:`--otype=am`).
Additionally, ORBKIT can directly create ready to use VMD_ scripts 
(:literal:`--otype=vmd`) and ZIBAmira_ (:literal:`--otype=hx`) networks for the 
easy depiction in VMD_ and ZIBAmira_, respectively. 
While a VMD_ script is a plain text script based on Gaussian_ cube files, the  
ZIBAmira_ network uses ZIBAmira_ Mesh files and a simple ZIBAmira_ colormap file 
as input data.

The Gaussian_ cube file (file extension: .cb) is a plain text file containing
the important informations of the calculation. For large systems,
the cube file requires a lot of space on the hard drive and is very 
time-consuming.

In contrast, the HDF5_ file format (file extension: .h5) which is a hierarchical 
data format can store and organize large amounts of numerical data. More 
informations about this file format can be found at:

  http://en.wikipedia.org/wiki/Hierarchical_Data_Format

Another advantage of this file format is the easy readability with Matlab_, Python
and many other languages. 
Example files for loading HDF5_ files with Matlab_ or Python are available in our 
program package (``examples/HDF5_Examples``). With the JAVA program, HDFVIEW_, 
the user can easily load and read the HDF5_ files. 

The last available plain text data format are the ZIBAmira_ Mesh files. Those 
can be directly opened in ZIBAmira_.

If you want to simply visualize the data created by ORBKIT, you can also use 
a simple Mayavi_ interface (:literal:`--otype=mayavi`). Here, no data will be
saved to disc, if no other output type is specified.

.. hint::

  If you want to use ORBKIT for visualization, you may want to read the 
  :doc:`../quick`.

.. _grid:

Grid Related Options
--------------------

There are two ways to specify the grid when calling ORBKIT via the Terminal. 
You can adapt the grid to the molecular geometry: 

.. code-block:: bash

        $ orbkit -i INPUT --adjust_grid=D X

Here, ORBKIT creates a grid with a grid spacing of ``X`` a\ :sub:`0` and the size
of the molecule plus ``D`` a\ :sub:`0` in each direction.

Alternatively you can also modify the grid via an external plain text file 
(:literal:`GRID_FILE`):

.. code-block:: bash

	$ orbkit -i INPUT --grid=GRID_FILE

This file can have two possible formats. It can be represented either by the boundary
conditions of an equidistant rectangular grid (**regular grid**) or by a list of 
data points (**vector grid**):


  +----------------------------------------------------------------+----------------------------------------------------------------+
  | **regular grid**                                               | **vector grid**                                                |
  +----------------------------------------------------------------+----------------------------------------------------------------+
  | .. literalinclude:: ../../examples/basic_examples/grid_reg.txt | .. literalinclude:: ../../examples/basic_examples/grid_vec.txt |
  |    :language: bash                                             |    :language: bash                                             |
  +----------------------------------------------------------------+----------------------------------------------------------------+

.. note:: A :literal:`#` at the beginning of a line implicates a comment line.

ORBKIT performs all computations internally on a **vector grid**. 
For this purpose it converts a **regular grid** beforehand to a **vector grid**.
After the computation the original grid is recreated.

By default, ORBKIT the 1-dimensional **vector grids** into 1-dimensional slices of equal length. 
The atomic orbitals, the molecular orbitals, and the density are calculated for 
each slice separately. At the end of the calculation, the data
is reassembled and stored in an output file. 
This enables an easy parallelization and requires a smaller amount of RAM.

The length of the 1-dimensional slices can be set with

.. code-block:: bash

    $ orbkit -i INPUT --slice_length=1e4

In the default setting, ORBKIT performs the density calculation by starting 
only one subprocess. The number of subprocesses, which are distributed over 
the existing CPUs, can be modified with the subsequent command:

.. code-block:: bash

    $ orbkit -i INPUT --numproc=4

.. _mo:

Molecular Orbital Selection
---------------------------

ORBKIT is capable of calculating a selected set of molecular orbitals. This set
can be specified either **inline** or by using an **external file**.

You can use the **MOLPRO-like nomenclature**, e.g., ``3.1`` for the third orbital 
in symmetry one, or you choose it by the 
**index within the input file** (counting from **one**!). 

.. hint:: 
  
  For Gaussian_ and Gamess-US_, the symmetry labels used are different, 
  e.g., ``3.A1`` for the third orbital in symmetry A1.

.. note:: 
  
  For unrestricted calculations, the symmetry labels are extended by ``_a`` 
  for alpha and by ``_b`` for beta molecular orbitals, e.g., ``3.A1_b``.

In the latter case, you can additionally use the keywords ``homo`` (highest occupied 
molecular orbital), ``lumo`` (lowest unoccupied molecular orbital),
as well as ``last_bound`` (highest-lying orbital with negative orbital energy) 
You can further select a range of orbitals, e.g., ``--calc_mo=1:homo-1``, which evokes the 
computation of the molecular orbitals 1, 2, 3, ..., and homo-2.

+-------------------+--------------------------------------------------------------+------------------------------------------------------------------+
|                   |  **MOLPRO-like Nomenclature**                                | **Index within the Input File**                                  |
+-------------------+--------------------------------------------------------------+------------------------------------------------------------------+
| **Inline**        |.. code-block:: bash                                          |.. code-block:: bash                                              |
|                   |                                                              |                                                                  |
|                   |    $ orbkit -i INP --calc_mo=1.1,1.3                         |    $ orbkit -i INP --calc_mo=3:lumo+3,1                          |
|                   |                                                              |                                                                  |
|                   |Hint: Multiple calls are possible.                            |Hint: Multiple calls are possible.                                |
+-------------------+--------------------------------------------------------------+------------------------------------------------------------------+
| **Ext. File**     |.. code-block:: bash                                          |.. code-block:: bash                                              |
|                   |                                                              |                                                                  |
|                   |    $ orbkit -i INP --calc_mo=MO_LIST                         |    $ orbkit -i INP --calc_mo=MO_LIST                             |
|                   |                                                              |                                                                  |
|                   |``MO_LIST``:                                                  |``MO_LIST``:                                                      |
|                   |                                                              |                                                                  |
|                   |.. literalinclude:: ../../examples/basic_examples/MO_List.tab |.. literalinclude:: ../../examples/basic_examples/MO_List_int.tab |
|                   |    :language: bash                                           |    :language: bash                                               |
|                   |                                                              |                                                                  |
+-------------------+--------------------------------------------------------------+------------------------------------------------------------------+

The computation and storage of all molecular orbitals can be called by 

.. code-block:: bash

    $ orbkit -i INPUT --calc_mo=all_mo

One special capability of ORBKIT is the computation of the density with a selected 
set of molecular orbitals. 

.. code-block:: bash

    $ orbkit -i INPUT --mo_set=MO_SET

The selection of molecular orbitals can be accomplished in the same manner as
described above for ``--calc_mo``. Although for ``--mo_set``, each line in the 
external file or each call of ``--mo_set`` corresponds to one density calculation.
    
Derivative Calculation
----------------------

ORBKIT can compute analytical spatial derivatives up to second order
with respect to :math:`x`, :math:`y`, or :math:`z` for the atomic and 
molecular orbitals, as well as for the electron density. 
For instance, a derivative of the density with 
respect to :math:`x` can be invoked as follows:

.. code-block:: bash

    $ orbkit -i INPUT --drv=x

Multiple calls of the option :literal:`--drv=DRV` are possible.

For second derivatives, you can specify the respective combinations, 
e.g., 'xx' or 'yz'.

The computation of the laplacian, i.e., 
:math:`\nabla^2 \rho = \nabla^2_x \rho + \nabla^2_y \rho + \nabla^2_z \rho`,
can be invoked by

.. code-block:: bash

    $ orbkit -i INPUT --laplacian

Spin-Density
------------

For unrestricted calculations, the spin density and related quantities 
(e.g. derivatives) may be calculated by 

.. code-block:: bash

    $ orbkit -i INPUT --spin=alpha

for alpha spin density and by

.. code-block:: bash

    $ orbkit -i INPUT --spin=beta

for beta spin density. 

.. note ::
  The usage of the ``--spin=`` keyword omits the reading
  of the molecular orbitals of the other spin.

Additional Options
------------------

In the following, two additional features are highlighted. 

Firstly, the atom-projected electron density can be computed by

.. code-block:: bash

    $ orbkit -i INPUT --atom_projected_density=INDEX

which is the integrand of the Mulliken charges, and secondly, ORBKIT 
is capable of calculating the molecular orbital transition electronic flux density 
(components :literal:`x`, :literal:`y`, and :literal:`z`) between the orbitals 
:literal:`I` and :literal:`J`:

.. code-block:: bash

    $ orbkit -i INPUT --mo_tefd=I J --drv=x --drv=y --drv=z

In order to compute and store **all** atomic orbitals or a derivative thereof,
you can call

.. code-block:: bash

    $ orbkit -i INPUT --calc_ao

Here, it is only possible to compute one derivative at a time.

.. _HDF5: http://www.hdfgroup.org/HDF5/
.. _HDFVIEW: http://www.hdfgroup.org/products/java/hdf-java-html/hdfview/
.. _MOLPRO: https://www.molpro.net/
.. _TURBOMOLE: http://www.turbomole.com/
.. _Gamess-US: http://www.msg.chem.iastate.edu/gamess/
.. _Gaussian: http://www.gaussian.com/
.. _ZIBAmira: http://amira.zib.de/
.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _Matlab: http://www.mathworks.de/products/matlab/
.. _Mayavi: http://docs.enthought.com/mayavi/mayavi/index.html
