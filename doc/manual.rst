Manual
======

orbkit can be easily operated via the terminal or within a Python program.  For 
advanced users, the latter type of usage is recommended.

Requirements for the Input Files
--------------------------------

orbkit requires Cartesian Harmonic Gaussian basis sets. Unless otherwise 
stated, it assumes the standard molden basis function order for the exponents
:math:`(l_x,l_y,l_z)`:

- s:
  (0,0,0)
- p:
  (1,0,0), (0,1,0), (0,0,1)
- d: 
  (2,0,0), (0,2,0), (0,0,2),
  (1,1,0), (1,0,1), (0,1,1)
- f: 
  (3,0,0), (0,3,0), (0,0,3), 
  (1,2,0), (2,1,0), (2,0,1), 
  (1,0,2), (0,1,2), (0,2,1), 
  (1,1,1)
- g: 
  (4,0,0), (0,4,0), (0,0,4), 
  (3,1,0), (3,0,1), (1,3,0), 
  (0,3,1), (1,0,3), (0,1,3), 
  (2,2,0), (2,0,2), (0,2,2), 
  (2,1,1), (1,2,1), (1,1,2)

**Molden File Format:**

- Contains Cartesian Harmonics by default
- Starts with :literal:`[Molden Format]`
- Contains :literal:`[Atoms]`, :literal:`[GTO]`, :literal:`[MO]`
- Generation:

  1) MOLPRO: http://www.molpro.net/info/2012.1/doc/manual/node102.html
  2) TURBOMOLE: tm2molden

**GAMESS-US Output File:**

- Please use Cartesian Harmonics (default)

**GAUSSIAN Formatted Checkpoint File:**

- Contains Cartesian Harmonics by default
- Not applicable for natural orbitals => occupation numbers are not printed
- Creation of a FChk file:

  1) Add :literal:`%Chk=chkpt-file` to your Gaussian_ input file
  2) Use :literal:`formchk` to convert the chk file:

.. code-block:: bash

    $ formchk chkpt-file formatted-file

**GAUSSIAN .log File:**

- Spherical Harmonics are chosen by default! You have to switch manually 
  to Cartesian Harmonics (:literal:`6D 10F`)!
- Use the following parameters in your root section

.. code-block:: bash

    6D 10F gfinput IOP(6/7=3)

- If more than one geometry/atomic orbitals/molecular orbitals sections are present
  in the .log file, orbkit provides an interactive selection.

Usage via the Terminal
----------------------

For an overview of all available features, a general help text can be
called in the terminal by using the command:

.. code-block:: bash

    $ orbkit -h

orbkit creates an logfile (:literal:`.oklog`) containing all information printed
during the calculation.

Input/Output
............
For the standard density calculation including all occupied molecular 
orbitals on the default grid, the following command must be entered in the
console: 

.. code-block:: bash

    $ orbkit -i INPUT --itype=ITYPE -o OUTPUTNAME --otype=OTYPE

The only option required by orbkit is the definition of the input file with 
:literal:`-i INPUT`. Unless otherwise stated, orbkit assumes 
:literal:`--itype=molden` and :literal:`--otype=h5`. For the 
:literal:`-o OUTPUTNAME`, it assumes the base name of :literal:`INPUT`. 
The file extension of the output will be automatically extended by the program. 
In order to produce different output types within one run, multiple calls 
of the option :literal:`--otype=OTYPE` are possible. You can choose between
standard Gaussian_ cube files (:literal:`--otype=cb`), 
HDF5_ files (:literal:`--otype=h5`) or ZIBAmira_ Mesh files (:literal:`--otype=am`).
Additionally, orbkit can directly create a ready to use ZIBAmira_ network for 
the easy depiction in ZIBAmira_ (:literal:`--otype=hx`). This network is 
based on a ZIBAmira_ Mesh file and includes a simple color map adapted to the
system.

The Gaussian_ cube file (file extension: .cb) is a normal text file containing
the important informations of the calculation. For large systems,
the cube file requires a lot of space on the hard drive and is very 
time-consuming.

In contrast, the HDF5_ file format (file extension: .h5) which is a hierarchical 
data format can store and organize large amounts of numerical data. More 
informations about this file format can be found at:

  http://en.wikipedia.org/wiki/Hierarchical_Data_Format

Another advantage of this file format is the easy readability with Matlab_, Python
and many other languages. 
Example files for loading HDF5_ files via Matlab_ or Python are available in the 
program package. With the JAVA program, HDFVIEW_, the user can easily load and 
read the HDF5_ files. 

The last available data format are the ZIBAmira_ Files. Those can be directly
loaded in ZIBAmira_.

Grid Related Options
....................

For the modification of the grid, one can specify a grid via an external 
plain text file (:code:`GRID_FILE`):

.. code-block:: bash

	$ orbkit -i INPUT --grid=GRID_FILE

This file can have two possible formats. (A :literal:`#` at the beginning of 
a line implicate a comment line) It can be represented either by the boundary
conditions of an equidistant rectangular grid (`regular grid`),

.. literalinclude:: ../examples/grid_reg.txt
   :language: bash

or by a list of data points
(`vector grid`),

.. literalinclude:: ../examples/grid_vec.txt
   :language: bash

By default, orbkit divides 3-dimensional `regular grids` into 2-dimensional 
slices or 1-dimensional `vector grids` into 1-dimensional slices of equal length. 
The atomic orbitals, the molecular orbitals, and the density are calculated for 
each slice separately. At the end of the calculation, the data
is reassembled and stored in an output file. 

For `vector grids`, the length of the 1-dimensional slices can be defined with

.. code-block:: bash

    $ orbkit -i INPUT --vector=1e4

In the default setting, orbkit 
performs the density calculation by starting only one subprocess.
The number of subprocesses, which are distributed over the existing CPUs, 
can be modified with the subsequent command:

.. code-block:: bash

    $ orbkit -i INPUT --numproc=4


Molecular Orbital Selection
...........................

One special capability of orbkit is the computation of the density with a selected 
set of molecular orbitals. This can be invoked by providing an :literal:`MO_SET`:

.. code-block:: bash

    $ orbkit -i INPUT --mo_set=MO_SET
    
Such an :literal:`MO_SET` is a plain text file, where every row signifies a new 
calculation of the density from the molecular orbitals specified in this row.

The molecular orbitals can be selected either via their index (As given in the 
output file, i.e., starting from one.)

.. literalinclude:: ../examples/MO_List_int.tab
   :language: bash

or by using the MOLPRO_-like nomenclature, 
e.g., 3.1 for the third orbital in symmetry one. 
(For Gaussian_ and Gamess-US_, the symmetry labels are used, e.g., 3.A1 for the third orbital in 
symmetry A1.) 

.. literalinclude:: ../examples/MO_List.tab
   :language: bash


In the same manner, the computation and storage of a selected set of molecular
orbitals can be accomplished by

.. code-block:: bash

    $ orbkit -i INPUT --calc_mo=CALC_MO
    
Here, the density computation is omitted. The computation and storage of all
molecular orbitals can be called by 

.. code-block:: bash

    $ orbkit -i INPUT --calc_mo=all_mo
    
Derivative Calculation
......................

orbkit can compute analytical spatial derivatives with respect to :literal:`x`,
:literal:`y`, and/or :literal:`z` for the atomic and molecular orbitals, as well
as for the electron density. E.g., a derivative of the density with respect to 
:literal:`x` can be invoked as follows:

.. code-block:: bash

    $ orbkit -i INPUT --drv=x

Multiple calls of the option :literal:`--drv=DRV` are possible.

Additional Options
..................

In the following, two additional features are highlighted. 
On the one hand, the atom-projected electron density 

.. code-block:: bash

    $ orbkit -i INPUT --atom_projected_density=INDEX

which is the integrand of the Mulliken charges, and on the other hand, the 
molecular orbital transition electronic flux density (components :literal:`x`,
:literal:`y`, and :literal:`z`) 
between the orbitals :literal:`I` and :literal:`J`:

.. code-block:: bash

    $ orbkit -i INPUT --mo_tefd=I J --drv=x --drv=y --drv=z


.. _HDF5: http://www.hdfgroup.org/HDF5/
.. _HDFVIEW: http://www.hdfgroup.org/products/java/hdf-java-html/hdfview/
.. _MOLPRO: https://www.molpro.net/
.. _TURBOMOLE: http://www.turbomole.com/
.. _Gamess-US: http://www.msg.chem.iastate.edu/gamess/
.. _Gaussian: http://www.gaussian.com/
.. _ZIBAmira: http://amira.zib.de/
.. _Matlab: http://www.mathworks.de/products/matlab/

Usage Within a Python Program
-----------------------------

The following examples show exemplary how to use orbkit within your Python 
programs. These and more examples can be found in the :literal:`orbkit/examples` 
folder. Please refer to the function references to get information about
all modules and functions available.

Using orbkit's main Function Interface
......................................

.. literalinclude:: ../examples/use_orbkit_function.py
   :language: python

Using orbkit's as a Module
..........................

.. literalinclude:: ../examples/use_as_module.py
   :language: python

Using orbkit to Calculate Analytical Derivatives of Molecular Orbitals
......................................................................

.. literalinclude:: ../examples/calculate_derivatives.py
   :language: python



Limitations of orbkit
---------------------

One limitation of orbkit is the ability to calculate only s, p, d, f and g atomic 
orbitals (Molden file limitation).

For the calculation of atomic and molecular orbitals, orbkit requires 
Cartesian Gaussians. The usage of spherical Gaussian functions is not yet implemented. 