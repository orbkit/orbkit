.. _`High-Level Interface`:

ORBKIT's High-Level Interface
=============================

This chapter should serve as an overview of how to use ORBKIT's high-level 
interface within your Python programs. For clarity, it is structured 
equivalently to :ref:`terminal_interface`.

.. The following examples show exemplary how to use ORBKIT from within your Python 
   programs. These and more examples can be found in the :literal:`orbkit/examples` 
   folder. 


.. contents:: Table of Contents:
  :local:
  :depth: 1

General Aspects
---------------

To access the main features of ORBKIT, you have to import it using::

  import orbkit

Now, you can specify all options for the computation analogously to 
:ref:`terminal_interface`, e.g.::

  orbkit.options.numproc = 4                       # number of processes

Within this chapter we will discuss the most important options. 
For an overview of all available options, please refer to the chapter 
:ref:`options`.

Finally, you have to run ORBKIT with::

  data = orbkit.run_orbkit()

This function will not only save the results to the requested output formats,
it will also return the results for further post-processing.

.. hint:: 

  If you want to run several calculations with different options, you have to
  reset all options. This can be accomplished by calling ``orbkit.init()``
  between the calculations, e.g.::

    import orbkit as ok
    # Set some options ...
    # run orbkit
    data_1 = orbkit.run_orbkit()
    
    # Reset all options  
    orbkit.init()          
    # Set some options ...   
    # run orbkit
    data_2 = orbkit.run_orbkit()            

Input/Output
------------

The input filename can be specified by::

  orbkit.options.filename = 'h2o.molden'

The available file types are **'molden'**, **'aomix'** (AOMix file),
**'gamess'** (GAMESS-US output file), **'gaussian.log'** (GAUSSIAN output file), 
**'gaussian.fchk'** (GAUSSIAN formatted checkpoint file),
**'wfn'** and **'wfx'** files. 

Concerning ORBKIT's output, you can choose between the following options:

  - **'h5'** (HDF5 file)
  - **'cb'** (Gaussian cube file)
  - **'vmd'** (VMD network) 
  - **'am'** (ZIBAmiraMesh file)
  - **'hx'** (ZIBAmira network)
  - **'vmd'** (VMD network) 
  - **'mayavi'** (opens a simple interactive Mayavi interface) 

You can also specify more than one output format at the same time::

  orbkit.options.otype = ['h5','vmd']

ORBKIT assumes the input file, if not otherwise stated, e.g.::

  orbkit.options.outputname = 'h2o' # output file (base) name

For more information on the different output types, refer to 
:ref:`io` (:doc:`./terminal`).

.. hint::

  You can omit the creation of an output file by either setting 
  ``orbkit.options.otype = []`` or by setting 
  ``orbkit.options.no_output = True``. Furthermore, you can disable the creation of 
  .oklog file can with ``orbkit.options.no_log = True``, and the terminal output can 
  be disabled with ``orbkit.options.quiet = True``.

Grid Related Options
--------------------

Although the default setting for a grid in ORBKIT is a regular grid, i.e., 
:math:`N_{\sf data points} = N_{\sf x} \times N_{\sf y} \times N_{\sf z}`,
ORBKIT carries out all computations on a a vector grid, i.e.,
:math:`N_{\sf data points} = N_{\sf x} = N_{\sf y} = N_{\sf z}`.
Thus for the regular case, the grid and all output data is converted 
automatically back and forth within the computational procedures.

To omit the back transformation to a regular grid, you may set::

  options.is_vector = True

For the main computational tasks, the grid is divided into slices,
the length of which, i.e., the number of grid of points per 
subprocess, can be set by::
 
  orbkit.options.slice_length = 1e4 

.. note::
  
  If you want to omit the slicing, please set::
    
    orbkit.options.no_slice = True
    

There are several ways to specify the grid in ORBKIT (in a.u.):
  
**Adjusting the grid to the geometry**::

  orbkit.options.adjust_grid = [5, 0.1]
  
Here, ORBKIT creates a grid with a grid spacing of 0.1 a\ :sub:`0` and the size
of the molecule plus 5 a\ :sub:`0` in each direction.

**Reading the grid parameters from a file**::

  orbkit.options.grid_file  = 'grid.txt'

This file can have two possible formats. It can be represented either by the boundary
conditions of an equidistant rectangular grid (**regular grid**) or by a list of 
grid points (**vector grid**). For more information, refer to
:ref:`grid` (:doc:`./terminal`).

**Specifying the boundary conditions manually**::

  orbkit.grid.N_   = [  201,   201,   101]   # grid points (regular grid)
  orbkit.grid.max_ = [ 10.0,  10.0,   5.0]   # maximum grid value
  orbkit.grid.min_ = [-10.0, -10.0,  -5.0]   # minimum grid value

.. note::

  If you prefer setting the grid spacing instead of the number of data points,
  you may set this parameters by::
    
    orbkit.grid.delta_ = [0.1, 0.2, 0.1]
  
**Specifying the grid manually**::

  import numpy
  orbkit.grid.x = numpy.linspace(-10,10,201)  
  orbkit.grid.y = numpy.array([0],dtype=float)   
  orbkit.grid.z = numpy.array([-1.0,1.1])    
  # We have already initialized a grid for orbkit:
  orbkit.grid.is_initialized = True

Here, x, y and z have to be one-dimensional ``numpy.array`` of type ``float``
(``numpy.float64``). 

.. attention::
  
  The last line is **mandatory**, i.e., we have to tell ORBKIT, that there is no 
  need to initialize the grid.

.. hint ::

  Please keep in mind that for a  **vector grid** the relation
  :math:`N_{\sf data points} = N_{\sf x} = N_{\sf y} = N_{\sf z}`
  has to hold.
  
  If you have initialized a **vector grid** manually, do not forget to  
  also set the variable ``grid.is_vector = True``. You can use this 
  standard variable as an input parameter in other ORBKIT functions.


.. _`mo high-level`:

Molecular Orbital Selection
---------------------------

ORBKIT is capable of calculating a selected set of molecular orbitals::

  orbkit.options.calc_mo = ['3.1','1.1','2.3']

and of calculating the density with a selected set of molecular orbitals::

  orbkit.options.mo_set = [[1,2,3],                   # first set
		       ['homo', 'lumo+2:lumo+4']] # second set

.. note::
  
  While the first example uses the **MOLPRO-like nomenclature**, e.g., ``3.1`` for 
  the third orbital in symmetry one, the second example uses the 
  **index within the input file** (counting from one). 

  For unrestricted calculations, the symmetry labels are extended by ``_a`` 
  for alpha and by ``_b`` for beta molecular orbitals, e.g., ``3.A1_b``.  

  For more information, refer to :ref:`mo` (Usage via the Terminal).

Derivative Calculation
----------------------

ORBKIT can compute analytical spatial derivatives with respect to :math:`x`,
:math:`y`, and :math:`z` for the atomic and molecular orbitals, as well
as for the electron density::

  orbkit.options.drv = ['x', 'z']

This invokes the computation of the derivatives with respect to :math:`x`
and the computation of the derivatives with respect to :math:`z`. 
For second derivatives, specify the respective combinations,e.g., 'xx' or 'yz'.

Spin-Density
------------

For unrestricted calculations, the spin density and related quantities 
(e.g. derivatives) may be calculated by::

  orbkit.options.spin = 'alpha'

or::

  orbkit.options.spin = 'beta'

The usage of this keyword omits the reading of the molecular orbitals of the other 
spin.

Return Values
-------------

Besides writing the requested output, the function ``run_orbkit()``,
returns all data computed::

  data = orbkit.run_orbkit()

Depending on your options, this data set has a different structure.

+---------------------------------+-------------------------------------------------------------------+
|**Computed Quantity**            | **Returned Data**                                                 |
+---------------------------------+-------------------------------------------------------------------+
|density                          | ``numpy.ndarray`` with ``shape=(N)``                              |
+---------------------------------+-------------------------------------------------------------------+
|derivative of density            |1. density ``numpy.ndarray`` with ``shape=(N)``                    |
|                                 |2. derivative of density ``numpy.ndarray`` with ``shape=(NDRV,N)`` |
+---------------------------------+-------------------------------------------------------------------+
|molecular orbitals               |1. ``numpy.ndarray`` with ``shape=((NMO,) + N)``                   |
|                                 |2. ``dict`` with information on selected molecular orbitals        |
+---------------------------------+-------------------------------------------------------------------+
|derivative of molecular orbitals |1. ``numpy.ndarray`` with ``shape=((NDRV,NMO,) + N)``              |
|                                 |2. ``dict`` with information on selected molecular orbitals        |
+---------------------------------+-------------------------------------------------------------------+
|density from a set of            |1. ``numpy.ndarray`` with ``shape=((NSET,) + N)``                  |
|                                 |                                                                   |
|molecular orbitals               |2. ``dict`` with information on selected molecular orbitals        |
+---------------------------------+-------------------------------------------------------------------+
|derivative of density from a     |1. ``numpy.ndarray`` with ``shape=((NSET,NDRV,) + N)``             |
|                                 |                                                                   |
|set of molecular orbitals        |2. ``dict`` with information on selected molecular orbitals        |
+---------------------------------+-------------------------------------------------------------------+

- ``N`` is shape as the grid.
- ``NDRV`` is the number derivatives requested.
- ``NMO`` is the number of molecular orbitals requested.
- ``NSET`` is the number of molecular orbital sets requested.

