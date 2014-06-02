Function Reference
==================

Central Variables
-----------------

The following variables are the central variables in orbkit and are global
within :mod:`orbkit.main`.


.. data:: geo_info

  Contains information about the atoms and is only required for the creation 
  of the output file.
  
  - Type: list of str
      | 
  - Shape: (number of atoms, 3)
      | 
  - Members:
      :1st element: atom symbol
      :2nd element: atom number (according to the input file)
      :3rd element: nuclear charge of atom (if present in the input file)

.. data:: geo_spec

  Contains the atom positions in Bohr.
  
  - Type: numpy.ndarray or list of floats
      | 
  - Shape: (number of atoms, 3)
      | 
 
.. data:: ao_spec
  
  Contains all information about the atomic orbitals.
  
  - Type: list of dictionaries 
      - number of dictionaries corresponds to number of contracted Gaussians
  - Members:    
      :pnum: - Number of number of primitives (integer)
      :atom: - Atom number (according to the input file, i.e., starting from 1) (integer)
      :coeffs: - numpy.ndarray

	       - Shape: (number of primitives, 2)

	       :1st element: exponent
	       :2nd element: coefficient



.. data:: mo_spec

  Contains all information about the molecular orbitals.

  - Type: list of python dictionaries 
	- number of dictionaries corresponds to number of molecular orbitals
  - Members: 
      :energy: - Energy of molecular orbital (float)
      :occ_num: - Occupation of molecular orbital (float)
      :sym: - Symmetry of molecular orbital (string, MOLPRO-like: e.g., 12.1 or 12.A1)
      :coeffs: - Molecular orbital coefficients (numpy.ndarray) 
	       - Shape: (number of atomic orbitals, )
  

.. data:: (delta_)ao_list

  numpy.ndarray containing (the derivatives of) all atomic orbitals on the 
  specified grid.

.. data:: (delta_)mo_list

  numpy.ndarray containing (the derivatives of) all molecular orbitals on the 
  specified grid.
  
.. data:: (delta_)rho

  numpy.ndarray containing (the derivatives of) the electron density on the 
  specified grid.


orbkit.options
--------------
.. automodule:: orbkit.options 

Input/Output Options
....................
.. autodata:: orbkit.options.filename
.. autodata:: orbkit.options.itype
.. autodata:: orbkit.options.outputname
.. autodata:: orbkit.options.otype

Grid-Related Options
....................

.. autodata:: orbkit.options.vector
.. autodata:: orbkit.options.grid_file
.. autodata:: orbkit.options.center_grid
.. autodata:: orbkit.options.random_grid

Computational Options
.....................

.. autodata:: orbkit.options.numproc
.. autodata:: orbkit.options.mo_set
.. autodata:: orbkit.options.calc_mo
.. autodata:: orbkit.options.all_mo
.. autodata:: orbkit.options.drv

Additional Options
..................

.. autodata:: orbkit.options.z_reduced_density
.. autodata:: orbkit.options.atom_projected_density
.. autodata:: orbkit.options.mo_tefd

Options for Advanced Users
..........................

.. autodata:: orbkit.options.quiet
.. autodata:: orbkit.options.no_log
.. autodata:: orbkit.options.no_output
.. autodata:: orbkit.options.no_slice
.. autodata:: orbkit.options.interactive

Options for developers
......................

.. autofunction:: orbkit.options.get_options
.. autofunction:: orbkit.options.check_options
.. autofunction:: orbkit.options.check_if_exists

Parameters
..........

.. autodata:: orbkit.options.itypes
.. autodata:: orbkit.options.otypes
.. autodata:: orbkit.options.drv_options

orbkit.main
-----------
.. automodule:: orbkit.main 
    :members:

orbkit.read
------------
.. automodule:: orbkit.read 
    :members:

orbkit.grid
-----------
.. automodule:: orbkit.grid 
    :members:

orbkit.core
-----------
.. automodule:: orbkit.core 
    :members:

orbkit.display
--------------
.. automodule:: orbkit.display 
    :members:

orbkit.output
-------------
.. automodule:: orbkit.output 

.. autofunction:: orbkit.output.main_output
.. autofunction:: orbkit.output.HDF5_creator
.. autofunction:: orbkit.output.amira_creator
.. autofunction:: orbkit.output.hx_network_creator
.. autofunction:: orbkit.output.cube_creator

orbkit.extras
-------------
.. automodule:: orbkit.extras 
    :members:
 
