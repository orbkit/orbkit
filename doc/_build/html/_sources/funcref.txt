Function Reference
==================

Central Variables
-----------------

The following variables are the central variables in orbkit. When/if 
they are created, they are global variables within :mod:`orbkit.main`


.. data:: geo_info

  Contains information about the atoms and is only required for the creation 
  of the output file.
  
  - Type: list of str
  - Shape: (number of atoms, 3)
  - Members:
  
    #. element: atom symbol
    #. element: atom number (according to the input file)
    #. element: nuclear charge of atom (if present in the input file)

.. data:: geo_spec

  Contains the atom positions in Bohr.
  
  - Type: numpy.ndarray or list of floats
  - Shape: (number of atoms, 3)
    
.. data:: ao_spec
  
  Contains all information on the atomic orbitals.
  
  - Type: list of dictionaries 
	  (number of dictionaries corresponds to number of contracted Gaussians)
  - Members:    
      :pnum: Number of number of primitives (integer)
      :atom: Atom number (according to the input file, i.e., starting from 1) (integer)
      :coeffs: numpy.ndarray

	Shape: (number of primitives, 2)

	#. element: coefficient
	#. element: exponent

.. data:: mo_spec

  Contains all information on the atomic orbitals.

  - Type: list of python dictionaries 
	  (number of dictionaries corresponds to number of molecular orbitals)
  - Members: 
      :energy: Energy of molecular orbital (float)
      :occ_num: Occupation of molecular orbital (float)
      :sym: Symmetry of molecular orbital (string, MOLPRO-like: e.g., 12.1 or 12.A1)
      :coeffs: Molecular orbital coefficients (numpy.ndarray) 

	Shape: (number of atomic orbitals, )
  

.. data:: (delta_)ao_list

  numpy.ndarray containing (the derivatives of) all atomic orbitals on the 
  specified grid.

.. data:: (delta_)mo_list

  numpy.ndarray containing (the derivatives of) all molecular orbitals on the 
  specified grid.
  
.. data:: (delta_)rho

  numpy.ndarray containing (the derivatives of) the electron density on the 
  specified grid.



orbkit.main
-----------
.. automodule:: orbkit.main 
    :members:

orbkit.core
-----------
.. automodule:: orbkit.core 
    :members:

orbkit.read
------------
.. automodule:: orbkit.read 
    :members:

orbkit.output
-------------
.. automodule:: orbkit.output 
    :members:

orbkit.extras
-------------
.. automodule:: orbkit.extras 
    :members:
 
