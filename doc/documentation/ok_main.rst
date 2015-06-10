orbkit's High-Level Interface
=============================

For an overview of all available features, refer to the chapter :ref:`options`.

The following examples show exemplary how to use orbkit within your Python 
programs. These and more examples can be found in the :literal:`orbkit/examples` 
folder. 

Input/Output
------------

For the standard density calculation including all occupied molecular 
orbitals on the default grid, the following command must be entered in the
console::

  import orbkit as ok
  ok.options.filename   = 'h2o.md'                # input file name
  ok.options.itype      = 'molden'                # input file type [default]
  ok.options.outputname = 'h2o'                   # output file (base) name
  ok.options.otype      = ['h5','vmd']            # output file type [default]
  ok.options.numproc    = 4                       # number of processes
  # run orbkit
  ok.run_orbkit()

.. hint:: 

  If you want to run several calculations with different options, you have to
  reset all options. This can be accomplished by calling :mod:`orbkit.main.init`
  between the calculations, e.g.::

    import orbkit as ok
    # Set some options ...
    # run orbkit
    ok.run_orbkit()
    
    # Reset all options  
    ok.main.init()          
    # Set some options ...   
    # run orbkit
    ok.run_orbkit()                   

Grid Related Options
--------------------

::

  ok.options.adjust_grid= [5, 0.1]                # adjust the grid to the geometry

::

  ok.options.grid_file  = 'grid_reg.txt'          # grid file to read from (regular grid)
  
::
  
  ok.options.grid_file  = 'grid_vec.txt'          # grid file to read from (vector grid)
  ok.options.vector     = 1e4                     # number of points per subprocess

::

  ok.grid.N_            = [  201,   201,   101]   # grid points (regular grid)
  ok.grid.max_          = [ 10.0,  10.0,   5.0]   # maximum grid value
  ok.grid.min_          = [-10.0, -10.0,  -5.0]   # minimum grid value

Molecular Orbital Selection
---------------------------

::

  ok.options.calc_mo    = ['2.1','1.1','2.3']     # list of MO labels (Molpro notation)

::

  ok.options.mo_set     = [[1,2,3],['homo', 
                           'lumo+2:lumo+4']]      # list of MO labels (Molden enumeration)

Derivative Calculation
----------------------

::

  ok.options.drv        = ['x', 'z']              # derivatives along x and z

.. Additional Options
.. ------------------

.. .. literalinclude:: ../../examples/use_orbkit_function.py
..    :language: python
