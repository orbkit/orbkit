.. _`Options`:

Options
=======

Whenever ORBKIT's main function (:mod:`orbkit.main.run_orbkit`) is used within a Python 
program, all actions are controlled via the ORBKIT module :mod:`orbkit.options`.
This module also controls the parser, which is executed when ORBKIT is called
via the Terminal.

To reset all options, simply reload the module or call :mod:`orbkit.main.init`::

  from orbkit import main, options
  # Set an option
  options.all_mo = True
  # Reset all options
  reload(options)
  
  # Set another option
  options.numproc = 4
  # Reset all options
  main.init()
  
The options itself are global variables within the :mod:`orbkit.options` module
and are listed and explained below.
A detailed description can also be found by calling:

.. code-block:: bash

    $ orbkit -h

Input/Output Options
....................

.. autodata:: orbkit.options.filename
.. autodata:: orbkit.options.itype
.. autodata:: orbkit.options.outputname
.. autodata:: orbkit.options.otype

Grid-Related Options
....................

.. autodata:: orbkit.options.vector
.. autodata:: orbkit.options.center_grid
.. autodata:: orbkit.options.grid_file
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

.. autodata:: orbkit.options.gross_atomic_density
.. autodata:: orbkit.options.mo_tefd

Options for Advanced Users
..........................

.. autodata:: orbkit.options.quiet
.. autodata:: orbkit.options.no_log
.. autodata:: orbkit.options.no_output
.. autodata:: orbkit.options.no_slice
.. autodata:: orbkit.options.interactive
