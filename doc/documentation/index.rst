Documentation
=============

.. orbkit can be easily operated via the terminal or within your own Python programs.
    For advanced users, the latter type of usage is recommended.

This documentation should serve as an overview of how to operate orbkit. 

First, we describe how to prepare inputs from several quantum chemistry programs
(:doc:`requirements`). Then, we show how to use orbkit as standalone program in 
the Terminal (:doc:`terminal`).

Subsequently, we show how to use those functionalities in your Python programs
with orbkit's :doc:`ok_main`. You can also use orbkit's functions directly 
in your Python programs. This is shown in :doc:`ok_modules`.

Besides computing quantities using arbitrary grids, orbkit is also able to
compute analytical integrals between atomic and molecular orbitals. Thus,
you can use orbkit to compute expectation values analytically. The basics are
described in :doc:`ok_analytical`.

Finally, we provide some chapters on orbkit's :doc:`cvars`  and :doc:`options`.
The chapter :doc:`funcref` contains the docstrings of all functions in orbkit.

.. hint::
  
  Please refer to :doc:`../tutorials/index` for practical applications.

**Contents:**

  .. toctree::
    :maxdepth: 4

    requirements
    terminal
    ok_main
    ok_modules
    ok_analytical
    cvars
    options
    funcref

.. note::
  
  The documentation is under construction. If you have any questions,
  do not hesitate to contact us.