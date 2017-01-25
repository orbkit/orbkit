.. image:: detci_orbkit.png
   :height: 200 px
   :width: 200 px
   :alt: detCI@ORBKIT Logo
   :align: right
   
.. _`detCI\@ORBKIT`:

detCI\@ORBKIT
=============

detCI\@ORBKIT is an extension of ORBKIT for post-processing **multi-determinantal**
wave function data from the output of quantum chemical calculations.
Our procedure, allows to apply different one-electron operators to Configuration Interaction (CI) wave functions.
The main scope of this procedure is to analyze and visualize correlated electron dynamics
by interpreting several fundamental one-electron quantities, such as the :ref:`detCI:Electron Density`,
the transient :ref:`detCI:Electronic Flux Density` (also called Current Density), and the :ref:`detCI:Electronic Dipole Moment`.

For the theoretical background, details on the implementation, and the robustness of these quantities, see 

  Vincent Pohl, Gunter Hermann, and Jean Christophe Tremblay,
  "An Open-Source Framework of Analyzing *N*-electron Dynamics: I. Multi-Determinantal Wave Functions", 
  `arXiv:1701.06885`__ (2017).

__ http://arxiv.org/abs/1701.06885

Please cite this publication, if you use detCI\@ORBKIT in your work. 

We provide well documented example files for both model systems described in this publication in the `examples folder`__.
The procedure of how to operate detCI\@ORBKIT is described in detail in the following chapters.

__ https://github.com/orbkit/orbkit/tree/cython/examples/detci

**Table of Contents:**

.. toctree::
   :maxdepth: 4

   detci

.. note::
  
  This section is under construction. If you have any questions,
  do not hesitate to contact us.
