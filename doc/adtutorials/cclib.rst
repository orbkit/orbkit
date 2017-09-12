.. _`Tutorial for Input Processing with cclib`:

Tutorial for Input Processing with cclib
========================================

This tutorial shows how to use the cclib_ Python program package to process 
quantum chemical output files together with ORBKIT.
It is prerequisite to install a current version cclib (>=1.3.1) beforehand by
following the installation guide:
  
  http://cclib.github.io/tutorial.html#manual-download-and-install

.. _cclib: http://cclib.github.io/

.. note::
  
  Currently, the cclib interface is only tested for Gaussian .log files.
  If you find any incompatibilities, please do not hesitate to contact us!

.. contents:: Table of Contents:
  :local:
  :depth: 1

General Aspects
---------------

In order to parse a quantum chemical output file with cclib, ORBKIT follows the 
standard procedure described in cclib's tutorial, e.g. for a Gaussian output 
file::

  from cclib.parser import Gaussian
  
  p = Gaussian('mycalc.log')
  ccData = p.parse()

Afterwards, ORBKIT converts the ``ccData`` data set to ORBKIT's ``QCinfo`` class 
(cf. :ref:`Central Variables`), i.e.::

  from orbkit import read
  
  qc = read.convert_cclib(ccData)

Usage via the Terminal
----------------------

For the standard ORBKIT terminal interface you have to specify two options: 
``--itype=cclib`` and ``--cclib_parser=CCLIB_PARSER``. The first tells ORBKIT 
to use cclib for parsing the output and the latter specifies which 
``cclib.parser`` has to be imported. (The available parser may be found on the
cclib_ website.)

Thus, to process the Gaussian output file ``mycalc.log``, you have to call:

.. code-block:: bash

  $ orbkit -i mycalc.log --itype=cclib --cclib_parser=Gaussian

.. attention::
  
  The option ``--cclib_parser`` is case sensitive! 

ORBKIT's High-Level Interface
-----------------------------

Again, to use the cclib paser, two options have to be called::

  import orbkit as ok

  ok.options.filename = 'mycalc.log'
  ok.options.itype = 'cclib'
  ok.options.cclib_parser = 'Gaussian'
  
  data = ok.run_orbkit()

.. attention::
  
  The option ``ok.options.cclib_parser`` is case sensitive! 

ORBKIT's Low-Level Interface
----------------------------

For ORBKIT's low-level interface, there are two possible ways to use cclib.
You can either tell ORBKIT, that it should use the ``cclib.parser`` ::

  from orbkit import read

  qc = read.main_read(filename,itype='cclib',all_mo=False,cclib_parser='Gaussian')

  # This is basically the same as
  qc = read.read_with_cclib(filename, cclib_parser='Gaussian', all_mo=False)

or you can do the parsing by yourself and let ORBKIT only convert the data::

  from orbkit import read
  from cclib.parser import Gaussian
  
  p = Gaussian('mycalc.log')
  ccData = p.parse()
  
  qc = read.convert_cclib(ccData)
