Tutorial for Atomic Population and Charge Analysis
==================================================

This short tutorial shows, how to perform a Mulliken and Löwdin population
analysis with ORBKIT analytically, and how to write the output to an ``xyz`` or
a ``PDB`` file.

Computation of the Atomic Populations and Charges
-------------------------------------------------

First, we have to import some modules and set some of ORBKITs options::

  from orbkit import read, atomic_populations

Then, we have to read the input file::

  # Open molden file and read parameters
  qc = read.main_read('h2o.molden')

Now, we have to call the respective functions:
  
For a **Mulliken Population Analysis**::
  
  pop = atomic_populations.mulliken(qc)

and for a **Löwdin Population Analysis**::
  
  pop = atomic_populations.lowdin(qc)

The return value ``pop`` is a dictionary, which contains information about the 
population analysis and has following members:

+-----------------+----------------------------------------+
| Variable        | Contents                               |
+-----------------+----------------------------------------+
|``population``   | Contains the population for each atom. |
+-----------------+----------------------------------------+
|``charge``       | Contains the charges for each atom.    |
+-----------------+----------------------------------------+

Creation of an Output File
--------------------------

ORBKIT provides two possible output file formats for the population analysis: 
the ``PDB`` and the ``xyz`` file format::

  from orbkit.output import pdb_creator,xyz_creator
  # Write a PDB file
  pdb_creator(qc.geo_info,qc.geo_spec,filename='charges',charges=pop['charge'],comments='')
  # Write an xyz file
  xyz_creator(qc.geo_info,qc.geo_spec,filename='charges',charges=pop['charge'],comments='')

This will write two files: ``charges.pdb`` and ``charges.xyz``, which can be 
processed for visualization by several programs. 

For instance, you can depict the PDB file containing the partial charges in VMD:

  1. Load the file:
  
    .. code-block:: bash  

      $ vmd -m charges.pbd

  2. Go to **Graphics->Representations** 
  3. Select ``CPK`` as **Drawing Method** 
  4. Select ``Beta`` as **Coloring Method** 
