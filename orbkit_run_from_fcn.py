# -*- coding: iso-8859-1 -*-

'''
orbkit
Axel Schild (axel.schild [at] fu-berlin.de)
Gunter Hermann
Vincent Pohl

Institut fuer Chemie und Biochemie, Freie Universitaet Berlin, 14195 Berlin, Germany

This file is part of orbkit.

orbkit is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or any later version.

orbkit is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with orbkit.  If not, see <http://www.gnu.org/licenses/>.
'''

'''
Example file for the execution of the MAGIC DENSITY CREATOR
as a python module in a python program.

The program is started for several scenarios. 

An attempt was made to cover all options implemented.
'''

#--- Import the orbkit_main module ---
import orbkit_main

#---- First scenario ----
'''
Calculate the electron density, 
work with 2 parallel processes (default = 4)
Set the grid by hand, and save the data to a HDF5 file.
Print out again the number of electrons externally.
'''
in_fid  = 'h2o.md'
out_fid = 'h2o'

##--- Initilize orbkit ---
orbkit_main.init()

#--- Set the options ---
orbkit_main.grid.N_   		= [  300,   450,   275]
orbkit_main.grid.max_ 		= [ 18.0,  10.0,   8.0]
orbkit_main.grid.min_ 		= [-18.0, -10.0,  -8.0]
orbkit_main.options.numproc	= 2
orbkit_main.options.hdf5	= True
orbkit_main.options.filename	= in_fid 
orbkit_main.options.outputname	= out_fid

#--- Run the orbkit ---
orbkit_main.main()

import numpy
print('We have around %.3f electrons.' % (numpy.sum(orbkit_main.rho)*orbkit_main.grid.d3r))

#-------------------------------------------------------

#---- Second scenario ----
'''
Calculate the reduced density (with respect to the z-axis) of MOs selected by symmetry,
Read the grid from a .cvs file, and
save the data to HDF5 files.
'''
in_fid  = 'h2o.md'
out_fid = 'h2o'

#--- Initilize orbkit ---
orbkit_main.init()

#--- Set the options ---
orbkit_main.grid.N_   		= [  300,   450,   275]
orbkit_main.grid.max_ 		= [ 18.0,  10.0,   8.0]
orbkit_main.grid.min_ 		= [-18.0, -10.0,  -8.0]
orbkit_main.options.quiet	= False
orbkit_main.options.hdf5	= True
orbkit_main.options.reduced_density = True
orbkit_main.options.csv_grid	= './Tab.csv'
orbkit_main.options.calc_mo	= './MO_List.tab'
orbkit_main.options.filename	= in_fid 
orbkit_main.options.outputname	= out_fid


#--- Run the orbkit ---
orbkit_main.main()

#-------------------------------------------------------

#---- Third Run ----
'''
Create an ZIBAmira Network from the density caclulated of MOs
selected by their specific number in the molden file
'''
in_fid  = 'h2o.md'
out_fid = 'h2o'

#--- Initilize orbkit ---
orbkit_main.init()

#--- Set the options ---
orbkit_main.options.amira	= True	# Has not necessarily to be set
orbkit_main.options.hx_network	= True
orbkit_main.options.mo_list	= './MO_List_int.tab'
orbkit_main.options.filename	= in_fid 
orbkit_main.options.outputname	= out_fid

#--- Run the orbkit ---
orbkit_main.main()

#-------------------------------------------------------

#---- Multiprocessing Speed Test ----
'''
Uncomment this section to test the speed dependence 
with respect to the number of processes run.
'''
#in_fid  = 'h2o.md'
#out_fid = 'h2o'

#for ii in [1, 2, 3, 4, 5, 6, 7, 8, 10, 20]:
  #t.append(time.time())
  ##--- Initilize orbkit ---
  #orbkit_main.init()

  ##--- Set the options ---
  #orbkit_main.options.no_log	= True
  #orbkit_main.options.quiet	= True
  #orbkit_main.options.hdf5	= True
  #orbkit_main.options.numproc	= ii
  #orbkit_main.options.csv_grid	= './Tab.csv'
  #orbkit_main.options.filename	= in_fid 
  #orbkit_main.options.outputname = out_fid


  ##--- Run the orbkit ---
  #orbkit_main.main()

  #t.append(time.time()) # Final time
  #print(orbkit_main.tForm('The test run with %d processes took' % ii, t[-1]-t[-2]))
