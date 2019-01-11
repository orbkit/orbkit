# -*- coding: iso-8859-1 -*-

'''
This file is part of orbkit. See the main program or documentation 
for information on the license.

Example file that shows how to use the orbkit main function to 
calculate electron densities.
'''

# import the functions of orbkit (import function 
# ok that imports all other functions)
import orbkit

#---- first scenario ----
'''
Calculate the electron density with 4 parallel subprocesses
on a user-defined grid, save data to a hdf5 file and print 
number of electrons
'''
print('\n \n First scenario: calculate electron density. \n \n')

in_fid  = 'h2o.molden'
out_fid = 'h2o_r1'

# initialize orbkit with default parameters and options
orbkit.init()

# set some options
orbkit.grid.N_            = [  201,   201,   101]   # grid points (regular grid)
orbkit.grid.max_          = [ 10.0,  10.0,   5.0]   # maximum grid value
orbkit.grid.min_          = [-10.0, -10.0,  -5.0]   # minimum grid value
orbkit.options.filename   = in_fid                  # input file name
orbkit.options.itype      = 'molden'                # input file type [default]
orbkit.options.outputname = out_fid                 # output file (base) name
orbkit.options.otype      = 'h5'                    # output file type [default]
orbkit.options.numproc    = 4                       # number of processes

# run orbkit
rho = orbkit.run_orbkit()

import numpy as np
print('We have %.3f electrons.' % (np.sum(rho)*orbkit.grid.d3r))

#---- Second scenario ----
'''
Calculate the density and its x- and z-derivative using some orbitals
defined in an external file (by symmetry notation), read the grid from a file,
and save the data to cube and hdf5 files.

Note that for the derivatives, there will be one cube file per MO per direction,
while there will only be one hdf5 file containing information of all MOs per 
direction.

Ommit the creation of the *.oklog file
'''
print('\n \n Second scenario: calculate derivatives of electron density. \n \n')

in_fid  = 'h2o.molden'
out_fid = 'h2o_r2'

# initialize orbkit with default parameters and options
orbkit.init()

# set some options
orbkit.options.grid_file  = 'grid_reg.txt'          # grid file to read from (regular grid)
orbkit.options.filename   = in_fid                  # input file name
orbkit.options.itype      = 'molden'                # input file type [default]
orbkit.options.outputname = out_fid                 # output file (base) name
orbkit.options.otype      = ['h5']                  # output file types
orbkit.options.calc_mo    = ['2.1','1.1','2.3']     # list of MO labels (Molpro notation)
orbkit.options.drv        = ['x', 'z']              # derivatives along x and z
orbkit.options.no_log     = True                    # do not write log this time

# run orbkit
orbkit.run_orbkit()


#---- third scenario ----
'''
Create an ZIBAmira Network from the density caclulated of MOs
selected by their number in the molden file
'''
print('\n \n Third scenario: calculate electron density for Amira and VMD. \n \n')

in_fid  = 'h2o.molden'
out_fid = 'h2o_r3'

# initialize orbkit with default parameters and options
orbkit.init()

# set some options
orbkit.options.adjust_grid= [5, 0.1]                # adjust the grid to the geometry
orbkit.options.otype      = ['vmd']           # output file types
orbkit.options.mo_set     = [[1,2,'homo-2'],
                         ['homo-1:lumo']]           # list of MO labels (Molden enumeration)
orbkit.options.filename   = in_fid                  # input file name
orbkit.options.itype      = 'molden'                # input file type [default]
orbkit.options.outputname = out_fid                 # output file (base) name
orbkit.options.numproc    = 2                       # number of processes

# run orbkit
orbkit.run_orbkit()
