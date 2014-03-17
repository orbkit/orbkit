# -*- coding: iso-8859-1 -*-

'''
This file is part of orbkit. See the main program or documentation 
for information on the license.

Example file that shows how to use the orbkit main function to 
calculate electron densities.
'''

# import the functions of orbkit (import function 
# ok that imports all other functions)
from orbkit import ok


#---- first scenario ----
'''
Calculate the electron density with 4 parallel subprocesses
on a user-defined grid, save data to a HDF5 file and print 
number of electrons
'''
print('\n \n First scenario: calculate electron density. \n \n')

in_fid  = 'h2o.md'
out_fid = 'h2o_r1'

# initialize orbkit with default parameters and options
ok.main.init()

# set some options
ok.grid.N_            = [  201,   201,   101]   # grid points (regular grid)
ok.grid.max_          = [ 10.0,  10.0,   5.0]   # maximum grid value
ok.grid.min_          = [-10.0, -10.0,  -5.0]   # minimum grid value
ok.options.filename   = in_fid                  # input file name
ok.options.itype      = 'molden'                # input file type [default]
ok.options.outputname = out_fid                 # output file (base) name
ok.options.otype      = 'h5'                    # output file type [default]
ok.options.numproc    = 4                       # number of processes

# run orbkit
ok.main.main()

import numpy as np
print('We have %.3f electrons.' % (np.sum(ok.main.rho)*ok.grid.d3r))


#---- Second scenario ----
'''
Calculate the density and its x- and z-derivative using some orbitals
defined in an external file (by symmetry notation), read the grid from a file,
and save the data to cube and HDF5 files.

Note that for the derivatives, there will be one cube file per MO per direction,
while there will only be one HDF5 file containing information of all MOs per 
direction.

Ommit the creation of the *.oklog file
'''
print('\n \n Second scenario: calculate derivatives of electron density. \n \n')

in_fid  = 'h2o.md'
out_fid = 'h2o_r2'

# initialize orbkit with default parameters and options
ok.main.init()

# set some options
ok.options.grid_file  = 'grid_reg.txt'   # grid file to read from (regular grid)
ok.options.filename   = in_fid           # input file name
ok.options.itype      = 'molden'         # input file type [default]
ok.options.outputname = out_fid          # output file (base) name
ok.options.otype      = ['h5', 'cb']     # output file types
ok.options.calc_mo    = './MO_List.tab'  # list of MOs (Molpro notation)
ok.options.drv        = ['x', 'z']       # derivatives along x and z
ok.options.no_log     = True             # do not write log this time

# run orbkit
ok.main.main()


#---- third scenario ----
'''
Create an ZIBAmira Network from the density caclulated of MOs
selected by their number in the molden file
'''
print('\n \n Third scenario: calculate electron density for Amira. \n \n')

in_fid  = 'h2o.md'
out_fid = 'h2o_r3'

# initialize orbkit with default parameters and options
ok.main.init()

# set some options
ok.options.otype      = ['am', 'hx']         # output file types
ok.options.mo_list    = './MO_List_int.tab'  # list of MOs (Molden enumeration)
ok.options.filename   = in_fid 
ok.options.outputname = out_fid

# run orbkit
ok.main.main()
