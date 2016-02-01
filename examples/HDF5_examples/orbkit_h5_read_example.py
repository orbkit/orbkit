# -*- coding: iso-8859-1 -*-

'''
This file is part of orbkit. See the main program or documentation 
for information on the license.

Example file that shows how to read an HDF5 orbkit output with python.
'''

# Import the HDF5 python module
import h5py

fid = 'h2o.h5'

# Open the HDF5 File created by orbkit
HDF5_file = h5py.File(fid,'r')

# Load the grid
x = HDF5_file['/x'][...]
y = HDF5_file['/y'][...]
z = HDF5_file['/z'][...]

# Load the density
rho = HDF5_file['/rho'][...]

# Do something with the data
import numpy
d3r = (x[0,1]-x[0,0]) * (y[0,1]-y[0,0]) * (z[0,1]-z[0,0])
print('We have %.3f electrons.' % (numpy.sum(rho)*d3r))
