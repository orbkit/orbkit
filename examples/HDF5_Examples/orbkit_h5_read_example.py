# -*- coding: iso-8859-1 -*-

'''
orbkit
Axel Schild (axel.schild [at] fu-berlin.de)
Gunter Hermann
Vincent Pohl


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

# Import the HDF5 python module
import h5py

fid = 'h2o.h5'

# Open the HDF5 File created by orbkit
HDF5_file = h5py.File(fid,'r')

# Load Grid
x = HDF5_file['/x']
x = x[...]
y = HDF5_file['/y']
y = y[...]
z = HDF5_file['/z']
z = z[...]

# Load the density calculated
rho = HDF5_file['/rho']
rho = rho[...]

# Do something with the data
import numpy
d3r = (x[0,1]-x[0,0]) * (y[0,1]-y[0,0]) * (z[0,1]-z[0,0])
print('We have %.3f electrons.' % (numpy.sum(rho)*d3r))
