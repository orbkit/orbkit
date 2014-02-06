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
Example file for the execution of orbkit
for calculating the derivative of the electron density.

'''
import time

import numpy
import pylab as plt

#--- Import the orbkit_main module ---
import orbkit_main

# Measurement of required execution time
t=[time.time()]
path = '/home/vincpohl/mnt/soroban/CalcMaster/H2+/'


out_fid = 'h2o'

#--- Initilize orbkit ---
orbkit_main.init()

#--- Set the options ---
#--- Grid for H2O
in_fid  = 'h2o.md'
orbkit_main.grid.N_   		= [  1,    50,    50]
orbkit_main.grid.max_ 		= [ 0.0,  5.0,   5.0]
orbkit_main.grid.min_ 		= [ 0.0, -5.0,  -5.0]
orbkit_main.options.hdf5	= True
orbkit_main.options.filename	= in_fid 
orbkit_main.options.outputname	= out_fid

orbkit_main.core.options = orbkit_main.options

geo_spec, geo_info, ao_spec, mo_spec = orbkit_main.main_read(orbkit_main.options.filename)

orbkit_main.grid.grid_init()
orbkit_main.grid.grid_display()

rho, delta_rho = orbkit_main.core.delta_rho_compute(geo_spec, geo_info, ao_spec, mo_spec)

t.append(time.time())

orbkit_main.orbkit_output.display(orbkit_main.tForm('The calculation took',t[-1]-t[0]))

#--- Plot the vector fields ---
#--- as a 2D Plot if len(x) or len(y) is 1 ---
#--- else show a 3D Plot ---

x = orbkit_main.grid.x
y = orbkit_main.grid.y
z = orbkit_main.grid.z

b2D = None
if len(x) == 1:
  r = y
  d_r = delta_rho[1,0,:,:]
  d_z = delta_rho[2,0,:,:]
  b2D = True
elif len(y) == 1:
  r = x
  d_r = delta_rho[0,:,0,:]
  d_z = delta_rho[2,:,0,:]
  b2D = True


if b2D:
  [Z,R] = numpy.meshgrid(z,r)
  fig = plt.figure(figsize = (6,5))

  Q = plt.quiver(Z, R, d_z, d_r, scale=10, units='width')
  qk = plt.quiverkey(Q, 0.5, 0.95, 1, r'$1\,\frac{1}{\rm a^2_0}$',
		labelpos='E',
		coordinates='figure',
		fontproperties={'weight': 'bold'})

  plt.xlabel(r'$y [\sf{a}_0]$')
  plt.ylabel(r'$z [\sf{a}_0]$')
  
  plt.show()
  
else:
  from enthought.mayavi import mlab
  
  mlab.pipeline.iso_surface(src, contours=[0.08, ], opacity=0.3, color=(0.8, 0.8, 0.8))
