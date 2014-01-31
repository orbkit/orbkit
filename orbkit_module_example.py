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
as ...
'''

#--- Import the orbkit_main module ---
import orbkit_main

#---- First scenario ----
'''
...
'''
in_fid  = 'h2o.md'
out_fid = 'h2o'

#--- Initilize orbkit ---
orbkit_main.init()

#--- Set the options ---
orbkit_main.grid.N_   		= [  100,   100,   100]
orbkit_main.grid.max_ 		= [ 10.0,  10.0,   10.0]
orbkit_main.grid.min_ 		= [-10.0, -10.0,  -10.0]
orbkit_main.options.hdf5	= True
orbkit_main.options.filename	= in_fid 
orbkit_main.options.outputname	= out_fid

orbkit_main.core.options = orbkit_main.options

geo_spec, geo_info, ao_spec, mo_spec = orbkit_main.main_read(orbkit_main.options.filename)

orbkit_main.grid.grid_init()
orbkit_main.grid.grid_display()

selected_MO = ['3.1']

Selected_mo_spec = []

for k in range(len(mo_spec)):
  if mo_spec[k]['sym'] in selected_MO:
    Selected_mo_spec.append(mo_spec[k])

ao_list = orbkit_main.core.ao_creator(geo_spec,ao_spec)
mo_list = orbkit_main.core.mo_creator(ao_list,Selected_mo_spec)

import numpy
from enthought.mayavi import mlab

x = orbkit_main.grid.x
y = orbkit_main.grid.y
z = orbkit_main.grid.z

src = mlab.pipeline.scalar_field(mo_list[0])

mlab.pipeline.iso_surface(src, contours=[0.001, ], opacity=0.3, color=(0, 0, 0.8))
mlab.pipeline.iso_surface(src, contours=[-0.001, ], opacity=0.3, color=(0.8, 0, 0))

mlab.show()