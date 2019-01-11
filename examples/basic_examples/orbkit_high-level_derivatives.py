# -*- coding: iso-8859-1 -*-

'''
This file is part of orbkit. See the main program or documentation 
for information on the license.

Example file that shows how orbkit can calculate gradients of 
molecular orbitals. Also, a user defined grid in spherical 
coordinates will be used for calculation.

If mayavi2 is installed, the vector field of the gradient of
one MO will be shown.
'''

import numpy, orbkit

# import the functions of orbkit (import function 
# ok that imports all other functions)
# defining some grid in spherical coordinates
r     = numpy.arange(0.1,       3.,              0.2)
theta = numpy.arange(0,         numpy.pi+0.1, numpy.pi/10.)
phi   = numpy.arange(-numpy.pi, numpy.pi+0.1, numpy.pi/10.)

# initialize orbkit with default parameters and options
orbkit.init()

# Setting up a spherical grid -> also sets option for vector grid
orbkit.grid.sph2cart_vector(r,theta,phi)

# orbkit options
orbkit.options.filename     = 'h2o.molden'       # input file name
orbkit.options.itype        = 'molden'       # input file type
orbkit.options.numproc      = 1              # number of processes for multiprocessing
orbkit.options.slice_length = 4e4            # number of points per process
orbkit.options.calc_mo      = 'MO_List.tab'  # list of molecular orbitals to be used

# first run: do not calculate derivatives
orbkit.options.drv       = ['None', 'x', 'y', 'z'] # calculate mos and their derivatives along x, y and z
orbkit.options.no_output = True            # we will create our own output

# run orbkit
mos = orbkit.run_orbkit()

# Write a npz file with all the ouput
orbkit.output.main_output(mos, outputname='h2o_MO', otype='npz') 

# plot derivative
x = orbkit.grid.x
y = orbkit.grid.y
z = orbkit.grid.z

# test if mayavi exists and possibly plot
maya = False
try:
    from enthought.mayavi import mlab
    maya = True
except ImportError:
    pass
try:
    from mayavi import mlab
    maya = True
except Exception:
    pass

if maya == False:
    print('mayavi module could not be loaded and vector field cannot be shown')
else:
    mo_num = 3
    mlab.quiver3d(x,y,z,mos[1,mo_num,:],mos[2,mo_num,:],\
                    mos[3,mo_num,:],line_width=1.5,scale_factor=0.1)
    mlab.show()


