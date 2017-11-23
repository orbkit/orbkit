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
orbkit.options.outputname   = 'h2o_MO'       # output file (base) name
orbkit.options.otype        = 'h5'           # output file type
orbkit.options.numproc      = 1              # number of processes for multiprocessing
orbkit.options.slice_length = 4e4            # number of points per process
orbkit.options.calc_mo      = 'MO_List.tab'  # list of molecular orbitals to be used

# first run: do not calculate derivatives
orbkit.options.drv          = None           # do not calculate derivative
orbkit.options.no_output    = True           # we will create our own output

# run orbkit
mo = orbkit.run_orbkit()

# create output: molecular orbital data
orbkit.output.hdf5_creator(mo,
    orbkit.options.outputname,
    orbkit.main.qc.geo_info,orbkit.main.qc.geo_spec,
    data_id='MO',               # name of data set
    append=None,                # create new file [default]
    data_only=False,            # include grid, structure, and MO data [default]
    is_mo_output=True,
    mo_spec=orbkit.main.qc.mo_spec)

# second run: calculate derivatives
orbkit.options.drv       = ['x', 'y', 'z'] # calculate derivatives along x, y and z
orbkit.options.no_output = True            # we will create our own output

# run orbkit
mo = orbkit.run_orbkit()

# append output: derivative data
orbkit.output.hdf5_creator(mo,orbkit.options.outputname,None,None,
  data_id='delta_MO',             # name of data set
  append='/',                     # where to append in file
  data_only=True,                 # do not include grid, structure, and MO data
  is_mo_output=False,
  mo_spec=orbkit.main.qc.mo_spec)

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
    mlab.quiver3d(x,y,z,mo[0,mo_num,:],mo[1,mo_num,:],\
                    mo[2,mo_num,:],line_width=1.5,scale_factor=0.1)
    mlab.show()


