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

import numpy as np

# import the functions of orbkit (import function 
# ok that imports all other functions)
import orbkit as ok

# defining some grid in spherical coordinates
r     = np.arange(0.1,    3.,              0.2)
theta = np.arange(0,      np.pi+0.1, np.pi/10.)
phi   = np.arange(-np.pi, np.pi+0.1, np.pi/10.)

# initialize orbkit with default parameters and options
ok.init()

# Setting up a spherical grid -> also sets option for vector grid
ok.grid.sph2cart_vector(r,theta,phi)

# orbkit options
ok.options.filename     = 'h2o.molden'       # input file name
ok.options.itype        = 'molden'       # input file type
ok.options.outputname   = 'h2o_MO'       # output file (base) name
ok.options.otype        = 'h5'           # output file type
ok.options.numproc      = 1              # number of processes for multiprocessing
ok.options.slice_length = 4e4            # number of points per process
ok.options.calc_mo      = 'MO_List.tab'  # list of molecular orbitals to be used

# first run: do not calculate derivatives
ok.options.drv          = None           # do not calculate derivative
ok.options.no_output    = True           # we will create our own output

# run orbkit
mo = ok.run_orbkit()

# create output: molecular orbital data
ok.output.HDF5_creator(mo,
    ok.options.outputname,
    ok.main.qc.geo_info,ok.main.qc.geo_spec,
    data_id='MO',               # name of data set
    append=None,                # create new file [default]
    data_only=False,            # include grid, structure, and MO data [default]
    is_mo_output=True,
    mo_spec=ok.main.qc.mo_spec)

# second run: calculate derivatives
ok.options.drv       = ['x', 'y', 'z'] # calculate derivatives along x, y and z
ok.options.no_output = True            # we will create our own output

# run orbkit
mo = ok.run_orbkit()

# append output: derivative data
ok.output.HDF5_creator(mo,ok.options.outputname,None,None,
  data_id='delta_MO',             # name of data set
  append='/',                     # where to append in file
  data_only=True,                 # do not include grid, structure, and MO data
  is_mo_output=False,
  mo_spec=ok.main.qc.mo_spec)

# plot derivative
x = ok.grid.x
y = ok.grid.y
z = ok.grid.z

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


