# -*- coding: iso-8859-1 -*-

'''
This file is part of orbkit. See the main program or documentation 
for information on the license.

Example file that shows how orbkit can be used as a module. Functions
of orbkit are used to read a molden file and to calculate a single 
molecular orbital. Finally, show cuts of the orbital in the xy-, 
xz- and yz-plane using matplotlib.

If mayavi2 is installed, you may set the variable mayavi_yes = True
to see an isosurface plot of the orbital.
'''

mayavi_yes = False

# import the functions of orbkit (import function 
# ok that imports all other functions)
from orbkit import ok

# name and type of input file
fid_in  = 'h2o.md'
itype   = 'molden'

# number of subprocesses
numproc = 4

# set grid parameters 
ok.grid.N_   = [  50,   52,   54]
ok.grid.max_ = [ 6.5,  6.5,   6.5]
ok.grid.min_ = [-6.5, -6.5,  -6.5]

# define variable to store quantum chemical information
qc = ok.QCinfo()

# open molden file and read parameters
qc = ok.read.main_read(fid_in,itype=itype)

# initialize grid
ok.grid.grid_init()

# print grid information
print(ok.grid.get_grid())

# define the molecular orbital to be calculated
selected_MO = ['3.1']

# get only the information of selected MO
Selected_mo_spec = []
for k in range(len(qc.mo_spec)):
  if qc.mo_spec[k]['sym'] in selected_MO:
    Selected_mo_spec.append(qc.mo_spec[k])

# calculate MO
mo_list = ok.core.rho_compute(qc.geo_spec,qc.ao_spec,
                              Selected_mo_spec,calc_mo=True,numproc=numproc)    

# plot the results
x = ok.grid.x
y = ok.grid.y
z = ok.grid.z

# if selected, use mayavi2 to make isosurface plot
if mayavi_yes == True:
    from enthought.mayavi import mlab
    src = mlab.pipeline.scalar_field(mo_list[0])
    mlab.pipeline.iso_surface(\
        src, contours=[0.001, ], opacity=0.3, color=(0, 0, 0.8))
    mlab.pipeline.iso_surface(\
        src, contours=[-0.001, ], opacity=0.3, color=(0.8, 0, 0))
    mlab.show()

# use matplotlib to show cuts of the molecular orbitals
import matplotlib.pyplot as plt
import numpy as np

# select cuts
xd = mo_list[0][ok.grid.N_[0]/2-1,:,:]
yd = mo_list[0][:,ok.grid.N_[1]/2-1,:]
zd = mo_list[0][:,:,ok.grid.N_[2]/2-1]

# plot cuts
f, (pic1, pic2, pic3) = \
            plt.subplots(3,1,sharex=True,sharey=True,figsize=(6,14))
pic1.contour(z,y,xd,50,linewidths=0.5,colors='k')
pic1.contourf(\
    z,y,xd,50,cmap=plt.cm.rainbow,vmax=abs(xd).max(),vmin=-abs(xd).max())
pic1.set_xlabel('z')
pic1.set_ylabel('y')

pic2.contour(z,x,yd,50,linewidths=0.5,colors='k')
pic2.contourf(\
    z,x,yd,50,cmap=plt.cm.rainbow,vmax=abs(yd).max(),vmin=-abs(yd).max())    
pic2.set_xlabel('z')
pic2.set_ylabel('x')    

pic3.contour(y,x,zd,50,linewidths=0.5,colors='k')
pic3.contourf(\
    y,x,zd,50,cmap=plt.cm.rainbow,vmax=abs(zd).max(),vmin=-abs(zd).max())  
pic3.set_xlabel('y')
pic3.set_ylabel('x')     

# following options applied for all subplots as they share x- and y-axis
pic1.xaxis.set_ticks(np.arange(-5,6,5))
pic1.yaxis.set_ticks(np.arange(-5,6,5))
pic1.set_aspect('equal')

# plot
f.subplots_adjust(left=0.15,bottom=0.05,top=0.95,right=0.95)
f.show()
