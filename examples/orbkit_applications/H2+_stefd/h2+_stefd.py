# Import the functions of orbkit
from orbkit import grid,read,core,output,display

# Import general modules
import numpy
try:
  from enthought.mayavi import mlab
except ImportError:
  from mayavi import mlab

# Name of input and output files
fid_in  = 'h2+.molden'

# Choose the molecular orbitals to be calculated
selected_MO = ['1.1_a','1.5_a']

# number of processes and points per process
numproc = 2
slice_length = 1e4

# Open molden file and read parameters
qc = read.main_read(fid_in, itype='molden',all_mo=True)

# Set grid parameters 
grid.adjust_to_geo(qc,extend=5.0,step=0.5)
# Initialize grid
grid.grid_init()
# Print grid information
display.display(grid.get_grid())

# Choose the molecular orbitals to be calculated
selected_MO = ['1.1 alpha','1.5 alpha']
qc.mo_spec = qc.mo_spec.select(selected_MO)

# Calculate molecular orbitals
mo_list = core.rho_compute(qc,calc_mo=True,slice_length=slice_length,drv=None,
                              numproc=numproc) 

# Calculate analytic derivatives of molecular orbitals 
mo_list_drv = core.rho_compute(qc,calc_mo=True,slice_length=slice_length,drv='xyz',
                                  numproc=numproc)

# Calculate Transition Electronic Flux Densities (time independent)

# Calculate the STEFD [in units of (i E_h/(hbar a_0^2))]
j_stefd = -0.5 * (  mo_list[numpy.newaxis,0] * mo_list_drv[:,1] 
                  - mo_list[numpy.newaxis,1] * mo_list_drv[:,0])

Z,Y,X = output.meshgrid2(*grid.tolist()[::-1])

mlab.quiver3d(X,Y,Z,*j_stefd)
mlab.show()
