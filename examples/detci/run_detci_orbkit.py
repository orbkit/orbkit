from orbkit import options, grid, read ,core, output
from orbkit.display import init_display, display, good_bye_message
import numpy, h5py,time

options.quiet = False

# Begin writing logfile
out_fid = 'testing_detci_orbkit'
init_display(out_fid)

# Start time measurement
overalltime = [time.time()]  

# Set some options
numproc = 4

# Load orbkit module
from orbkit.detci import ci_read

# Read CIS states
fid_qc = 'water_CIS_cc-pVDZ.log'
qc = read.main_read(fid_qc,itype='gamess',all_mo=True)
select_state = range(2)
ci = ci_read.gamess_cis(fid_qc,select_state=select_state,threshold=0.0)
moocc = [i['occ_num'] for i in qc.mo_spec]
moocc = numpy.array(moocc,dtype=numpy.intc)

## Read MCSCF states
#fid_molpro = 'ch2o.out'
#fid_molden = 'ch2o.molden'
#qc = read.main_read(fid_molden,all_mo=True)
#ci = ci_read.molpro_mcscf(fid_molpro,select_run=0,threshold=0.0)
#closed,active,external = ci_read.molpro_mo_order_ci(ci[0].info['occ_info'],qc.mo_spec)
#qc.mo_spec = closed+active+external
#moocc = numpy.zeros(len(closed),dtype=numpy.intc) + 2

# Initialize grid
grid.adjust_to_geo(qc,extend=5.0,step=0.5)
display('Setting up the grid...')
grid.grid_init(is_vector=True)
display(grid.get_grid())
slice_length = grid.N_[1]*grid.N_[2]/2

# Calculate molist, molistdrv, moom, omr, omv
display('Computing molecular orbitals...')
molist = core.rho_compute(qc,
                          calc_mo=True,
                          slice_length=slice_length,
                          drv=None,
                          numproc=numproc)
display('Computing derivatives of molecular orbitals...')
molistdrv = core.rho_compute(qc,
                             calc_mo=True,
                             slice_length=slice_length,
                             drv='xyz',
                             numproc=numproc)

from orbkit.analytical_integrals import get_ao_overlap,get_mo_overlap_matrix
from orbkit.analytical_integrals import get_ao_dipole_matrix

display('Computing analytical overlaps of molecular orbitals...')
aoom = get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec)
moom = get_mo_overlap_matrix(qc.mo_spec,qc.mo_spec,aoom,numproc=numproc)

delta_moom = []
dm_moom = []
for component in 'xyz':
  delta_moom.append(
    get_mo_overlap_matrix(qc.mo_spec,
                          qc.mo_spec,
                          get_ao_overlap(qc.geo_spec,
                                         qc.geo_spec,
                                         qc.ao_spec,
                                         drv=component),
                          numproc=numproc)
                          )
  dm_moom.append(
    get_mo_overlap_matrix(qc.mo_spec,
                          qc.mo_spec,
                          get_ao_dipole_matrix(qc,component=component),
                          numproc=numproc)
                          )

omr = numpy.array(dm_moom)
omv = numpy.array(delta_moom)

# CI calculations
from orbkit.detci import occ_check,ci_core
overalltime.append(time.time())
# Create HDF5 File
output.hdf5_write('%s.h5' % out_fid,mode='w',gname='general_info',
                  x=grid.x,y=grid.y,z=grid.z,N=grid.N_,
                  geo_info=qc.geo_info,geo_spec=qc.geo_spec,
                  nstates=len(ci))
st_comb = (len(ci)*len(ci) + len(ci))/2
f = h5py.File('%s.h5' % out_fid, 'a')
dset = f.create_dataset('Rhonm',[st_comb,len(grid.x)])
dset = f.create_dataset('Rhon',[len(ci),len(grid.x)])
dset = f.create_dataset('J',[st_comb,3,len(grid.x)])
dset = f.create_dataset('en',[len(ci)])
dset = f.create_dataset('mur',[st_comb,3]) 
dset = f.create_dataset('muv',[st_comb,3]) 
dset = f.create_dataset('mur_v',[st_comb,3]) 

count = 0
states = []
for a in range(len(ci)):
  display('\nState %s:' % ci[a].info['state'])
  zero,sing = occ_check.compare(ci[a],ci[a],moocc,numproc=numproc)
  display('\n\tElectron density:')
  rhon = ci_core.rho(zero,sing,molist,slice_length=slice_length,numproc=numproc)
  en = ci_core.enum(zero,sing,moom)
  display('\n\tNorm of the state: %f' %en)
  # Save data to HDF5 File
  f['Rhon'][a,:] = rhon
  f['en'][a] = en
  for b in range(a+1,len(ci)):
    display('\nState %s to %s:'% (ci[a].info['state'], ci[b].info['state']))
    zero,sing = occ_check.compare(ci[a],ci[b],moocc,numproc=numproc)
    display('\n\tTransition electron density:')
    rhonm = ci_core.rho(zero,sing,molist,slice_length=slice_length,
                        numproc=numproc)
    display('\n\tTransition electronic flux density:')
    jnm = ci_core.jab(zero,sing,molist,molistdrv,slice_length=slice_length,
                      numproc=numproc)    
    mu = ci_core.mu(ci[a],ci[b],qc,zero,sing,omr,omv)
    display('\n\tTransition dipole moment in length gauge:')
    display('\t%+.04f\t%+.04f\t%+.04f' % tuple(mu[0]))
    # Save data to HDF5 File
    f['Rhonm'][count,:] = rhonm
    f['J'][count,:,:] = jnm
    f['mur'][count,:] = mu[0]
    f['muv'][count,:] = mu[1]
    f['mur_v'][count,:] = mu[2]
    states.append([a,b])
    count += 1

f['states'] = states
f.close()

overalltime.append(time.time())

good_bye_message(overalltime)

