'''
Partial Charges and Gross Atomic Density

Example for the orbkit's non-grid based capabilities and for the computation
of the gross atomic density.

Hint: Use the VMD script file (ch2o.vmd) for depicting the results.
'''

from orbkit import read, grid, extras, output, atomic_populations, display

# Open molden file and read parameters
qc = read.main_read('ch2o.molden',itype='molden')

# Perform a Mulliken population analysis and write the output to a PDB file
pop = atomic_populations.mulliken(qc)
output.pdb_creator(qc.geo_info,qc.geo_spec,filename='ch2o',charges=pop['charge'])

# Initialize the grid
display.display('Setting up the grid...')
grid.adjust_to_geo(qc,extend=5.0,step=0.1)
grid.grid_init()
display.display(grid.get_grid())

# Compute and save the gross atomic density of the C atom
rho_atom = extras.gross_atomic_density(1,qc)
output.main_output(rho_atom[0],qc,outputname='ch2o',
                   otype='cb')

