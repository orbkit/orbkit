# -*- coding: iso-8859-1 -*-
'''
This file is part of orbkit. See the main program or documentation 
for information on the license.

Example file that shows how to use orbkit for ordering molecular orbitals
using analytical integrals.

Moreover, it shows...
  - how to read a list of files
  - how to order the molecular orbital coefficients
  - how to save and read the information obtained to and from an hdf5-file
  - how to depict several molecular orbitals
  - how to perform a standard orbkit computation for one molecular structures

Please note that the input files are compressed in .tar.gz file in the examples 
folder and need to be decompressed before running this example.
'''
import os
from time import time
from orbkit import Multi
from orbkit.read.high_level import main_read
from orbkit.output.high_level import main_output
from orbkit.display import init_display,display,tForm

create_plots = False # Specifies, if plots shall be created

# Options for the depiction of contour plots of selected molecular orbitals
selected_mos = ['24.1','23.2'] # Specifies, which MOs to be plotted
r0 = 1                         # Specifies the starting structure geo_spec_all[r0]
steps = 5                      # Specifies, how many steps to printed in one graph

select_slice = 'xz'            # Selects which plane to be plotted (Here, xz-plane)
where = 0.0                    # Selects where to place the plane (Here, y=0)

# A list containing all steps to be performed. Is not mandatory for the execution
run_id = ['Reading',
          'Plotting before ordering',
          'Ordering routine',
          'Plotting after ordering',
          'Plotting selected MOs',
          'Test the save and read routines',
          'Running orbkit for one selected structure'
          ]

init_display(name = 'mo_ordering')  # Specify a filename for the oklog file 

t = [time()]

#How are input files formatted?
fid_list = 'NaCl_molden_files.tar.gz'
#fid = 'nacl.%03d.molden'
#fid_list = []
#for i in range(0,16,1):
#  f = os.path.join(path,fid % i)
#  if not os.path.exists(f):
#    raise IOError('%s does not exist!' % f)
#  fid_list.append(f)

# Read all input files
mult = Multi()
mult.read(fid_list,itype='molden',all_mo=True,nosym=False)

t.append(time())

if create_plots:
  # Plot molecular orbital coefficients before the ordering routine starts
  for j,i in mult.sym.items():
    mo = mult.mo_coeff_all[i]
    mult.plot(mo,                   # Specifies the 3D matrix to be plotted
              symmetry='%s NaCl MO - Before Ordering' %j,
              output_format='png',
              ylim=[-10,10],
              thresh=1e-3,          # Curves smaller than this value are omitted
              plt_dir='Plots'       # Specifies the output folder
              )

t.append(time())

# Run the ordering routine using analytical overlap integrals
# Input argument None has been used because input files have been read already
index_list, mo_overlap = mult.order_using_analytical_overlap(None)

t.append(time())

if create_plots:
  # Plot molecular orbitals and overlap after the ordering routine
  cs = 0 
  for j,i in mult.sym.items():
    mo = mult.mo_coeff_all[i]
    mult.plot(mo,                   # Specifies the 3D matrix to be plotted
              symmetry='%s NaCl MO - Ordered' %j,
              output_format='png',
              ylim=[-10,10],
              thresh=1e-3,          # Curves smaller than this value are omitted
              plt_dir='Plots'       # Specifies the output folder
              )
    if 0:
      # Plot overlap matrix
      mult.plot(mo_overlap[cs],     # Specifies the 3D matrix to be plotted
                symmetry='%s NaCl MO Overlap Matrix'%j,
                output_format='png',
                ylim=[-1.1,1.1],
                plt_dir='Plots',    # Specifies the output folder
                x0=1                # Start x-axis from 1
                )
    cs += 1 

t.append(time())

if create_plots:
  # Plot selected molecular orbitals for some molecular structures
  mult.show_selected_mos(selected_mos,r0=r0,steps=steps)

t.append(time())

npz_fid = 'nacl.npz' # Specifies the filename of the npz file

# Save the results to a hdf5 file
main_output(mult, outputname=npz_fid, otype='native', ftype='numpy')


# Reset the `orbkit.multiple_files` module
try:
  from importlib import reload # >= Python3.4
except ImportError:
  try: 
    from imp import reload # <= Python3.3
  except ImportError:
    pass # Python2.X

# Read in the results
mult = main_read(npz_fid,gname='multi')

t.append(time())

# Perform a standard orbkit computation
mult.construct_qc() # Construct the qc_info class for every structure
r = 0                    # Index to be calculated
outputname = 'nacl_r%d' % r # Specifies the name of the output file

display('Running orbkit for the structure %d' % r)
import orbkit

# Initialize orbkit with default parameters and options
orbkit.init(reset_display=False)

# Set some options
orbkit.options.adjust_grid= [5, 0.1]                # adjust the grid to the geometry
orbkit.options.outputname = outputname              # output file (base) name
orbkit.options.otype      = 'h5'                    # output file type [default]
orbkit.options.numproc    = 4                       # number of processes

# Run orbkit with qc as input
data = orbkit.run_orbkit(mult.QC[r])

t.append(time())

# Print information about the timing
string = '\n%s\nRequired time:\n\t' % (80*'-')
for i,j in enumerate(run_id):
  if t[i+1]-t[i] > 1e-1:
    string += tForm(j,t[i+1]-t[i],extra='')[1:] + '\t'

display(string)
