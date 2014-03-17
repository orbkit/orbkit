# -*- coding: iso-8859-1 -*-
'''Contains and processes all orbkit options.'''
import os
import sys
thismodule = sys.modules[__name__]

from orbkit import grid

available = ['filename','itype','outputname','otype',
	     'vector','grid_file','center_grid',
	     'numproc','mo_set','all_mo','drv',
	     'z_reduced_density','atom_projected_density','mo_tefd',
	     'quiet','no_log','no_output','no_slice','random_grid',
	     'interactive']

__all__ = ['get_options','check_options','reset_grid',].extend(available)

itypes = ['molden', 'gamess', 'gaussian.log', 'gaussian.fchk']	#:
otypes = ['h5', 'cb', 'am', 'hx']	#:
drv_options = ['x','y','z']	#:

def get_options():
  '''Return all possible options and their value'''
  opts = [(i,globals()[i]) for i in available]
  return dict(opts)

def check_options(error=sys.stdout.write,display=sys.stdout.write,
		  interactive=False,
		  info=True):
  '''Check the options for errors
  
  :Parameters:
  options : instance
	All orbkit options
  error : function, optional
	A function or something equivalent to handle the errors
  interactive : bool, optional
	If True, the user is ask for the correction of 
	filenames which do not exist
  info : bool, optional
	If True, some additional information is printed
  
  Default Error and Exception Handling: 
	print the error and continue
  '''
  
  #--- Input/Output Options ---#
  
  # Look for the input file
  setattr(thismodule,'filename',check_if_exists(filename, 
				     what='filename for the input file', 
				     interactive=interactive,
				     error=error))
  
  # Check the input type for correctness
  if itype not in itypes:
    error('Invalid input file format (choose from "%s")\n' % 
	  '", "'.join(itypes))
  
  
  fid_base = os.path.splitext(filename)[0]
    
  if outputname is None:
    setattr(thismodule,'outputname',fid_base)
  
  # Check the output types for correctness
  if otype is None:
    setattr(thismodule,'otype',['h5'])
  elif not isinstance(otype,list):
    setattr(thismodule,'otype',[otype])
  if not all(i in otypes for i in otype):
    error('Invalid output file formats (choose from "%s")\n' % 
	  '", "'.join(otypes))
  
  # Check if h5py is installed
  if 'h5' in otype:
    try: 
      import h5py
    except ImportError:
      error('ImportError: The module h5py is not installed!\n')  

  #--- Grid-Related Options ---#  
  
  # Check for the compability of the vector option and the output formats
  if vector is not None and ('cb' in otype or 
	      'am' in otype or 'hx' in otype):
    error('For a vectorized grid, only HDF5 \n' +
    'is available as output format. Type: --otype=h5\n')
    #'The Gaussian cube file output format does not allow\n\t'+
    #'a vectorized grid, i.e., "-v" and "-t cb" are mutually exclusive.\n')
  
  # Look for the grid input file
  if grid_file is not None: 
    setattr(thismodule,'grid_file',check_if_exists(grid_file, 
					what='filename for the grid input file', 
					interactive=interactive,
					error=error))

  #--- Computational Options ---#
  
  if not isinstance(numproc,int):
    error('The number of processes (--numproc) has to be an integer value.\n')
  
  # Check the files specified by --calc_mo or --mo_set for existance
  if calc_mo != False and mo_set != False:
    error('Please choose --calc_mo OR --mo_set, not both. \n'+
	  '--calc_mo will be done.\n')
  if calc_mo != False:
    if calc_mo.lower() != 'all_mo':
      setattr(thismodule,'calc_mo',
	      check_if_exists(calc_mo,
		       what='filename for the MO list',
		       interactive=interactive))
  if mo_set != False:
    if mo_set.lower() != 'all_mo':
      setattr(thismodule,'mo_set',
	      check_if_exists(mo_set, 
		       what='filename for the MO list', 
		       interactive=interactive))
  
  if not isinstance(all_mo,bool):
    error('The option --all_mo has to be a boolean.\n')
  
  if (drv is not None) and not all(i in drv_options for i in drv):
    error('Invalid derivative option (choose from "%s")\n' % 
	  '", "'.join(drv_options))
  
  #--- Additional Options ---#
  if atom_projected_density is not None and drv is not None:
    error('The derivative of the atom_projected density is not implemented.\n')
    
  if mo_tefd is not None:
    setattr(thismodule,'all_mo',True)  
  
  if mo_tefd is not None and drv is None:
    error('The computation of molecular orbital transition electronic \n' + 
	  'flux density between two orbitals requires the selection of \n' +
	  'the component (e.g. --drv=x)\n')
  
  # The following options cannot be checked before running the main program
  if info:
    string = 'The option %s--center cannot be checked before %s...\n'
    if center_grid is not None:
      display(string % ('--center','reading\nthe input file'))
    if z_reduced_density:
      display(string % ('--z_reduced_density',
				      'creating\nthe grid'))
    if atom_projected_density is not None:
      display(string % ('--atom_projected_density',
				      'reading\nthe input file'))
    if mo_tefd is not None:
      display(string % ('--mo_tefd',
				      'reading\nthe input file'))
  
  return True

def check_if_exists(fid, what='',interactive=False,error=IOError,
		    display=sys.stdout.write):
  '''Checks if a file exists
  
  :Parameters:
  
  fid : string
	Filename of the requested file.
  what : string, optional
        Description of the file.
  interactive : string, optional
        If True, ask for user interaction.
  error : exception, optional
        Specifies the error which shall be raised.
  
  :Returns:
  
  fid : string
	Filename of the requested file.
  '''
  
  while not (isinstance(fid,str) and os.path.exists(fid)): 
    if fid != '':
      display('%s does not exist!\n' % fid)
    if interactive:
      fid = raw_input('Please insert a correct %s: ' % what)
    else:
      error('Insert a correct %s!\n' % what)	
      break 
  return fid


def reset_grid():
  '''Reset the grid parameters.'''
  grid.is_initialized = False
  grid.min_ = [-8.0, -8.0, -8.0]
  grid.max_ = [ 8.0,  8.0,  8.0]
  grid.N_   = [ 101,  101,  101]

# initiating the parser variables
# the names are chosen according to core.init_parser()

#--- Input/Output Options ---
filename	= ''			#:
itype		= 'molden'		#:
outputname	= None			#:
otype		= 'h5'			#:
#--- Grid-Related Options ---
vector		= None			#:
grid_file	= None			#:
center_grid	= None			#:
#--- Computational Options ---
numproc		= 1				#:
mo_set		= False			#:
calc_ao		= False			#:
calc_mo		= False			#:
all_mo		= False			#:
drv			= None			#:
#--- Additional Options ---
z_reduced_density = False	#:
atom_projected_density = None #:
mo_tefd		= None			#:
#--- Options for developers ---
quiet		= False			#:
no_log		= False			#:
no_output	= False			#:
no_slice	= False			#:
random_grid	= False			#:
interactive	= False			#:
#--- Default values for the grid parameters ---
reset_grid()
