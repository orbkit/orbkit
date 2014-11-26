# -*- coding: iso-8859-1 -*-
'''Module containing and processing all orbkit options.'''

lgpl = '''
orbkit
Gunter Hermann, Vincent Pohl, and Axel Schild

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

lgpl_short = '''This is orbkit.
Copyright (C) 2014 Gunter Hermann, Vincent Pohl, and Axel Schild. 
This program comes with ABSOLUTELY NO WARRANTY. 
This is free software, and you are welcome to redistribute it 
under certain conditions. Type '-l' for details.
'''

import os
import sys
from copy import deepcopy
thismodule = sys.modules[__name__]

from orbkit import grid

available = ['filename','itype','outputname','otype',
             'vector','grid_file','center_grid',
             'numproc','mo_set','all_mo','drv',
             'z_reduced_density','atom_projected_density','mo_tefd',
             'quiet','no_log','no_output','no_slice','random_grid',
             'interactive']

#__all__ = ['get_options','check_options','reset_grid',].extend(available)

itypes = ['molden', 
          'gamess', 
          'gaussian.log', 
          'gaussian.fchk']                        #: Specifies possible input types.
otypes = ['h5', 'cb', 'am', 'hx', 'vmd']       #: Specifies possible output types.
drv_options = ['x','y','z']             #: Specifies possible derivative variables.

def get_options():
  '''Returns all possible options and their value.'''
  opts = [(i,globals()[i]) for i in available]
  return dict(opts)

def init_parser():
  '''Initializes parser and processes the options.
  '''
  import optparse
  global parser
  
  def default_if_called(option, opt, value, parser, default=1e4):
    try:
      arg = parser.rargs[0]
      if ((arg[:2] == "--" and len(arg) > 2) or
        (arg[:1] == "-" and len(arg) > 1 and arg[1] != "-")):
        raise ValueError
      value = int(float(arg))
    except (IndexError, ValueError):
      value = int(default)
    setattr(parser.values, option.dest, value)
  
  #optparse.Option.STORE_ACTIONS += ('call_back',)
  usage = 'Usage: %prog [options] -i INPUT'
  parser = optparse.OptionParser(usage=usage,description=lgpl_short) 
  
  parser.add_option("-l", dest="show_lgpl",
                      default=False,action="store_true", 
                      help="show license information and exit")
  parser.add_option("--quiet",dest="quiet",
                      default=False,action="store_true", 
                      help="suppress terminal output")  
  parser.add_option("--no_log",dest="no_log",
                      default=False,action="store_true", 
                      help="suppress output of a INPUT.oklog logfile")
  group = optparse.OptionGroup(parser, "Input/Output Options", 
        '''Comment: So far, orbkit can only handle Cartesian Gaussian basis
        functions. See the manual for details.''')
  group.add_option("-i", "--input", dest="filename",metavar="INPUT",
                      default='', type="string",nargs=1,
                      help="input file")
  group.add_option("--itype", dest="itype",
                      default='molden', type="choice",choices=itypes,
                      help='''input type: ''' + ', '.join(itypes) + 
                      " [default: '%default']")
  group.add_option("-o", "--output",dest="outputname",
                      type="string",
                      help='''name of the output file 
                      [default: base name of INPUT]''')
  group.add_option("-t", "--otype", dest="otype",
                      type="choice", action="append", choices=otypes,
                      help='''output formats (multiple calls possible):  
                      '%s' (HDF5 file), '%s' (Gaussian cube file), 
                      '%s' (ZIBAmiraMesh file), '%s' (ZIBAmira network), 
                      '%s' (VMD network) '''
                      % tuple(otypes) + "[default: 'h5']")
  parser.add_option_group(group)
  group = optparse.OptionGroup(parser, "Computational Options")
  group.add_option("-p", "--numproc",dest="numproc",
                      default=1, type="int",
                      help='''number of subprocesses to be started 
                      during the execution [default: %default]''')
  group.add_option("--mo_set",dest="mo_set",
                      default=False, type="string", 
                      help='''read the plain text file MO_SET containing row 
                      vectors of molecular orbital indeces (delimiter=' ', 
                      Integer numbering or MOLPRO's symmetry numbering) 
                      and compute the electron density 
                      using exclusively those orbitals'''.replace('  ','').replace('\n',''))  
  group.add_option("--calc_ao",dest="calc_ao",
                      default=False, type="string", 
                      help=optparse.SUPPRESS_HELP)
                      #="calculate and save the AOs specified by the indices 
                      #in their selected file (delimiter=' ')") #INCLUDEME  
  group.add_option("--calc_mo",dest="calc_mo",
                      default=False, type="string", 
                      help=('''calculate and save the MOs specified in the 
                      plain text file CALC_MO by the indices (delimiter=' ') 
                      (Type 'all_mo' to store all occupied and virtual
                      orbitals)''').replace('  ','').replace('\n','')) 
  group.add_option("--all_mo",dest="all_mo",
                      default=False, action="store_true", 
                      help='''take into account all (occupied and virtual) MOs 
                      for all computations''')
  group.add_option("-d", "--drv",dest="drv",choices=drv_options,
                      type="choice",action="append",
                      help=('''compute the analytical derivative of the requested
                      quantities with respect to DRV, i.e., 'x', 'y', and/or 'z' 
                      (multiple calls possible)'''
                      ).replace('  ','').replace('\n',''))
  parser.add_option_group(group)
  group = optparse.OptionGroup(parser, "Grid-Related Options")      
  group.add_option("-v", "--vector",dest="vector",
                      action="callback",callback=default_if_called,
                      callback_kwargs={'default': dvec},
                      help=('''perform the computations for a vectorized grid, 
                      i.e., with x, y, and z as vectors. Compute successively 
                      VECTOR points at once per subprocess
                      [default: --vector=%0.0e]''' % dvec
                      ).replace('  ','').replace('\n',''))   
  group.add_option("--grid", dest="grid_file",
                      type="string",
                      help='''read the grid from the plain text file GRID_FILE''')    
  group.add_option("--adjust_grid",dest="adjust_grid", 
                      type="float",nargs=2,
                      help=('''create a grid using a spacing of X a_0 and having 
                      the size of the molecule plus D a_0 in each direction, 
                      e.g., --adjust_grid=D X'''
                      ).replace('  ','').replace('\n',''))
  group.add_option("--random_grid", dest="random_grid",
                      default=False, action="store_true",  
                      help=optparse.SUPPRESS_HELP)
  group.add_option("--center", dest="center_grid",
                      metavar="ATOM",type="int",
                      help='''center with respect to the origin and the 
                      atom number ATOM (input order)''')
  parser.add_option_group(group)
  group = optparse.OptionGroup(parser, "Additional Options")
  group.add_option("--z_reduced_density",dest="z_reduced_density",
                      default=False, action="store_true", 
                      help="reduce the density with respect to the z-axis")
  group.add_option("--atom_projected_density",dest="atom_projected_density",
                      metavar="INDEX",action="append",type="int",
                      help='''compute the atom-projected electron density with
                      respect to atom INDEX (multiple calls possible)''')
  group.add_option("--mo_tefd",dest="mo_tefd", 
                      type="int",nargs=2,action="append",
                      help=('''compute the molecular orbital transition electronic 
                      flux density between the orbitals I and J specify the 
                      requested component with "--drv", e.g., 
                      --mo_tefd=I J --drv=x (multiple calls possible)'''
                      ).replace('  ','').replace('\n',''))
                      
  # The following parser options are hidden 
  group.add_option("--no_slice",dest="no_slice",
                      default=False, action="store_true",
                      help=optparse.SUPPRESS_HELP)
  group.add_option("--no_output",dest="no_output",
                      default=False, action="store_true",
                      help=optparse.SUPPRESS_HELP)
  group.add_option("--not_interactive",dest="interactive",
                      default=True, action="store_false",
                      help=optparse.SUPPRESS_HELP)
  parser.add_option_group(group)

  (kwargs, args) = parser.parse_args()
  
  # Print the licence, if requested
  if kwargs.show_lgpl:
    print(lgpl.replace('\nThis file is part of orbkit.\n',''))
    sys.exit(0)
  
  # Print help if no input file has been set
  if kwargs.filename == '':
    parser.print_help()
    sys.exit(1)
  
  for i,j in vars(kwargs).iteritems():
    setattr(thismodule,i,j)
    #setattr(options,i,j)
  
  # Check the options for compatibility and correctness
  check_options(error=parser.error,
                        interactive=interactive,
                        info=False)
  
  return
  # init_parser 


def check_options(error=sys.stdout.write,display=sys.stdout.write,
            interactive=False,
            info=True):
  '''Checks options for errors.
  
  **Parameters:**
  
    error : function, optional
      Handles the errors.
    display :  function, optional
      Handles the print commands.
    interactive : bool, optional
      If True and a file does not exist, asks the user to insert name of 
      existing file.
    info : bool, optional
      If True, some additional information is printed.

    
  :Default Error and Exception Handling: 
    Prints the errors and continues.
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
  if vector is not None and ('cb' in otype or 'vmd' in otype or 
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
    if not (calc_mo.lower() == 'all_mo' or ',' in calc_mo.lower() or ':' in  calc_mo.lower()):
      try:
        i = deepcopy(calc_mo)
	if i != 'homo' and i != 'lumo':
	  for r in ['homo','lumo','-','+']:
            i = i.replace(r,'')
          int(i.split('.')[0])
      except ValueError: 
        setattr(thismodule,'calc_mo',
                check_if_exists(calc_mo,
                what='filename for the MO list',
                interactive=interactive))
  if mo_set != False:
    if not (mo_set.lower() == 'all_mo' or ',' in mo_set.lower() or ':' in  mo_set.lower()):
      try:
        i = deepcopy(mo_set)
        for r in ['homo','lumo','-','+']:
          i = i.replace(r,'')
        int(i.split('.')[0])
      except ValueError: 
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

def check_if_exists(fid, what='',error=IOError,display=sys.stdout.write,
                interactive=False):
  '''Checks the existence of a file.
  
  **Parameters:**
  
    fid : string
      Specifies filename of the requested file.
    what : string, optional
      Describes the file.
    error : function, optional
      Handles the errors.
    display :  function, optional
      Handles the print commands.
    interactive : bool, optional
      If True and a file does not exist, asks the user to insert name of 
      existing file.
  
  **Returns:**
  
    fid : string
      Specifies filename of the requested file.
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

# initiating the parser variables
# the names are chosen according to core.init_parser()

#--- Input/Output Options ---
filename        = ''            #: Specifies input file name. (str)
itype           = 'molden'      #: Specifies input file type. See :data:`itypes` for details. (str) 
outputname      = None          #: Specifies output file (base) name. (str)
otype           = 'h5'          #: Specifies output file type. See :data:`otypes` for details. (str or list of str)
#--- Computational Options ---
numproc         = 1             #: Specifies number of subprocesses for multiprocessing. (int)
mo_set          = False         #: Specifies molecular orbitals used for density calculation. (filename)
calc_ao         = False         #
calc_mo         = False         #: Specifies which molecular orbitals will be calculated. (filename)
all_mo          = False         #: If True, all molecular orbitals will be computed. (bool)
drv             = None          #: Specifies derivative variables. (list of str)
#--- Grid-Related Options ---
vector          = None          #: If not None, vector grid is used. Specifies number of points per subprocess. (int)
dvec            = 1e4           #(No Option) Specifies the standard value for the points per subprocess. (int)
grid_file       = None          #: Specifies file to read grid from. (filename)
adjust_grid     = None          #: If not None, create a grid using a spacing of X a_0 and having the size of the molecule plus D a_0 in each direction. (list: [D, x])
center_grid     = None          #: If not None, grid is centered to specified atom and origin. (int) 
random_grid     = False         #: If True, creates random grid around atom positions. (bool)
#--- Additional Options ---
z_reduced_density = False       #: If True, reduces the density with respect to the z-axis. (bool)
atom_projected_density = None   #: Computes the atom-projected electron density with respect to specified atom. (int or list of int)
mo_tefd         = None          #: Computes the molecular orbital transition electronic flux density between the orbitals I and J specify the requested component with :data:`orbkit.options.drv`. (list of [I, J])
#--- Options for Advanced Users ---
quiet           = False         #: If True, omits terminal output. (bool)
no_log          = False         #: If True, omits logfile output. (bool)
no_output       = False         #: If True, omits creation of output. (bool)
no_slice        = False         #: If True, omits slicing of the grid. (bool)
interactive     = False         #: If True, asks user to select unclarified options. (bool)
#--- Default values for the grid parameters ---
grid.reset_grid()
