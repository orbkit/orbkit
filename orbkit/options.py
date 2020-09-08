# -*- coding: iso-8859-1 -*-
'''Module containing and processing all orbkit options.'''

lgpl = '''ORBKIT
Gunter Hermann, Vincent Pohl, Lukas Eugen Marsoner Steinkasserer, Axel Schild, and Jean Christophe Tremblay

Institut fuer Chemie und Biochemie, Freie Universitaet Berlin, 14195 Berlin, Germany

This file is part of ORBKIT.

ORBKIT is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or any later version.

ORBKIT is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with ORBKIT.  If not, see <http://www.gnu.org/licenses/>.
'''

lgpl_short = '''This is ORBKIT.
Copyright (C) 2017 Gunter Hermann, Vincent Pohl, Lukas Eugen Marsoner Steinkasserer, Axel Schild, and Jean Christophe Tremblay. 
This program comes with ABSOLUTELY NO WARRANTY. 
This is free software, and you are welcome to redistribute it 
under certain conditions. Type '-l' for details.
'''

import os
import sys
from copy import deepcopy
thismodule = sys.modules[__name__]

from orbkit import grid

available = [
  'filename','itype','cclib_parser','outputname','otype',
  'numproc','mo_set','calc_ao','all_mo','calc_mo','spin','drv','laplacian',
  'slice_length','is_vector','grid_file','adjust_grid','center_grid','random_grid',
  'gross_atomic_density',
  'quiet','no_log','no_output','no_slice','interactive', 'test'
  ]

itypes = ['auto',
          'molden',
          'aomix',
          'gamess', 
          'gaussian.log', 
          'gaussian.fchk',
          'wfn',
          'wfx',
          'cclib',
          'native',
          'orbkit.dump']                        #: Specifies possible input types.

niotypes = ['npz',
            'hdf5']                         #: Specifies file format for native io

otypes = ['h5', 'hdf5', 
          'npz', 
          'cube', 'cb', 
          'cube.gz', 'cb.gz', 
          'am', 
          'hx', 
          'vmd', 
          'mayavi',
          'native',
          'auto'] #: Specifies possible output types.

drv_options = ['None','x','y','z',
               'xx','yy','zz','x2','y2','z2',
               'xy','yx','xz','zx','yz','zy']     #: Specifies possible derivative variables.

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
  group = optparse.OptionGroup(parser, "Input/Output Options")
  group.add_option("-i", "--input", dest="filename",metavar="INPUT",
                      default='', type="string",nargs=1,
                      help="input file")
  group.add_option("-e", "--itype", dest="itype",
                      default='auto', type="choice",choices=itypes,
                      help="input type: '" + "', '".join(itypes) + 
                      "' [default: '%default']")
  group.add_option("--niotype", dest="niotype",
                      default='npz', type="choice",choices=niotypes,
                      help="input type: '" + "', '".join(niotypes) + 
                      "' [default: '%default']")
  group.add_option("--cclib_parser",dest="cclib_parser",
                      type="string",
                      help='''if '--itype=cclib', this argument determines what
                      cclib.parser will be used, e.g., 'Gaussian' or 'GAMESS'.''')
  group.add_option("-o", "--output",dest="outputname",
                      type="string",
                      help='''name of the output file 
                      [default: base name of INPUT]''')
  group.add_option("-t", "--otype", dest="otype",
                      type="choice", action="append", choices=otypes,
                      help='''output formats (multiple calls possible):  
                      '{0}' or '{1}' (HDF5 file),                        
                      '{2}' (Compressed numpy file),                             
                      '{3}' or '{5}' (Gaussian cube file),                     
                      '{7}' (VMD network),                           
                      '{8}' (ZIBAmiraMesh file),                     
                      '{9}' (ZIBAmira network),                      
                      '{10}' (simple interactive Mayavi interface) 
                      '{11}' (determine from OUTPUTNAME) 
                      [default: '{11}' if OUTPUTNAME has file extension 
                      else '{0}']'''.format(*otypes))
  parser.add_option_group(group)
  
  group = optparse.OptionGroup(parser, "Computational Options")
  group.add_option("-p", "--numproc",dest="numproc",
                      default=1, type="int",
                      help='''number of subprocesses to be started 
                      during the execution [default: %default]''')
  group.add_option("--mo_set",dest="mo_set",
                      default=[], type="string",action="append",
                      help='''read the plain text file MO_SET containing row 
                      vectors of molecular orbital indeces (delimiter=' ', 
                      Integer numbering or MOLPRO's symmetry numbering) 
                      and compute the electron density 
                      using exclusively those orbitals'''.replace('  ','').replace('\n',''))  
  group.add_option("--calc_ao",dest="calc_ao",
                      default=False,action="store_true", 
                      help="calculate and save all AOs.")  
  group.add_option("--calc_mo",dest="calc_mo",
                      default=[], type="string", action="append",
                      help=('''calculate and save the MOs specified in the 
                      plain text file CALC_MO by the indices (delimiter=' ') 
                      (Type 'all_mo' to store all occupied and virtual
                      orbitals)''').replace('  ','').replace('\n','')) 
  group.add_option("--all_mo",dest="all_mo",
                      default=False, action="store_true", 
                      help='''take into account all (occupied and virtual) MOs 
                      for all computations''')
  group.add_option("--spin",dest="spin",
                      default=None, type=spin, choices=['alpha','beta'],
                      help='''consider only `alpha` or `beta` molecular orbitals
                      for the computations. Only available for unrestricted
                      calculations.'''.replace('  ','').replace('\n',''))
  group.add_option("-d", "--drv",dest="drv",choices=drv_options,
                      type="choice",action="append",
                      help=('''compute the analytical derivative of the requested
                      quantities with respect to DRV, i.e., 'x', 'y', and/or 'z'. 
                      For 2nd derivatives, specify the respective combinations
                      , e.g., 'xx' or 'yz'. (multiple calls possible)'''
                      ).replace('  ','').replace('\n',''))
  group.add_option("--laplacian",dest="laplacian",
                      default=False, action="store_true",
                      help='''compute the analytical laplacian of the density
                      or the specified mo_set, respectively.
                      ''')
  parser.add_option_group(group)
  
  group = optparse.OptionGroup(parser, "Grid-Related Options")
  group.add_option("--adjust_grid",dest="adjust_grid", 
                      type="float",nargs=2,default=[5,0.5],
                      help=('''create a grid using a spacing of X a_0 and having 
                      the size of the molecule plus D a_0 in each direction, 
                      e.g., --adjust_grid=D X [default: --adjust_grid=5 0.5]'''
                      ).replace('  ','').replace('\n',''))
  group.add_option("--grid", dest="grid_file",
                      type="string",
                      help='''read the grid from the plain text file GRID_FILE''')    
  group.add_option("--random_grid", dest="random_grid",
                      default=False, action="store_true",  
                      help=optparse.SUPPRESS_HELP)
  group.add_option("--center", dest="center_grid",
                      metavar="ATOM",type="int",
                      help='''center with respect to the origin and the 
                      atom number ATOM (input order)''')
  group.add_option("-s", "--slice_length",dest="slice_length",
                      default=1e4, type="int",
                      help=('''specify how many grid points are computed at once 
                      (per subprocess).''').replace('  ','').replace('\n',''))
  group.add_option("-v", "--vector",dest="is_vector",
                      default=False, action="store_true", 
                      help=('''store the output in a vector format.''')
                      )   
  parser.add_option_group(group)
  group = optparse.OptionGroup(parser, "Additional Options")
  group.add_option("--gross_atomic_density",dest="gross_atomic_density",
                      metavar="INDEX",action="append",type="int",
                      help='''compute the atom-projected electron density with
                      respect to atom INDEX (multiple calls possible)''')
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
  
  if kwargs.otype is None:
    kwargs.otype = ['auto']
  
  for i,j in vars(kwargs).items():
    setattr(thismodule,i,j)
  
  # Check the options for compatibility and correctness
  check_options(error=parser.error,
                interactive=interactive,
                info=False,
                check_io=(not len(args)))
  
  
  if len(args) and args[0] == 'test':
    from orbkit.test import test
    test()

  return
  # init_parser 

def raise_error(string,error=IOError):
  if hasattr(thismodule,'parser'):
    error = parser.error
  raise error(string)

def print_message(string):
  print(string)

def check_options(error=raise_error,display=print_message,
                  interactive=False,info=True,check_io=True):
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
  if check_io:
    # Look for the input file

    setattr(thismodule,'filename',check_if_exists(filename, 
                          what='filename for the input file', 
                          interactive=interactive,
                          error=error))
    # Check the input type for correctness
    if itype not in itypes:
      error('Invalid input file format (choose from "%s")\n' % 
          '", "'.join(itypes))
    if itype == 'cclib' and cclib_parser is None:
      error('The input type cclib requires the specification of parser, ' + 
            'e.g., --cclib_parser=Gaussian')

    if niotype not in niotypes:
      error('Unsupported format for native io (choose from "%s")\n' %
            '", "'.join(niotypes))
    if niotype == 'hdf5':
      try:
        __import__('h5py')
      except ImportError:
        error('External IO to HDF5 file was requested but no\n' +
              'HDF5 module could be found.')
    
    fid_base,ext = os.path.splitext(filename)
    
    if outputname is None:
      setattr(thismodule,'outputname',fid_base)
    else:
      outpath = os.path.dirname(outputname.split('@')[0])
      if not (outpath == '' or os.path.exists(outpath)):
        error('Output path "%s" does not exist!' % outpath)
      
  
    # Check the output types for correctness
    
    if otype is None:
      setattr(thismodule,'otype',[])
    elif ('auto' in otype and os.path.splitext(outputname)[1] is not None 
          and os.path.splitext(outputname)[1][1:] not in otypes):
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
  
  # Look for the grid input file
  if grid_file is not None: 
    setattr(thismodule,'grid_file',check_if_exists(grid_file, 
                        what='filename for the grid input file', 
                        interactive=interactive,
                        error=error))
  
  if adjust_grid is not None:
    if (not isinstance(adjust_grid,(list,tuple)) or 
       (len(adjust_grid) != 2) or 
       (not isinstance(adjust_grid[0],(int,float))) or 
       (not isinstance(adjust_grid[1],(int,float)))
       ):
      error('The grid parameter (--adjust_grid), has to be a list containing '
            'containing two floats.\n')
    elif adjust_grid[1] == 0:
      error('The grid spacing (second value in --adjust_grid) cannot be zero.\n')
  
  #--- Computational Options ---#
  
  if not isinstance(numproc,int):
    error('The number of processes (--numproc) has to be an integer value.\n')
  
  # Check the files specified by --calc_mo or --mo_set for existance
  def check_mo(attr):
    data = getattr(thismodule,attr)
    if not data:
      setattr(thismodule,attr,False)
      return False
    if isinstance(data,int):
      data = str(data)
    if isinstance(data,str):
      data = [data]
    try:
      for d in data:
        d = str(d)
        if not (',' in d.lower() or ':' in  d.lower()):        
          i = deepcopy(d)
          if i != 'homo' and i != 'lumo':
            for r in ['homo','lumo','-','+']:
              i = i.replace(r,'')
            int(i.split('.')[0])
    except ValueError: 
      if len(data) == 1:
        data = data[0]
        if not any([i != data.lower() for i in ['all_mo','occupied','unoccupied','virtual']]):   
          setattr(thismodule,attr,
                  check_if_exists(data,
                  what='filename for the MO list',
                  interactive=interactive))
        else:
          setattr(thismodule,attr,data)
      else:
        display('You have called `%s` multiple times. So, you have\n' % attr + 
                'to give the molecular orbital labels explicitly, i.e.,\n' + 
                'no filenames and no usage of the keyword `all_mo`.\n\n')
        error('Entry `%s` is not valid!' % d)
    return True
  
  i = check_mo('calc_mo')
  j = check_mo('mo_set')
  
  if i and j:
    error('Please choose --calc_mo OR --mo_set, not both. \n'+
        '--calc_mo will be done.\n')
  
  if not isinstance(all_mo,bool):
    error('The option --all_mo has to be a boolean.\n')
  
  if spin is not None and not (spin == 'alpha' or spin == 'beta'):
    error('The option --spin has to be `alpha` or `beta`.\n')
  
  if (drv is not None) and not all(i in drv_options for i in drv):
    error('Invalid derivative option (choose from "%s")\n' % 
        '", "'.join(drv_options))
  
  if laplacian:
    if not (drv is None or drv == ['xx','yy','zz'] or drv == ['x2','y2','z2']):
      display('Note: You have set the option --laplacian and specified values\n' +
              'for --drv. Both options are not compatible.\n\n' +
              'The options have been changed to -dxx -dyy -dzz.\n')
    setattr(thismodule,'drv', ['xx','yy','zz'])
  
  #--- Additional Options ---#
  if gross_atomic_density is not None and drv is not None:
    error('The derivative of the gross atomic density is not implemented.\n')
      
  # The following options cannot be checked before running the main program
  if info:
    string = 'The option %s--center cannot be checked before %s...\n'
    if center_grid is not None:
      display(string % ('--center','reading\nthe input file'))
    if gross_atomic_density is not None:
      display(string % ('--gross_atomic_density',
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
      fid = input('Please insert a correct %s: ' % what)
    else:
      error('Insert a correct %s!\n' % what)
      break 
  return fid

def check_grid_output_compatibilty(error=raise_error): 
  if not grid.is_regular and ('cube'   in otype or 
                              'cb'     in otype or
                              'vmd'    in otype or
                              'am'     in otype or 
                              'hx'     in otype or 
                              'mayavi' in otype):
    error('For a non-regular vector grid, only HDF5 ' +
    'is available as output format. Choose: --otype=h5\n')

# initiating the parser variables
# the names are chosen according to core.init_parser()

#--- Input/Output Options ---
filename        = ''            #: Specifies input file name. (str)
itype           = 'auto'        #: Specifies input file type. See :data:`itypes` for details. (str) 
niotype         = 'npz'         #: Specifies output filetype for native io
cclib_parser    = None          #: If itype is 'cclib', specifies the cclib.parser. (str)
outputname      = ''            #: Specifies output file (base) name. (str)
otype           = 'auto'        #: Specifies output file type. See :data:`otypes` for details. (str or list of str or None)
#--- Computational Options ---
numproc         = 1             #: Specifies number of subprocesses for multiprocessing. (int)
mo_set          = False         #: Specifies molecular orbitals used for density calculation. (filename or list of indices)
calc_ao         = False         #: If True, all atomic orbitals will be computed and saved.
calc_mo         = False         #: Specifies which molecular orbitals will be calculated. (filename or list of indices)
all_mo          = False         #: If True, all molecular orbitals will be computed. (bool)
spin            = None          #: If not None, exclusively 'alpha' or 'beta' molecular orbitals are taken into account. (None,'alpha', or 'beta')
drv             = None          #: Specifies derivative variables. (list of str)
laplacian       = False         #: If True, computes the laplacian of the density or of the mo_set. (bool)
#--- Grid-Related Options ---
slice_length    = 1e4           #: Specifies the number of points per subprocess. (int)
vector          = None          #  This option is only present because of backward compatibility
is_vector       = False         #: If True, vector grid is used for the output. (bool)
grid_file       = None          #: Specifies file to read grid from. (filename)
adjust_grid     = None          #: If not None, create a grid using a spacing of X a_0 and having the size of the molecule plus D a_0 in each direction. (list: [D, x])
center_grid     = None          #: If not None, grid is centered to specified atom and origin. (int) 
random_grid     = False         #: If True, creates random grid around atom positions. (bool)
#--- Additional Options ---
gross_atomic_density = None     #: Computes the gross atomic electron density with respect to specified atom. (int or list of int)
#--- Options for Advanced Users ---
quiet           = False         #: If True, omits terminal output. (bool)
no_log          = False         #: If True, omits logfile output. (bool)
no_output       = False         #: If True, omits creation of output. (bool)
no_slice        = False         #: If True, omits slicing of the grid. (bool)
interactive     = False         #: If True, asks user to select unclarified options. (bool)
#--- Default values for the grid parameters ---
grid.reset_grid()
