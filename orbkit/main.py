#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
'''Module for controlling all computational tasks.'''
'''
ORBKIT
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
  Copyright (C) 2017 Gunter Hermann, Vincent Pohl, Lukas Eugen Marsoner Steinkasserer, Axel Schild, and Jean Christophe Tremblay
  This program comes with ABSOLUTELY NO WARRANTY.
  This is free software, and you are welcome to redistribute it
  under certain conditions. Type '-l' for details.
'''

# Import general modules
import time

# Import orbkit modules
from orbkit import core, grid, extras, read
from orbkit import options, output
import orbkit.display as display_module
from orbkit.display import display,good_bye_message

def run_orbkit(use_qc=None,check_options=True,standalone=False):
  '''Controls the execution of all computational tasks.
  
  **Parameters:**
  
  use_qc : QCinfo, optional
    If not None, the reading of a quantum chemistry output is omitted and
    the given QCinfo class is used for all computational tasks. 
    (See :ref:`Central Variables` in the manual for details on QCinfo.) 
  check_options : bool, optional
    If True, the specified options will be validated. 
  
  **Returns:**
  
  data : type and shape depend on the options.
    Contains orbkit's output. 
    See :ref:`High-Level Interface` in the manual for details.
  '''
  # Set some global variables
  global qc
  
  # Display program information
  display(lgpl_short)
  
  # Check for the correctness of orbkit.options
  if check_options:
    display('Checking orbkit.options...\n')
    options.check_options(display=display,
                          interactive=False,
                          info=True,check_io=(use_qc is None))
  
  # Measurement of required execution time
  t=[time.time()]
  
  # Do we need to read out the info of all MOs?  
  if (options.mo_set or options.calc_mo) is not False: 
    options.all_mo = True
  
  if use_qc is None:  
    # Read the input file
    qc = read.main_read(options.filename,
                        itype=options.itype,
                        all_mo=options.all_mo,
                        spin=options.spin,
                        cclib_parser=options.cclib_parser)
  else:
    # Use a user defined QCinfo class.
    qc = use_qc

  if 'native' in options.otype:
    output.main_output(qc, outputname=options.outputname, otype='native', ftype=options.niotype)
    options.otype.remove('native')
    if not len(options.otype):
      t.append(time.time()) # Final time
      good_bye_message(t)
      return
  
  display('\nSetting up the grid...')
  if options.grid_file is not None: 
    # Read the grid from an external file
    grid.read(options.grid_file)
  elif options.adjust_grid is not None:
    # Adjust the grid to geo_spec 
    extend,step = options.adjust_grid      
    grid.adjust_to_geo(qc,extend=extend,step=step)  
  elif options.random_grid:
    # Create a randomized grid
    grid.random_grid(qc.geo_spec)
  
  # Initialize grid
  grid.grid_init(is_vector=options.is_vector)
  if options.is_vector:
    grid.is_regular = False
  display(grid.get_grid())   # Display the grid
  
  if not grid.is_regular and options.center_grid is not None:
    raise IOError('The option --center is only supported for regular grids.')
  elif options.center_grid is not None:
    atom = grid.check_atom_select(options.center_grid,qc.geo_info,qc.geo_spec,
                                  interactive=True,display=display)
    # Center the grid to a specific atom and (0,0,0) if requested
    grid.center_grid(qc.geo_spec[atom-1],display=display)
  
  if check_options or standalone:
    options.check_grid_output_compatibilty()
    
  t.append(time.time()) # A new time step
  
  # The calculation of all AOs (--calc_ao)
  if options.calc_ao != False:
    data = extras.calc_ao(qc,
                          drv=options.drv,
                          otype=options.otype)
    t.append(time.time()) # Final time
    good_bye_message(t)
    return data
  
  # The calculation of selected MOs (--calc_mo) or 
  # the density formed by selected MOs (--mo_set)
  if (options.mo_set or options.calc_mo) != False: 
    # What should the program do?
    if options.calc_mo != False:
      fid_mo_list = options.calc_mo
    elif options.mo_set != False:
      fid_mo_list = options.mo_set
    
    # Call the function for actual calculation 
    if options.calc_mo != False:
      data = extras.calc_mo(qc,
                            fid_mo_list,
                            drv=options.drv,
                            otype=options.otype)
    elif options.mo_set != False: 
      data = extras.mo_set(qc,
                           fid_mo_list,
                           drv=options.drv,
                           laplacian=options.laplacian,
                           otype=options.otype)
    
    t.append(time.time()) # Final time
    good_bye_message(t)
    return data


  if options.gross_atomic_density is not None:
    rho_atom = extras.gross_atomic_density(options.gross_atomic_density,qc,
                                           drv=options.drv)

    if not options.no_output:
      output_written = output.main_output(rho_atom,
                                          qc,
                                          outputname=options.outputname,
                                          otype=options.otype)
    t.append(time.time())
    good_bye_message(t)
    return rho_atom
  
  t.append(time.time()) # A new time step
  
  # Compute the (derivative of the) electron density 
  if options.no_slice:
    data = core.rho_compute_no_slice(qc,
                                     drv=options.drv,
                                     laplacian=options.laplacian,
                                     return_components = False)
  
  else:
    data = core.rho_compute(qc,
                            drv=options.drv,
                            slice_length=options.slice_length,
                            laplacian=options.laplacian,
                            numproc=options.numproc)
  if options.drv is None:
    rho = data
  elif options.laplacian:
    rho, delta_rho, laplacian_rho = data    
  else:
    rho, delta_rho = data
    
  t.append(time.time()) # A new time step

  # Generate the output requested 
  if not options.no_output:
    if not (options.drv is not None or options.laplacian):
      plt_data = rho
      datalabels = 'rho'
    else:
      plt_data = [rho]
      datalabels = ['rho']
      if options.drv is not None:
        plt_data.extend(delta_rho)    
        datalabels.extend(['d/d%s of %s' % (ii_d,'rho') for ii_d in options.drv])
      if options.laplacian:
        plt_data.append(laplacian_rho)
        datalabels.append('laplacian of rho')
    output.main_output(plt_data,qc,outputname=options.outputname,
                       otype=options.otype,datalabels=datalabels)
    
  t.append(time.time()) # Final time
  
  good_bye_message(t)
  
  # Return the computed data, i.e., rho for standard, and (rho,delta_rho)  
  # for derivative calculations. For laplacian (rho,delta_rho,laplacian_rho) 
  return data
  # run_orbkit 

def init(reset_display=True):  
  ''' Resets all :mod:`orbkit.options` and :mod:`orbkit.display`. 
  '''
  try:
    from importlib import reload # >= Python3.4
  except ImportError:
    try: 
      from imp import reload # <= Python3.3
    except ImportError:
      pass # Python2.X
    
  reload(options)
  if reset_display:
    reload(display_module)

def run_standalone():
  '''Starts orbkit as a standalone program using parser options (:mod:`orbkit.core.init_parser`).
  '''  
  # Call the parser
  options.init_parser()
  
  # Call the main loop
  return run_orbkit(check_options=False,standalone=True)
