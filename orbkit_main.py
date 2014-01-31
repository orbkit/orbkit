#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

'''
orbkit
Axel Schild (axel.schild [at] fu-berlin.de)
Gunter Hermann
Vincent Pohl

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
  Copyright (C) 2014 Axel Schild, Gunter Hermann, and Vincent Pohl
  This program comes with ABSOLUTELY NO WARRANTY.
  This is free software, and you are welcome to redistribute it
  under certain conditions. Type '-l' for details.
'''


# Import general modules
import os
import resource
import time
from scipy import integrate

# Import orbkit modules
import orbkit_core as core
import orbkit_grid as grid
import orbkit_extras as extras 
import orbkit_output 
import orbkit_integrate

class cOptions:
  #--- Class initiating the several options if ---------------------------------
  #--- orbkit_main is called as a function ----------------------------------------
  def __init__(self):
    # initiating the parser variables
    # the names are chosen according to core.init_parser()
    self.filename	= None
    self.outputname	= None
    self.no_display	= True
    self.numproc	= 4
    self.all_MO		= False
    self.hx_network	= None
    self.hdf5		= None
    self.calc_mo	= False
    self.mo_list	= False
    self.gaussian	= None
    self.gamess		= None
    self.csv_grid	= False
    self.vector_grid	= False
    self.delta_slice	= 4e4
    self.amira		= False
    self.a_num		= None
    self.discard_ao	= None
    self.save_ao	= None
    self.rho_creator	= False
    self.center		= None 
    self.quiet		= False
    self.no_log		= False
    self.reduced_density= False
    
    # Options for developers    
    self.no_slice	= False
    self.no_output	= False
    
    #--- Default values for the grid parameters ---
    self.reset_grid()
    
  def reset_grid(self):
    grid.min_ = [-8.0, -8.0, -8.0]
    grid.max_ = [ 8.0,  8.0,  8.0]
    grid.N_   = [ 100,  100,  100]

def tForm(string,T,extra=''):
  #--- Function for the formating of the time output string --------------------
  t_diff = int(round(T))
  tF = {}
  tF['str'] = string
  tF['extra'] = extra
  tF['min'] = t_diff/60
  tF['sec'] = t_diff%60
  tF['h']   = tF['min']/60
  tF['min'] = tF['min']%60
  if tF['h'] == 0:
    if tF['min'] == 0: 
      return ('\n%(str)s took %(sec).3fs%(extra)s.\n' % 
	      {'str':string, 'sec': T, 'extra':extra})
    else: 
      return ('\n%(str)s took %(min)dmin and %(sec)ds%(extra)s.\n' % tF)
  else: return ('\n%(str)s took %(h)dh, %(min)dmin and %(sec)ds%(extra)s.\n' % tF)
  #--- tForm ---

def main_read(filename):
  #--- Function to read the input file -----------------------------------------
  orbkit_output.display('Opened \n\t%s\n' % filename)
  
  # Do we need to read out the info of all MOs?  
  if (core.options.mo_list or core.options.calc_mo) != False: 
    core.options.all_MO = True
  
  # What kind of input file has to be read?
  if core.options.gaussian:
    orbkit_output.display('Loading gaussian file...')
    if ('.log' in filename[-4:]):
      geo_spec, geo_info, ao_spec, mo_spec = core.read_gaussian_log(filename)
    elif ('.FChk' in filename[-5:]):
      geo_spec, geo_info, ao_spec, mo_spec = core.read_gaussian_fchk(filename)
    else:
      core.parser.error('Please choose a valid gaussian output file!' + 
					     ' (File extension: .log or .FChk)')
  elif core.options.gamess:
    orbkit_output.display('Loading Gamess file...')
    geo_spec, geo_info, ao_spec, mo_spec = core.read_gamess(filename, 
						      all_MO=core.options.all_MO)
  else:
    orbkit_output.display('Loading molden file...')
    geo_spec, geo_info, ao_spec, mo_spec = core.read_molden(filename, 
						      all_MO=core.options.all_MO)

  # Return required data
  return geo_spec, geo_info, ao_spec, mo_spec 
  #--- main_read ---

def main_calc():
  #--- Main program for all computational tasks --------------------------------
  # Set some global variables
  global cancel, geo_spec, geo_info, ao_spec, mo_spec, rho
  
  cancel = False	# Initialize variable for canceling the program
  
  orbkit_output.display(lgpl_short)
  
  # Measurement of required execution time
  t=[time.time()]

  # Read the input file
  geo_spec, geo_info, ao_spec, mo_spec = main_read(core.options.filename)

  # Read grid from selected *.csv file if requested
  if core.options.csv_grid != False: 
    grid.grid_csv(core.options.csv_grid)
  
  # Initialize grid
  grid.grid_init()
  orbkit_output.display('\nSetting up the grid...')
  grid.grid_display()	  # Display the grid
  # Center the grid to a specific atom and (0,0,0) if requested
  grid.call_center_grid(geo_info,geo_spec,core.options)
  
  # Initialize vector grid
  if core.options.vector_grid:
    grid.grid2vector()
  #if core.options.random:
  #grid.random_grid(geo_spec)
  
  t.append(time.time())	# A new time step
  
  if cancel: return 1	# Cancel the program if requested
  
  # The calculation of selected MOs (--calc_mo) or 
  # the density formed by selected MOs (--mo_list)
  if (core.options.mo_list or core.options.calc_mo) != False: 
    # What should the program do?
    if core.options.calc_mo != False:
      fid_List = core.options.calc_mo
      calc_mo  = True
    elif core.options.mo_list != False:
      fid_List = core.options.mo_list
      calc_mo  = False
      
    # Does the requested file exist?
    if fid_List != 'ALL_MO':
      while (os.path.exists(fid_List)) == False: 
	fid_List = raw_input(fid_List + ' does not exist: ')    
    
    # Call the function for actual calculation
    extras.mo_select(geo_spec, geo_info, ao_spec, 
		       mo_spec, fid_List, calc_mo=calc_mo)  
    
    t.append(time.time())	# Final time
    orbkit_output.display(tForm('The calculation took',t[-1]-t[0]))
    orbkit_output.display('Thank you. Good bye.')
    return 0
  
  t.append(time.time())	# A new time step
  
  # Compute the electron density
  if core.options.no_slice:
    rho = core.rho_compute_no_slice(geo_spec, geo_info, ao_spec, mo_spec,
				    vector=core.options.vector_grid)
  else:
    if not core.options.vector_grid:
      rho = core.rho_compute(geo_spec, geo_info, ao_spec, mo_spec)
    else:
      rho = core.rho_compute_vector(geo_spec, geo_info, ao_spec, mo_spec, 
				   delta_slice=core.options.delta_slice)
  
  # Compute the reduced electron density if requested
  if core.options.reduced_density:
    orbkit_output.display('\nReducing the density with respect to the z-axis.\n')
    rho = integrate.simps(rho, grid.x, dx=grid.delta_[0], axis=0, even='avg')
    rho = integrate.simps(rho, grid.y, dx=grid.delta_[1], axis=0, even='avg')
  
  t.append(time.time())	# A new time step

  if cancel: return 1	# Cancel the program if requested

  # Generate the output requested
  if not core.options.no_output and not core.options.vector_grid:
    orbkit_output.main_output(rho,geo_info,geo_spec,mo_spec,core.options.outputname)
  if core.options.vector_grid:
    orbkit_output.amira_creator_vector(rho,core.options.outputname)
  t.append(time.time())	# Final time
  
  ram_requirement = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
  ram_req = '\nand required %.2f MB of RAM' % (ram_requirement/1000.)
  orbkit_output.display(tForm('The calculation',t[-1]-t[1],extra=ram_req))
  orbkit_output.display('Thank you. Good bye.')
    
  return rho
  #--- main_calc ---

def init():
  #--- Initiate orbkit_main for the case if ---------------------------------------
  #--- orbkit_main is called as a function ----------------------------------------
  # Set some global variables  
  global options
  
  # Initiate the parser options
  options = cOptions()

  return 0
  #--- init ---
  
def main():
  #---- MAIN LOOP ----
  #---- This function calls all funtions ----
  
  # Remove old .log file
  orbkit_output.rm_log()
  
  # Set the options for the execution
  core.options = options
  
  # Call the main function performing all computations 
  main_calc()
    
  return 0
  #--- main loop ---

if __name__ == '__main__':
  #--- The following comands are executed if the program is --------------------
  #--- started as a standalone program -----------------------------------------
  
  # Initiate the console output
  textbuffer = ''
  
  # Call the parser
  core.init_parser()
  
  # Set the parser options
  options = core.options
  
  # Call the main loop
  main()
