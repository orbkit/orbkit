#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
'''Module for controlling all computational tasks.'''
'''
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
  Copyright (C) 2014 Gunter Hermann, Vincent Pohl, and Axel Schild
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
from orbkit import core, grid, extras, read, output, integrate
from orbkit import options, display
from orbkit.qcinfo import QCinfo

def tForm(string,T,extra=''):
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
  # tForm 

def good_bye_message(t):
  ram_requirement = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
  ram_req = '\nand required %.2f MB of RAM' % (ram_requirement/1000.)
  display.display(tForm('The calculation',t[-1]-t[0],extra=ram_req))
  display.display('Thank you. Good bye.')

def main():
  ''' Controls the execution of all computational tasks.
  '''
  # Set some global variables
  global qc
  global rho, delta_rho
  
  # initialize quantum chemical information variable
  qc = QCinfo()
  
  # Display program information
  display.display(lgpl_short)
  
  # Measurement of required execution time
  t=[time.time()]
  
  # Do we need to read out the info of all MOs?  
  if (options.mo_set or options.calc_mo) is not False: 
    options.all_mo = True
  
  # Read the input file
  qc = read.main_read(options.filename,
                      itype=options.itype,
                      all_mo=options.all_mo)

  display.display('\nSetting up the grid...')
  if options.grid_file is not None: 
    if grid.read(options.grid_file) and (options.vector is None):
      options.vector = core.dvec
  
  if options.random_grid:
    grid.random_grid(qc.geo_spec)
    if (options.vector is None):
      options.vector = core.dvec
  else:
    # Initialize grid
    grid.grid_init(is_vector=(options.vector is not None))
  
  display.display(grid.get_grid())   # Display the grid

  if options.vector and options.center_grid is not None:
    raise IOError('The option --center is only supported for regular grids.')
  elif options.center_grid is not None:
    atom = options.center_grid
    if not((isinstance(atom, int)) and (0 < atom <= len(qc.geo_spec))):
      display.display('Not a Valid atom number for centering the grid')
      display.display('Coose a valid index:')
      for i,j in enumerate(qc.geo_info): display.display('\t%s\t%d' % (i[0], j+1))
      if options.interactive:
        while not((isinstance(atom, int)) and (0 < atom <= len(qc.geo_spec))):
          atom = raw_input('Please insert a correct index: ')
      else: raise IOError('Insert a correct filename for the MO list!')

    # Center the grid to a specific atom and (0,0,0) if requested
    grid.center_grid(qc.geo_spec[atom])

  if options.vector is not None:
    info = 'vectorized'
  else:
    info = 'regular'
  display.display('The computations will be carried out applying ' +
                'a %s grid...' % info)
  
  t.append(time.time()) # A new time step
  
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
      data = extras.calc_mo(qc.geo_spec, qc.geo_info, qc.ao_spec, qc.mo_spec, 
                fid_mo_list, drv=options.drv, vector=options.vector, 
                otype=options.otype)
    elif options.mo_set != False: 
      data = extras.mo_set(qc.geo_spec, qc.geo_info, qc.ao_spec, qc.mo_spec, 
                fid_mo_list, drv=options.drv, vector=options.vector, 
                otype=options.otype)
    
    t.append(time.time()) # Final time
    good_bye_message(t)
    return data


  if options.atom_projected_density is not None:
    atom = options.atom_projected_density
    rho_atom, mulliken_charge = extras.compute_mulliken_charges(atom,
                        qc.geo_info,qc.geo_spec,qc.ao_spec,qc.mo_spec,
                        is_vector=(options.vector is not None))
    
    if not options.no_output:
      from numpy import array
      fid = '%s.h5' % options.outputname
      display.display('\nSaving to Hierarchical Data Format file (HDF5)...' +
                '\n\t%(o)s' % {'o': fid})
      HDF5_File = output.hdf5_open(fid,mode='w')
      data = {'geo_info':array(qc.geo_info), 'geo_spec':array(qc.geo_spec),
            'atom_projected_density':rho_atom, 'atom':array(atom),
            'mulliken_charge':mulliken_charge,  
            'x':grid.x, 'y':grid.y, 'z':grid.z}
      output.hdf5_append(data,HDF5_File,name='')
      HDF5_File.close()
    
    t.append(time.time())
    good_bye_message(t)
    return rho_atom
  
  
  if options.mo_tefd is not None:
    mos = options.mo_tefd
    ao_list = core.ao_creator(qc.geo_spec,qc.ao_spec,
                              is_vector=(options.vector is not None))
    mo_tefd = []
    index = []
    for i,j in mos:
      mo_tefd.append([])
      index.append([])
      for ii_d in options.drv:
        display.display('\nMO-TEFD: %s->%s %s-component'%(i,j,ii_d))
        tefd = extras.mo_transition_flux_density(i,j,
                                        qc.geo_spec,qc.ao_spec,qc.mo_spec,
                                        drv=ii_d,ao_list=ao_list,
                                        is_vector=(options.vector is not None))
        mo_tefd[-1].append(tefd)
        index[-1].append('%s->%s:%s'%(i,j,ii_d))
    
    if not options.no_output:
      from numpy import array
      fid = '%s.h5' % options.outputname
      display.display('\nSaving to Hierarchical Data Format file (HDF5)...' +
                      '\n\t%(o)s' % {'o': fid})
      HDF5_File = output.hdf5_open(fid,mode='w')
      data = {'geo_info':array(qc.geo_info), 'geo_spec':array(qc.geo_spec),
              'mo_tefd:info':array(index), 'mo_tefd':array(mo_tefd), 
              'x':grid.x, 'y':grid.y, 'z':grid.z}
      output.hdf5_append(data,HDF5_File,name='')
      HDF5_File.close()
    
    t.append(time.time())
    good_bye_message(t)
    return mo_tefd
  
  t.append(time.time()) # A new time step
  
  # Compute the (derivative of the) electron density 
  if options.no_slice:
    data = core.rho_compute_no_slice(qc.geo_spec, qc.ao_spec, qc.mo_spec, 
                                     drv=options.drv,
                                     is_vector=(options.vector is not None),
                                     return_components = False)
  
  else:
    data = core.rho_compute(qc.geo_spec, qc.ao_spec, qc.mo_spec, 
                            drv=options.drv,
                            vector=options.vector,
                            numproc=options.numproc)
  if options.drv is None:
    rho = data
  else:
    rho, delta_rho = data
  
  # Compute the reduced electron density if requested 
  if options.z_reduced_density:
    if options.vector is not None:
      display.display(
      '\nSo far, reducing the density is not supported for vectorized grids.\n'+ 
      'Skipping the reduction...\n')
    elif options.drv is not None:
      display.display(
      '\nSo far, reducing the density is not supported for the derivative of the density.\n'+ 
      'Skipping the reduction...\n')
    else:
      display.display('\nReducing the density with respect to the z-axis.\n')
      rho = integrate.simps(rho, grid.x, dx=grid.delta_[0], axis=0, even='avg')
      rho = integrate.simps(rho, grid.y, dx=grid.delta_[1], axis=0, even='avg')
  
  t.append(time.time()) # A new time step

  # Generate the output requested 
  if not options.no_output:
    output.main_output(rho,qc.geo_info,qc.geo_spec,
                    outputname=options.outputname,
                    otype=options.otype,
                    data_id='rho',mo_spec=qc.mo_spec,
                    is_vector=(options.vector is not None))
    if options.drv is not None:
      output.main_output(delta_rho,qc.geo_info,qc.geo_spec,
                    outputname=options.outputname,
                    otype=options.otype,
                    data_id='delta_rho',mo_spec=qc.mo_spec,
                    drv=options.drv,
                    is_vector=(options.vector is not None))
      
  t.append(time.time()) # Final time
  
  good_bye_message(t)
    
  # Return the computed data, i.e., rho for standard, and (rho,delta_rho)  
  # for derivative calculations 
  return data
  # main 

def init():  
  ''' Resets all :mod:`orbkit.options` and :mod:`orbkit.display`. 
  '''
  reload(options)
  display.is_initiated = False
  
  return 
  # init 
  

#if __name__ == '__main__':
def run_standalone():
  '''Starts orbkit as a standalone program using parser options (:mod:`orbkit.core.init_parser`).
  '''  
  # Call the parser
  options.init_parser()
  
  # Reset orbkit.display
  display.is_initiated = False

  # Call the main loop
  main()
