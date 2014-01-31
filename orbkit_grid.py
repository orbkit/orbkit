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

# Import general modules
import csv 
import os

import numpy
from scipy import weave

# Import orbkit modules
import orbkit_core as core
import orbkit_output

def grid_init():
  #--- FUNCTION grid_init ------------------------------------------------------
  #--- Set up the x-, y-, z-grid specified by the global lists: ----------------
  #--- min_ (minimum values), max_ (maximum values), and N_ (# of grid points) -
  #-----------------------------------------------------------------------------
  
  #--- All grid related variables should be globals ---
  global x, y, z, d3r, min_, max_, N_, delta_, grid
  
  #--- Initialize a list for the grid ---
  grid = [[],[],[]]
  delta_ = numpy.zeros((3,1))
  
  #--- Loop over the three dimensions ---
  for ii in range(3):
    if max_[ii] == min_[ii]:
      #--- If min-value is equal to max-value, write only min-value to grid --- 
      grid[ii]   = numpy.array([min_[ii]])
      delta_[ii] = 1
    else:
      #--- Calculate the grid using the input parameters ---
      delta_[ii] = (max_[ii]-min_[ii]) / float(N_[ii] - 1)
      grid[ii] = min_[ii] + numpy.arange(N_[ii]) * delta_[ii]
  
  #--- Write grid ---
  x = grid[0]  
  y = grid[1]  
  z = grid[2]
  d3r = numpy.product(delta_)
  
  return 0
  #--- grid_init ---

def grid2vector():
  #--- FUNCTION grid2vector ----------------------------------------------------
  #--- Convert the regular grid characterized by x-, y-, z-vectors -------------
  #--- to a 3x(Nx*Ny*Nz) grid matrix -------------------------------------------
  #-----------------------------------------------------------------------------
  
  #--- All grid related variables should be globals ---
  global x, y, z
  
  #--- Initialize a list for the grid ---
  grid = numpy.zeros((3,numpy.product(N_)))
  
  grid_code = """
  int count=0;
  
  for (int i=0; i<Nx[0]; i++)
  {
    for (int j=0; j<Ny[0]; j++)
    {
      for (int k=0; k<Nz[0]; k++)
      {
        GRID2(0,count) = x[i];
        GRID2(1,count) = y[j];
        GRID2(2,count) = z[k];
	count += 1;
      }
    }
  }
  """
  weave.inline(grid_code, ['x','y','z','grid'], verbose = 1, support_code = core.cSupportCode.math)
  
  #--- Write grid ---
  x = grid[0,:]  
  y = grid[1,:]  
  z = grid[2,:]
  
  return 0
  #--- grid2vector ---

def sph2cart_vector(r,theta,phi):
  #--- All grid related variables should be globals ---
  global x, y, z
  
  grid = numpy.zeros((3,numpy.product([len(r),len(theta),len(phi)])))
  grid_code = """
  int count=0;

  for (int i=0; i<Nr[0]; i++)
  {
    for (int j=0; j<Ntheta[0]; j++)
    {
      for (int k=0; k<Nphi[0]; k++)
      {
	GRID2(0,count) = r[i] * sin(theta[j]) * cos(phi[k]);
	GRID2(1,count) = r[i] * sin(theta[j]) * sin(phi[k]);
	GRID2(2,count) = r[i] * cos(theta[j]);
	count += 1;
      }
    }
  }
  """
  weave.inline(grid_code, ['r','theta','phi','grid'], verbose = 1, support_code = core.cSupportCode.math)

  #--- Write grid ---
  x = grid[0,:]  
  y = grid[1,:]  
  z = grid[2,:]
  
  return 0
  #--- grid2vector ---

def vector2grid(matrix):
  #--- FUNCTION grid2vector ----------------------------------------------------
  #--- Convert the regular grid characterized by x-, y-, z-vectors -------------
  #--- to a 3x(Nx*Ny*Nz) grid matrix -------------------------------------------
  #-----------------------------------------------------------------------------
  
  #--- All grid related variables should be globals ---
  global x, y, z
  
  #--- Initialize a list for the grid ---
  grid = numpy.zeros((3,numpy.product(N_)))
  
  grid_code = """
  int count=0;
  
  for (int i=0; i<Nx[0]; i++)
  {
    for (int j=0; j<Ny[0]; j++)
    {
      for (int k=0; k<Nz[0]; k++)
      {
        GRID2(0,count) = x[i];
        GRID2(1,count) = y[j];
        GRID2(2,count) = z[k];
	count += 1;
      }
    }
  }
  """
  weave.inline(grid_code, ['x','y','z','grid'], verbose = 1, support_code = core.cSupportCode.math)
  
  #--- Write grid ---
  x = grid[0,:]  
  y = grid[1,:]  
  z = grid[2,:]
  
  return 0
  #--- grid_init ---

def random_grid(geo_spec,N=1e6,scale=0.5):
  
  #--- All grid related variables should be globals ---
  global x, y, z
  geo_spec = numpy.array(geo_spec)
  at_num = len(geo_spec)
  #--- Initialize a list for the grid ---
  grid = numpy.zeros((3,at_num,N))
  
  #--- Loop over the three dimensions ---
  for ii_d in range(3):
    for ii_a in range(at_num):
      grid[ii_d,ii_a,:] = numpy.random.normal(loc=geo_spec[ii_a,ii_d],scale=0.5,size=N)
  
  grid = numpy.reshape(grid,(3,N*at_num))

  #--- Write grid ---
  x = grid[0]  
  y = grid[1]  
  z = grid[2]
  
  return 0
  #--- random_grid ---

def grid_display(quiet=False,start='\t'):
  #--- FUNCTION grid_display ---------------------------------------------------
  #--- Display the current x-, y-, z-grid --------------------------------------
  #-----------------------------------------------------------------------------
  coord = ['x', 'y', 'z']
  display = ''
  for ii in range(3):
    display += ('%(s)s%(c)smin = %(min).2f %(c)smax = %(max).2f N%(c)s = %(N)d ' % 
      {'s': start, 'c': coord[ii], 'min': min_[ii], 'max': max_[ii], 'N': N_[ii]})
    if max_[ii] != min_[ii]:
      # Print the delta values only if min-value is not equal to max-value ---
      display += 'd%(c)s = %(d).3f' % {'c': coord[ii], 'd': delta_[ii]}
    display += '\n'
  
  if not quiet:
    #--- Display the string ---
    orbkit_output.display(display)
  
  return display
  #--- grid_display ---

def grid_csv(csv_grid):
  #--- FUNCTION grid_csv -------------------------------------------------------
  #--- Load the grid paramters from a csv-file (delimiter ",") -----------------
  #-----------------------------------------------------------------------------
  
  #--- All grid related variables should be globals ---
  global  min_, max_, N_
  
  if os.path.exists(csv_grid):
    #--- Read the csv-file ---
    try:
      with open(csv_grid, 'rb') as f:
	  reader = csv.reader(f, 'excel')
	  ii=0
	  for row in reader:
	      min_[ii]=float(row[1])
	      max_[ii]=float(row[3])
	      N_[ii]=int(row[5])
	      ii=ii+1
    except:
      #--- If the file cannot be loaded correctly, exit the program ---
      core.parser.error("Invalid .csv file!\n" +
	'Example:\n' +
	'\txmin=,-8.5,xmax=,8.5,Nx=,100' + '\n'
	'\tymin=,-8.5,ymax=,8.5,Ny=,100' + '\n'
	'\tzmin=,-8.5,zmax=,8.5,Nz=,100' + '\n')
  else: 
    #--- If the file does not exist, exit the program ---
    core.parser.error(csv_grid + ' does not exist!')  
  
  return 0
  #--- grid_csv ---

def center_grid(ac):
  #--- FUNCTION center_grid ----------------------------------------------------
  #--- Center the grid to the point ac and to the origin (0,0,0) ---------------
  #-----------------------------------------------------------------------------
  
  #--- All grid related variables should be globals ---
  global x, y, z, d3r, min_, max_, N_, delta_
  
  P=[numpy.zeros((3,1)), numpy.reshape(ac,(3,1))]
  
  d_tilde = numpy.abs(P[0] - P[1])
  N_tilde = numpy.round(numpy.abs(d_tilde / delta_))
  
  for ii in range(3): 
    if N_tilde[ii] != 0:
      delta_[ii] = d_tilde[ii] / N_tilde[ii]
  
  grid = [x, y, z]
  
  for ii in range(3):
    position = numpy.nonzero(ac[ii] <= grid[ii])[0][0]
    g = numpy.abs(grid[ii][position] - ac[ii]);
    c = 1/2.*delta_[ii] - g;
    grid[ii] += c;
  
  x = grid[0]  
  y = grid[1]  
  z = grid[2]
  d3r = numpy.product(delta_)
  
  min_ = [min(grid[0]), min(grid[1]), min(grid[2])]
  max_ = [max(grid[0]), max(grid[1]), max(grid[2])]
  N_   = [len(grid[0]), len(grid[1]), len(grid[2])]
  
  orbkit_output.display('Centered Grid to (%.2f %.2f %.2f): ' % (ac[0], ac[1], ac[2]))
  grid_display()
  
  for ii in range(3):
    if len(numpy.nonzero(0. == numpy.round(grid[ii]*10000))[0])!= 0: orbkit_output.display('Warning!\n\tAt least one grid point is equal to zero.\n')
  
  return 0
  #--- center_grid ---

def call_center_grid(geo_info,geo_spec,options):
  #--- FUNCTION call_center_grid -----------------------------------------------
  #--- Call the function center_grid centering the grid to one atom and (0,0,0)
  #-----------------------------------------------------------------------------
  
  #--- All grid related variables should be globals ---
  global x, y, z, d3r, min_, max_, N_, delta_
  
  if options.center: 
    try: 
      center_grid(geo_spec[options.a_num-1])
    except: 
      falsch = True;
      while falsch: ###FIXME
	orbkit_output.display('Failure: Valid atom numbers:')
	for ii in range(len(geo_info)): orbkit_output.display('\t' + str(geo_info[ii][0:2]))
	try:
	  options.a_num = int(raw_input('Atom number: '))
	  if(options.a_num > 0) and (options.a_num <= len(geo_info)):
	    falsch = False
	except: orbkit_output.display('Not a number!\n')
      
      center_grid(geo_spec[options.a_num-1])
  return 0
  #--- call_center_grid ---

#--- Default values for the grid parameters ---
min_ = [-8.0, -8.0, -8.0]
max_ = [ 8.0,  8.0,  8.0]
N_   = [ 101,  101,  101]

#--- Initialize some lists ---
x = [0]
y = [0]
z = [0]
delta_ = numpy.zeros((3,1))