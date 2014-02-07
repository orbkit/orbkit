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
import os
import numpy

# Import orbkit modules
import orbkit_grid as grid
import orbkit_core as core
import orbkit_extras as extras 

def display(string):
  try:
    if not core.options.quiet:
      print(string)
    if not core.options.no_log:
      fid = '%s.log' % core.options.outputname
      f = open(fid, 'a')
      f.write('%s\n' % string)
      f.close()
  except AttributeError:
    print(string)

def rm_log(no_log=False):
  if not no_log:
    try:
	os.remove('%s.log' % core.options.outputname)
    except OSError:
	pass
  

def main_output(rho,geo_info,geo_spec,mo_spec,outputname):
  #--- Function to create the desired output -----------------------------------
  if core.options.amira:
    display('Saving to ZIBAmiraMesh file...' +
		    '\n\t%(o)s.am' % {'o': outputname})
    amira_creator(rho,outputname)
  elif core.options.hdf5:
    display('Saving to Hierarchical Data Format file (HDF5)...' +
		    '\n\t%(o)s.h5' % {'o': outputname})
    HDF5_creator(rho,outputname,geo_info,geo_spec,mo_spec)
  else:
    display('Saving to .cb file...' +
		    '\n\t%(o)s.cb' % {'o': outputname})
    cube_creator(rho,outputname,geo_info,geo_spec)
    #output_creator(rho,outputname,geo_info,geo_spec)	# Axel's cube files
  
  if core.options.hx_network:
    display('\nCreating ZIBAmira network file...')
    if core.options.amira == False:    
      display('\tSaving to ZIBAmiraMesh file...' +
		      '\n\t\t%(o)s.am' % {'o': outputname})
      amira_creator(rho,core.options.outputname)
    #--- Create Amira network incl. Alphamap ---
    extras.hx_network_creator(rho,outputname)
  return 0
  #--- main_output ---

def output_creator(rho,filename,geo_info,geo_spec):
  #--- FUNCTION output_creator -------------------------------------------------
  #--- Creates a standard plain text output ------------------------------------
  #-----------------------------------------------------------------------------
  
  #--- Open an empty file ---
  fid = open('%(f)s.cb' % {'f': filename}, 'w')
  
  #--- Write the type and the position of the atoms in the header ---
  fid.write('Position of the nuclei (in a_0)\n')
  string = 'Number of nuclei:\t%d\n' % len(geo_info)
  for ii in range(len(geo_info)):
    for geo in geo_info[ii]:
      string += '%s\t' % geo
    for geo in geo_spec[ii]:
      string += '%0.8f\t' % geo
    string += '\n' 
  
  fid.write(string)
  
  # Write grid information for grid reconstruction
  fid.write('Grid data (in a_0)\n')
  fid.write(grid.grid_display(quiet=True,start=''))
  
  # Write density information
  fid.write('Data section begins (z runs fastest)  \n')
  
  string = ''
  for rr in range(len(grid.x)):
    for ss in range(len(grid.y)):
      for tt in range(len(grid.z)):
        string += '%+g\n' % rho[rr,ss,tt]
      string += '\n'
  
  fid.write(string)
  
  #--- Close the file ---
  fid.close()
  
  return 0
  #--- output_creator ---

def cube_creator(rho,filename,geo_info,geo_spec):
  #--- FUNCTION output_creator -------------------------------------------------
  #--- Creates a standard plain text output ------------------------------------
  #-----------------------------------------------------------------------------
  
  #--- Open an empty file ---
  fid = open('%(f)s.cb' % {'f': filename}, 'w')
  
  #--- Write the type and the position of the atoms in the header ---
  string = ' orbkit Density calculation\n'
  string += ' %(f)s\n'  % {'f': filename}
  #--- How many atoms ---
  string += ('%(at)d' % {'at': len(geo_info)}).rjust(6)
  #--- Minima
  for ii in range(3):
    string += ('%(min)0.6f' % {'min': grid.min_[ii]}).rjust(12)

  for ii in range(3):
    string += '\n'
    string += ('%(N)d' % {'N': grid.N_[ii]}).rjust(6)
    for jj in range(3):
      if jj == ii:
	string += ('%(dr)0.6f' % {'dr': grid.delta_[ii]}).rjust(12)
      else:
	string += ('%(dr)0.6f' % {'dr': 0}).rjust(12)
  
  for ii in range(len(geo_info)):
    string += '\n'
    string += ('%(N)s' % {'N': geo_info[ii][1]}).rjust(6)
    string += ('%(ch)0.6f' % {'ch': float(geo_info[ii][2])}).rjust(12)
    for jj in range(3):
      string += ('%(r)0.6f' % {'r': geo_spec[ii][jj]}).rjust(12)
  string += '\n'
  for rr in range(len(grid.x)):
    for ss in range(len(grid.y)):
      for tt in range(len(grid.z)):
	string += ('%(rho).6E' % {'rho': rho[rr,ss,tt]}).rjust(13)
	if (tt % 6 == 5):
	  string += '\n'
      string += '\n'
  
  
  fid.write(string)
  
  #--- Close the file ---
  fid.close()
  
  return 0
  #--- output_creator ---

def amira_creator(rho,filename):
  #--- FUNCTION amira_creator --------------------------------------------------
  #--- Creates a ZIBAmira mesh file --------------------------------------------
  #-----------------------------------------------------------------------------
  
  #--- Open an empty file ---
  fid = open('%(f)s.am' % {'f': filename},'w')

  # usage:
  #     - open Amira
  #     - left-click File -> Open Data
  #     - choose the sys.argv[2].am
  #     - Press OK
  #     - right-click on the new Data Object -> Compute -> Arithmetic
  #     - choose the Arithmetic object
  #     - select  Expr -> A and Result type -> regular and Apply
  #     - use the new data object to display your density as usual,
  
  #--- Write Header ---
  fid.write('# AmiraMesh 3D ASCII 2.0\n\n\n')
  fid.write('define Lattice %(Nx)d %(Ny)d %(Nz)d\n' % 
	      {'Nx': grid.N_[0],'Ny': grid.N_[1],'Nz': grid.N_[2]})
  fid.write('define Coordinates %(N)d\n\n' % {'N': numpy.sum(grid.N_)})
  fid.write('Parameters {\n')
  fid.write('    Content "%(Nx)dx%(Ny)dx%(Nz)d float, uniform coordinates",\n' %
	      {'Nx': grid.N_[0],'Ny': grid.N_[1],'Nz': grid.N_[2]})
  fid.write('    BoundingBox %(xmin)f %(xmax)f %(ymin)f %(ymax)f %(zmin)f %(zmax)f,\n' %
	      {'xmin': grid.min_[0],'xmax': grid.max_[0],
	       'ymin': grid.min_[1],'ymax': grid.max_[1],
	       'zmin': grid.min_[2],'zmax': grid.max_[2]})
  fid.write('    CoordType "uniform"\n}\n\n')
  fid.write('Lattice { float Data } @1\n')
  fid.write('# Data section follows\n@1\n')
  
  # Write density information to .am file
  string = ''
  for tt in range(len(grid.z)):
    for ss in range(len(grid.y)):
  	  for rr in range(len(grid.x)): 
	    string += '%g\n' % rho[rr,ss,tt]
  
  fid.write(string)
  
  #--- Close the file ---
  fid.close()
  
  return 0
  #--- amira_creator ---

def amira_creator_vector(rho,filename):
  #--- FUNCTION amira_creator --------------------------------------------------
  #--- Creates a ZIBAmira mesh file using a recti --------------------------------------------
  #-----------------------------------------------------------------------------
  
  #--- Open an empty file ---
  fid = open('%(f)s.am' % {'f': filename},'w')

  # usage:
  #     - open Amira
  #     - left-click File -> Open Data
  #     - choose the sys.argv[2].am
  #     - Press OK
  #     - right-click on the new Data Object -> Compute -> Arithmetic
  #     - choose the Arithmetic object
  #     - select  Expr -> A and Result type -> regular and Apply
  #     - use the new data object to display your density as usual,
  
  N = len(rho)
  
  #--- Write Header ---
  fid.write('# AmiraMesh 3D ASCII 2.0\n\n\n')
  fid.write('define Lattice %(Nx)d %(Ny)d %(Nz)d\n' % 
	      {'Nx': N,'Ny': N,'Nz': N})
  fid.write('define Coordinates %(N)d\n\n' % {'N': N})
  fid.write('Parameters {\n')
  fid.write('    Content "%(Nx)dx%(Ny)dx%(Nz)d float, rectilinear coordinates",\n' %
	      {'Nx': N,'Ny': N,'Nz': N})
  fid.write('    CoordType "rectilinear"\n}\n\n')
  fid.write('Lattice { float Data } @1\n')
  fid.write('Coordinates { float xyz } @2\n\n')
  fid.write('# Data section follows\n@1\n')
  
  # Write density information to .am file
  string = ''
  #for tt in range(len(grid.z)):
    #for ss in range(len(grid.y)):
  	  #for rr in range(len(grid.x)): 
	    #string += '%0.8f\n' % rho[rr,ss,tt]
  for ii in range(len(rho)):
    string += '%0.8f\n' % rho[ii]
  string += '\n@2\n'
  for xx in grid.x: 
    string += '%0.8f\n' % xx
  for yy in grid.y:  
    string += '%0.8f\n' % yy
  for zz in grid.z:  
    string += '%0.8f\n' % zz
  
  fid.write(string)
  
  #--- Close the file ---
  fid.close()
  
  return 0
  #--- amira_creator_rectilinear_coord ---

def amira_creator_rectilinear_coord(rho,filename):
  #--- FUNCTION amira_creator --------------------------------------------------
  #--- Creates a ZIBAmira mesh file using a recti --------------------------------------------
  #-----------------------------------------------------------------------------
  
  #--- Open an empty file ---
  fid = open('%(f)s.am' % {'f': filename},'w')

  # usage:
  #     - open Amira
  #     - left-click File -> Open Data
  #     - choose the sys.argv[2].am
  #     - Press OK
  #     - right-click on the new Data Object -> Compute -> Arithmetic
  #     - choose the Arithmetic object
  #     - select  Expr -> A and Result type -> regular and Apply
  #     - use the new data object to display your density as usual,
  
  #--- Write Header ---
  fid.write('# AmiraMesh 3D ASCII 2.0\n\n\n')
  fid.write('define Lattice %(Nx)d %(Ny)d %(Nz)d\n' % 
	      {'Nx': grid.N_[0],'Ny': grid.N_[1],'Nz': grid.N_[2]})
  fid.write('define Coordinates %(N)d\n\n' % {'N': numpy.sum(grid.N_)})
  fid.write('Parameters {\n')
  fid.write('    Content "%(Nx)dx%(Ny)dx%(Nz)d float, rectilinear coordinates",\n' %
	      {'Nx': grid.N_[0],'Ny': grid.N_[1],'Nz': grid.N_[2]})
  fid.write('    CoordType "rectilinear"\n}\n\n')
  fid.write('Lattice { float Data } @1\n')
  fid.write('Coordinates { float xyz } @2\n\n')
  fid.write('# Data section follows\n@1\n')
  
  # Write density information to .am file
  string = ''
  for tt in range(len(grid.z)):
    for ss in range(len(grid.y)):
  	  for rr in range(len(grid.x)): 
	    string += '%0.8f\n' % rho[rr,ss,tt]
  string += '\n@2\n'
  for xx in grid.x: 
    string += '%0.8f\n' % xx
  for yy in grid.y:  
    string += '%0.8f\n' % yy
  for zz in grid.z:  
    string += '%0.8f\n' % zz
  
  fid.write(string)
  
  #--- Close the file ---
  fid.close()
  
  return 0
  #--- amira_creator_rectilinear_coord ---

def HDF5_creator(rho,filename,geo_info,geo_spec,mo_spec):
  # --- FUNCTION HDF5_creator --- creates HDF5 file (Hierarchical Data Format) ------------------------
  import h5py
  f = h5py.File(filename + '.h5', 'w')
  if core.options.reduced_density:
    dset = f.create_dataset('z',(1,len(grid.z)),data=grid.z)
    dset = f.create_dataset('rho',((1,len(grid.z))),data=rho)
  else:
    dset = f.create_dataset('x',(1,len(grid.x)),data=grid.x)
    dset = f.create_dataset('y',(1,len(grid.y)),data=grid.y)
    dset = f.create_dataset('z',(1,len(grid.z)),data=grid.z)
    dset = f.create_dataset('rho',((len(grid.x),len(grid.y),len(grid.z))),data=rho)
  MO_info = f.create_group('MO_info')
  occ_num=[]
  energy=[]
  sym=[]
  for ii in range(len(mo_spec)):
    occ_num.append(mo_spec[ii]['occ_num'])
    energy.append(mo_spec[ii]['energy'])
    sym.append(mo_spec[ii]['sym'])
  dset = MO_info.create_dataset('occ_num',((1,len(mo_spec))),data=occ_num)
  dset = MO_info.create_dataset('energy',((1,len(mo_spec))),data=energy)
  dset = MO_info.create_dataset('sym',((1,len(mo_spec))),data=sym)
  
  dset = f.create_dataset('geo_info',(numpy.shape(geo_info)),data=numpy.array(geo_info))
  dset = f.create_dataset('geo_spec',(numpy.shape(geo_spec)),data=geo_spec)
  f.close()
  
  return 0
  #--- HDF5_creator ---

def HDF5_creator_geo(filename,geo_info,geo_spec):
  import h5py
  fid = filename if filename.endswith('.h5') else '%s.h5' % filename
  f = h5py.File(fid, 'a')
  dset = f.create_dataset('geo_info',(numpy.shape(geo_info)),data=numpy.array(geo_info))
  dset = f.create_dataset('geo_spec',(numpy.shape(geo_spec)),data=geo_spec)
  f.close()
  
  return 0

def HDF5_creator_MO(MO,filename,geo_info,geo_spec,mo_spec):
  # --- FUNCTION HDF5_creator --- creates HDF5 file (Hierarchical Data Format) ------------------------
  import h5py
  f = h5py.File(filename + '.h5', 'w')
  
  dset = f.create_dataset('x',(1,len(grid.x)),data=grid.x)
  dset = f.create_dataset('y',(1,len(grid.y)),data=grid.y)
  dset = f.create_dataset('z',(1,len(grid.z)),data=grid.z)
  
  mo_name = []
  for ij in range(len(MO)):
    mo_name.append(mo_spec[ij]['sym'])
    dID = 'MO:%s' % mo_spec[ij]['sym']
    
    dset = f.create_dataset(dID,((len(grid.x),len(grid.y),len(grid.z))),data=MO[ij])
  
  dset = f.create_dataset('MO:Content',data=numpy.array(mo_name))
  
  MO_info = f.create_group('MO_info')
  occ_num=[]
  energy=[]
  sym=[]
  for ii in range(len(mo_spec)):
    occ_num.append(mo_spec[ii]['occ_num'])
    energy.append(mo_spec[ii]['energy'])
    sym.append(mo_spec[ii]['sym'])
  dset = MO_info.create_dataset('occ_num',((1,len(mo_spec))),data=occ_num)
  dset = MO_info.create_dataset('energy',((1,len(mo_spec))),data=energy)
  dset = MO_info.create_dataset('sym',((1,len(mo_spec))),data=sym)
  
  dset = f.create_dataset('geo_info',(numpy.shape(geo_info)),data=numpy.array(geo_info))
  dset = f.create_dataset('geo_spec',(numpy.shape(geo_spec)),data=geo_spec)
  f.close()
  
  return 0
  #--- HDF5_creator ---

def h5am_creator(rho,filename):
  # --- FUNCTION HDF5_creator --- creates HDF5 file (Hierarchical Data Format) ------------------------
  import h5py

  f = h5py.File(filename + '.h5am', 'w')

  amira = f.create_group('amira')
  amira.attrs['version'] 	= 2
  amira.attrs['contenttype'] 	= 'HDF5amiralattice'
  amira.attrs['numdims'] 	= 3
  amira.attrs['dims'] 		= [len(grid.x),len(grid.y),len(grid.z)]  
  amira.attrs['boundingbox'] 	= [grid.x[0],grid.x[-1],grid.y[0],grid.y[-1],grid.z[0],grid.z[-1]]
  amira.attrs['latticetype'] 	= 'http://amira.zib.de/latticetypes#uniform'
  #amira.attrs['latticetype'] 	= 'http://amira.zib.de/latticetypes#rectilinear'
  #amira.attrs['coordX'] 	= grid.x
  #amira.attrs['coordY'] 	= grid.y  
  #amira.attrs['coordZ'] 	= grid.z
  amira.attrs['ndatasets'] 	= 1

  dataset0 = amira.create_group('dataset:0')
  dataset0.attrs['ndatavar'] 	= 0
  dataset0.attrs['datatype'] 	= 'http://amira.zib.de/types#double'
  a=h5py.h5t.array_create(h5py.h5t.IEEE_F64LE, (len(grid.x),len(grid.y),len(grid.z)))
  dset = dataset0.create_dataset('timestep:0',(len(grid.x),len(grid.y),len(grid.z)),data=rho)
  
  #print dset
  #dset = dataset0.create_dataset('timestep:0',(len(grid.x),len(grid.y),len(grid.z)),data=rho)
    
  #with h5py.File(filename + '.h5am', 'w') as f:
    #f['/amira/Attribute/contenttype'] = 'HDF5amiralattice';
    #f['/amira/Attribute/version'] = 2
    #f['/amira/Attribute/numdims'] = 3
    #f['/amira/Attribute/dims'] = [len(grid.x),len(grid.y),len(grid.z)]
    #f['/amira/Attribute/boundingbox'] = [grid.x[0],grid.x[-1],grid.y[0],grid.y[-1],grid.z[0],grid.z[-1]]
    #f['/amira/Attribute/latticetype'] = 'http://amira.zib.de/latticetypes#uniform'
    #f['/amira/Attribute/ndatasets'] = 1
    #f['/amira/Attribute/ndatavar'] = 1
    #f['/amira/Attribute/datatype'] = 'http://amira.zib.de/types#float'
    #f['/amira/dataset:0/timestep:0'] =rho

