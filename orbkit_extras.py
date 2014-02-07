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
import sys

import numpy
from scipy import integrate

# Import orbkit modules
import orbkit_core as core
import orbkit_grid as grid
import orbkit_output

def mo_select(geo_spec, geo_info, ao_spec, mo_spec, fid_mo_list, calc_mo=False):
  '''FUNCTION mo_select
  Calculate and save selected MOs (calc_mo = True)
  or Calculate and save the density only for selected MOs
  '''
  File_MO_List = []
  all_MO = []
  Selected_mo_spec = []
  Selected_mo_ii = [] 
  sym_select = False
  
  if core.options.calc_mo=='ALL_MO':
    selected_MO = numpy.array(numpy.arange(len(mo_spec))+1, dtype=numpy.str)
    File_MO_List = [selected_MO]
    Selected_mo_spec = mo_spec
  else:  
    try:
      fid=open(fid_mo_list,'rb')
      flines = fid.readlines()
      fid.close()
      for stuff in flines:
	stuff_int = stuff.split()
	File_MO_List.append(stuff_int)
	all_MO = all_MO+stuff_int
	
    except:
      orbkit_output.display('The selected MO-list (%(m)s) is not valid!' % 
		      {'m': fid_mo_list} + '\ne.g.\n\t1\t3\n\t2\t7\t9\n')
      sys.exit(1)
    
    #--- Check if the MOs are specified ---
    #--- by symmetry (e.g. 1.1 in MOLPRO nomenclature) or ---
    #--- by the number in the molden file (e.g. 1) ---
    selected_MO=list(set(all_MO))
    if '.' in selected_MO[0]: sym_select = True
    
    if sym_select:
      for k in range(len(mo_spec)):
	if mo_spec[k]['sym'] in selected_MO:
	  Selected_mo_spec.append(mo_spec[k])
	  Selected_mo_ii.append(mo_spec[k]['sym'])
      Selected_mo_ii = numpy.array(Selected_mo_ii)
    else:
      selected_MO = map(int, selected_MO)            
      selected_MO.sort()
      selected_MO = map(str, selected_MO)
      for k in range(len(mo_spec)):
	if str(k+1) in selected_MO:
	  Selected_mo_spec.append(mo_spec[k])
	  
  if calc_mo and core.options.hdf5:
    save_mo_hdf5(core.options.outputname,geo_info,geo_spec,ao_spec,
		  Selected_mo_spec)
    return 0
  elif calc_mo:
    #--- Calculate the AOs and MOs ---
    if options.save_ao:
      ao_count = core.ao_creator2(geo_spec,ao_spec,
				  fid=core.options.outputname)
      mo_list = core.mo_creator2(ao_count,Selected_mo_spec,
				  fid=core.options.outputname)
    elif options.discard_ao:
      mo_list = core.mo_creator3(geo_spec,ao_spec,Selected_mo_spec)
    else:
      ao_list = core.ao_creator(geo_spec,ao_spec,vector=False)
      mo_list = core.mo_creator(ao_list,Selected_mo_spec,vector=False)

    orbkit_output.display("Norm of the MOs:")
    MO=[]
    for ij in range(len(selected_MO)):
      if selected_MO[ij] in File_MO_List[ii]: 
	if sym_select: ii_MO = numpy.argwhere(Selected_mo_ii==selected_MO[ij])
	else: ii_MO = ij
	MO.append(mo_list[ii_MO])
	orbkit_output.display("\t%(m).6f\tMO %(n)s" 
	  % {'m':numpy.sum(mo_list[ii_MO]**2)*grid.d3r, 'n':mo_spec[ij]['sym']})
    #--- Create Output ---
    for ij in range(len(MO)):
      outputname = (core.options.outputname + '_MO_' + mo_spec[ij]['sym'])#selected_MO[ij])
      orbkit_output.main_output(MO[ij],geo_info,geo_spec,Selected_mo_spec,outputname)
    
    return MO
  else:
    for ii in range(len(File_MO_List)):
      orbkit_output.display('\nStarting with the ' + str(ii+1) + '. element of the MO-List (' 
		      + fid_mo_list + ')...\n\t' + str(File_MO_List[ii]) + 
		      '\n\t(Only regarding existing and occupied MOs.)\n')
      
      Spec = []
      for ij in range(len(selected_MO)):
	if selected_MO[ij] in File_MO_List[ii]: 
	  if sym_select: ii_MO = numpy.argwhere(Selected_mo_ii==selected_MO[ij])
	  else: ii_MO = ij
	  Spec.append(Selected_mo_spec[ii_MO])
      rho = core.rho_compute(geo_spec, geo_info, ao_spec, Spec)
      
      if core.options.reduced_density:
	rho = integrate.simps(rho, grid.x, dx=grid.delta_[0], axis=0, even='avg')
	rho = integrate.simps(rho, grid.y, dx=grid.delta_[1], axis=0, even='avg')
	
      outputname = (core.options.outputname + '_' + str(ii+1))
      orbkit_output.main_output(rho,geo_info,geo_spec,Spec,outputname)
    return 0
  #--- mo_select ---

def save_mo_hdf5(filename,geo_info,geo_spec,ao_spec,mo_spec,
		  x=None,y=None,z=None,N=None,only_a_few=False):
  import h5py
  if x == None: x = grid.x
  if y == None: y = grid.y
  if z == None: z = grid.z
  if N == None: N = grid.N_
  
  #--- Initialize HDF5_File ---
  fid = filename if filename.endswith('.h5') else '%s.h5' % filename
  f = h5py.File(fid, 'w')
  
  f.create_dataset('x',(1,len(x)),data=x)
  f.create_dataset('y',(1,len(y)),data=y)
  f.create_dataset('z',(1,len(z)),data=z)
  
  f.create_dataset('geo_info',(numpy.shape(geo_info)),data=numpy.array(geo_info))
  f.create_dataset('geo_spec',(numpy.shape(geo_spec)),data=geo_spec)
  
  mo_name = []
  occ_num=[]
  energy=[]
  sym=[]
  
  
  for mo in mo_spec:
    mo_name.append(mo['sym'])
    occ_num.append(mo['occ_num'])
    energy.append(mo['energy'])
    sym.append(mo['sym'])
    
    dID = 'MO:%s' % mo['sym']
    dset = f.create_dataset(dID,N)
  
  
  f.create_dataset('MO:Content',data=numpy.array(mo_name))
  
  MO_info = f.create_group('MO_info')
  MO_info.create_dataset('occ_num',((1,len(mo_spec))),data=occ_num)
  MO_info.create_dataset('energy',((1,len(mo_spec))),data=energy)
  MO_info.create_dataset('sym',((1,len(mo_spec))),data=sym)
  
  f.close()
  
  if only_a_few:
    mo_list = core.rho_compute(geo_spec,geo_info,ao_spec,mo_spec,calc_mo=True)
    
    orbkit_output.display('Writing %s ...' % fid)
    
    f = h5py.File(fid, 'a')
    
    for ii in range(len(mo_spec)):
      dID = 'MO:%s' % mo_spec[ii]['sym']
      f[dID][:,:,:] = mo_list[ii,:,:,:]
    
    f.close()
  
  else:
    #--- Slice the grid ---
    sDim = 0
    sNum = N[sDim]
    N[sDim] = 1
    xyz = [x,y,z]
    zz = [x,y,z]
    zz[sDim] = numpy.zeros(1)
    for ii_z in range(len(xyz[sDim])):
      zz[sDim][:] = xyz[sDim][ii_z]
      
      ao_list = core.ao_creator(geo_spec,ao_spec,x=zz[0],y=zz[1],z=zz[2],N=N)
      core.mo_creator(ao_list,mo_spec,x=zz[0],y=zz[1],z=zz[2],N=N,
		      vector=False,HDF5_save=fid,h5py=h5py,s=ii_z)
  
  return 0
  #--- mo_select ---

def atom2index(atom,geo_info=None):
  if not (isinstance(atom,list) or isinstance(atom,numpy.ndarray)):
    atom = [atom]
    
  index = []
  
  if geo_info != None:
    geo_info = numpy.array(geo_info)[:,1]
    for a in atom:
      i = numpy.argwhere(geo_info.astype(int) == a)
      if len(i) != 1:
	orbkit_output.display('Error: No or multiple occurence ' + 
			      'of the atom number %d in geo_info!' % a)
	return 1
      index.append(int(i))
    
    index = numpy.array(index, dtype=int)
  
  else:
    try:
      index = numpy.array(atom,dtype=int)
    except ValueError:
      orbkit_output.display('Error: Cannot convert atom to integer array!')
      return 1
  
  return atom, index

def atom_projected_density(atom,geo_spec,ao_spec,mo_spec,geo_info=None,
			    bReturnMO=False,ao_list=None,mo_list=None,
			    x=None,y=None,z=None,N=None):
  '''Compute the projected electron density with respect of the selected atoms.
  
  Rho^a = \sum_i occ_i * mo_i^a * mo_i
  
  NOT TESTED YET!
  
  Parameters
  ----------
  atom: int or list of int
	    
  all_MO:   bool, optional
	    If True, all molecular orbitals are returned.
  
  Returns
  -------
  geo_spec, geo_info, ao_spec, mo_spec	  
	    See the manual for details.
	    
  '''  
  orbkit_output.display('Attention: ' + 
	'The function atom_projected_density is not fully tested yet!\n')
  
  if x == None: x = grid.x
  if y == None: y = grid.y
  if z == None: z = grid.z
  if N == None: N = numpy.array(grid.N_)
  
  atom, index = atom2index(atom,geo_info=geo_info)
  
  orbkit_output.display('Computing the atom-projected density with respect to '+ 
			'the atom(s)')
  outp = '\t['
  for i,a in enumerate(atom):
    if not i % 10 and i != 0:
      outp += '\n\t'
    outp += '%d,\t' % a
  orbkit_output.display('%s]\n' % outp[:-2])
  
  orbkit_output.display('\tCalculating ao_list & mo_list')
  if ao_list == None:
    ao_list = core.ao_creator(geo_spec,ao_spec,x=x,y=y,z=z,N=N)
  if mo_list == None:
    mo_list = core.mo_creator(ao_list,mo_spec,x=x,y=y,z=z,N=N)
    
  orbkit_output.display('\tCalculating the atom-projected density')
  
  if bReturnMO: mo_atom = [[] for a in index]
  rho_atom = [numpy.zeros(N) for a in index]
  for i,a in enumerate(index):
    orbkit_output.display('\t\tFinished %d of %d' % (i+1,len(index)))
    ao_index = []
    ao = []
    ll = 0
    for ii in ao_spec:
      for l in range(core.l_deg(l=ii['type'])):
	if ii['atom'] == a:
	  ao_index.append(ll)
	  ao.append(ao_list[ll])
	ll += 1
    
    for ii_mo,spec in enumerate(mo_spec):
      mo = numpy.zeros(N)
      for jj in range(len(ao_index)):
	mo += spec['coeffs'][ao_index[jj]] * ao[jj]
      rho_atom[i] += spec['occ_num'] * mo_list[ii_mo] * mo
      if bReturnMO: mo_atom[i].append(mo)
  
  if bReturnMO:
    orbkit_output.display('Returning the atom-projected density and' +
    '\n\tthe atom-projected molecular orbitals')
    return rho_atom, mo_atom
  else:
    orbkit_output.display('Returning the atom-projected density')
    return rho_atom
  #--- atom_projected_density ---

def compute_mulliken_charges(atom,geo_info,geo_spec,ao_spec,mo_spec,
			    ao_list=None,mo_list=None,
			    x=None,y=None,z=None,N=None):
  '''Compute the Mulliken charges of the selected atoms using
  the respective projected electron densities.
  
  Rho^a = \sum_i occ_i * mo_i^a * mo_i
  
  NOT TESTED YET!
  
  Parameters
  ----------
  atom: int or list of int
	    
  all_MO:   bool, optional
	    If True, all molecular orbitals are returned.
  
  Returns
  -------
  geo_spec, geo_info, ao_spec, mo_spec	  
	    See the manual for details.
	    
  '''
  orbkit_output.display('Attention: ' + 
	'The function compute_mulliken_charges is not fully tested yet!\n')
  
  if x == None: x = grid.x
  if y == None: y = grid.y
  if z == None: z = grid.z
  if N == None: N = numpy.array(grid.N_)
  
  atom, index = atom2index(atom,geo_info=geo_info)
  
  rho_atom = atom_projected_density(index,geo_spec,ao_spec,mo_spec,
			    ao_list=ao_list,mo_list=mo_list,
			    x=x,y=y,z=z,N=N)
  
  orbkit_output.display('Mulliken Charges:')
  mulliken_charge = []  
  for i,a in enumerate(index):
    mulliken_charge.append(core.integration(rho_atom[i]))
    #--- Print the Mulliken Charges ---
    orbkit_output.display('\tAtom %s (%s):\t%+0.4f' % 
		  (geo_info[a][1],geo_info[a][0],mulliken_charge[-1]))
  
  return mulliken_charge


def mo_transition_flux_density(i,j,geo_spec,ao_spec,mo_spec,drv='x',
			    ao_list=None,mo_list=None,
			    delta_ao_list=None,delta_mo_list=None,
			    x=None,y=None,z=None,N=None):
  '''Calculate one component (e.g. drv='x') of the 
  transition electoronic flux density between the 
  molecular orbitals i and j.
  
  MOTEFD_{i,j}(r) = <mo_i|\delta(r-r')\nabla_x|mo_j>
  
  NOT TESTED YET!
  
  Parameters
  ----------
  i: int
    index of the first MO (Bra)
  j: int
    index of the second MO (Ket)
  drv: {0,1,2,'x','y','z'}
    The desired component of the vector field of the 
    transition electoronic flux density
	    
  '''
  orbkit_output.display('Attention: ' + 
	'The function mo_transition_flux_density is not fully tested yet!\n')
  
  if x == None: x = grid.x
  if y == None: y = grid.y
  if z == None: z = grid.z
  if N == None: N = numpy.array(grid.N_)
  
  
  if mo_list == None:
    if ao_list == None:
      orbkit_output.display('\tComputing ao_list and ' +
			    'MO #%d since it is not given.' % i)
      ao_list = core.ao_creator(geo_spec,ao_spec,x=x,y=y,z=z,N=N)  
    else:
      orbkit_output.display('\tComputing MO #%d since it is not given.' % i)
    mo = core.mo_creator(ao_list,[mo_spec[i]])[0]
  else:
    mo = mo_list[i]
  if delta_mo_list == None:
    if delta_ao_list == None:
      orbkit_output.display('\tComputing delta_ao_list and the derivative of ' +
			    'MO #%d since it is not given.' % j)
      delta_ao_list = core.delta_ao_creator(geo_spec,ao_spec,drv=drv,
					  x=x,y=y,z=z,N=N)
    else:
      orbkit_output.display('\tComputing MO #%d since it is not given.' % j)
    delta_mo = core.mo_creator(delta_ao_list,[mo_spec[j]])[0]
  else:
    delta_mo = delta_mo_list[i]
  
  return mo*delta_mo
  #--- mo_transition_flux_density ---


def colormap_creator(rho,filename,n_peaks=5,start=0.01,stop=0.999,peak_width=0.1):
  #--- FUNCTION colormap_creator -----------------------------------------------
  #--- Create a .cmap colormap for ZIBAmira adjusted to the density ------------
  #--- Default: Isosurface values between 1% and 99.9% of the total density ----
  #-----------------------------------------------------------------------------
  
  #--- Sort the density values ---
  a=numpy.reshape(rho,-1)
  a=a[numpy.argsort(a)]
  
  #--- Where do we have the start and stop percentage? ---
  n1=int(start*len(a))
  n2=int(stop*len(a))
  rho_min=a[n1]
  rho_max=a[n2]
  
  #--- Compute the distance between two isosurface values ---
  delta_peak =(rho_max-rho_min)/(n_peaks-1)
  
  #--- Open a cmap file ---
  fid = open('%(f)s.cmap' % {'f': filename}, 'w')
  
  #--- Write the header ---
  fid.write('<!DOCTYPE Colormap>\n')
  fid.write('<ColormapVisage2.0 Name="%(f)s">\n' % {'f': filename})
  fid.write('  <Graph Active="1" Type="0" Name="">\n')
  
  #--- Initialize a counter for the contour values ---
  counter = 0
  
  #--- Initialize a string for the contours ---
  c_str = ('    <Control Opacity="%(o)f" Number="%(c)d" Blue="%(v)f"' + 
	      ' Red="%(v)f" Green="%(v)f" Value="%(p)f"/>\n' )
  
  #--- Write the initial value at zero with zero opacity ---
  fid.write(c_str % {'o': 0, 'c': counter, 'v': 0, 'p': 0})
  counter += 1
  
  #--- Loop over the contour values ---
  for ii in range(n_peaks):
    #--- Calculate the value for the isosurface peak and two values ---
    #--- next to the peak with zero opacity ---
    peak = rho_min+delta_peak*(ii+1)
    peak_minus = peak * (1 - peak_width/2.)
    peak_plus = peak * (1 + peak_width/2.)
    
    #--- Calculate a value for the opacity and the color---
    value = 1-(float(ii+1)/(n_peaks+1)*0.9)
    
    #--- Write the peak ---
    #--- Left edge of the peak ---
    fid.write(c_str % {'o': 0, 'c': counter, 'v': value, 'p': peak_minus})
    counter += 1
    #--- The peak ---
    fid.write(c_str % {'o': 1-value, 'c': counter, 'v': value, 'p': peak})
    counter += 1
    #--- Right edge of the peak ---
    fid.write(c_str % {'o': 0, 'c': counter, 'v': value, 'p': peak_plus})
    counter += 1
  
  #--- Write a final value at 1.5 * (final peak value) with zero opacity ---
  fid.write(c_str % {'o': 0, 'c': counter, 'v': 0, 'p': peak_plus*1.5})
  
  #--- Finalize the file ---
  fid.write('  </Graph>\n')
  fid.write('</ColormapVisage2.0>')
  
  #--- Close the file ---
  fid.close()
  
  return 0
  #--- colormap_creator ---

def hx_network_creator(rho,filename):
  #--- FUNCTION hx_network_creator ---------------------------------------------
  #--- Create a ZIBAmira hx network including a .cmap colormap file ------------
  #--- adjusted to the density for the easy depiction of the density -----------
  #-----------------------------------------------------------------------------
  
  #--- Create a .cmap colormap file using the default values ---
  orbkit_output.display('\tCreating ZIBAmira colormap file...\n\t\t%(f)s.cmap' % 
		      {'f': filename})
  colormap_creator(rho,filename)
  
  #--- Create a .hx network file based on the file orbkit_hx_network_draft.hx ---
  orbkit_output.display('\tCreating ZIBAmira network file...\n\t\t%(f)s.hx' % 
		      {'f': filename})
  #--- Open an empty file
  fid = open('%(f)s.hx' % {'f': filename},'w')
  
  filename = filename.split('/')[-1]
  #--- Copy the content of the draft file and replace the keywords ---
  for line in open('orbkit_hx_network_draft.hx'):
    line = line.replace("FILENAME",filename)
    fid.write(line) 
  
  #--- Close the file ---
  fid.close()  
  
  return 0

#--- Class for the creation of a MO amira network ---
#### NOT FINISHED YET ###
class cMO_display:
  def __init__(self):
    self.filename = core.options.outputname + '_MO'
    orbkit_output.display("\tCreating ZIBAmira network file for \n\t  the depiction of the MOs...\n\t\t" + self.filename + ".hx")    
    fid = open(self.filename + '.hx','w') 
    name = self.filename.split('/')[-1] 
    fid.close() 
    
  def colormap(self):
    orbkit_output.display("\tCreating ZIBAmira colormap file for \n\t  the depiction of the MOs...\n\t\t" + self.filename + ".hx")    
    fid = open(self.filename + '.col.am','w') 
    col_am = '''
# AmiraMesh 3D ASCII 2.0\n\n
define Lattice 4\n
Parameters {
    ContentType "Colormap",
    MinMax -50 50,
    Interpolate 0
}\n
Lattice { float[4] Data } @1\n
# Data section follows\n@1
1.000000000000000e+00 0.000000000000000e+00 0.000000000000000e+00 1.000000000000000e+00 
1.000000000000000e+00 0.000000000000000e+00 0.000000000000000e+00 1.000000000000000e+00 
0.000000000000000e+00 0.000000000000000e+00 1.000000000000000e+00 1.000000000000000e+00 
0.000000000000000e+00 0.000000000000000e+00 1.000000000000000e+00 1.000000000000000e+00 
    '''
    fid.write(col_am)
    fid.close()
#def hx_network_MO_diplay_creator
  ## Read line by line
  #counter = 0
  #factor=abs(float(flines[-1].split(" ")[0])-float(flines[0].split(" ")[0]))/100;
  #for ii, line in enumerate(flines):
    #counter += 1;
    #line = line.replace('\n', '')
    #column = line.split(" ")
    #if column[1] != str(0):
	#fid.write('    <Control Opacity="'+ str(0) + '" Number="'+ str(counter) + '" Blue="0" Red="0" Value="' + str(float(column[0])-factor) +'" Green="0"/>\n')
	#counter += 1;
	#fid.write('    <Control Opacity="'+ str(column[1]) + '" Number="'+ str(counter) + '" Blue="0" Red="0" Value="' + str(column[0]) +'" Green="0"/>\n')
	#counter += 1;
	#fid.write('    <Control Opacity="'+ str(0) + '" Number="'+ str(counter) + '" Blue="0" Red="0" Value="' + str(float(column[0])+factor) +'" Green="0"/>\n')
	#counter += 1;
    #else:
      #fid.write('    <Control Opacity="'+ str(column[1]) + '" Number="'+ str(counter) + '" Blue="0" Red="0" Value="' + str(column[0]) +'" Green="0"/>\n')

  #fid.write('  </Graph>\n')
  #fid.write('</ColormapVisage2.0>')

  #fid.close()
