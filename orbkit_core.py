# -*- coding: iso-8859-1 -*-

lgpl = '''
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
import copy
import os
import optparse
import re
import string
import sys
import time

import numpy
from scipy import weave
from scipy import integrate

import multiprocessing as mp
from multiprocessing import Process
from multiprocessing import Pool

# Import orbkit modules
import orbkit_grid as grid
import orbkit_output 

def init_parser():
  #--- FUNCTION init_parser ----------------------------------------------------
  #--- Initialize parser and process the options -------------------------------
  #-----------------------------------------------------------------------------
  
  global parser, options, args
  default_fid='';default_oid='';
  
  usage = 'Usage: %prog [options] <input> <outputname>'
  parser = optparse.OptionParser(usage) 
  parser.add_option("-l", dest="show_lgpl",
                      default=False,action="store_true", 
                      help="show license information and exit")
  parser.add_option("-i", "--input", dest="filename",
                      default=default_fid, type="string",
                      help="molden input file")
  parser.add_option("-o", "--output",dest="outputname",
                      default=default_oid, type="string",
                      help="name of the output file")
  #parser.add_option("-g", "--gaussian",dest="gaussian",
                      #action="store_true", 
                      #help="read gaussian output file (*.log, *.FChk)") 
  parser.add_option("-g", "--gamess",dest="gamess",
                      action="store_true", 
                      help="read gamess output file (*.log)")  
  parser.add_option("--quiet",dest="quiet",
                      default=False,action="store_true", 
                      help="suppress terminal output")  
  parser.add_option("--no_log",dest="no_log",
                      default=False,action="store_true", 
                      help="suppress output of a <outputname>.log file")  
  parser.add_option("-p", "--numProc",dest="numproc",
                      default=1, type="int",
                      help="number of subprocesses started during the execution (default=1)")
  parser.add_option("--hdf5",dest="hdf5",
                      action="store_true", 
                      help="store Hierarchical Data Format file (.h5)")  
  parser.add_option("--amira",dest="amira",
                      default=False,action="store_true", 
                      help="store AmiraMesh file (.am)")   
  parser.add_option("--hx_network",dest="hx_network",
                      action="store_true", 
                      help="create ZIBAmira network (.hx)") 
  parser.add_option("--mo_list",dest="mo_list",
                      default=False, type="string", 
                      help="Read a file containing row vectors with the index of the selcted MOs (delimiter=' ')")  
  parser.add_option("--calc_mo",dest="calc_mo",
                      default=False, type="string", 
                      help="Calculate and save the MOs specified by the indices in their selected file (delimiter=' ')"+
                      ". Type 'ALL_MO' to store all the orbitals")
  parser.add_option("--reduced_density",dest="reduced_density",
                      default=False, action="store_true", 
                      help="Reduce the density with respect to the z-axis")             
  #parser.add_option("--save_ao",dest="save_ao",
                      #action="store_true", 
                      #help="write and read AOs (ao_<output_file>_0.p) ")     
  parser.add_option("--discard_ao",dest="discard_ao",
                      action="store_true", 
                      help="compute and discard AOs during the MO calculation") 
  parser.add_option("--all_MO",dest="all_MO",
                      default=False, action="store_true", 
                      help="Calculate all (occupied and virtual) MOs")
  parser.add_option("-c",dest="center",
                      action="store_true", help="center grid...") 
  parser.add_option("--atom",dest="a_num",
                      type="int", help="\tto atom number and the origin")  
  parser.add_option("--read_grid",dest="csv_grid",
                      default=False, type="string", 
                      help="Read a .csv file containing information about the grid (delimiter=',')")           
  parser.add_option("--vector",dest="vector_grid",
                      default=False,action="store_true",
                      help="Under construction...")
  (options, args) = parser.parse_args()                                                  
  
  #--- The following parser options are temporarily hidden: 
  options.save_ao 	= False
  options.gaussian 	= None
  options.no_slice 	= False
  options.no_output	= False
  options.delta_slice	= 1e4
  
  if options.show_lgpl:
    print(lgpl.replace('\nThis file is part of orbkit.\n',''))
    sys.exit(0)
  
  if options.save_ao and options.discard_ao:
    parser.error('\noptions --save_ao and --discard_ao are mutually exclusive\n')
  if options.amira and options.hdf5:
    parser.error('\nYou can only create one output file at one time.\nThus, the options --amira and --hdf5 are mutually exclusive.\n')
  
  while (os.path.exists(options.filename))==False:
    options.filename = raw_input('Insert a valid molden input file: ')

  if len(args) == 2:
    options.filename = args[0]
    options.outputname = args[1]
  else:
    if options.filename=='':
      options.filename = raw_input('Insert molden input file: ')

    if options.outputname=='':
      options.outputname = raw_input('Insert name of the output file: ')

  # Check if h5py is installed
  if options.hdf5 or options.save_ao:
    global h5py
    try: 
      import h5py
    except ImportError:
      parser.error('\nImportError: The module h5py is not installed!\n')
  
  if (options.mo_list or options.calc_mo) != False: options.all_MO = True
  return 0
  #--- init_parser ---
  
def read_molden(filename, all_MO=False):
  #--- FUNCTION read_molden ----------------------------------------------------
  #--- Load all information desired from a molden file -------------------------
  #-----------------------------------------------------------------------------
  
  fid    = open(filename,'r')		# Open the file
  flines = fid.readlines()		# Read the WHOLE file into RAM
  fid.close()				# Close the file
  
  #--- Is this really a molden file? ---
  if not '[Molden Format]\n' in flines:
    orbkit_output.display('The input file %s is no valid molden file!\nIt does'+
		 ' not contain the keyword: [Molden Format]\nEXIT\n' % filename)
    sys.exit(42)
  
  #--- Set a counter for the AOs ---
  basis_count = 0

  #--- Declare synonyms for molden keywords ---
  synonyms = {'Sym': 'sym',
              'Ene': 'energy',
              'Occup': 'occ_num',
             }
  MO_keys = synonyms.keys()
  
  #--- Go through the file line by line ---
  for il in range(len(flines)):
    line = flines[il]			# The current line as string
    thisline = line.split()		# The current line split into segments
    
    #--- Check the file for keywords ---
    if '[Molden Format]' in line:
      #--- A new file begins ---
      #--- Initialize the variables ---
      geo_info = []			# Information about atoms
      geo_spec = []			# Atom postitions
      ao_spec = []			# Information about atomic orbitals
      mo_spec = []			# Information about molecular orbitals
      sec_flag = False			# A Flag specifying the current section 
    elif '[Atoms]' in line:
      #--- The section containing information about ---
      #--- the molecular geometry begins ---
      sec_flag = 1
      if 'Angs' in line:
	#--- The length are given in Angstroem ---
	#--- and have to be converted to Bohr radii --
        aa_to_au = 1/0.52917720859
      else:
	#--- The length are given in Bohr radii ---
        aa_to_au = 1.0
    elif '[GTO]' in line:
      #--- The section containing information about ---
      #--- the atomic orbitals begins ---
      sec_flag = 2
      bNew = True			# Indication for start of new AO section
    elif '[MO]' in line:
      #--- The section containing information about ---
      #--- the molecular orbitals begins ---
      sec_flag = 3
      bNew = True			# Indication for start of new MO section
    elif '[STO]' in line:
      #--- The orbkit does not support Slater type orbitals ---
      orbkit_output.display('orbkit does not work for STOs!\nEXIT\n');
      sys.exit(42)
    else:
      #--- Check if we are in a specific section ---
      if sec_flag == 1:
	#--- Geometry section ---
	geo_info.append(thisline[0:3])
	geo_spec.append([float(ii)*aa_to_au for ii in thisline[3:]])
      if sec_flag == 2:
	#--- Atomic orbital section ---
	if thisline == []:
	  #--- There is a blank line after every AO ---
	  bNew = True
	elif bNew:
	  #--- The following AOs are for which atom? ---
	  bNew = False
	  at_num = int(thisline[0]) - 1
	  ao_num = 0
	elif len(thisline) == 3:
	  #--- AO information section ---
	  #--- Initialize a new dict for this AO ---
	  ao_num = 0			# Initialize number of atomic orbiatls 
	  ao_type = thisline[0]		# Which type of atomic orbital do we have
	  pnum = int(thisline[1])	# Number of primatives
	  #--- Calculate the degeneracy of this AO and increase basis_count ---
	  basis_count += l_deg(lquant[ao_type])
	  ao_spec.append({'atom': at_num,
			  'type': ao_type,
			  'pnum': pnum,
			  'coeffs': numpy.zeros((pnum, 2))
			  })
	else:
	  #--- Append the AO coefficients ---
	  coeffs = numpy.array(line.replace('D','e').split(), dtype=numpy.float64)
	  ao_spec[-1]['coeffs'][ao_num,:] = coeffs
	  ao_num += 1
      if sec_flag == 3:
	#--- Molecular orbital section ---
	if '=' in line:
	  #--- MO information section ---	  
	  if bNew:
	    #--- Create a numpy array for the MO coefficients ---
	    mo_spec.append({'coeffs': numpy.zeros(basis_count)})
	    bNew = False
	  #--- Append information to dict of this MO ---
	  info = line.replace('\n','').replace(' ','')
	  info = info.split('=')
	  if info[0] in MO_keys: 
	    if info[0] != 'Sym':
	      info[1] = float(info[1])
	    mo_spec[-1][synonyms[info[0]]] = info[1]
	else:
	  if ('[' or ']') in line:
	    sec_flag = 0 
	  else:
	    #--- Append the MO coefficients ---
	    bNew = True			# Reset bNew
	    index = int(thisline[0])-1
	    try: 
	      #--- Try to convert coefficient to float ---
	      mo_spec[-1]['coeffs'][index] = float(thisline[1])
	    except ValueError:
	      #--- If it cannot be converted print error message ---
	      orbkit_output.display(
		'Error in coefficient %d of MO %s!' % (index, mo_spec[-1]['sym']) + 
		'\nSetting this coefficient to zero...')
  
  #--- Are all MOs requested for the calculation? ---
  if not all_MO:
    mo_range = copy.deepcopy(mo_spec)
    mo_spec = []
    for ii_mo in mo_range:
      if ii_mo['occ_num'] > 0.0000001:
        mo_spec.append(ii_mo)
  
  return geo_spec, geo_info, ao_spec, mo_spec
  #--- read_molden ---

def read_gamess(filename, all_MO=False, loc='', ex_atom=''):
  #--- FUNCTION read_molden ----------------------------------------------------
  #--- Load all information desired from a molden file -------------------------
  #-----------------------------------------------------------------------------
  #--- Initialize the variables ---
  geo_info = []			# Information about atoms
  geo_spec = []			# Atom postitions
  ao_spec = []			# Information about atomic orbitals
  mo_spec = []			# Information about molecular orbitals
  occ = []
  AO = []
  sec_flag = False			# A Flag specifying the current section
  element_list = {}
  element_count = 0
  geo_skip = 0
  ao_skip = 0
  basis_count = 0
  blast = False

  #--- Go through the file line by line ---
  with open(filename) as fileobject:
    for line in fileobject:
      thisline = line.split()		# The current line split into segments
      #--- Check the file for keywords ---
      if ' ATOM      ATOMIC                      COORDINATES' in line:
	#--- The section containing information about ---
	#--- the molecular geometry begins ---
	sec_flag = 1
	geo_skip = 1
	atom_count = 1
	if '(BOHR)' in line:
	  #--- The length are given in Angstroem ---
	  #--- and have to be converted to Bohr radii --
	  aa_to_au = 1.0
	else:
	  #--- The length are given in Bohr radii ---
	  aa_to_au = 1/0.52917720859
	  
      elif 'ATOMIC BASIS SET' in line:
	#--- The section containing information about ---
	#--- the atomic orbitals begins ---
	sec_flag = 2
	ao_skip = 6
	at_type = []
      ##elif '          EIGENVECTORS' in line:
      elif LOC[loc.lower()] in line:
	#--- The section containing information about ---
	#--- the molecular orbitals begins ---
	sec_flag = 3
	if loc == '':
	  mo_skip = 1
	else:
	  mo_skip = 0
	init_mo = False			# Initialize new MO section
	mo_new = False			# Indication for start of new MO section
	ene = False
	len_mo = 0
      elif ' NUMBER OF OCCUPIED ORBITALS (ALPHA)          =' in line:
	occ.append(int(thisline[-1]))
      elif ' NUMBER OF OCCUPIED ORBITALS (BETA )          =' in line:
	occ.append(int(thisline[-1]))
      elif '          ECP POTENTIALS' in line:
	occ = []
      elif ' NUMBER OF OCCUPIED ORBITALS (ALPHA) KEPT IS    =' in line:
	occ.append(int(thisline[-1]))
      elif ' NUMBER OF OCCUPIED ORBITALS (BETA ) KEPT IS    =' in line:
	occ.append(int(thisline[-1]))
      else:
	#--- Check if we are in a specific section ---
	if sec_flag == 1:
	  if not geo_skip:
	    if line == '\n':
	      sec_flag = 0
	      a = numpy.array(geo_info)
	      #for jj in range(len(element_list)):
		#element_list[element_list.items()[jj][0]] = numpy.sum(a==element_list.items()[jj][0])
	    else:
	      if thisline[0] != ex_atom: 
	      #element_list[thisline[0]] = 0
		geo_info.append([thisline[0],atom_count,thisline[1]])
		geo_spec.append([float(ii)*aa_to_au for ii in thisline[2:]])
		atom_count += 1
	  elif geo_skip:
	    geo_skip -= 1
	    
	if sec_flag == 2:
	  if not ao_skip:
	    if ' TOTAL NUMBER OF BASIS SET SHELLS' in line:
	      sec_flag = 0
	    else:
	      if len(thisline) == 1:
		#--- Read atom type ---#
		at_type = thisline[0]
		if thisline[0] != ex_atom:
		  AO.append([])
		  new_ao = False
	      elif len(thisline) == 0 and new_ao == False and at_type != ex_atom:
		new_ao = True
	      elif len(thisline) == 5 and new_ao == True and at_type != ex_atom:
		AO[-1].append({'atom_type': at_type,
		'type': '',
		'pnum': 0,
		'coeffs': []})
		new_ao = False
		AO[-1][-1]['coeffs'].append([float(ii) for ii in thisline[3:]])
		AO[-1][-1]['type'] = thisline[1].lower()
		AO[-1][-1]['pnum'] += 1
	      elif len(thisline) == 5 and new_ao == False and at_type != ex_atom:
		AO[-1][-1]['coeffs'].append([float(ii) for ii in thisline[3:]])
		AO[-1][-1]['type'] = thisline[1].lower()
		AO[-1][-1]['pnum'] += 1
		
	  elif ao_skip:
	    ao_skip -= 1
	
	if sec_flag == 3:
	  if not mo_skip:
	    if END_LOC[loc.lower()] in line:
	      sec_flag = 0
	    else:
	      if thisline == []:
		if blast:
		  sec_flag = 0
		  blast = False
		blast = True
		init_mo = True
		mo_new = False
		ene = False
	      elif init_mo and not mo_new:
		init_len = len(thisline)
		for ii in range(len(thisline)):
		  #len_mo += 1
		  mo_spec.append({'coeffs': [],
		    'energy': 0.0,
		    'occ_num': 0.0,
		    'sym': ''
		  })
		  if loc != '':
		    #mo_spec[-1]['sym'] = '%d.%s' % (len_mo-1, flines[il+2].split()[ii])
		    #
		    #mo_skip = 2
		    mo_spec[-1]['sym'] = '%d.A1' % (len_mo-1)
		    mo_skip = 1
		init_mo = False
		mo_new = True
		blast = False
	      elif len(thisline) == init_len and ene == False:
		for ii in range(init_len,0,-1):
		  mo_spec[-ii]['energy'] = float(thisline[init_len-ii])
		ene = True
	      elif len(thisline) == init_len and ene == True:
		for ii in range(init_len,0,-1):
		  len_mo += 1
		  mo_spec[-ii]['sym'] = '%d.%s' % (len_mo-1, thisline[init_len-ii])
	      elif thisline != [] and not init_mo and mo_new:
		if ex_atom not in line[:16] or ex_atom == '':
		  for ii in range(init_len,0,-1):
		    mo_spec[-ii]['coeffs'].append(float(line[16:].split()[init_len-ii]))
	  elif mo_skip:
	    mo_skip -= 1


  #--- Check usage of same atomic basis sets ---#
  basis_set = {}
  for ii in range(len(AO)):
    if not AO[ii][0]['atom_type'] in basis_set.keys():
      basis_set[AO[ii][0]['atom_type']] = AO[ii]
    else:
      for jj in range(len(AO[ii])):
	if AO[ii][jj]['coeffs'] !=  basis_set[AO[ii][0]['atom_type']][jj]['coeffs']:
	  raise IOError('Different basis sets for the same atom.')

  #--- Numpy array ---#
  for ii in basis_set.iterkeys():
    for jj in range(len(basis_set[ii])):
      basis_set[ii][jj]['coeffs'] = numpy.array(basis_set[ii][jj]['coeffs'])

  for kk in range(len(mo_spec)):
    mo_spec[kk]['coeffs'] = numpy.array(mo_spec[kk]['coeffs'])

  #--- Complement atomic basis sets ---#
  for kk in range(len(geo_info)):
    for ll in range(len(basis_set[geo_info[kk][0]])):
      ao_spec.append({'atom': geo_info[kk][1]-1,
	'type': basis_set[geo_info[kk][0]][ll]['type'],
	'pnum': basis_set[geo_info[kk][0]][ll]['pnum'],
	'coeffs': basis_set[geo_info[kk][0]][ll]['coeffs']
	})

  for ii in range(len(mo_spec)):
    if occ[0] and occ[1]:
      mo_spec[ii]['occ_num'] += 2.0
      occ[0] -= 1
      occ[1] -= 1
    if not occ[0] and occ[1]:
      mo_spec[ii]['occ_num'] += 1.0
      occ[1] -= 1
      print 'occ_1'
    if not occ[1] and occ[0]:
      print 'occ_0'
      mo_spec[ii]['occ_num'] += 1.0
      occ[0] -= 1
    
  return geo_spec, geo_info, ao_spec, mo_spec


def read_gaussian_fchk(filename, all_MO=False):
  #--- FUNCTION read_gaussian_fchk ---------------------------------------------
  #--- Load all information desired from a gaussian FChk file ------------------
  #-----------------------------------------------------------------------------

  fid    = open(filename,'r')		# Open the file
  flines = fid.readlines()		# Read the WHOLE file into RAM
  fid.close()				# Close the file
  
  sec_flag = 0
  
  el_num = [0,0]
  index = 0
  at_num = 0
  ao_num = 0 
  switch = 0
  geo_info = []
  ao_spec = []
  mo_spec = []
  count_mo = 0
  
  #--- Set a counter for the AOs ---
  basis_count = 0
  
  #--- Go through the file line by line ---
  for il in range(len(flines)):
    line = flines[il]			# The current line as string
    thisline = line.split()		# The current line split into segments
    
    #--- Check the file for keywords ---
    if 'Number of alpha electrons' in line:
      el_num[0] = int(thisline[5]) 
    elif 'Number of beta electrons' in line:
      el_num[1] = int(thisline[5])
    #elif 'Number of basis functions' in line:
      #basis_count = int(thisline[5])
    elif 'Atomic numbers'  in line:
      sec_flag = 1
      index = 0
      at_num = int(thisline[-1])
      count = 0
      if geo_info == []:
	for ii in range(at_num):
	  geo_info.append(['',ii,''])
    elif 'Integer atomic weights' in line:
      sec_flag = 1
      index = 2
      at_num = int(thisline[-1])
      count = 0
      if geo_info == []:
	for ii in range(at_num):
	  geo_info.append(['',ii,''])	
    elif 'Current cartesian coordinates' in line:
      at_num = int(thisline[-1])/3
      sec_flag = 2
      geo_spec = []
      count = 0
      xyz = []
    elif 'Shell types' in line:
      sec_flag = 3
      index = 'type'
      ao_num = int(thisline[-1])
      count = 0
      if ao_spec == []:
	for ii in range(ao_num):
	  ao_spec.append({})
    elif 'Number of primitives per shell' in line:
      sec_flag = 3
      index = 'pnum'
      ao_num = int(thisline[-1])
      count = 0
      if ao_spec == []:
	for ii in range(ao_num):
	  ao_spec.append({})
    elif 'Shell to atom map' in line:
      sec_flag = 3
      index = 'atom'
      ao_num = int(thisline[-1])
      count = 0
      if ao_spec == []:
	for ii in range(ao_num):
	  ao_spec.append({})
    elif 'Primitive exponents' in line:
      sec_flag = 4
      ao_num = int(thisline[-1])
      count = 0
      switch = 0
      index = 0
      if ao_spec == []:
	raise IOError('Shell types need to be defined before the AO exponents!')
      if not 'coeffs' in ao_spec[0].keys():
	for ii in range(len(ao_spec)):
	  pnum = ao_spec[ii]['pnum']
	  ao_spec[ii]['coeffs'] = numpy.zeros((pnum, 2))
    elif 'Contraction coefficients' in line:
      sec_flag = 4
      ao_num = int(thisline[-1])
      count = 0
      switch = 1
      index = 0
      if ao_spec == []:
	raise IOError('Shell types need to be defined before the AO exponents!')
      if not 'coeffs' in ao_spec[0].keys():
	for ii in range(len(ao_spec)):
	  pnum = ao_spec[ii]['pnum']
	  ao_spec[ii]['coeffs'] = numpy.zeros((pnum, 2))
    elif 'Orbital Energies' in line:
      sec_flag = 5
      mo_num = int(thisline[-1])
      index = 0
      count = len(mo_spec)
      if el_num[0] == el_num[1]:
	i = el_num[0]
	occ = 2
      else:
	i = el_num[0 if 'Alpha' in line else 1]
	occ = 1
      for ii in range(mo_num):
	mo_spec.append({'coeffs': numpy.zeros(basis_count),
	  'energy': 0.0,
	  'occ_num': float(occ if ii < i else 0)
	})
	mo_spec[-1]['sym'] = '%i.1' % len(mo_spec)
    elif 'MO coefficients' in line:
      sec_flag = 6
      count = 0
      mo_num = int(thisline[-1])
    else:
      #--- Check if we are in a specific section ---
      if sec_flag == 1:
	for ii in thisline:
	  geo_info[count][index] = ii
	  count += 1
	  if count == at_num:
	    sec_flag = 0
      elif sec_flag == 2:
	for ii in thisline:
	  xyz.append(float(ii))
	  if len(xyz) == 3:
	    geo_spec.append(xyz)
	    xyz = []
	    count += 1
	    if count == at_num:
	      sec_flag = 0
      elif sec_flag == 3:
	for ii in thisline:
	  ii = int(ii)
	  if index is 'type':
	    ii = orbit[abs(ii)]
	    basis_count += l_deg(lquant[ii])
	  elif index is 'atom':
	    ii -= 1
	  ao_spec[count][index] = ii
	  count += 1
	  if count == ao_num:
	    sec_flag = 0
      elif sec_flag == 4:
	for ii in thisline:
	  ao_spec[index]['coeffs'][count,switch] = float(ii)
	  count += 1
	  ao_num -= 1
	  if count == ao_spec[index]['pnum']:
	    index += 1
	    count = 0
	if not ao_num:
	  sec_flag = 0
      elif sec_flag == 5:
	for ii in thisline:
	  mo_spec[count]['energy'] = float(ii)
	  count += 1
	  if index != 0 and not count % basis_count:
	    sec_flag = 0
      elif sec_flag == 6:
	for ii in thisline:
	  mo_spec[index]['coeffs'][count] = float(ii)
	  count += 1
	  if count == basis_count:
	    count = 0
	    index += 1
	  if index != 0 and not index % basis_count:
	    sec_flag = 0
	
	
  #--- Are all MOs requested for the calculation? ---
  if not all_MO:
    mo_range = copy.deepcopy(mo_spec)
    mo_spec = []
    for ii_mo in mo_range:
      if ii_mo['occ_num'] > 0.0000001:
        mo_spec.append(ii_mo)

  return geo_spec, geo_info, ao_spec, mo_spec
  #--- read_gaussian_fchk ---

def read_gaussian_log(filename):
  #--- FUNCTION read_gaussian_log ----------------------------------------------
  #--- Load all information desired from a gaussian log file -------------------
  #-----------------------------------------------------------------------------
  
  fid    = open(filename,'r')		# Open the file
  flines = fid.readlines()		# Read the WHOLE file into RAM
  fid.close()				# Close the file

  aa_to_au = 1/0.52917720859  # conversion factor for Angstroem to bohr radii

  sec_flag = 0;skip=0;orb_skip=0;
  ij=0;kk=0;
  geo_info = [];geo_spec = []; ao_spec = []; mo_spec = [];

  dummy_array = []; ao = [];

  ao_type=[];ao_pnum=[];mo_coeff=[]
  el_num=[0,0];bf_num=0
  dummy_num=0;geo_info_n=0;found_MO=True
  for stuff in flines:
    thisline = stuff.split() 
    if ('basis functions,' in stuff) == 1: bf_num = int(thisline[0])
    elif ('alpha electrons' in stuff) == 1: el_num = [int(thisline[0]), int(thisline[3])]
    
    #elif ('                          Input orientation:' in stuff) == 1: sec_flag = 1; skip=5;geo_info = []; geo_spec = []
    elif ('                         Standard orientation:' in stuff) == 1: sec_flag = 1; skip=5;geo_info = []; geo_spec = []
    
    elif ('AO basis set in the form of general basis input (Overlap normalization):' in stuff) == 1: sec_flag = 2; skip=1;dummy_num=0;
    elif ('     Alpha Molecular Orbital Coefficients:' in stuff) == 1: found_MO=False;sec_flag = 3; skip=1;ij=0; start=3; mo_coeff=numpy.array(numpy.zeros( (bf_num,bf_num) ))#FIXME MO-BETA?
    elif (found_MO and '     Molecular Orbital Coefficients:' in stuff) == 1: found_MO=False;sec_flag = 3; skip=1;ij=0; start=3; mo_coeff=numpy.array(numpy.zeros( (bf_num,bf_num) ))
    
    if skip == 0:
      if ('---------------------------------------------------------------------' in stuff) == 1: sec_flag = 0;
      elif   sec_flag == 1: 
	geo_info.append([thisline[1],thisline[0],thisline[1]]); 
	geo_spec.append([aa_to_au*float(ij) for ij in thisline[3:]]);
      elif   sec_flag == 2:  
	dummy_array.append(thisline)
	if(' ****' in stuff) == 1: 
	  ao.append(dummy_array[:-1]);
	  dummy_array=[];
	  dummy_num+=1;
	  if dummy_num == len(geo_info): sec_flag = 0; 
	  ###FIXME
	  
      elif   sec_flag == 3:
	if start == 0:
	  tmp_line=stuff[20:].split();
	  for jj in range(len(tmp_line)):
	    mo_coeff[kk+jj][ij]=float(tmp_line[jj])
	  ij+=1
	  if ij == bf_num:start=3;kk+=len(tmp_line);ij=0
	  if kk == bf_num: sec_flag = 0;
	else:start-=1	
    else: skip -= 1;

  # ao_spec
  k_flag = 1; kk = 0;
  for stuff in ao: 
    at_num = int(stuff[0][0]) - 1
    for ii in stuff[1:]:
      if k_flag == 1:
	ao_type = ii[0]
	ao_pnum = int(ii[1])  
	dummy_dict = {}
	dummy_dict['coeffs'] = numpy.array(numpy.zeros( (ao_pnum,2) )) 
	k_flag = 0
	kk=0     
      else:
	dummy_dict['coeffs'][kk][0] = float(ii[0].replace('D','e'))
	dummy_dict['coeffs'][kk][1] = float(ii[1].replace('D','e'))
	kk += 1
      if kk == ao_pnum:
	dummy_dict['pnum'] = ao_pnum
	dummy_dict['type'] = ao_type.lower()
	dummy_dict['atom'] = at_num
	ao_spec.append(dummy_dict)
	k_flag = 1
      
  #mo_spec
  for ij in range(bf_num):
    dummy_dict={}
    dummy_dict['coeffs'] = numpy.array(mo_coeff[ij])
    if ij<el_num[0] and ij<el_num[1]: occ_num=2.0
    elif ij<el_num[0] or ij<el_num[1]: occ_num=1.0
    else: occ_num = 0.0
    dummy_dict['occ_num'] = occ_num
    if occ_num != 0: mo_spec.append(dummy_dict)
    
  return geo_spec, geo_info, ao_spec, mo_spec
  #--- read_gaussian_log ---

def l_creator(geo_spec,ao_spec,sel_ao,exp_list=None,coeff_list=None,
	      x=None,y=None,z=None,N=None,vector=False):
  #--- FUNCTION l_creator ------------------------------------------------------
  #--- Calculate the contracted atomic orbitals of q.n. l ----------------------
  #--- for the AO: ao_spec[sel_ao] ---------------------------------------------
  #-----------------------------------------------------------------------------
  #http://docs.scipy.org/doc/numpy/user/c-info.python-as-glue.html#inline-c-code
  
  if x == None: x = grid.x
  if y == None: y = grid.y
  if z == None: z = grid.z
  
  if not vector:
    if N == None: N = numpy.array(grid.N_)
  else:
    if len(x) != len(y) or len(x) != len(z):
      orbkit_output.display("Dimensions of x-, y-, and z- coordinate differ!")
      return 0
    else:
      N = (len(x),)
  
  if exp_list == None:
    l = lquant[ao_spec[sel_ao]['type']]
    exp_list = numpy.array(exp[l])
  if coeff_list == None:
    coeff_list = numpy.array(ao_spec[sel_ao]['coeffs'])
  at_pos = numpy.array(geo_spec[ao_spec[sel_ao]['atom']])
  ao_list = numpy.zeros(((len(exp_list),) + tuple(N)))
  
  ao_num = numpy.shape(coeff_list)[0]
  if not vector:
    ao_code = '''
    double Norm[ao_num][Nao_list[0]];
    double X, Y, Z;
    int lx[Nao_list[0]], ly[Nao_list[0]], lz[Nao_list[0]];
    double rr, ao_l0[ao_num], ao_xyz[Nao_list[0]];
    

    for (int il=0; il<Nao_list[0]; il++)
    {
      lx[il] = EXP_LIST2(il,0);
      ly[il] = EXP_LIST2(il,1);
      lz[il] = EXP_LIST2(il,2);
      
      for (int ii=0; ii<ao_num; ii++)
      {
	Norm[ii][il] = ao_norm(lx[il],ly[il],lz[il],&COEFF_LIST2(ii,0));
      }
    }
    
    for (int i=0; i<Nx[0]; i++)
    {
      for (int j=0; j<Ny[0]; j++)
      {
	for (int k=0; k<Nz[0]; k++)
	{
	  X = x[i]-at_pos[0];
	  Y = y[j]-at_pos[1];
	  Z = z[k]-at_pos[2];
	  
	  rr = pow(X,2)+pow(Y,2)+pow(Z,2);
	  
	  for (int il=0; il<Nao_list[0]; il++)
	  {	
	    ao_xyz[il] = xyz(X, Y, Z, lx[il], ly[il], lz[il]);
	  }
	  
	  for (int ii=0; ii<ao_num; ii++)
	  {
	    ao_l0[ii] = COEFF_LIST2(ii,1) * exp(-COEFF_LIST2(ii,0) * rr);
	    
	    for (int il=0; il<Nao_list[0]; il++)
	    {
	      AO_LIST4(il,i,j,k) += Norm[ii][il] * ao_xyz[il] * ao_l0[ii];
	    }
	  }
	  
	}
      }
    }
    '''
  else:
    ao_code = '''
    double Norm[ao_num][Nao_list[0]];
    double X, Y, Z;
    int lx[Nao_list[0]], ly[Nao_list[0]], lz[Nao_list[0]];
    double rr, ao_l0[ao_num], ao_xyz[Nao_list[0]];
    

    for (int il=0; il<Nao_list[0]; il++)
    {
      lx[il] = EXP_LIST2(il,0);
      ly[il] = EXP_LIST2(il,1);
      lz[il] = EXP_LIST2(il,2);
      
      for (int ii=0; ii<ao_num; ii++)
      {
	Norm[ii][il] = ao_norm(lx[il],ly[il],lz[il],&COEFF_LIST2(ii,0));
      }
    }
    
    for (int i=0; i<Nx[0]; i++)
    {
      X = x[i]-at_pos[0];
      Y = y[i]-at_pos[1];
      Z = z[i]-at_pos[2];
      
      rr = pow(X,2)+pow(Y,2)+pow(Z,2);
      
      for (int il=0; il<Nao_list[0]; il++)
      {	
	ao_xyz[il] = xyz(X, Y, Z, lx[il], ly[il], lz[il]);
      }
      
      for (int ii=0; ii<ao_num; ii++)
      {
	ao_l0[ii] = COEFF_LIST2(ii,1) * exp(-COEFF_LIST2(ii,0) * rr);
	
	for (int il=0; il<Nao_list[0]; il++)
	{
	  AO_LIST2(il,i) += Norm[ii][il] * ao_xyz[il] * ao_l0[ii];
	}
      }
	  

    }
    '''
  weave.inline(ao_code, ['x','y','z','ao_num','exp_list',
	      'coeff_list','at_pos','ao_list'], 
	      support_code = cSupportCode.norm + cSupportCode.xyz,verbose = 1)

  return ao_list
  #--- l_creator ---

def l_creator_atom(at_pos=None,exp_list=None,coeff_list=None,
	      pnum_list=None,x=None,y=None,z=None,N=None):
  #--- FUNCTION l_creator_atom -------------------------------------------------
  #--- Calculate all contracted atomic orbitals of ONE selected atom -----------
  #-----------------------------------------------------------------------------
  #http://docs.scipy.org/doc/numpy/user/c-info.python-as-glue.html#inline-c-code
  
  if x == None: x = grid.x
  if y == None: y = grid.y
  if z == None: z = grid.z
  if N == None: N = numpy.array(grid.N_)
  
  #N = numpy.array(grid.N_)
  
  #--- All AOs from one atom
  if exp_list == None or coeff_list == None or numpy.shape(at_pos) != (3,):
    orbkit_output.display(
      "If all AOs for one atom are requested, the variables " +
      "exp_list and coeff_list and pnum_list \nhave to be specify explicitly.\n" +
      "Moreover, geo_spec has to contain only one set of [x,y,z] coordinates:\n\t" +
      "l_creator(at_pos=[x,y,z],exp_list=exp,coeff_list=coeff,pnum_list=pnum)")
    return 255
  else:
    at_pos = numpy.array(at_pos)
    x = x - at_pos[0]
    y = y - at_pos[1]
    z = z - at_pos[2]
    exp_list = numpy.array(exp_list)
    coeff_list = numpy.array(coeff_list)
    pnum = numpy.shape(coeff_list)[0]
    pnum_list = numpy.array(pnum_list)
    which_ao = numpy.zeros((pnum,))

    ii_p = 0
    for ii in range(len(pnum_list)):
      for ij in range(pnum_list[ii]):
	which_ao[ii_p] = ii
	ii_p += 1
  
  ao_list = numpy.zeros((len(exp_list),len(x),len(y),len(z)))  
  
  ao_code = '''
  double Norm[pnum][Nao_list[0]];
  double X, Y, Z;
  int il; 
  int lx[Nao_list[0]], ly[Nao_list[0]], lz[Nao_list[0]];
  double rr, ao_l0, ao_xyz;
  
  for (int ii=0; ii<pnum; ii++)
  {
    il = which_ao[ii];
    lx[il] = EXP_LIST2(il,0);
    ly[il] = EXP_LIST2(il,1);
    lz[il] = EXP_LIST2(il,2);
    Norm[ii][il] = ao_norm(lx[il],ly[il],lz[il],&COEFF_LIST2(ii,0));
  }
  
  for (int i=0; i<Nao_list[1]; i++)
  {
    for (int j=0; j<Nao_list[2]; j++)
    {
      for (int k=0; k<Nao_list[3]; k++)
      {
	X = x[i];
	Y = y[j];
	Z = z[k];
	
	rr = pow(X,2)+pow(Y,2)+pow(Z,2);
	
	
	for (int ii=0; ii<pnum; ii++)
	{
	  ao_l0 = COEFF_LIST2(ii,1) * exp(-COEFF_LIST2(ii,0) * rr);
	  
	  il = which_ao[ii];
	  ao_xyz = xyz(X, Y, Z, lx[il], ly[il], lz[il]);
	  AO_LIST4(il,i,j,k) += Norm[ii][il] * ao_xyz * ao_l0;

	}
	
      }
    }
  }
  '''

  weave.inline(ao_code, ['x','y','z','N','exp_list','coeff_list',
			 'pnum','which_ao','ao_list'], 
	      support_code = cSupportCode.norm + cSupportCode.xyz,verbose = 1)
  
  return ao_list
  #--- l_creator_atom ---
  

###################################################
# Derivative of AO

def delta_l_creator(geo_spec,ao_spec,sel_ao,exp_list=None,coeff_list=None,
		    drv='x',x=None,y=None,z=None,N=None):
  #--- FUNCTION delta_l_creator ------------------------------------------------
  #--- Calculate the derivative of the contracted atomic orbitals --------------
  #--- of q.n. l for the AO: ao_spec[sel_ao] -----------------------------------
  #--- with respect to a specific variable (e.g. drv = 'x' or drv = 0) ---------
  #-----------------------------------------------------------------------------
  #http://docs.scipy.org/doc/numpy/user/c-info.python-as-glue.html#inline-c-code
  
  if x == None: x = grid.x
  if y == None: y = grid.y
  if z == None: z = grid.z
  if N == None: N = numpy.array(grid.N_)
  
  #N = numpy.array(grid.N_)
  
  if exp_list == None:
    l = lquant[ao_spec[sel_ao]['type']]
    exp_list = numpy.array(exp[l])
  if coeff_list == None:
    coeff_list = numpy.array(ao_spec[sel_ao]['coeffs'])
  at_pos = numpy.array(geo_spec[ao_spec[sel_ao]['atom']])
  
  ao_list = numpy.zeros(((len(exp_list),) + tuple(N)))
  
  ao_num = numpy.shape(coeff_list)[0]

  #l = lquant[ao_spec[sel_ao]['type']]
  #l_deg = (l+1)*(l+2)/2
  #coeff_list = numpy.array(ao_spec[sel_ao]['coeffs'])
  #exp_list = numpy.array(exp[l])
  #at_pos = numpy.array(geo_spec[ao_spec[sel_ao]['atom']])
  #ao_list = numpy.zeros((l_deg,len(x),len(y),len(z)))
  #ao_num = numpy.shape(coeff_list)[0]

  # Calculate the derivative with respect to which variable
  #try: 
    #drv = int(drv)
  #except ValueError:
  if not isinstance(drv, (int, long)):
    drv = 'xyz'.find(drv)
  # Was it a valid selection? If not calculate derivative with respect to x...
  if drv == -1:
    drv = 0
    orbkit_output.display("The selection of the derivative variable was not valid!" +
		    " (drv = 'x' or 'y' or 'z')")
    orbkit_output.display("Calculating the derivative with respect to x...")	

  ao_code = '''
  double Norm[ao_num][Nao_list[0]];
  double X, Y, Z;
  int lx[Nao_list[0]], ly[Nao_list[0]], lz[Nao_list[0]];
  double rr, ao_l0[ao_num], ao_xyz;
  

  for (int il=0; il<Nao_list[0]; il++)
  {
    lx[il] = EXP_LIST2(il,0);
    ly[il] = EXP_LIST2(il,1);
    lz[il] = EXP_LIST2(il,2);
    
    for (int ii=0; ii<ao_num; ii++)
    {
      Norm[ii][il] = ao_norm(lx[il],ly[il],lz[il],&COEFF_LIST2(ii,0));
    }
  }
  
  for (int i=0; i<Nao_list[1]; i++)
  {
    for (int j=0; j<Nao_list[2]; j++)
    {
      for (int k=0; k<Nao_list[3]; k++)
      {
	X = x[i]-at_pos[0];
	Y = y[j]-at_pos[1];
	Z = z[k]-at_pos[2];
	
	rr = pow(X,2)+pow(Y,2)+pow(Z,2);
	
	for (int ii=0; ii<ao_num; ii++)
	{
	  ao_l0[ii] = COEFF_LIST2(ii,1) * exp(-COEFF_LIST2(ii,0) * rr);
	  
	  for (int il=0; il<Nao_list[0]; il++)
	  {
	    switch(drv)
	    {
	      case 0:
	      {
		if (lx[il] == 0)
		{
		  ao_xyz = - 2 * COEFF_LIST2(ii,0) * xyz(X, Y, Z, lx[il]+1, ly[il], lz[il]);
		}
		else
		{
		  ao_xyz = lx[il] * xyz(X, Y, Z, lx[il]-1, ly[il], lz[il]) - 2 * COEFF_LIST2(ii,0) * xyz(X, Y, Z, lx[il]+1, ly[il], lz[il]);	    
		}
	      } break;
	      
	      case 1:
	      {
		if (ly[il] == 0)
		{
		  ao_xyz = - 2 * COEFF_LIST2(ii,0) * xyz(X, Y, Z, lx[il], ly[il]+1, lz[il]);
		}
		else
		{
		  ao_xyz = ly[il] * xyz(X, Y, Z, lx[il], ly[il]-1, lz[il]) - 2 * COEFF_LIST2(ii,0) * xyz(X, Y, Z, lx[il], ly[il]+1, lz[il]);	    
		}
	      } break;
	      
	      case 2:
	      {
		if (lz[il] == 0)
		{
		  ao_xyz = - 2 * COEFF_LIST2(ii,0) * xyz(X, Y, Z, lx[il], ly[il], lz[il]+1);
		}
		else
		{
		  ao_xyz = lz[il] * xyz(X, Y, Z, lx[il], ly[il], lz[il]-1) - 2 * COEFF_LIST2(ii,0) * xyz(X, Y, Z, lx[il], ly[il], lz[il]+1);	    
		}
	      } break;
	      
	      default:
	      {
		std::cout << "False statement for derivative variable!" << std::endl;
	      }
	    }
	    AO_LIST4(il,i,j,k) += Norm[ii][il] * ao_xyz * ao_l0[ii];
	  }
	}
	
      }
    }
  }
  
  '''

  weave.inline(ao_code, ['x','y','z','N','ao_num','l_deg','exp_list',
	      'coeff_list','at_pos','ao_list','drv'], 
	      support_code = cSupportCode.norm + cSupportCode.xyz,verbose = 1)

  return ao_list
  #--- delta_l_creator ---

def delta_l_creator_atom(at_pos=None,exp_list=None,coeff_list=None,
		    pnum_list=None,drv='x',x=None,y=None,z=None,N=None):
  #--- FUNCTION delta_l_creator_atom -------------------------------------------
  #--- Calculate the derivative of ALL contracted atomic orbitals --------------
  #--- of ONE selected atom ----------------------------------------------------
  #--- with respect to a specific variable (e.g. drv = 'x' or drv = 0) ---------
  #-----------------------------------------------------------------------------
  #http://docs.scipy.org/doc/numpy/user/c-info.python-as-glue.html#inline-c-code
  
  if x == None: x = grid.x
  if y == None: y = grid.y
  if z == None: z = grid.z
  if N == None: N = numpy.array(grid.N_)
  
  #N = numpy.array(grid.N_)

  # Calculate the derivative with respect to which variable
  #try: 
    #drv = int(drv)
  #except ValueError:
  if not isinstance(drv, (int, long)):
    drv = 'xyz'.find(drv)
  # Was it a valid selection? If not calculate derivative with respect to x...
  if drv == -1:
    drv = 0
    orbkit_output.display("The selection of the derivative variable was not valid!" +
		    " (drv = 'x' or 'y' or 'z')")
    orbkit_output.display("Calculating the derivative with respect to x...")	


  #--- All AOs from one atom
  if exp_list == None or coeff_list == None or numpy.shape(at_pos) != (3,):
    orbkit_output.display(
      "If all AOs for one atom are requested, the variables " +
      "exp_list and coeff_list and pnum_list \nhave to be specify explicitly.\n" +
      "Moreover, geo_spec has to contain only one set of [x,y,z] coordinates:\n\t" +
      "l_creator(at_pos=[x,y,z],exp_list=exp,coeff_list=coeff,pnum_list=pnum)")
    return 255
  else:
    at_pos = numpy.array(at_pos)
    x = x - at_pos[0]
    y = y - at_pos[1]
    z = z - at_pos[2]
    exp_list = numpy.array(exp_list)
    coeff_list = numpy.array(coeff_list)
    pnum = numpy.shape(coeff_list)[0]
    pnum_list = numpy.array(pnum_list)
    which_ao = numpy.zeros((pnum,))

    ii_p = 0
    for ii in range(len(pnum_list)):
      for ij in range(pnum_list[ii]):
	which_ao[ii_p] = ii
	ii_p += 1

  ao_list = numpy.zeros((len(exp_list),len(x),len(y),len(z)))  


  ao_code = '''
  double Norm[pnum][Nao_list[0]];
  double X, Y, Z;
  int il; 
  int lx[Nao_list[0]], ly[Nao_list[0]], lz[Nao_list[0]];
  double rr, ao_l0, ao_xyz;
  
  for (int ii=0; ii<pnum; ii++)
  {
    il = which_ao[ii];
    lx[il] = EXP_LIST2(il,0);
    ly[il] = EXP_LIST2(il,1);
    lz[il] = EXP_LIST2(il,2);
    Norm[ii][il] = ao_norm(lx[il],ly[il],lz[il],&COEFF_LIST2(ii,0));
  }
  
  for (int i=0; i<Nao_list[1]; i++)
  {
    for (int j=0; j<Nao_list[2]; j++)
    {
      for (int k=0; k<Nao_list[3]; k++)
      {
	X = x[i];
	Y = y[j];
	Z = z[k];
	
	rr = pow(X,2)+pow(Y,2)+pow(Z,2);
	
	for (int ii=0; ii<pnum; ii++)
	{
	  ao_l0 = COEFF_LIST2(ii,1) * exp(-COEFF_LIST2(ii,0) * rr);
	  
	  il = which_ao[ii];
	  switch(drv)
	  {
	    case 0:
	    {
	      if (lx[il] == 0)
	      {
		ao_xyz = - 2 * COEFF_LIST2(ii,0) * xyz(X, Y, Z, lx[il]+1, ly[il], lz[il]);
	      }
	      else
	      {
		ao_xyz = lx[il] * xyz(X, Y, Z, lx[il]-1, ly[il], lz[il])  - 2 * COEFF_LIST2(ii,0) * xyz(X, Y, Z, lx[il]+1, ly[il], lz[il]);	    
	      }
	    } break;
	    
	    case 1:
	    {
	      if (ly[il] == 0)
	      {
		ao_xyz = - 2 * COEFF_LIST2(ii,0) * xyz(X, Y, Z, lx[il], ly[il]+1, lz[il]);
	      }
	      else
	      {
		ao_xyz = ly[il] * xyz(X, Y, Z, lx[il], ly[il]-1, lz[il])  - 2 * COEFF_LIST2(ii,0) * xyz(X, Y, Z, lx[il], ly[il]+1, lz[il]);	    
	      }
	    } break;
	    
	    case 2:
	    {
	      if (lz[il] == 0)
	      {
		ao_xyz = - 2 * COEFF_LIST2(ii,0) * xyz(X, Y, Z, lx[il], ly[il], lz[il]+1);
	      }
	      else
	      {
		ao_xyz = lz[il] * xyz(X, Y, Z, lx[il], ly[il], lz[il]-1)  - 2 * COEFF_LIST2(ii,0) * xyz(X, Y, Z, lx[il], ly[il], lz[il]+1);	    
	      }
	    } break;
	    
	    default:
	    {
	      std::cout << "False statement for derivative variable!" << std::endl;
	    }
	  }
	  AO_LIST4(il,i,j,k) += Norm[ii][il] * ao_xyz * ao_l0;
	
	}
	
      }
    }
  }
  
  '''

  weave.inline(ao_code, ['x','y','z','N','at_pos','exp_list','coeff_list',
			 'pnum','which_ao','ao_list','drv'], 
	      support_code = cSupportCode.norm + cSupportCode.xyz,verbose = 1)

  return ao_list
  #--- delta_l_creator_atom ---

###################################################
def ao_creator(geo_spec,ao_spec,x=None,y=None,z=None,N=None,vector=False,exp_list=False):
  #--- FUNCTION ao_creator -----------------------------------------------------
  #--- Calculate the contracted atomic orbitals --------------------------------
  if x == None: x = grid.x
  if y == None: y = grid.y
  if z == None: z = grid.z
  if N == None: N = numpy.array(grid.N_)
  
  if not exp_list:
    ii_exp = None
  
  # Generalized AO creator
  for ii in range(len(ao_spec)):
    if exp_list:
      ii_exp = numpy.array([ao_spec[ii]['Exponents']])
    
    if not ii:
      ao_list = l_creator(geo_spec,ao_spec,ii,
			      x=x,y=y,z=z,N=N,vector=vector,exp_list=ii_exp)
    else:
      ao_list = numpy.append(ao_list,l_creator(geo_spec,ao_spec,ii,
				x=x,y=y,z=z,N=N,vector=vector,exp_list=ii_exp)	
				, axis = 0)
    
  return ao_list
  #--- ao_creator ---

def delta_ao_creator(geo_spec,ao_spec,drv='x',x=None,y=None,z=None,exp_list=None):
  #--- FUNCTION delta_ao_creator -----------------------------------------------
  #--- Calculate the derivative of the contracted atomic orbitals --------------
  #--- with respect to a specific variable (e.g. drv = 'x' or drv = 0) ---------
  if x == None: x = grid.x
  if y == None: y = grid.y
  if z == None: z = grid.z
  
  if not exp_list:
    ii_exp = None
  
  # Generalized delta AO creator
  for ii in range(len(ao_spec)):
    if exp_list:
      ii_exp = numpy.array([ao_spec[ii]['Exponents']])
    if not ii:
      ao_list = delta_l_creator(geo_spec,ao_spec,ii,
				drv=drv,x=x,y=y,z=z,exp_list=ii_exp)
    else:
      ao_list = numpy.append(ao_list,delta_l_creator(geo_spec,ao_spec,ii,
			      drv=drv,x=x,y=y,z=z,exp_list=ii_exp), axis = 0)
  
  return ao_list
  #--- delta_ao_creator ---

def calc_single_mo(zz):
  try:
    ao_list = Spec['ao_list']
    mo_spec = Spec['mo_spec']
    N = Spec['N']
    mo = numpy.zeros(N)
    for jj in range(len(ao_list)):
      mo += mo_spec[zz]['coeffs'][jj] * ao_list[jj]
    return mo
  except KeyboardInterrupt:
    #--- Catch keybord interrupt signal to prevent a hangup of the worker processes ---
    return 0

def mo_creator(ao_list,mo_spec,x=None,y=None,z=None,N=None,vector=False,HDF5_save=False,h5py=False,s=0):
  #--- FUNCTION mo_creator------------------------------------------------------
  #--- Calculate the molecular orbitals ----------------------------------------
  #-----------------------------------------------------------------------------
  #--- If for HDF5_save a string (filenmae) is given the slice s of each -------
  #--- molecular orbital will be saved to the disk. This module is used by -----
  #--- orbkit_extras. Here the HDF5 file is initalized ----------------------------
  #-----------------------------------------------------------------------------
  
  if x == None: x = grid.x
  if y == None: y = grid.y
  if z == None: z = grid.z
  
  if not vector:
    if N == None: N = tuple(grid.N_)
  else:
    if len(x) != len(y) or len(x) != len(z):
      orbkit_output.display("Dimensions of x-, y-, and z- coordinate differ!")
      return 0
    else:
      N = (len(x),)
  
  if not HDF5_save:
    #--- Standard mo_creator ---
    mo_list = []
    for ii in range(len(mo_spec)):
      mo_list.append(numpy.zeros(N))
      for jj in range(len(ao_list)):
	mo_list[ii] += mo_spec[ii]['coeffs'][jj] * ao_list[jj]
    return mo_list
  else:
    
    f = h5py.File(HDF5_save, 'a')
    
    global Spec    
    Spec = {'ao_list': ao_list, 'mo_spec': mo_spec, 'N': N}
    
    if options.numproc > len(mo_spec): options.numproc = len(mo_spec)
    
    #--- Start the worker processes --
    pool = Pool(processes=options.numproc)
    
    #--- Write the slices in x to an array zz ---
    zz=[]
    for ii in range(len(mo_spec)):
      zz.append(ii)
    
    it = pool.imap(calc_single_mo, zz)
    
    
    
    for ii in range(len(mo_spec)):
      dID = 'MO:%s' % mo_spec[ii]['sym']
      a = it.next()
      f[dID][s,:,:] = a
    
    #--- Close the worker processes --
    pool.close()
    
    f.close()
    return 0
  #--- mo_creator ---

def mo_creator_matrop(ao_list,mo_coeff,AO_Nr,x=None,y=None,z=None,vector=False):
  #--- FUNCTION mo_creator------------------------------------------------------
  #--- Calculate the molecular orbitals ----------------------------------------
  if x == None: x = grid.x
  if y == None: y = grid.y
  if z == None: z = grid.z
  
  if not vector:
    N = tuple(grid.N_)
  else:
    if len(x) != len(y) or len(x) != len(z):
      orbkit_output.display("Dimensions of x-, y-, and z- coordinate differ!")
      return 0
    else:
      N = (len(x),)
  
  mo_list = []
  i = 0
  for ii_s in range(len(mo_coeff)):
    shape = numpy.shape(mo_coeff[ii_s])
    for ii in range(shape[0]):
      mo_list.append(numpy.zeros(N))
      for jj in range(shape[1]):
        for jj_ao in range(len(ao_list)):
	  if AO_Nr[jj_ao] == (jj+1,ii_s+1):
	    mo_list[i] += mo_coeff[ii_s][ii,jj] * ao_list[jj_ao]
      i += 1
  return mo_list
  #--- mo_creator_matrop ---

###################################################
def ao_creator2(geo_spec,ao_spec,fid='tmp',x=None,y=None,z=None):
  #--- FUNCTION ao_creator2 ----------------------------------------------------
  #--- The atomic orbital creator for save_ao ----------------------------------
  #-----------------------------------------------------------------------------
  #--- save_ao: Save all AOs to disk and reload for every MO calculation -------
  if x == None: x = grid.x
  if y == None: y = grid.y
  if z == None: z = grid.z
  
  HDF5_f = h5py.File('%(f)s_ao.h5' % {'f': fid},'w')
  shape = (len(x),len(y),len(z))
  ii_ao = 0
  # Generalized AO creator
  for ii in range(len(ao_spec)):
    ao_list = l_creator(geo_spec,ao_spec,ii,x=x,y=y,z=z)    
    for kk in range(len(ao_list)):
      #--- Save the atomic orbitals as hdf5 file ---
      dset_id = 'ao:%(i)03d' % {'i': ii_ao}
      dset = HDF5_f.create_dataset(dset_id, shape, data=ao_list[kk])
      #with h5py.File(filename) as f: 
	#f['/ao'] = ao_list[kk]  
      ii_ao += 1
  HDF5_f.close()
  #--- Store and return the number of atomic orbitals ---
  ao_num = ii_ao
  
  return ao_num
  #--- ao_creator2 ---

def mo_creator2(ao_num,mo_spec,fid='tmp',x=None,y=None,z=None):
  #--- FUNCTION mo_creator2 ----------------------------------------------------
  #--- The molecular orbital creator for save_ao -------------------------------
  #-----------------------------------------------------------------------------
  #--- save_ao: Save all AOs to disk and reload for every MO calculation -------
  if x == None: x = grid.x
  if y == None: y = grid.y
  if z == None: z = grid.z
  
  HDF5_f = h5py.File('%(f)s_ao.h5' % {'f': fid},'r')
  mo_list = []
  for ii in range(len(mo_spec)):
    mo_list.append(numpy.zeros((len(x),len(y),len(z))));
  
  mo_count = 0
  for ii_ao in range(ao_num): 
    dset_id = 'ao:%(i)03d' % {'i': ii_ao}
    ao_list = HDF5_f[dset_id]
    #filename = '%(f)s_ao_%(i)03d.h5' % {'f': fid, 'i': ii_ao}
    #with h5py.File(filename) as f: 
      #ao_list = f['/ao'].value 
    for ii in range(len(mo_spec)):
      mo_list[ii] += mo_spec[ii]['coeffs'][mo_count] * ao_list
    mo_count += 1
  HDF5_f.close()
  return mo_list
  #--- mo_creator2 ---

###################################################  
def mo_creator3(geo_spec,ao_spec,mo_spec,x=None,y=None,z=None):
  #--- FUNCTION mo_creator3 ----------------------------------------------------
  #--- The molecular orbital creator for discard_ao ----------------------------
  #-----------------------------------------------------------------------------
  #--- discard_ao: Calculate MOs without keeping the AOs in the RAM ------------
  if x == None: x = grid.x
  if y == None: y = grid.y
  if z == None: z = grid.z
  
  mo_list = []
  for ii in range(len(mo_spec)):
    mo_list.append(numpy.zeros((len(x),len(y),len(z))))

  mo_count = 0 
  # Generalized AO creator
  for ii in range(len(ao_spec)):
    ao_list = l_creator(geo_spec,ao_spec,ii,x=x,y=y,z=z)
    
    for kk in range(len(ao_list)):
      for jj in range(len(mo_spec)):
	mo_list[jj] += mo_spec[jj]['coeffs'][mo_count] * ao_list[kk]
      mo_count += 1
              
  return mo_list
  #--- mo_creator3 ---

def slice_rho(x):
  #--- FUNCTION slice_rho ------------------------------------------------------
  #--- Calculate the density or the derivative of the density with respect -----
  #--- to Spec['Derivative'] for one slice in x-direction ----------------------
  #-----------------------------------------------------------------------------
  try:
    #--- All desired information is stored in the Global variable Spec ---
    geo_spec = Spec['geo_spec']
    ao_spec = Spec['ao_spec']
    mo_spec = Spec['mo_spec']
    drv = Spec['Derivative']
    mo = 'mo' in Spec.keys()
    
    #--- Calculate the MOs and AOs for this slice ---
    if options.save_ao:
      fid = '%(f)s_x_%(i)0.4f' % {'f': options.outputname, 'i': x} 
      ao_count = ao_creator2(geo_spec,ao_spec,fid=fid,x=x)
      mo_list = mo_creator2(ao_count,mo_spec,fid=fid,x=x)
    elif options.discard_ao:
      mo_list = mo_creator3(geo_spec,ao_spec,mo_spec,x=x)
    else:
      ao_list = ao_creator(geo_spec,ao_spec,x=x)
      mo_list = mo_creator(ao_list,mo_spec,x=x)
    
    if mo:
      return numpy.array(mo_list)
    
    #--- Initialize a numpy array for the density ---
    rho = numpy.zeros((len(x),len(grid.y),len(grid.z)))
    
    #--- Initialize a numpy array for the norm of the MOs ---
    mo_norm = numpy.zeros((len(mo_list)))
    
    #--- Calculate the density and the norm ---
    for ii_mo in range(len(mo_list)): 
      mo_norm[ii_mo] = numpy.sum(numpy.square(mo_list[ii_mo]))
      rho += mo_spec[ii_mo]['occ_num'] * numpy.square(numpy.abs(mo_list[ii_mo]))
      
    if drv == None:
      #--- Return the density and the norm ---
      return rho, mo_norm
    else:
      #--- Initialize a numpy array for the derivative of the density ---
      delta_rho = numpy.zeros((3,len(x),len(grid.y),len(grid.z)))
      for ii_d in drv:
	#--- Calculate the derivatives of the AOs and MOs for this slice ---
	delta_ao_list = delta_ao_creator(geo_spec,ao_spec,drv=ii_d,x=x)
	delta_mo_list = mo_creator(delta_ao_list,mo_spec,x=x)
	
	#--- Calculate the derivative of the density
	for ii_mo in range(len(delta_mo_list)): 
	  delta_rho[ii_d,:,:,:] += (mo_spec[ii_mo]['occ_num'] * 
			    2 * delta_mo_list[ii_mo]*mo_list[ii_mo])
      #--- Return the derivative of the density ---
      #print numpy.shape(rho),numpy.shape(mo_norm),numpy.shape(delta_rho)
      return rho, mo_norm, delta_rho
  except KeyboardInterrupt:
    #--- Catch keybord interrupt signal to prevent a hangup of the worker processes ---
    return 0
  #--- slice_rho ---

def rho_compute(geo_spec,geo_info,ao_spec,mo_spec,calc_mo=False):
  #--- FUNCTION rho_compute ----------------------------------------------------
  #--- Function for the computation of the density -----------------------------
  #-----------------------------------------------------------------------------
  #--- The 3D grid is divided into Slices, and the computational tasks ---------
  #--- are distributed to the worker processes ---------------------------------
  #-----------------------------------------------------------------------------
  
  #--- Specify the global variable containing all desired information needed ---
  #--- by the function slice_rho ---
  global Spec
  
  Spec = {'geo_spec': geo_spec, 'ao_spec': ao_spec, 
	      'mo_spec': mo_spec, 'Derivative': None}
  
    
  
  #Spec = {'geo_spec': geo_spec, 'ao_spec': ao_spec, 'mo_spec': mo_spec}
  
  #--- Save the complete initial grid ---
  x = grid.x
  y = grid.y
  z = grid.z
  N = copy.deepcopy(grid.N_)
  Min = copy.deepcopy(grid.min_)
  Max = copy.deepcopy(grid.max_)
  
  #--- Make slices ---
  rho = numpy.zeros((len(grid.x),len(grid.y),len(grid.z)))
  sDim = 0
  sNum = grid.N_[sDim]
  grid.N_[sDim] = 1
  
  #--- The number of worker processes is capped to the number of ---
  #--- grid points in x-direction. --- 
  if options.numproc > sNum: options.numproc = sNum
  
  #--- Print information regarding the density calculation ---
  orbkit_output.display("\nStarting the density calculation...")
  orbkit_output.display("The grid has been seperated in 2d slices and the")
  if options.numproc == 1:
    orbkit_output.display("calculation will be carried out with 1 subprocess.\n" + 
    "\n\tThe number of subprocesses can be changed with -p\n")
  else:
    orbkit_output.display("calculation will be carried out with %d subprocesses." 
		  % options.numproc)
  orbkit_output.display("\nThere are %d contracted AOs and %d MOs to be calculated."
		  % (len(mo_spec[0]['coeffs']), len(mo_spec)))
  
  #--- Initialize some additional user information ---
  status_old = 0
  s_old = 0
  t = [time.time()]
  
  #--- Initialize an array to store the norm of the MOs ---
  mo_norm = numpy.zeros((len(mo_spec)))
  
  if calc_mo:
    Spec['mo'] = True
    mo_list = numpy.zeros((len(mo_spec),) + tuple(N))
  
  #--- Start the worker processes --
  pool = Pool(processes=options.numproc)
  
  #--- Write the slices in x to an array zz ---
  zz=[]
  for s in range(sNum):
    zz.append((numpy.array([x[s]])))
  
  #--- Compute the density slice by slice ---
  it = pool.imap(slice_rho, zz)
  for s in range(sNum):
    if calc_mo:
      mo_list[:,s,:,:] = it.next()[:,0,:,:]
	
    else:
      rho[s,:,:], norm = it.next()
      mo_norm += norm
    
    #--- Print out the progress of the computation ---
    status = round(s*100/float(sNum))
    if not status % 20 and status != status_old: 
      t.append(time.time())
      orbkit_output.display("\tFinished %(f)d%% (%(s)d slices in %(t).3fs)" 
		      % {'f': status,
			 's': s + 1 - s_old,
			 't': t[-1]-t[-2]})
      status_old = status
      s_old = s + 1
  
  #--- Close the worker processes --
  pool.close()
  
  #--- Restore the complete initial grid ---
  grid.N_ = N
  grid.min_ = Min
  grid.max_ = Max
  grid.grid_init()
  
  if calc_mo:
    return mo_list
  
  #--- Print the norm of the MOs ---
  orbkit_output.display("\nNorm of the MOs:")
  for ii_mo in range(len(mo_norm)):
    orbkit_output.display("\t%(m).6f\tMO %(n)s" 
		    % {'m':mo_norm[ii_mo]*grid.d3r, 'n':mo_spec[ii_mo]['sym']})
  
  #--- Print the number of electrons ---
  orbkit_output.display("We have " + str(numpy.sum(rho)*grid.d3r) + " electrons.")
  
  return rho
  #--- rho_compute ---

####
def slice_rho_vector(zz):
  #--- FUNCTION slice_rho ------------------------------------------------------
  #--- Calculate the density or the derivative of the density with respect -----
  #--- to Spec['Derivative'] for one slice in x-direction ----------------------
  #-----------------------------------------------------------------------------
  try:
    #--- All desired information is stored in the Global variable Spec ---
    geo_spec = Spec['geo_spec']
    ao_spec = Spec['ao_spec']
    mo_spec = Spec['mo_spec']
    drv = Spec['Derivative']
    
    # Set up Grid
    x = grid.x[zz[0]:zz[1]]
    y = grid.y[zz[0]:zz[1]]
    z = grid.z[zz[0]:zz[1]]
    
    
    #--- Calculate the MOs and AOs for this slice ---
    if options.save_ao:
      fid = '%(f)s_sl_%(i)03d' % {'f': options.outputname, 'i': zz[0]} 
      ao_count = ao_creator2(geo_spec,ao_spec,fid=fid,x=x,y=y,z=z,vector=True)
      mo_list = mo_creator2(ao_count,mo_spec,fid=fid,x=x,y=y,z=z,vector=True)
    elif options.discard_ao:
      mo_list = mo_creator3(geo_spec,ao_spec,mo_spec,x=x,y=y,z=z,vector=True)
    else:
      ao_list = ao_creator(geo_spec,ao_spec,x=x,y=y,z=z,vector=True)
      mo_list = mo_creator(ao_list,mo_spec,x=x,y=y,z=z,vector=True)

    #--- Initialize a numpy array for the density ---
    rho = numpy.zeros((len(x),))
    
    #--- Initialize a numpy array for the norm of the MOs ---
    mo_norm = numpy.zeros((len(mo_list)))
    
    #--- Calculate the density and the norm ---
    for ii_mo in range(len(mo_list)): 
      mo_norm[ii_mo] = numpy.sum(numpy.square(mo_list[ii_mo]))
      rho += mo_spec[ii_mo]['occ_num'] * numpy.square(numpy.abs(mo_list[ii_mo]))
      
    if drv == None:
      #--- Return the density and the norm ---
      return rho, mo_norm
    else:
      #--- Initialize a numpy array for the derivative of the density ---
      delta_rho = numpy.zeros((3,len(x)))
      for ii_d in drv:
	#--- Calculate the derivatives of the AOs and MOs for this slice ---
	delta_ao_list = delta_ao_creator(geo_spec,ao_spec,drv=ii_d,x=x,y=y,z=z,vector=True)
	delta_mo_list = mo_creator(delta_ao_list,mo_spec,x=x,y=y,z=z,vector=True)
	
	#--- Calculate the derivative of the density
	for ii_mo in range(len(delta_mo_list)): 
	  delta_rho[ii_d,:] += (mo_spec[ii_mo]['occ_num'] * 
			    2 * delta_mo_list[ii_mo]*mo_list[ii_mo])
      #--- Return the derivative of the density ---
      #print numpy.shape(rho),numpy.shape(mo_norm),numpy.shape(delta_rho)
      return rho, mo_norm, delta_rho
  except KeyboardInterrupt:
    #--- Catch keybord interrupt signal to prevent a hangup of the worker processes ---
    return 0
  #--- slice_rho ---

def rho_compute_vector(geo_spec,geo_info,ao_spec,mo_spec,delta_slice=1e4):
  #--- FUNCTION rho_compute ----------------------------------------------------
  #--- Function for the computation of the density -----------------------------
  #-----------------------------------------------------------------------------
  #--- The 3D grid is divided into Slices, and the computational tasks ---------
  #--- are distributed to the worker processes ---------------------------------
  #-----------------------------------------------------------------------------
  
  #--- Specify the global variable containing all desired information needed ---
  #--- by the function slice_rho ---
  global Spec
  
  Spec = {'geo_spec': geo_spec, 'ao_spec': ao_spec, 
	      'mo_spec': mo_spec, 'Derivative': None}
  #Spec = {'geo_spec': geo_spec, 'ao_spec': ao_spec, 'mo_spec': mo_spec}
  
  if len(grid.x) != len(grid.y) or len(grid.x) != len(grid.z):
      orbkit_output.display("Dimensions of x-, y-, and z- coordinate differ!")
      sys.exit(99)
  
  N = len(grid.x)
  sNum = int(numpy.floor(N/delta_slice)+1)
  #--- Make slices ---
  
  rho = numpy.zeros((N,))
  
  #--- Print information regarding the density calculation ---
  orbkit_output.display("\nStarting the density calculation...")
  
  #--- Initialize some additional user information ---
  status_old = 0
  s_old = 0
  t = [time.time()]
  
  #--- Initialize an array to store the norm of the MOs ---
  mo_norm = numpy.zeros((len(mo_spec)))
  
  #--- Start the worker processes --
  pool = Pool(processes=options.numproc)
  
  #--- Write the slices in x to an array zz ---
  zz=[]
  i = 0 
  
  for s in range(sNum):
    if i+delta_slice >= N:
      zz.append((numpy.array([i,N+1],dtype=int)))      
    else:
      zz.append((numpy.array([i,i+delta_slice],dtype=int)))
    i += delta_slice 
  
  #--- Compute the density slice by slice ---
  it = pool.imap(slice_rho_vector, zz)
  for s in range(sNum):
    rho[zz[s][0]:zz[s][1]], norm = it.next()
    mo_norm += norm
    
    #--- Print out the progress of the computation ---
    status = round(s*100/float(sNum))
    if not status % 20 and status != status_old: 
      t.append(time.time())
      orbkit_output.display("\tFinished %(f)d%% (%(s)d slices in %(t).3fs)" 
		      % {'f': status,
			 's': s + 1 - s_old,
			 't': t[-1]-t[-2]})
      status_old = status
      s_old = s + 1
  
  #--- Close the worker processes --
  pool.close()
  
  #--- Print the norm of the MOs ---
  orbkit_output.display("\nNorm of the MOs:")
  for ii_mo in range(len(mo_norm)):
    orbkit_output.display("\t%(m).6f\tMO %(n)s" 
		    % {'m':mo_norm[ii_mo]*grid.d3r, 'n':mo_spec[ii_mo]['sym']})
  
  #--- Print the number of electrons ---
  orbkit_output.display("We have " + str(numpy.sum(rho)*grid.d3r) + " electrons.")
  
  return rho
  #--- rho_compute ---


####

def rho_compute_no_slice(geo_spec,geo_info,ao_spec,mo_spec,vector=False):
  #--- FUNCTION rho_compute_no_slice -------------------------------------------
  #--- Function for the computation of the density without slicing the grid ----
  #-----------------------------------------------------------------------------
  
  #--- Calculate the AOs and MOs ---
  if options.save_ao:
    ao_count = ao_creator2(geo_spec,ao_spec,fid=options.outputname)
    mo_list = mo_creator2(ao_count,mo_spec,fid=options.outputname)
  elif options.discard_ao:
    mo_list = mo_creator3(geo_spec,ao_spec,mo_spec)
  else:
    ao_list = ao_creator(geo_spec,ao_spec,vector=vector)
    mo_list = mo_creator(ao_list,mo_spec,vector=vector)
  
  if not vector:
    N = tuple(grid.N_)
  else:
    if len(grid.x) != len(grid.y) or len(grid.x) != len(grid.z):
      orbkit_output.display("Dimensions of x-, y-, and z- coordinate differ!")
      return 0
    else:
      N = (len(grid.x),)
  
  #--- Initialize a numpy array for the density
  rho = numpy.zeros(N)
  
  #--- Calculate the density and the norm
  for ii_mo in range(len(mo_list)): 
    rho += numpy.square(numpy.abs(mo_list[ii_mo])) * mo_spec[ii_mo]['occ_num']
  
  #--- Print the norm of the MOs ---
  orbkit_output.display("Norm of the MOs:")
  for ii_mo in range(len(mo_list)): 
    orbkit_output.display("\t%(m).6f\tMO %(n)s" 
	    % {'m':numpy.sum(mo_list[ii_mo]**2)*grid.d3r, 'n':mo_spec[ii_mo]['sym']})
  
  #--- Print the number of electrons ---
  orbkit_output.display("We have " + str(numpy.sum(rho)*grid.d3r) + " electrons.")
  
  return rho
  #--- rho_compute_no_slice ---

def delta_rho_compute(geo_spec, geo_info, ao_spec, mo_spec, drv=[0,1,2]):
  #--- FUNCTION rho_compute ----------------------------------------------------
  #--- Function for the computation of the density -----------------------------
  #-----------------------------------------------------------------------------
  #--- The 3D grid is divided into Slices, and the computational tasks ---------
  #--- are distributed to the worker processes ---------------------------------
  #-----------------------------------------------------------------------------
  
  #--- Specify the global variable containing all desired information needed ---
  #--- by the function slice_rho ---
  global Spec
  Spec = {'geo_spec': geo_spec, 'ao_spec': ao_spec, 
	      'mo_spec': mo_spec, 'Derivative': drv}
  
  #--- Save the complete initial grid ---
  x = grid.x
  y = grid.y
  z = grid.z
  N = copy.deepcopy(grid.N_)
  Min = copy.deepcopy(grid.min_)
  Max = copy.deepcopy(grid.max_)
  
  #--- Make slices ---
  rho = numpy.zeros((len(grid.x),len(grid.y),len(grid.z)))
  delta_rho = numpy.zeros((3,len(grid.x),len(grid.y),len(grid.z)))
  sDim = 0
  sNum = grid.N_[sDim]
  grid.N_[sDim] = 1
  
  #--- The number of worker processes is capped to the number of ---
  #--- grid points in x-direction. --- 
  if options.numproc > sNum: options.numproc = sNum
  
  #--- Print information regarding the density calculation ---
  orbkit_output.display("\nStarting the density calculation...")
  orbkit_output.display("The grid has been seperated in 2d slices and the")
  orbkit_output.display("calculation will be carried out with %d subprocesses." 
		  % options.numproc)
  orbkit_output.display("\nThere are %d contracted AOs and %d MOs to be calculated."
		  % (len(mo_spec[0]['coeffs']), len(mo_spec)))
  
  #--- Initialize some additional user information ---
  status_old = 0
  s_old = 0
  t = [time.time()]
  
  #--- Initialize an array to store the norm of the MOs ---
  mo_norm = numpy.zeros((len(mo_spec)))
    
  #--- Start the worker processes --
  pool = Pool(processes=options.numproc)
  
  #--- Write the slices in x to an array zz ---
  zz=[]
  for s in range(sNum):
    zz.append((numpy.array([x[s]])))
  
  #--- Compute the derivative of the density slice by slice ---

  it = pool.imap(slice_rho, zz)
  for s in range(sNum):
    rho[s,:,:], norm, d_rho = it.next()
    for ii_d in drv:
      delta_rho[ii_d,s,:,:] = d_rho[ii_d,:,:,:]
    mo_norm += norm
    
    #--- Print out the progress of the computation ---
    status = round(s*100/float(sNum))
    if not status % 20 and status != status_old: 
      t.append(time.time())
      orbkit_output.display('\tFinished %(f)d%% (%(s)d slices in %(t).3fs)'
		      % {'f': status,
			 's': s + 1 - s_old,
			 't': t[-1]-t[-2]})
      status_old = status
      s_old = s + 1
  
  #--- Close the worker processes --
  pool.close()
  
  #--- Restore the complete initial grid ---
  grid.N_ = N
  grid.min_ = Min
  grid.max_ = Max
  grid.grid_init()

  #--- Print the norm of the MOs ---
  orbkit_output.display('\nNorm of the MOs:')
  for ii_mo in range(len(mo_norm)):
    orbkit_output.display('\t%(m).6f\tMO %(n)s' 
		    % {'m':mo_norm[ii_mo]*grid.d3r, 'n':mo_spec[ii_mo]['sym']})
  
  #--- Print the number of electrons ---
  orbkit_output.display('We have %f electrons.' % integration(rho))
  
  return rho, delta_rho
  #--- delta_rho_compute ---

def delta_rho_compute_no_slice(geo_spec, geo_info, ao_spec, mo_spec):
  #--- FUNCTION rho_compute_no_slice -------------------------------------------
  #--- Function for the computation of the density without slicing the grid ----
  #-----------------------------------------------------------------------------
  
  #--- Initialize a numpy array for the derivative of the density ---
  ao_list = ao_creator(geo_spec,ao_spec)
  mo_list = mo_creator(ao_list,mo_spec)
  
  #--- Initialize a numpy array for the density
  rho = numpy.zeros((len(grid.x),len(grid.y),len(grid.z)))
  
  #--- Calculate the density and the norm
  for ii_mo in range(len(mo_list)): 
    rho += numpy.square(numpy.abs(mo_list[ii_mo])) * mo_spec[ii_mo]['occ_num']
  
  #--- Print the norm of the MOs ---
  orbkit_output.display('Norm of the MOs:')
  for ii_mo in range(len(mo_list)): 
    orbkit_output.display('\t%(m).6f\tMO %(n)s' % 
		    {'m':numpy.sum(mo_list[ii_mo]**2)*grid.d3r, 
		     'n':mo_spec[ii_mo]['sym']})
  
  #--- Print the number of electrons ---
  orbkit_output.display('We have %f electrons.' % integration(rho))
  
  #--- Print information ---
  orbkit_output.display('\nCalculating the derivative of the density...')
  delta_rho = numpy.zeros((3,len(grid.x),len(grid.y),len(grid.z)))
  
  #--- Loop over spatial directions ---
  string = ['x', 'y', 'z']
  for d_ii in range(3):
    orbkit_output.display('\t...with respect to %s' % string[d_ii])
    #--- Calculate the derivatives of the AOs and MOs ---
    delta_ao_list = delta_ao_creator(geo_spec,ao_spec,drv=d_ii)
    delta_mo_list = mo_creator(delta_ao_list,mo_spec)
    
    
    #--- Calculate the derivative of the density
    for ii_mo in range(len(delta_mo_list)): 
      delta_rho[d_ii,:,:,:] += (mo_spec[ii_mo]['occ_num'] * 
				2 * delta_mo_list[ii_mo]*mo_list[ii_mo])
  
  return rho, delta_rho
  #--- delta_rho_compute_no_slice ---


#####################################
#--- Support Code ---
#####################################
#--- Assign to every AO symbol (s,p,d,etc.) the quantum number l ---------------
allTheLetters = string.lowercase
orbit = 'spd%s' % allTheLetters[5:]
lquant = {}
for ii in range(len(orbit))[::-1]:
  lquant[orbit[ii]] = ii

#--- MOLPRO AO order ---
exp = []
exp.append([(0,0,0)])			# s orbitals

exp.append([(1,0,0), (0,1,0), (0,0,1)])	# p orbitals

exp.append([(2,0,0), (0,2,0), (0,0,2),
	    (1,1,0), (1,0,1), (0,1,1)])	# d orbitals
	    
exp.append([(3,0,0), (0,3,0), (0,0,3),
	    (1,2,0), (2,1,0), (2,0,1),
	    (1,0,2), (0,1,2), (0,2,1),
	    (1,1,1)])			# f orbitals
	    
exp.append([(4,0,0), (0,4,0), (0,0,4),
	    (3,1,0), (3,0,1), (1,3,0),
	    (0,3,1), (1,0,3), (0,1,3),
	    (2,2,0), (2,0,2), (0,2,2),
	    (2,1,1), (1,2,1), (1,1,2)])	# g orbitals

def l_deg(l=None,ao=None):
  #--- FUNCTION l_deg ----------------------------------------------------------
  #--- Calculate the degeneracy of a given atomic orbital ----------------------
  #-----------------------------------------------------------------------------
  #--- Works with the molpro out nomenclature for cartesian Harmonics: ---------
  #---      s->'s', p->['x','y','z'], d-> ['xx','yy', etc.], etc.      ---------
  #---      e.g., l_deg(ao='xxy')                                      ---------
  #-----------------------------------------------------------------------------
  #--- Works with quantum number l for the cartesian Harmonic:         ---------
  #---      e.g., l_deg(l=1)                                           ---------
  #-----------------------------------------------------------------------------  
  if ao != None:
    if ao == 's':
      return 1
    else:
      l = len(ao)
  return (l+1)*(l+2)/2
  #--- l_deg ---

def integration(matrix,x=None,y=None,z=None):
  if x == None: x = grid.x
  if y == None: y = grid.y
  if z == None: z = grid.z  
  
  if matrix.squeeze().ndim == 3:
    integral = integrate.simps(matrix, x, axis=0, even='avg')
    integral = integrate.simps(integral, y, axis=0, even='avg')
    integral = integrate.simps(integral, z, axis=0, even='avg')
  elif matrix.squeeze().ndim == 2:
    if len(x) == 1:
      r = y
      matrix = matrix[0,:,:]
    elif len(y) == 1:
      r = x
      matrix = matrix[:,0,:]
    else:
      print 'dim(z) = 1! No cylindrical coordinates...'
      return 255
    [Z,R] = numpy.meshgrid(z,r)
    integral = 2*numpy.pi*integrate.simps(R*matrix, r, axis=0, even='avg')
    integral = integrate.simps(integral, z, axis=0, even='avg')
  else: 
    #print 'Error in integration! ndim is not 3 or 2...'
    return numpy.sum(matrix)
  return integral

#--- C++ support code ---
class cSupportCode:
  math = '''
    #include <math.h>
    #define _USE_MATH_DEFINES
    '''
  norm = '''
    #include <math.h>
    #define _USE_MATH_DEFINES
    
    double ao_norm(int l,int m,int n,double *alpha);
    int doublefactorial(int n);
    
    double ao_norm(int l,int m,int n,double *alpha)
    {
      // --- FUNCTION ao_norm --- calculate normalization of uncontracted AOs --
      double norm;
      norm = pow((2./M_PI),(3./4.)) * (pow(2.,(l+m+n)) * 
	     pow(*alpha,((2.*l+2.*m+2.*n+3.)/4.))) / (pow((doublefactorial(2.*l-1.)*
	     doublefactorial(2.*m-1.) * doublefactorial(2.*n-1.)),0.5));
      return norm;
    }

    int doublefactorial(int n)
      // --- FUNCTION doublefactorial --- does what the name suggests ----------
    {
    if (n <= 0)
      {
	return 1;
      }
      else
      {
	return n * doublefactorial(n-2);
      }
    }
  '''
  
  xyz = '''
  
    double xyz(double x,double y,double z,int lx,int ly,int lz);
    
    double xyz(double x,double y,double z,int lx,int ly,int lz)
    {
      // --- FUNCTION xyz --- calculate x^lx * y^ly * z^lz ---------------------
      double Xl, Yl, Zl;
      
      if (lx == 0)
	Xl = 1;
      else if (lx == 1)
	Xl = x;
      else
	Xl = pow(x, lx);
	
      if (ly == 0)
	Yl = 1;
      else if (ly == 1)
	Yl = y;
      else
	Yl = pow(y, ly);
	
      if (lz == 0)
	Zl = 1;
      else if (lz == 1)
	Zl = z;
      else
	Zl = pow(z, lz);

      return Xl * Yl * Zl;
    }

  '''
