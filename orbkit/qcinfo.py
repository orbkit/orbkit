# -*- coding: iso-8859-1 -*-
'''Module for processing the data read from the output files of quantum chemical 
software. '''
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
#from scipy.constants import value as physical_constants
import numpy
from os import path

#: Contains the mass conversion factor to atomic units
u_to_me = 1822.88839
nist_mass = None
#: Standard atomic masses as "Linearized ASCII Output", see http://physics.nist.gov
nist_file = path.join(path.dirname(path.realpath(__file__)),
                      'supporting_data/Atomic_Weights_NIST.html')
# see http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii2&isotype=some

def read_nist():
  global nist_mass
  
  f = open(nist_file,'r')
  flines = f.readlines()
  f.close()
  
  nist_mass = []
  index = None
  new = True
  
  def rm_brackets(text,rm=['(',')','[',']']):
    for i in rm:
      text = text.replace(i,'')
    return text
  
  for line in flines:
    thisline = line.split()
    if 'Atomic Number =' in line:
      i = int(thisline[-1]) - 1
      new = (i != index)
      if new:
        nist_mass.append(['',0])
      index = i
    elif 'Atomic Symbol =' in line and new:
      nist_mass[index][0] = thisline[-1]
    elif 'Standard Atomic Weight =' in line and new:
      nist_mass[index][1] = float(rm_brackets(thisline[-1]))

def standard_mass(atom):
  if nist_mass is None:
    read_nist()
  
  try:
    atom = int(atom) - 1
    return nist_mass[atom][1] * u_to_me
  except ValueError:
    return dict(nist_mass)[atom.title()] * u_to_me

class QCinfo:
  def __init__(self):
    self.geo_spec = []
    self.geo_info = []
    self.ao_spec  = []
    self.mo_spec  = []
    self.etot     = 0.
    self.com      = None
    self.coc      = None
    # transition dipole information
    self.tdm_states      = {'multiplicity' : None,
                            'energy'       : None,
                            'dipole'       : None}
    self.tdm_transitions = {'dipole'       : None}
#  self.mo_coeff = None
#  self.mo_occup = None
#  self.mo_energ = None
#  self.mo_sym   = None

  def get_com(self,nuc_list=[]):
    self.com   = numpy.zeros(3)
    total_mass = 0.
    if not nuc_list:
      nuc_list = range(len(self.geo_spec)) # iterate over all nuclei
    for ii in nuc_list:
      nuc_mass    = standard_mass(self.geo_info[ii][0])
      self.com   += numpy.multiply(self.geo_spec[ii],nuc_mass)
      total_mass += nuc_mass
    self.com = self.com/total_mass
    return self.com

  def get_coc(self):
    if self.coc is None: self.coc = 0
    return self.coc

