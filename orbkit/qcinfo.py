# -*- coding: iso-8859-1 -*-
from scipy.constants import value as physical_constants
import numpy as np

# standard atomic masses from http://physics.nist.gov
standard_mass = {
  '1' : 1.00794     ,  'H' : 1.00794    ,
  '2' : 4.002602    , 'He' : 4.002602    , 'HE' : 4.002602    ,
  '3' : 6.941       , 'Li' : 6.941       , 'LI' : 6.941       ,
  '4' : 9.012182    , 'Be' : 9.012182    ,
  '5' : 10.811      ,  'B' : 10.811      ,
  '6' : 12.0107     ,  'C' : 12.0107     ,
  '7' : 14.0067     ,  'N' : 14.0067     ,
  '8' : 15.9994     ,  'O' : 15.9994     ,
  '9' : 18.9984032  ,  'F' : 18.9984032  ,
 '10' : 20.1797     , 'Ne' : 20.1797     , 'NE' : 20.1797     ,
 '11' : 22.98976928 , 'Na' : 22.98976928 , 'NA' : 22.98976928 ,
 '12' : 24.3050     , 'Mg' : 24.3050     , 'MG' : 24.3050     ,
 '13' : 26.9815386  , 'Al' : 26.9815386  , 'AL' : 26.9815386  ,
 '14' : 28.0855     , 'Si' : 28.0855     , 'SI' : 28.0855     ,
 '15' : 30.973762   ,  'P' : 30.973762   ,
 '16' : 32.065      ,  'S' : 32.065      ,
 '17' : 35.453      , 'Cl' : 35.453      , 'CL' : 35.453      ,
 '18' : 39.948      , 'Ar' : 39.948      , 'AR' : 39.948      ,
 '19' : 39.0983     ,  'K' : 39.0983     , 
 '20' : 40.078      , 'Ca' : 40.078      , 'CA' : 40.078      ,
 '21' : 44.955912   , 'Sc' : 44.955912   , 'SC' : 44.955912   ,
 '22' : 47.867      , 'Ti' : 47.867      , 'TI' : 47.867      ,
 '23' : 50.9415     ,  'V' : 50.9415     ,
 '24' : 51.9961     , 'Cr' : 51.9961     , 'CR' : 51.9961     ,
 '25' : 54.938045   , 'Mn' : 54.938045   , 'MN' : 54.938045   ,
 '79' : 196.966569  , 'Au' : 196.966569  , 'AU' : 196.966569
 }
# convert to atomic units
u_to_me      = 1822.88839
dict_entries = [elements for elements, masses in standard_mass.items()]
for dict_entry in dict_entries: standard_mass[dict_entry] *= u_to_me

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
    self.com   = np.zeros(3)
    total_mass = 0.
    if not nuc_list:
      # iterate over all nuclei
      for ii in range(len(self.geo_spec)):
        nuc_mass    = standard_mass[self.geo_info[ii][0]]
        self.com   += np.multiply(self.geo_spec[ii],nuc_mass)
        total_mass += nuc_mass
    else:
      # iterate over given nuclei
      for ii in nuc_list:
        nuc_mass    = standard_mass[self.geo_info[ii][0]]
        self.com   += np.multiply(self.geo_spec[ii],nuc_mass)
        total_mass += nuc_mass
    self.com = self.com/total_mass
    return self.com

  def get_coc(self):
    if self.coc is None: self.coc = 0
    return self.coc

