# -*- coding: iso-8859-1 -*-
from scipy.constants import value as physical_constants

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

  def get_com(self):
    if self.com is None: self.com = 0
    return self.com

  def get_coc(self):
    if self.coc is None: self.coc = 0
    return self.coc
