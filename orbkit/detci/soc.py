from __future__ import print_function
'''Interface to MolSOC code by Sandro Giuseppe Chiodo and Monica Leopoldini

The code can be obtained here

http://cpc.cs.qub.ac.uk/summaries/AERK_v1_0.html

If you use the code the required reference is

Comput. Phys. Commun. 185(2014)676

https://doi.org/10.1016/j.cpc.2013.10.014
'''

import os
import numpy

class SOC:
  '''Class used to calculate SO-coupling between singlet and tripplet states.

  The implementation is based on PySOC by Gao et al.

  J. Chem. Theory Comput., 2017, 13 (2), pp 515-524

  https://doi.org/10.1021/acs.jctc.6b00915

  Currently only singlet and tripplet states are supported.
  '''
  def __init__(self, qc_bra, qc_ket, ci_bra, ci_ket):
    self.qc_bra = qc_bra
    self.qc_ket = qc_ket
    self.ci_bra = ci_bra
    self.ci_ket = ci_ket

    if not self.qc_ket:
      self.qc_ket = self.qc_bra

    if self.ci_bra.info['spin'] not in ['Singlet', 'Triplet']\
    or self.ci_ket.info['spin'] not in ['Singlet', 'Triplet']:
      raise ValueError('Bra and Ket states must be either singlets or triplets!')

    if not self.qc_bra.charge == self.qc_ket.charge:
      raise ValueError('Bra and Ket states must have the same total charge!')

    if not self.qc_bra.ao_spec == self.qc_ket.ao_spec:
      raise ValueError('Bra and Ket states must have the same basis set!')

    if not self.qc_bra.comp_geo_info(self.qc_ket.geo_info) or not numpy.allclose(self.qc_bra.geo_spec, self.qc_ket.geo_spec):
      raise ValueError('Bra and Ket states must have the same atomic positions!')


  def write_molsoc_inp(self, path='.'):
    '''Writes molsoc.inp file.
    
    **Parameters:**
    
      path : str
        Specifies where file are written.
    '''

    multiplicities = {'Singlet': 1, 'Triplet': 3}

    charge = self.qc_bra.charge
    m1 = multiplicities[self.ci_bra.info['spin']]
    m2 = multiplicities[self.ci_ket.info['spin']]

    #Use one-electron SO coupling calculations with screened-nuclear charge
    #J. Comput. Chem., 30 (2009), p. 83
    #https://doi.org/10.1002/jcc.20982
    #Positions given in Angstrom
    #Calculate dipole moments
    header = '''#input for MolSOC
ANG  ZEF DIP CAR
#\n'''
    
    geom = ''

    info = self.qc_ket.geo_info
    pos = self.qc_ket.geo_spec
    for ia in range(len(self.qc_ket.geo_spec)):
      geom += str(info[0,0]) + ' '
      for ip in range(3):
        geom += str(pos[ia,0]) + ' '
      geom += str(info[ia,-1])
      geom += '\n'

    geom += 'End'
    out = header + geom

    with open(os.path.join(path, 'molsoc.inp'), 'wb') as fd:
      print(out, file=fd)


  def write_molsoc_basis(self, path='.'):
    '''Writes molsoc_basis file.
    
    **Parameters:**
    
      path : str
        Specifies where file are written.
    '''

    names = {0: 'S', 1: 'P', 2: 'D', 3: 'F', 4: 'G'}
    out = ''
    for ia, atom in enumerate(self.qc_bra.geo_info):
      out += atom[0] + ' ' + str(0) + '\n'
      for ic in numpy.argwhere(self.qc_bra.ao_spec.get_cont2atoms() == ia):
        if self.qc_bra.ao_spec[ic]['type']:
          out += self.qc_bra.ao_spec[ic]['type'].upper()
        else:
          out += names[self.qc_bra.ao_spec[ic]['ao_spherical'][0][0]]
        out += ' ' + str(len(self.qc_bra.ao_spec[ic])) + ' '+ str(1.0) + '\n'
        for prim in self.qc_bra.ao_spec[ic]['coeffs']:
          out += str(prim[0]) + ' ' + str(prim[1]) + '\n'
      if ia < len(self.qc_bra.geo_info)-1:
        out += ' ****\n'

    out += 'END'

    with open(os.path.join(path, 'molsoc_basis'), 'wb') as fd:
      print(out, file=fd)



  def write_mos(self, path='.'):
    '''Writes mos1 and mos2 files in Turbomole format.
    
    **Parameters:**
    
      path : str
        Specifies where file are written.
    '''

    qc_data = {0: self.qc_bra, 1: self.qc_ket}
    ci_data = {0: self.ci_bra, 1: self.ci_ket}
    spins = {0: 'Alpha', 1: 'Beta'}
#    for fname in range(2):
#      for spin in 











