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
  def __init__(self, qc_bra, qc_ket, ci_bra, ci_ket):
    self.qc_bra = qc_bra
    self.qc_ket = qc_ket
    self.ci_bra = ci_bra
    self.ci_ket = ci_ket

    if not self.qc_ket:
      self.qc_ket = self.qc_bra

    if not self.qc_bra.charge == self.qc_ket.charge:
      raise ValueError('Bra and Ket states must have the same total charge!')

    if not self.qc_bra.ao_spec == self.qc_ket.ao_spec:
      raise ValueError('Bra and Ket states must have the same basis set!')

    if not self.qc_bra.comp_geo_info(self.qc_ket.geo_info) or not numpy.allclose(self.qc_bra.geo_spec, self.qc_ket.geo_spec):
      raise ValueError('Bra and Ket states must have the same atomic positions!')


  def write_molsoc_inp(self, path='.', **args):
    '''Writes molsoc.inp file.
    
    **Parameters:**
    
      path : str
        Specifies where file are written.

    **Note:**
  
      All arguments supported by MolSOC can be passed as optinal extra arguments.
      Note that most parameters are determied automatically.
      Important options include:
      
        - ZE: Screened-nuclear charge method is used to calculate SO coupling.
        - Dip/Qua/Oct: Dipole, quadrupole, and octupole moment components are calculated.
    '''

    header = '''#input for MolSOC
Tur Ang'''

    assert self.qc_bra.ao_spec.spherical == self.qc_ket.ao_spec.spherical

    if self.qc_bra.ao_spec.spherical:
      header += 'Sph'
    else:
      header += 'Car'

    assert self.qc_bra.mo_spec.spinpolarized == self.qc_ket.mo_spec.spinpolarized

    if self.qc_bra.mo_spec.spinpolarized:
      header += 'UKS'
    else:
      header += 'RKS'

    for arg in **args:
      header += ' ' + arg

    header += '\n'

    multiplicities = {'Singlet': 1, 'Doublet': 2, 'Triplet': 3, 'Quartet': 4}
    m1 = multiplicities[self.ci_bra.info['spin']]
    m2 = multiplicities[self.ci_ket.info['spin']]
    
    header += str(self.qc_bra.charge) + '  ' + str(m1)

    

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


  def unnormalize_orbitals(self, qc):
    '''Apply the inverse of of the normalization procedure used to
       normalize molecular orbitals during reading.
    
    **Parameters:**
    
      qc : QCinfo instance

    **Returns:**
    
      qc : QCinfo instance
    '''

    # Convert MO coefficients
    def dfact(n):
      if n <= 0:
        return 1
      else:
        return n * dfact(n-2)

    mo = qc.mo_spec.get_coeff()
    for i,j in enumerate(qc.ao_spec.get_lxlylz()):
      norm = (dfact(2*j[0] - 1) * dfact(2*j[1] - 1) * dfact(2*j[2] - 1))
      j = sum(j)
      if j >1: 
        mo[:,i] /= numpy.sqrt(norm)   
    for ii in range(len(qc.mo_spec)):
      qc.mo_spec[ii]['coeffs'] = mo[ii]
  
    return qc

  def write_mo(self, eigen, coeff):
    tmp = ''
    for imo in range(len(eigen)):
      tmp += '     {0}  a      eigenvalue={1}   nsaos={2}\n'.format(imo+1,eigen[imo],len(coeff[imo]))
      iao = 0
      while iao <= len(coeff[imo]):
        if iao < len(coeff[imo]) - 4:
          for c in coeff[imo,iao:iao+4]:
            tmp += str(c) + '  '
        else:
          for c in coeff[imo,iao:]:
            tmp += str(c) + '  '
        iao += 4
        tmp += '\n'
    return tmp

  def write_molecular_orbitals(self, path='.'):
    '''Writes mos1 and mos2 files in Turbomole format.
    
    **Parameters:**
    
      path : str
        Specifies where file are written.
    '''

    qc_data = {0: self.qc_bra, 1: self.qc_ket}
    #This is not normal Turbomole behavior but its what MolSOC wants
    spin_seperators = {0: 'alpha', 1: 'beta '}
    files = {0: 'mos1', 1: 'mos2'}
    for iqc in range(2):

      header = '''$scfmo    scfconv=None   format(4d20.14)                                           
# SCF total energy is      None
#
'''
      out = header
      qc = self.unnormalize_orbitals(qc_data[iqc])
      eigen = qc_data[iqc].mo_spec.get_eig()
      coeff = qc_data[iqc].mo_spec.get_coeff()
      if qc.mo_spec.spinpolarized:
        for spin in ['alpha', 'beta']:
          out += spin + '\n'
          spin_index = qc.mo_spec.get_spin(spin)
          eigen = eigen[spin_index]
          coeff = coeff[spin_index]
          out += write_mo(eigen, coeff)
      else:
        out += write_mo(eigen, coeff)

      out += '$end'
      with open(os.path.join(path, files[iqc]), 'wb') as fd:
        print(out, file=fd)
      











