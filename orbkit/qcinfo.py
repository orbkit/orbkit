# -*- coding: iso-8859-1 -*-
'''Module for processing the data read from the output files of quantum chemical
software. '''
'''
orbkit
Gunter Hermann, Vincent Pohl, Lukas Eugen Marsoner Steinkasserer, Axel Schild, and Jean Christophe Tremblay

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
from copy import copy,deepcopy

from .display import display
from .units import u_to_me, aa_to_a0
from .tools import get_atom_symbol, standard_mass
from .orbitals import AOClass, MOClass

class QCinfo:
  '''Class managing all information from the from the output
  files of quantum chemical software.

  See :ref:`Central Variables` in the manual for details.
  '''
  def __init__(self, data=None):
    self.geo_info = []
    self.geo_spec = []
    
    if data:
      self.geo_spec = data['geo_spec']
      self.geo_info = data['geo_info']
      self.format_geo()
      if isinstance(data['ao_spec'], numpy.ndarray):
        ao_spec = data['ao_spec'][numpy.newaxis][0]
      else:
        ao_spec = data['ao_spec']
      if isinstance(data['mo_spec'], numpy.ndarray):
        mo_spec = data['mo_spec'][numpy.newaxis][0]
      else:
        mo_spec = data['mo_spec']
      self.ao_spec = AOClass(restart=ao_spec)
      self.mo_spec = MOClass(restart=mo_spec)
    else:
      self.geo_info = []
      self.geo_spec = []
      self.ao_spec = AOClass()
      self.mo_spec = MOClass()

  def __eq__(self, other):
    if not isinstance(other, QCinfo):
      raise TypeError('Comaring of QCinfo to non QCinfo object not defined')
    same = [
            self.comp_geo_info(other.geo_info),
            numpy.allclose(self.geo_spec, other.geo_spec),
            self.ao_spec == other.ao_spec,
            self.mo_spec == other.mo_spec
            ]
    return all(same)

  def update(self):
    self.ao_spec.update()
    self.mo_spec.update()

  def comp_geo_info(self, geo2):
    same = True
    for atom1, atom2 in zip(self.geo_info, geo2):
      if not len(atom1) == len(atom2):
        raise ValueError('Atom object are of different length!')
      for i in range(len(atom1)):
        if atom1[i] != atom2[i]:
          same = False
    return same

  def copy(self):
    from copy import deepcopy
    qcinfo = deepcopy(self)
    qcinfo.update()
    return qcinfo

  def format_geo(self, is_angstrom=False):
    '''Converts geo_info and geo_spec to a universal format.
    **Parameters:**

    is_angstrom : bool, optional
      If True, input is assumed to be in Angstrom and positions are converted to Bohr radii.
    '''
    for i in self.geo_info:
      i[0] = get_atom_symbol(i[0])
      i[2] = float(i[-1])
    self.geo_info = numpy.array(self.geo_info)
    self.geo_spec = numpy.array(self.geo_spec,dtype=float)
    if is_angstrom:
      self.geo_spec *= aa_to_a0

  def get_com(self,nuc_list=None):
    '''Computes the center of mass.
    '''
    com   = numpy.zeros(3)
    total_mass = 0.
    if nuc_list is None:
      nuc_list = list(range(len(self.geo_spec))) # iterate over all nuclei
    for ii in nuc_list:
      nuc_mass    = standard_mass(self.geo_info[ii][0])
      com   += numpy.multiply(self.geo_spec[ii],nuc_mass)
      total_mass += nuc_mass
    com = com/total_mass
    return com


  def get_charge(self, nuclear=True, electron=True):
    '''Computes total charge of the system.
    '''
    charge = 0.
    if electron:
      charge -= sum(self.mo_spec.get_occ())
    if nuclear:
      for ii in range(len(self.geo_info)):
        charge += float(self.geo_info[ii][2])
    return charge

  def get_coc(self):
    '''Computes the center of charge.
    '''
    coc     = numpy.zeros(3)
    for ii in range(len(self.geo_info)):
      nuc_charge    = float(self.geo_info[ii][2])
      coc += numpy.multiply(self.geo_spec[ii],nuc_charge)
    coc = coc / self.get_charge(nuclear=True)
    return coc

  def get_bc(self,matrix=None,is_vector=False):
    '''Calculates Barycenter for scalar field
    '''
    # Initialize variable
    self.bc = numpy.zeros(3)
    # Calculation of barycenter
    from orbkit import grid
    if not is_vector:
      grid.grid2vector()
    xyz = grid.tolist()
    for i in range(3):
      self.bc[i] = (matrix.reshape((-1,))*xyz[i]).sum()
    self.bc /= matrix.sum()
    if not is_vector:
      grid.vector2grid(*grid.N_)
    return self.bc

  def select_spin(self,restricted,spin=None):
    '''For an unrestricted calculation, the name of the MO
    ('sym' keyword in ``qc.mo_spec``) is modified, e.g., 3.1_b for MO 3.1 with
    beta spin and 3.1_a for MO 3.1 with alpha spin.
    For restricted calculation, the 'spin' keyword from ``qc.mo_spec`` is
    removed.

    **Parameters:**

    restricted : bool
      If True, removes the 'spin' keyword from ``qc.mo_spec``.
    spin : {None, 'alpha', or 'beta'}, optional
      If not None, returns exclusively 'alpha' or 'beta' molecular orbitals.
    '''
    # Only molecular orbitals of one spin requested?
    if spin is not None:
      for i in range(len(self.mo_spec))[::-1]:
        if self.mo_spec[i]['spin'] != spin:
          del self.mo_spec[i]

    if restricted:
      # Closed shell calculation
      for mo in self.mo_spec:
        del mo['spin']
    else:
      # Rename MOs according to spin
      for mo in self.mo_spec:
        mo['sym'] += '_%s' % mo['spin'][0]
      if not isinstance(self.mo_spec, MOClass):
        self.mo_spec = MOClass(self.mo_spec)
        self.mo_spec.get_spinstate()

  def todict(self):
    '''Returns the dictionary that is used to save QCinfo instance
    '''
    data = {}
    data['ao_spec'] = self.ao_spec.todict()
    data['mo_spec'] = self.mo_spec.todict()
    data['geo_spec'] = self.geo_spec
    data['geo_info'] = self.geo_info
    data['parent_class_name'] = self.__module__ + '.' + self.__class__.__name__
    return data

  def get_ase_atoms(self,bbox=None,**kwargs):
    '''Create an ASE atoms object.
    (cf. https://wiki.fysik.dtu.dk/ase/ase/atoms.html )

    **Parameters:**

    bbox : list of floats (bbox=[xmin,xmax,ymin,ymax,zmin,zmax]), optional
      If not None, sets the unit cell to the grid boundaries and moves the
      molecule in its center.

    **Returns:**

    atoms : Atoms object
      See https://wiki.fysik.dtu.dk/ase/ase/atoms.html for details

    .. Note::

      ASE has to be in the PYTHONPATH
    '''
    from ase import Atoms
    from ase.units import Bohr

    atoms = Atoms("".join(self.geo_info[:,0]),
                  positions=self.geo_spec*Bohr,
                  **kwargs)
    if bbox is not None:
      if len(bbox) != 6:
        raise ValueError("bbox has to have 6 elements")
      bbox = numpy.array(bbox)
      atoms.translate(-bbox[::2]*Bohr)
      atoms.cell = numpy.eye(3) * (bbox[1::2] - bbox[::2])*Bohr

    return atoms
  # Synonym
  atoms = get_ase_atoms

  def view(self,select=slice(None,None,None),bbox=None,**kwargs):
    '''Opens ase-gui with the atoms of the QCinfo class.
    (cf. https://wiki.fysik.dtu.dk/ase/ase/visualize/visualize.html )

    **Parameters:**

    select : slice or (array of int), default: all atoms
      Specifies the atoms to be shown.
    bbox : list of floats (bbox=[xmin,xmax,ymin,ymax,zmin,zmax]), optional
      If not None, sets the unit cell to the grid boundaries and moves the
      molecule in its center.

    .. Note::

      ASE has to be in the PYTHONPATH
    '''
    from ase import visualize
    visualize.view(self.get_ase_atoms(bbox=bbox,**kwargs)[select])

  @property
  def nuclear_repulsion(self):
    '''Calculates nuclear repulsion energy.'''
    from numpy.linalg import norm
    Vnn = 0
    Natoms = self.geo_spec.shape[0]
    for a in range(Natoms):
      Za = float(self.geo_info[a,2])
      Ra = self.geo_spec[a,:].astype(float)
      for b in range(a+1, Natoms):
        Zb = float(self.geo_info[b,2])
        Rb = self.geo_spec[b,:].astype(float)
        Vnn += Za*Zb / norm(Ra-Rb)
    return Vnn

class CIinfo():
  '''Class managing all information from the from the output
  files of quantum chemical software for CI calculations.

  The CI related features are in ongoing development.
  '''
  def __init__(self,method='ci'):
    self.method = method
    self.info   = None
    self.coeffs = []
    self.occ    = []
    self.moocc  = None

  def __str__(self):
    string = '%s' % self.method.upper()
    if self.info is not None:
      string += ' State %(state)s' % self.info
      if 'spin' in self.info.keys() and self.info['spin'] != 'Unknown':
         string += ' (%(spin)s)' % self.info
    if numpy.shape(self.coeffs) != (0,):
      string += ':\tNorm = %0.8f (%d Coefficients)' %(self.get_norm(),
                                                      len(self.coeffs))
    return  string
  def __eq__(self, other):
    try:
      return self.__dict__ == other.__dict__
    except ValueError:
      return False
  def get_norm(self):
    return sum(self.coeffs**2)
  def renormalize(self):
    self.coeffs /= self.get_norm()
  def apply_threshold(self,threshold,keep_length=False):
    i = numpy.abs(self.coeffs) > threshold
    if keep_length:
      self.coeffs[numpy.invert(i)] = 0.0
    else:
      if self.info['state'] == '0':
       self.coeffs = self.coeffs[i]
       self.occ = self.occ[:1]
      else:
       self.coeffs = self.coeffs[i]
       self.occ = self.occ[i]
      
  def copy(self):
    ciinfo = deepcopy(self)
    return ciinfo

  def get_moocc(self):
    if self.moocc is None:
      raise ValueError('ci.set_moocc(qc) has to be called first! (ci.moocc is not initialized)')
    return self.moocc

  def set_moocc(self,moocc):
    assert (moocc.dtype == numpy.intc), 'moocc has to be numpy.intc'
    self.moocc = moocc

  def hdf5_save(self,fid='out.h5',group='/ci:0',mode='w'):
    from orbkit.output import hdf5_open,hdf5_append
    from copy import copy
    for hdf5_file in hdf5_open(fid,mode=mode):
      dct = copy(self.todict())
      dct['info'] = numpy.array(dct['info'].items(),dtype=str)
      hdf5_append(dct,hdf5_file,name=group)

  def hdf5_read(self,fid='out.h5',group='/ci:0'):
    from orbkit.output import hdf5_open,hdf52dict
    for hdf5_file in hdf5_open(fid,mode='r'):
      for key in self.__dict__.keys():
        try:
          self.__dict__[key] = hdf52dict('%s/%s' % (group,key),hdf5_file)
        except KeyError:
          self.__dict__[key] = hdf5_file['%s' % group].attrs[key]
      self.__dict__['info'] = dict(self.__dict__['info'])

