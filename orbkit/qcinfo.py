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
from copy import copy

from orbkit.display import display
from .units import u2me, aa2a0
from .tools import get_lxlylz
from orbkit.read.tools import get_atom_symbol, standard_mass

class QCinfo:
  '''Class managing all information from the from the output 
  files of quantum chemical software.
  
  See :ref:`Central Variables` in the manual for details.
  '''
  def __init__(self, filename=None):
    self.geo_info = []
    self.geo_spec = []
    self.ao_spec  = []
    self.ao_spherical = None
    self.mo_spec  = []
    self.etot     = 0.
    self.com      = 'Center of mass can be calculated with self.get_com().'
    self.coc      = 'Center of charge can be calculated with self.get_coc().'
    self.bc      =  'Barycenter of scalar field can be calculated with self.get_bc().'
    self.pop_ana  = {}
    # transition dipole information
    self.states         = {'multiplicity' : None,
                           'energy'       : None}
    self.dipole_moments = None

    if filename:
      self.read(filename)

  def save(self, filename=None):
    import hashlib, time
    ao_lxlylz, cont2atoms, pg_expcont, cont2prim, mo_occ, mo_eig, mo_sym, mo_coeff = self.old2new_formats()
    date = time.strftime("%d/%m/%Y") 
    time = time.strftime("%H:%M:%S")
    date_time = date + time
    fingerprint = hashlib.md5(filename + date_time).hexdigest()
    if not filename:
      filename = 'default_output'
    if '.npz' not in filename:
      filename += '.npz'

    numpy.savez_compressed(filename,
                           fingerprint=fingerprint,
                           date=date,
                           time=time,
                           cont2atoms=cont2atoms,
                           cont2prim=cont2prim,
                           pg_expcont=pg_expcont,
                           ao_lxlylz=ao_lxlylz,
                           geo_spec=self.geo_spec,
                           geo_info=self.geo_info,
                           ao_spherical=self.ao_spherical,
                           mo_occ=mo_occ,
                           mo_eig=mo_eig,
                           mo_sym=mo_sym,
                           mo_coeff=mo_coeff)
    return

  def read(self, filename):
    if isinstance(filename, str):
      if '.npz' not in filename:
        filename += '.npz'
      fname = open(filename, 'r')

    data = numpy.load(fname)
    cont2atoms = data['cont2atoms']
    cont2prim = data['cont2prim']
    pg_expcont = data['pg_expcont']
    ao_lxlylz = data['ao_lxlylz']
    mo_occ = data['mo_occ']
    mo_eig =  data['mo_eig']
    mo_sym = data['mo_sym']
    mo_coeff = data['mo_coeff']
    self.geo_spec = data['geo_spec']
    self.geo_info = data['geo_info']
    self.ao_spherical = None

    if data['ao_spherical']:
      ao_spherical = list(map(tuple, numpy.reshape(data['ao_spherical'], (2, -1))))

    self.new2old_formats(ao_lxlylz, cont2atoms, pg_expcont, cont2prim, mo_occ, mo_eig, mo_sym, mo_coeff)
    display('Loaded QCInfo class from file {0}'.format(filename))
    display('File was created on {0} at {1}'.format(data['date'], data['time']))

    return

  def new2old_formats(self, ao_lxlylz, cont2atoms, pg_expcont, cont2prim, mo_occ, mo_eig, mo_sym, mo_coeff):
    self.ao_spec = [{'atom': None,
                       'coeffs': None,
                       'pnum': None,
                       'exp_list': None}
                        for i in range(len(cont2prim))]

    for ic in range(len(cont2prim)):
      self.ao_spec[ic]['atom'] = cont2atoms[ic]
      self.ao_spec[ic]['coeffs'] = pg_expcont[cont2prim[ic]]
      self.ao_spec[ic]['pnum'] = len(cont2prim[ic])
      self.ao_spec[ic]['exp_list'] = list(map(tuple, numpy.reshape(ao_lxlylz[ic] , (3, -1))))

    self.mo_spec = [{'coeffs': None,
                     'energy': None,
                     'occ_num': None,
                     'sym': None}
                     for i in range(len(mo_occ))]

    for imo in range(len(mo_occ)):
      self.mo_spec[imo]['coeffs'] = mo_coeff[imo]
      self.mo_spec[imo]['energy'] = mo_eig[imo]
      self.mo_spec[imo]['sym'] = mo_sym[imo]
      self.mo_spec[imo]['occ_num'] = mo_occ[imo]
    return

  def copy(self):
    from copy import deepcopy
    qcinfo = deepcopy(self)
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
      self.geo_spec *= aa2a0
  
  def sort_mo_sym(self):
    '''Sorts mo_spec by symmetry.
    '''
    keys = []
    for i_mo in self.mo_spec:
      keys.append(i_mo['sym'].split('.'))
    keys = numpy.array(keys,dtype=int)
    self.mo_spec = list(numpy.array(self.mo_spec)[numpy.lexsort(keys.T)])

  def old2new_formats(self):
    ao_lxlylz = get_lxlylz(self.ao_spec)

    cont2atoms = numpy.zeros(len(self.ao_spec), dtype=numpy.intc)
    pg_expcont = []
    for ic, contracted in enumerate(self.ao_spec):
      cont2atoms[ic] = contracted['atom']
      for prim in contracted['coeffs']:
        pg_expcont.append(prim)

    pg_expcont = numpy.array(pg_expcont, dtype=numpy.float64)

    cont2prim = numpy.zeros((len(self.ao_spec), len(pg_expcont)), dtype=bool)
    i_prim = 0
    for i_cont, contracted in enumerate(self.ao_spec):
      max_i = i_prim + len(contracted['coeffs'])
      cont2prim[i_cont, i_prim:max_i] = True
      i_prim += len(contracted['coeffs'])

    mo_occ = numpy.zeros(len(self.mo_spec), dtype=numpy.float64)
    mo_eig = numpy.zeros(len(self.mo_spec), dtype=numpy.float64)
    mo_sym = []
    mo_coeff = numpy.zeros((len(self.mo_spec), len(self.mo_spec[0]['coeffs'])), dtype=numpy.float64)

    for i_mo, mo in enumerate(self.mo_spec):
      mo_occ[i_mo] = mo['occ_num']
      mo_eig[i_mo] = mo['energy']
      mo_sym.append(mo['sym'])
      for i_cont, coeff in enumerate(mo['coeffs']):
        mo_coeff[i_mo,i_cont] = coeff
    mo_sym = numpy.array(mo_sym)
    return ao_lxlylz, cont2atoms, pg_expcont, cont2prim, mo_occ, mo_eig, mo_sym, mo_coeff 

  def get_mo_labels(self):
    return ['MO %(sym)s, Occ=%(occ_num).2f, E=%(energy)+.4f E_h' % 
                  i for i in self.mo_spec]
  
  def get_mo_energies(self):
    return copy(self.mo_eig)
  
  def get_mo_occ(self):
    return copy(self.mo_occ)
  
  def get_nmoocc(self):
    return sum(self.get_mo_occ())
  
  def get_com(self,nuc_list=None):
    '''Computes the center of mass.
    '''
    self.com   = numpy.zeros(3)
    total_mass = 0.
    if nuc_list is None:
      nuc_list = list(range(len(self.geo_spec))) # iterate over all nuclei
    for ii in nuc_list:
      nuc_mass    = standard_mass(self.geo_info[ii][0])
      self.com   += numpy.multiply(self.geo_spec[ii],nuc_mass)
      total_mass += nuc_mass
    self.com = self.com/total_mass
    return self.com

  def get_coc(self):
    '''Computes the center of charge.
    '''
    self.coc     = numpy.zeros(3)
    total_charge = 0.
    for ii in range(len(self.geo_info)):
      nuc_charge    = float(self.geo_info[ii][2])
      self.coc     += numpy.multiply(self.geo_spec[ii],nuc_charge)
      total_charge += nuc_charge
    self.coc = self.coc/total_charge
    return self.coc
    
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
  
  def todict(self):
    '''Converts all essential variables into a dictionary.
    '''
    dct = {}
    keys = ['geo_spec',
            'geo_info',
            'ao_spec',
            'ao_spherical',
            'mo_spec']
    for key in keys:
      dct[key] = getattr(self,key)
    return dct
  
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
    

class CIinfo:
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
      self.coeffs = self.coeffs[i]
      self.occ = self.occ[i]
  def copy(self):
    ciinfo = self.__class__(method=self.method)
    if self.coeffs != []:
      ciinfo.coeffs = numpy.copy(self.coeffs)
    if self.occ != []:
      ciinfo.occ = numpy.copy(self.occ)
    if self.info is not None:
      ciinfo.info = self.info.copy()    
    if self.moocc is not None:
      ciinfo.moocc = self.moocc.copy()    
    return ciinfo
  def todict(self):
    return self.__dict__
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

