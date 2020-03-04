import numpy
from os import path
from copy import copy
import sys
import warnings

try:
  from UserList import UserList
except ImportError:
  from collections import UserList

from .tools import *
from .display import display

ao_dict_synonyms = {'atom':             'atom',
                    'pnum':             'pnum',
                    'type':             'type',
                    'coeffs':           'coeffs',
                    'lxlylz':           'lxlylz',
                    'exp_list':         'lxlylz',
                    'lm':               'lm',
                    'ao_spherical':     'lm',
                    }

class AOClass(UserList):
  '''AO base class which contains all information on atomic orbitals.
  Two types of dataformats are available:

  1. Numpy-style data:

      cont2atoms : numpy.ndarray, dtype=numpy.intc, shape = (NAO)
        Transformation matrix between contracted GTO's and atoms.
      _assign_prim_to_cont : numpy.ndarray, dtype=numpy.intc, shape = (NPAO)
        Transformation matrix between contracted GTO's and primitive GTO's.
      pg_expcont : numpy.ndarray, dtype=float64, shape = (NPAO, 2)
        Information on primitive GTO's:
          1st element exponent and 2nd element contraction
      contspher : numpy.ndarray, dtype=numpy.intc, shape = (NAO, 2)
        Same information as ao_spherical as numpy.ndarray
      lxlylz : numpy.ndarray, dtype=numpy.intc, shape = (NAO, 3)
        Contains the expontents lx, ly, lz for the Cartesian Gaussians.

  2. Lists of dictionaries / list of tuples:

    Member of dict          Content
    --------------------    -------------------------------------
    'atom'                  Index of atom
    'pnum'                  Number of primitives
    'type'                  Type of AO
    'coeffs'                AO coefficients
    'lxlylz'                Exponents of Cartesian Gaussians
    'lm' (if spherical)     Quantum number of Spherical Gaussians

    See :ref:`Central Variables` in the manual for details.
  '''
  def __init__(self, data = None, seq = [], restart=None):
    if isinstance(data, list):
      seq = data
    elif isinstance(data, dict):
      restart = data
    
    UserList.__init__(self, seq)
    self._up_to_date = False
    self.normalized = False
    self.spherical = False

    # prim -> primitives
    # cont -> contracted

    self._cont_types = None
    self._nprim_per_cont = None
    self._prim_coeffs = None

    self._assign_prim_to_cont = None
    self._assign_cont_to_atoms = None

    self._lxlylz = None
    self._assign_lxlylz_to_cont = None
    self._nlxlylz_per_cont = None

    self._lm = None
    self._assign_lm_to_cont = None


    #self.atom_indices = None
    #self.type_list = None
    #self.pnum_list = None
    #self.ao_coeffs = None
    #self.prim2cont = None

    #self.lxlylz = None
    #self.assign_lxlylz = None
    #self.bincount_lxlylz = None
    #self.lm = None
    #self.assign_lm = None

    if restart is not None:
      self._up_to_date = True
      self.normalized           = restart['normalized']
      self.spherical            = restart['spherical']

      self._assign_cont_to_atoms         = restart['_assign_cont_to_atoms']
      self._cont_types            = restart['_cont_types']
      self._nprim_per_cont            = restart['_nprim_per_cont']
      self._prim_coeffs            = restart['_prim_coeffs']
      self._assign_prim_to_cont            = restart['_assign_prim_to_cont']

      self._lxlylz               = restart['_lxlylz']
      self._assign_lxlylz_to_cont        = restart['_assign_lxlylz_to_cont']
      self._nlxlylz_per_cont      = restart['_nlxlylz_per_cont']

      self._lm                   = restart.get('_lm', None)
      self._assign_lm_to_cont            = restart.get('_assign_lm_to_cont', None)

      self.internal_to_dict()

  def todict(self):
    self.update()
    data = {'normalized':       self.normalized,
            'spherical':        self.spherical,

            '_assign_cont_to_atoms':     self._assign_cont_to_atoms,
            '_cont_types':        self._cont_types,
            '_nprim_per_cont':        self._nprim_per_cont,
            '_prim_coeffs':        self._prim_coeffs,
            '_assign_prim_to_cont':        self._assign_prim_to_cont,

            '_lxlylz':           self._lxlylz,
            '_assign_lxlylz_to_cont':    self._assign_lxlylz_to_cont,
            '_nlxlylz_per_cont':  self._nlxlylz_per_cont,
            '_lm':               self._lm,
            '_assign_lm_to_cont':        self._assign_lm_to_cont,
            'parent_class_name' : self.__module__ + '.' + self.__class__.__name__
            }
    return data
  #def __repr__(self):
  def __str__(self):
    return '\n'.join(self.get_labels())
  def __getitem__(self, item):
    if isinstance(item,int):
      return UserList.__getitem__(self, item)
    else:
      ao_out = AOClass(seq=UserList.__getitem__(self, item))
      ao_out.update()
      return ao_out

  def __eq__(self, other):
    cases = [isinstance(other, AOClass), other == [], other is None]
    if not any(cases):
      raise TypeError('Comaring of AOClass to non AOClass object not defined')
    if cases[0]:
      self.update()
      if self.spherical:
        angular = all([numpy.allclose(self._lm, other._lm),
                       numpy.allclose(self._assign_lm_to_cont, other._assign_lm_to_cont)])
      else:
        angular = all([numpy.allclose(self._lxlylz, other._lxlylz),
                       numpy.allclose(self._assign_lxlylz_to_cont, other._assign_lxlylz_to_cont)])

      cont_types = all([i==j for i,j in zip(self._cont_types,other._cont_types)])

      same = [self.spherical == other.spherical,
              self.normalized == other.normalized,
              angular,
              cont_types,
              numpy.allclose(self._assign_cont_to_atoms, other._assign_cont_to_atoms),
              numpy.allclose(self._nprim_per_cont, other._nprim_per_cont),
              numpy.allclose(self._prim_coeffs, other._prim_coeffs),
              numpy.allclose(self._assign_prim_to_cont, other._assign_prim_to_cont),
              numpy.allclose(self._nprim_per_cont, other._nprim_per_cont),
              numpy.allclose(self._prim_coeffs , other._prim_coeffs )]
      return all(same)
    else:
      if self.data is None or len(self.data) == 0:
        return True
      else:
        return False

  def __setitem__(self, i, item):
    self.data[i] = item
    self._up_to_date = False
  def __delitem__(self, i):
    del self.data[i]
    self._up_to_date = False

  def append(self, item):
    UserList.append(self, item)
    self._up_to_date = False

  def extend(self, item):
    UserList.extend(self, item)

    self._up_to_date = False
  def remove(self, item):
    UserList.remove(self, item)
    self._up_to_date = False

  def update(self):
    '''Transfers UserList data dictionary to internal variables
    '''
    self._up_to_date = False
    self.check_members()
    self.update_ao_data()
    self.update_lxlylz()
    self.update_lm()

    self.is_normlized(force=True)

    self._up_to_date = True

  def check_members(self):
    if self.data == []:
      raise ValueError('ao_spec not initialized')
    for i,j in enumerate(self.data):
      if not isinstance(j,dict):
        raise ValueError('ao_spec[{0}] has to be a dictionary'.format(i))
      missing = []
      for key in ['atom','type','pnum','coeffs']:
        if key not in j.keys():
          missing.append(key)
      if self.spherical and 'lm' not in j.keys():
        missing.append('lm')
      if missing:
        raise ValueError('ao_spec[{0}] misses {1}'.format(i,str(missing)))

  def update_ao_data(self):
    self._assign_cont_to_atoms = []
    self._cont_types = []
    self._nprim_per_cont = []
    self._prim_coeffs = numpy.zeros((0,2))
    self._assign_prim_to_cont = []
    for i,cont in enumerate(self.data):
      self._assign_cont_to_atoms.append(cont['atom'])
      self._cont_types.append(cont['type'])
      cont_coeffs = cont['coeffs']
      self._nprim_per_cont.append(len(cont_coeffs))
      self._prim_coeffs = numpy.append(self._prim_coeffs,cont_coeffs,axis=0)
      self._assign_prim_to_cont.extend([i]*len(cont_coeffs))

    self._assign_cont_to_atoms = require(self._assign_cont_to_atoms, dtype='i')
    self._nprim_per_cont = require(self._nprim_per_cont, dtype='i')
    self._prim_coeffs = require(self._prim_coeffs, dtype='f')
    self._assign_prim_to_cont = require(self._assign_prim_to_cont, dtype='i')

  #def get_conts_are_prenormalized
  def is_normlized(self,force=False):
    '''Check if orbitals in AOClass are normalized.
    '''
    if not self._up_to_date or force:
      self.normalized = False
      conts_are_norm = []
      for cont in self.data:
        conts_are_norm.append(cont['pnum'] < 0)
      if all(conts_are_norm) != any(conts_are_norm):
        raise ValueError('Either all or none of the atomic orbitals have to be normalized!')
      self.normalized = all(conts_are_norm)
    return copy(self.normalized)

  def update_lxlylz(self):
    '''Extracts the exponents lx, ly, lz for the Cartesian Gaussians.

    **Parameters:**

    get_assign : bool, optional
      Specifies, if the index of the atomic orbital shall be returned as well.

    **Returns:**

    _lxlylz : numpy.ndarray, dtype=numpy.intc, shape = (NAO,3)
      Contains the expontents lx, ly, lz for the Cartesian Gaussians.
    _assign_lxlylz_to_cont : list of int, optional
      Contains the index of the atomic orbital in ao_spec.
    '''

    self._lxlylz = []
    self._assign_lxlylz_to_cont = []
    self._nlxlylz_per_cont = None
    for sel_ao in range(len(self.data)):
      if 'lxlylz' in self.data[sel_ao].keys():
        l = self.data[sel_ao]['lxlylz']
      else:
        l = exp[lquant[self.data[sel_ao]['type']]]
      self._lxlylz.extend(l)
      self._assign_lxlylz_to_cont.extend([sel_ao]*len(l))

    self._lxlylz = require(self._lxlylz,dtype='i')
    self._assign_lxlylz_to_cont = require(self._assign_lxlylz_to_cont, dtype='i')
    self._nlxlylz_per_cont = require(numpy.bincount(self._assign_lxlylz_to_cont), dtype='i')

    # Get label
    #if get_label: return copy(1000*self._assign_lxlylz_to_cont + (self._lxlylz * numpy.array([100,10,1])).sum(axis=1,dtype=numpy.intc))

  def update_lm(self):
    if self.spherical:
      self._lm = []
      self._assign_lm_to_cont = []
      for sel_ao in range(len(self.data)):
        for lm in self.data[sel_ao]['lm']:
          self._lm.append(lm)
          self._assign_lm_to_cont.append(sel_ao)
      self._assign_lm_to_cont = require(self._assign_lm_to_cont, dtype='i')

  def set_lm_dict(self,p=[1,0]):
    '''Sets the l,m quantum numbers for the contracted spherical harmonic
    Gaussian basis set.
    '''
    for sel_ao in range(len(self.data)):
      self.data[sel_ao]['lm'] = []
      ii = self.data[sel_ao]['type']
      l = lquant[ii]
      for m in (range(0,l+1) if l != 1 else p):
        self.data[sel_ao]['lm'].append((l,m))
        if m != 0:
          self.data[sel_ao]['lm'].append((l,-m))
    self.spherical = True
    self._up_to_date = False

  def internal_to_dict(self):
    '''Transforms Numpy-style data to lists of dictionary style data
      for compatability.
    '''
    ao_spec = []
    for ic in range(len(self._assign_cont_to_atoms)):
      ao_spec.append({'atom': self._assign_cont_to_atoms[ic],
                      'type': self._cont_types[ic],
                      'pnum': self._nprim_per_cont[ic],
                      'coeffs': self._prim_coeffs[self._assign_prim_to_cont == ic],
                      'lxlylz': self._lxlylz[self._assign_lxlylz_to_cont == ic]
                      })
      if self.spherical:
        ao_spec[-1]['lm'] = self._lm[self._assign_lm_to_cont == ic]
    self.data = ao_spec


  #def get_atom_indices(self):
  def get_assign_cont_to_atoms(self):
    if not self._up_to_date: self.update()
    return copy(self._assign_cont_to_atoms)

  #def get_type_list(self):
  def get_cont_types(self):
    if not self._up_to_date: self.update()
    return copy(self._cont_types)

  #def get_pnum_list(self):
  def get_nprim_per_cont(self):
    if not self._up_to_date: self.update()
    return copy(self._nprim_per_cont)

  #def get_ao_coeffs(self):
  def get_prim_coeffs(self):
    if not self._up_to_date: self.update()
    return copy(self._prim_coeffs)

  #def get_prim2cont(self):
  def get_assign_prim_to_cont(self):
    if not self._up_to_date: self.update()
    return copy(self._assign_prim_to_cont)

  def get_lxlylz(self):
    if not self._up_to_date: self.update()
    return copy(self._lxlylz)

  #def get_assign_lxlylz(self):
  def get_assign_lxlylz_to_cont(self):
    if not self._up_to_date: self.update()
    return copy(self._assign_lxlylz_to_cont)

  #def get_bincount_lxlylz(self):
  def get_nlxlylz_per_cont(self):
    if not self._up_to_date: self.update()
    return copy(self._nlxlylz_per_cont)

  def get_lm(self):
    if not self._up_to_date: self.update()
    return copy(self._lm)

  #def get_assign_lm(self):
  def get_assign_lm_to_cont(self):
    if not self._up_to_date: self.update()
    return copy(self._assign_lm_to_cont)

  def get_normalized(self):
    if not self._up_to_date: self.update()
    return int(self.normalized)


#This should really die sooner or later...
#For now there are still a few calls to it
#and I don't know how to solve those easily
  def get_old_ao_spherical(self):
    '''Compatability funtction to allow access to old-style ao_spherical.
    '''
    if not self._up_to_date: self.update()
    return list(zip(self.get_assign_lm_to_cont(),self.get_lm())) if self.spherical else []

  def get_labels(self):
    if not self._up_to_date: self.update()
    if self.spherical:
      labels = ['l,m=%s,atom=%d' % (self.get_lm()[i],self.get_assign_cont_to_atoms()[j])
                for i,j in enumerate(self.get_assign_lm_to_cont())]
    else:
      labels = ['lxlylz=%s,atom=%d' % (self.get_lxlylz()[i],self.get_assign_cont_to_atoms()[j])
                for i,j in enumerate(self.get_assign_lxlylz_to_cont())]
    return labels

  def get_ao_num(self):
    if not self._up_to_date: self.update()
    return len(self.get_lm()) if self.spherical else len(self.get_lxlylz())

class MOClass(UserList):
  '''MO base class which contains all information on atomic orbitals.
  Two types of dataformats are available:

  1. Numpy-style data:

      coeffs : numpy.ndarray, dtype=float64, shape = (NMO, NAO)
        Molecular orbital coefficients.
      occ : numpy.ndarray, dtype=float64, shape = (NMO)
        Occupation numbers for molecular orbitals.
      eig : numpy.ndarray, dtype=float64, shape = (NMO)
        Eigenvalues for molecular orbitals.
      sym : numpy.ndarray, dtype=str, shape = (NMO)
        MOLPRO-like symmetry label of molecular orbitals.

  2. Lists of dictionaries:

    mo_spec

    See :ref:`Central Variables` in the manual for details.
  '''
  def __init__(self, seq = [], restart=None, data = None):
    if isinstance(data, list):
      seq = data
    elif isinstance(data, dict):
      restart = data
    UserList.__init__(self, seq)
    self._up_to_date = False
    self.coeffs = None
    self.occ = None
    self.sym = None
    self.spin = None
    self.eig = None
    self.spinpolarized = False
    self.selection_string = None
    self.selected_mo = None
    self.alpha_index = None
    self.beta_index = None
    if restart is not None:
      self._up_to_date = True
      self.coeffs = restart['coeffs']
      self.occ = restart['occ']
      self.sym = restart['sym']
      self.eig = restart['eig']
      self.spinpolarized = restart['spinpolarized']
      self.selection_string = restart.get('selection_string', None)
      self.selected_mo = restart.get('selected_mo', None)
      self.alpha_index = restart['alpha_index']
      self.beta_index = restart['beta_index']
      try:
        self.spin = restart['spin']
      except KeyError:
        warnings.warn('"spin" not found in restart file', UserWarning)
        self.spin = numpy.zeros_like(self.sym,dtype=str)
        self.spin[self.alpha_index] = 'alpha'
        self.spin[self.beta_index] = 'beta'
      self.internal_to_dict()

  def todict(self):
    self.update()
    data = {'coeffs': self.coeffs,
            'selection_string': self.selection_string,
            'alpha_index': self.alpha_index,
            'beta_index': self.beta_index,
            'selected_mo': self.selected_mo,
            'spinpolarized': self.spinpolarized,
            'occ': self.occ,
            'eig': self.eig,
            'parent_class_name' : self.__module__ + '.' + self.__class__.__name__,
            'sym': self.sym,
            'spin': self.spin}
    return data

  def __getitem__(self, item):
    parse_directly = False
    if isinstance(item, (int, numpy.int64)):
      return UserList.__getitem__(self, item)
    elif isinstance(item, slice):
      return MOClass(UserList.__getitem__(self, item))
    elif isinstance(item, (list, numpy.ndarray)) or \
         (sys.version_info.major == 3 and isinstance(item, range)):
      if isinstance(item, numpy.ndarray):
        parse_directly = item.dtype in [int, numpy.int_, numpy.intc, bool, numpy.bool_]
      else:
        parse_directly = True
        for it in item:
          if not isinstance(it, (int, numpy.int_, numpy.intc, bool, numpy.bool_)):
            parse_directly = False
            break
    if parse_directly:
      item = numpy.array(item)
      if item.ndim > 1:
        raise ValueError('Only 1D arrays can be used for indexing!')
      data_out = []
      for i, c in enumerate(item):
        if isinstance(c, numpy.bool_):
          if c:
            data_out.append(self.data[i])
        else:
          data_out.append(self.data[c])
      mo_out = MOClass(data_out)
      mo_out.selected_mo = item
      mo_out.update()
      del data_out
    else:
      mo_out = self.select(item)
    return mo_out

  def __setitem__(self, i, item):
    self.data[i] = item
    self._up_to_date = False

  def __delitem__(self, i):
    del self.data[i]
    self._up_to_date = False

  def __eq__(self, other):
    cases = [isinstance(other, MOClass), other == [], other is None]
    if not any(cases):
      raise TypeError('Comaring of MOClass to non MOClass object not defined')
    if cases[0]:
      self.update()
      same = [numpy.allclose(self.coeffs, other.coeffs),
      numpy.allclose(self.eig, other.eig),
      numpy.allclose(self.occ, other.occ),
      self.compare_sym(other.sym)]
      return all(same)
    else:
      if self.data is None or len(self.data) == 0:
        return True
      else:
        return False

  def compare_sym(self, sym2):
    same = True
    for atom1, atom2 in zip(self.sym, sym2):
      if not len(atom1) == len(atom2):
        raise ValueError('sym object are of different length!')
      for i in range(len(self.sym)):
        if self.sym[i] != sym2[i]:
          same = False
    return same

  def __str__(self):
    return '\n'.join(self.get_labels())

  def splinsplit_array(self, array):
    array_alpha = array[self.alpha_index]
    array_beta = array[self.beta_index]
    return array_alpha, array_beta

  @property
  def is_energy_sorted(self):
    return all(numpy.diff(self.get_eig())>=0)
  
  def get_homo(self, tol=1e-5, sort=True):
    '''Returns index of highest occupied MO.
    '''
    if not self._up_to_date:
      self.update()
    if sort:
      self.sort_by_energy()
    
    ihomo = (self.get_occ() > tol).nonzero()[0]
    
    return  None if not len(ihomo) else ihomo[-1]

  def get_lumo(self, tol=1e-5, sort=True):
    '''Returns index of lowest unoccupied MO.
    '''
    if not self._up_to_date:
      self.update()
    if sort:    
      self.sort_by_energy()
    
    ilumo = (self.get_occ() < tol).nonzero()[0]
    
    return None if not len(ilumo) else ilumo[0]

  def get_lastbound(self, sort=True):
    '''Returns index of highest bound MO.
    '''
    if not self._up_to_date:
      self.update()
    if sort:
      self.sort_by_energy()
    
    imaxbound = (self.get_eig() <= 0.).nonzero()[0]
    
    return None if not len(imaxbound) else imaxbound[-1]

  def sort_by_sym(self):
    '''Sorts mo_spec by symmetry.
    '''
    self.data.sort()
    self.update()

  def append(self, item):
    UserList.append(self, item)
    self._up_to_date = False

  def extend(self, item):
    UserList.extend(self, item)
    self._up_to_date = False

  def remove(self, item):
    UserList.remove(self, item)
    self._up_to_date = False

  def mo_template(self):
    template = {'coeffs': None,
                'energy': None,
                'occ_num': None,
                'sym': None,
                'spin': None}
    return template

  def internal_to_dict(self):
    self.data = []
    for imo in range(len(self.occ)):
      self.data.append(self.mo_template())
      self.data[-1]['coeffs'] = self.coeffs[imo]
      self.data[-1]['energy'] = self.eig[imo]
      self.data[-1]['occ_num'] = self.occ[imo]
      self.data[-1]['sym'] = self.sym[imo]
      self.data[-1]['spin'] = self.spin[imo]
    return

  def update(self):
    self._up_to_date = False
    self.get_coeffs()
    self.get_occ()
    self.get_eig()
    self.get_sym()
    self.get_spinstate()
    self._up_to_date = True
    return

  def sort_by_energy(self):
    if self.is_energy_sorted:
      return
    warnings.warn('MOs are not sorted by energy. Sorting them...', UserWarning)
    tmp_data = []
    sort = numpy.argsort(self.get_eig())
    for s in sort:
      tmp_data.append(self.data[s])
    self.data = tmp_data
    self.update()
    return

  def get_spinstate(self):
    '''Determines whether the MOClass has alpha and beta spins and removes the _a/_b spin labels.
    For spin-paired calculations all spins are set to 'alpha'.
    '''
    if not self._up_to_date:
      self.get_sym()
      
    self.spinpolarized = False
    self.alpha_index = []
    self.beta_index = []
    spindic = {'a': self.alpha_index, 
               'b': self.beta_index, 
               'u': self.alpha_index} # unknown is counted as alpha_index
    spinlabel = {'a':'alpha','b':'beta'}
    for isym,jsym in enumerate(self.sym):
      if jsym.endswith('_a') or jsym.endswith('_b'):
        spindic[jsym[-1]].append(isym)
        self.sym[isym] = jsym[:-2]
        self.spin[isym] = spinlabel[jsym[-1]]
      else:
        spindic[self.spin[isym][0]].append(isym)
    
    self.internal_to_dict()
    if len(self.beta_index) != 0:
      self.spinpolarized = True

  def get_labels(self,format='default'):
    if format=='short':
      return ['%(sym)s' %
                  i for i in self.data]
    elif format=='print':
      return ['%(sym)s (Occ = %(occ_num).2f, E = %(energy)+.4f E_h)' %
                  i for i in self.data]
      
    elif format in ['cb','cube','vmd']:
      return ['%(sym)s,Occ=%(occ_num).1f,E=%(energy)+.2f' %
                  i for i in self.data]
    else:
      return ['%(sym)s, Occ=%(occ_num).2f, E=%(energy)+.4f E_h' %
                    i for i in self.data]

  def get_indices(self):
    if self.selected_mo is not None:
      return self.selected_mo
    else:
      return list(range(len(self.data)))
  
  def set_coeffs(self, item):
    '''Set function for numpy array version of molecular orbital coefficients.

       **Returns:**

        coeff: numpy.ndarray, dtype=numpy.float64, shape = (NMO, NAO)
    '''
    require(item, dtype=numpy.float64)
    if not self.coeffs.shape == item.shape:
      raise ValueError('Old and new arrays need to be of the same size!')
    self.coeffs = item
    self.internal_to_dict()

  def set_occ(self, item):
    '''Set function for numpy array version of molecular orbital occupancies.

       **Parameters:**

        occ: numpy.ndarray, dtype=numpy.float64, shape = (NMO)
    '''
    require(item, dtype=numpy.float64)
    if not self.occ.shape == item.shape:
      raise ValueError('Old and new arrays need to be of the same size!')
    self.occ = item
    self.internal_to_dict()

  def set_eig(self, item):
    '''Set function for numpy array version of molecular orbital energies.

       **Parameters:**

        eig: numpy.ndarray, dtype=numpy.float64, shape = (NMO)
    '''
    require(item, dtype=numpy.float64)
    if not self.eig.shape == item.shape:
      raise ValueError('Old and new arrays need to be of the same size!')
    self.eig = item
    self.internal_to_dict()

  def set_sym(self, item):
    '''Set function for numpy array version of molecular orbital symmetries.

       **Parameters:**

        sym: numpy.ndarray, dtype=numpy.str, shape = (NMO)
    '''
    require(item, dtype=str)
    if not self.sym.shape == item.shape:
      raise ValueError('Old and new arrays need to be of the same size!')
    self.sym = item
    self.internal_to_dict()

  def get_coeffs(self):
    '''Get function for numpy array version of molecular orbital coefficients.

       **Returns:**

        coeffs: numpy.ndarray, dtype=numpy.float64, shape = (NMO, NAO)
    '''
    if not self._up_to_date:
      self.coeffs = numpy.zeros(shape=(len(self.data), len(self.data[0]['coeffs'])), dtype=numpy.float64)
      for imo, mo in enumerate(self.data):
        self.coeffs[imo] = mo['coeffs']
    return copy(self.coeffs)

  def get_eig(self):
    '''Get function for numpy array version of molecular orbital energies.

       **Returns:**

        eig: numpy.ndarray, dtype=numpy.float64, shape = (NMO)
    '''
    if not self._up_to_date:
      self.eig = numpy.zeros(shape=(len(self.data)), dtype=numpy.float64)
      for imo, mo in enumerate(self.data):
        self.eig[imo] = mo['energy']
    return copy(self.eig)

  def get_occ(self, return_alpha_beta=False, return_int=False, tol_int=1e-5, sum_occ=False):
    '''Get function for numpy array version of molecular orbital occupancies.

       **Parameters:**

        return_alpha_beta: bool, optional, determines whether occupancies should be
        returned as a shape = (NMO) array (False) or a shape = (2,NMO/2) array (True)
        return_int: bool, optional, determines whether occupancies should be
        returned as a numpy.intc array. Raises an error if there are partially occupied orbitals.
        tol_int: float, optional, defines the tolerance for when an occupation number is considered an integer
        sum_occ: bool, optional, determines whether occupancies should be summed over

       **Note:**

        return_alpha_beta=True is only supported for spin-polarized calculations.

       **Returns:**

        occ: numpy.ndarray, dtype=numpy.float64, shape = (NMO)
    '''
    if not self._up_to_date:
      self.occ = numpy.zeros(shape=(len(self.data)), dtype=numpy.float64)
      for imo, mo in enumerate(self.data):
        self.occ[imo] = mo['occ_num']
    if return_int:
      for f in self.occ:
        if self.spinpolarized and 0. > f < 1.-tol_int:
          raise ValueError('Occupation numbers are not integers')
        elif not self.spinpolarized and 0. > f < 2.-tol_int:
          raise ValueError('Occupation numbers are not integers')
    if sum_occ:
      return sum(self.get_occ(return_int=True))
    else:
      if not return_alpha_beta or not self.spinpolarized:
        if not return_int:
          return copy(self.occ)
        else:
          return copy(numpy.array(self.occ, dtype=numpy.intc))
      else:
        occ_alpha = self.occ[self.alpha_index]
        occ_beta = self.occ[self.beta_index]
        if not return_int:
          return copy(numpy.array([occ_alpha,occ_beta]))
        else:
          occ_alpha = numpy.array(occ_alpha, dtype=numpy.intc)
          occ_beta = numpy.array(occ_beta, dtype=numpy.intc)
          return copy(numpy.array([occ_alpha,occ_beta]))

  def get_sym(self):
    '''Get function for numpy array version of molecular orbital symmetries.

       **Returns:**

        sym: numpy.ndarray, dtype=numpy.str, shape = (NMO)
    '''
    if not self._up_to_date:
      self.sym = []
      self.spin = []
      for imo, mo in enumerate(self.data):
        self.sym.append(mo['sym'])
        if 'spin' in mo.keys():
          self.spin.append(mo['spin'])
        else:
          self.spin.append('unknown')
      
      self.sym = numpy.array(self.sym, dtype=str)
      
    return copy(self.sym)

  def get_spin_index(self, spin):
    '''Function used to select MO's by spin. A numpy.ndarray is returned
       which contains the indexes of MO's of the selected spin.

      **Parameters:**

      spin : 'str', can be either 'alpha' or 'beta'

      **Returns:**

      alpha : numpy.ndarray, dtype=numpy.intc

    '''
    spindic = {'alpha': self.alpha_index, 'beta': self.beta_index}
    return numpy.array(spindic[spin], dtype=numpy.intc)

  def combine_restricted_orbitals(self):
    '''Combines unrestricted orbitals with same MO coefficients to restricted orbitals.
    (For example Psi4 creates a molden file with unrestricted orbitals after a ROHF calculation.)
    '''

    if not self.spinpolarized:
      # orbitals are already restricted
      return

    # assert same energy
    ener_a, ener_b = self.eig[self.alpha_index], self.eig[self.beta_index]
    assert numpy.all(ener_a == ener_b), "MO eigenvalues do not match"

    # assert same MO coefficients
    coeff_a, coeff_b = self.coeffs[self.alpha_index], self.coeffs[self.beta_index]
    assert numpy.all(coeff_a == coeff_b), "MO coeffs do not match"

    # add occupation and select alpha orbitals only
    occ = self.occ[self.alpha_index] + self.occ[self.beta_index]
    self.data = self[self.alpha_index].data
    self.spinpolarized = False
    self.update()
    self.set_occ(occ)
    self.update()

  def select(self, fid_mo_list, flatten_input=True, sort_indices=True):
    '''Selects molecular orbitals from an external file or a list of molecular
       orbital labels.

    **Parameters:**

      mo_spec :
        See :ref:`Central Variables` for details.
      fid_mo_list : str, `'all_mo'`, or list
        | If fid_mo_list is a str, specifies the filename of the molecular orbitals list.
        | If fid_mo_list is 'all_mo', creates a list containing all molecular orbitals.
        | If fid_mo_list is a list, provides a list (or a list of lists) of molecular
          orbital labels.

      flatten_input : boolean, optional
        Specifies wheter lists of lists should be flattened so a single MOClass instance can be returned rather than a list of MOClass instances.

      sort_indices : boolean, optional
        Specifies wheter list of indexes should be sorted before it is returned. This is only supported if flatten_input is set to ``True``.


    **Supported Formats:**

      Integer List (Counting from **Zero**!)::

        1       2       3
        5       4
        homo    lumo+2:lumo+4

      List with Symmetry Labels::

        1.1     2.1     1.3
        1.1     4.1
        4.1     2.3     2.1

      ``alpha`` and ``beta`` can be used together with symmetry labels to restrict the selection to orbitals of that symmetry.
      This option is not supported for integer lists. Note also that ``alpha`` and ``beta`` only restrict selection within one
      string of inputs. If you which to spin-restrict orbitlas given as a list of strings please use ``all_alpha`` or ``all_beta``.
      You may also select all MO's of a given symmetry by specifying only the symmetry label
      with a preceeding asterisc, e.g. ``*.A1`` to get all orbitals of A1 symmetry. 

    **Returns:**

      List of MOClass instances containing the selected orbitals as well as further information on the selection criteria used
      If a sinlge list is used as in input and/or flatten_input=True, an MOClass instance is returned instead.
    '''
    import re
    display('\nProcessing molecular orbital list...')
    if flatten_input and isinstance(list(fid_mo_list)[0], (list, numpy.ndarray)):
      display('\nWarning! Flattening of input lists requested!')
      display('\nIf this was not intended please set ``flatten_input=False``')
    
    mo_in_file = []
    selected_mo = []
    sym_select = False

    def ordered_set(inlist):
      outlist = []
      for i in inlist:
        if i not in outlist:
          outlist.append(i)
      return outlist

    def eval_mp(i, r):
      if r not in i:
        return i
      else:
        tmp = i.split(r)
        if r == '-':
          i = int(tmp[0]) - int(tmp[1])
        elif r == '+':
          i = int(tmp[0]) + int(tmp[1])
        else:
          raise ArithmeticError('Unknown Operation in input string.')
        return str(i)

    def remove_empty(items):
      out = []
      for i, item in enumerate(items):
        if item and not item == 'None':
          out.append(item)
      return out

    def range_separators(items):
      pos_range = []
      for i, item in enumerate(items):
        if item == '#':
          pos_range.append(i)
      pos_range = numpy.array(pos_range)
      return pos_range

    def parse_nosym(item):
      keys = {'homo': str(self.get_homo()),
              'lumo': str(self.get_lumo()),
              'last_bound': str(self.get_lastbound()),
              'lastbound': str(self.get_lastbound())}
            
      warnings = {'homo': 'occupied',
                  'lumo': 'unoccupied',
                  'last_bound': 'bound'}
      
      for key in keys:
        if key in item and key in warnings.keys() and keys[key] == 'None':
          print('No {0} orbitals are present. Ignoring respective entries...'.format(warnings[key]))
        item = item.replace(key, keys[key])

      #This seems quite insane to me... I don't really know what else to do though...
      pos_range = range_separators(remove_empty(item.replace(':',':#:').split(':')))
      item = remove_empty(item.split(':'))
      if isinstance(item, list):
        for i in range(len(item)):
          for operation in ['-','+']:
            item[i] = eval_mp(item[i], operation)
      else:
        for operation in ['-','+']:
          item = eval_mp(item, operation)
      
      if len(pos_range) > 0:
        s = 1
        if numpy.allclose(pos_range,numpy.array([0])) and len(item) == 1: #Stuff like :b
          a = 0
          b = int(item[0])
        elif numpy.allclose(pos_range,numpy.array([1])) and len(item) == 1: #Stuff like a:
          a = int(item[0])
          b = len(self.data)
        elif numpy.allclose(pos_range,numpy.array([1])) and len(item) == 2: #Stuff like a:b
          a = int(item[0])
          b = int(item[1])
        elif numpy.allclose(pos_range,numpy.array([1,3])) and len(item) == 3: #Stuff like a:b:s
          a = int(item[0])
          b = int(item[1])
          s = int(item[2])
        elif numpy.allclose(pos_range,numpy.array([0,1])) and len(item) == 1: #Stuff like ::s
          a = 0
          b = len(self.data)
          s = int(item[0])
        else:
          raise ValueError('Format not recognized')
        if a < 0:
          raise ValueError('Lower index or range out of bounds')
        if b > len(self.data):
          raise ValueError('Upper index or range out of bounds')
        item = [k for k in range(a, b, s)]
      else:
        item = [int(k) for k in item]
      return item

    def parse_sym(item):
      error = False
      if any([operation in item for operation in ['+', '-', ':']]):
        raise ValueError('Combining (mathematical) operators and symmetry' +
                         'labels is not supported')
      elif '*' in item:
        return [i for i,s in enumerate(self.get_sym()) if item.split('*.')[-1] in s]
      else:
        return [i for i,s in enumerate(self.get_sym()) if s == item]

    def parse_spin(item, all_alpha_beta):
      spindic = {0: 'all_alpha', 1: 'all_beta'}
      if isinstance(item, str):
        for i_s in range(2):
          if spindic[i_s] in item:
            all_alpha_beta[i_s] = True
            item = item.replace(spindic[i_s], '')
        if numpy.any(all_alpha_beta):
          return all_alpha_beta, item, None
        else:
          if 'alpha' in item:
            return all_alpha_beta, item.replace('alpha', ''), self.get_spin_index('alpha')
          elif 'beta' in item:
            return all_alpha_beta, item.replace('beta', ''), self.get_spin_index('beta')
          else:
            return all_alpha_beta, item, None
      elif isinstance(item, list):
        for i_s in range(2):
          if spindic[i_s] in item:
            all_alpha_beta[i_s] = True
            item = remove_from_list(item, spindic[i_s])
        if numpy.any(all_alpha_beta):
          return all_alpha_beta, item, None
        else:
          if 'alpha' in item:
            return all_alpha_beta, remove_from_list(item, 'alpha'), self.get_spin_index('alpha')
          elif 'beta' in item:
            return all_alpha_beta, remove_from_list(item, 'beta'), self.get_spin_index('beta')
          else:
            return all_alpha_beta, item, None

    
    regsplit = re.compile(r"[\s,;]")
    
    # We set these variables here for later reference
    all_alpha_beta = [False, False]

    if isinstance(fid_mo_list,str) and 'all_mo' in fid_mo_list.lower():
      all_alpha_beta, fid_mo_list, srec = parse_spin(fid_mo_list, all_alpha_beta)
      spinrestructions = [srec]
      mo_in_file_new = [[i for i in range(len(self.data))]]
      
    else:
      if isinstance(fid_mo_list,str) and not path.exists(fid_mo_list):
        if fid_mo_list == 'occupied':
          display('All occupied orbitals have been selected:')
          fid_mo_list = [map(str,(self.get_occ() > 1e-5).nonzero()[0])]
        elif fid_mo_list == 'unoccupied' or fid_mo_list == 'virtual':
          display('All unoccupied/virtual orbitals have been selected:')
          fid_mo_list = [map(str,(self.get_occ() <= 1e-5).nonzero()[0])]
        elif ',' in fid_mo_list:
          fid_mo_list = fid_mo_list.split(',')
        else:
          fid_mo_list = [fid_mo_list]
      if isinstance(fid_mo_list, list):
        for i in fid_mo_list:
          if not isinstance(i, list):
            if isinstance(i, int) or isinstance(i, float):
              i = str(int(i))
            i = re.sub(' +',' ', i)
            i = regsplit.split(i) if isinstance(i,str) else [i]
          mo_in_file.append(list(map(str,i)))
      else:
        try:
          fid=open(fid_mo_list,'r')
          flines = fid.readlines()
          fid.close()
          for line in flines:
            integer = line.replace(',',' ').split()
            mo_in_file.append(integer)
        except:
          raise IOError('The selected mo_list (%(m)s) is not valid!' %
                        {'m': fid_mo_list} + '\ne.g.\n\t1\t3\n\t2\t7\t9\n')
        fid_mo_list = mo_in_file

      # Print some information
      for i,j in enumerate(mo_in_file):
        display('\tLine %d: %s' % (i+1,', '.join(j)))
      display('')
      
      # Check if the molecular orbitals are specified by symmetry
      # (e.g. 1.1 in MOLPRO nomenclature) or
      # by the number in the input file (e.g. 1)
      mo_in_file_new = []
      spinrestructions = []
      for sublist in mo_in_file:
        all_alpha_beta, sublist, srec = parse_spin(sublist, all_alpha_beta)
        spinrestructions.append(srec)
        tmp = []
        for item in sublist:
          if '.' not in item and '*' not in item:
            if isinstance(item, str):
              tmp.extend(parse_nosym(item))
            else:
              tmp.extend(item)
          else:
            tmp.extend(parse_sym(item))
        mo_in_file_new.append(tmp)
    
    mo_in_file = []
    for isub, sublist in enumerate(mo_in_file_new):
      if numpy.any(all_alpha_beta):
        spindic = {0: 'alpha', 1: 'beta'}
        for i_s, all_spin in enumerate(all_alpha_beta):
          if all_spin:
            selected = []
            for isel in range(len(sublist)):
              if sublist[isel] in self.get_spin_index(spindic[i_s]):
                selected.append(isel)
      else:
        if spinrestructions[isub] is not None:
          selected = []
          for isel in range(len(sublist)):
            if sublist[isel] in spinrestructions[isub]:
              selected.append(isel)
        else:
          selected = range(len(sublist))
      mo_in_file.append([sublist[i] for i in selected])

    if flatten_input:
      mo_in_file = [[item for sublist in mo_in_file for item in sublist]]
      if sort_indices:
        mo_in_file = [numpy.sort([item for sublist in mo_in_file for item in sublist])]
    
    if len(mo_in_file) == 1 and flatten_input:
      mo_spec = MOClass([])
      mo_spec.spinpolarized = self.spinpolarized
      for index in mo_in_file[0]:
        mo_spec.append(self.data[index])
      mo_spec.selected_mo = mo_in_file[0]
      mo_spec.selection_string = str(fid_mo_list[0])
      # Update alpha_index and beta_index in newly created mo_spec
      mo_spec.alpha_index = require([i for i,j in enumerate(mo_in_file[0]) if j in self.alpha_index], dtype='i')
      mo_spec.beta_index = require([i for i,j in enumerate(mo_in_file[0]) if j in self.beta_index], dtype='i')

    else:
      mo_spec = []
      for i, sublist in enumerate(mo_in_file):
        mo_sub = MOClass([])
        mo_sub.spinpolarized = self.spinpolarized
        for index in sublist:
          mo_sub.append(self.data[index])
        mo_sub.selected_mo = sublist
        mo_sub.selection_string = str(fid_mo_list[i])
        # Update alpha_index and beta_index in newly created mo_spec
        mo_sub.alpha_index = require([i for i,j in enumerate(sublist) if j in self.alpha_index], dtype='i')
        mo_sub.beta_index = require([i for i,j in enumerate(sublist) if j in self.beta_index], dtype='i')
        mo_spec.append(mo_sub)

    # Print some information
    display('\nThe following orbitals will be considered...')
    for i, sublist in enumerate(mo_in_file):
      display('\tLine %d: %s' % (i+1, str(sublist)))

    display('')
    return mo_spec
