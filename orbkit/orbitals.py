import numpy
from os import path
from copy import copy
import sys

try:
  from UserList import UserList
except ImportError:
  from collections import UserList

from .tools import *
from .display import display

class AOClass(UserList):
  '''AO base class which contains all information on atomic orbitals.
  Two types of dataformats are available:

  1. Numpy-style data:

      cont2atoms : numpy.ndarray, dtype=numpy.intc, shape = (NAO) 
        Transformation matrix between contracted GTO's and atoms.
      prim2cont : numpy.ndarray, dtype=numpy.intc, shape = (NPAO)
        Transformation matrix between contracted GTO's and primitive GTO's.
      pg_expcont : numpy.ndarray, dtype=float64, shape = (NPAO, 2) 
        Information on primitive GTO's:
          1st element exponent and 2nd element contraction 
      contspher : numpy.ndarray, dtype=numpy.intc, shape = (NAO, 2)
        Same information as ao_spherical as numpy.ndarray
      lxlylz : numpy.ndarray, dtype=numpy.intc, shape = (NAO, 3)
        Contains the expontents lx, ly, lz for the Cartesian Gaussians.

  2. Lists of dictionaries / list of tuples:

    ao_spec, ao_spherical

    See :ref:`Central Variables` in the manual for details.
  '''
  def __init__(self, seq = [], restart=None):
    UserList.__init__(self, seq)
    self.up2date = False
    self.spherical = False
    self.cont2atoms = None
    self.prim2cont = None
    self.contspher = None
    self.pao = None
    self.lxlylz = None
    self.assign_lxlylz = None
    self.lmpao = None
    self.lmprim2cont = None
    self.normalized = False
    if restart is not None:
      self.up2date = True
      self.spherical = restart['spherical']
      self.cont2atoms = restart['cont2atoms']
      self.prim2cont = restart['prim2cont']
      self.contspher = restart['contspher']
      self.pao = restart['pao']
      self.lxlylz = restart['lxlylz']
      self.assign_lxlylz = restart['assign_lxlylz']
      self.lmpao = None
      self.lmprim2cont = None
      self.normalized = False
      self.new2old()

  def todict(self):
    self.update()
    data = {'cont2atoms': self.cont2atoms,
            'prim2cont': self.prim2cont,
            'pao': self.pao,
            'spherical': self.spherical,
            'contspher': self.contspher,
            'lxlylz': self.lxlylz,
            'class_name': self.__module__ + '.' + self.__class__.__name__,
            'assign_lxlylz': self.assign_lxlylz}
    return data

  def __getitem__(self, item):
    if isinstance(item, int):
      return UserList.__getitem__(self, item)
    elif isinstance(item, (list,numpy.ndarray)):
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
      ao_out = AOClass(data_out)
      ao_out.update()
      del data_out
      return ao_out
    else:
      raise NotImplementedError('Only lists and arrays of integers and/or booleans are supported for array indexing!')

  def __eq__(self, other):
    cases = [isinstance(other, AOClass), other == [], other is None]
    if not any(cases):
      raise TypeError('Comaring of AOClass to non AOClass object not defined')
    if cases[0]:
      self.update()
      same = [self.spherical == other.spherical,
              numpy.allclose(self.cont2atoms, other.cont2atoms),
              numpy.allclose(self.prim2cont, other.prim2cont),
              numpy.allclose(self.contspher, other.contspher),
              numpy.allclose(self.pao, other.pao),
              numpy.allclose(self.lxlylz, other.lxlylz)]
      return all(same)
    else:
      if self.data is None or len(self.data) == 0:
        return True
      else:
        return False
  def update(self):
    self.get_contspher()
    self.get_cont2atoms()
    self.get_pao()
    self.get_lmpao()
    self.get_prim2cont()
    self.get_lmprim2cont()
    self.get_lxlylz()
    self.is_normlized()
    self.up2date = True
    return

  def __setitem__(self, i, item):
    self.data[i] = item
    self.up2date = False
  def __delitem__(self, i):
    del self.data[i]
    self.up2date = False

  def append(self, item):
    if 'ao_spherical' not in item:
      item['ao_spherical'] = []
    UserList.append(self, item)
    self.up2date = False

  def extend(self, item):
    if 'ao_spherical' not in item:
      item['ao_spherical'] = []
    UserList.extend(self, item)
    self.up2date = False

  def remove(self, item):
    UserList.remove(self, item)
    self.up2date = False

  def ao_template(self):
    template = {'atom': None,
                'pnum': None,
                'coeffs': None,
                'exp_list': None,
                'ao_spherical': None}
    return template
  def new2old(self):
    '''Transforms Numpy-style data to lists of dictionary style data
      for compatability.
    '''
    if self.cont2atoms is not None:
      ao_spec = []
      ispher = 0
      cont2prim = self.get_cont2prim()
      for ic in range(len(self.cont2atoms)):
        ao_spec.append(self.ao_template())
        ao_spec[-1]['atom'] = self.cont2atoms[ic]
        ao_spec[-1]['coeffs'] = self.pao[cont2prim[ic]]
        ao_spec[-1]['pnum'] = len(ao_spec[-1]['coeffs'])
        ao_spec[-1]['exp_list'] = tuple(self.lxlylz[ic])
        if self.contspher is not None and len(self.contspher) > 0:
          lic = self.contspher[ic]['ao_spherical'][0]
          sphericals = [tuple(cont['ao_spherical']) for cont in self.contspher[ispher:ispher+lic]]
          for lm in self.contspher[ic]['ao_spherical']:
            sphericals.append(tuple(lm))
          ao_spec[-1]['ao_spherical'] = self.contspher[ic]
      self.data = ao_spec

  def set_contspher(self, item):
    '''Set function for numpy array version of ao_spherical.

       **Parameters:**

        contspher: numpy.ndarray, dtype=numpy.intc
    '''
    require(item, dtype=numpy.intc)
    if not self.contspher.shape == item.shape:
      raise ValueError('Old and new arrays need to be of the same size!')
    self.contspher = item
    self.new2old()

  def set_cont2atoms(self, item):
    '''Set mapping between contracted GTO's and atoms.

       **Parameters:**

        cont2atoms: numpy.ndarray, dtype=numpy.intc, shape = (NAO)
    '''
    require(item, dtype=numpy.intc)
    if not self.cont2atoms.shape == item.shape:
      raise ValueError('Old and new arrays need to be of the same size!')
    self.cont2atoms = item
    self.new2old()

  def set_prim2cont(self, item):
    '''Set mapping between primitive and contracted GTO's.

       **Parameters:**

        prim2cont: numpy.ndarray, dtype=numpy.intc, shape = (NPAO)
    '''
    require(item, dtype=numpy.intc)
    if not self.prim2cont.shape == item.shape:
      raise ValueError('Old and new arrays need to be of the same size!')
    self.prim2cont = item
    self.new2old()

  def set_lxlylz(self, item):
    '''Set function for the exponents lx, ly, lz for the Cartesian Gaussians.

    **Parameters:**

    lxlylz : numpy.ndarray, dtype=numpy.intc, shape = (NAO,3)
      Contains the expontents lx, ly, lz for the Cartesian Gaussians.
    '''
    require(item, dtype=numpy.intc)
    if not self.lxlylz.shape == item.shape:
      raise ValueError('Old and new arrays need to be of the same size!')
    self.lxlylz = item
    self.new2old()

  def set_pao(self, item):
    '''Set exponents and contraction coefficients for primitive GTO's.

       **Parameters:**

        pao: numpy.ndarray, dtype=numpy.float64, shape = (NPAO, 2)
    '''
    require(item, dtype=numpy.float64)
    if not self.pao.shape == item.shape:
      raise ValueError('Old and new arrays need to be of the same size!')
    self.pao = item
    self.new2old()

  def is_normlized(self):
    '''Check if orbitals in AOClass are normalized.
    '''
    if not self.up2date:
      self.normalized = False
      conts_are_norm = []
      for cont in self.data:
        conts_are_norm.append(cont['pnum'] < 0)
      if all(conts_are_norm) != any(conts_are_norm):
        raise ValueError('Either all or none of the atomic orbitals have to be normalized!')
      self.normalized = all(conts_are_norm)
    return copy(self.normalized)

#This should really die sooner or later...
#For now there are still a few calls to it
#and I don't know how to solve those easily
  def get_old_ao_spherical(self):
    '''Compatability funtction to allow access to old-style ao_spherical.
    '''
    old_ao_spherical = []
    li = -1
    for cont in self.data:
      if cont['ao_spherical'] is not None:
        lm1 = None
        for (l,m) in cont['ao_spherical']:
          if lm1 != l:
            lm1 = l
            li += 1
          spher = [li, (l,m)]
          old_ao_spherical.append(spher)
    return old_ao_spherical

  def get_contspher(self):
    '''Get function for numpy array version of ao_spherical.

       **Returns:**

        contspher: numpy.ndarray, dtype=numpy.intc, shape = (NSPH)
    '''
    if not self.up2date:
      self.contspher = []
      for cont in self.data:
        if cont['ao_spherical'] is not None:
          for lm in cont['ao_spherical']:
            self.contspher.append(lm)
    self.contspher = numpy.array(self.contspher, dtype=numpy.intc)
    return copy(self.contspher)
  def get_cont2atoms(self):
    '''Get mapping between contracted GTO's and atoms.

       **Returns:**

        cont2atoms: numpy.ndarray, dtype=numpy.intc, shape = (NAO)
    '''
    if not self.up2date:
      self.cont2atoms = numpy.zeros(shape=(len(self.data)), dtype=numpy.intc)
      for ic, contracted in enumerate(self.data):
        self.cont2atoms[ic] = contracted['atom']
    return copy(self.cont2atoms)
  def get_pao(self):
    '''Get exponents and contraction coefficients for primitive GTO's.

       **Returns:**

        pao: numpy.ndarray, dtype=numpy.float64, shape = (NPAO, 2)
    '''
    if not self.up2date:
      self.pao = []
      for ic, contracted in enumerate(self.data):
        for ip, prim in enumerate(contracted['coeffs']):
          self.pao.append(prim)
    self.pao = numpy.array(self.pao, dtype=numpy.float64)
    return copy(self.pao)
  def get_l(self, contracted):
    if 'exp_list' in contracted:
      l = contracted['exp_list']
    else:
      l = exp[lquant[contracted['type']]]
    return l
  def get_lmpao(self):
    '''Get exponents and contraction coefficients for each spherical harmonic basis function.

       **Returns:**

        lmpao: numpy.ndarray, dtype=numpy.float64, shape = (NSPH)
    '''
    if not self.up2date:
      self.lmpao = []
      for ic, contracted in enumerate(self.data):
        for l in self.get_l(contracted):
          self.lmpao.extend(contracted['coeffs'])
    self.lmpao = numpy.array(self.lmpao, dtype=numpy.float64)
    return copy(self.lmpao)
  def get_prim2cont(self):
    '''Get mapping between primitive and contracted GTO's.

       **Returns:**

        prim2cont: numpy.ndarray, dtype=numpy.intc, shape = (NPAO)
    '''
    if not self.up2date:
      self.prim2cont = []
      for ic, contracted in enumerate(self.data):
        for ipc in range(len(contracted['coeffs'])):
          self.prim2cont.append(ic)
    self.prim2cont = numpy.array(self.prim2cont, dtype=numpy.intc)
    return copy(self.prim2cont)
  #This array is nice to have becaus it allows for efficient indexing of pao's
  #We do not set it as an attribute though as it might bocome very large e.g.
  #in the case of CI/TD-DFT calculations
  def get_cont2prim(self):
    '''Get 2D-matrix mapping between primitive and contracted GTO's.

       **Returns:**

        cont2prim: numpy.ndarray, dtype=bool, shape = (NAO, NPAO)
    '''
    if not self.up2date:
      self.get_prim2cont()
      self.get_cont2atoms()
    nc = self.cont2atoms.shape[0]
    np = self.prim2cont.shape[0]
    cont2prim = numpy.zeros(shape=(nc, np), dtype=bool)
    for ic in range(nc):
      cont2prim[ic, [numpy.argwhere(self.prim2cont == ic)]] = True
    return cont2prim
  def get_lmprim2cont(self, return_l=False):
    '''Get mapping between spherical harmonic basis functions and contracted GTO's.

      **Parameters:**
     
        return_l : bool, optional
          Return also l-index of spherical harmonic basis functions.

       **Returns:**

        return_l = False:
          lmprim2cont: numpy.ndarray, dtype=numpy.intc, shape = (NSPH)
        return_l = True:
          lmprim2cont: numpy.ndarray, dtype=numpy.intc, shape = (NSPH, NSPH)
    '''
    if not self.up2date or return_l:
      self.lmprim2cont = []
      if return_l:
        ipl = []
        il = 0
      for ic, contracted in enumerate(self.data):
        for l in self.get_l(contracted):
          self.lmprim2cont.extend([[ic]*len(contracted['coeffs'])][0])
          if return_l:
            ipl.extend([[il]*len(contracted['coeffs'])][0])
            il += 1
      self.lmprim2cont = numpy.array(self.lmprim2cont, dtype=numpy.intc)
      if return_l:
        ipl = numpy.array(ipl, dtype=numpy.intc)
        return numpy.array([[i,j] for i, j in zip(self.lmprim2cont, ipl)], dtype=numpy.intc)
    else:
      return copy(self.lmprim2cont)

  def get_lxlylz(self, get_assign=False, bincount=False, get_label=False):
    '''Extracts the exponents lx, ly, lz for the Cartesian Gaussians.

    **Parameters:**

    get_assign : bool, optional
      Specifies, if the index of the atomic orbital shall be returned as well.

    **Returns:**

    lxlylz : numpy.ndarray, dtype=numpy.intc, shape = (NAO,3)
      Contains the expontents lx, ly, lz for the Cartesian Gaussians.
    assign_lxlylz : list of int, optional
      Contains the index of the atomic orbital in ao_spec.
    '''
    if get_label:
      get_assign = True
    if not self.up2date:
      self.lxlylz = []
      self.assign_lxlylz = []
      for sel_ao in range(len(self.data)):
        if 'exp_list' in self.data[sel_ao].keys():
          l = self.data[sel_ao]['exp_list']
        else:
          l = exp[lquant[self.data[sel_ao]['type']]]
        self.lxlylz.extend(l)
        self.assign_lxlylz.extend([sel_ao]*len(l))

      self.lxlylz = numpy.array(self.lxlylz,dtype=numpy.intc,order='C')
      self.assign_lxlylz = numpy.array(self.assign_lxlylz,dtype=numpy.intc,order='C')

    if get_assign and get_label:
      return copy(1000*self.assign_lxlylz + (self.lxlylz * numpy.array([100,10,1])).sum(axis=1,dtype=numpy.intc))
    elif get_assign and not get_label:
      if bincount:
        return (copy(self.lxlylz), copy(numpy.bincount(self.assign_lxlylz)))
      else:
        return (copy(self.lxlylz), copy(self.assign_lxlylz))
    else:
      return copy(self.lxlylz)


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
  def __init__(self, seq = [], restart=None):
    UserList.__init__(self, seq)
    self.up2date = False
    self.coeffs = None
    self.occ = None
    self.sym = None
    self.eig = None
    self.spinpolarized = False
    self.selection_string = None
    self.selected_mo = None
    self.alpha_index = None
    self.beta_index = None
    if restart is not None:
      self.up2date = True
      self.coeffs = require(restart['coeffs'], dtype=numpy.float64)
      self.occ = require(restart['occ'], dtype=numpy.float64)
      self.sym = require(restart['sym'], dtype=str)
      self.eig = require(restart['eig'], dtype=numpy.float64)
      self.spinpolarized = bool(restart['spinpolarized'])
      self.selection_string = str(restart['selection_string'])
      self.alpha_index = require(restart['alpha_index'], dtype=numpy.float64)
      self.beta_index = require(restart['beta_index'], dtype=numpy.intc)
      if restart['selected_mo']:
        self.selected_mo = require(restart['selected_mo'], dtype=numpy.intc)
      else:
        self.selected_mo = None
      self.new2old()

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
            'class_name': self.__module__ + '.' + self.__class__.__name__,
            'sym': self.sym}
    return data

  def __getitem__(self, item):
    parse_directly = False
    if isinstance(item, int):
      return UserList.__getitem__(self, item)
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
    self.up2date = False

  def __delitem__(self, i):
    del self.data[i]
    self.up2date = False

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

  def splinsplit_array(self, array):
    array_alpha = array[:len(self.data)//2]
    array_beta = array[len(self.data)//2:]
    return array_alpha, array_beta

  def get_homo(self, tol=1e-5):
    '''Returns index of highest occupied MO.
    '''
    if not self.up2date:
      self.update()
    if not self.spinpolarized:
      return (self.get_occ() > tol).nonzero()[0][-1]
    else:
      occ_alpha, occ_beta = self.splinsplit_array(self.get_occ())
      return min([(occ_alpha > tol).nonzero()[0][-1],
                  (occ_beta > tol).nonzero()[0][-1]])

  def get_lumo(self, tol=1e-5):
    '''Returns index of lowest unoccupied MO.
    '''
    if not self.up2date:
      self.update()
    if not self.spinpolarized:
      ilumo = (self.get_occ() > tol).nonzero()[0][-1]+1
      if ilumo > len(self.data):
        raise ValueError('No unoccupied orbitals present!')
      else:
        return ilumo

    else:
      occ_alpha, occ_beta = self.splinsplit_array(self.get_occ())
      ilumo = max([(occ_alpha > tol).nonzero()[0][-1]+1,
                   (occ_beta > tol).nonzero()[0][-1]+1])
      if ilumo > len(self.data):
        raise ValueError('No unoccupied orbitals present!')
      else:
        return ilumo

  def get_lastbound(self):
    '''Returns index of highest bound MO.
    ''' 
    if not self.up2date:
      self.update()
    if not self.spinpolarized:
      imaxbound = (self.get_eig() <= 0.).nonzero()[0][-1]
      if imaxbound > len(self.data):
        raise ValueError('No unoccupied orbitals present!')
      else:
        return imaxbound
    else:
      eigen_alpha, eigen_beta = self.splinsplit_array(self.get_eig())
      imaxbound = max([(eigen_alpha <= 0.).nonzero()[0][-1],
                       (eigen_beta <= 0.).nonzero()[0][-1]])
      if imaxbound > len(self.data):
        raise ValueError('No unoccupied orbitals present!')
      else:
        return imaxbound   

  def sort_by_sym(self):
    '''Sorts mo_spec by symmetry.
    '''
    self.data.sort()
    self.update()

  def append(self, item):
    UserList.append(self, item)
    self.up2date = False

  def extend(self, item):
    UserList.extend(self, item)
    self.up2date = False

  def remove(self, item):
    UserList.remove(self, item)
    self.up2date = False

  def mo_template(self):
    template = {'coeffs': None,
                'energy': None,
                'occ_num': None,
                'sym': None}
    return template

  def new2old(self):
    self.data = []
    for imo in range(len(self.occ)):
      self.data.append(self.mo_template())
      self.data[-1]['coeffs'] = self.coeffs[imo]
      self.data[-1]['energy'] = self.eig[imo]
      self.data[-1]['occ_num'] = self.occ[imo]
      self.data[-1]['sym'] = self.sym[imo]
    return

  def update(self):
    self.get_coeffs()
    self.get_occ()
    self.get_eig()
    self.get_sym()
    self.up2date = True
    if self.alpha_index is None and self.beta_index is None:
      self.get_spinstate()
    return

  def get_spinstate(self):
    '''Determines whether the MOClass has alpha and beta spins and removes the _a/_b spin labels.
    For spin-paired calculations all spins are set to 'alpha'.
    '''
    if not self.up2date:
      self.get_sym()
    self.spinpolarized = False
    self.alpha_index = []
    self.beta_index = []
    spindic = {'a': self.alpha_index, 'b': self.beta_index}
    for isym in range(len(self.sym)):
      split_label = self.sym[isym].split('_')
      if len(split_label) == 2:
        spindic[split_label[-1]].append(isym)
        self.sym[isym] = split_label[0]
      else:
        spindic['a'].append(isym)
    if len(self.beta_index) != 0:
      self.spinpolarized = True

  def get_labels(self):
    return ['MO %(sym)s, Occ=%(occ_num).2f, E=%(energy)+.4f E_h' % 
                  i for i in self.data]

  def set_coeffs(self, item):
    '''Set function for numpy array version of molecular orbital coefficients.

       **Returns:**

        coeff: numpy.ndarray, dtype=numpy.float64, shape = (NMO, NAO)
    '''
    require(item, dtype=numpy.float64)
    if not self.coeffs.shape == item.shape:
      raise ValueError('Old and new arrays need to be of the same size!')
    self.coeffs = item
    self.new2old()

  def set_occ(self, item):
    '''Set function for numpy array version of molecular orbital occupancies.

       **Parameters:**

        occ: numpy.ndarray, dtype=numpy.float64, shape = (NMO)
    '''
    require(item, dtype=numpy.float64)
    if not self.occ.shape == item.shape:
      raise ValueError('Old and new arrays need to be of the same size!')
    self.occ = item
    self.new2old()

  def set_eig(self, item):
    '''Set function for numpy array version of molecular orbital energies.

       **Parameters:**

        eig: numpy.ndarray, dtype=numpy.float64, shape = (NMO)
    '''
    require(item, dtype=numpy.float64)
    if not self.eig.shape == item.shape:
      raise ValueError('Old and new arrays need to be of the same size!')
    self.eig = item
    self.new2old()

  def set_sym(self, item):
    '''Set function for numpy array version of molecular orbital symmetries.

       **Parameters:**

        sym: numpy.ndarray, dtype=numpy.str, shape = (NMO)
    '''
    require(item, dtype=str)
    if not self.sym.shape == item.shape:
      raise ValueError('Old and new arrays need to be of the same size!')
    self.sym = item
    self.new2old()

  def get_coeffs(self):
    '''Get function for numpy array version of molecular orbital coefficients.

       **Returns:**

        coeffs: numpy.ndarray, dtype=numpy.float64, shape = (NMO, NAO)
    '''
    if not self.up2date:
      self.coeffs = numpy.zeros(shape=(len(self.data), len(self.data[0]['coeffs'])), dtype=numpy.float64)
      for imo, mo in enumerate(self.data):
        self.coeffs[imo] = mo['coeffs']
    return copy(self.coeffs)

  def get_eig(self):
    '''Get function for numpy array version of molecular orbital energies.

       **Returns:**

        eig: numpy.ndarray, dtype=numpy.float64, shape = (NMO)
    '''
    if not self.up2date:
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
    if not self.up2date:
      self.occ = numpy.zeros(shape=(len(self.data)), dtype=numpy.float64)
      for imo, mo in enumerate(self.data):
        self.occ[imo] = mo['occ_num']
    if return_int:
      for f in self.occ:
        self.get_spinstate()
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
    if not self.up2date:
      self.sym = []
      for imo, mo in enumerate(self.data):
        self.sym.append(mo['sym'])
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
    spindic = {'alpha': self.alpha_index, 'beta': self.alpha_index}
    return numpy.array(spindic[spin], dtype=numpy.intc)

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
    
    **Returns:**
    
      List of MOClass instances containing the selected orbitals as well as further information on the selection criteria used
      If a sinlge list is used as in input and/or flatten_input=True, an MOClass instance is returned instead.
    '''
    import re
    display('\nProcessing molecular orbital list...')
    if flatten_input and isinstance(list(fid_mo_list)[0], (list, numpy.ndarray)):
       display('\nWarning! Flattening of input lists requested!')
    
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
        if item:
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

      for key in keys:
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
      tmp = []
      if any([operation in item for operation in ['+', '-', ':']]) \
         or '.' not in item:
         raise IOError('{0} is not a valid label according '.format(item) +
                      'to the MOLPRO nomenclature, e.g., `5.1` or `5.A1`.' +
                      '\n\tHint: You cannot mix integer numbering and MOLPRO\'s ' +
                      'symmetry labels')
      for i in numpy.argwhere(self.get_sym() == item):
        tmp.extend(i)
      return tmp

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
        if ',' in fid_mo_list:
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
          raise IOError('The selected mo-list (%(m)s) is not valid!' % 
                        {'m': fid_mo_list} + '\ne.g.\n\t1\t3\n\t2\t7\t9\n')
        fid_mo_list = mo_in_file

      # Print some information
      for i,j in enumerate(mo_in_file):
        display('\tLine %d: %s' % (i+1,', '.join(j)))
      
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
          if '.' not in item:
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

    if len(mo_in_file) == 1:
      mo_spec = MOClass([])
      for index in mo_in_file[0]:
        mo_spec.append(self.data[index])
      mo_spec.selected_mo = mo_in_file[0]
      mo_spec.selection_string = str(fid_mo_list[0])
    else:
      mo_spec = []
      for i, sublist in enumerate(mo_in_file):
        mo_sub = MOClass([])
        for index in sublist:
          mo_sub.append(self.data[index])
        mo_sub.selected_mo = sublist
        mo_sub.selection_string = str(fid_mo_list[i])
        mo_spec.append(mo_sub)

    # Print some information
    display('\nThe following orbitals will be considered...')
    for i, sublist in enumerate(mo_in_file):
      display('\tLine %d: %s' % (i+1, str(sublist)))
    
    display('')
    return mo_spec

















