import numpy
try:
  from UserList import UserList
except ImportError:
  from collections import UserList

from .tools import *

class AOClass(UserList):
  '''
  AO base class which contains all information on atomic orbitals.
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
  def __init__(self, restart=None, seq = []):
    UserList.__init__(self, seq)
    self.up2date = False
    self.spherical = False
    self.cont2atoms = None
    self.prim2cont = None
    self.contspher = None
    self.pao = None
    self.lxlylz = None
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
            'lxlylz': self.lxlylz}
    return data
  def __getitem__(self, item):
    return UserList.__getitem__(self, item)
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
  def append(self, item):
    if 'ao_spherical' not in item:
      item['ao_spherical'] = None
    UserList.append(self, item)
    self.up2date = False
  def extend(self, item):
    if 'ao_spherical' not in item:
      item['ao_spherical'] = None
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

  def set_template(self, array, item):
    if not numpy.allclose(array, item):
      raise ValueError('Old and new arrays need to be of the same size!')
    array = item
    self.new2old()
  def set_contspher(self, item):
    require(item, dtype=numpy.intc)
    self.set_template(self.contspher, item)
  def set_cont2atoms(self, item):
    require(item, dtype=numpy.intc)
    self.set_template(self.cont2atoms, item)
  def set_prim2cont(self, item):
    require(item, dtype=numpy.intc)
    self.set_template(self.prim2cont, item)
  def set_lxlylz(self, item):
    require(item, dtype=numpy.intc)
    self.set_template(self.lxlylz, item)
  def set_pao(self, item):
    require(item, dtype=numpy.float64)
    self.set_template(self.pao, item)

  def is_normlized(self):
    if not self.up2date:
      self.normalized = False
      conts_are_norm = []
      for cont in self.data:
        conts_are_norm.append(cont['pnum'] < 0)
      if all(conts_are_norm) != any(conts_are_norm):
        raise ValueError('Either all or none of the atomic orbitals have to be normalized!')
      self.normalized = all(conts_are_norm)
    return self.normalized

#This should really die sooner or later...
#For now there are still a few calls to it
#and I don't know how to solve those easily
  def get_old_ao_spherical(self):
    old_ao_spherical = []
    for cont in self.data:
      if cont['ao_spherical'] is not None:
        old_ao_spherical.extend(cont['ao_spherical'])
    return old_ao_spherical

  def get_contspher(self):
    if not self.up2date:
      self.contspher = []
      for cont in self.data:
        if cont['ao_spherical'] is not None:
          for lm in cont['ao_spherical']:
            self.contspher.append(lm)
    self.contspher = numpy.array(self.contspher, dtype=numpy.intc)
    return self.contspher
  def get_cont2atoms(self):
    if not self.up2date:
      self.cont2atoms = numpy.zeros(shape=(len(self.data)), dtype=numpy.intc)
      for ic, contracted in enumerate(self.data):
        self.cont2atoms[ic] = contracted['atom']
    return self.cont2atoms
  def get_pao(self):
    if not self.up2date:
      self.pao = []
      for ic, contracted in enumerate(self.data):
        for ip, prim in enumerate(contracted['coeffs']):
          self.pao.append(prim)
    self.pao = numpy.array(self.pao, dtype=numpy.float64)
    return self.pao
  def get_l(self, contracted):
    if 'exp_list' in contracted:
      l = contracted['exp_list']
    else:
      l = exp[lquant[contracted['type']]]
    return l
  def get_lmpao(self):
    if not self.up2date:
      self.lmpao = []
      for ic, contracted in enumerate(self.data):
        for l in self.get_l(contracted):
          self.lmpao.extend(contracted['coeffs'])
    self.lmpao = numpy.array(self.lmpao, dtype=numpy.float64)
    return self.lmpao
  def get_prim2cont(self):
    if not self.up2date:
      self.prim2cont = []
      for ic, contracted in enumerate(self.data):
        for ipc in range(len(contracted['coeffs'])):
          self.prim2cont.append(ic)
    self.prim2cont = numpy.array(self.prim2cont, dtype=numpy.intc)
    return self.prim2cont
  #This array is nice to have becaus it allows for efficient indexing of pao's
  #We do not set it as an attribute though as it might bocome very large e.g.
  #in the case of CI/TD-DFT calculations
  def get_cont2prim(self):
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
      return self.lmprim2cont
  def get_lxlylz(self, get_assign=False, bincount=False, get_label=False):
    '''
    Extracts the exponents lx, ly, lz for the Cartesian Gaussians.

    **Parameters:**

    get_assign : bool, optional
      Specifies, if the index of the atomic orbital shall be returned as well.

    **Returns:**

    lxlylz : numpy.ndarray, dtype=numpy.intc, shape = (NAO,3)
      Contains the expontents lx, ly, lz for the Cartesian Gaussians.
    assign : list of int, optional
      Contains the index of the atomic orbital in ao_spec.
    '''
    if not self.up2date:
      self.lxlylz = []
    assign = []
    for sel_ao in range(len(self.data)):
      if 'exp_list' in self.data[sel_ao].keys():
        l = self.data[sel_ao]['exp_list']
      else:
        l = exp[lquant[self.data[sel_ao]['type']]]
      if not self.up2date:
        self.lxlylz.extend(l)
      assign.extend([sel_ao]*len(l))
    if not self.up2date:
      self.lxlylz = numpy.array(self.lxlylz,dtype=numpy.intc,order='C')
    assign = numpy.array(assign,dtype=numpy.intc,order='C')
    if get_label:
      return 1000*assign + (self.lxlylz * numpy.array([100,10,1])).sum(axis=1,dtype=numpy.intc)
    elif get_assign:
      if bincount:
        assign = numpy.bincount(assign)
      return (self.lxlylz,assign)
    return self.lxlylz



class MOClass(UserList):
  '''
  MO base class which contains all information on atomic orbitals.
  Two types of dataformats are available:

  1. Numpy-style data:

      mo_coeff : numpy.ndarray, dtype=float64, shape = (NMO, NAO) 
        Molecular orbital coefficients.
      mo_occ : numpy.ndarray, dtype=float64, shape = (NMO)
        Occupation numbers for molecular orbitals.
      mo_eig : numpy.ndarray, dtype=float64, shape = (NMO)
        Eigenvalues for molecular orbitals.
      mo_sym : numpy.ndarray, dtype=str, shape = (NMO)
        MOLPRO-like symmetry label of molecular orbitals.

  2. Lists of dictionaries:

    mo_spec

    See :ref:`Central Variables` in the manual for details.
  '''
  def __init__(self, restart=None, seq = []):
    UserList.__init__(self, seq)
    self.up2date = False
    self.mo_coeff = None
    self.mo_occ = None
    self.mo_sym = None
    self.mo_eig = None
    if restart is not None:
      self.up2date = True
      self.mo_coeff = restart['mo_coeff']
      self.mo_occ = restart['mo_occ']
      self.mo_sym = restart['mo_sym']
      self.mo_eig = restart['mo_eig']
      self.new2old()
  def todict(self):
    self.update()
    data = {'mo_coeff': self.mo_coeff,
            'mo_occ': self.mo_occ,
            'mo_eig': self.mo_eig,
            'mo_sym': self.mo_sym}
    return data
  def __getitem__(self, item):
    return UserList.__getitem__(self, item)
  def __setitem__(self, i, item):
    self.data[i] = item
    self.up2date = False
  def __eq__(self, other):
    cases = [isinstance(other, MOClass), other == [], other is None]
    if not any(cases):
      raise TypeError('Comaring of MOClass to non MOClass object not defined')
    if cases[0]:
      self.update()
      same = [numpy.allclose(self.mo_coeff, other.mo_coeff),
      numpy.allclose(self.mo_eig, other.mo_eig),
      numpy.allclose(self.mo_occ, other.mo_occ),
      self.compare_mo_sym(other.mo_sym)]
      return all(same)
    else:
      if self.data is None or len(self.data) == 0:
        return True
      else:
        return False
  def compare_mo_sym(self, sym2):
    same = True
    for atom1, atom2 in zip(self.mo_sym, sym2):
      if not len(atom1) == len(atom2):
        raise ValueError('mo_sym object are of different length!')
      for i in range(len(self.mo_sym)):
        if self.mo_sym[i] != sym2[i]:
          same = False
    return same
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
    for imo in range(len(self.mo_occ)):
      self.data.append(self.mo_template())
      self.data[-1]['coeffs'] = self.mo_coeff[imo]
      self.data[-1]['energy'] = self.mo_eig[imo]
      self.data[-1]['occ_num'] = self.mo_occ[imo]
      self.data[-1]['sym'] = self.mo_sym[imo]
    return
  def update(self):
    self.get_mo_coeff()
    self.get_mo_occ()
    self.get_mo_eig()
    self.get_mo_sym()
    self.up2date = True
    return
  def set_template(self, array, item):
    if not numpy.allclose(array, item):
      raise ValueError('Old and new arrays need to be of the same size!')
    array = item
    self.new2old()
  def set_mo_coeff(self, item):
    require(item, dtype=numpy.float64)
    self.set_template(self.mo_coeff, item)
  def set_mo_eig(self, item):
    require(item, dtype=numpy.float64)
    self.set_template(self.mo_eig, item)
  def set_mo_sym(self, item):
    require(item, dtype=str)
    self.set_template(self.mo_sym, item)

  def get_mo_coeff(self):
    if not self.up2date:
      self.mo_coeff = numpy.zeros(shape=(len(self.data), len(self.data[0]['coeffs'])), dtype=numpy.float64)
      for imo, mo in enumerate(self.data):
        self.mo_coeff[imo] = mo['coeffs']
    return self.mo_coeff
  def get_mo_eig(self):
    if not self.up2date:
      self.mo_eig = numpy.zeros(shape=(len(self.data)), dtype=numpy.float64)
      for imo, mo in enumerate(self.data):
        self.mo_eig[imo] = mo['energy']
    return self.mo_eig
  def get_mo_occ(self):
    if not self.up2date:
      self.mo_occ = numpy.zeros(shape=(len(self.data)), dtype=numpy.float64)
      for imo, mo in enumerate(self.data):
        self.mo_occ[imo] = mo['occ_num']
    return self.mo_occ
  def get_mo_sym(self):
    if not self.up2date:
      self.mo_sym = numpy.zeros(shape=(len(self.data)), dtype=str)
      for imo, mo in enumerate(self.data):
        self.mo_sym[imo] = mo['sym']
    return self.mo_sym











