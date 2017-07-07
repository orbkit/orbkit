import numpy
from os import path
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

  def set_template(self, array, item):
    '''Template for updating Numpy-style data.
    '''
    if not numpy.allclose(array, item):
      raise ValueError('Old and new arrays need to be of the same size!')
    array = item
    self.new2old()
  def set_contspher(self, item):
    '''Set function for numpy array version of ao_spherical.

       **Parameters:**

        contspher: numpy.ndarray, dtype=numpy.intc
    '''
    require(item, dtype=numpy.intc)
    self.set_template(self.contspher, item)
  def set_cont2atoms(self, item):
    '''Set mapping between contracted GTO's and atoms.

       **Parameters:**

        cont2atoms: numpy.ndarray, dtype=numpy.intc, shape = (NAO)
    '''
    require(item, dtype=numpy.intc)
    self.set_template(self.cont2atoms, item)
  def set_prim2cont(self, item):
    '''Set mapping between primitive and contracted GTO's.

       **Parameters:**

        prim2cont: numpy.ndarray, dtype=numpy.intc, shape = (NPAO)
    '''
    require(item, dtype=numpy.intc)
    self.set_template(self.prim2cont, item)
  def set_lxlylz(self, item):
    '''Set function for the exponents lx, ly, lz for the Cartesian Gaussians.

    **Parameters:**

    lxlylz : numpy.ndarray, dtype=numpy.intc, shape = (NAO,3)
      Contains the expontents lx, ly, lz for the Cartesian Gaussians.
    '''
    require(item, dtype=numpy.intc)
    self.set_template(self.lxlylz, item)
  def set_pao(self, item):
    '''Set exponents and contraction coefficients for primitive GTO's.

       **Parameters:**

        pao: numpy.ndarray, dtype=numpy.float64, shape = (NPAO, 2)
    '''
    require(item, dtype=numpy.float64)
    self.set_template(self.pao, item)

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
    return self.normalized

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
    return self.contspher
  def get_cont2atoms(self):
    '''Get mapping between contracted GTO's and atoms.

       **Returns:**

        cont2atoms: numpy.ndarray, dtype=numpy.intc, shape = (NAO)
    '''
    if not self.up2date:
      self.cont2atoms = numpy.zeros(shape=(len(self.data)), dtype=numpy.intc)
      for ic, contracted in enumerate(self.data):
        self.cont2atoms[ic] = contracted['atom']
    return self.cont2atoms
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
    return self.pao
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
    return self.lmpao
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
    return self.prim2cont
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
      return self.lmprim2cont
  def get_lxlylz(self, get_assign=False, bincount=False, get_label=False):
    '''Extracts the exponents lx, ly, lz for the Cartesian Gaussians.

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
  '''MO base class which contains all information on atomic orbitals.
  Two types of dataformats are available:

  1. Numpy-style data:

      coeff : numpy.ndarray, dtype=float64, shape = (NMO, NAO) 
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
    self.coeff = None
    self.occ = None
    self.sym = None
    self.eig = None
    self.selection_string = None
    self.selected_mo = None
    if restart is not None:
      self.up2date = True
      self.coeff = restart['coeff']
      self.occ = restart['occ']
      self.sym = restart['sym']
      self.eig = restart['eig']
      self.selection_string = restart['selection_string']
      self.selected_mo = restart['selected_mo']
      self.new2old()
  def todict(self):
    self.update()
    data = {'coeff': self.coeff,
            'selection_string': self.selection_string,
            'selected_mo': self.selected_mo,
            'occ': self.occ,
            'eig': self.eig,
            'sym': self.sym}
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
      same = [numpy.allclose(self.coeff, other.coeff),
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
  def get_homo(self, tol=1e-5):
    '''Returns index of highest occupied MO.
    '''
    return (self.get_occ() > tol).nonzero()[0][-1]
  def get_lumo(self, tol=1e-5):
    '''Returns index of lowest unoccupied MO.
    '''
    ilumo = (self.get_occ() > tol).nonzero()[0][-1]+1
    if ilumo > len(self.data):
      raise ValueError('No unoccupied orbitals present!')
    else:
      return ilumo
  def get_lastbound(self):
    '''Returns index of highest bound MO.
    ''' 
    imaxbound = (self.get_eig() <= 0.).nonzero()[0][-1]
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
      self.data[-1]['coeffs'] = self.coeff[imo]
      self.data[-1]['energy'] = self.eig[imo]
      self.data[-1]['occ_num'] = self.occ[imo]
      self.data[-1]['sym'] = self.sym[imo]
    return
  def update(self):
    self.get_coeff()
    self.get_occ()
    self.get_eig()
    self.get_sym()
    self.up2date = True
    return
  def set_template(self, array, item):
    '''Template for updating Numpy-style data.
    '''
    if not numpy.allclose(array, item):
      raise ValueError('Old and new arrays need to be of the same size!')
    array = item
    self.new2old()
  def set_coeff(self, item):
    '''Set function for numpy array version of molecular orbital symmetries.

       **Parameters:**

        sym: numpy.ndarray, dtype=numpy.str, shape = (NMO)
    '''
    require(item, dtype=numpy.float64)
    self.set_template(self.coeff, item)
  def set_occ(self, item):
    '''Set function for numpy array version of molecular orbital occupancies.

       **Parameters:**

        eig: numpy.ndarray, dtype=numpy.float64, shape = (NMO)
    '''
    require(item, dtype=numpy.float64)
    self.set_template(self.occ, item)
  def set_eig(self, item):
    '''Set function for numpy array version of molecular orbital energies.

       **Parameters:**

        eig: numpy.ndarray, dtype=numpy.float64, shape = (NMO)
    '''
    require(item, dtype=numpy.float64)
    self.set_template(self.eig, item)
  def set_sym(self, item):
    '''Set function for numpy array version of molecular orbital coefficients.

       **Parameters:**

        coeff: numpy.ndarray, dtype=numpy.float64, shape = (NMO, NAO)
    '''
    require(item, dtype=str)
    self.set_template(self.sym, item)

  def get_coeff(self):
    '''Get function for numpy array version of molecular orbital coefficients.

       **Returns:**

        coeff: numpy.ndarray, dtype=numpy.float64, shape = (NMO, NAO)
    '''
    if not self.up2date:
      self.coeff = numpy.zeros(shape=(len(self.data), len(self.data[0]['coeffs'])), dtype=numpy.float64)
      for imo, mo in enumerate(self.data):
        self.coeff[imo] = mo['coeffs']
    return self.coeff
  def get_eig(self):
    '''Get function for numpy array version of molecular orbital energies.

       **Returns:**

        eig: numpy.ndarray, dtype=numpy.float64, shape = (NMO)
    '''
    if not self.up2date:
      self.eig = numpy.zeros(shape=(len(self.data)), dtype=numpy.float64)
      for imo, mo in enumerate(self.data):
        self.eig[imo] = mo['energy']
    return self.eig
  def get_occ(self):
    '''Get function for numpy array version of molecular orbital occupancies.

       **Returns:**

        occ: numpy.ndarray, dtype=numpy.float64, shape = (NMO)
    '''
    if not self.up2date:
      self.occ = numpy.zeros(shape=(len(self.data)), dtype=numpy.float64)
      for imo, mo in enumerate(self.data):
        self.occ[imo] = mo['occ_num']
    return self.occ
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
    return self.sym

  def select(self, fid_mo_list, flatten_input=True):
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
        Specifies wheter lists of lists should be flattened so a single MOClass instance can be returned rather than a list of MOClass instances
        

    **Supported Formats:**
    
      Integer List (Counting from **ONE**!)::
      
        1       2       3
        5       4
        homo    lumo+2:lumo+4
      
      List with Symmetry Labels::
      
        1.1     2.1     1.3
        1.1     4.1
        4.1     2.3     2.1
    
    **Returns:**
    
      List of MOClass instances containing the selected orbitals as well as further information on the selection criteria used
      If a sinlge list is used as in input and/or flatten_input=True, an MOClass instance is returned instead
    
    ..attention:
      
      For **unrestricted** calculations, orbkit adds `_a` (alpha) or `_b` (beta) to
      the symmetry labels, e.g., `1.1_a`. 
      If you have specified the option `spin=alpha` or `spin=beta`, only the 
      alpha or the beta orbitals are taken into account for the counting 
      within the Integer List.
    '''
    import re
    display('\nProcessing molecular orbital list...')
    if flatten_input:
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
      tmp.extend(list(numpy.argwhere(self.get_sym() == item)[0]))
      return tmp

    regsplit = re.compile(r"[\s,;]")

    if isinstance(fid_mo_list,str) and fid_mo_list.lower() == 'all_mo':
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
            i = re.sub(' +',' ', i)
            i = regsplit.split(i) if isinstance(i,str) else [i]
          mo_in_file.append(map(str,i))
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
      for sublist in mo_in_file:
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

    mo_in_file = mo_in_file_new
    if flatten_input:
      mo_in_file = [[item for sublist in mo_in_file for item in sublist]]

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

















