import numpy
from copy import copy

from orbkit.display import display
from orbkit.qcinfo import QCinfo, CIinfo
from orbkit.orbitals import MOClass
from orbkit.units import ev_to_ha

from .tools import molpro_mo_order_ci
from .psi4 import psi4_detci
from .tmol import tmol_escf, tmol_tddft
from .gamess import gamess_cis, gamess_tddft
from .molpro import molpro_mcscf
from .gaussian import gaussian_tddft

def main_ci_read(qc,fname,itype='psi4_detci',threshold=0.0,
                 select=None,nforbs=0,bortho=False,
                 **kwargs):
  '''Reads determinant CI calculation. 
  
  Supported input files (``itype``):
  
  itype           QC-Program  Level of Theory
  ==============  ==========  =====================
  'psi4_detci'    PSI4        All CI calculations
  'gamess_cis'    GAMESS-US   CIS
  'gamess_tddft'  GAMESS-US   TD-DFT
  'gauss16_tddft' Gaussian16  TD-DFT
  'tmol_tddft'    TURBOMOLE   TD-DFT
  'molpro_mcscf'  MOLPRO      MCSCF
  ==============  ==========  =====================
  
  **Parameters:**
  
    qc : class QCinfo
      See :ref:`Central Variables` for details.
    fname: str, file descriptor
      Specifies the filename for the input file.
      fname can also be used with a file descriptor instad of a filename.
    itype : str, choices={'tar', 'psi4_detci', 'gamess_cis', 'tmol_tddft', 'molpro_mcscf'}
      Specifies the type of the input file.
    threshold : float, optional
      Specifies a read threshold for the CI coefficients.
    select : None or (list of) int (counting from zero), optional
      if None, reads all states, else ...
      | TURBOMOLE & GAMESS-US: Specifies the STATE to be read (0 -> ground state)
      | PSI4 & MOLPRO: Specifies the CALCULATION to be read
    nforbs : int, optional, TURBOMOLE-specific
      Specifies the number of frozen orbitals.
    bortho : bool, optional, TURBOMOLE-specific
      If True, an orthonormalization of  the TD-DFT coefficients is perfomed
      with Gram-Schmidt.
  
  **Returns:**
  
    qc : class QCinfo
      New instance of QCinfo. See :ref:`Central Variables` for details.
      qc.mo_spec is possibly reordered.
    ci : list of CIinfo class instances 
      See :ref:`Central Variables` for details.
  '''

  if isinstance(fname, str):
    filename = fname
  else:
    filename = fname.name

  display('Opened \n\t%s\n' % filename)  
  assert isinstance(qc,QCinfo), '`qc` has to be an instance of the QCinfo class'
  
  # What kind of input file has to be read?
  reader = {'psi4_detci': psi4_detci,
            'gamess_cis': gamess_cis,
            'gamess_tddft': gamess_tddft,
            'gaussian_tddft': gaussian_tddft,
            'tmol_tddft': tmol_tddft,
            'molpro_mcscf': molpro_mcscf}
  
  display('Loading %s file...' % itype)
  if itype not in reader.keys():
    display('Available reader (`itype`) are:\n  ' + ', '.join(reader.keys()))
    raise NotImplementedError("itype='%s' not implemented!"%itype)
  
  kwargs['nmoocc'] = qc.mo_spec.get_lumo()
  
  ci = reader[itype](fname,select_state=select,threshold=threshold, 
                     select_run=select,                 # PSI4/MOLPRO specific
                     nforbs=nforbs,bortho=bortho,       # TURBOMOLE specific
                     **kwargs)
  
  # Get a copy of qc
  qc = qc.copy()
  if itype in ['tmol_tddft','gamess_cis','gamess_tddft','gaussian_tddft']: # CIS-like
    moocc = qc.mo_spec.get_occ(return_int=True)
  elif itype in ['psi4_detci','molpro_mcscf']: # detCI-like
    # Reorder qc.mo_spec
    irreps=None if 'irreps' not in ci[0].info.keys() else ci[0].info['irreps']
    closed,active,external = molpro_mo_order_ci(ci[0].info['occ_info'],
                                                qc.mo_spec,
                                                irreps=irreps
                                                )
    qc.mo_spec = closed+active+external
    qc.mo_spec = MOClass(qc.mo_spec)
    moocc = numpy.zeros(len(closed),dtype=numpy.intc) + 2
  
  # Add moocc to CI class instances
  for i in ci:
    i.set_moocc(moocc)
  
  return qc,ci
