import numpy, h5py
from orbkit.display import display
from orbkit.units import pi, c_in_au,ev2ha
from orbkit.analytical_integrals import get_ao_dipole_matrix,get_ao_overlap
from orbkit.analytical_integrals import get_mo_overlap,get_mo_overlap_matrix
#from scipy.cluster.vq import kmeans,vq
from scipy.signal import argrelextrema

class RPAClass:
  '''
  '''
  def __init__(self,qc):
    
    self.qc = qc
      
  def spectrum(self,unit_cell=numpy.identity(3),pbc=[False, False, False],omega=numpy.linspace(0.0,6.0,401)*ev2ha,eta=0.1*ev2ha,hdf5_fid=None,return_props=False):
    
    display('Calculating RPA absorption spectrum.\n')
    # Calculate volume, area or length of periodic axes
    pbc = numpy.array(pbc)
    unit_cell[numpy.invert(pbc)] = numpy.identity(3)[numpy.invert(pbc)]
    volume = numpy.linalg.det(unit_cell)
    
    # Get MO energies
    qc = self.qc
    moene = qc.mo_spec.get_eig()
    HOMO = qc.mo_spec.get_homo()
    LUMO = HOMO + 1
    q = 'xyz'
    
    # Calculate AO overlap moments
    if True in pbc:
      daoom = []
      for c in q:
        daoom.append(get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec,drv=c))
      daoom = numpy.array(daoom)
    # Calculate AO dipole matrix
    else:
      aodm = []
      for c in q:
        aodm.append(get_ao_dipole_matrix(qc,component=c))
      aodm = numpy.array(aodm)
    
    # Calculate RPA signal
    xx = omega
    moocc = qc.mo_spec.get_occ()
    # Calculating spectra
    signal = numpy.zeros((len(q),len(xx)))
    oscstr = []
    moexene = []
    mos = []
    for jj in range(LUMO):
      for kk in range(LUMO,len(qc.mo_spec)):
        docc = moocc[jj] - moocc[kk]
        tmpexene = abs(moene[jj] - moene[kk])
        if abs(tmpexene) >= omega[0] and abs(tmpexene) <= omega[-1]:
          moexene.append(tmpexene)
          oscstr.append([])
          mos.append([jj,kk])
          for xyz in range(len(q)):
            if True in pbc:
              dm = get_mo_overlap(qc.mo_spec[kk],qc.mo_spec[jj],daoom[xyz])/(moene[kk] - moene[jj])
            else:
              dm = get_mo_overlap(qc.mo_spec[kk],qc.mo_spec[jj],aodm[xyz])
            oscstr[-1].append(dm)
            signal[xyz] += ((4*pi*eta*docc)/volume)*(1/((xx-moene[kk]+moene[jj])**2+(eta**2)))*(dm)**2
            signal[xyz] += (4*pi*eta*(-docc)/volume)*(1/((xx-moene[jj]+moene[kk])**2+(eta**2)))*(dm)**2
    
    oscstr = numpy.array(oscstr)
    moexene = numpy.array(moexene)
    mos = numpy.array(mos)
    # Calculate absorption spectra
    for xyz in range(len(q)):
      signal[xyz] *= xx
      
    if hdf5_fid:
      fid = hdf5_fid if hdf5_fid.endswith('.h5') else '%s.h5' % hdf5_fid
      f = h5py.File(fid, 'w')
      f['omega'] = xx
      f['signal'] = signal
      f['volume'] = volume
      if return_props:
        f['oscstr'] = oscstr
        f['moexene'] = moexene
        f['mos'] = mos
      f.close()
    if return_props: return signal, oscstr, moexene
    else: return signal
  
  def lorentz(self,x,eta):
    '''Returns Lotentz distribution'''
    return (eta*0.5)/(numpy.pi * ((x)**2. + (eta*0.5)**2.))

  def get_Tij_ditos(self,dE,pbc=[False, False, False],unit_cell=numpy.identity(3),eta=0.1*ev2ha):
    
    display('Calculating Transition Density Matrix for DITOs.\n')

    # Define MO variables
    qc = self.qc
    HOMO = qc.mo_spec.get_homo()
    LUMO = HOMO + 1
    q = 'xyz'

    # Get MO energies and excitation energy matrix
    moene = qc.mo_spec.get_eig()
    exene_mat = numpy.array([moene for i in moene]) 
    exene_mat -= exene_mat.T
    exene_mat = numpy.abs(exene_mat[:LUMO,LUMO:])
    
    
    # Calculate derivatives of AO overlap matrix
    if True in pbc:
      daoom = []
      for c in q:
        daoom.append(get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec,drv=c))
      daoom = numpy.array(daoom)
    # Calculate AO dipole matrix
    else:
      aodm = []
      for c in q:
        aodm.append(get_ao_dipole_matrix(qc,component=c))
      aodm = numpy.array(aodm)
    
    # Calculate maxima in RPA spectra
    omega = numpy.linspace(dE[0],dE[1],401)
    signal = self.spectrum(unit_cell=unit_cell,pbc=pbc,omega=omega,eta=eta)
    s = (signal).sum(axis=0)#signal[2]#
    m = argrelextrema(s, numpy.greater)[0]
    maxima = {'omega': omega[m], 'signal': s[m]} 

    # Calculating transition dipole matrix
    Tij = numpy.zeros((len(m),len(q),LUMO,len(qc.mo_spec)-LUMO))
    qc = self.qc
    for pp in range(len(m)):
      for jj in range(LUMO):
        for kk in range(LUMO,len(qc.mo_spec)):
          delta = (abs(exene_mat[jj,kk-LUMO]) - maxima['omega'][pp])
          if abs(exene_mat[jj,kk-LUMO]) >= dE[0] and abs(exene_mat[jj,kk-LUMO]) <= dE[1]:
            for xyz in range(len(q)):
              if True in pbc:
                # From velocity gauge to length gauge
                dm = get_mo_overlap(qc.mo_spec[jj],qc.mo_spec[kk],daoom[xyz])/(moene[jj] - moene[kk])
              else:
                # Length gauge
                dm = get_mo_overlap(qc.mo_spec[jj],qc.mo_spec[kk],aodm[xyz])
              Tij[pp,xyz,jj,kk-LUMO] += 2*(numpy.abs(dm)**2)*self.lorentz(delta,eta)
      norm = numpy.linalg.norm(Tij[pp])
      Tij[pp] /= max(norm,1e-14)
      
    return Tij,maxima
