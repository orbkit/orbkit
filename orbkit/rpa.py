import numpy, h5py
from orbkit.display import display
from orbkit.units import pi, c_in_au,ev2ha
from orbkit.analytical_integrals import get_ao_dipole_matrix,get_ao_overlap
from orbkit.analytical_integrals import get_mo_overlap_matrix

class RPAClass:
  '''
  '''
  def __init__(self,qc,hdf5_input=None):
    
    self.qc = qc
    self.hdf5_input = hdf5_input
    if hdf5_input:
      fid = hdf5_input if hdf5_input.endswith('.h5') else '%s.h5' % hdf5_input
      HDF5_file = h5py.File(fid,'r')
      self.dm = HDF5_file['dm'][...]
    else:
      self.dm = []
      

  def spectrum(self,unit_cell=numpy.identity(3),pbc=[False, False, False],omega=numpy.linspace(0.0,6.0,401)*ev2ha,eta=0.1,hdf5_fid=None,numproc=4):
    
    display('Calculating RPA absorption spectrum.\n')
    # Calculate volume, area or length of periodic axes
    pbc = numpy.array(pbc)
    unit_cell[numpy.invert(pbc)] = numpy.identity(3)[numpy.invert(pbc)]
    volume = numpy.linalg.det(unit_cell)
    
    # Get MO energies and coefficients
    qc = self.qc
    mo_ene = qc.mo_spec.get_eig()
    HOMO = qc.mo_spec.get_homo()
    LUMO = HOMO + 1
    q = 'xyz'
    
    # Calculate dipole moments in velocity gauge
    dm = self.dm
    if dm == []:
      if True in pbc:
        display('Calculating dipole moments in velocity gauge.\n')
        for c in q:
          aoom = get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec,drv=c)
          dm.append(get_mo_overlap_matrix(qc.mo_spec,qc.mo_spec,aoom,numproc=numproc))
        dm = numpy.array(dm)

      # Calculate dipole moments in length gauge
      else:
        display('Calculating dipole moments in length gauge.\n')
        for c in q:
          dm_aoom = get_ao_dipole_matrix(qc,component=c)
          dm.append(get_mo_overlap_matrix(qc.mo_spec,qc.mo_spec,dm_aoom,numproc=numproc))
        dm = numpy.array(dm)
    
    # Calculate RPA signal
    xx = omega
    moocc = qc.mo_spec.get_occ()
    # Calculating spectra
    signal = numpy.zeros((len(q),len(xx)))
    for xyz in range(len(q)):
      for jj in range(LUMO):
        for kk in range(LUMO,len(mo_ene)):
          if jj != kk:
            docc = moocc[jj] - moocc[kk]
            signal[xyz] += ((4*pi*eta*docc)/volume)*(1/((xx-mo_ene[kk]+mo_ene[jj])**2+(eta**2)))*(dm[xyz,kk,jj]/(mo_ene[kk]-mo_ene[jj]))**2
            signal[xyz] += (4*pi*eta*(-docc)/volume)*(1/((xx-mo_ene[jj]+mo_ene[kk])**2+(eta**2)))*(dm[xyz,jj,kk]/(mo_ene[jj]-mo_ene[kk]))**2
  
    # Calculate absorption spectra
    if True in pbc:
      for xyz in range(len(q)):
        signal[xyz] *= volume*xx/c_in_au
      
    if hdf5_fid and not self.hdf5_input:
      fid = hdf5_fid if hdf5_fid.endswith('.h5') else '%s.h5' % hdf5_fid
      f = h5py.File(fid, 'w')
      f['dm'] = dm
      f['omega'] = xx
      f['signal'] = signal
      f['volume'] = volume
      f.close()

    return signal
  
  def get_Tij_otos(self,dE,pbc=[False, False, False],q = 'xyz',numproc=4):
    
    # Get MO energies and coefficients
    qc = self.qc
    mo_ene = qc.mo_spec.get_eig()
    HOMO = qc.mo_spec.get_homo()
    LUMO = HOMO + 1
    #q = 'xyz'
    
    # Calculate dipole moments in velocity gauge
    dm = self.dm
    if not dm:
      if True in pbc:
        display('Calculating dipole moments in velocity gauge.\n')
        for c in q:
          aoom = get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec,drv=c)
          dm.append(get_mo_overlap_matrix(qc.mo_spec,qc.mo_spec,aoom,numproc=numproc))
        dm = numpy.array(dm)

      # Calculate dipole moments in length gauge
      else:
        display('Calculating dipole moments in length gauge.\n')
        for c in q:
          dm_aoom = get_ao_dipole_matrix(qc,component=c)
          dm.append(get_mo_overlap_matrix(qc.mo_spec,qc.mo_spec,dm_aoom,numproc=numproc))
        dm = numpy.array(dm)
    
    # Calculating transition dipole matrix
    display('Calculating density matrix for OTOs.\n')
    Tij = numpy.zeros((len(dE),LUMO,len(qc.mo_spec)-LUMO))
    qc = self.qc
    for thr in range(len(dE)):
      for xyz in range(len(q)):
        for jj in range(LUMO):
          for kk in range(LUMO,len(qc.mo_spec)):
            if abs(mo_ene[jj] - mo_ene[kk]) >= dE[thr,0] and abs(mo_ene[jj] - mo_ene[kk]) <= dE[thr,1]:
              Tij[thr,jj,kk-LUMO] += dm[xyz,jj,kk]
            
      norm = numpy.linalg.norm(Tij[thr])
      Tij[thr] /= norm
    
    return Tij
