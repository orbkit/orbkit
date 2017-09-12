import numpy
from orbkit.display import display
from orbkit.units import pi
from orbkit.analytical_integrals import get_ao_dipole_matrix,get_ao_overlap
from orbkit.analytical_integrals import get_mo_overlap_matrix

class RPAClass:
  '''
  '''
  def spectrum(self,qc,periodic_axis='xyz',q='xyz',broadening=,numproc=4):
    
    # Components of spectrum
    q = ['%s' % i for i in q]
    periodic_axis = ['%s' % i for i in periodic_axis]
    non_perdiodic_axis = list(set(q) - set(periodic_axis))
    
    
    # Get MO energies and coefficients
    mo_ene = qc.mo_spec.get_eig()
    HOMO = qc.mo_spec.get_homo()
    LUMO = HOMO + 1
    
    # Calculate dipole moments in length gauge
    if non_perdiodic_axis:
      display('Calculating dipole moments in length gauge.\n')
      omr = []
      for c in non_perdiodic_axis:
        dm_aoom = get_ao_dipole_matrix(qc,component=c)
        omr.append(get_mo_overlap_matrix(coeffs,coeffs,dm_aoom,numproc=numproc))
      omr = numpy.array(omr)
    
    # Calculate dipole moments in velocity gauge
    if periodic_axis:
      display('Calculating dipole moments in velocity gauge.\n')
      omv = []
      for c in periodic_axis:
        aoom = get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec,drv=component)
        omv.append(get_mo_overlap_matrix(qc.mo_spec,qc.mo_spec,aoom,numproc=numproc))
      omv = numpy.array(omv)
      
    # 
    
    return 'lalala'
  
  def 
    
    

#L = 18.885*aa2au
#alpha = 60/180.0*pi
#A = L**2*numpy.sin(alpha)


  ## Arrange excitation energies
  #mo_ene = numpy.zeros(len(qc.mo_spec))
  ## Get HOMO and energy of MOs
  #for ii in range(len(qc.mo_spec)):
    #mo_ene[ii] = qc.mo_spec[ii]['energy']
    #if qc.mo_spec[ii]['occ_num'] != 0:
      #HOMO = ii

  #fid_dm = '/home/jochen19/PHD/sens_mos2/calc_data/mu_v.h5'
  #HDF5_file = h5py.File(fid_dm,'r')
  #omv = HDF5_file['%s/omv' % n][...] 

  #xmin = 0.0
  #xmax = 8.2
  #sigma = 0.1

  ## Read and calculate absorption spectra
  #xx = numpy.linspace(xmin,xmax,401)/hartree2eV
  #sigma /= hartree2eV
  #display('Calculate absorption spectra for the %s phenol on MoS2\n' % n)
  ## Calculating spectra
  #signal = numpy.zeros((len(component),len(xx)))
  #for xyz in range(len(component)):
    #for jj in range(HOMO+1):
      #for kk in range(HOMO+1,len(mo_ene)):
        #if jj != kk:
          #signal[xyz] += ((sigma*2)/A)*(1/((xx-mo_ene[kk]+mo_ene[jj])**2+(sigma**2)))*
          #signal[xyz] += (sigma*(-2)/A)*(1/((xx-mo_ene[jj]+mo_ene[kk])**2+(sigma**2)))*(omv[xyz,jj,kk]/(mo_ene[jj]-mo_ene[kk]))**2
  
  ## Calculate absorption spectra
  #for xyz in range(len(component)):
    #signal[xyz] *= 4*pi*xx/137

  #def get_oto(self, qc)
