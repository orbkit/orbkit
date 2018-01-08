import numpy, h5py
from orbkit.display import display
from orbkit.units import pi, c_in_au,ev2ha
from orbkit.analytical_integrals import get_ao_dipole_matrix,get_ao_overlap
from orbkit.analytical_integrals import get_mo_overlap,get_mo_overlap_matrix
from scipy.cluster.vq import kmeans,vq

class RPAClass:
  '''
  '''
  def __init__(self,qc):
    
    self.qc = qc
      
  def spectrum(self,unit_cell=numpy.identity(3),pbc=[False, False, False],omega=numpy.linspace(0.0,6.0,401)*ev2ha,eta=0.1*ev2ha,hdf5_fid=None):
    
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
      aoom = []
      for c in q:
        aoom.append(get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec,drv=c))
      aoom = numpy.array(aoom)
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
    for jj in range(LUMO):
      for kk in range(LUMO,len(qc.mo_spec)):
        docc = moocc[jj] - moocc[kk]
        moexene = abs(moene[jj] - moene[kk])
        if abs(moexene) >= omega[0] and abs(moexene) <= omega[-1]:
          for xyz in range(len(q)):
            if True in pbc:
              dm = get_mo_overlap(qc.mo_spec[kk],qc.mo_spec[jj],aoom[xyz])/(moene[kk] - moene[jj])
            else:
              dm = get_mo_overlap(qc.mo_spec[kk],qc.mo_spec[jj],aodm[xyz])
            signal[xyz] += ((4*pi*eta*docc)/volume)*(1/((xx-moene[kk]+moene[jj])**2+(eta**2)))*(dm)**2
            signal[xyz] += (4*pi*eta*(-docc)/volume)*(1/((xx-moene[jj]+moene[kk])**2+(eta**2)))*(dm)**2
  
    # Calculate absorption spectra
    for xyz in range(len(q)):
      signal[xyz] *= xx
      
    if hdf5_fid:
      fid = hdf5_fid if hdf5_fid.endswith('.h5') else '%s.h5' % hdf5_fid
      f = h5py.File(fid, 'w')
      f['omega'] = xx
      f['signal'] = signal
      f['volume'] = volume
      f.close()

    return signal

  def weighted_kmeans_clustering(self,points,centers,npeaks,it_max=100):
    
    # Assign points to nearest center
    for i in range(len(points['ene'])):
      distances = abs(points['ene'][i] - centers['ene'])
      idx = numpy.argmin(distances)
      points['k'][i] = idx
      centers['n'][idx] += 1
      centers['w'][idx] += points['w'][i]
    centers['ene'] = numpy.zeros(len(centers['ene']))

    #Average the points in each cluster to get a new cluster center.
    for i in range(len(points['ene'])):
      centers['ene'][points['k'][i]] += points['ene'][i] * points['w'][i]
    for i in range(len(centers['ene'])):
      centers['ene'][i] /= centers['w'][i]

    it_num = 0
    distsq = numpy.zeros(npeaks)
    while ( it_num < it_max ):
      it_num +=1
      swap = 0
      for i in range(len(points['ene'])):
        ci=points['k'][i]
        if centers['n'][ci] <= 1:
          continue

        for j in range(len(centers['ene'])):
          distpc = (abs(points['ene'][i]-centers['ene'][j])*centers['w'][j])
          if ci == j:
            distsq[j]= distpc/(centers['w'][j] - points['w'][i])
          elif centers['n'][j] == 0:
            centers['ene'][j] = numpy.copy(points['ene'][i])
            distsq[j] = 0 
          else:
            distsq[j]= distpc/(centers['w'][j] + points['w'][i])

        # Find the index of the minimum value of DISTSQ.
        nearest_cluster = numpy.argmin(distsq)

        # If that is not the cluster to which point I now belongs, move it there.
        if nearest_cluster == ci:
          continue

        j = nearest_cluster
        centers['ene'][ci] = (centers['w'][ci]*centers['ene'][ci] - points['w'][i] * points['ene'][i] ) / ( centers['w'][ci] - points['w'][i])
        centers['ene'][j] = (centers['w'][j]*centers['ene'][j] + points['w'][i]* points['ene'][i]) / ( centers['w'][j] + points['w'][i])
        centers['n'][ci] -= 1
        centers['n'][j] += 1
        centers['w'][ci] -= points['w'][i]
        centers['w'][j] += points['w'][i]

        # assign the point its new home
        points['k'][i] = j

        swap += 1
      # Exit if no reassignments were made during this iteration.
      if swap==0: 
        break
    
    return points,centers,it_num

  def get_Tij_ditos(self,dE,pbc=[False, False, False],npeaks=0):
    
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
    
    
    # Calculate AO overlap moments
    if True in pbc:
      aoom = []
      for c in q:
        aoom.append(get_ao_overlap(qc.geo_spec,qc.geo_spec,qc.ao_spec,drv=c))
      aoom = numpy.array(aoom)
    # Calculate AO dipole matrix
    else:
      aodm = []
      for c in q:
        aodm.append(get_ao_dipole_matrix(qc,component=c))
      aodm = numpy.array(aodm)
      
    if npeaks > 0 and len(dE) > 1:
      raise IOError('This combination of selected energy ranges and k-means clustering is not possible.\n')
    else:
      # Calculating transition dipole matrix
      Tij = numpy.zeros((len(dE),len(q),LUMO,len(qc.mo_spec)-LUMO))
      qc = self.qc
      for thr in range(len(dE)):
        for jj in range(LUMO):
          for kk in range(LUMO,len(qc.mo_spec)):
            if abs(exene_mat[jj,kk-LUMO]) >= dE[thr,0] and abs(exene_mat[jj,kk-LUMO]) <= dE[thr,1]:
              for xyz in range(len(q)):
                if True in pbc:
                  dm = get_mo_overlap(qc.mo_spec[jj],qc.mo_spec[kk],aoom[xyz])/(moene[jj] - moene[kk])
                else:
                  dm = get_mo_overlap(qc.mo_spec[jj],qc.mo_spec[kk],aodm[xyz])
                Tij[thr,xyz,jj,kk-LUMO] += 2*dm
        norm = numpy.linalg.norm(Tij[thr])
        Tij[thr] /= max(norm,1e-14)
      
    # k-means Clustering
    if npeaks > 0 and len(dE) == 1:
      display('K-Means Clustering for Transition Density Matrix.\n')
      
      # Variable initialization
      exene = (exene_mat).reshape((-1,))
      kmcluster = numpy.zeros((npeaks,len(q),len(exene)))
      
      # Sorting excitation energies
      bsort = numpy.argsort((exene))
      bmo = numpy.logical_and(numpy.abs(exene) >= dE[0,0],numpy.abs(exene) <= dE[0,1])
      exene *= bmo
      exene = exene[bsort]
      exnonzero = numpy.nonzero((exene))[0]
      exene = exene[exnonzero]
      
      # Sorting boolian with indices
      btmp = (numpy.arange(len(bmo))+1)*bmo
      btmp = btmp[bsort]
      btmp = btmp[exnonzero]-1
      
      # k-means Clustering
      c,tmp = kmeans(exene,npeaks)
      centroids = []
      xyzpoints = []
      
      # Sorting of transition dipole matrix
      for xyz in range(len(q)):
        # Define non-zero transition matrix
        xyzTij = (Tij[0,xyz]).reshape((-1,))
        xyzTij *= bmo
        xyzTij = xyzTij[bsort]
        xyzTij = xyzTij[exnonzero]
        # Initialize points and centers
        points = {'ene': exene,
                  'w': numpy.abs(xyzTij),
                  'k': numpy.zeros(len(exene),dtype=int)}
        centers = {'ene': c,
                  'n': numpy.zeros(len(c),dtype=int),
                  'w': numpy.zeros(len(c))}
        # Weighted k means clustering
        points,centers,it_num = self.weighted_kmeans_clustering(points,centers,npeaks,it_max=100)
        centroids.append(centers)
        xyzpoints.append(points)

        for u,v in enumerate(points['k']):
          kmcluster[v,xyz,btmp[u]] = xyzTij[u]

      kmcluster = kmcluster.reshape((npeaks,len(q),LUMO,len(qc.mo_spec)-LUMO))
      # Normalize new transition density matrix
      for thr in range(npeaks):
        norm = numpy.linalg.norm(kmcluster[thr])
        kmcluster[thr] /= max(norm,1e-14)

      return kmcluster,centroids,xyzpoints
    
    else:
      return Tij
