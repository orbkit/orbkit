def molden_writer(qc,filename='new'):
  '''Creates a molden file with information concerning the molecular
  structure, basis set, and molecular expansion coefficients.

  **Parameters:**

  qc : class or dict
    QCinfo class or dictionary containing the following attributes/keys.
    See :ref:`Central Variables` for details.
  filename : str
    Contains the base name of the output file.
  '''

  if not filename.endswith('.molden'):
    filename += '.molden'

  # Open an empty file
  fid = open('%(f)s' % {'f': filename},'w')

  #Write HEADER
  fid.write('[Molden Format]\n')

  #Write GEOMETRY SPECIFICATION
  fid.write('[Atoms] AU\n')
  string = ''
  for i in range(len(qc.geo_spec)):
    string += '%3s %5d %5d %12.8f %12.8f %12.8f\n' %  (qc.geo_info[i][0],int(qc.geo_info[i][1]),
                                                     int(float(qc.geo_info[i][2])),qc.geo_spec[i][0],
                                                     qc.geo_spec[i][1],qc.geo_spec[i][2])
  fid.write(str(string))

  #Write BASIS SET
  fid.write('[GTO]')
  anum = 0
  string = ''
  for i in range(len(qc.ao_spec)):
    if anum != (qc.ao_spec[i]['atom']+1):
      string += '\n%5i 0\n' % (qc.ao_spec[i]['atom']+1)
    string += '%2s %6i 1.00\n' % (qc.ao_spec[i]['type'],qc.ao_spec[i]['pnum'])
    for j in range(len(qc.ao_spec[i]['coeffs'])):
      string += '%18.10E %18.10E\n' % (qc.ao_spec[i]['coeffs'][j][0],qc.ao_spec[i]['coeffs'][j][1])
    anum = (qc.ao_spec[i]['atom']+1)

  fid.write(str(string) + '\n')

  #Write MOs
  fid.write('[MO]\n')
  string = ''
  sym = qc.mo_spec.get_sym()
  ene = qc.mo_spec.get_eig()
  occ = qc.mo_spec.get_occ()
  coeff = qc.mo_spec.get_coeffs()
  spindic = {0: 'alpha', 1: 'beta'}
  alpha = qc.mo_spec.get_spin_index('alpha')
  beta = qc.mo_spec.get_spin_index('beta')
  for i_s, spin in enumerate([alpha, beta]):
    for i in spin:
      string += ' Sym= %s\n'              % sym[i]
      string += ' Ene= %10.8f\n'          % ene[i]
      string += ' Spin= %s\n'             % spindic[i_s]
      string += ' Occup= %10.8f\n'        % occ[i]
      for j, c in enumerate(coeff[i]):
        string += '%10i %18.10E\n'        % (j+1, c)

  fid.write(str(string))

  #Close the file
  fid.close()
