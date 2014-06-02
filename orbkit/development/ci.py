#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
'''Module for calculating the one-electron density for ground and excited states 
and for transition flux density between various states. 
Required data is obtained from CI calculation of a QC program.
'''

import numpy
import os
from orbkit import read
from orbkit.display import display

def cis_read_gamess(filename,select_state=None):
  select_state = range(1,11)
  cis_threshold = 0.0#1e-4
  init_state = False
  with open(filename) as fileobject:
    for line in fileobject:
      thisline = line.split()             # The current line split into segments
      #--- Check the file for keywords ---
      if ' PRINTING CIS COEFFICIENTS LARGER THAN' in line:
        cis_exc = []
        cis_coeffs = []
        cis_state = []
        if float(thisline[-1]) > 0.0:
          #print('Attention: The threshold for printing of the cis coefficients is %.2f.' % (float(thisline[-1])))
          diplay('Attention: The threshold for printing of the CIS coefficients is %.2f.' % (float(thisline[-1])))
      elif ' EXCITED STATE ' in line and 'ENERGY=' and 'SPACE SYM' in line:
        if int(thisline[2]) in select_state:
          init_state = True  # # Initialize new excited state
          cis_skip = 6
          cis_state.append([thisline[2],thisline[-1]])
          cis_exc.append([])
          cis_coeffs.append([])
      if init_state == True:
        if not cis_skip:
          if '----------------------------------------------' in line:
            init_state = False
          else:
            if abs(float(thisline[2])) >= cis_threshold:
              cis_exc[-1].append(thisline[:2])
              cis_coeffs[-1].append(thisline[2])
        elif cis_skip:
          cis_skip -= 1
  for i,j in enumerate(cis_coeffs):
    j = numpy.array(j,dtype=float)
    norm = numpy.sum(j**2)
    cis_coeffs[i] = j
    
  cis_coeffs = [numpy.array(i,dtype=float) for i in cis_coeffs]
  
  return cis_exc, cis_coeffs, cis_state
  
def cis_rho_compute(filename,select_state=None):

moocc = numpy.array([numpy.array(i['occ_num'], dtype=float) for i in mo_spec])
def cis_rho(xx):
  ao_list = mdc_main.core.ao_creator(geo_spec,ao_spec,x=xx,N=(len(xx),len(y),len(z)))
  molist = numpy.array(mdc_main.core.mo_creator(ao_list,mo_spec,x=xx,N=(len(xx),len(y),len(z))))
  
  #--------------------#
  #--- Which State? ---#
  #--------------------#
  Rho = []
  
  for ii_s in range(nStates):

    ci = CIS_coeffs[ii_s]
    exc = CIS_exc[ii_s]
    #-----------------------------#
    #--- Calculate the density ---#
    #-----------------------------#

    #print 'Calculating the density'

    ##--------------------------------------------------------------#
    ##--- Calculation of CIS Density for ii as bra and jj as ket ---#
    ##--------------------------------------------------------------#
    
    rho = numpy.zeros((len(xx),len(y),len(z)))
    weave.inline(code, ['ci','exc', 'rho','molist','moocc'],verbose = 1, support_code = cSupportCode.math)
    Rho.append(rho)
  return Rho

#print 'Saving to Hierarchical Data Format file (HDF5)...'

import h5py
filename = '/home/jochen19/Master/Computational_Sections/CIS_Density/adp_pyr_dens_151_151_151'

#--- Initialize HDF5_File ---
fid = filename if filename.endswith('.h5') else '%s.h5' % filename
f = h5py.File(fid, 'w')

f.create_dataset('x',(1,len(x)),data=x)
f.create_dataset('y',(1,len(y)),data=y)
f.create_dataset('z',(1,len(z)),data=z)

f.create_dataset('geo_info',(numpy.shape(geo_info)),data=numpy.array(geo_info))
f.create_dataset('geo_spec',(numpy.shape(geo_spec)),data=geo_spec)

MO_info = f.create_group('MO_info')
occ_num=[]
energy=[]
sym=[]
for ii in range(len(mo_spec)):
  occ_num.append(mo_spec[ii]['occ_num'])
  energy.append(mo_spec[ii]['energy'])
  sym.append(mo_spec[ii]['sym'])
dset = MO_info.create_dataset('occ_num',((1,len(mo_spec))),data=occ_num)
dset = MO_info.create_dataset('energy',((1,len(mo_spec))),data=energy)
dset = MO_info.create_dataset('sym',((1,len(mo_spec))),data=sym)

for ii_s in range(nStates):  
  dID = 'rho:%d' % ii_s
  dset = f.create_dataset(dID,N)

f.close()

sNum = len(x)

#--- The number of worker processes is capped to the number of ---
#--- grid points in x-direction. --- 
if mdc_main.core.options.numproc > sNum: mdc_main.core.options.numproc = sNum

#--- Initialize some additional user information ---
status_old = 0
s_old = 0


#--- Initialize an array to store the norm of the MOs ---
mo_norm = numpy.zeros((len(mo_spec)))

#--- Start the worker processes --
pool = Pool(processes=mdc_main.core.options.numproc)

#--- Write the slices in x to an array zz ---
zz=[]
for s in range(sNum):
  zz.append((numpy.array([x[s]])))



#--- Compute the density slice by slice ---
it = pool.imap(cis_rho, zz)

write2log("\n\nThe calculation will be carried out with %d subprocesses.\n" 
                  % mdc_main.core.options.numproc)

t.append(time.time())
for s in range(sNum):
  Rho = it.next()
  f = h5py.File(fid, 'a')
  for ii_s in range(nStates):  
    dID = 'rho:%d' % ii_s
    f[dID][s,:,:] = Rho[ii_s]
  f.close()
  #--- Print out the progress of the computation ---
  
  status = s*100/float(sNum)
  if not (s+1) % mdc_main.core.options.numproc and status != status_old:
    t.append(time.time())
    write2log("\tFinished %(f).3f%% (%(s)d slices in %(t).3fs)\n" 
                    % {'f': status,
                        's': s + 1 - s_old,
                        't': t[-1]-t[-2]})
    status_old = status
    s_old = s + 1

#--- Close the worker processes --
pool.close()

code = '''
float occ;
bool bOne=false,bTwo=false;
double i_ci, j_ci;
int index=0;
for (int i=0; i<Nci[0]; i++)
{
  i_ci = ci[i];
  if (i_ci != 0.0)
  {
    for (int j=0; j<Nci[0]; j++)
    {
      j_ci = ci[j];
      if (j_ci != 0.0)
      {
        if (i == j)
        {
          for (int mm=0; mm<Nmolist[0]; mm++)
          {
            if (EXC2(i,0) == mm || EXC2(i,1) == mm) 
            {
              occ = 1.0;
            }
            else
            {
              occ = moocc[mm];
            }
            if (occ != 0.0)
            {
              for (int xx=0; xx<Nmolist[1]; xx++)
              {
                for (int yy=0; yy<Nmolist[2]; yy++)
                {
                  for (int zz=0; zz<Nmolist[3]; zz++)
                  {
                    RHO3(xx,yy,zz) += occ*pow(i_ci*MOLIST4(mm,xx,yy,zz),2);
                  }
                
                }
              }
            }
          }
        }
        else
        {
          bOne = (EXC2(i,0)-EXC2(j,0) == 0);
          bTwo =  (EXC2(i,1)-EXC2(j,1) == 0);
          if (bOne || bTwo)
          {
            if (bOne)
            {
              index = 1;
            }
            else if (bTwo)
            {
              index = 0;
            }
            else
            {
              std::cout << "False statement!" << std::endl;
            }
            for (int xx=0; xx<Nmolist[1]; xx++)
            {
              for (int yy=0; yy<Nmolist[2]; yy++)
              {
                for (int zz=0; zz<Nmolist[3]; zz++)
                {
                  RHO3(xx,yy,zz) += i_ci*j_ci*MOLIST4(EXC2(i,index),xx,yy,zz)*MOLIST4(EXC2(j,index),xx,yy,zz);
                }
              }
            }
            
          }
        }
      }
    }
  }  
}

'''