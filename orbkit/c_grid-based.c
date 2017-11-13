#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>

#include "c_support.h"
#include "c_grid-based.h"

void c_lcreator(double* ao_list, int* lxlylz, double* coeff_list, 
                double* at_pos, double* x, double* y, double* z, 
                int npts, int ao_num , int pnum, int drv, int is_normalized)
{
  // vector grid only!
  double *norm;
  double X, Y, Z;
  int *lx,*ly,*lz;
  double rr, *ao_l0;
  double sum;
  int i,il,ii;

  norm = (double*) malloc(ao_num * pnum * sizeof(double));
  lx = (int*) malloc(ao_num * sizeof(int));
  ly = (int*) malloc(ao_num * sizeof(int));
  lz = (int*) malloc(ao_num * sizeof(int));
  ao_l0 = (double*) malloc(pnum * sizeof(double));

  for (il=0; il<ao_num; il++)
  {
    lx[il] = lxlylz[3*il];
    ly[il] = lxlylz[3*il+1];
    lz[il] = lxlylz[3*il+2];
    for (ii=0; ii<pnum; ii++) 
    {
      norm[il*pnum+ii] = ao_norm(lx[il],ly[il],lz[il],coeff_list[2*ii],is_normalized);
    }
  }
  
  for (i=0; i<npts; i++)
  {
    X = x[i]-at_pos[0];
    Y = y[i]-at_pos[1];
    Z = z[i]-at_pos[2];
    rr = X*X+Y*Y+Z*Z;
    
    for (ii=0; ii<pnum; ii++)
    {
      ao_l0[ii] = coeff_list[2*ii+1] * exp(-coeff_list[2*ii] * rr);
    }
    
    for (il=0; il<ao_num; il++)
    {
      sum = 0.;
      if (drv == 0)
      {
        for (ii=0; ii<pnum; ii++)
        {
          sum += norm[il*pnum + ii] * ao_l0[ii];
        }
        sum *= xyz(X, Y, Z, lx[il], ly[il], lz[il]);
      }
      else
      {
        for (ii=0; ii<pnum; ii++)
        {
          sum += norm[il*pnum + ii] * ao_l0[ii] *
             get_ao_xyz(X, Y, Z, lx[il], ly[il], lz[il], coeff_list[2*ii], drv);
        }        
      }
      ao_list[il*npts+i] = sum;
    }
  }

  
  free(norm);
  free(lx);
  free(ly);
  free(lz);
  free(ao_l0);
}
  
