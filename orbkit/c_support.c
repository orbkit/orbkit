#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include <stdlib.h>

double xyz(double x,double y,double z,int lx,int ly,int lz);
double get_ao_xyz(double X, double Y, double Z, int lx, int ly, int lz, 
                  double alpha, int drv);
double ao_norm(int l,int m,int n,double alpha, int is_normalized);
int doublefactorial(int n);

void c_lcreator(double* ao_list, int* exp_list, double* coeff_list,
                double* at_pos, double* x, double* y, double* z, 
                int npts, int ao_num, int pnum, int drv, int is_normalized);

void c_ao_overlap();

void c_ao_overlap()
{
  
}

void c_lcreator(double* ao_list, int* exp_list, double* coeff_list, 
                double* at_pos, double* x, double* y, double* z, 
                int npts, int ao_num , int pnum, int drv, int is_normalized)
{
  // vector grid only!
  double *Norm;
  double X, Y, Z;
  int *lx,*ly,*lz;
  double rr, *ao_l0;
  double sum;
  int i,il,ii;

  Norm = (double*) malloc(ao_num * pnum * sizeof(double));
  lx = (int*) malloc(ao_num * sizeof(int));
  ly = (int*) malloc(ao_num * sizeof(int));
  lz = (int*) malloc(ao_num * sizeof(int));
  ao_l0 = (double*) malloc(pnum * sizeof(double));

  for (il=0; il<ao_num; il++)
  {
    lx[il] = exp_list[3*il];
    ly[il] = exp_list[3*il+1];
    lz[il] = exp_list[3*il+2];
    for (ii=0; ii<pnum; ii++) 
    {
      Norm[il*pnum+ii] = ao_norm(lx[il],ly[il],lz[il],coeff_list[2*ii],is_normalized);
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
          sum += Norm[il*pnum + ii] * ao_l0[ii];
        }
        sum *= xyz(X, Y, Z, lx[il], ly[il], lz[il]);
      }
      else
      {
        for (ii=0; ii<pnum; ii++)
        {
          sum += Norm[il*pnum + ii] * ao_l0[ii] *
             get_ao_xyz(X, Y, Z, lx[il], ly[il], lz[il], coeff_list[2*ii], drv);
        }        
      }
      ao_list[il*npts+i] = sum;
    }
  }

  
  free(Norm);
  free(lx);
  free(ly);
  free(lz);
  free(ao_l0);
}
  
double xyz(double x,double y,double z,int lx,int ly,int lz)
{
  //  FUNCTION xyz  calculate x^lx * y^ly * z^lz 
  double Xl, Yl, Zl;
  
  if (lx == 0)
    Xl = 1;
  else if (lx == 1)
    Xl = x;
  else if (lx == 2)
    Xl = x*x;
  else
    Xl = pow(x, lx);
    
  if (ly == 0)
    Yl = 1;
  else if (ly == 1)
    Yl = y;
  else if (ly == 2)
    Yl = y*y;
  else
    Yl = pow(y, ly);
    
  if (lz == 0)
    Zl = 1;
  else if (lz == 1)
    Zl = z;
  else if (lz == 2)
    Zl = z*z;
  else
    Zl = pow(z, lz);

  return Xl * Yl * Zl;
}

double get_ao_xyz(double X, double Y, double Z, int lx, int ly, int lz, double alpha, int drv)
{
  /*FUNCTION get_ao_xyz  
  drv : type int
  
  0 - No derivative
  
  1 2 3 
  1 - d/dx
  2 - d/dy
  3 - d/dz
  
  4 7 8
  7 5 9
  8 9 6
  4 - d2/dx2
  5 - d2/dy2
  6 - d2/dz2
  
  7 - d2/(dx dy)
  8 - d2/(dx dz)
  9 - d2/(dy dz)
  */
  double ao_xyz = 0;
  switch(drv)      
  {
    case 0: // No derivative
    {
      ao_xyz = xyz(X, Y, Z, lx, ly, lz);
    } break;
    case 1: // d/dx
    {
      if (lx == 0)
      {
          ao_xyz = - 2 * alpha * xyz(X, Y, Z, lx+1, ly, lz);
      }
      else
      {
          ao_xyz = lx * xyz(X, Y, Z, lx-1, ly, lz) 
              - 2 * alpha * xyz(X, Y, Z, lx+1, ly, lz);    
      }
    } break;
    case 2: // d/dy
    {
      if (ly == 0)
        {
            ao_xyz = - 2 * alpha * xyz(X, Y, Z, lx, ly+1, lz);
        }
        else
        {
            ao_xyz = ly * xyz(X, Y, Z, lx, ly-1, lz) 
                - 2 * alpha * xyz(X, Y, Z, lx, ly+1, lz);   
        }
    } break;
    case 3: // d/dz
    {
      if (lz == 0)
        {
            ao_xyz = - 2 * alpha * xyz(X, Y, Z, lx, ly, lz+1);
        }
        else
        {
            ao_xyz = lz * xyz(X, Y, Z, lx, ly, lz-1) 
                - 2 * alpha * xyz(X, Y, Z, lx, ly, lz+1);   
        }
    } break;
    case 4: // d2/dx2
    {
      ao_xyz = 2 * alpha * xyz(X, Y, Z, lx, ly, lz) 
                * (2 * alpha * pow(X,2) - (2*lx+1));
      if (lx >= 2)
      { 
          ao_xyz += (pow(lx,2)-lx) * xyz(X, Y, Z, lx-2, ly, lz);
      }
    } break;
    case 5: // d2/dy2
    {
      ao_xyz = 2 * alpha * xyz(X, Y, Z, lx, ly, lz) 
                * (2 * alpha * pow(Y,2) - (2*ly+1));
      if (ly >= 2)
      {
          ao_xyz += (pow(ly,2)-ly) * xyz(X, Y, Z, lx, ly-2, lz);
      }
    } break;
    case 6: // d2/dz2
    {
      ao_xyz = 2 * alpha * xyz(X, Y, Z, lx, ly, lz) 
                * (2 * alpha * pow(Z,2) - (2*lz+1));
      if (lz >= 2)
      {
          ao_xyz += (pow(lz,2)-lz) * xyz(X, Y, Z, lx, ly, lz-2);
      }
    } break;
    case 7: // d/dx d/dy
    {
      ao_xyz = 4 * pow(alpha,2) * xyz(X, Y, Z, lx+1, ly+1, lz);        
      if (lx == 0 && ly > 0)
      {
        ao_xyz += ly * xyz(X, Y, Z, lx, ly-1, lz);
      }
      else if (lx > 0 && ly == 0)
      {
        ao_xyz += lx * xyz(X, Y, Z, lx-1, ly, lz);
      }
      else if (lx > 0 && ly > 0)
      {
        ao_xyz += lx * ly * xyz(X, Y, Z, lx-1, ly-1, lz);
      }
    } break;
    case 8: // d/dx d/dz
    {
      ao_xyz = 4 * pow(alpha,2) * xyz(X, Y, Z, lx+1, ly, lz+1);        
      if (lx == 0 && lz > 0)
      {
        ao_xyz += lz * xyz(X, Y, Z, lx, ly, lz-1);
      }
      else if (lx > 0 && lz == 0)
      {
        ao_xyz += lx * xyz(X, Y, Z, lx-1, ly, lz);
      }
      else if (lx > 0 && lz > 0)
      {
        ao_xyz += lx * lz * xyz(X, Y, Z, lx-1, ly, lz-1);
      }
    } break;
    case 9: // d/dy d/dz
    {
      ao_xyz = 4 * pow(alpha,2) * xyz(X, Y, Z, lx, ly+1, lz+1);        
      if (ly == 0 && lz > 0)
      {
        ao_xyz += lz * xyz(X, Y, Z, lx, ly, lz-1);
      }
      else if (ly > 0 && lz == 0)
      {
        ao_xyz += ly * xyz(X, Y, Z, lx, ly-1, lz);
      }
      else if (ly > 0 && lz > 0)
      {
        ao_xyz += ly * lz * xyz(X, Y, Z, lx, ly-1, lz-1);
      }
    } break;
    default:
    {
      printf("False statement for derivative variable!");
    }
  }
  return ao_xyz;
}

double ao_norm(int l,int m,int n,double alpha, int is_normalized)
{
  //  FUNCTION ao_norm  calculate normalization of uncontracted AOs --
  double norm;
  
  if (is_normalized > 0) return 1.0;
  
  norm = pow((2./M_PI),(3./4.)) * (pow(2.,(l+m+n)) * 
      pow(alpha,((2.*l+2.*m+2.*n+3.)/4.))) / (pow((doublefactorial(2*l-1)*
      doublefactorial(2*m-1) * doublefactorial(2*n-1)),0.5));
  return norm;
}

int doublefactorial(int n)
  //  FUNCTION doublefactorial  does what the name suggests -
{
  if (n <= 0) return 1;
  else return n * doublefactorial(n-2);
}
