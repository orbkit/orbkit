#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>

#include "c_support.h"

double ipow(double base, int exp)
{
    double result = 1.;
    while (exp>0)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }

    return result;
}

double xyz(double x,double y,double z,int lx,int ly,int lz)
{
  //  FUNCTION xyz  calculate x^lx * y^ly * z^lz 
  return ipow(x,lx)*ipow(y,ly)*ipow(z,lz);
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
