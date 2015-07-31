# -*- coding: iso-8859-1 -*-
'''C++ Support Code '''

math = '''
  #include <math.h>
  #define _USE_MATH_DEFINES
  '''
norm = '''
  #include <math.h>
  #define _USE_MATH_DEFINES
  
  double ao_norm(int l,int m,int n,double *alpha);
  int doublefactorial(int n);
  
  double ao_norm(int l,int m,int n,double *alpha)
  {
    //  FUNCTION ao_norm  calculate normalization of uncontracted AOs --
    double norm;
    norm = pow((2./M_PI),(3./4.)) * (pow(2.,(l+m+n)) * 
        pow(*alpha,((2.*l+2.*m+2.*n+3.)/4.))) / (pow((doublefactorial(2.*l-1.)*
        doublefactorial(2.*m-1.) * doublefactorial(2.*n-1.)),0.5));
    return norm;
  }

  int doublefactorial(int n)
    //  FUNCTION doublefactorial  does what the name suggests -
  {
  if (n <= 0)
    {
      return 1;
    }
    else
    {
      return n * doublefactorial(n-2);
    }
  }
'''
xyz = '''

  double xyz(double x,double y,double z,int lx,int ly,int lz);
  
  double xyz(double x,double y,double z,int lx,int ly,int lz)
  {
    //  FUNCTION xyz  calculate x^lx * y^ly * z^lz 
    double Xl, Yl, Zl;
    
    if (lx == 0)
      Xl = 1;
    else if (lx == 1)
      Xl = x;
    else
      Xl = pow(x, lx);
      
    if (ly == 0)
      Yl = 1;
    else if (ly == 1)
      Yl = y;
    else
      Yl = pow(y, ly);
      
    if (lz == 0)
      Zl = 1;
    else if (lz == 1)
      Zl = z;
    else
      Zl = pow(z, lz);

    return Xl * Yl * Zl;
  }

'''

ao_xyz = '''
  double get_ao_xyz(double X, double Y, double Z, int lx, int ly, int lz, double alpha, int drv);

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
        std::cout << "False statement for derivative variable!" 
                  << std::endl;
      }
    }
    return ao_xyz;
  }
'''

overlap = '''/* Code adapted from 
M. Hô, J. M. Hernandez-Perez: "Evaluation of Gaussian Molecular Integrals",
DOI:10.3888/tmj.14-3
*/
  #include <math.h>
  #define _USE_MATH_DEFINES
  
  // Structures
  struct S_Primitive 
  {
    double alpha; // Primitve exponent
    int l[3]; // Exponents lx, ly, and lz;
    double R[3]; // Center of the primitve orbital
  };
  
  // Prototypes  
  double get_overlap(S_Primitive *pA, S_Primitive *pB);
  
  double s(int i, int a, int b, S_Primitive *pA, S_Primitive *pB);
  
  double get_overlap(S_Primitive *pA, S_Primitive *pB)
  {
    
    double EAB, overlap;
    double rr = 0;
    
    for (int i=0; i<3; i++)
    {
      rr += pow((pA->R[i] - pB->R[i]), 2);
    }
    
    EAB = exp( -((pA->alpha * pB->alpha) / 
          (pA->alpha + pB->alpha)) * rr);
    overlap = EAB * pow((M_PI/(pA->alpha + pB->alpha)), 3./2.);
    
    for (int i=0; i<3; i++)
    {
      overlap *= s(i, pA->l[i], pB->l[i], pA, pB);     
    }
    
    return overlap;
  }
  
  double s(int i, int a, int b, S_Primitive *pA, S_Primitive *pB) 
  {
      // Initial Conditions
      if ((a == 0) && (b == 0))
        return 1.;
      else if ((a == 1) && (b == 0))
        return -(pA->R[i] - ((pA->alpha * pA->R[i] + pB->alpha * pB->R[i]) / 
               (pA->alpha + pB->alpha)));
      else if (b == 0) // Recurrence Index
        return -(pA->R[i] - (pA->alpha * pA->R[i] + pB->alpha * pB->R[i]) / 
               (pA->alpha + pB->alpha)) * s(i, a-1, 0, pA, pB) + (( a - 1 ) / 
               (2. * ( pA->alpha + pB->alpha ))) * s(i, a-2, 0, pA, pB);
      else // Transfer Equation
        return s(i, a+1, b-1, pA, pB) + ( pA->R[i] - pB->R[i] ) * 
               s(i, a, b-1, pA, pB);
  }'''