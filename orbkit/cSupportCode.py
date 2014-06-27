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

overlap = '''/* Code adapted from 
M. Hô, J. M. Hernandez-Perez - "Evaluation of Gaussian Molecular Integrals"
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