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