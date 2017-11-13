#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>

#include "c_support.h"
#include "c_non-grid-based.h"

double get_overlap(S_Primitive *pA, S_Primitive *pB)
{
  /* Code adapted from 
  M. Hô, J. M. Hernandez-Perez: "Evaluation of Gaussian Molecular Integrals",
  DOI:10.3888/tmj.14-3
  */
  double EAB, overlap = 0.;
  double rr = 0;
  int i;
  
  for (i=0; i<3; i++)
  {
    rr += (pA->R[i] - pB->R[i]) * (pA->R[i] - pB->R[i]);
  }
  
  EAB = exp( -((pA->alpha * pB->alpha) / 
        (pA->alpha + pB->alpha)) * rr);
  overlap = EAB * pow((M_PI/(pA->alpha + pB->alpha)), 3./2.);  
  
  for (i=0; i<3; i++)
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
}
