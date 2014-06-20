# -*- coding: iso-8859-1 -*-
from orbkit import cSupportCode
from scipy import weave
import numpy

def run_water_example():
  code = '''
  // overlap between primitive basis funciton 0,0 and 6,0

  S_Primitive A, B;

  int i = 0;
  int j = 6;

  A.alpha = ORBCOEFF2(i,0);
  B.alpha = ORBCOEFF2(j,0);
  for (int rr=0; rr<3; rr++)
  {
    A.R[rr] = R2(FCENTER1(i),rr);
    B.R[rr] = R2(FCENTER1(j),rr);
    A.l[rr] = CARTANG2(i,rr);
    B.l[rr] = CARTANG2(j,rr);
    //std::cout << A.R[rr] << i << j << std::endl;
  }

  std::cout << "Overlap between primitive basis funciton 0,0 and 6,0: ";
  std::cout << get_overlap(&A, &B) << std::endl;

  // overlap between MO 0 and MO 6

  double sum = 0;

  for (int rr=0; rr<3; rr++)
  {
    A.R[rr] = R2(FCENTER1(i),rr);
    B.R[rr] = R2(FCENTER1(j),rr);
    A.l[rr] = CARTANG2(i,rr);
    B.l[rr] = CARTANG2(j,rr);
  }
  for (int k=0; k<Norbcoeff[1]; k++)
  {
    for (int l=0; l<Norbcoeff[1]; l++)
    {
      A.alpha = ORBCOEFF2(i,k);
      B.alpha = ORBCOEFF2(j,l);
      
      sum += ao_norm(A.l[0], A.l[1], A.l[2], &A.alpha) * 
            ao_norm(B.l[0], B.l[1], B.l[2], &B.alpha) *
            PRIMCOEFF2(i,k) * PRIMCOEFF2(j,l) *
            get_overlap(&A, &B);
    }
  }
  std::cout << "Overlap between between MO 0 and MO 6: ";
  std::cout << sum << std::endl;

  // Construct the full overlap matrix

  for (int i=0; i<Norbcoeff[0]; i++)
  {
    
    for (int j=0; j<Norbcoeff[0]; j++)
    {
      sum = 0;
      for (int rr=0; rr<3; rr++)
      {
        A.R[rr] = R2(FCENTER1(i),rr);
        B.R[rr] = R2(FCENTER1(j),rr);
        A.l[rr] = CARTANG2(i,rr);
        B.l[rr] = CARTANG2(j,rr);
      }
      for (int k=0; k<Norbcoeff[1]; k++)
      {
        for (int l=0; l<Norbcoeff[1]; l++)
        {
          A.alpha = ORBCOEFF2(i,k);
          B.alpha = ORBCOEFF2(j,l);
          
          sum += ao_norm(A.l[0], A.l[1], A.l[2], &A.alpha) * 
                ao_norm(B.l[0], B.l[1], B.l[2], &B.alpha) *
                PRIMCOEFF2(i,k) * PRIMCOEFF2(j,l) *
                get_overlap(&A, &B);
        }
      }
      std::cout << sum << "\t\t";
      OVERLAP_MATRIX2(i,j) = sum;
    }
    std::cout << std::endl;
  }

  '''
  r = numpy.array([[0., 1.43233673, - 0.96104039],
      [0., - 1.43233673, - 0.96104039],
      [0., 0., 0.24026010]])

  primcoeff = numpy.array([[0.1543289673, 0.5353281423, 0.4446345422],
    [0.1543289673, 0.5353281423, 0.4446345422],
    [0.1543289673, 0.5353281423, 0.4446345422],
    [-0.09996722919, 0.3995128261, 0.7001154689],
    [0.155916275, 0.6076837186, 0.3919573931],
    [0.155916275, 0.6076837186, 0.3919573931],
    [0.155916275, 0.6076837186, 0.3919573931]])

  orbcoeff = numpy.array([[3.425250914, 0.6239137298, 0.168855404],
    [3.425250914, 0.6239137298, 0.168855404] ,
    [130.7093214, 23.80886605, 6.443608313],
    [5.033151319, 1.169596125, 0.38038896] ,
    [5.033151319, 1.169596125, 0.38038896],
    [5.033151319, 1.169596125, 0.38038896],
    [5.033151319, 1.169596125, 0.38038896]])

  fcenter = numpy.array([0,1,2,2,2,2,2])

  from orbkit.core import exp
  s = exp[0][0]
  p = exp[1] 

  cartang = numpy.array([s,s,s,s,p[0],p[1],p[2]])

  overlap_matrix = numpy.zeros((len(orbcoeff),)*2)

  arg_names = ['r', 'primcoeff', 'orbcoeff',  'fcenter', 'cartang', 'overlap_matrix']

  weave.inline(code, arg_names = arg_names, 
          support_code = cSupportCode.overlap + cSupportCode.norm,verbose = 1)
  