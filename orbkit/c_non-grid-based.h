// Structures
typedef struct prim 
{
  double alpha; // Primitve exponent
  int l[3]; // Exponents lx, ly, and lz;
  double R[3]; // Center of the primitve orbital
} S_Primitive;

// Prototypes  
double get_overlap(S_Primitive *pA, S_Primitive *pB);

double s(int i, int a, int b, S_Primitive *pA, S_Primitive *pB);
