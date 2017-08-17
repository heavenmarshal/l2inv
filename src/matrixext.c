#include<math.h>
#include "matrix.h"

void onev(double* v, unsigned int n)
{
  unsigned int i;
  for(i=0; i<n; ++i) v[i] = 1.0;
}
double* new_one_vector(unsigned int n)
{
  double *v;
  v = new_vector(n);
  onev(v,n);
  return v;
}
void sum_scalar(double* vin, double* vout, double scalar, unsigned int n)
{
  unsigned int i;
  for(i=0; i<n; ++i)
    vout[i] = vin[i] + scalar;
}
void sum_scalar_inplace(double* v, double scalar, unsigned int n)
{
  unsigned int i;
  for(i=0; i<n; ++i)
    v[i] += scalar;
}
double* new_reduce_row(double **mat, unsigned int n, unsigned int m)
{
  int i,j;
  double *v = new_zero_vector(m);
  for(i = 0; i < m; ++i)
    for(j = 0; j < n; ++j)
      v[i] += mat[j][i];
  return v;
}

void vector_exp_check(double *v, int n)
{
  int i = 0;
  for (i = 0; i < n; i++) {
    if (v[i] < -500) v[i] = 0;
    else v[i] = exp(v[i]);
  }
}

void vector_log(double *v, int n)
{
  int i = 0;
  for (i = 0; i < n; i++) {
    v[i] = log(v[i]);
  }
}
