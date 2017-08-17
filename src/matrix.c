/********************************************************************************
 *
 * Bayesian Regression and Adaptive Sampling with Gaussian Process Trees
 * Copyright (C) 2005, University of California
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Questions? Contact Robert B. Gramacy (rbgramacy@ams.ucsc.edu)
 *
 ********************************************************************************/


#include "rhelp.h"

#include "matrix.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>


#define DEBUG

/*
 * replace matrix with zeros
 */

void zero(double **M, unsigned int n1, unsigned int n2)
{
  unsigned int i, j;
  for(i=0; i<n1; i++) for(j=0; j<n2; j++) M[i][j] = 0;
}


/*
 * replace square matrix with identitiy
 */

void id(double **M, unsigned int n)
{
  unsigned int i;
  zero(M, n, n);
  for(i=0; i<n; i++) M[i][i] = 1.0;
}


/*
 * same as new_matrix below, but for creating
 * n x n identity matrices
 */

double ** new_id_matrix(unsigned int n)
{
  unsigned int i;
  double** m = new_zero_matrix(n, n);
  for(i=0; i<n; i++) m[i][i] = 1.0;
  return m;
}


/*
 * same as new_matrix below, but zeros out the matrix
 */

double ** new_zero_matrix(unsigned int n1, unsigned int n2)
{
  unsigned int i, j;
  double **m = new_matrix(n1, n2);
  for(i=0; i<n1; i++) for(j=0; j<n2; j++) m[i][j] = 0.0;
  return m;
}


/*
 * create a new n1 x n2 matrix which is allocated like
 * and n1*n2 array, but can be referenced as a 2-d array
 */

double ** new_matrix(unsigned int n1, unsigned int n2)
{
  int i;
  double **m;

  if(n1 == 0 || n2 == 0) return NULL;

  m = (double**) malloc(sizeof(double*) * n1);
  assert(m);
  m[0] = (double*) malloc(sizeof(double) * (n1*n2));
  assert(m[0]);

  for(i=1; i<n1; i++) m[i] = m[i-1] + n2;

  return m;
}


/*
 * create a double ** Matrix from a double * vector
 * should be freed with the free command, rather than
 * delete_matrix
 */

double ** new_matrix_bones(double *v, unsigned int n1, unsigned int n2)
{
  double **M;
  int i;
  M = (double **) malloc(sizeof(double*) * n1);
  M[0] = v;
  for(i=1; i<n1; i++) M[i] = M[i-1] + n2;
  return(M);
}


/*
 * create a new n1 x n2 matrix which is allocated like
 * and n1*n2 array, and copy the of n1 x n2 M into it.
 */

double ** new_dup_matrix(double** M, unsigned int n1, unsigned int n2)
{
  double **m;

  if(n1 <= 0 || n2 <= 0) {
    /* assert(M == NULL); */
    return NULL;
  }

  m = new_matrix(n1, n2);
  dup_matrix(m, M, n1, n2);
  return m;
}


/*
 * copy M2 to M1
 */

void dup_matrix(double** M1, double **M2, unsigned int n1, unsigned int n2)
{
  unsigned int i;
  if(n1 == 0 || n2 == 0) return;
  assert(M1 && M2);
  for(i=0; i<n1; i++) dupv(M1[i], M2[i], n2);
}


/*
 * delete a matrix allocated as above
 */

void delete_matrix(double** m)
{
  if(m == NULL) return;
  assert(*m);
  free(*m);
  assert(m);
  free(m);
}

/*
 * print an n x col matrix allocated as above out an opened outfile.
 * actually, this routine can print any double**
 */

void printMatrix(double **M, unsigned int n, unsigned int col, FILE *outfile)
{
  int i,j;
  assert(outfile);
  if(n > 0 && col > 0) assert(M);
  for(i=0; i<n; i++) {
    for(j=0; j<col; j++) {
#ifdef DEBUG
      if(j==col-1) MYprintf(outfile, "%.20f\n", M[i][j]);
      else MYprintf(outfile, "%.20f ", M[i][j]);
#else
      if(j==col-1) MYprintf(outfile, "%g\n", M[i][j]);
      else MYprintf(outfile, "%g ", M[i][j]);
#endif
    }
  }
}

/*
 * allocate and return an array of length n with scale*1 at
 * each entry
 */

double* ones(unsigned int n, double scale)
{
  double *o;
  unsigned int i;
  /* o = (double*) malloc(sizeof(double) * n); */
  o = new_vector(n);
  /* assert(o); */
  for(i=0; i<n; i++) o[i] = scale;
  return o;
}


/*
 * send a vector
 * of the matrix M out to a file
 */

void vector_to_file(const char* file_str, double* vector, unsigned int n)
{
  FILE* VOUT;
  unsigned int i;

  VOUT = fopen(file_str, "w");
  assert(VOUT);
  for(i=0; i<n; i++) MYprintf(VOUT, "%g\n", vector[i]);
  fclose(VOUT);
}


/*
 * open file with the given name
 * and print the passed matrix to it
 */

void matrix_to_file(const char* file_str, double** matrix, unsigned int n1, unsigned int n2)
{
  FILE* MOUT;

  MOUT = fopen(file_str, "w");
  assert(MOUT);
  printMatrix(matrix, n1, n2, MOUT);
  fclose(MOUT);
}

/*
 * allocates a new double array of size n1
 */

double* new_vector(unsigned int n)
{
  double *v;
  if(n == 0) return NULL;
  v = (double*) malloc(sizeof(double) * n);
  return v;
}


/*
 * allocates a new double array of size n1
 * and fills it with the contents of vold
 */

double* new_dup_vector(double* vold, unsigned int n)
{
  double *v;
  v = new_vector(n);
  dupv(v, vold, n);
  return v;
}

/*
 * zeros out v
 * (assumes that it has already been allocated)
 */

void zerov(double*v, unsigned int n)
{
  unsigned int i;
  for(i=0; i<n; i++) v[i] = 0;
}

/*
 * allocates a new double array of size n1
 * and fills it with zeros
 */

double* new_zero_vector(unsigned int n)
{
  double *v;
  v = new_vector(n);
  zerov(v, n);
  return v;
}

/*
 * copies vold to v
 * (assumes v has already been allcocated)
 */

void dupv(double *v, double* vold, unsigned int n)
{
  unsigned int i;
  for(i=0; i<n; i++) v[i] = vold[i];
}



/*
 * sumv:
 *
 * return the sum of the contents of the vector
 */

double sumv(double *v, unsigned int n)
{
  unsigned int i;
  double s;
  if(n==0) return 0;
  assert(v);
  s = 0;
  for(i=0; i<n; i++) s += v[i];
  return(s);
}

/*
 * printing a vector out to outfile
 */

void printVector(double *v, unsigned int n, FILE *outfile, PRINT_PREC type)
{
  unsigned int i;
  if(type==HUMAN) for(i=0; i<n; i++) MYprintf(outfile, "%g ", v[i]);
  else if(type==MACHINE) for(i=0; i<n; i++) MYprintf(outfile, "%.20f ", v[i]);
  else error("bad PRINT_PREC type");
  MYprintf(outfile, "\n");
}


/*
 * return the minimum element in the vector.
 * pass back the index of the minimum through
 * the which pointer
 */

double min(double *v, unsigned int n, unsigned int *which)
{
  unsigned int i;
  double min;

  *which = 0;
  min = v[0];

  for(i=1; i<n; i++) {
    if(v[i] < min)  {
      min = v[i];
      *which = i;
    }
  }

  return min;
}


/*
 * return the maximum element in the vector.  pass back the index of
 * the maximum through the which pointer
 */

double max(double *v, unsigned int n, unsigned int *which)
{
  unsigned int i;
  double max;

  *which = 0;
  max = v[0];

  for(i=1; i<n; i++) {
    if(v[i] > max)  {
      max = v[i];
      *which = i;
    }
  }

  return max;
}


/*
 * new vector of integers of length n
 */

int *new_ivector(unsigned int n)
{
  int *iv;
  if(n == 0) return NULL;
  iv = (int*)  malloc(sizeof(int) * n);
  assert(iv);
  return iv;
}

/*
 * sq:
 *
 * calculate the square of x
 */

double sq(double x)
{
  return x*x;
}
