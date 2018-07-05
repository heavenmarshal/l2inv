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


#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <stdio.h>

typedef enum PRINT_PREC {HUMAN=1001, MACHINE=1002} PRINT_PREC;

void zero(double **M, unsigned int n1, unsigned int n2);
void id(double **M, unsigned int n);
double ** new_id_matrix(unsigned int n);
double ** new_zero_matrix(unsigned int n1, unsigned int n2);
double ** new_matrix(unsigned int m, unsigned int n);
double ** new_matrix_bones(double *v, unsigned int n1, unsigned int n2);
double ** new_dup_matrix(double** M, unsigned int n1, unsigned int n2);
void dup_matrix(double** M1, double **M2, unsigned int n1, unsigned int n2);
void delete_matrix(double** m);
void printMatrix(double **M, unsigned int n, unsigned int col, FILE *outfile);
double* ones(unsigned int n, double scale);
void vector_to_file(const char* file_str, double *quantiles, unsigned int n);
void matrix_to_file(const char* file_str, double** matrix, unsigned int n1, unsigned int n2);
double* new_vector(unsigned int n);
double* new_dup_vector(double* vold, unsigned int n);
void zerov(double*v, unsigned int n);
double* new_zero_vector(unsigned int n);
void dupv(double *v, double* vold, unsigned int n);
double sumv(double *v, unsigned int n);
void printVector(double *v, unsigned int n, FILE *outfile, PRINT_PREC type);
double max(double *v, unsigned int n, unsigned int *which);
double min(double *v, unsigned int n, unsigned int *which);
int *new_ivector(unsigned int n);
double sq(double x);
#endif

