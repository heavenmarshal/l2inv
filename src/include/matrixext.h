#ifndef __MATRIX_EXT_H__
#define __MATRIX_EXT_H__

void onev(double*, unsigned int);

double* new_one_vector(unsigned int);

void sum_scalar(double*, double*, double, unsigned int);

void sum_scalar_inplace(double*, double, unsigned int);

double* new_reduce_row(double**, unsigned int, unsigned int);

void vector_exp_check(double *, int);

void vector_log(double*, int);
#endif
