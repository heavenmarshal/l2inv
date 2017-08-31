#ifndef NEWTONSOLVER_H
#define NEWTONSOLVER_H
#include<gsl/gsl_errno.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_roots.h>

typedef double (*targfun) (double, void*);
typedef double (*targderv) (double, void*);
typedef void (*targfunderv) (double, void*, double*, double*);

enum Stat
{
  success,
  failure,
  exceed
};
int newtonsolver(double, targfun, targderv, targfunderv, void*,
		 double*, int, double, double);

#endif
