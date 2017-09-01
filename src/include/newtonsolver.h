#ifndef NEWTONSOLVER_H
#define NEWTONSOLVER_H

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
