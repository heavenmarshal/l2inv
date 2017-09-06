#ifndef MVCONAPPROX_H
#define MVCONAPPROX_H

typedef struct
{
  int p;
  const double *sig2;
  const double *mu2;
  double barval;
  double upb;
} parDkaps;
double dkappaSeq(double, void*);
double dkappa2(double, void*);
void dkappadd(double, void*, double*, double*);
void kappafs(double, const double*, const double*,
	     int, double*, double*, double*);
double posapprox(double, double, double, double, double);
double negapprox(double, double, double, double, double, double);
void saddleapprox(const double*, const double*,
		  const double*, const double*,
		  const double*, int, int, double*);
#endif
