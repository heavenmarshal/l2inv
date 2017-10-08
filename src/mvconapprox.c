#include<Rmath.h>
#include"newtonsolver.h"
#include"mvconapprox.h"
#include"nleqslv.h"

#define SQRPI2 2.50662827463
#define SGN(x) ((double)((0.0<x)-(0.0>x)))

double dkappaSeq(double x, void *param)
{
  parDkaps* ppar = (parDkaps*) param;
  int p = ppar -> p;
  const double *sig2 = ppar -> sig2;
  const double *mu2 = ppar -> mu2;
  double barval = ppar -> barval;
  double denom, dkappa = 0.0;
  int i;
  for(i=0; i!=p; ++i)
  {
    denom = 1.0 - 2.0*sig2[i]*x;
    dkappa += sig2[i]/denom;
    dkappa += 4.0*sig2[i]*mu2[i]*(1.0-sig2[i]*x)*x/denom/denom;
    dkappa += mu2[i];
  }
  return dkappa - barval;
}
double logdkappaSeq(double x, void *param)
{
  parDkaps* ppar = (parDkaps*) param;
  double dkapseq = dkappaSeq(x,param);
  double barval = ppar -> barval;
  double logdkap = log(dkapseq+barval)-log(barval);
  return logdkap;
}
double dkappa2(double x, void* param)
{
  parDkaps* ppar = (parDkaps*) param;
  int p = ppar -> p;
  const double *sig2 = ppar -> sig2;
  const double *mu2 = ppar -> mu2;
  double dkappa2 = 0.0;
  double denom, denom2;
  int i;
  for(i=0; i!=p; ++i)
  {
    denom = 1.0 - 2.0*sig2[i]*x;
    denom2 = denom*denom;
    dkappa2 += 2.0*sig2[i]*sig2[i]/denom2;
    dkappa2 += 4.0*sig2[i]*mu2[i]/denom2/denom;
  }
  return dkappa2;
}
double logdkappa2(double x, void* param)
{
  parDkaps* ppar = (parDkaps*) param;
  double dkaseq = dkappaSeq(x,param);
  double barval = ppar -> barval;
  double dkap2 = dkappa2(x,param);
  double dkap = dkaseq+barval;
  return dkap2/dkap;
}
void dkappadd(double x, void *param,
	      double* fv, double* dfv)
{
  parDkaps* ppar = (parDkaps*) param;
  int p = ppar -> p;
  const double *sig2 = ppar -> sig2;
  const double *mu2 = ppar -> mu2;
  double barval = ppar -> barval;
  double denom, denom2, sigfmu;
  double dkappa = 0.0, dkappa2 = 0.0;
  int i;
  for(i = 0; i != p; ++i)
  {
    denom = 1.0 - 2.0*sig2[i]*x;
    denom2 = denom*denom;
    sigfmu = 4.0*sig2[i]*mu2[i];
    dkappa += sig2[i]/denom;
    dkappa += sigfmu * (1.0-sig2[i]*x)*x/denom2;
    dkappa += mu2[i];
    dkappa2 += 2.0*sig2[i]*sig2[i]/denom2;
    dkappa2 += sigfmu/denom2/denom;
  }
  *fv = dkappa - barval;
  *dfv = dkappa2;
}
void logdkappadd(double x, void *param,
		 double *fv, double *dfv)
{
  double dkapseq, dkap, dkap2;
  double barval;
  parDkaps* ppar = (parDkaps*) param;
  dkappadd(x,param,&dkapseq,&dkap2);
  barval = ppar -> barval;
  dkap = dkapseq + barval;
  *fv = log(dkap)-log(barval);
  *dfv = dkap2/dkap;
}
double transfun(double x, double upb)
{
  double y = x>0? upb - exp(-x): x+upb-1.0;
  return y;
}
double dtransfun(double x)
{
  double y = x>0? exp(-x): 1.0;
}
double transdkappaSeq(double x, void *param)
{
  int i, p;
  parDkaps* ppar = (parDkaps*) param;
  p = ppar -> p;
  const double *sig2 = ppar -> sig2;
  const double *mu2 = ppar -> mu2;
  double barval = ppar -> barval;
  double upb = ppar -> upb;
  double xt = transfun(x,upb);
  double denom, dkappa = 0.0;
  for(i=0; i!=p; ++i)
  {
    denom = 1.0 - 2.0*sig2[i]*xt;
    dkappa += sig2[i]/denom;
    dkappa += 4.0*sig2[i]*mu2[i]*(1.0-sig2[i]*xt)*xt/denom/denom;
    dkappa += mu2[i];
  }
  return dkappa - barval;
}
double transdkappa2(double x, void* param)
{
  int i, p;
  parDkaps* ppar = (parDkaps*) param;
  p = ppar -> p;
  const double *sig2 = ppar -> sig2;
  const double *mu2 = ppar -> mu2;
  double upb = ppar -> upb;
  double dkappa2 = 0.0;
  double xt = transfun(x,upb);
  double dxt = dtransfun(x);
  double denom, denom2;
  for(i=0; i!=p; ++i)
  {
    denom = 1.0 - 2.0*sig2[i]*xt;
    denom2 = denom * denom;
    dkappa2 += 2.0*sig2[i]*sig2[i]/denom2;
    dkappa2 += 4.0*sig2[i]*mu2[i]/denom2/denom;
  }
  dkappa2 *= dxt;
  return dkappa2;
}
/* void transdkappadd(double x, void *param, */
/* 		   double* fv, double* dfv) */
/* { */
/*   int i, p; */
/*   parDkaps* ppar = (parDkaps*) param; */
/*   p = ppar -> p; */
/*   const double *sig2 = ppar -> sig2; */
/*   const double *mu2 = ppar -> mu2; */
/*   double barval = ppar -> barval; */
/*   double upb = ppar -> upb; */
/*   double xt = transfun(x,upb); */
/*   double dxt = dtransfun(x); */
/*   double denom, denom2, sigfmu; */
/*   double dkappa = 0.0, dkappa2 = 0.0; */
/*   for(i = 0; i !=p; ++i) */
/*   { */
/*     denom = 1.0 - 2.0*sig2[i]*xt; */
/*     denom2 = denom*denom; */
/*     sigfmu = 4.0*sig2[i]*mu2[i]; */
/*     dkappa += sig2[i]/denom; */
/*     dkappa += sigfmu * (1.0-sig2[i]*xt)*xt/denom2; */
/*     dkappa += mu2[i]; */
/*     dkappa2 += 2.0*sig2[i]*sig2[i]/denom2; */
/*     dkappa2 += sigfmu/denom2/denom; */
/*   } */
/*   *fv = dkappa - barval; */
/*   *dfv = dkappa2*dxt; */
/* } */
void kappafs(double x, const double* sig2,
	     const double* mu2, int p,
	     double* kappa, double* kappa2,
	     double* kappa3)
{
  double k=0.0, k2=0.0, k3=0.0;
  double denom, denom2, denom3, denom4;
  int i;
  for(i = 0; i !=p; ++i)
  {
    denom = 1.0 - 2.0 * sig2[i] * x;
    denom2 = denom * denom;
    denom3 = denom2 *denom;
    denom4 = denom2 * denom2;
    k -= 0.5 * log(denom);
    k += 2.0 * sig2[i]*mu2[i]*x*x/denom + mu2[i]*x;
    k2 += 2.0*sig2[i]*sig2[i]/denom2;
    k2 += 4.0*sig2[i]*mu2[i]/denom3;
    k3 += 8.0*sig2[i]*sig2[i]*sig2[i]/denom3;
    k3 += 24.0*sig2[i]*sig2[i]*mu2[i]/denom4;
  }
  *kappa = k;
  *kappa2 = k2;
  *kappa3 = k3;
}

double posapprox(double tval, double zz, double ww,
		 double dk2, double lambda)
{
  double zz2 = zz*zz;
  double eww = exp(-0.5*ww*ww);
  double pnzz = 1 - pnorm(zz,0.0,1.0,1,0);
  double ezz = exp(0.5*zz2);
  double sqrdk2 = sqrt(dk2);
  double c1 = eww*(sqrdk2/SQRPI2-tval*dk2*ezz*pnzz);
  double c2 = ezz*eww*sqrdk2*lambda/6.0*(pnzz*zz2*(zz2+3.0)-dnorm(zz,0.0,1.0,0)*zz*(zz2+2));
  return c1 + c2;
}

double negapprox(double tval, double zz, double ww,
		 double dk2, double lambda, double mumk)
{
  double zz2 = zz*zz;
  double eww = exp(-0.5*ww*ww);
  double pzz = pnorm(zz, 0.0, 1.0, 1, 0);
  double ezz = exp(0.5 * zz2);
  double sqrdk2 = sqrt(dk2);
  double c1 = mumk + eww*(sqrdk2/SQRPI2+tval*dk2*ezz*pzz);
  double c2 = ezz*eww*sqrdk2*lambda/6.0*
    (pzz*zz2*(zz2+3)+dnorm(zz,0.0,1.0,0)*zz*(zz2+2));
  return c1 - c2;
}

void saddleapprox(const double* sig2m, const double *mu2m,
		  const double *barval, const double *mumk,
		  const double *upb, int n, int p, double *info)
{
  double tval, kp, dk2, dk3, sgnt;
  double zz, lambda, ww;
  int i, j, stat;
  for(i = 0, j = 0; i != n; ++i, j+=p)
  {
    parDkaps param = {p, sig2m+j, mu2m+j, barval[i], upb[i]};

    /* stat = newtonsolver(0.5*upb[i], &dkappaSeq, &dkappa2, */
    /* 			&dkappadd, (void*) &param, &tval, */
    /* 			100, 1E-8, 1E-10); */

    /* if(stat != success) */
    /* { */
      /* call more expensive solution */
    stat = nleqslv(0.0, &transdkappaSeq, &transdkappa2,
		   (void*) &param, &tval, 100, 1E-8, 1E-8);
    if(stat != success)
    {
      info[i] = NAN; 		// na
      continue;
    }
    tval = transfun(tval,upb[i]);
    /* } */
    kappafs(tval,&sig2m[j],&mu2m[j],p,
	    &kp, &dk2, &dk3);
    sgnt = SGN(tval);
    zz = tval * sqrt(dk2);
    lambda = dk3/pow(dk2,1.5);
    ww = sgnt * sqrt(2*barval[i]*tval-2.0*kp);
    info[i] = tval>0.0? posapprox(tval,zz,ww,dk2,lambda):
      negapprox(tval,zz,ww,dk2,lambda,mumk[i]);
  }
}
void saddleapprox_R(const double *sig2_, const double *mu2_,
		    const double *barval_, const double *mumk_,
		    const double *upb_, const int* n_,
		    const int* p_, double *info_)
{
  saddleapprox(sig2_, mu2_, barval_, mumk_, upb_, *n_, *p_, info_);
}
