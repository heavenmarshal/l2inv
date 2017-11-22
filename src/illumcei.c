#include<R.h>
#include<Rmath.h>
#include "matrix.h"
#include "linalg.h"

void illumcei(double **coeffmean, double **coeffsd,
	      double* basis, double *xi, double barval,
	      double sdhat, int nmc, int nfea, int tlen,
	      int numbas, double* ei)
{
    int i, j, k;
    double u, info, sum, *cmean, *csd;
    double *cvec, *dev;
    cvec = new_vector(numbas);
    dev = new_vector(tlen);
    GetRNGstate();
    for(i = 0; i < nfea; i++)
    {
	cmean = coeffmean[i];
	csd = coeffsd[i];
	sum = 0.0;
	for(j = 0; j < nmc; ++j)
	{
	    for(k = 0; k < numbas; ++k)
	    {
		u = norm_rand();
		cvec[k] = u*csd[k]+cmean[k];
	    }
	    for(k = 0; k < tlen; ++k)
	    {
		u = norm_rand();
		dev[k] = u*sdhat - xi[k];
	    }
	    linalg_dgemv(CblasNoTrans,tlen,numbas,1.0,&basis,tlen,cvec,1,1.0,dev,1);
	    info = barval - linalg_ddot(tlen,dev,1,dev,1);
	    if(info<0.0) info = 0.0;
	    sum += info;
	}
	ei[i] = sum/nmc;
    }
    free(cvec);
    free(dev);
    PutRNGstate();
}
void illumcei_R(double *coeffmean_, double *coeffsd_,
		double *basis_, double* xi_, double *barval_,
		double *sdhat_, int *nmc_, int *nfea_, int *tlen_,
		int *numbas_, double *ei_)
{
  int nfea = *nfea_, numbas = *numbas_;
  double **coeffmean, **coeffsd;
  coeffmean = new_matrix_bones(coeffmean_, nfea, numbas);
  coeffsd = new_matrix_bones(coeffsd_, nfea, numbas);
  illumcei(coeffmean, coeffsd, basis_, xi_, *barval_,
	   *sdhat_, *nmc_, nfea, *tlen_, numbas, ei_);
}
