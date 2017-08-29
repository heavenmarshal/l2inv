#include<R.h>
#include<Rmath.h>
#include"matrix.h"
#define SQ(x) ((x)*(x))
void mvcon(double** coefmean, double** coefsd,
	   double* barrier, double* chatstar,
	   double* dsq, int nfea, int nsim,
	   int numbas, double* criterion)
{
  int i, j, k;
  double sum, fv, u;
  GetRNGstate();
  for(i = 0; i < nfea; ++i)
  {
    sum = 0.0;
    for(j = 0; j < nsim; ++j)
    {
      fv = barrier[i];
      for(k = 0; k < numbas; ++k)
      {
	u   = norm_rand();
	fv -= dsq[k]*SQ(u*coefsd[i][k]+coefmean[i][k]-chatstar[i]);
      }
      if(fv>0)
	sum += fv;
    }
    criterion[i] = sum/nsim;
  }
  PutRNGstate();
}

void mvcon_R(double* coefmean_, double* coefsd_,
	     double* barrier_, double* chatstar_,
	     double* dsq_, int* nfea_, int* nsim_,
	     int* numbas_, double* criterion_)
{
  int nfea = *nfea_, numbas = *numbas_;
  double** coefmean = new_matrix_bones(coefmean_,nfea,numbas);
  double** coefsd = new_matrix_bones(coefsd_,nfea,numbas);
  mvcon(coefmean,coefsd,barrier_,chatstar_,dsq_,
	nfea,*nsim_,numbas,criterion_);
  free(coefmean);
  free(coefsd);
}
