#include <stdlib.h>
#include "matrix.h"
#include "gp_sep.h"
#include "gpseplm.h"

extern unsigned int  NGPsep;
extern GPsep **gpseps;

extern unsigned int NGPsepLm;
extern GPsepLm **gplms;

void esl2dZmean(unsigned int *nbas, unsigned int *gpidces,
		double *x, double *reddsq, double *cht,
		double *esl2d)
{
  int i, p;
  double mean, sig2, df;
  GPsep* gpsep;
  p = *nbas;
  *esl2d=0.0;
  for(i=0; i < p; ++i)
  {
    gpsep = gpseps[gpidces[i]];
    predGPsep_lite(gpsep, 1, &x, &mean, &sig2, &df, NULL);
    *esl2d += reddsq[i]*sq(mean-cht[i])+reddsq[i]*sig2;
  }
}
void esl2dCmean(unsigned int *nbas, unsigned int *gpidces,
		double *x, double *reddsq, double *cht,
		double *esl2d)
{
  int i, p;
  double mean, sig2, df;
  double *hh;
  GPsepLm* gplm;
  hh = new_vector(1);
  *hh = 1.0;
  p = *nbas;
  *esl2d = 0.0;
  for(i=0; i < p; ++i)
  {
    gplm = gplms[gpidces[i]];
    predGPsepLm_lite(gplm,1,&x,&hh,&mean,&sig2,&df,NULL);
    *esl2d += reddsq[i]*sq(mean-cht[i])+reddsq[i]*sig2;
  }
  free(hh);
}

void esl2dLmean(unsigned int *nbas, unsigned int *ndim, unsigned int *gpidces,
		double *x, double *reddsq, double *cht, double *esl2d)
{
  int i, p;
  double mean, sig2, df;
  double *hh;
  GPsepLm* gplm;
  hh = new_vector(*ndim+1);
  hh[0] = 1.0;
  dupv(hh+1,x,*ndim);
  p = *nbas;
  *esl2d = 0.0;
  for(i=0; i < p; ++i)
  {
    gplm = gplms[gpidces[i]];
    predGPsepLm_lite(gplm,1,&x,&hh,&mean,&sig2,&df,NULL);
    *esl2d += reddsq[i]*sq(mean-cht[i])+reddsq[i]*sig2;
  }
  free(hh);
}
