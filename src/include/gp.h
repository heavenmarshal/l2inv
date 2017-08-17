#ifndef __GP_H__
#define __GP_H__
void pred_generic(const unsigned int n, const double phidf, double *Z, 
  double **Ki, const unsigned int nn, double **k, double *mean, 
  double **Sigma);
void new_predutil_generic_lite(const unsigned int n, double **Ki, 
  const unsigned int nn, double **k, double ***ktKi, double **ktKik);
#endif
