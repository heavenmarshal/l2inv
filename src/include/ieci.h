#ifndef __IECI_H__
#define __IECI_H__

typedef enum STYPE {UL=3001, MEAN=3002, NORM=3004} Stype;

double calc_alc(const int m, double *ktKik, double *s2p, const double phi,
		double *badj, const double tdf, double *w);
void calc_ktKikx(double *ktKik, const int m, double **k, const int n,
		 double *g, const double mui, double *kxy, double **Gmui_util,
		 double *ktGmui_util, double *ktKikx);
#endif

