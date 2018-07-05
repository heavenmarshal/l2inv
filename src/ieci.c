/****************************************************************************
 *
 * Local Approximate Gaussian Process Regression
 * Copyright (C) 2013, The University of Chicago
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301  USA
 *
 * Questions? Contact Robert B. Gramacy (rbg@vt.edu)
 *
 ****************************************************************************/


#include "matrix.h"
#include "linalg.h"
#include "rhelp.h"
#include "ieci.h"
#include "util.h"
#include <stdlib.h>
#include <assert.h>
#ifdef RPRINT
#include <R.h>
#include <Rmath.h>
#else
#include <math.h>
#endif



/*
 * calc_alc:
 *
 * function that iterates over the m Xref locations, and the
 * stats calculated by previous calc_* function in order to
 * calculate the reduction in variance
 */

double calc_alc(const int m, double *ktKik, double *s2p, const double phi,
		double *badj, const double tdf, double *w)
{
  int i;
  double zphi, ts2, alc, dfrat;

  dfrat = tdf/(tdf - 2.0);
  alc = 0.0;
  for(i=0; i<m; i++) {
    zphi = (s2p[1] + phi)*ktKik[i];
    if(badj) ts2 = badj[i] * zphi / (s2p[0] + tdf);
    else ts2 = zphi / (s2p[0] + tdf);
    if(w) alc += w[i]*dfrat*ts2;
    else alc += ts2*dfrat;
  }

  return (alc/m);
}


/*
 * calc_ktKikx:
 *
 * function for calculating the ktKikx vector used in the
 * IECI calculation -- writes over the KtKik input --
 * R interface (calc_ktKikx_R) available in plgp source tree
 */

void calc_ktKikx(double *ktKik, const int m, double **k, const int n,
		 double *g, const double mui, double *kxy, double **Gmui_util,
		 double *ktGmui_util, double *ktKikx)
{
  int i;
  double **Gmui;
  double *ktGmui;

  /* first calculate Gmui = g %*% t(g)/mu */
  if(!Gmui_util) Gmui = new_matrix(n, n);
  else Gmui = Gmui_util;
  linalg_dgemm(CblasNoTrans,CblasTrans,n,n,1,
               mui,&g,n,&g,n,0.0,Gmui,n);

  /* used in the for loop below */
  if(!ktGmui_util) ktGmui = new_vector(n);
  else ktGmui = ktGmui_util;

  /* loop over all of the m candidates */
  for(i=0; i<m; i++) {

    /* ktGmui = drop(t(k) %*% Gmui) */
    /* zerov(ktGmui, n); */
    linalg_dsymv(n,1.0,Gmui,n,k[i],1,0.0,ktGmui,1);

    /* ktKik += diag(t(k) %*% (g %*% t(g) * mui) %*% k) */
    if(ktKik) ktKikx[i] = ktKik[i] + linalg_ddot(n, ktGmui, 1, k[i], 1);
    else ktKikx[i] = linalg_ddot(n, ktGmui, 1, k[i], 1);

    /* ktKik.x += + 2*diag(kxy %*% t(g) %*% k) */
    ktKikx[i] += 2.0*linalg_ddot(n, k[i], 1, g, 1)*kxy[i];

    /* ktKik.x + kxy^2/mui */
    ktKikx[i] += sq(kxy[i])/mui;
  }

  /* clean up */
  if(!ktGmui_util) free(ktGmui);
  if(!Gmui_util) delete_matrix(Gmui);
}
