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
#include "gp.h"
#include "linalg.h"
#include "rhelp.h"
/* 
 * pred_generic:
 *
 * a function that captures the predictive equations without reference
 * struct gp objects.  Created so that code can be shared between GP
 * and GPsep objects, and beyond
 */

void pred_generic(const unsigned int n, const double phidf, double *Z, 
  double **Ki, const unsigned int nn, double **k, double *mean, double **Sigma)
{
  int i, j;
  double **ktKi, **ktKik;

  /* ktKi <- t(k) %*% util$Ki */
  ktKi = new_matrix(n, nn);
  linalg_dsymm(CblasRight,nn,n,1.0,Ki,n,k,nn,0.0,ktKi,nn);
  /* ktKik <- ktKi %*% k */
  ktKik = new_matrix(nn, nn);
  linalg_dgemm(CblasNoTrans,CblasTrans,nn,nn,n,
               1.0,k,nn,ktKi,nn,0.0,ktKik,nn);

  /* mean <- ktKi %*% Z */
  linalg_dgemv(CblasNoTrans,nn,n,1.0,ktKi,nn,Z,1,0.0,mean,1);

  /* Sigma <- phi*(Sigma - ktKik)/df */
  for(i=0; i<nn; i++) {
    Sigma[i][i] = phidf * (Sigma[i][i] - ktKik[i][i]);
    for(j=0; j<i; j++)
      Sigma[j][i] = Sigma[i][j] = phidf * (Sigma[i][j] - ktKik[i][j]);
  }

  /* clean up */
  delete_matrix(ktKi);
  delete_matrix(ktKik);
}

void new_predutil_generic_lite(const unsigned int n, double ** Ki, 
  const unsigned int nn, double **k, double ***ktKi, double **ktKik)
{
  unsigned int i, j;

  /* ktKi <- t(k) %*% util$Ki */
  *ktKi = new_matrix(n, nn);
  linalg_dsymm(CblasRight,nn,n,1.0,Ki,n,k,nn,0.0,*ktKi,nn);
  /* ktKik <- diag(ktKi %*% k) */
  *ktKik = new_zero_vector(nn); 
  for(i=0; i<nn; i++) for(j=0; j<n; j++) (*ktKik)[i] += (*ktKi)[j][i]*k[j][i];
}
