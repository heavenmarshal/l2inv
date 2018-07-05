/******************************************************************************** 
 *
 * Bayesian Regression and Adaptive Sampling with Gaussian Process Trees
 * Copyright (C) 2005, University of California
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Questions? Contact Robert B. Gramacy (rbgramacy@ams.ucsc.edu)
 *
 ********************************************************************************/


#include <stdlib.h>
#include <assert.h>
#include "linalg.h"
#include "matrix.h"
#include "rhelp.h"

#ifdef FORTPACK
char uplo = 'U';
#endif
/* #define DEBUG */

/*
 * linalg_ddot:
 *
 * analog of ddot in cblas nad blas
 */

double linalg_ddot(n, X, ldx, Y, ldy)
int n, ldx, ldy;
double *X, *Y;
{
  double result;

#ifdef FORTBLAS
  size_t n64,ldx64,ldy64;
  n64 = n; ldx64 = ldx; ldy64=ldy;
  result = ddot(&n64,X,&ldx64,Y,&ldy64);
#else
  result = cblas_ddot(n, X, ldx, Y, ldy);
#endif
  return result;
}

/*
 * linalg_dgemm:
 *
 * analog of dgemm in cblas nad blas
 * assumed column major representation
 */

void linalg_dgemm(TA, TB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
const enum CBLAS_TRANSPOSE TA, TB;
int m, n, k, lda, ldb, ldc;
double alpha, beta;
double **A, **B, **C;
{
#ifdef FORTBLAS
  size_t m64, n64, k64, lda64, ldb64, ldc64;
  char ta, tb;
  m64 = m; n64 = n; k64 = k; lda64 = lda; ldb64 = ldb; ldc64 = ldc;
  if(TA == CblasTrans) ta = 'T'; else ta = 'N';
  if(TB == CblasTrans) tb = 'T'; else tb = 'N';
  /* dgemm(&ta,&tb,&m,&n,&k,&alpha,*A,&lda,*B,&ldb,&beta,*C,&ldc); */
  dgemm(&ta,&tb,&m64,&n64,&k64,&alpha,*A,&lda64,*B,&ldb64,&beta,*C,&ldc64);
#else
  cblas_dgemm(CblasColMajor,TA,TB,m,n,k,alpha,*A,lda,*B,ldb,beta,*C,ldc);
#endif
}

/*
 * linalg_dgemv:
 *
 * analog of dgemv in cblas nad blas
 * assumed column major representation
 */

void linalg_dgemv(TA, m, n, alpha, A, lda, X, ldx, beta, Y, ldy)
const enum CBLAS_TRANSPOSE TA;
int m, n, lda, ldx, ldy;
double alpha, beta;
double **A;
double *X, *Y;
{
#ifdef FORTBLAS
  size_t m64, n64, lda64, ldx64, ldy64;
  char ta;
  m64 = m; n64 = n, lda64 = lda; ldx64 = ldx; ldy64 = ldy;
  if(TA == CblasTrans) ta = 'T'; else ta = 'N';
  /* dgemv(&ta,&m,&n,&alpha,*A,&lda,X,&ldx,&beta,Y,&ldy); */
  dgemv(&ta,&m64,&n64,&alpha,*A,&lda64,X,&ldx64,&beta,Y,&ldy64);
#else
  cblas_dgemv(CblasColMajor,TA,m,n,alpha,*A,lda,X,ldx,beta,Y,ldy);
#endif
}

/*
 * linalg_dsymm:
 *
 * analog of dsymm in cblas nad blas
 * assumed column major and upper-triangluar representation
 */


void linalg_dsymm(SIDE, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
const enum CBLAS_SIDE SIDE;
int m, n, lda, ldb, ldc;
double alpha, beta;
double **A, **B, **C;
{
#ifdef FORTBLAS
  size_t m64, n64, lda64, ldb64, ldc64;
  char side;
  m64 = m; n64 = n; lda64 = lda; ldb64 = ldb; ldc64 = ldc;
  if(SIDE == CblasRight) side = 'R'; else side = 'L';
  /* dsymm(&side,&uplo,&m,&n,&alpha,*A,&lda,*B,&ldb,&beta,*C,&ldc); */
  dsymm(&side,&uplo,&m64,&n64,&alpha,*A,&lda64,*B,&ldb64,&beta,*C,&ldc64);
#else
  cblas_dsymm(CblasColMajor,SIDE,CblasUpper,m,n,alpha,*A,lda,*B,ldb,beta,*C,ldc);
#endif
}

/*
 * linalg_dsymv:
 *
 * analog of dsymv in cblas and blas
 * assumed column major representation
 */


void linalg_dsymv(n, alpha, A, lda, X, ldx, beta, Y, ldy)
int n, lda, ldx, ldy;
double alpha, beta;
double **A;
double *X, *Y;
{
#ifdef FORTBLAS
  size_t n64, lda64, ldy64, ldx64;
  n64 = n; lda64 = lda; ldx64 = ldx; ldy64 = ldy;
  /* dsymv(&uplo,&n,&alpha,*A,&lda,X,&ldx,&beta,Y,&ldy); */

  dsymv(&uplo,&n64,&alpha,*A,&lda64,X,&ldx64,&beta,Y,&ldy64);
#else
  cblas_dsymv(CblasColMajor,CblasUpper,n,alpha,*A,lda,X,ldx,beta,Y,ldy);
#endif
}

/*
 * linalg_dposv:
 *
 * analog of dposv in clapack and lapack where 
 * Mutil is with colmajor and uppertri or rowmajor
 * and lowertri
 */

int linalg_dposv(n, Mutil, Mi)
int n;
double **Mutil, **Mi;
{
  long info;
	
  /* then use LAPACK */
#ifdef FORTPACK
  size_t n64;
  n64 = n;
  dposv(&uplo,&n64,&n64,*Mutil,&n64,*Mi,&n64,&info);
#else
  /*info = clapack_dposv(CblasColMajor,CblasUpper,n,n,*Mutil,n,*Mi,n);*/
  info = clapack_dposv(CblasRowMajor,CblasLower,n,n,*Mutil,n,*Mi,n);
#endif
  
#ifdef DEBUG
  if(info != 0) {
    matrix_to_file("M.dump", Mutil, n, n);
    error("offending matrix dumped into matrix.dump");
  }
#endif
  
  return (int) info;
}
