#include "linalg.h"
#include "linalgext.h"

#ifdef FORTPACK
extern char uplo;
#endif

int linalgext_dposv(n, m, Mutil, Mi)
  int n, m;
  double **Mutil, **Mi;
{
  long info;
	
  /* then use LAPACK */
#ifdef FORTPACK
  size_t n64;
  size_t m64;
  n64 = n;
  m64 = m;
  dposv(&uplo,&n64,&m64,*Mutil,&n64,*Mi,&n64,&info);
#else
  /*info = clapack_dposv(CblasColMajor,CblasUpper,n,n,*Mutil,n,*Mi,n);*/
  info = clapack_dposv(CblasRowMajor,CblasLower,n,m,*Mutil,n,*Mi,n);
#endif
  
  return (int) info;
}
