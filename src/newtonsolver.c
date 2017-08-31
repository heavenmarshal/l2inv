#include"newtonsolver.h"

int newtonsolver(double x0, targfun fun, targderv dfun, targfunderv fdf,
		 void *param, double* root, int maxiter,
		 double xabs, double xrel)
{
  gsl_function_fdf FDF;
  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;
  double xp;
  int iter, status;
  FDF.f = fun;
  FDF.df = dfun;
  FDF.fdf = fdf;
  FDF.params = param;
  T = gsl_root_fdfsolver_newton;
  s = gsl_root_fdfsolver_alloc (T);
  gsl_root_fdfsolver_set (s, &FDF, x0);

  status = GSL_CONTINUE;
  for(iter = 0; iter < maxiter && status==GSL_CONTINUE; ++iter)
  {
    status = gsl_root_fdfsolver_iterate (s);
    if(status != GSL_SUCCESS)
    {
      *root = x0;
      return failure;
    }
    xp = x0;
    x0 = gsl_root_fdfsolver_root (s);
    status = gsl_root_test_delta (x0, xp, xabs, xrel);
  }
  gsl_root_fdfsolver_free (s);
  *root = x0;
  if(status == GSL_CONTINUE)
    return exceed;
  return success;
}
