#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

gsl_interp *intp;
#pragma omp threadprivate(intp)

  double interp_( double * yarray,
                           double * xarray,
                           double *  val_interp,
                           int    *  size) 
/* take in intp so it is only allocated/deallocated */
/* once per time step */
{
  double ret;

  (void) gsl_interp_init(intp,xarray,yarray,*size);
  ret = gsl_interp_eval(intp, xarray,yarray,*val_interp,NULL);

  return ret;

}

void interp_alloc_(int *size)
{
#pragma omp critical
   { 
   intp = gsl_interp_alloc(gsl_interp_akima, *size); 
   }
}

void interp_free_(void)
{
#pragma omp critical
   {
   gsl_interp_free(intp);
   }
}
