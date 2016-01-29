#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

  double interp_( double * yarray,
                           double * xarray,
                           double *  val_interp,
                           int    *  size ) 
{
  gsl_interp * intp = gsl_interp_alloc(gsl_interp_akima, *size);
  double ret;

  (void) gsl_interp_init(intp,xarray,yarray,*size);
  ret = gsl_interp_eval(intp, xarray,yarray,*val_interp,NULL);

  gsl_interp_free(intp);

  return ret;

}

