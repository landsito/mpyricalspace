/*
 * drift_calc.h
 * Patrick Alken
 */

#ifndef INCLUDED_drift_calc_h
#define INCLUDED_drift_calc_h

#include <gsl/gsl_vector.h>
#include "ndlinear/gsl_multifit_ndlinear.h"
#include "drift_basis.h"

typedef struct
{
  gsl_vector *mean_coeffs;
  gsl_vector *stddev_coeffs;
  gsl_multifit_ndlinear_workspace *ndlinear_p;
  drift_basis_workspace *drift_basis_workspace_p;

  int (*u[N_DRIFT_FIT])(double x, double y[], void *p);
} drift_calc_workspace;

/* location of JULIA radar is lat: 11.95 S, lon: 76.87 W */
#define DRIFT_LONGITUDE    (-76.87)
#define DRIFT_LATITUDE     (-11.95)

/*
 * Prototypes
 */

drift_calc_workspace *drift_calc_alloc(void);
void drift_calc_free(drift_calc_workspace *w);
double drift_calc_mean(double t, double s, double f,
                       drift_calc_workspace *w);
double drift_calc_stddev(double t, double s, double f,
                         drift_calc_workspace *w);

#endif /* INCLUDED_drift_calc_h */
