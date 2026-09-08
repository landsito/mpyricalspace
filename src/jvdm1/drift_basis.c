/*
 * drift_basis.h
 * Patrick Alken
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_bspline.h>

#include "drift_basis.h"

size_t n_drift_nt = N_DRIFT_SUM_DEFAULT_T;
size_t n_drift_ns = N_DRIFT_SUM_DEFAULT_S;
size_t n_drift_ne = N_DRIFT_SUM_DEFAULT_EUVAC;

drift_basis_workspace *
drift_basis_alloc(const size_t N[])
{
  drift_basis_workspace *w;

  w = (drift_basis_workspace *) calloc(1, sizeof(drift_basis_workspace));
  if (w == 0)
    {
      perror("malloc");
      return 0;
    }

  if (N[IDX_DRIFT_T] >= 4)
    {
      const size_t nbreak_t = N[IDX_DRIFT_T] - 4 + 2;

      w->bspline_workspace_t = gsl_bspline_alloc(4, nbreak_t);
      if (w->bspline_workspace_t == 0)
        {
          drift_basis_free(w);
          fprintf(stderr, "drift_basis_alloc: malloc failed\n");
          return 0;
        }

      gsl_bspline_knots_uniform(DRIFT_LT_MIN,
                                DRIFT_LT_MAX,
                                w->bspline_workspace_t);
    }

  return (w);
} /* drift_basis_alloc() */

void
drift_basis_free(drift_basis_workspace *w)
{
  if (!w)
    return;

  if (w->bspline_workspace_t)
    gsl_bspline_free(w->bspline_workspace_t);

  free(w);
} /* drift_basis_free() */

/*
get_drift_basis_t()
  Return basis function for local time variable

u_{i}(t) = cos(i * t)

Inputs: t - local time in hours
*/

int
get_drift_basis_t(double t, double y[], void *params)

{
  drift_basis_workspace *w = (drift_basis_workspace *) params;
  gsl_vector_view v = gsl_vector_view_array(y, N_DRIFT_SUM_T);

  if (t < DRIFT_LT_MIN || t > DRIFT_LT_MAX)
    return GSL_FAILURE;

  gsl_bspline_eval(t, &v.vector, w->bspline_workspace_t);

  return GSL_SUCCESS;
} /* get_drift_basis_t() */

/*
get_drift_basis_s()
  Return basis function for season variable

u_{i}(s) = { cos(i*s/2),     i even
             sin((i+1)*s/2), i odd  }

with s \in [0, 365]
*/

int
get_drift_basis_s(double s, double y[], void *params)

{
  size_t i;

  for (i = 0; i < N_DRIFT_SUM_S; ++i)
    {
      if ((i % 2) == 0)
        y[i] = cos(i * M_PI * s / 365.25);
      else
        y[i] = sin((i + 1) * M_PI * s / 365.25);
    }

  return GSL_SUCCESS;
} /* get_drift_basis_s() */

int
get_drift_basis_euvac(double euvac, double y[], void *params)
{
  size_t i;

  if (euvac < DRIFT_EUVAC_MIN || euvac > DRIFT_EUVAC_MAX)
    return GSL_FAILURE;

  for (i = 0; i < N_DRIFT_SUM_EUVAC; ++i)
    y[i] = pow(euvac, (double) i);

  return GSL_SUCCESS;
} /* get_drift_basis_euvac() */
