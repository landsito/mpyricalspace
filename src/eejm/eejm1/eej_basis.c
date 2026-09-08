/* eej_basis.c
 * 
 * Copyright (C) 2006, 2007 Patrick Alken
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
 * This file provides the basis functions used for each parameter
 * in the statistical model of the EEJ.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_bspline.h>

#include "eej_basis.h"

/*
eej_basis_alloc()
  Allocate a eej_basis workspace

Inputs: N        - array containing number of terms in each sum
        satidx   - index of satellite (EEJ_BASIS_satname)
        fit_type - type of fit (mean/stddev) (EEJ_BASIS_type)

Return: pointer to workspace
*/

eej_basis_workspace *
eej_basis_alloc(const size_t N[], size_t satidx, size_t fit_type)
{
  eej_basis_workspace *w;
  size_t nbreak;

  w = (eej_basis_workspace *) calloc(1, sizeof(eej_basis_workspace));
  if (w == 0)
    {
      fprintf(stderr, "eej_basis_alloc: malloc failed\n");
      return 0;
    }

  w->satidx = satidx;
  w->fit_type = fit_type;

  if (N[EEJ_IDX_T] < 4)
    return (w); /* don't initialize bspline stuff */

#ifdef EEJ_BSPLINE_FORCE_ZERO
  if (fit_type == EEJ_BASIS_MEAN)
    {
      /*
       * To make sure the mean goes to 0 at the endpoints, add
       * (k-1) splines to both ends but they won't be used in
       * the fit, since our basis function will ignore them.
       * Therefore the coeffs of these are always 0, guaranteeing
       * the mean drops to 0.
       *
       * nbreak = n - k + 2 + 2(k-1) = n + k
       */
      nbreak = N[EEJ_IDX_T] + EEJ_BSPLINE_ORDER_T;
    }
  else
#endif /* EEJ_BSPLINE_FORCE_ZERO */
    nbreak = N[EEJ_IDX_T] - EEJ_BSPLINE_ORDER_T + 2;

  w->bspline_workspace_p = gsl_bspline_alloc(EEJ_BSPLINE_ORDER_T,
                                             nbreak);
  if (w->bspline_workspace_p == 0)
    {
      eej_basis_free(w);
      fprintf(stderr, "eej_basis_alloc: malloc failed\n");
      return 0;
    }

  w->B = gsl_vector_alloc(gsl_bspline_ncoeffs(w->bspline_workspace_p));
  if (w->B == 0)
    {
      eej_basis_free(w);
      fprintf(stderr, "eej_basis_alloc: malloc failed\n");
      return 0;
    }

  gsl_bspline_knots_uniform(EEJ_BASIS_LT_MIN,
                            EEJ_BASIS_LT_MAX,
                            w->bspline_workspace_p);

  return (w);
} /* eej_basis_alloc() */

void
eej_basis_free(eej_basis_workspace *w)
{
  if (!w)
    return;

  if (w->bspline_workspace_p)
    gsl_bspline_free(w->bspline_workspace_p);

  if (w->B)
    gsl_vector_free(w->B);

  free(w);
} /* eej_basis_free() */

/*
eej_basis_phi()
  Return basis function for longitude variable

u_{i}(phi) = { cos(i*phi/2),     i even
               sin((i+1)*phi/2), i odd  }

with phi \in [-pi, pi]
*/

int
eej_basis_phi(double phi, double y[], void *params)

{
  size_t i;

  for (i = 0; i < EEJ_MEAN_FIT_N_PHI_CHAMP; ++i)
    {
      if ((i % 2) == 0)
        y[i] = cos((double)(i/2) * phi);
      else
        y[i] = sin((double)((i+1)/2) * phi);
    }

  return GSL_SUCCESS;
} /* eej_basis_phi() */

/*
eej_basis_t()
  Return basis function for local time variable

u_{i}(t) = B_i(t)

Inputs: i - basis index
        t - local time (hours) (t \in [EEJ_BASIS_LT_MIN, EEJ_BASIS_LT_MAX])
*/

int
eej_basis_t(double t, double y[], void *params)

{
  eej_basis_workspace *w = (eej_basis_workspace *) params;

  if (w->satidx == EEJ_BASIS_SAC)
    {
      y[0] = 0.5;
    }
  else
    {
      gsl_vector_view v = gsl_vector_view_array(y, gsl_bspline_ncoeffs(w->bspline_workspace_p));

      gsl_bspline_eval(t, &v.vector, w->bspline_workspace_p);

#if 0
#ifdef EEJ_BSPLINE_FORCE_ZERO
      if (w->fit_type == EEJ_BASIS_MEAN)
        return (gsl_vector_get(w->B, i + EEJ_BSPLINE_ORDER_T - 1));
      else
#endif
        return (gsl_vector_get(w->B, i));
#endif
    }

  return GSL_SUCCESS;
} /* eej_basis_t() */

/*
eej_basis_s()
  Return basis function for season variable

u_{i}(s) = { cos(2i*pi*(s-s0)/365.25),  i even
             sin(2i*pi*(s-s0)/365.25),  i odd  }

with s \in [0, 365]
*/

int
eej_basis_s(double s, double y[], void *params)

{
  size_t i;
  double s0 = 79.0; /* spring equinox: 20 March */

  for (i = 0; i < EEJ_MEAN_FIT_N_S_CHAMP; ++i)
    {
      if ((i % 2) == 0)
        y[i] = cos(2.0*i*M_PI*(s - s0)/365.25);
      else
        y[i] = sin(2.0*i*M_PI*(s - s0)/365.25);
    }

  return GSL_SUCCESS;
} /* eej_basis_s() */

/*
eej_basis_euvac()
  Return basis function for EUVAC variable

u_{i}(f) = f^i
*/

int
eej_basis_euvac(double euvac, double y[], void *params)

{
  size_t i;

  for (i = 0; i < EEJ_MEAN_FIT_N_EUVAC_CHAMP; ++i)
    y[i] = pow(euvac, (double) i);

  return GSL_SUCCESS;
} /* eej_basis_euvac() */
