/* eefm_basis.c
 * 
 * Copyright (C) 2008 Patrick Alken
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
 * in the statistical model of the EEF.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_bspline.h>

#include "eefm_basis.h"

size_t n_eefm_np = EEF_MEAN_FIT_N_PHI_DEFAULT;
size_t n_eefm_nt = EEF_MEAN_FIT_N_T_DEFAULT;
size_t n_eefm_ns = EEF_MEAN_FIT_N_S_DEFAULT;
size_t n_eefm_ne = EEF_MEAN_FIT_N_EUVAC_DEFAULT;
size_t n_eefm_nl = EEF_MEAN_FIT_N_LUNAR_DEFAULT;

/*
eefm_basis_alloc()
  Allocate a eefm_basis workspace

Inputs: N        - array containing number of terms in each sum
        satidx   - index of satellite (EEF_BASIS_satname)
        fit_type - type of fit (mean/stddev) (EEF_BASIS_type)

Return: pointer to workspace
*/

eefm_basis_workspace *
eefm_basis_alloc(const size_t N[], size_t satidx, size_t fit_type)
{
  size_t i;
  eefm_basis_workspace *w;
  size_t nbreak;

  w = (eefm_basis_workspace *) calloc(1, sizeof(eefm_basis_workspace));
  if (w == 0)
    {
      fprintf(stderr, "eefm_basis_alloc: malloc failed\n");
      return 0;
    }

  w->satidx = satidx;
  w->fit_type = fit_type;

  for (i = 0; i < EEF_FIT_DIM; ++i)
    w->N[i] = N[i];

  if (N[EEF_IDX_T] < 4)
    {
      fprintf(stderr, "warning: N_t < 4, bspline workspace not initialized\n");
      return (w); /* don't initialize bspline stuff */
    }

#ifdef EEF_BSPLINE_FORCE_ZERO
  if (fit_type == EEF_BASIS_MEAN)
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
      nbreak = N[EEF_IDX_T] + EEF_BSPLINE_ORDER_T;
    }
  else
#endif /* EEF_BSPLINE_FORCE_ZERO */
    nbreak = N[EEF_IDX_T] - EEF_BSPLINE_ORDER_T + 2;

  w->bspline_workspace_p = gsl_bspline_alloc(EEF_BSPLINE_ORDER_T,
                                             nbreak);
  if (w->bspline_workspace_p == 0)
    {
      eefm_basis_free(w);
      fprintf(stderr, "eefm_basis_alloc: malloc failed\n");
      return 0;
    }

  gsl_bspline_knots_uniform(EEF_BASIS_LT_MIN,
                            EEF_BASIS_LT_MAX,
                            w->bspline_workspace_p);

#if 0
  if (fit_type == EEF_BASIS_STDDEV)
    {
      if (N[EEF_IDX_PHI] < 4)
        {
          fprintf(stderr, "warning: stddev: N_phi < 4, bspline workspace not initialized\n");
          exit(1);
        }

      nbreak = N[EEF_IDX_PHI] - EEF_BSPLINE_ORDER_T + 2;
      w->bspline_workspace_sigma = gsl_bspline_alloc(EEF_BSPLINE_ORDER_T,
                                                     nbreak);
      gsl_bspline_knots_uniform(-M_PI,
                                M_PI,
                                w->bspline_workspace_sigma);
    }
#endif

  return (w);
} /* eefm_basis_alloc() */

void
eefm_basis_free(eefm_basis_workspace *w)
{
  if (!w)
    return;

  if (w->bspline_workspace_p)
    gsl_bspline_free(w->bspline_workspace_p);

  if (w->bspline_workspace_sigma)
    gsl_bspline_free(w->bspline_workspace_sigma);

  free(w);
} /* eefm_basis_free() */

/*
eefm_basis_phi()
  Return basis function for longitude variable

u_{i}(phi) = { cos(i*phi/2),     i even
               sin((i+1)*phi/2), i odd  }

with phi \in [-pi, pi]
*/

int
eefm_basis_phi(double phi, double y[], void *params)

{
  eefm_basis_workspace *w = params;
  size_t i;

#if 0
  if (w->fit_type == EEF_BASIS_STDDEV)
    {
      size_t n = gsl_bspline_ncoeffs(w->bspline_workspace_sigma);
      gsl_vector_view v = gsl_vector_view_array(y, n);
      int s;

      s = gsl_bspline_eval(phi, &v.vector, w->bspline_workspace_sigma);
      if (s)
        return s;
    }
  else
#endif
    {
      for (i = 0; i < w->N[EEF_IDX_PHI]; ++i)
        {
          if ((i % 2) == 0)
            y[i] = cos((double)(i/2) * phi);
          else
            y[i] = sin((double)((i+1)/2) * phi);
        }
    }

  return GSL_SUCCESS;
} /* eefm_basis_phi() */

/*
eefm_basis_t()
  Return basis function for local time variable

u_{i}(t) = B_i(t)

Inputs: i - basis index
        t - local time (hours) (t \in [EEF_BASIS_LT_MIN, EEF_BASIS_LT_MAX])
*/

int
eefm_basis_t(double t, double y[], void *params)

{
  eefm_basis_workspace *w = (eefm_basis_workspace *) params;

  if (w->satidx == EEF_BASIS_SAC)
    {
      y[0] = 0.5;
    }
  else
    {
      size_t n = gsl_bspline_ncoeffs(w->bspline_workspace_p);
      gsl_vector_view v = gsl_vector_view_array(y, n);
      int s;

      s = gsl_bspline_eval(t, &v.vector, w->bspline_workspace_p);
      if (s)
        return s;
    }

  return GSL_SUCCESS;
} /* eefm_basis_t() */

/*
eefm_basis_s()
  Return basis function for season variable

u_{i}(s) = { cos(i*pi*s/365.25),  i even
             sin(i*pi*s/365.25),  i odd  }

with s \in [0, 365]
*/

int
eefm_basis_s(double s, double y[], void *params)

{
  eefm_basis_workspace *w = params;
  size_t i;

  for (i = 0; i < w->N[EEF_IDX_S]; ++i)
    {
      if ((i % 2) == 0)
        y[i] = cos(i * M_PI * s / 365.25);
      else
        y[i] = sin((i + 1) * M_PI * s / 365.25);
    }

  return GSL_SUCCESS;
} /* eefm_basis_s() */

/*
eefm_basis_euvac()
  Return basis function for EUVAC variable

u_{i}(f) = f^i
*/

int
eefm_basis_euvac(double euvac, double y[], void *params)

{
  eefm_basis_workspace *w = params;
  size_t i;

  for (i = 0; i < w->N[EEF_IDX_EUVAC]; ++i)
    y[i] = pow(euvac, (double) i);

  return GSL_SUCCESS;
} /* eefm_basis_euvac() */

/*
eefm_basis_lunar()
  Return basis function for lunar local time variable. tau is in
radians \in [0, 2pi]

u_{i}(tau) = { cos(i * tau),     i even
               sin((i+1) * tau), i odd  }
*/

int
eefm_basis_lunar(double tau, double y[], void *params)

{
  eefm_basis_workspace *w = params;
  size_t i;

  for (i = 0; i < w->N[EEF_IDX_LUNAR]; ++i)
    {
      if ((i % 2) == 0)
        y[i] = cos(i * tau);
      else
        y[i] = sin((i + 1) * tau);
    }

  return GSL_SUCCESS;
} /* eefm_basis_lunar() */
