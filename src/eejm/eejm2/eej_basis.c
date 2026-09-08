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

size_t n_eejm_np = EEJ_MEAN_FIT_N_PHI_CHAMP;
size_t n_eejm_nt = EEJ_MEAN_FIT_N_T_CHAMP;
size_t n_eejm_ns = EEJ_MEAN_FIT_N_S_CHAMP;
size_t n_eejm_ne = EEJ_MEAN_FIT_N_EUVAC_CHAMP;
size_t n_eejm_nl = EEJ_MEAN_FIT_N_LUNAR_CHAMP;

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
  size_t i;
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

  for (i = 0; i < EEJ_FIT_DIM; ++i)
    w->N[i] = N[i];

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
  eej_basis_workspace *w = params;
  size_t i;

  for (i = 0; i < w->N[EEJ_IDX_PHI]; ++i)
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
      size_t n = gsl_bspline_ncoeffs(w->bspline_workspace_p);
      gsl_vector_view v = gsl_vector_view_array(y, n);
      int s;

      s = gsl_bspline_eval(t, &v.vector, w->bspline_workspace_p);
      if (s)
        return s;
    }

  return GSL_SUCCESS;
} /* eej_basis_t() */

/*
eej_basis_s()
  Return basis function for season variable

u_{i}(s) = { cos(i*pi*s/365.25),  i even
             sin(i*pi*s/365.25),  i odd  }

with s \in [0, 365]
*/

int
eej_basis_s(double s, double y[], void *params)

{
  eej_basis_workspace *w = params;
  size_t i;

  for (i = 0; i < w->N[EEJ_IDX_S]; ++i)
    {
      if ((i % 2) == 0)
        y[i] = cos(i * M_PI * s / 365.25);
      else
        y[i] = sin((i + 1) * M_PI * s / 365.25);
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
  eej_basis_workspace *w = params;
  size_t i;

  for (i = 0; i < w->N[EEJ_IDX_EUVAC]; ++i)
    y[i] = pow(euvac, (double) i);

  return GSL_SUCCESS;
} /* eej_basis_euvac() */

/*
eej_basis_lunar()
  Return basis function for lunar local time variable, tau is in
radians \in [0, 2pi]

u_{i}(phi_l) = { cos(i * tau),      i even
                 sin((i+1) * tau),  i odd  }
*/

int
eej_basis_lunar(double tau, double y[], void *params)

{
  eej_basis_workspace *w = params;
  size_t i;

  for (i = 0; i < w->N[EEJ_IDX_LUNAR]; ++i)
    {
      if ((i % 2) == 0)
        y[i] = cos(i * tau);
      else
        y[i] = sin((i + 1) * tau);
    }

  return GSL_SUCCESS;
} /* eej_basis_lunar() */
