/* eefm_basis.h
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

#ifndef INCLUDED_eefm_basis_h
#define INCLUDED_eefm_basis_h

#include <gsl/gsl_bspline.h>

/* dimension of fit function (phi, t, s, euvac, phi_lunar) */
#define EEF_FIT_DIM                       5

extern size_t n_eefm_np;
extern size_t n_eefm_nt;
extern size_t n_eefm_ns;
extern size_t n_eefm_ne;
extern size_t n_eefm_nl;

/* number of terms in each sum of mean fit functions */

#define EEF_MEAN_FIT_N_PHI_DEFAULT        9
#define EEF_MEAN_FIT_N_T_DEFAULT          6
#define EEF_MEAN_FIT_N_S_DEFAULT          5
#define EEF_MEAN_FIT_N_EUVAC_DEFAULT      2
#define EEF_MEAN_FIT_N_LUNAR_DEFAULT      3

/* number of terms in each sum of stddev fit functions */

#define EEF_STDDEV_FIT_N_PHI_DEFAULT      3
#define EEF_STDDEV_FIT_N_T_DEFAULT        6
#define EEF_STDDEV_FIT_N_S_DEFAULT        5
#define EEF_STDDEV_FIT_N_EUVAC_DEFAULT    2
#define EEF_STDDEV_FIT_N_LUNAR_DEFAULT    3

/* B-spline order for time variable (cubic) */
#define EEF_BSPLINE_ORDER_T               4

/* force bspline fit to zero at endpoints? */
#undef EEF_BSPLINE_FORCE_ZERO

/* order of variables in fit */
#define EEF_IDX_PHI                       0
#define EEF_IDX_T                         1
#define EEF_IDX_S                         2
#define EEF_IDX_EUVAC                     3
#define EEF_IDX_LUNAR                     4

/* satellite indices for eefm_basis_alloc() */
#define EEF_BASIS_CHAMP                   0
#define EEF_BASIS_OERSTED                 1
#define EEF_BASIS_SAC                     2

/* fit types for eefm_basis_alloc() */
#define EEF_BASIS_MEAN                    0
#define EEF_BASIS_STDDEV                  1

/* minimum/maximum local times for fit */
#define EEF_BASIS_LT_MIN                  (7.0)
#define EEF_BASIS_LT_MAX                  (17.0)

typedef struct
{
  gsl_bspline_workspace *bspline_workspace_p;
  gsl_bspline_workspace *bspline_workspace_sigma;
  size_t satidx;   /* EEF_BASIS_satname */
  size_t fit_type; /* EEF_BASIS_fittype */
  size_t N[EEF_FIT_DIM];
} eefm_basis_workspace;

/*
 * Prototypes
 */

eefm_basis_workspace *eefm_basis_alloc(const size_t N[], size_t satidx,
                                       size_t fit_type);
void eefm_basis_free(eefm_basis_workspace *w);

int eefm_basis_phi(double phi, double y[], void *params);
int eefm_basis_t(double t, double y[], void *params);
int eefm_basis_s(double t, double y[], void *params);
int eefm_basis_euvac(double t, double y[], void *params);
int eefm_basis_lunar(double phi, double y[], void *params);

#endif /* INCLUDED_eefm_basis_h */
