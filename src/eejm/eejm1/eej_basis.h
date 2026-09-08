/* eej_basis.h
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

#ifndef INCLUDED_eej_basis_h
#define INCLUDED_eej_basis_h

#include <gsl/gsl_bspline.h>

/* dimension of fit function (phi, t, s, euvac) */
#define EEJ_FIT_DIM                  4

/* number of terms in each sum of mean fit functions */

#define EEJ_MEAN_FIT_N_PHI_CHAMP     9
#define EEJ_MEAN_FIT_N_T_CHAMP       8
#define EEJ_MEAN_FIT_N_S_CHAMP       3
#define EEJ_MEAN_FIT_N_EUVAC_CHAMP   2

#define EEJ_MEAN_FIT_N_PHI_OER       EEJ_MEAN_FIT_N_PHI_CHAMP
#define EEJ_MEAN_FIT_N_T_OER         EEJ_MEAN_FIT_N_T_CHAMP
#define EEJ_MEAN_FIT_N_S_OER         EEJ_MEAN_FIT_N_S_CHAMP
#define EEJ_MEAN_FIT_N_EUVAC_OER     EEJ_MEAN_FIT_N_EUVAC_CHAMP

#define EEJ_MEAN_FIT_N_PHI_SAC       EEJ_MEAN_FIT_N_PHI_CHAMP
#define EEJ_MEAN_FIT_N_T_SAC         1
#define EEJ_MEAN_FIT_N_S_SAC         EEJ_MEAN_FIT_N_S_CHAMP
#define EEJ_MEAN_FIT_N_EUVAC_SAC     EEJ_MEAN_FIT_N_EUVAC_CHAMP

/* number of terms in each sum of stddev fit functions */

#define EEJ_STDDEV_FIT_N_PHI_CHAMP   9
#define EEJ_STDDEV_FIT_N_T_CHAMP     7
#define EEJ_STDDEV_FIT_N_S_CHAMP     3
#define EEJ_STDDEV_FIT_N_EUVAC_CHAMP 2

#define EEJ_STDDEV_FIT_N_PHI_OER     EEJ_STDDEV_FIT_N_PHI_CHAMP
#define EEJ_STDDEV_FIT_N_T_OER       EEJ_STDDEV_FIT_N_T_CHAMP
#define EEJ_STDDEV_FIT_N_S_OER       EEJ_STDDEV_FIT_N_S_CHAMP
#define EEJ_STDDEV_FIT_N_EUVAC_OER   EEJ_STDDEV_FIT_N_EUVAC_CHAMP

#define EEJ_STDDEV_FIT_N_PHI_SAC     EEJ_STDDEV_FIT_N_PHI_CHAMP
#define EEJ_STDDEV_FIT_N_T_SAC       1
#define EEJ_STDDEV_FIT_N_S_SAC       EEJ_STDDEV_FIT_N_S_CHAMP
#define EEJ_STDDEV_FIT_N_EUVAC_SAC   EEJ_STDDEV_FIT_N_EUVAC_CHAMP

/* B-spline order for time variable (cubic) */
#define EEJ_BSPLINE_ORDER_T          4

/* force bspline fit to zero at endpoints? */
#undef EEJ_BSPLINE_FORCE_ZERO

/* order of variables in fit */
#define EEJ_IDX_PHI                  0
#define EEJ_IDX_T                    1
#define EEJ_IDX_S                    2
#define EEJ_IDX_EUVAC                3

/* satellite indices for eej_basis_alloc() */
#define EEJ_BASIS_CHAMP              0
#define EEJ_BASIS_OERSTED            1
#define EEJ_BASIS_SAC                2

/* fit types for eej_basis_alloc() */
#define EEJ_BASIS_MEAN               0
#define EEJ_BASIS_STDDEV             1

/* minimum/maximum local times for fit */
#define EEJ_BASIS_LT_MIN             (5.0)
#define EEJ_BASIS_LT_MAX             (19.0)

typedef struct
{
  gsl_vector *B;   /* spline values vector */
  gsl_bspline_workspace *bspline_workspace_p;
  size_t satidx;   /* EEJ_BASIS_satname */
  size_t fit_type; /* EEJ_BASIS_fittype */
} eej_basis_workspace;

/*
 * Prototypes
 */

eej_basis_workspace *eej_basis_alloc(const size_t N[], size_t satidx,
                                     size_t fit_type);
void eej_basis_free(eej_basis_workspace *w);

int eej_basis_phi(double phi, double y[], void *params);
int eej_basis_t(double t, double y[], void *params);
int eej_basis_s(double t, double y[], void *params);
int eej_basis_euvac(double t, double y[], void *params);

#endif /* INCLUDED_eej_basis_h */
