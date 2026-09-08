/* eefm_calc.h
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

#ifndef INCLUDED_eefm_calc_h
#define INCLUDED_eefm_calc_h

#include <gsl/gsl_vector.h>
#include <ndlinear/gsl_multifit_ndlinear.h>
#include "eefm_basis.h"

typedef struct
{
  gsl_vector *mean_coeffs;
  gsl_vector *stddev_coeffs;
  gsl_multifit_ndlinear_workspace *fit_mean_p;
  gsl_multifit_ndlinear_workspace *fit_stddev_p;
  eefm_basis_workspace *eefm_basis_mean_p;
  eefm_basis_workspace *eefm_basis_stddev_p;
  gsl_vector *x;

  int (*u[EEF_FIT_DIM])(double x, double y[], void *p);
} eefm_calc_workspace;

/* different possible models to use */
#define EEF_CALC_CHAMP       0
#define EEF_CALC_OERSTED     1
#define EEF_CALC_SAC         2
#define EEF_CALC_CHAMP_OER   3

/* default model */
#define EEF_CALC_MODEL       EEF_CALC_CHAMP

/*
 * Prototypes
 */

eefm_calc_workspace *eefm_calc_alloc(size_t model);
void eefm_calc_free(eefm_calc_workspace *w);
double eefm_calc_mean(double phi, double t, double s, double e, double tau,
                      eefm_calc_workspace *w);
double eefm_calc_stddev(double phi, double t, double s, double e,
                        double tau, eefm_calc_workspace *w);

#endif /* INCLUDED_eefm_calc_h */
