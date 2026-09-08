/* eej_calc.h
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

#ifndef INCLUDED_eej_calc_h
#define INCLUDED_eej_calc_h

#include <gsl/gsl_vector.h>
#include "fit.h"
#include "eej_basis.h"

typedef struct
{
  gsl_vector *mean_coeffs;
  gsl_vector *stddev_coeffs;
  fit_workspace *fit_mean_p;
  fit_workspace *fit_stddev_p;
  eej_basis_workspace *eej_basis_mean_p;
  eej_basis_workspace *eej_basis_stddev_p;
  gsl_vector *x;

  int (*u[EEJ_FIT_DIM])(double x, double y[], void *p);
} eej_calc_workspace;

/* different possible models to use */
#define EEJ_CALC_CHAMP       0
#define EEJ_CALC_OERSTED     1
#define EEJ_CALC_SAC         2
#define EEJ_CALC_CHAMP_OER   3

/* default model */
#define EEJ_CALC_MODEL       EEJ_CALC_CHAMP

/*
 * Prototypes
 */

eej_calc_workspace *eej_calc_alloc(size_t model);
void eej_calc_free(eej_calc_workspace *w);
double eej_calc_mean(double phi, double t, double s, double e,
                     eej_calc_workspace *w);
double eej_calc_stddev(double phi, double t, double s, double e,
                       eej_calc_workspace *w);

#endif /* INCLUDED_eej_calc_h */
