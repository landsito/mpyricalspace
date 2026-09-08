/* fit.h
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

#ifndef INCLUDED_fit_h
#define INCLUDED_fit_h

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

typedef struct
{
  size_t n;        /* dimension of fit function */
  size_t *N;       /* number of terms in fit sums N[i] = N_i */
  size_t n_coeffs; /* number of fit coefficients */
  size_t *r;       /* indices of basis functions of fit (n of these) */
  gsl_vector *work; /* scratch array of size n_coeffs */
  gsl_vector *work2; /* scratch array */

  /*
   * Views into the 'work' array which will be used to store
   * the results of calling the basis functions, so that
   * (v[i])_j = u^{(i)}_j(x_i)
   */
  gsl_vector_view *v;

  /* pointer to basis functions and parameters */
  int (**u)(double x, double y[], void *p);
  void *params;
} fit_workspace;

/*
 * Prototypes
 */

fit_workspace *fit_ndim_alloc(size_t n, size_t N[],
                              int (**u)(double x, double y[], void *p),
                              void *params);
void fit_ndim_free(fit_workspace *w);
int fit_ndim_design(const gsl_matrix *data, gsl_matrix *X,
                    fit_workspace *w);
int fit_ndim_calc(const gsl_vector *x, const gsl_vector *c,
                  const gsl_matrix *cov, double *y, double *y_err,
                  fit_workspace *w);
double fit_ndim_calc2(const gsl_vector *x, const gsl_vector *c,
                      fit_workspace *w);
int fit_ndim_residuals(gsl_matrix *X, gsl_vector *b, gsl_vector *x,
                       gsl_vector *r);

#endif /* INCLUDED_fit_h */
