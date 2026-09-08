/* fit.c
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

#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>

#include "fit.h"

static inline int fit_ndim_design_row(const gsl_vector *d,
                                      gsl_vector *x,
                                      fit_workspace *w);

/*
fit_ndim_alloc()
  Allocate a fit workspace

Inputs: n      - dimension of fit function
        N      - number of terms in each sum; N[i] = N_i, 0 <= i < n
        u      - basis functions to call
                 u[j] = u^{(j)}, 0 <= j < n
        params - parameters to pass to basis functions

Return: pointer to new workspace

Notes: the supplied basis functions 'u[j]' must accept three
       arguments:

int uj(double x, double y[], void *params)

and fill the y[] vector so that y[i] = u_{i}^{(j)}(x) (the ith
basis function for the jth parameter evaluated at x)
*/

fit_workspace *
fit_ndim_alloc(size_t n, size_t N[],
               int (**u)(double x, double y[], void *p), void *params)
{
  fit_workspace *w;
  size_t n_coeffs; /* total number of fit coefficients */
  size_t i, idx;
  size_t sum_N;

  if (n == 0)
    {
      GSL_ERROR_NULL("n must be at least 1", GSL_EINVAL);
    }

  w = (fit_workspace *) calloc(1, sizeof(fit_workspace));
  if (!w)
    {
      GSL_ERROR_NULL("failed to allocate space for workspace", GSL_ENOMEM);
    }

  w->r = (size_t *) calloc(n, sizeof(size_t));
  if (!w->r)
    {
      fit_ndim_free(w);
      GSL_ERROR_NULL("failed to allocate space for index vector", GSL_ENOMEM);
    }

  w->N = (size_t *) calloc(n, sizeof(size_t));
  if (!w->N)
    {
      fit_ndim_free(w);
      GSL_ERROR_NULL("failed to allocate space for N vector", GSL_ENOMEM);
    }

  n_coeffs = 1;
  sum_N = 0;
  for (i = 0; i < n; ++i)
    {
      if (N[i] == 0)
        {
          fit_ndim_free(w);
          GSL_ERROR_NULL("one of the sums is empty", GSL_EINVAL);
        }

      /*
       * The total number of coefficients is: N_1 * N_2 * ... * N_n
       */
      n_coeffs *= N[i];
      w->N[i] = N[i];
      sum_N += N[i];
    }

  w->n = n;
  w->n_coeffs = n_coeffs;

  w->work = gsl_vector_alloc(n_coeffs);
  w->work2 = gsl_vector_alloc(sum_N);
  if (!w->work || !w->work2)
    {
      fit_ndim_free(w);
      GSL_ERROR_NULL("failed to allocate space for basis vector", GSL_ENOMEM);
    }

  w->v = (gsl_vector_view *) calloc(n, sizeof(gsl_vector_view));
  if (!w->v)
    {
      fit_ndim_free(w);
      GSL_ERROR_NULL("failed to allocate space for basis vector", GSL_ENOMEM);
    }

  idx = 0;
  for (i = 0; i < n; ++i)
    {
      w->v[i] = gsl_vector_subvector(w->work2, idx, N[i]);
      idx += N[i];
    }

  w->u = u;
  w->params = params;

  return (w);
} /* fit_ndim_alloc() */

/*
fit_ndim_free()
  Free workspace w
*/

void
fit_ndim_free(fit_workspace *w)
{
  if (!w)
    return;

  if (w->r)
    free(w->r);

  if (w->N)
    free(w->N);

  if (w->work)
    gsl_vector_free(w->work);

  if (w->work2)
    gsl_vector_free(w->work2);

  if (w->v)
    free(w->v);

  free(w);
} /* fit_ndim_free() */

/*
fit_ndim_design()
  This function constructs the coefficient design matrix 'X'

Inputs: data  - data vectors for matrix X
                data is a ndata-by-n matrix with the ith row
                given by x_i, so that
                data_{ij} = (x_i)_j, the jth element of the
                ith data vector
        X     - (output) design matrix (must be ndata by w->n_coeffs)
        w     - workspace

Return: success or error
*/

int
fit_ndim_design(const gsl_matrix *data, gsl_matrix *X, fit_workspace *w)
{
  const size_t ndata = data->size1;
  size_t i; /* looping */
  int s;

  if ((X->size1 != ndata) || (X->size2 != w->n_coeffs))
    {
      GSL_ERROR("X matrix has wrong dimensions", GSL_EBADLEN);
    }
  else
    {
      for (i = 0; i < ndata; ++i)
        {
          gsl_vector_const_view d = gsl_matrix_const_row(data, i);
          gsl_vector_view xv = gsl_matrix_row(X, i);

          s = fit_ndim_design_row(&d.vector, &xv.vector, w);
          if (s != GSL_SUCCESS)
            return s;
        }

      return GSL_SUCCESS;
    }
} /* fit_ndim_design() */

/*
fit_ndim_calc()
  Compute the fit function at a given data point

Inputs: x     - data point (w->n elements)
        c     - coefficient vector
        cov   - covariance matrix
        y     - where to store fit function result
        y_err - standard deviation of fit
        w     - workspace

Return: success or error
*/

int
fit_ndim_calc(const gsl_vector *x, const gsl_vector *c,
              const gsl_matrix *cov, double *y, double *y_err,
              fit_workspace *w)

{
  if (c->size != w->n_coeffs)
    {
      GSL_ERROR("c vector has wrong size", GSL_EBADLEN);
    }
  else
    {
      int s;

      s = fit_ndim_design_row(x, w->work, w);
      if (s != GSL_SUCCESS)
        return s;

      /*
       * Now w->work contains the appropriate basis functions
       * evaluated at the given point - compute the function value
       */
      s = gsl_multifit_linear_est(w->work, c, cov, y, y_err);

      return s;
    }
} /* fit_ndim_calc() */

/*
fit_ndim_calc2()
  Compute the fit function at a given data point

Inputs: x - data point (w->n elements)
        c - coefficient vector
        w - fit workspace

Return: value of fit function f(x)
*/

double
fit_ndim_calc2(const gsl_vector *x, const gsl_vector *c,
               fit_workspace *w)
{
  double y;
  int s;

  s = fit_ndim_design_row(x, w->work, w);
  if (s != GSL_SUCCESS)
    return s;

  gsl_blas_ddot(w->work, c, &y);

  return y;
}

/*
fit_ndim_residuals()
  Compute vector of residuals from fit

Inputs: X - design matrix
        b - rhs vector
        x - fit coefficients
        r - (output) where to store residuals
*/

int
fit_ndim_residuals(gsl_matrix *X, gsl_vector *b, gsl_vector *x,
                   gsl_vector *r)
{
  const size_t n = b->size;
  size_t i;

  for (i = 0; i < n; ++i)
    {
      double bi = gsl_vector_get(b, i);
      gsl_vector_const_view row = gsl_matrix_const_row(X, i);
      double b_est, ri;

      gsl_blas_ddot(&row.vector, x, &b_est);
      ri = bi - b_est;

      gsl_vector_set(r, i, ri);
    }

  return (GSL_SUCCESS);
} /* fit_ndim_residuals() */

/******************************************
 *         INTERNAL ROUTINES              *
 ******************************************/

/*
fit_ndim_design_row()

  Compute a row of the design matrix X:

X(:,j) =  A_{:,j} = u_{r_0}^{(0)}(d_0) * u_{r_1)^{(1)}(d_1) * ... *
                    u_{r_{n-1}}^{(n-1)}(d_{n-1})

where 'd' is the corresponding data vector for that row

Inputs: d - data vector of length n
        x - (output) where to store row of design matrix X
        w - workspace

Return: success or error
*/

static inline int
fit_ndim_design_row(const gsl_vector *d, gsl_vector *x,
                    fit_workspace *w)
{
  size_t j;
  int k;
  size_t denom;
  double melement;

  /* compute basis functions for this data point */
  for (j = 0; j < w->n; ++j)
    {
      int s = w->u[j](gsl_vector_get(d, j),
                      w->v[j].vector.data,
                      w->params);

      if (s != GSL_SUCCESS)
        return s;
    }

  for (j = 0; j < w->n_coeffs; ++j)
    {
      /*
       * The (:,j) element of the matrix X will be:
       *
       * X_{:,j} = u_{r_0}^{(0)}(d0) *
       *           u_{r_1)^{(1)}(d1) *
       *           ... *
       *           u_{r_{n-1}}^{(n-1)}(d_{n-1})
       *
       * with the basis function indices r_k given by
       *
       * r_k = floor(j / Prod_{i=(k+1)..(n-1)} [ N_i ]) (mod N_k)
       *
       * In the case where N_i = N for all i,
       *
       * r_k = floor(j / N^{n - k - 1}) (mod N)
       *
       * n: dimension of fit function
       * N_i: number of terms in sum i of fit function
       */

      /* calculate the r_k and the matrix element X_{:,j} */

      denom = 1;
      melement = 1.0;
      for (k = (int)(w->n - 1); k >= 0; --k)
        {
          w->r[k] = (j / denom) % w->N[k];
          denom *= w->N[k];

          melement *= gsl_vector_get(&(w->v[k]).vector, w->r[k]);
        }

      /* set the matrix element */
      gsl_vector_set(x, j, melement);
    }

  return GSL_SUCCESS;
} /* fit_ndim_design_row() */
