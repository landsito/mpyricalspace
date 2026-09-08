/*
 * drift_basis.h
 * Patrick Alken
 */

#ifndef INCLUDED_drift_basis_h
#define INCLUDED_drift_basis_h

#include <gsl/gsl_bspline.h>

/* dimension of fit function (t, s, euvac) */
#define N_DRIFT_FIT                3

/* number of terms in each sum of statistical fit functions */

extern size_t n_drift_nt;
extern size_t n_drift_ns;
extern size_t n_drift_ne;

#define N_DRIFT_SUM_DEFAULT_T      7
#define N_DRIFT_SUM_DEFAULT_S      11
#define N_DRIFT_SUM_DEFAULT_EUVAC  2

#define N_DRIFT_SUM_T              (n_drift_nt)
#define N_DRIFT_SUM_S              (n_drift_ns)
#define N_DRIFT_SUM_EUVAC          (n_drift_ne)

/* order of variables in fit */
#define IDX_DRIFT_T                0
#define IDX_DRIFT_S                1
#define IDX_DRIFT_EUVAC            2

#define DRIFT_LT_MIN               (8.0)
#define DRIFT_LT_MAX               (16.0)

#define DRIFT_EUVAC_MIN            (60.0)
#define DRIFT_EUVAC_MAX            (260.0)

typedef struct
{
  gsl_bspline_workspace *bspline_workspace_t;
} drift_basis_workspace;

/*
 * Prototypes
 */

drift_basis_workspace *drift_basis_alloc(const size_t N[]);
void drift_basis_free(drift_basis_workspace *w);

int get_drift_basis_t(double t, double y[], void *params);
int get_drift_basis_s(double s, double y[], void *params);
int get_drift_basis_euvac(double e, double y[], void *params);

#endif /* INCLUDED_drift_basis_h */
