/*
 * rteef.h
 */

#ifndef INCLUDED_rteef_h
#define INCLUDED_rteef_h

#include <time.h>

#include "f107.h"
#include "eefm_calc.h"

typedef struct
{
  double a[50];
  double b[50];
  size_t na;
  size_t nb;
  double *data_ief;     /* IEF (ACE) data */
  size_t ndata_ief;     /* total number of IEF data */
  size_t ndata_ief_max; /* maximum number of IEF data */
  double *eef_rt;       /* RT model output */
  double *eef_mean;     /* Mean model output */
  double *dbuffer;

  eefm_calc_workspace *eefm_calc_p;
  f107_workspace *f107_p;
} rteef_workspace;

/*
 * Prototypes
 */

rteef_workspace *rteef_alloc(void);
void rteef_free(rteef_workspace *w);

/*
 * rteef_compute(): if ief_override is non-NULL it is used as the IEF Ey input
 * (nief samples, 5-min cadence, starting one hour before t_start) instead of
 * reading ace<YYYY>.dat, and the climatological term is skipped (eef_mean stays
 * 0) -- a pure transfer-function evaluation. Pass NULL, 0 for the normal path.
 */
int rteef_compute(time_t t_start, time_t t_end, double longitude,
                  rteef_workspace *w, const double *ief_override, size_t nief);

#endif /* INCLUDED_rteef_h */
