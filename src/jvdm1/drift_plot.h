/*
 * drift_plot.h
 * Patrick Alken
 */

#ifndef INCLUDED_drift_plot_h
#define INCLUDED_drift_plot_h

#include "drift_basis.h"
#include "drift_calc.h"

typedef struct
{
  double *t;            /* local time in hours */
  double *s;            /* season in days */
  double *euvac;        /* EUVAC */
  double *v;            /* drift in m/s */

  size_t n;             /* number of data points read */

  double plot_lt;       /* local time value for plot */
  double plot_euvac;    /* euvac value for plot */
  double plot_season;   /* season value for plot */

  size_t plot_idx;      /* index of parameter to plot */
  int plot_data;        /* plot data points? */

  drift_calc_workspace *drift_calc_workspace_p;
} drift_plot_workspace;

#define MAX_FILENAME                100

#define MAX_DATA                    80000

#define MAX_BUFFER                  5000

/*
 * Prototypes
 */

drift_plot_workspace *drift_plot_alloc(void);
void drift_plot_free(drift_plot_workspace *w);
void drift_plot_proc(drift_plot_workspace *w);

#endif /* INCLUDED_drift_plot_h */
