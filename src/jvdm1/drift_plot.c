/*
 * drift_plot.c
 * Patrick Alken
 *
 * Output data plotting the model for each variable while holding
 * the other parameters constant.
 *
 * ./drift_plot -h for options
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <getopt.h>
#include <assert.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>

#include "drift_basis.h"
#include "drift_calc.h"
#include "drift_files.h"
#include "drift_plot.h"

const char *param_str[] = {
  "time (hours)",
  "season (day of year)",
  "euvac (W/m^2)"
};

/*
drift_plot_alloc()
  Initialize a drift workspace
*/

drift_plot_workspace *
drift_plot_alloc(void)

{
  drift_plot_workspace *w;

  w = (drift_plot_workspace *) calloc(1, sizeof(drift_plot_workspace));
  if (!w)
    {
      perror("calloc");
      return (0);
    }

  w->t = (double *) malloc(sizeof(double) * MAX_DATA);
  w->s = (double *) malloc(sizeof(double) * MAX_DATA);
  w->euvac = (double *) malloc(sizeof(double) * MAX_DATA);
  w->v = (double *) malloc(sizeof(double) * MAX_DATA);

  w->drift_calc_workspace_p = drift_calc_alloc();

  w->plot_lt = 10.5;
  w->plot_euvac = 100.0;
  w->plot_season = 79.0;

  w->plot_idx = IDX_DRIFT_T;
  w->plot_data = 0;

  return (w);
} /* drift_plot_alloc() */

/*
drift_plot_free()
  Terminate a drift workspace
*/

void
drift_plot_free(drift_plot_workspace *w)

{
  if (!w)
    return;

  if (w->t)
    free(w->t);

  if (w->s)
    free(w->s);

  if (w->euvac)
    free(w->euvac);

  if (w->v)
    free(w->v);

  if (w->drift_calc_workspace_p)
    drift_calc_free(w->drift_calc_workspace_p);

  free(w);
} /* drift_plot_free() */

/*
drift_plot_proc()
  Plot drift data

Inputs: w - drift workspace
*/

void
drift_plot_proc(drift_plot_workspace *w)

{
  size_t i, j, k;
  FILE *fp;
  double *param_data[N_DRIFT_FIT];
  double param_plot[N_DRIFT_FIT];
  double param_min[N_DRIFT_FIT];
  double param_max[N_DRIFT_FIT];
  double param_step[N_DRIFT_FIT];
  double param_dplot[] = { 1.0,
                           5.0,
                           20.0 };
  const char *param_files[] = { FILE_MODEL_T_OUT,
                                FILE_MODEL_S_OUT,
                                FILE_MODEL_EUVAC_OUT };
  double *x;
  double y;
  double v;        /* model mean value */
  double sigma;    /* model sigma value */
  double plotval;
  double dplotval;
  size_t cnt_data;

  /* parameter data locations */
  param_data[IDX_DRIFT_T] = w->t;
  param_data[IDX_DRIFT_S] = w->s;
  param_data[IDX_DRIFT_EUVAC] = w->euvac;

  /* parameter constant plot values */
  param_plot[IDX_DRIFT_T] = w->plot_lt;
  param_plot[IDX_DRIFT_S] = w->plot_season;
  param_plot[IDX_DRIFT_EUVAC] = w->plot_euvac;

  /* parameter minimum values */
  param_min[IDX_DRIFT_T] = DRIFT_LT_MIN;
  param_min[IDX_DRIFT_S] = 0.0;
  param_min[IDX_DRIFT_EUVAC] = DRIFT_EUVAC_MIN;

  /* parameter maximum values */
  param_max[IDX_DRIFT_T] = DRIFT_LT_MAX;
  param_max[IDX_DRIFT_S] = 365.0;
  param_max[IDX_DRIFT_EUVAC] = DRIFT_EUVAC_MAX;

  /* parameter step values */
  param_step[IDX_DRIFT_T] = 0.01;
  param_step[IDX_DRIFT_S] = 1.0;
  param_step[IDX_DRIFT_EUVAC] = 1.0;

  /* plot_idx contains the variable we want to plot */
  i = w->plot_idx;

  fp = fopen(param_files[i], "w");
  if (!fp)
    {
      perror("fopen");
      return;
    }

//  fprintf(stderr, "Writing %s...", param_files[i]);

  /* now print out model */

  j = 1;
  fprintf(fp, "# Index 0:\n");
  fprintf(fp, "# Field %zu: %s\n", j++, param_str[i]);
  fprintf(fp, "# Field %zu: v (m/s)\n", j++);
  fprintf(fp, "# Field %zu: sigma (m/s)\n", j++);

  fprintf(fp, "#\n# Parameter values:\n");
  for (j = 0; j < N_DRIFT_FIT; ++j)
    {
      if (j == i)
        continue;

      fprintf(fp, "# %s = %f\n", param_str[j], param_plot[j]);
    }

  if (w->plot_data)
    {
      j = 1;
      fprintf(fp, "#\n# Index 1:\n");
      fprintf(fp, "# Field %zu: time (hours)\n", j++);
      fprintf(fp, "# Field %zu: season (days)\n", j++);
      fprintf(fp, "# Field %zu: EUVAC (W/m^2)\n", j++);
      fprintf(fp, "# Field %zu: v (m/s)\n", j++);
    }

  plotval = param_plot[i]; /* save constant plot value */

  for (y = param_min[i]; y < param_max[i]; y += param_step[i])
    {
      param_plot[i] = y;
      v = drift_calc_mean(param_plot[IDX_DRIFT_T],
                          param_plot[IDX_DRIFT_S],
                          param_plot[IDX_DRIFT_EUVAC],
                          w->drift_calc_workspace_p);

      sigma = drift_calc_stddev(param_plot[IDX_DRIFT_T],
                                param_plot[IDX_DRIFT_S],
                                param_plot[IDX_DRIFT_EUVAC],
                                w->drift_calc_workspace_p);

      fprintf(fp, "%f %f %f\n", y, v, sigma);
    }

  param_plot[i] = plotval; /* restore constant plot value */

  if (w->plot_data)
    {
      fprintf(fp, "\n\n");

      /*
       * The 'j' loop goes through all data points and each point
       * is tested via the 'k' loop to make sure it fits within
       * the constant values of the other parameters
       */

      /* print out selected data */
      cnt_data = 0;
      for (j = 0; j < w->n; ++j)
        {
          int foundpoint = 1;

          for (k = 0; k < N_DRIFT_FIT; ++k)
            {
              /*
               * we are plotting vs param_i so ignore that when
               * we get to it
               */
              if (k == i)
                continue;

              /*
               * Test current parameter and data point to see
               * if it is within the desired plot limits:
               *
               * paramval_j \in [plotval_k - eps, plotval_k + eps]
               */
              x = param_data[k];
              plotval = param_plot[k];
              dplotval = param_dplot[k];

              if ((x[j] > (plotval + dplotval)) || 
                  (x[j] < (plotval - dplotval)))
                {
                  foundpoint = 0;
                  break;
                }
            }

          if (!foundpoint)
            continue;

          fprintf(fp,
                  "%f %f %f %f\n",
                  w->t[j],
                  w->s[j],
                  w->euvac[j],
                  w->v[j]);

          ++cnt_data;
        }

      fprintf(stderr, "found %zu data points...", cnt_data);
    } /* if (w->plot_data) */

  fclose(fp);

//  fprintf(stderr, "done\n");
} /* drift_plot_proc() */

int
main(int argc, char *argv[])

{
  drift_plot_workspace *drift_plot_workspace_p;
  int c;

  drift_plot_workspace_p = drift_plot_alloc();
  if (!drift_plot_workspace_p)
  {
    fprintf(stderr, "main: drift_plot_alloc failed\n");
    exit(1);
  }

  while ((c = getopt(argc, argv, "t:s:e:dp:")) != (-1))
  {
    switch (c)
    {
      case 't':
        drift_plot_workspace_p->plot_lt = strtod(optarg, NULL);
        break;

      case 's':
        drift_plot_workspace_p->plot_season = strtod(optarg, NULL);
        break;

      case 'e':
        drift_plot_workspace_p->plot_euvac = strtod(optarg, NULL);
        break;

      case 'd':
        drift_plot_workspace_p->plot_data = 1;
        break;

      case 'p':
        {
          switch (*optarg)
            {
              case 't':
                drift_plot_workspace_p->plot_idx = IDX_DRIFT_T;
                break;

              case 's':
                drift_plot_workspace_p->plot_idx = IDX_DRIFT_S;
                break;

              case 'e':
                drift_plot_workspace_p->plot_idx = IDX_DRIFT_EUVAC;
                break;

              default:
                break;
            }
          break;
        }

      case '?':
      default:
        printf("usage: %s [options]\n", argv[0]);
        printf("\n\nOptions:\n\n");
        printf("  -t <local time> : Specifies value of time parameter in hours [default: %f]\n",
               drift_plot_workspace_p->plot_lt);
        printf("  -s <season>     : Specifies value of season parameter in days [default: %f]\n",
               drift_plot_workspace_p->plot_season);
        printf("  -e <euvac>      : Specifies value of EUVAC index [default: %f]\n",
               drift_plot_workspace_p->plot_euvac);
        printf("  -p <t|s|e>    : Specifies which variable to plot against (t = local time, s = season, e = EUVAC) [default: local time]\n");
        exit(1);
        break;
    }
  }

  drift_plot_proc(drift_plot_workspace_p);

  drift_plot_free(drift_plot_workspace_p);

  return (0);
} /* main() */
