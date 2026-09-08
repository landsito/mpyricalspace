/* eefm_plot.c
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

/*
 * Output data plotting the model for each variable while holding
 * the other parameters constant.
 *
 * ./eefm_plot -h for options
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

#include "eefm_basis.h"
#include "eefm_calc.h"
#include "eefm_files.h"
#include "eefm_plot.h"

const char *param_str[] = {
  "longitude",
  "local time",
  "season",
  "euvac",
  "lunar local time"
};

const char *param_units[] = {
  "degrees",
  "hours",
  "day of year",
  "W/m^2",
  "lunar hours"
};

const double param_stepsize[] = {
  1.0 * M_PI / 180.0, /* 1 degree step for phi */
  0.1,                /* hours */
  1.0,                /* days */
  1.0,                /* W/m^2 */
  0.1                 /* lunar hours */
};

/*
eefm_plot_alloc()
  Initialize an eef workspace
*/

eefm_plot_workspace *
eefm_plot_alloc(eefm_plot_parameters *params)

{
  eefm_plot_workspace *w;

  w = (eefm_plot_workspace *) calloc(1, sizeof(eefm_plot_workspace));
  if (!w)
    {
      perror("calloc");
      return (0);
    }

  w->eefm_calc_workspace_p = eefm_calc_alloc(params->model);
  if (!w->eefm_calc_workspace_p)
    {
      eefm_plot_free(w);
      return 0;
    }

  w->params = *params;

  return (w);
} /* eefm_plot_alloc() */

/*
eefm_plot_free()
  Terminate a eef workspace
*/

void
eefm_plot_free(eefm_plot_workspace *w)

{
  if (!w)
    return;

  if (w->eefm_calc_workspace_p)
    eefm_calc_free(w->eefm_calc_workspace_p);

  free(w);
} /* eefm_plot_free() */

/*
eefm_plot_proc()
  Perform fitting of mean/stddev data

Inputs: w - eef workspace
*/

void
eefm_plot_proc(eefm_plot_workspace *w)

{
  size_t i, j;
  FILE *fp;
  double param_plot[EEF_FIT_DIM];
  double param_min[EEF_FIT_DIM];
  double param_max[EEF_FIT_DIM];
  const char *param_files[] = { FILE_MODEL_PHI_OUT,
                                FILE_MODEL_T_OUT,
                                FILE_MODEL_S_OUT,
                                FILE_MODEL_EUVAC_OUT,
                                FILE_MODEL_LUNAR_OUT };
  const char *sat_names[4];
  double y;
  double E;        /* model mean value */
  double sigma;    /* model standard deviation */
  double plotval;

  /* parameter constant plot values */
  param_plot[EEF_IDX_PHI] = w->params.plot_phi;
  param_plot[EEF_IDX_T] = w->params.plot_lt;
  param_plot[EEF_IDX_S] = w->params.plot_season;
  param_plot[EEF_IDX_EUVAC] = w->params.plot_euvac;
  param_plot[EEF_IDX_LUNAR] = w->params.plot_lunar;

  /* parameter minimum values */
  param_min[EEF_IDX_PHI] = -M_PI;
  param_min[EEF_IDX_T] = EEF_BASIS_LT_MIN;
  param_min[EEF_IDX_S] = 0.0;
  param_min[EEF_IDX_EUVAC] = 73.0;
  param_min[EEF_IDX_LUNAR] = 0.0;

  /* parameter maximum values */
  param_max[EEF_IDX_PHI] = M_PI;
  param_max[EEF_IDX_T] = EEF_BASIS_LT_MAX;
  param_max[EEF_IDX_S] = 365.0;
  param_max[EEF_IDX_EUVAC] = 250.0;
  param_max[EEF_IDX_LUNAR] = 24.0 + 5.0 / 6.0;

  /* plot_idx contains the variable we want to plot */
  i = w->params.plot_idx;

  if (!w->params.filename)
    w->params.filename = param_files[i];

  fp = fopen(w->params.filename, "w");
  if (!fp)
    {
      perror("fopen");
      return;
    }

  sat_names[EEF_CALC_CHAMP] = "CHAMP";
  sat_names[EEF_CALC_OERSTED] = "Oersted";
  sat_names[EEF_CALC_SAC] = "SAC-C";

  fprintf(stderr,
          "Plotting EEF vs %s using %s model with parameters:\n",
          param_str[i],
          sat_names[w->params.model]);
  for (j = 0; j < EEF_FIT_DIM; ++j)
    {
      double param;

      if (j == i)
        continue;

      if (j == EEF_IDX_PHI)
        param = param_plot[j] * 180.0 / M_PI;
      else
        param = param_plot[j];

      fprintf(stderr, "%s (%s) = %f\n", param_str[j], param_units[j],
              param);
    }

  fprintf(stderr, "Writing %s...", w->params.filename);

  /* now print out model */

  j = 1;
  fprintf(fp, "# Field %zu: %s (%s)\n", j++,
          param_str[i], param_units[i]);
  fprintf(fp, "# Field %zu: J (A/m)\n", j++);
  fprintf(fp, "# Field %zu: sigma (A/m)\n", j++);

  fprintf(fp, "#\n# Parameter values:\n");
  for (j = 0; j < EEF_FIT_DIM; ++j)
    {
      double param;

      if (j == i)
        continue;

      if (j == EEF_IDX_PHI)
        param = param_plot[j] * 180.0 / M_PI;
      else
        param = param_plot[j];

      fprintf(fp, "# %s (%s) = %f\n", param_str[j], param_units[j],
              param);
    }

  plotval = param_plot[i]; /* save constant plot value */

  for (y = param_min[i]; y < param_max[i]; y += param_stepsize[i])
    {
      param_plot[i] = y;
      E = eefm_calc_mean(param_plot[EEF_IDX_PHI],
                         param_plot[EEF_IDX_T],
                         param_plot[EEF_IDX_S],
                         param_plot[EEF_IDX_EUVAC],
                         param_plot[EEF_IDX_LUNAR],
                         w->eefm_calc_workspace_p);

      sigma = eefm_calc_stddev(param_plot[EEF_IDX_PHI],
                               param_plot[EEF_IDX_T],
                               param_plot[EEF_IDX_S],
                               param_plot[EEF_IDX_EUVAC],
                               param_plot[EEF_IDX_LUNAR],
                               w->eefm_calc_workspace_p);

      if (i == EEF_IDX_PHI)
        fprintf(fp, "%f %f %f\n", y * 180.0 / M_PI, E, sigma);
      else
        fprintf(fp, "%f %f %f\n", y, E, sigma);
    }

  param_plot[i] = plotval; /* restore constant plot value */

  fclose(fp);

  fprintf(stderr, "done\n");
} /* eefm_plot_proc() */

int
main(int argc, char *argv[])

{
  eefm_plot_workspace *eefm_plot_workspace_p;
  eefm_plot_parameters eefm_plot_params;
  int c;

  /* defaults */
  eefm_plot_params.plot_phi = 90.0 * M_PI / 180.0;
  eefm_plot_params.plot_lt = 10.5;
  eefm_plot_params.plot_season = 266.0;
  eefm_plot_params.plot_euvac = 100.0;
  eefm_plot_params.plot_lunar = 12.5;
  eefm_plot_params.plot_idx = EEF_IDX_PHI;
  eefm_plot_params.filename = 0;
  eefm_plot_params.model = EEF_CALC_MODEL;

  while ((c = getopt(argc, argv, "l:t:s:e:m:p:o:hy:")) != (-1))
  {
    switch (c)
    {
      case 'l':
        eefm_plot_params.plot_phi = strtod(optarg, NULL) * M_PI / 180.0;
        break;

      case 't':
        eefm_plot_params.plot_lt = strtod(optarg, NULL);
        break;

      case 's':
        eefm_plot_params.plot_season = strtod(optarg, NULL);
        break;

      case 'e':
        eefm_plot_params.plot_euvac = strtod(optarg, NULL);
        break;

      case 'm':
        eefm_plot_params.plot_lunar = strtod(optarg, NULL);
        break;

      case 'p':
        {
          switch (*optarg)
            {
              case 'l':
                eefm_plot_params.plot_idx = EEF_IDX_PHI;
                break;

              case 't':
                eefm_plot_params.plot_idx = EEF_IDX_T;
                break;

              case 's':
                eefm_plot_params.plot_idx = EEF_IDX_S;
                break;

              case 'e':
                eefm_plot_params.plot_idx = EEF_IDX_EUVAC;
                break;

              case 'm':
                eefm_plot_params.plot_idx = EEF_IDX_LUNAR;
                break;

              default:
                break;
            }
          break;
        }

      case 'y':
        {
          switch (*optarg)
            {
              case 'c':
                eefm_plot_params.model = EEF_CALC_CHAMP;
                break;

              case 'o':
                eefm_plot_params.model = EEF_CALC_OERSTED;
                break;

              case 's':
                eefm_plot_params.model = EEF_CALC_SAC;
                break;

              default:
                break;
            }
          break;
        }

      case 'o':
        eefm_plot_params.filename = optarg;
        break;

      case '?':
      case 'h':
      default:
        printf("usage: %s [options]\n", argv[0]);
        printf("\n\nOptions:\n\n");
        printf("  -o <filename>   : Specifies output filename [default: plot_<var>_data]\n");
        printf("  -l <longitude>  : Specifies value of longitude parameter in degrees [default: %f]\n",
               eefm_plot_params.plot_phi * 180.0 / M_PI);
        printf("  -t <local time> : Specifies value of time parameter in hours [default: %f]\n",
               eefm_plot_params.plot_lt);
        printf("  -s <season>     : Specifies value of season parameter in day of year (0-365) [default: %f (spring equinox)]\n",
               eefm_plot_params.plot_season);
        printf("  -e <euvac>      : Specifies value of euvac in W/m^2 [default: %f]\n",
               eefm_plot_params.plot_euvac);
        printf("  -m <lunar local time> : Specifies value of lunar local time in hours [default: %f]\n",
               eefm_plot_params.plot_lunar);
        printf("  -y <c|o|s>      : Specifies which satellite model to use for plot (c = CHAMP, o = Oersted, s = SAC-C) [default: CHAMP]\n");
        printf("  -p <l|t|s|e|m>    : Specifies which variable to plot against (l = longitude, t = local time, s = season, e = euvac, m = lunar local time) [default: longitude]\n");
        exit(1);
        break;
    }
  }

  eefm_plot_workspace_p = eefm_plot_alloc(&eefm_plot_params);
  if (!eefm_plot_workspace_p)
  {
    fprintf(stderr, "main: eefm_plot_alloc failed\n");
    exit(1);
  }

  fprintf(stderr, "To change parameter values, see %s -h\n\n", argv[0]);

  eefm_plot_proc(eefm_plot_workspace_p);

  eefm_plot_free(eefm_plot_workspace_p);

  return (0);
} /* main() */
