/* eej_plot.c
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
 * ./eej_plot -h for options
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

#include "common.h"

#include "eej_basis.h"
#include "eej_calc.h"
#include "eej_files.h"
#include "eej_plot.h"

const char *param_str[] = {
  "longitude",
  "local time",
  "season",
  "euvac"
};

const char *param_units[] = {
  "degrees",
  "hours",
  "day of year",
  "W/m^2"
};

const double param_stepsize[] = {
  1.0 * M_PI / 180.0, /* 1 degree step for phi */
  0.1,                /* hours */
  1.0,                /* days */
  1.0                 /* W/m^2 */
};

/*
eej_plot_alloc()
  Initialize a eej workspace
*/

eej_plot_workspace *
eej_plot_alloc(eej_plot_parameters *params)

{
  eej_plot_workspace *w;

  w = (eej_plot_workspace *) calloc(1, sizeof(eej_plot_workspace));
  if (!w)
    {
      perror("calloc");
      return (0);
    }

  w->eej_calc_workspace_p = eej_calc_alloc(params->model);
  if (!w->eej_calc_workspace_p)
    {
      eej_plot_free(w);
      return 0;
    }

  w->params = *params;

  return (w);
} /* eej_plot_alloc() */

/*
eej_plot_free()
  Terminate a eej workspace
*/

void
eej_plot_free(eej_plot_workspace *w)

{
  if (!w)
    return;

  if (w->eej_calc_workspace_p)
    eej_calc_free(w->eej_calc_workspace_p);

  free(w);
} /* eej_plot_free() */

/*
eej_plot_proc()
  Perform fitting of mean/stddev data

Inputs: w - eej workspace
*/

void
eej_plot_proc(eej_plot_workspace *w)

{
  size_t i, j;
  FILE *fp;
  double param_plot[EEJ_FIT_DIM];
  double param_min[EEJ_FIT_DIM];
  double param_max[EEJ_FIT_DIM];
  const char *param_files[] = { FILE_MODEL_PHI_OUT,
                                FILE_MODEL_T_OUT,
                                FILE_MODEL_S_OUT,
                                FILE_MODEL_EUVAC_OUT };
  const char *sat_names[4];
  double y;
  double J;        /* model mean value */
  double sigma;    /* model standard deviation */
  double plotval;

  /* parameter constant plot values */
  param_plot[EEJ_IDX_PHI] = w->params.plot_phi;
  param_plot[EEJ_IDX_T] = w->params.plot_lt;
  param_plot[EEJ_IDX_S] = w->params.plot_season;
  param_plot[EEJ_IDX_EUVAC] = w->params.plot_euvac;

  /* parameter minimum values */
  param_min[EEJ_IDX_PHI] = -M_PI;
  param_min[EEJ_IDX_T] = 5.0;
  param_min[EEJ_IDX_S] = 0.0;
  param_min[EEJ_IDX_EUVAC] = 71.2;

  /* parameter maximum values */
  param_max[EEJ_IDX_PHI] = M_PI;
  param_max[EEJ_IDX_T] = 19.0;
  param_max[EEJ_IDX_S] = 365.0;
  param_max[EEJ_IDX_EUVAC] = 314.6;

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

  sat_names[EEJ_CALC_CHAMP] = "CHAMP";
  sat_names[EEJ_CALC_OERSTED] = "Oersted";
  sat_names[EEJ_CALC_SAC] = "SAC-C";

  fprintf(stderr,
          "Plotting EEJ vs %s using %s model with parameters:\n",
          param_str[i],
          sat_names[w->params.model]);
  for (j = 0; j < EEJ_FIT_DIM; ++j)
    {
      double param;

      if (j == i)
        continue;

      if (j == EEJ_IDX_PHI)
        param = param_plot[j] * 180.0 / M_PI;
      else
        param = param_plot[j];

      fprintf(stderr, "%s (%s) = %f\n", param_str[j], param_units[j],
              param);
    }

  fprintf(stderr, "Writing %s...", w->params.filename);

  /* now print out model */

  j = 1;
  fprintf(fp, "# Field %u: %s (%s)\n", j++,
          param_str[i], param_units[i]);
  fprintf(fp, "# Field %u: J (A/m)\n", j++);
  fprintf(fp, "# Field %u: sigma (A/m)\n", j++);

  fprintf(fp, "#\n# Parameter values:\n");
  for (j = 0; j < EEJ_FIT_DIM; ++j)
    {
      double param;

      if (j == i)
        continue;

      if (j == EEJ_IDX_PHI)
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
      J = eej_calc_mean(param_plot[EEJ_IDX_PHI],
                        param_plot[EEJ_IDX_T],
                        param_plot[EEJ_IDX_S],
                        param_plot[EEJ_IDX_EUVAC],
                        w->eej_calc_workspace_p);

      sigma = eej_calc_stddev(param_plot[EEJ_IDX_PHI],
                              param_plot[EEJ_IDX_T],
                              param_plot[EEJ_IDX_S],
                              param_plot[EEJ_IDX_EUVAC],
                              w->eej_calc_workspace_p);

      if (i == EEJ_IDX_PHI)
        fprintf(fp, "%f %f %f\n", y * 180.0 / M_PI, J, sigma);
      else
        fprintf(fp, "%f %f %f\n", y, J, sigma);
    }

  param_plot[i] = plotval; /* restore constant plot value */

  fclose(fp);

  fprintf(stderr, "done\n");
} /* eej_plot_proc() */

int
main(int argc, char *argv[])

{
  eej_plot_workspace *eej_plot_workspace_p;
  eej_plot_parameters eej_plot_params;
  int c;

  /* defaults */
  eej_plot_params.plot_phi = 90.0 * M_PI / 180.0;
  eej_plot_params.plot_lt = 10.5;
  eej_plot_params.plot_euvac = 180.0;
  eej_plot_params.plot_season = 79.0;
  eej_plot_params.plot_idx = EEJ_IDX_PHI;
  eej_plot_params.filename = 0;
  eej_plot_params.model = EEJ_CALC_MODEL;

  while ((c = getopt(argc, argv, "l:t:s:e:p:o:hy:")) != (-1))
  {
    switch (c)
    {
      case 'l':
        eej_plot_params.plot_phi = strtod(optarg, NULL) * M_PI / 180.0;
        break;

      case 't':
        eej_plot_params.plot_lt = strtod(optarg, NULL);
        break;

      case 's':
        eej_plot_params.plot_season = strtod(optarg, NULL);
        break;

      case 'e':
        eej_plot_params.plot_euvac = strtod(optarg, NULL);
        break;

      case 'p':
        {
          switch (*optarg)
            {
              case 'l':
                eej_plot_params.plot_idx = EEJ_IDX_PHI;
                break;

              case 't':
                eej_plot_params.plot_idx = EEJ_IDX_T;
                break;

              case 's':
                eej_plot_params.plot_idx = EEJ_IDX_S;
                break;

              case 'e':
                eej_plot_params.plot_idx = EEJ_IDX_EUVAC;
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
                eej_plot_params.model = EEJ_CALC_CHAMP;
                break;

              case 'o':
                eej_plot_params.model = EEJ_CALC_OERSTED;
                break;

              case 's':
                eej_plot_params.model = EEJ_CALC_SAC;
                break;

              default:
                break;
            }
          break;
        }

      case 'o':
        eej_plot_params.filename = optarg;
        break;

      case '?':
      case 'h':
      default:
        printf("usage: %s [options]\n", argv[0]);
        printf("\n\nOptions:\n\n");
        printf("  -o <filename>   : Specifies output filename [default: plot_<var>_data]\n");
        printf("  -l <longitude>  : Specifies value of longitude parameter in degrees [default: %f]\n",
               eej_plot_params.plot_phi * 180.0 / M_PI);
        printf("  -t <local time> : Specifies value of time parameter in hours [default: %f]\n",
               eej_plot_params.plot_lt);
        printf("  -s <season>     : Specifies value of season parameter in day of year (0-365) [default: %f (spring equinox)]\n",
               eej_plot_params.plot_season);
        printf("  -e <euvac>      : Specifies value of euvac in W/m^2 [default: %f]\n",
               eej_plot_params.plot_euvac);
        printf("  -y <c|o|s>      : Specifies which satellite model to use for plot (c = CHAMP, o = Oersted, s = SAC-C) [default: CHAMP]\n");
        printf("  -p <l|t|s|e>    : Specifies which variable to plot against (l = longitude, t = local time, s = season, e = euvac) [default: longitude]\n");
        exit(1);
        break;
    }
  }

  eej_plot_workspace_p = eej_plot_alloc(&eej_plot_params);
  if (!eej_plot_workspace_p)
  {
    fprintf(stderr, "main: eej_plot_alloc failed\n");
    exit(1);
  }

  fprintf(stderr, "To change parameter values, see %s -h\n\n", argv[0]);

  eej_plot_proc(eej_plot_workspace_p);

  eej_plot_free(eej_plot_workspace_p);

  return (0);
} /* main() */
