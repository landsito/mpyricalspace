/* eefm_plot.h
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

#ifndef INCLUDED_eefm_plot_h
#define INCLUDED_eefm_plot_h

#include "eefm_basis.h"
#include "eefm_calc.h"

#include <gsl/gsl_vector.h>

typedef struct
{
  double plot_phi;      /* longitude value for plot */
  double plot_lt;       /* local time value for plot */
  double plot_season;   /* season value for plot */
  double plot_euvac;    /* euvac value for plot */
  double plot_lunar;    /* local lunar time value for plot */

  size_t plot_idx;      /* index of parameter to plot */

  size_t model;         /* which satellite model to plot */

  const char *filename; /* output filename */
} eefm_plot_parameters;

typedef struct
{
  eefm_plot_parameters params;

  eefm_calc_workspace *eefm_calc_workspace_p;
} eefm_plot_workspace;

/*
 * Prototypes
 */

eefm_plot_workspace *eefm_plot_alloc(eefm_plot_parameters *params);
void eefm_plot_free(eefm_plot_workspace *w);
void eefm_plot_proc(eefm_plot_workspace *w);

#endif /* INCLUDED_eefm_plot_h */
