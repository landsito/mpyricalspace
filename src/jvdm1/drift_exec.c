/*
 * drift_exec.c
 * Developed by L. Navarro.
 * Based on drift_plot.c developed by Patrick Alken
 *
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

double* getJULIAverticaldrift(double slt, double doy, double f107avg)
{
	static double values[2];
	drift_calc_workspace *w= drift_calc_alloc();
	values[0] = drift_calc_mean(slt, doy,f107avg,w);
	values[1] = drift_calc_stddev(slt, doy,f107avg,w);
	drift_calc_free(w);
	return values;
}

int main(int argc, char *argv[])

{
	double*values;
	double doy;

	for( doy = 1; doy < 366; doy = doy + 1 ){
		values=getJULIAverticaldrift(11.,doy,80.);
		printf("%f %f %f\n",doy,values[0],values[1]);
	}
	return (0);
} /* main() */
