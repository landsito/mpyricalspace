/*
 * drift_exec.c
 * Developed by L. Navarro.
 * Based on drift_plot.c developed by Patrick Alken
 *
 */
#include <Python.h>
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
	if (!w)
		return NULL;                /* the coefficient files could not be read */
	values[0] = drift_calc_mean(slt, doy,f107avg,w);
	values[1] = drift_calc_stddev(slt, doy,f107avg,w);
	drift_calc_free(w);
	return values;
}

static PyObject *jvdm1(PyObject *self, PyObject *args, PyObject *kwargs)
{
	double slt,doy,f107avg;
	static char * keywords[] = {"slt","doy","f107avg",NULL};

	double * result;
	PyObject * ret;
	PyObject * mean;
	PyObject * std;

	// parse arguments
	if (!PyArg_ParseTupleAndKeywords(args,
									kwargs,
									"ddd:jvdm1",
									keywords,
									&slt,
									&doy,
									&f107avg)) {
		return NULL;
	}

	// run the actual function
	result = getJULIAverticaldrift( slt, doy, f107avg);
	if (!result) {
		PyErr_SetString(PyExc_RuntimeError,
			"jvdm1: could not read the coefficient files drift_mean_coeffs / drift_stddev_coeffs "
			"from the current directory");
		return NULL;
	}

	// build the resulting string into a Python object.
	mean = Py_BuildValue("f", result[0]);
	std = Py_BuildValue("f",result[1]);
	ret =PyTuple_Pack(2, mean, std);

	return ret;
}

static PyMethodDef jvdm1_methods[] =
{
     {"jvdm1", (PyCFunction) jvdm1, METH_VARARGS | METH_KEYWORDS, "evaluate Alken's JULIA Vertical Drift Model"},
     {NULL, NULL, 0, NULL},
};

static struct PyModuleDef jvdm1module = {
    PyModuleDef_HEAD_INIT,
    "jvdm1",   /* name of module */
    "Module for Alken's JULIA Vertical Drift Model", /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    jvdm1_methods
};

PyMODINIT_FUNC PyInit_jvdm1(void)
{
  /* Create the module and add the functions */
  return PyModule_Create(&jvdm1module);
}
