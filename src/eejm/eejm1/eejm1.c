/*
 * eejm1 -- Python/NumPy wrapper for P. Alken's Equatorial Electrojet model, v1.
 * Developed by L. Navarro.  Based on eej_plot.c developed by Patrick Alken.
 *
 * eejm1.eej(phi, t, s, e) -> (n, 2) float64 array : [ mean J (A/m), sigma (A/m) ]
 *   phi : longitude    [radians]
 *   t   : local time   [hours,  ~5 .. 19]
 *   s   : season       [day of year, 0 .. 365]
 *   e   : F10.7 / EUVAC solar-activity proxy [sfu]
 *
 * The workspace reads champ_mean_coeffs / champ_stddev_coeffs from the current
 * working directory, so it is built once for the whole input batch (not per point).
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "eej_basis.h"
#include "eej_calc.h"

static PyObject *
eej(PyObject *self, PyObject *args)
{
    PyObject *a0, *a1, *a2, *a3;
    PyArrayObject *phi = NULL, *t = NULL, *s = NULL, *e = NULL, *out = NULL;
    eej_calc_workspace *w = NULL;
    npy_intp n, i, dims[2];

    if (!PyArg_ParseTuple(args, "OOOO", &a0, &a1, &a2, &a3))
        return NULL;

    phi = (PyArrayObject *) PyArray_FROM_OTF(a0, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    t   = (PyArrayObject *) PyArray_FROM_OTF(a1, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    s   = (PyArrayObject *) PyArray_FROM_OTF(a2, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    e   = (PyArrayObject *) PyArray_FROM_OTF(a3, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (!phi || !t || !s || !e)
        goto fail;

    n = PyArray_SIZE(phi);
    if (PyArray_SIZE(t) != n || PyArray_SIZE(s) != n || PyArray_SIZE(e) != n) {
        PyErr_SetString(PyExc_ValueError, "eejm1.eej: phi, t, s, e must have equal length");
        goto fail;
    }

    dims[0] = n; dims[1] = 2;
    out = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (!out)
        goto fail;

    w = eej_calc_alloc(EEJ_CALC_MODEL);
    if (!w) {
        PyErr_SetString(PyExc_RuntimeError,
                        "eej_calc_alloc failed (champ_*_coeffs not found in cwd?)");
        goto fail;
    }

    for (i = 0; i < n; ++i) {
        double P = *(double *) PyArray_GETPTR1(phi, i);
        double T = *(double *) PyArray_GETPTR1(t, i);
        double S = *(double *) PyArray_GETPTR1(s, i);
        double E = *(double *) PyArray_GETPTR1(e, i);
        *(double *) PyArray_GETPTR2(out, i, 0) = eej_calc_mean(P, T, S, E, w);
        *(double *) PyArray_GETPTR2(out, i, 1) = eej_calc_stddev(P, T, S, E, w);
    }

    eej_calc_free(w);
    Py_DECREF(phi); Py_DECREF(t); Py_DECREF(s); Py_DECREF(e);
    return (PyObject *) out;

fail:
    if (w) eej_calc_free(w);
    Py_XDECREF(phi); Py_XDECREF(t); Py_XDECREF(s); Py_XDECREF(e); Py_XDECREF(out);
    return NULL;
}

static PyMethodDef eejm1_methods[] = {
    {"eej", (PyCFunction) eej, METH_VARARGS,
     "Alken Equatorial Electrojet model v1: eej(phi, t, s, e) -> (n,2) [J, sigma]"},
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef eejm1module = {
    PyModuleDef_HEAD_INIT, "eejm1",
    "P. Alken's climatological Equatorial Electrojet model, version 1", -1, eejm1_methods
};

PyMODINIT_FUNC
PyInit_eejm1(void)
{
    PyObject *m = PyModule_Create(&eejm1module);
    if (!m)
        return NULL;
    import_array();
    return m;
}
