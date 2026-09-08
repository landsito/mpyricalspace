/*
 * eejm2 -- Python/NumPy wrapper for P. Alken's Equatorial Electrojet model, v2.
 * Developed by L. Navarro.  Based on eej_plot.c developed by Patrick Alken.
 *
 * Ref: Alken, P., and S. Maus (2007), Spatio-temporal characterization of the
 *      equatorial electrojet from CHAMP, Orsted, and SAC-C satellite magnetic
 *      measurements, J. Geophys. Res., 112, A09305, doi:10.1029/2007JA012524
 *
 * eejm2.eej(phi, t, s, e, tau [, model]) -> (n, 2) float64 : [ mean J (A/m), sigma (A/m) ]
 *   phi : longitude          [radians]
 *   t   : local time         [hours,  ~5 .. 19]
 *   s   : season             [day of year, 0 .. 365]
 *   e   : EUVAC / F10.7 proxy [sfu]
 *   tau : lunar local time    [hours, 0 .. 24.833]
 *   model : 0 = CHAMP (default), 1 = Orsted, 2 = SAC-C  -- which satellite the
 *           coefficients were fit from. SAC-C flew a fixed local time, so its J
 *           has no local-time dependence.
 *
 * The workspace reads <sat>_mean_coeffs / <sat>_stddev_coeffs from the current
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
    PyObject *a0, *a1, *a2, *a3, *a4;
    PyArrayObject *phi = NULL, *t = NULL, *s = NULL, *e = NULL, *tau = NULL, *out = NULL;
    eej_calc_workspace *w = NULL;
    npy_intp n, i, dims[2];
    int model = EEJ_CALC_CHAMP;

    if (!PyArg_ParseTuple(args, "OOOOO|i", &a0, &a1, &a2, &a3, &a4, &model))
        return NULL;
    if (model != EEJ_CALC_CHAMP && model != EEJ_CALC_OERSTED && model != EEJ_CALC_SAC) {
        PyErr_SetString(PyExc_ValueError,
                        "eejm2.eej: model must be 0 (CHAMP), 1 (Orsted) or 2 (SAC-C)");
        return NULL;
    }

    phi = (PyArrayObject *) PyArray_FROM_OTF(a0, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    t   = (PyArrayObject *) PyArray_FROM_OTF(a1, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    s   = (PyArrayObject *) PyArray_FROM_OTF(a2, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    e   = (PyArrayObject *) PyArray_FROM_OTF(a3, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    tau = (PyArrayObject *) PyArray_FROM_OTF(a4, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (!phi || !t || !s || !e || !tau)
        goto fail;

    n = PyArray_SIZE(phi);
    if (PyArray_SIZE(t) != n || PyArray_SIZE(s) != n ||
        PyArray_SIZE(e) != n || PyArray_SIZE(tau) != n) {
        PyErr_SetString(PyExc_ValueError, "eejm2.eej: phi, t, s, e, tau must have equal length");
        goto fail;
    }

    dims[0] = n; dims[1] = 2;
    out = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (!out)
        goto fail;

    w = eej_calc_alloc((size_t) model);
    if (!w) {
        PyErr_SetString(PyExc_RuntimeError,
                        "eej_calc_alloc failed (<sat>_*_coeffs not found in cwd?)");
        goto fail;
    }

    for (i = 0; i < n; ++i) {
        double P = *(double *) PyArray_GETPTR1(phi, i);
        double T = *(double *) PyArray_GETPTR1(t, i);
        double S = *(double *) PyArray_GETPTR1(s, i);
        double E = *(double *) PyArray_GETPTR1(e, i);
        double U = *(double *) PyArray_GETPTR1(tau, i);
        *(double *) PyArray_GETPTR2(out, i, 0) = eej_calc_mean(P, T, S, E, U, w);
        *(double *) PyArray_GETPTR2(out, i, 1) = eej_calc_stddev(P, T, S, E, U, w);
    }

    eej_calc_free(w);
    Py_DECREF(phi); Py_DECREF(t); Py_DECREF(s); Py_DECREF(e); Py_DECREF(tau);
    return (PyObject *) out;

fail:
    if (w) eej_calc_free(w);
    Py_XDECREF(phi); Py_XDECREF(t); Py_XDECREF(s); Py_XDECREF(e); Py_XDECREF(tau); Py_XDECREF(out);
    return NULL;
}

static PyMethodDef eejm2_methods[] = {
    {"eej", (PyCFunction) eej, METH_VARARGS,
     "Alken Equatorial Electrojet model v2: eej(phi, t, s, e, tau[, model]) -> (n,2) [J, sigma]"},
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef eejm2module = {
    PyModuleDef_HEAD_INIT, "eejm2",
    "P. Alken's climatological Equatorial Electrojet model, version 2", -1, eejm2_methods
};

PyMODINIT_FUNC
PyInit_eejm2(void)
{
    PyObject *m = PyModule_Create(&eejm2module);
    if (!m)
        return NULL;
    import_array();
    return m;
}
