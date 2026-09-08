/*
 * eefm1 -- Python/NumPy wrapper for P. Alken's Equatorial Electric Field model.
 * Developed by L. Navarro.  Based on eefm_plot.c developed by Patrick Alken.
 *
 * Ref: Alken, P., and S. Maus (2010), Electric fields in the equatorial ionosphere
 *      derived from CHAMP satellite magnetic field measurements,
 *      J. Atmos. Sol.-Terr. Phys., 72(4), 319-326, doi:10.1016/j.jastp.2009.02.006
 *
 * eefm1.eef(phi, t, s, e, tau) -> (n, 2) float64 array : [ mean E, sigma ]  (mV/m)
 *   phi : longitude          [radians]
 *   t   : local time         [hours,  7 .. 17]
 *   s   : season             [day of year, 0 .. 365]
 *   e   : solar flux (F10.7 / EUVAC proxy) [sfu]
 *   tau : lunar local time    [hours, 0 .. 24.833]
 *
 * The workspace reads champ_E_mean_coeffs / champ_E_stddev_coeffs from the current
 * working directory, so it is built once for the whole input batch (not per point).
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "eefm_basis.h"
#include "eefm_calc.h"

static PyObject *
eef(PyObject *self, PyObject *args)
{
    PyObject *a0, *a1, *a2, *a3, *a4;
    PyArrayObject *phi = NULL, *t = NULL, *s = NULL, *e = NULL, *tau = NULL, *out = NULL;
    eefm_calc_workspace *w = NULL;
    npy_intp n, i, dims[2];

    if (!PyArg_ParseTuple(args, "OOOOO", &a0, &a1, &a2, &a3, &a4))
        return NULL;

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
        PyErr_SetString(PyExc_ValueError, "eefm1.eef: phi, t, s, e, tau must have equal length");
        goto fail;
    }

    dims[0] = n; dims[1] = 2;
    out = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (!out)
        goto fail;

    w = eefm_calc_alloc(EEF_CALC_MODEL);
    if (!w) {
        PyErr_SetString(PyExc_RuntimeError,
                        "eefm_calc_alloc failed (champ_E_*_coeffs not found in cwd?)");
        goto fail;
    }

    for (i = 0; i < n; ++i) {
        double P = *(double *) PyArray_GETPTR1(phi, i);
        double T = *(double *) PyArray_GETPTR1(t, i);
        double S = *(double *) PyArray_GETPTR1(s, i);
        double E = *(double *) PyArray_GETPTR1(e, i);
        double U = *(double *) PyArray_GETPTR1(tau, i);
        *(double *) PyArray_GETPTR2(out, i, 0) = eefm_calc_mean(P, T, S, E, U, w);
        *(double *) PyArray_GETPTR2(out, i, 1) = eefm_calc_stddev(P, T, S, E, U, w);
    }

    eefm_calc_free(w);
    Py_DECREF(phi); Py_DECREF(t); Py_DECREF(s); Py_DECREF(e); Py_DECREF(tau);
    return (PyObject *) out;

fail:
    if (w) eefm_calc_free(w);
    Py_XDECREF(phi); Py_XDECREF(t); Py_XDECREF(s); Py_XDECREF(e); Py_XDECREF(tau); Py_XDECREF(out);
    return NULL;
}

static PyMethodDef eefm1_methods[] = {
    {"eef", (PyCFunction) eef, METH_VARARGS,
     "Alken Equatorial Electric Field model: eef(phi, t, s, e, tau) -> (n,2) [E, sigma]"},
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef eefm1module = {
    PyModuleDef_HEAD_INIT, "eefm1",
    "P. Alken's climatological Equatorial Electric Field model", -1, eefm1_methods
};

PyMODINIT_FUNC
PyInit_eefm1(void)
{
    PyObject *m = PyModule_Create(&eefm1module);
    if (!m)
        return NULL;
    import_array();
    return m;
}
