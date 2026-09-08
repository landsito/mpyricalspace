/*
 * ppeefm1 -- Python/NumPy wrapper for the Real-Time Equatorial Electric Field
 * model (RTEEF) of Manoj, Maus & Alken -- the same model behind the Manoj & Maus
 * PPEF web service, run locally from the bundled ACE solar-wind data (2001-2007).
 * Developed by L. Navarro.  Based on rteef.c developed by Patrick Alken.
 *
 * ppeefm1.rteef(t_start, t_end, longitude [, ief]) -> (n, 3) float64 array :
 *   [ unix time [s, UT],  climatological EEF [mV/m],  prompt-penetration EEF [mV/m] ]
 * one row every 5 min over [t_start, t_end] (span <= 5 days, one calendar year).
 * Needs ace<YYYY>.dat / TF.COF / F107.txt / champ_E_*_coeffs in the cwd.
 *
 * ief (optional): 1-D float64 IEF Ey [mV/m] at 5-min cadence covering
 *   [t_start - 1h, t_end] -- length (t_end - t_start)/300 + 13.  When given, the
 *   ACE file is not read: this synthetic input is pushed through the same TF.COF
 *   prompt-penetration filter, and column 1 (climatology) is left 0.
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "rteef.h"

static PyObject *
rteef(PyObject *self, PyObject *args)
{
    long t0, t1;
    double lon;
    PyObject *ief_obj = NULL;
    PyArrayObject *ief_arr = NULL;
    const double *ief_data = NULL;
    npy_intp nief = 0;
    rteef_workspace *w = NULL;
    PyArrayObject *out = NULL;
    npy_intp n, i, dims[2];
    int s;

    if (!PyArg_ParseTuple(args, "lld|O", &t0, &t1, &lon, &ief_obj))
        return NULL;
    if (t1 < t0) {
        PyErr_SetString(PyExc_ValueError, "ppeefm1.rteef: t_end precedes t_start");
        return NULL;
    }
    if (t1 - t0 > 5L * 86400L) {
        PyErr_SetString(PyExc_ValueError, "ppeefm1.rteef: span exceeds the model's 5-day limit");
        return NULL;
    }

    n = (t1 - t0) / 300 + 1;

    if (ief_obj && ief_obj != Py_None) {
        ief_arr = (PyArrayObject *) PyArray_FROMANY(ief_obj, NPY_DOUBLE, 1, 1,
                                                    NPY_ARRAY_IN_ARRAY);
        if (!ief_arr)
            return NULL;
        nief = PyArray_SIZE(ief_arr);
        if (nief != n + 12) {
            PyErr_Format(PyExc_ValueError,
                "ppeefm1.rteef: ief must have (t_end-t_start)/300 + 13 = %ld samples "
                "(5-min cadence over [t_start-1h, t_end]), got %ld",
                (long) (n + 12), (long) nief);
            Py_DECREF(ief_arr);
            return NULL;
        }
        ief_data = (const double *) PyArray_DATA(ief_arr);
    }

    putenv("TZ=GMT");
    tzset();

    w = rteef_alloc();
    if (!w) {
        Py_XDECREF(ief_arr);
        PyErr_SetString(PyExc_RuntimeError,
                        "rteef_alloc failed (TF.COF / F107.txt / champ_E_*_coeffs not found in cwd?)");
        return NULL;
    }

    s = rteef_compute((time_t) t0, (time_t) t1, lon, w, ief_data, (size_t) nief);
    Py_XDECREF(ief_arr);
    if (s) {
        rteef_free(w);
        PyErr_SetString(PyExc_RuntimeError,
                        "rteef_compute failed (no ACE data for that year? only 2001-2007 are bundled)");
        return NULL;
    }

    dims[0] = n; dims[1] = 3;
    out = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (!out) {
        rteef_free(w);
        return NULL;
    }

    /* rteef reads 1 h (12 samples) of ACE data before t_start to prime its filter */
    for (i = 0; i < n; ++i) {
        size_t k = (size_t) i + 12;
        int ok = k < w->ndata_ief;
        *(double *) PyArray_GETPTR2(out, i, 0) = (double) (t0 + i * 300);
        *(double *) PyArray_GETPTR2(out, i, 1) = ok ? w->eef_mean[k] : NAN;
        *(double *) PyArray_GETPTR2(out, i, 2) = ok ? w->eef_rt[k]   : NAN;
    }

    rteef_free(w);
    return (PyObject *) out;
}

static PyMethodDef ppeefm1_methods[] = {
    {"rteef", (PyCFunction) rteef, METH_VARARGS,
     "RTEEF (Manoj/Maus/Alken): rteef(t_start, t_end, longitude[, ief]) -> (n,3) [ut, mean, prompt] mV/m"},
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef ppeefm1module = {
    PyModuleDef_HEAD_INIT, "ppeefm1",
    "Real-Time Equatorial Electric Field model (Manoj, Maus & Alken)", -1, ppeefm1_methods
};

PyMODINIT_FUNC
PyInit_ppeefm1(void)
{
    PyObject *m = PyModule_Create(&ppeefm1module);
    if (!m)
        return NULL;
    import_array();
    return m;
}
