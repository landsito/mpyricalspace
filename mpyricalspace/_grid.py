'''
Shared machinery for the empirical models: the n-d grid builder, the cwd guard the
bundled fortran needs, and the xr.Dataset packer for the "time series" / "hypercube"
cases.
'''
import os
import contextlib
import numpy as np
import xarray as xr
from datetime import timedelta


def flatten_ndgrid(add_slt=False, add_doy=False, add_uthour=False, add_year=False, **kwargs):
    '''
    Combine every input into one flat (n_points, n_params) array, plus derived columns.

    If all inputs share a length they are treated as aligned samples; otherwise a full
    meshgrid is taken. add_year/doy/uthour/slt append columns computed from 'time'
    (and 'lon' for slt). Inspired by pymsis.create_input.

    Returns (keys, arr, original_shape).
    '''
    for k in kwargs:
        kwargs[k] = np.atleast_1d(kwargs[k])

    hypercube, dimension = True, 0
    for v in kwargs.values():
        if dimension == 0:
            dimension = len(v)
        if dimension != v.size:
            hypercube = False
            break

    if "time" in kwargs:
        kwargs["time"] = np.atleast_1d(np.array(kwargs["time"], dtype='datetime64[s]'))
    if "lon" in kwargs:
        kwargs["lon"] = kwargs["lon"] % 360

    keys, _vars = zip(*kwargs.items())
    keys = list(keys)

    if not hypercube:
        _vars = np.meshgrid(*_vars, indexing='ij')

    vars = [var.astype(object) for var in _vars]           # object dtype handles datetime64
    original_shape = vars[0].shape
    arr = np.stack(vars, -1).reshape(-1, len(vars))

    it = keys.index("time") if "time" in keys else None
    if (add_year or add_doy or add_uthour or add_slt) and it is None:
        raise ValueError("add_year/add_doy/add_uthour/add_slt require 'time' among the gridded inputs")

    if add_year:
        years = arr[:, it].astype('datetime64[Y]').astype(int) + 1970.
        arr = np.concatenate((arr, years[:, None]), axis=1); keys += ["year"]
    if add_doy:
        doys = (arr[:, it].astype('datetime64[D]') - arr[:, it].astype('datetime64[Y]')).astype(float) + 1
        arr = np.concatenate((arr, doys[:, None]), axis=1); keys += ["doy"]
    if add_uthour:
        uthrs = (arr[:, it].astype('datetime64[s]') - arr[:, it].astype('datetime64[D]')).astype(float)
        arr = np.concatenate((arr, uthrs[:, None]), axis=1); keys += ["ut"]
    if add_slt:
        il = keys.index("lon")
        dt2hr = lambda dt: dt.hour + dt.minute / 60. + dt.second / 3600.
        lon = (arr[:, il] + 180) % 360 - 180
        offsets = np.array([timedelta(hours=x / 15.) for x in lon])
        slts = np.array([dt2hr(s) for s in arr[:, it] + offsets])
        arr = np.concatenate((arr, slts[:, None]), axis=1); keys += ["slt"]

    return keys, arr, original_shape


@contextlib.contextmanager
def fortran_cwd(path):
    'the bundled fortran models read their coefficient files relative to cwd'
    pwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(pwd)


def pack(vals, time, keys=None, shape=None, grid_coords=None):
    '''
    Build the output Dataset: a plain time series when `shape` is 1-D (or absent),
    otherwise an N-D cube whose dims are the first len(shape) entries of `keys`.
    '''
    if shape is None or len(shape) == 1:
        return xr.Dataset({k: ('time', np.asarray(v, float)) for k, v in vals.items()},
                          coords={'time': np.atleast_1d(time)})
    dims = list(keys[:len(shape)])
    return xr.Dataset({k: (dims, np.reshape(np.asarray(v, float), shape)) for k, v in vals.items()},
                      coords=grid_coords)
