'''
The empirical models as pure functions: each takes time / location / (optional) indices
and returns an xr.Dataset. Missing indices are pulled from mpyricalspace.DataManager.

Empirical.run_* in Predictor.py is a thin wrapper around these.

Layout: shared internals first, then one section per model family --
  equatorial drifts | equatorial E-field & electrojet | thermospheric winds |
  geomagnetic field & ionosphere | Manoj & Maus PPEF.
Each model's private helpers sit directly above the model that uses them.
'''
import os
import shutil
import warnings
import numpy as np
import xarray as xr
from datetime import datetime, timedelta

from mpyricalspace import rocsatf as _rocsat
from mpyricalspace import scherliessfejerf as _sf
from mpyricalspace import hwm14f as _hwm14, hwm07f as _hwm07, hwm93f as _hwm93
from mpyricalspace import hltwimf as _hltwim
from mpyricalspace import jvdm1 as _jvdm
from mpyricalspace import eejm1 as _eejm1, eejm2 as _eejm2, eefm1 as _eefm1
from mpyricalspace import ppeefm1 as _ppeefm1
from mpyricalspace import igrf14f, igrf13f, igrf12f, igrf11f, igrf10f, igrf09f, iri26f, iri20f, iri16f, iri12f, iri07f, iri01f
from mpyricalspace._grid import flatten_ndgrid, fortran_cwd, pack
from mpyricalspace.DataManager import manager


# ===========================================================================
#  Shared internals
# ===========================================================================
DIR = os.path.dirname(__file__)
JRO = (-11.958256, -76.859012, 250)          # lat, lon, alt -- the only location the drift models are valid for
_CACHE = os.path.expanduser(os.environ.get('MPYRICALSPACE_DATA_DIR', '~/.cache/mpyricalspace'))

_t2dt = lambda t: np.atleast_1d(np.asarray(t, dtype='datetime64[s]')).astype('datetime64[s]').astype(datetime)
_bcast = lambda v, n: list(v) if np.ndim(v) else [v] * n


def _datadir(name):
    'coefficient-file directory for <name>: the installed package data, or the source tree (editable install)'
    for c in (os.path.join(DIR, '_data', name), os.path.join(DIR, '..', 'src', name),
              os.path.join(DIR, '..', 'src', 'hwm', name), os.path.join(DIR, '..', 'src', 'eejm', name)):
        if os.path.isdir(c):
            return os.path.abspath(c)
    raise FileNotFoundError("mpyricalspace: data files for %r not found (looked in %s/_data/ and ../src/)" % (name, DIR))


def _cwd_call(path, fn):
    'run fn() with the process cwd temporarily at `path` (several bundled models read their coeffs from cwd)'
    with fortran_cwd(path):
        return fn()


# ===========================================================================
#  Equatorial drifts -- Scherliess-Fejer, Alken JVDM, ROCSAT-1
# ===========================================================================
_SF_HISTORY_H    = 28                    # hours of AE history the storm model reads back
_SF_F107_DEFAULT = 120.                  # fallback when the store has no F10.7 for a date
_SF_AE_MAX       = 366 * 24 * 4          # the storm model's fixed AE buffer (a 15-min year)


def _fill_gaps(a):
    'forward- then leading-fill NaNs in a 1-D array (last-known value); all-NaN -> zeros'
    a = np.asarray(a, float)
    m = np.isfinite(a)
    if not m.any():
        return np.zeros_like(a)
    i = np.where(m, np.arange(a.size), -1)
    np.maximum.accumulate(i, out=i)
    i[i < 0] = int(np.argmax(m))         # leading NaNs -> first finite sample
    return a[i]


def _sf_ae_series(t0, t1, step_min, indices):
    '''AE on a regular `step_min`-minute grid over [t0, t1] (t0 step-aligned).
       indices=None -> the store; else an (n, >=2) [datetime, AE, ...] array resampled
       by last-known value. Returns (n_grid, ae_array, finite_fraction).'''
    n = int(round((t1 - t0).total_seconds() / 60. / step_min)) + 1
    grid = [t0 + timedelta(minutes=step_min * k) for k in range(n)]
    if indices is None:
        ae = manager().get_values(grid, 'ae_hourly' if step_min == 60 else 'ae_15min')
    else:
        idx = np.asarray(indices, dtype=object)
        it, ia = idx[:, 0], idx[:, 1].astype(float)
        order  = np.argsort(it)
        pos = np.clip(np.searchsorted(it[order], np.array(grid, dtype=object), side='right') - 1,
                      0, it.size - 1)
        ae = ia[order][pos]
    ae = np.asarray(ae, float)
    return n, ae, float(np.isfinite(ae).mean())


def scherliess_fejer(time, freq=None, lon=None, hour_resolution=False, F107=None,
                     indices=None, SLT=None, half_resolution=False):
    '''
    Scherliess & Fejer (1999) quiet + storm-time equatorial F-region vertical drift,
    evaluated at every datetime in `time`.

    The storm part reads 28 h of AE before the earliest time -- from the store
    (auto-fetched) or from `indices`; gaps are filled with the last-known value.

    Returns an xr.Dataset [m/s]:
        qvdrift  quiet-time drift (Scherliess-Fejer climatology)
        prompt   prompt-penetration storm drift
        dynamo   disturbance-dynamo storm drift
        dvdrift  total storm drift  (prompt + dynamo)
        vdrift   grand total        (qvdrift + dvdrift)

    lon      geographic longitude [deg] -- a scalar, or a per-time array (a track);
             defaults to Jicamarca. The quiet drift is defined around the magnetic
             equator at all longitudes; the storm part was derived at Jicamarca.
    F107     override F10.7 for every point (else per-time from the store / `indices`).
    SLT      evaluate at this fixed solar local time [h] rather than the one implied
             by `time` + `lon` (used by manoj_maus).
    indices  (n, >=2) array of [datetime, AE(, F10.7)] covering the window.
    hour_resolution  force the hourly AE model (default: 15-min, auto-falling back to
             hourly only when 15-min AE is essentially absent).

    `freq` and `half_resolution` are deprecated and ignored -- pass the exact time
    grid you want as `time`.
    '''
    if half_resolution:
        warnings.warn("scherliess_fejer: 'half_resolution' is deprecated and ignored; "
                      "pass the desired time grid as 'time'", DeprecationWarning, stacklevel=2)

    dts   = _t2dt(time)
    if dts.size == 0:
        raise ValueError("scherliess_fejer: 'time' is empty")
    lon_a = np.broadcast_to(np.asarray(JRO[1] if lon is None else lon, float), dts.size).astype(float)
    lon_m = lon_a % 360.
    d0, dn = min(dts), max(dts)

    step = 60 if hour_resolution else 15
    t0   = (d0 - timedelta(hours=_SF_HISTORY_H)).replace(second=0, microsecond=0)
    t0   = t0.replace(minute=t0.minute - t0.minute % step)
    t1   = dn.replace(second=0, microsecond=0)

    n_ae, ae, frac = _sf_ae_series(t0, t1, step, indices)
    if not hour_resolution and indices is None and frac < 0.10:
        t0h = t0.replace(minute=0)
        n60, ae60, frac60 = _sf_ae_series(t0h, t1, 60, None)
        if frac60 > frac:
            warnings.warn("15-min AE is %.0f%% complete for %s..%s; using the hourly AE model"
                          % (100 * frac, t0, t1))
            step, t0, n_ae, ae, frac = 60, t0h, n60, ae60, frac60
    if n_ae > _SF_AE_MAX:
        raise ValueError("scherliess_fejer: the time span exceeds ~1 year (the storm model's AE buffer)")
    if frac == 0.:
        warnings.warn("no AE data for %s..%s; storm drift will be 0 -- run: python -m mpyricalspace fetch"
                      % (t0, t1))
    elif frac < 0.95:
        warnings.warn("AE is %.0f%% complete for %s..%s; gaps filled with the last-known value"
                      % (100 * frac, t0, t1))

    hflag = 1 if step == 60 else 0
    aes   = np.zeros(_SF_AE_MAX)
    ae_f  = _fill_gaps(ae)
    aes[:ae_f.size] = ae_f

    if F107 is not None:
        f107s = np.full(dts.size, float(F107))
    elif (indices is not None and np.asarray(indices, dtype=object).ndim == 2
          and np.asarray(indices, dtype=object).shape[1] >= 3):
        idx   = np.asarray(indices, dtype=object)
        it    = idx[:, 0]; order = np.argsort(it)
        pos   = np.clip(np.searchsorted(it[order], np.array(dts, dtype=object), side='right') - 1,
                        0, it.size - 1)
        f107s = np.nan_to_num(idx[order, 2].astype(float)[pos], nan=_SF_F107_DEFAULT)
    else:
        f107s = np.nan_to_num(manager().get_values(dts, 'f107'), nan=_SF_F107_DEFAULT)

    min_iP = int(_SF_HISTORY_H * 60 / step) + 1
    q  = np.empty(dts.size); pr = np.empty(dts.size)
    dy = np.empty(dts.size); dv = np.empty(dts.size)
    for k, t in enumerate(dts):
        iP = int(round((t - t0).total_seconds() / 60. / step)) + 1
        if iP < min_iP:                                       # unreachable while t >= d0, kept as a guard
            raise ValueError("scherliess_fejer: need %d h of AE history before %s" % (_SF_HISTORY_H, t))
        if SLT is None:
            eer = t + timedelta(hours=lon_a[k] / 15.)
            slt = eer.hour + eer.minute / 60. + eer.second / 3600.
        else:
            slt = float(SLT)
        q[k] = _sf.vdrift_quiet_model(slt, lon_m[k], [t.timetuple().tm_yday, f107s[k]])
        pr[k], dy[k], dv[k] = _sf.vdrift_storm_model(hflag, iP, aes, slt)

    return xr.Dataset(
        data_vars={"qvdrift": ('time', q), "prompt": ('time', pr), "dynamo": ('time', dy),
                   "dvdrift": ('time', dv), "vdrift": ('time', q + pr + dy)},
        coords={'time': np.array(dts, dtype='datetime64[s]'), 'lon': ('time', lon_a)})


def jvdm1_drift(time, f107=None, f107a=None, slt=None, doy=None):
    dts = _t2dt(time)
    if (slt is None) != (doy is None):
        raise ValueError("slt and doy must be given together, or neither (then time is used)")

    if f107 is None or f107a is None:
        p = manager().get(dts, ['f107', 'f107a'])
        f107  = p[:, 1] if f107 is None else f107
        f107a = p[:, 2] if f107a is None else f107a
    f107, f107a = np.atleast_1d(f107), np.atleast_1d(f107a)

    gk = {'f107': f107, 'f107a': f107a, 'add_slt': True, 'add_doy': True}
    gk['lon'] = np.ones(len(f107)) * 283.140988 if len(f107) == len(f107a) == len(dts) else 283.140988
    if slt is not None:
        gk['slt'], gk['add_slt'] = slt, False
    if doy is not None:
        gk['doy'], gk['add_doy'] = doy, False
    if slt is None and doy is None:
        gk['time'] = dts
    keys, p, shape = flatten_ndgrid(**gk)

    favg = (p[:, keys.index("f107")] + p[:, keys.index("f107a")]) / 2.
    _p = np.concatenate((p[:, [keys.index("slt"), keys.index("doy")]], favg[:, None]), axis=1)

    day = np.logical_and(_p[:, 0] >= 8, _p[:, 0] <= 16)
    if np.count_nonzero(~day):
        warnings.warn("Alken model: %i points outside 08-16 LT set to nan" % np.count_nonzero(~day))
    DATA = np.full((_p.shape[0], 2), np.nan)
    if np.count_nonzero(day):
        DATA[day, :] = np.array([_jvdm.jvdm1(*row) for row in _p[day, :]])

    if len(shape) == 1:
        coords = {'time': p[:, keys.index("time")],
                  "lon": 283.140988,
                  "f107": ('time', p[:, keys.index("f107")]),
                  "f107a": ('time', p[:, keys.index("f107a")]),
                  "slt": ('time', p[:, keys.index("slt")])}
        data_vars = {"qvdrift_150km": ('time', DATA[:, 0]), "qvdrift_150km_error": ('time', DATA[:, 1])}
    else:
        coords = {"f107": np.atleast_1d(f107), "f107a": np.atleast_1d(f107a), "lon": np.atleast_1d([283.140988])}
        if "time" in keys:
            coords['time'] = np.atleast_1d(dts)
        if slt is not None:
            coords['slt'] = np.atleast_1d(slt)
        if doy is not None:
            coords['doy'] = np.atleast_1d(doy)
        dims = keys[:len(shape)]                     # add_slt/add_doy may have appended extra keys
        data_vars = {"qvdrift_150km": (dims, np.reshape(DATA[:, 0], shape)),
                     "qvdrift_150km_error": (dims, np.reshape(DATA[:, 1], shape))}
    return xr.Dataset(data_vars=data_vars, coords=coords).squeeze()


def rocsat_drift(time, lon=None, f107s=None, doys=None, slts=None):
    dts = _t2dt(time)
    series = doys is None and slts is None                # plain time series vs an (f107, lon, doy, slt) grid
    gk = {'add_slt': True, 'add_doy': True, 'add_uthour': False}
    if f107s is None:
        gk['f107'] = manager().get(dts, ['f107'])[:, 1]
    else:
        gk['f107'] = _bcast(f107s, dts.size) if series else f107s
    if lon is not None:
        gk['lon'] = _bcast(lon, dts.size)
    if doys is not None:
        gk['doy'], gk['add_doy'] = doys, False
    if series:
        gk['time'] = time
    if slts is not None:
        gk['slt'], gk['add_slt'] = slts, False

    keys, p, shape = flatten_ndgrid(**gk)
    idx = [keys.index(k) for k in ("lon", "doy", "f107", "slt")]
    _f  = np.asarray(p[:, idx], dtype=float)
    DATA = _rocsat.getverticaldrift_batch(_f[:, 0], _f[:, 1].astype(int), _f[:, 2], _f[:, 3])

    vals = {"qvdrift_rocsat": DATA, "slt": np.asarray(p[:, keys.index("slt")], float)}
    grid = {"f107": np.atleast_1d(gk['f107']), "lon": np.atleast_1d(lon),
            "doy": np.atleast_1d(doys) if doys is not None else None}
    grid = {k: v for k, v in grid.items() if v is not None and k in keys[:len(shape)]}
    return pack(vals, time, keys, shape, grid).squeeze()


# ===========================================================================
#  Equatorial electric field & electrojet -- P. Alken (EEF, EEJ v1 / v2)
# ===========================================================================
_NEW_MOON  = datetime(2000, 1, 6, 18, 14)     # a reference new moon (UT)
_SYNODIC   = 29.530588853                      # mean synodic month [days]
_LUNAR_DAY = 24.0 + 5.0 / 6.0                  # lunar day, 24 h 50 m

_EEJ = {1: (_eejm1, 'eejm1', (5.0, 19.0), False),
        2: (_eejm2, 'eejm2', (5.0, 19.0), True)}
_EEJ_SAT = {"champ": 0, "oersted": 1, "orsted": 1, "sac-c": 2, "sac": 2}


def _lunar_local_time(dts, slt):
    '''First-order lunar local time [h, 0..24.833] = solar local time - Moon-Sun elongation.

    The EEJ-v2 / EEF models carry the full lunar-tide dependence internally; this
    only supplies the *input* value tau. We approximate the Moon-Sun elongation as
    growing linearly at the mean synodic rate (360 deg / 29.53 d) from a reference
    new moon -- i.e. we ignore the eccentric, solar-perturbed real lunar motion.
    The true elongation leads/lags the mean by up to ~6 deg, so tau here can be off
    by ~20-25 min. That is negligible against a ~12.4 h tidal period for climatology
    but not for lunar-tide studies -- pass `lunars=` from a real ephemeris then.
    '''
    age = np.array([((d - _NEW_MOON).total_seconds() / 86400.) % _SYNODIC for d in np.atleast_1d(dts)])
    return np.mod(np.asarray(slt, float) - _LUNAR_DAY * age / _SYNODIC, _LUNAR_DAY)


def _alken_ef(time, lon, flux, slts, doys, lunars, mod, name, lt_lo, lt_hi, needs_lunar, vnames,
              call_extra=()):
    dts = _t2dt(time)
    gk = {'add_slt': slts is None, 'add_doy': doys is None}
    gk['flux'] = np.nan_to_num(manager().get_values(dts, 'f107') if flux is None
                               else np.atleast_1d(np.asarray(flux, float)), nan=120.)
    gk['lon'] = _bcast(JRO[1] if lon is None else lon, dts.size)
    if doys is not None:
        gk['doy'] = doys
    if slts is not None:
        gk['slt'] = slts
    if slts is None or doys is None:
        gk['time'] = time
    keys, p, shape = flatten_ndgrid(**gk)

    col = lambda k: np.asarray(p[:, keys.index(k)], dtype=float)
    slt, doy, flx = col('slt'), col('doy').astype(int), col('flux')
    phi = np.deg2rad((col('lon') + 180.) % 360. - 180.)

    out_lt = (slt < lt_lo) | (slt > lt_hi)
    if out_lt.any():
        warnings.warn("%s: %i point(s) outside the %g-%g h local-time window -> NaN"
                      % (name, int(out_lt.sum()), lt_lo, lt_hi))
    slt_c = np.clip(slt, lt_lo + 1e-4, lt_hi - 1e-4)   # keep the local-time bspline inside its knots

    args = [phi, slt_c, doy.astype(float), flx]
    if needs_lunar:
        if lunars is not None:
            args.append(np.asarray(_bcast(lunars, len(slt)), float))
        else:
            tt = _t2dt(p[:, keys.index('time')]) if 'time' in keys else dts
            args.append(np.broadcast_to(_lunar_local_time(tt, slt_c), slt_c.shape).copy())

    DATA = _cwd_call(_datadir(name),
                     lambda: getattr(mod, vnames[2])(
                         *(np.ascontiguousarray(a, float) for a in args), *call_extra))
    mean, sigma = DATA[:, 0].copy(), DATA[:, 1].copy()
    mean[out_lt] = sigma[out_lt] = np.nan

    grid = {"lon": np.atleast_1d(JRO[1] if lon is None else lon), 'time': np.atleast_1d(time)}
    return pack({vnames[0]: mean, vnames[1]: sigma}, time, keys, shape, grid).squeeze()


def eej(time, lon=None, flux=None, slts=None, doys=None, lunars=None, version=2, model="champ"):
    '''Equatorial Electrojet climatology (P. Alken). version = 2 (default) | 1.

    v1 is a function of longitude, local time, season (day of year) and F10.7;
    v2 replaces F10.7 with EUVAC and adds the lunar local time. Returns the
    height-integrated eastward current density and its 1-sigma [A/m], NaN outside
    the model's 5-19 h local-time window. F10.7 / EUVAC comes from the store when
    `flux` is not given.

    `model` (v2 only): which satellite the coefficients were fit from --
    "champ" (default), "oersted" or "sac-c" (Alken & Maus 2007). SAC-C flew a
    fixed local time, so its J has no local-time dependence.

    `lunars` (v2 only): if omitted, the lunar local time fed to the model is
    computed with a first-order formula (constant mean synodic rate), accurate to
    ~20-25 min -- fine for climatology, but pass your own values from a lunar
    ephemeris for lunar-tide work. See `_lunar_local_time`.
    '''
    if version not in _EEJ:
        raise ValueError("eej version must be 1 or 2 (got %s)" % version)
    m, name, (lo, hi), needs_lunar = _EEJ[version]
    extra = ()
    if version == 2:
        key = str(model).lower()
        if key not in _EEJ_SAT:
            raise ValueError("eej model must be one of %s (got %r)"
                             % (sorted(set(_EEJ_SAT)), model))
        extra = (_EEJ_SAT[key],)
    elif str(model).lower() != "champ":
        raise ValueError("eej version 1 has no satellite choice; model must be 'champ'")
    return _alken_ef(time, lon, flux, slts, doys, lunars, m, name, lo, hi, needs_lunar,
                     ('eej', 'eej_sigma', 'eej'), call_extra=extra)


def eef(time, lon=None, flux=None, slts=None, doys=None, lunars=None):
    '''Equatorial Electric Field climatology (P. Alken et al.): a function of
    longitude, local time, season, solar flux and lunar local time. Returns the
    zonal electric field and its 1-sigma [mV/m], NaN outside the 7-17 h local-time
    window. `flux` comes from the store (F10.7) when not given.

    `lunars`: if omitted, the lunar local time fed to the model is computed with a
    first-order formula (constant mean synodic rate), accurate to ~20-25 min --
    fine for climatology, but pass your own values from a lunar ephemeris for
    lunar-tide work. See `_lunar_local_time`.
    '''
    return _alken_ef(time, lon, flux, slts, doys, lunars, _eefm1, 'eefm1', 7.0, 17.0, True,
                     ('eef', 'eef_sigma', 'eef'))


# ===========================================================================
#  Thermospheric winds -- HWM 2014 / 2007 / 1993, HL-TWiM
# ===========================================================================
# HWM14 (hwm14.f90::findandopen) resolves its .bin/.dat via $HWMPATH -> no chdir needed.
# HWM07 opens its .dat files from cwd, so it uses fortran_cwd(); HWM93 has no data files.
try:
    os.environ.setdefault('HWMPATH', _datadir('hwm14'))
except FileNotFoundError:
    pass          # hwm14() itself will fail loudly if the data really is missing

_HWM = {2014: (_hwm14, 'hwm14_batch', None),
        2007: (_hwm07, 'hwm07_batch', 'hwm07'),        # 3rd item: data dir needing chdir
        1993: (_hwm93, 'gws5_batch',  None)}


def hwm(time, lat, lon, alt, ap=None, version=2014, ut=None, doy=None):
    '''
    Horizontal Wind Model. version = 2014 (default) | 2007 | 1993.
    Returns u, v [m/s] (all versions) plus du, dv (disturbance part; NaN for HWM93).
    '''
    if version not in _HWM:
        raise ValueError("HWM version must be 2014, 2007 or 1993 (got %s)" % version)
    _mod, _batch, _dd = _HWM[version]

    dts = _t2dt(time); n = dts.size
    gk = {'lat': lat, 'lon': lon, 'alt': alt, 'add_slt': False, 'add_doy': True, 'add_uthour': True}
    if ut is not None:
        gk['ut'], gk['add_uthour'] = ut, False
    if doy is not None:
        gk['doy'], gk['add_doy'] = doy, False
    per_time = doy is None and ut is None and n > 1           # else it's a single-time (grid) call
    if doy is None and ut is None:
        gk['time'] = time

    # ap / f10.7 are per-time quantities: length-n arrays for a time series, scalars for a grid
    def _idx(val, name):
        a = np.atleast_1d(np.asarray(val, float) if val is not None else manager().get_values(dts, name))
        if name[0] == 'f':
            a = np.nan_to_num(a, nan=150.)
        if not per_time:
            return float(a.reshape(-1)[0])
        return a if a.size == n else np.full(n, float(a.reshape(-1)[0]))

    gk['ap'] = _idx(ap, 'ap')
    if version != 2014:                                       # HWM07/93 also take f10.7
        gk['f107a'], gk['f107'] = _idx(None, 'f107a'), _idx(None, 'f107')
    # time series + scalar location -> broadcast so the grid stays 1-D (no meshgrid blow-up)
    if per_time and np.array(lat).size == np.array(lon).size == np.array(alt).size == 1:
        gk['lat'], gk['lon'], gk['alt'] = (np.full(n, x) for x in (lat, lon, alt))

    keys, p, shape = flatten_ndgrid(**gk)
    yy = (p[:, keys.index('time')].astype('datetime64[Y]').astype(str).astype(float) % 100
          if 'time' in keys else 0.)                          # HWM's iyd = yyddd; only ddd is used
    p[:, keys.index('doy')] = (yy * 1000 + p[:, keys.index('doy')]).astype(int)

    col = lambda *ks: np.asarray(p[:, [keys.index(k) for k in ks]], dtype=float)
    if version == 2014:
        f = col('doy', 'ut', 'alt', 'lat', 'lon', 'ap')
        call = lambda: _mod.hwm14_batch(f[:, 0].astype(int), f[:, 1], f[:, 2], f[:, 3], f[:, 4], f[:, 5],
                                        np.float32(np.nan))
    else:
        f = col('doy', 'ut', 'alt', 'lat', 'lon', 'f107a', 'f107', 'ap')
        call = lambda: getattr(_mod, _batch)(f[:, 0].astype(int), f[:, 1], f[:, 2], f[:, 3], f[:, 4],
                                             f[:, 5], f[:, 6], f[:, 7], np.float32(np.nan))
    DATA = call() if _dd is None else _cwd_call(_datadir(_dd), call)

    vals = {"u": DATA[:, 0], "v": DATA[:, 1], "du": DATA[:, 2], "dv": DATA[:, 3]}
    grid = {"ap": np.atleast_1d(ap), "lon": np.atleast_1d(lon), "lat": np.atleast_1d(lat),
            "alt": np.atleast_1d(alt), 'time': np.atleast_1d(time)}
    return pack(vals, time, keys, shape, grid).squeeze()


hwm14 = hwm          # backwards-compatible alias


def hltwim(time, lat, lon, kp=None, ut=None, doy=None):
    '''
    High-Latitude Thermospheric Wind Model (Dhadly et al., 2019) -- an NRL companion
    to HWM14 for the auroral / polar-cap winds. Height-independent (F-region).

    Valid only where |quasi-dipole magnetic latitude| > 40 deg; u, v, mu, mv are NaN
    elsewhere. Returns geographic u (zonal, + east) / v (meridional, + north), the
    same winds in quasi-dipole coordinates (mu, mv) [m/s], and the location in that
    frame: mlat (QD latitude, deg) and mlt (magnetic local time, h). kp is the
    3-hourly index (pulled from the store when not given). Give `time`, or both `ut` (hours) and
    `doy`.
    '''
    dts = _t2dt(time); n = dts.size
    gk = {'lat': lat, 'lon': lon, 'add_doy': True, 'add_uthour': True}
    if ut is not None:
        gk['ut'], gk['add_uthour'] = ut, False
    if doy is not None:
        gk['doy'], gk['add_doy'] = doy, False
    per_time = doy is None and ut is None and n > 1
    if doy is None and ut is None:
        gk['time'] = time

    a = np.atleast_1d(np.asarray(kp, float) if kp is not None else manager().get_values(dts, 'kp'))
    a = np.clip(np.nan_to_num(a, nan=3.0), 0., 10.)

    try:
        ns = int(np.broadcast_shapes(np.shape(lat), np.shape(lon))[0]) if (np.ndim(lat) or np.ndim(lon)) else 1
    except (IndexError, ValueError):
        ns = 1
    if not per_time and ns > 1:
        # many locations at one instant: broadcast the shared scalars so flatten_ndgrid
        # stays element-wise instead of taking a full meshgrid.
        gk['lat'] = np.broadcast_to(np.asarray(lat, float), ns)
        gk['lon'] = np.broadcast_to(np.asarray(lon, float), ns)
        if doy is not None: gk['doy'] = np.broadcast_to(np.asarray(doy, float), ns)
        if ut  is not None: gk['ut']  = np.broadcast_to(np.asarray(ut,  float), ns)
        gk['kp'] = np.full(ns, float(a.reshape(-1)[0]))
    elif not per_time:
        gk['kp'] = float(a.reshape(-1)[0])
    else:
        gk['kp'] = a if a.size == n else np.full(n, float(a.reshape(-1)[0]))
    if per_time and np.array(lat).size == np.array(lon).size == 1:
        gk['lat'], gk['lon'] = (np.full(n, x) for x in (lat, lon))

    keys, p, shape = flatten_ndgrid(**gk)
    col = lambda *ks: np.asarray(p[:, [keys.index(k) for k in ks]], dtype=float)
    day = col('doy')[:, 0].astype(int)
    uth = col('ut')[:, 0] * (1. if ut is not None else 1. / 3600.)     # add_uthour column is seconds
    uth = np.clip(uth, 0., 23.9999)
    lla = col('lat', 'lon', 'kp')
    DATA = _cwd_call(_datadir('hltwim'),
                     lambda: _hltwim.hltwim_batch(day, uth.astype('f4'), lla[:, 0].astype('f4'),
                                                  lla[:, 1].astype('f4'), lla[:, 2].astype('f4'),
                                                  np.float32(np.nan)))
    vals = {"u": DATA[:, 0], "v": DATA[:, 1], "mu": DATA[:, 2], "mv": DATA[:, 3],
            "mlat": DATA[:, 4], "mlt": np.mod(DATA[:, 5], 24.)}    # quasi-dipole lat, magnetic LT

    if not per_time and ns > 1:                                   # many locations at one instant
        return xr.Dataset(
            {k: (('points',), np.asarray(v, float)) for k, v in vals.items()},
            coords={'points': np.arange(ns), 'lat': ('points', col('lat')[:, 0]),
                    'lon': ('points', col('lon')[:, 0]),
                    'time': np.atleast_1d(np.asarray(time, dtype='datetime64[s]'))[0]})

    grid = {"kp": np.atleast_1d(kp), "lon": np.atleast_1d(lon), "lat": np.atleast_1d(lat),
            'time': np.atleast_1d(time)}
    return pack(vals, time, keys, shape, grid).squeeze()


# ===========================================================================
#  Geomagnetic field & ionosphere -- IGRF, NRLMSIS, IRI
# ===========================================================================
def igrf(time, lat, lon, alt, version=14):
    '''
    xr.Dataset with Bx, By, Bz, B [T], dip, dec, inc [deg]. version = 14 (default) | 13 | 12
    | 11 | 10 | 9. Source lives under src/igrf/ -- each version is one self-contained
    igrf<N>.f (coefficients baked in as DATA statements, no external files at all, unlike
    IRI); see src/igrf/sync_igrf.py (keeps it in sync with ngdc.noaa.gov) and
    src/igrf/build_pyf.py (regenerates a version's f2py .pyf, injecting the intent(in)/
    intent(out) split the bare Fortran source itself doesn't declare).
    Bx->east, By->north, Bz->up (package convention; raw IGRF is x->north, y->east, z->down).
    IGRF-14 (Nov 2024): same BGS synthesis routine as IGRF-13, one more 5-year coefficient
    epoch (2020-2025 now definitive, 2025-2030 secular-variation), valid to 2035.0.
    '''
    t64  = np.atleast_1d(np.asarray(time, dtype='datetime64[s]'))
    n    = t64.size
    lats, lons, alts = _bcast(lat, n), _bcast(lon, n), _bcast(alt, n)
    years = t64.astype('datetime64[Y]').astype(int) + 1970.
    if not (len(lons) == len(lats) == len(alts) == n):
        years = np.unique(years)

    _syn = {14: igrf14f.igrf14syn, 13: igrf13f.igrf13syn, 99: igrf13f.igrf13syn,
            12: igrf12f.igrf12syn, 11: igrf11f.igrf11syn,
            10: igrf10f.igrf10syn, 9: igrf09f.igrf9syn}
    if version not in _syn:
        raise ValueError("IGRF version must be 14, 13, 12, 11, 10 or 9 (got %s)" % version)
    syn = _syn[version]
    IGRF = lambda yr, al, la, lo: syn(0, yr, 1, al, 90 - la, lo % 360)
    if version == 99:
        import pyIGRF
        IGRF = lambda yr, al, la, lo: pyIGRF.igrf_value(la, lo % 360, al, yr)

    keys, p, shape = flatten_ndgrid(lon=lons, lat=lats, alt=alts, year=years)
    idx = [keys.index(k) for k in ("year", "alt", "lat", "lon")]
    DATA = np.array([IGRF(*row) for row in p[:, idx]])

    if version == 99:
        dec, inc, h, x, y, z, f = (DATA[:, i] for i in range(7))
        dip = np.degrees(np.arctan2(np.radians(z), 2 * np.radians(h)))
    else:
        x, y, z, f = DATA[:, 0], DATA[:, 1], DATA[:, 2], DATA[:, 3]
        h   = np.hypot(x, y)
        inc = np.degrees(np.arctan2(z, h))
        dip = np.degrees(np.arctan2(z, 2 * h))
        dec = np.degrees(np.arctan2(y, x))

    vals = {"Bx": y / 1e9, "By": x / 1e9, "Bz": -z / 1e9, "B": f / 1e9,
            "dip": dip, "dec": dec, "inc": inc}
    grid = {"lat": np.atleast_1d(lats), "lon": np.atleast_1d(lons),
            "alt": np.atleast_1d(alts), 'time': np.atleast_1d(time)}
    return pack(vals, time, keys, shape, grid)


def msis(time, lat, lon, alt, f107s=None, ap=None, nativelypackaged_indices=False, rho_m3=True):
    '''NRLMSIS 2.1 (via pymsis).

       Indices: by default (nativelypackaged_indices=False) F10.7 / F10.7a / ap come from
       DataManager, whichever of f107s / ap is not given explicitly. Pass
       nativelypackaged_indices=True to instead pull them from pymsis's own bundled/cached
       index file (pymsis.utils.get_f107_ap, source: CelesTrak) -- independent of
       DataManager; an explicit f107s / ap still overrides even then.'''
    from pymsis import msis as _msis
    dts = _t2dt(time)

    nat_f107 = nat_f107a = nat_aps = None
    if nativelypackaged_indices and (f107s is None or ap is None):
        from pymsis import utils as _msis_utils
        nat_f107, nat_f107a, nat_aps = _msis_utils.get_f107_ap(dts)

    if f107s is None:
        if nativelypackaged_indices:
            f107s, f107as = nat_f107, nat_f107a
        else:
            f107s  = manager().get_values(dts, 'f107')
            f107as = manager().get_values(dts, 'f107a')
    else:
        f107s  = f107s * np.ones(dts.size)
        f107as = f107s[:]
    bad = ~np.isfinite(f107as)
    if np.count_nonzero(bad):
        f107as[bad] = f107s[bad]

    if ap is None:
        if nativelypackaged_indices:
            aps = nat_aps
            ap  = aps[:, 1]                  # current 3-hr ap, matching the DataManager path below
        else:
            apdaily = manager().get_values(dts, 'ap_daily')
            hist = manager().get_history(dts, 'ap', [-20, 0], timedelta(hours=3))
            row1 = np.nanmean(hist[:, -12:-4], axis=1)
            row2 = np.nanmean(hist[:, -20:-12], axis=1)
            ap3  = hist[:, -4:]
            aps  = np.concatenate((row2[:, None], row1[:, None], ap3, apdaily[:, None]), axis=1)
            ap   = ap3[:, -1]
    else:
        try:    aps = np.ones((len(dts), 7)) * ap
        except Exception: aps = ap

    d = _msis.run(dts, lon, lat, alt, f107s, f107as, aps, version=2.1, geomagnetic_activity=-1)
    names = ["rho", "rho_N2", "rho_O2", "rho_O", "rho_He", "rho_H", "rho_Ar", "rho_N",
             "rho_Anomolous_O", "rho_NO", "T"]

    if d.ndim == 2:
        if not rho_m3:
            d[:, :10] = d[:, :10] * 1e-6
        return xr.Dataset({n: ('time', d[:, i]) for i, n in enumerate(names)}, coords={'time': dts}).squeeze()
    if not rho_m3:
        d[:, :, :, :, :10] = d[:, :, :, :, :10] * 1e-6
    dims = ('time', 'lon', 'lat', 'alt')
    coords = {"ap": np.atleast_1d(ap), "lon": np.atleast_1d(lon), "lat": np.atleast_1d(lat),
              "alt": np.atleast_1d(alt), 'time': np.atleast_1d(dts)}
    return xr.Dataset({n: (dims, d[:, :, :, :, i]) for i, n in enumerate(names)}, coords=coords).squeeze()


# --- IRI index files: ig_rz.dat (Rz12 / IG12) and apf107.dat (F10.7 / Ap) ---
_IRI_NAME = {2026: 'iri26', 2020: 'iri20', 2016: 'iri16', 2012: 'iri12', 2007: 'iri07', 2001: 'iri01'}


def _igrz_end_year(path):
    'ig_rz.dat line 1 = build info, line 2 (after a blank) = start_mo,start_yr,end_mo,end_yr'
    try:
        with open(path) as f:
            hdr = [ln for ln in f if ln.strip()][:2]
        return int(hdr[1].split(',')[3])
    except Exception:
        return 0


def _iri_srcdirs(name):
    '''Source directories for IRI version <name>: the installed package data (meson's
       install_data() already flattens iri/common/ + iri/<name>/ together there, so it's
       just one dir) or, for an editable/source-tree install, BOTH iri/common/ (CCIR/URSI/
       apf107/ig_rz, shared across every version) and iri/<name>/ (that version's own
       .for/.dat) -- kept as two separate directories on disk (no symlinks/copies between
       them) so iri/<name>/ only ever holds files that actually belong to that version;
       _iri_workdir() below is what flattens them, into the cache, not the source tree.'''
    installed = os.path.join(DIR, '_data', name)
    if os.path.isdir(installed):
        return [os.path.abspath(installed)]
    dirs = [os.path.join(DIR, '..', 'src', 'iri', 'common'), os.path.join(DIR, '..', 'src', 'iri', name)]
    dirs = [os.path.abspath(d) for d in dirs if os.path.isdir(d)]
    if not dirs:
        raise FileNotFoundError("mpyricalspace: data files for %r not found (looked in "
                                 "%s/_data/ and ../src/iri/{common,%s}/)" % (name, DIR, name))
    return dirs


def _iri_workdir(version, need_year=0):
    '''A writable copy of IRI-<version>'s data dir under the cache. Coefficient files are
       hard-linked from the package; ig_rz.dat / apf107.dat are real copies, refreshed from
       the CHAIN mirror when they do not reach need_year.'''
    dst = os.path.join(_CACHE, _IRI_NAME[version])
    os.makedirs(dst, exist_ok=True)
    for src in _iri_srcdirs(_IRI_NAME[version]):
        for f in os.listdir(src):
            d, s = os.path.join(dst, f), os.path.join(src, f)
            if os.path.exists(d):
                continue
            if f in ('ig_rz.dat', 'apf107.dat'):
                shutil.copy2(s, d)
            else:
                try:    os.link(s, d)
                except OSError: shutil.copy2(s, d)
    if _igrz_end_year(os.path.join(dst, 'ig_rz.dat')) < need_year:
        from mpyricalspace.DataManager import fetch_iri_indices
        try:
            fetch_iri_indices(dst)
        except Exception as e:
            warnings.warn("could not refresh IRI index files (%s); IRI may return -1 past %s" % (e, need_year))
    return dst


def refresh_iri_indices():
    'download the current ig_rz.dat / apf107.dat into every IRI working dir (used by `python -m mpyricalspace update`)'
    from mpyricalspace.DataManager import fetch_iri_indices
    dsts = [_iri_workdir(v) for v in _IRI_NAME]
    fetch_iri_indices(*dsts)
    return dsts


def iri(time, lat, lon, alt, quiet=True, F107=None, F107a=None, nativelypackaged_indices=False,
        NmF2=None, hmF2=None, version=2026,
        compute_Ne=True, compute_Te_Ti=True, compute_Ni=True, rho_m3=True):
    '''International Reference Ionosphere. version = 2026 (default) | 2020 | 2016 | 2012
       | 2007 | 2001. Source lives under src/iri/ -- see src/iri/README.md for the
       layout (common/ vs per-version), src/iri/sync_iri.py (keeps a version's .for/.dat
       set in sync with irimodel.org) and src/iri/build_pyf.py (regenerates a version's
       f2py .pyf when its .for set changes).

       Indices: by default (nativelypackaged_indices=False) F10.7 / F10.7a come from
       DataManager (same store every other model uses), whether auto-fetched or given
       as an explicit F107 / F107a (scalar or per-time array; the omitted one is then
       pulled from the store). Pass nativelypackaged_indices=True to instead have IRI
       read its own bundled apf107.dat (copied to the cache and auto-refreshed from
       the CHAIN mirror when a date runs past it, independently of DataManager) --
       an explicit F107 / F107a still overrides even then. Rz12 / IG12 are always
       read from IRI's own ig_rz.dat regardless of this flag.
       IRI-2012 caps its IGRF at 2015; IRI-2007/2001 at 2010/igrf10 respectively.

       2012/2016/2020/2026 need read_ig_rz()/readapf107() called before iri_sub (see
       the hasattr check below) -- irisub.for/irifun.for/even IRI's own iritest.for
       never call them internally; skip the call and every output is silently IRI's
       own -1 "not computed" sentinel, no error raised. 2007/2001 read both files
       inline on first use instead and don't expose these as separate subroutines
       (hasattr is False for them) -- see src/iri/sync_iri.py's VERSIONS comment for
       the full story, including a COMMON-block-size pitfall that looks identical to
       this one but isn't.'''
    keys, params, shape = flatten_ndgrid(time=time, lat=lat, lon=lon, alt=alt,
                                         add_slt=True, add_doy=True, add_year=True)
    params[:, keys.index('year')] = params[:, keys.index('year')].astype(int)
    row_dt = _t2dt(params[:, keys.index('time')])
    idx = [keys.index(k) for k in ("year", "doy", "lat", "lon", "alt", "slt")]
    params = params[:, idx]

    def get_options():
        jf, oarr = np.ones(50), np.zeros(100)
        if not compute_Ne:    jf[0] = 0
        if not compute_Te_Ti: jf[1] = 0
        if not compute_Ni:    jf[2] = 0
        jf[3], jf[30] = 0, 1            # B0,B1: ABT-2009
        jf[4]  = 1                      # foF2: CCIR
        jf[5]  = 0                      # Ni: RBV-2010 & TBT-2015
        if NmF2 is not None:
            jf[7], oarr[0] = 0, NmF2 * 100.**3
        if hmF2 is not None:
            jf[8], oarr[1] = 0, hmF2
        jf[21] = 0                      # ion densities in m-3
        jf[22] = 0                      # Te topside: newest model (TBPS-2026 / TBT-2012)
        jf[23] = 0                      # D-region FT-2001 / DRS-1995
        if want:      jf[24] = jf[31] = 0   # F10.7 daily -> oarr(41), 81-day -> oarr(46)
        if quiet:
            jf[25] = jf[32] = jf[34] = 0
            jf[35] = jf[36] = 1
        else:
            jf[25] = jf[34] = jf[32] = 1
            jf[35] = jf[36] = 0
        jf[27] = 0                      # spread-F off
        jf[28] = jf[29] = 0             # NeQuick topside
        jf[33] = 0                      # messages off
        jf[37] = 1                      # IRIFLIP writes off (<=2020) / IBP-2023 bubble model on (2026)
        jf[38] = jf[39] = 0             # hmF2: Shubin-COSMIC
        jf[40] = jf[41] = 1
        jf[46] = 0                      # CGM off
        jf[47] = 1                      # Ti: Tru-2021
        return jf, oarr

    try:
        mod = {2026: iri26f, 2020: iri20f, 2016: iri16f, 2012: iri12f,
               2007: iri07f, 2001: iri01f}[version]
    except KeyError:
        raise ValueError("IRI version must be 2026, 2020, 2016, 2012, 2007 or 2001 (got %s)" % version)

    n = params.shape[0]
    _NONE = np.full(n, -1., dtype='f4')                       # sentinel: "IRI, read your own file"

    def _series(val, name):
        a = np.asarray(val, dtype=float) if val is not None else manager().get_values(row_dt, name)
        a = np.atleast_1d(a)
        return np.nan_to_num(a if a.size == n else np.full(n, float(a.reshape(-1)[0])), nan=150.).astype('f4')

    want = (not nativelypackaged_indices) or F107 is not None or F107a is not None
    if want:
        f107x  = _series(F107, 'f107')
        f107ax = _series(F107 if F107a is None and F107 is not None else F107a, 'f107a')
    else:
        f107x = f107ax = _NONE

    jf, oarr0 = get_options()                                 # constant across the whole grid
    yr, dy, la, lo, al, slt = (np.asarray(params[:, i], dtype=t)
                               for i, t in ((0, int), (1, int), (2, float), (3, float), (4, float), (5, float)))
    wd = _iri_workdir(version, int(max(yr)))
    with fortran_cwd(wd):
        if hasattr(mod, 'read_ig_rz'):                        # IRI-2016/2020; IRI-2012 reads on first call
            mod.read_ig_rz()
            mod.readapf107()
        DATA, OARR = mod.iri_batch(jf, oarr0, yr, -dy, slt, la, lo, al,
                                   f107x, f107ax)                 # mmdd = -doy

    ne = DATA[:, 0].copy()
    if not compute_Ne:
        ne[:] = np.nan
    Tn, Ti, Te = DATA[:, 1], DATA[:, 2].copy(), DATA[:, 3].copy()
    if not compute_Te_Ti:
        Ti[:] = Te[:] = np.nan
    niO, niH, niHe, niO2, niNO = (DATA[:, i] for i in range(4, 9))
    NMF2, HMF2, NMF1, HMF1, TEC = OARR[:, 0], OARR[:, 1], OARR[:, 2], OARR[:, 3], OARR[:, 36]

    if not rho_m3:
        ne, niO, niH, niHe, niO2, niNO, NMF2, NMF1 = (a * 1e-6 for a in (ne, niO, niH, niHe, niO2, niNO, NMF2, NMF1))

    vals = {"ne": ne, "Te": Te, "Ti": Ti, "Tn": Tn, "ni_Oplus": niO, "ni_Hplus": niH,
            "ni_Heplus": niHe, "ni_O2plus": niO2, "ni_NOplus": niNO,
            "NMF2": NMF2, "HMF2": HMF2, "NMF1": NMF1, "HMF1": HMF1, "TEC": TEC}
    grid = {"lat": np.atleast_1d(lat), "lon": np.atleast_1d(lon),
            "alt": np.atleast_1d(alt), 'time': np.atleast_1d(time)}
    return pack(vals, time, keys, shape, grid)


# ===========================================================================
#  Manoj & Maus prompt-penetration EEF -- web service or local RTEEF
# ===========================================================================
_MANOJ_URL = 'https://electric-field-dot-gbp-oit-rc-hdgmrt-app-s8g.appspot.com/ppefm'
_RTEEF_YEARS = range(2001, 2008)                    # ACE solar-wind data bundled with ppeefm1


def _write_datamanager_f107txt(dst, d0, dn, pad_days=95):
    '''SPIDR-format F107.txt (the format ppeefm1's f107.c parses) built from DataManager's
       daily F10.7, covering [d0-pad_days, dn+pad_days] -- padding for f107.c's 81-day
       running average.'''
    day0 = (d0 - timedelta(days=pad_days)).date()
    dayn = (dn + timedelta(days=pad_days)).date()
    days = [day0 + timedelta(days=k) for k in range((dayn - day0).days + 1)]
    dts  = [datetime(d.year, d.month, d.day) for d in days]
    f107 = np.nan_to_num(manager().get_values(dts, 'f107'), nan=120.)
    with open(dst, 'w') as fh:
        fh.write("#F10.7 daily flux from mpyricalspace.DataManager (SPIDR-compatible format)\n")
        fh.write("#yyyy-MM-dd HH:mm value qualifier description\n")
        for d, v in zip(days, f107):
            fh.write('%s 00:00 %.1f   ""   ""\n' % (d.isoformat(), float(v)))


def _rteef_workdir(d0, dn, year=None, nativelypackaged_indices=False):
    '''A cwd with TF.COF / F107.txt / champ_E_* (+ ace<year>.dat, gunzipped from the package,
       when year is given). F107.txt: DataManager by default (regenerated per call, covering
       d0..dn), or the bundled SPIDR file when nativelypackaged_indices=True.'''
    src = _datadir('ppeefm1')
    dst = os.path.join(_CACHE, 'ppeefm1')
    os.makedirs(dst, exist_ok=True)
    for f in ('TF.COF', 'champ_E_mean_coeffs', 'champ_E_stddev_coeffs'):
        d = os.path.join(dst, f)
        if not os.path.exists(d):
            try:    os.link(os.path.join(src, f), d)
            except OSError: shutil.copy2(os.path.join(src, f), d)

    f107_dst = os.path.join(dst, 'F107.txt')
    if nativelypackaged_indices:
        shutil.copy2(os.path.join(src, 'F107.txt'), f107_dst)
    else:
        _write_datamanager_f107txt(f107_dst, d0, dn)

    if year is None:                              # synthetic-IEF path: no ACE file needed
        return dst
    ace = os.path.join(dst, 'ace%d.dat' % year)
    if not os.path.exists(ace):
        plain, gz = os.path.join(src, 'ace%d.dat' % year), os.path.join(src, 'ace%d.dat.gz' % year)
        if os.path.exists(plain):
            try:    os.link(plain, ace)
            except OSError: shutil.copy2(plain, ace)
        elif os.path.exists(gz):
            import gzip
            with gzip.open(gz, 'rb') as fi, open(ace + '.tmp', 'wb') as fo:
                shutil.copyfileobj(fi, fo)
            os.replace(ace + '.tmp', ace)
        else:
            raise FileNotFoundError("ppeefm1: no ACE data for %d (bundled range is 2001-2007)" % year)
    return dst


def _rteef(d0, dn, lon=-77., ief=None, nativelypackaged_indices=False):
    '''Local RTEEF (Manoj/Maus/Alken) -> (datetime, qef, tef, ppef) array, same layout
       as _manoj_web. qef = climatology, ppef = prompt penetration, tef = qef + ppef.

       ief: optional synthetic IEF Ey [mV/m], 5-min cadence, length
       (dn - d0)/300 s + 13 (one hour of priming before d0, through dn). When
       given, the bundled ACE data is not read and qef is 0 -- a pure
       transfer-function evaluation; the year is then unrestricted.

       nativelypackaged_indices=False (default): the F10.7 the climatology (qef) is built
       from comes from DataManager. =True: from the bundled SPIDR F107.txt.'''
    import calendar
    t0, t1 = calendar.timegm(d0.timetuple()), calendar.timegm(dn.timetuple())
    if ief is not None:
        ief = np.ascontiguousarray(ief, dtype=float)
        wd = _rteef_workdir(d0, dn, nativelypackaged_indices=nativelypackaged_indices)
        DATA = _cwd_call(wd, lambda: _ppeefm1.rteef(int(t0), int(t1), float(lon), ief))
    else:
        if d0.year not in _RTEEF_YEARS or dn.year not in _RTEEF_YEARS:
            raise ValueError("ppeefm1 has ACE data only for 2001-2007 (got %s .. %s)" % (d0.date(), dn.date()))
        if d0.year != dn.year:
            raise ValueError("ppeefm1: start and end must fall in the same calendar year")
        wd = _rteef_workdir(d0, dn, d0.year, nativelypackaged_indices=nativelypackaged_indices)
        DATA = _cwd_call(wd, lambda: _ppeefm1.rteef(int(t0), int(t1), float(lon)))
    ts   = np.array([datetime(1970, 1, 1) + timedelta(seconds=int(x)) for x in DATA[:, 0]], dtype=object)
    mean, rt = DATA[:, 1], DATA[:, 2]
    return np.column_stack((ts, mean, mean + rt, rt))          # dt, qef, tef, ppef


def _manoj_web(d0, dn, lon=-77.):
    'call the Manoj & Maus PPEF web model -> (datetime, qef, tef, ppef) array'
    import re
    import requests
    ndays = max((dn.date() - d0.date()).days, 1)
    p = {'year': '%04i' % d0.year, 'month': '%02i' % d0.month, 'day': '%02i' % d0.day,
         'utc': '%02i' % d0.hour, 'nDays': str(ndays), 'long': str(int(lon)), 'sat': 'auto', 'download': 'false'}
    txt = requests.get(_MANOJ_URL, params=p, timeout=120).text
    rows = re.findall(r'\[\s*new Date\s*\((.*?)\)\s*,\s*([\d.+\-eE]+)\s*,\s*([\d.+\-eE]+)\s*,\s*([\d.+\-eE]+)\s*\]', txt)
    out = []
    for expr, ppef, qef, tef in rows:                         # site column order: prompt, quiet, total (mV/m)
        expr = re.sub(r'new Date\([^)]+\)\.getTimezoneOffset\(\)', '0', expr).strip()
        if not re.fullmatch(r'[\d\s.+\-*/eE()]+', expr):      # arithmetic only, no names/calls
            continue
        ms = eval(expr, {"__builtins__": {}}, {})
        out.append([datetime.utcfromtimestamp(ms / 1000), float(qef), float(tef), float(ppef)])
    if not out:
        raise RuntimeError("the Manoj PPEF service returned no rows for %s..%s" % (d0, dn))
    return np.array(out)


def manoj_maus(time, d0, dn, freq, lon=-77., data=None, nativelypackaged_code=False, ief=None,
               nativelypackaged_indices=False):
    '''
    Manoj & Maus prompt-penetration / quiet equatorial electric field expressed as a
    vertical drift, by scaling the model onto the Scherliess-Fejer quiet drift.

    nativelypackaged_code=False (default) fetches the live Manoj & Maus PPEF web service.
    nativelypackaged_code=True runs the same model (RTEEF, the bundled ppeefm1 C code)
    locally from the ACE data bundled with the package -- no network, but valid only
    for 2001-2007, within one calendar year, and spans of at most 5 days. Pass a
    (datetime, qef, tef, ppef) array as `data` to use your own.

    ief: a synthetic IEF Ey time series [mV/m] (5-min cadence, length
    (dn - d0)/300 s + 13, covering one hour before d0 through dn). It overrides the
    model's internal ACE solar-wind buffer, so the *same* TF.COF prompt-penetration
    filter runs on your input. In this mode the return is the raw model output in
    mV/m -- ``qef`` (climatology; 0 here), ``ppef`` (prompt penetration), ``tef``
    (total) -- not the drift-scaled fields. It does not affect ``ppef``.

    nativelypackaged_indices (nativelypackaged_code=True / ief modes only): the F10.7
    the climatology (``qef``) is built from -- False (default) uses DataManager, True
    the bundled SPIDR ``F107.txt``. Only sets the climatological baseline either way;
    it does not affect ``ppef``. No effect when nativelypackaged_code=False (the web
    service has its own).
    '''
    if ief is not None:
        raw  = _rteef(d0, dn, lon, ief=np.asarray(ief, float),
                      nativelypackaged_indices=nativelypackaged_indices)
        dts  = _t2dt(time)
        hrs  = np.array([(x - d0).total_seconds() / 60. for x in dts])
        dhrs = np.array([(x - d0).total_seconds() / 60. for x in raw[:, 0]])
        return xr.Dataset(
            data_vars={"qef":  ('time', np.interp(hrs, dhrs, raw[:, 1].astype(float))),
                       "tef":  ('time', np.interp(hrs, dhrs, raw[:, 2].astype(float))),
                       "ppef": ('time', np.interp(hrs, dhrs, raw[:, 3].astype(float)))},
            coords={'time': dts})

    if data is None:
        data = (_rteef(d0, dn, lon, nativelypackaged_indices=nativelypackaged_indices)
                if nativelypackaged_code else _manoj_web(d0, dn, lon))
    else:
        data = np.asarray(data, dtype=object)
    qref = scherliess_fejer(time, lon=lon).qvdrift.data      # aligned 1:1 with `time`
    dts  = _t2dt(time)
    hrs  = np.array([(d - d0).total_seconds() / 60. for d in dts])
    dhrs = np.array([(d - d0).total_seconds() / 60. for d in data[:, 0]])
    qef  = np.interp(hrs, dhrs, data[:, 1].astype(float))
    tef  = np.interp(hrs, dhrs, data[:, 2].astype(float))
    ppef = np.interp(hrs, dhrs, data[:, 3].astype(float))

    factor = qref / qef
    bad = np.abs(factor) >= 70
    if np.count_nonzero(bad):
        factor = np.interp(hrs.astype(float), hrs[~bad].astype(float), factor[~bad].astype(float))

    return xr.Dataset(
        data_vars={"mV_to_ms": ('time', factor), "vdrift_manoj": ('time', factor * tef),
                   "prompt_manoj": ('time', factor * ppef), "qvdrift": ('time', qref)},
        coords={'time': dts})
